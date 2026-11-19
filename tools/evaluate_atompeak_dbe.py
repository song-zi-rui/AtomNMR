"""Compare AtomPeak decoding with and without the DBE budget on a frozen split."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import os
from pathlib import Path
import random
import sys
import time

os.environ['NMRPEAK_LIGHT_IMPORT'] = '1'
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from rdkit import Chem, RDLogger, rdBase
import torch

from nmrpeak.atompeak.candidate_policy_partial import heavy_atom_budget
from nmrpeak.atompeak.data_partial import PartialDataset
from nmrpeak.atompeak.formula_constraints_partial import formula_dbe, make_dbe_processor
from nmrpeak.atompeak.generation_metrics_partial import input_mode, size_bucket
from nmrpeak.atompeak.infer_constrained_partial import file_sha256, infer_constrained
from nmrpeak.atompeak.model_partial import PartialAtomPeakModel


def select_samples(dataset, limit, seed):
    strata = defaultdict(list)
    excluded = Counter()
    for index in range(len(dataset)):
        row = dataset.raw(index)
        formula = row['molecular_formula']
        mode = input_mode(row['branch_present'])
        if formula_dbe(formula) is None:
            excluded['unsupported_formula'] += 1
        elif '.' in row['smiles_canonical']:
            excluded['disconnected_reference'] += 1
        elif mode == 'no_input':
            excluded['missing_spectrum'] += 1
        else:
            size = size_bucket(sum(heavy_atom_budget(formula).values()))
            strata[size, mode].append((index, row['id']))

    eligible = sum(map(len, strata.values()))
    if not eligible:
        raise ValueError('No eligible records for DBE evaluation')
    selected_count = min(limit, eligible) if limit else eligible
    # Largest-remainder allocation keeps the split's eligible size/mode proportions.
    quotas = {key: divmod(len(rows) * selected_count, eligible) for key, rows in strata.items()}
    order = sorted(quotas, key=lambda key: (-quotas[key][1], key))
    extra = set(order[:selected_count - sum(q[0] for q in quotas.values())])
    rng = random.Random(seed)
    selected, counts = [], []
    for size, mode in sorted(strata):
        key = (size, mode)
        rows = strata[key]
        count = quotas[key][0] + int(key in extra)
        counts.append({'size_bin': size, 'input_mode': mode, 'eligible': len(rows), 'selected': count})
        sample = rng.sample(rows, count)
        selected.extend((index, mid, size, mode) for index, mid in sample)
    rng.shuffle(selected)
    return selected, {'eligible': eligible, 'excluded': dict(excluded), 'strata': counts}


def reference_blocked(record, dataset):
    vocab = dataset.target_vocab
    ids = torch.tensor([vocab.encode_smiles(record['smiles_canonical'])])
    scores = torch.zeros(1, len(vocab))
    processor = make_dbe_processor(vocab, record['molecular_formula'], True)
    return any(processor(ids[:, :i], scores)[0, ids[0, i]].isneginf().item()
               for i in range(1, ids.size(1)))


def topology(smiles):
    mol = Chem.MolFromSmiles(smiles)
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def measure_arm(generate, record, dbe, device):
    if device.type == 'cuda':
        torch.cuda.synchronize(device)
    start = time.perf_counter()
    beams, diagnostics = generate(record, dbe)
    if device.type == 'cuda':
        torch.cuda.synchronize(device)
    seconds = time.perf_counter() - start
    if (diagnostics['prefix_dbe_budget'] is not None) != dbe:
        raise ValueError('This comparison requires both arms to honor the requested DBE setting')

    truth = Chem.CanonSmiles(record['smiles_canonical'])
    skeleton = topology(truth)
    top = [candidate['canonical_smiles'] for candidate in beams]
    pool = [candidate['canonical_smiles'] for candidate in diagnostics['candidate_pool']]
    skeletons = {smiles: topology(smiles) for smiles in pool}
    return {'seconds': round(seconds, 6), 'top1_exact': truth in top[:1], 'top10_exact': truth in top,
            'pool_exact': truth in pool, 'top10_topology': any(skeletons[s] == skeleton for s in top),
            'pool_topology': skeleton in skeletons.values(), 'shortfall': bool(diagnostics['shortfall']),
            'empty_pool': not pool, 'pool_size': len(pool), 'top_candidates': top,
            'raw_returns': sum(a['raw_returned'] for a in diagnostics['attempts']),
            'masked_returns': sum(a.get('masked_returns', 0) for a in diagnostics['attempts']),
            'rejections': diagnostics['rejection_counts']}


def summarize_pairs(pairs, seed):
    result = {'samples': len(pairs), 'blocked_reference_sequences': sum(p['reference_blocked'] for p in pairs)}
    rng = np.random.default_rng(seed)
    # Resample whole connectivity groups to keep stereoisomer observations together.
    clusters = defaultdict(list)
    for pair in pairs:
        clusters[pair['connectivity_key']].append(pair)
    groups = list(clusters.values())
    draws = rng.integers(len(groups), size=(1000, len(groups))) if len(groups) > 1 else None
    group_sizes = np.array([len(group) for group in groups])
    resampled_size = group_sizes[draws].sum(axis=1) if draws is not None else None
    for metric in ('top1_exact', 'top10_exact', 'pool_exact', 'top10_topology', 'pool_topology', 'shortfall', 'empty_pool'):
        off = sum(p['off'][metric] for p in pairs)
        on = sum(p['on'][metric] for p in pairs)
        delta = (on - off) / len(pairs)
        interval = None
        if draws is not None:
            differences = np.array([sum(int(p['on'][metric]) - int(p['off'][metric]) for p in group) for group in groups])
            # A bootstrap cannot estimate unobserved switch effects with zero discordant pairs.
            if any(p['on'][metric] != p['off'][metric] for p in pairs):
                interval = np.quantile(differences[draws].sum(axis=1) / resampled_size, [.025, .975]).tolist()
        result[metric] = {'off': off, 'on': on, 'delta_rate': delta, 'paired_bootstrap_95': interval,
                          'on_only': sum(p['on'][metric] and not p['off'][metric] for p in pairs),
                          'off_only': sum(p['off'][metric] and not p['on'][metric] for p in pairs)}
    for mode in ('off', 'on'):
        arms = [p[mode] for p in pairs]
        raw = sum(arm['raw_returns'] for arm in arms)
        rejections = Counter()
        for arm in arms:
            rejections.update(arm['rejections'])
        result[mode] = {'mean_pool_size': sum(arm['pool_size'] for arm in arms) / len(arms),
                        'median_seconds': float(np.median([arm['seconds'] for arm in arms])),
                        'total_seconds': sum(arm['seconds'] for arm in arms), 'raw_returns': raw,
                        'masked_returns': sum(arm['masked_returns'] for arm in arms),
                        'rejections': dict(rejections), 'rejection_rates': {k: v / raw for k, v in rejections.items()}}
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', required=True, type=Path)
    parser.add_argument('--checkpoint', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--split', choices=('valid', 'test'), default='valid')
    parser.add_argument('--samples', type=int, default=300, help='Use 0 for all DBE-eligible records.')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--beams', type=int, default=100)
    parser.add_argument('--max-length', type=int)
    parser.add_argument('--device', default='cuda' if torch.cuda.is_available() else 'cpu')
    args = parser.parse_args(argv)
    if args.samples < 0 or args.beams < 10 or (args.max_length is not None and args.max_length < 2):
        parser.error('Require samples >= 0, beams >= 10, and max-length >= 2')
    if args.output_dir.exists():
        parser.error('Choose a new output directory')
    device = torch.device(args.device)
    RDLogger.DisableLog('rdApp.error')
    RDLogger.DisableLog('rdApp.warning')
    dataset = PartialDataset(args.data, args.split)
    try:
        selected, selection = select_samples(dataset, args.samples, args.seed)
        checkpoint = torch.load(args.checkpoint, map_location='cpu', weights_only=False)
        if (checkpoint['source_tokens'] != dataset.source_vocab.tokens
                or checkpoint['target_tokens'] != dataset.target_vocab.tokens):
            raise ValueError('Checkpoint/export vocabularies differ')
        tokenization = dataset.manifest['tokenization']
        fingerprint = dataset.manifest['source_snapshot_fingerprint']
        if checkpoint['tokenization'] != tokenization or checkpoint['dataset_fingerprint'] != fingerprint:
            raise ValueError('Checkpoint/export tokenization or snapshot differs')
        max_length = args.max_length or tokenization['max_target_len']
        model = PartialAtomPeakModel(**checkpoint['model_config']).to(device)
        model.load_state_dict(checkpoint['model'])
        model.eval()
        del checkpoint

        def generate(record, dbe):
            return infer_constrained(model, record, dataset.source_vocab, dataset.target_vocab, tokenization, device,
                                     num_candidates=10, search_beams=args.beams, max_search_beams=args.beams,
                                     max_length=max_length, dbe_budget=dbe)

        run = {'split': args.split, 'seed': args.seed, 'sample_ids': [mid for _, mid, _, _ in selected],
               'checkpoint': str(args.checkpoint.resolve()), 'checkpoint_sha256': file_sha256(args.checkpoint),
               'dataset_fingerprint': fingerprint, 'data': str(args.data.resolve()), 'fixed_beams': args.beams,
               'max_length': max_length, 'device': str(device), 'torch': str(torch.__version__),
               'rdkit': rdBase.rdkitVersion, 'selection': selection}
        args.output_dir.mkdir(parents=True, exist_ok=False)
        (args.output_dir / 'run.json').write_text(json.dumps(run, indent=2), encoding='utf-8')
        print(json.dumps({'selected': len(selected), 'eligible': selection['eligible'], 'beams': args.beams}), flush=True)
        if device.type == 'cuda':
            warmup = dataset.raw(selected[0][0])
            generate(warmup, False)
            generate(warmup, True)

        pairs = []
        with (args.output_dir / 'pairs.jsonl').open('x', encoding='utf-8') as handle:
            for i, (index, mid, size, mode) in enumerate(selected):
                record = dataset.raw(index)
                pair = {'id': mid, 'connectivity_key': record.get('connectivity_key') or mid,
                        'formula': record['molecular_formula'], 'target': record['smiles_canonical'],
                        'size_bin': size, 'input_mode': mode, 'reference_blocked': reference_blocked(record, dataset)}
                order = ('off', 'on') if i % 2 == 0 else ('on', 'off')
                for setting in order:
                    pair[setting] = measure_arm(generate, record, setting == 'on', device)
                pairs.append(pair)
                handle.write(json.dumps(pair) + '\n')
                handle.flush()
                if len(pairs) % 10 == 0:
                    print(json.dumps({'paired': len(pairs), 'total': len(selected)}), flush=True)

        report = {'run': run, 'overall': summarize_pairs(pairs, args.seed),
                  'by_size': {}, 'by_input_mode': {}}
        for field, group in (('size_bin', 'by_size'), ('input_mode', 'by_input_mode')):
            for value in sorted({p[field] for p in pairs}):
                report[group][value] = summarize_pairs([p for p in pairs if p[field] == value], args.seed)
        (args.output_dir / 'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        overall = report['overall']
        print(json.dumps({'samples': overall['samples'], 'top10_exact': overall['top10_exact'],
                          'pool_exact': overall['pool_exact'], 'blocked_reference_sequences': overall['blocked_reference_sequences'],
                          'report': str(args.output_dir / 'report.json')}), flush=True)
    finally:
        dataset.close()


if __name__ == '__main__':
    main()
