"""Compare paired HF/CF results and replay the generation baseline on each saved pool."""
from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
import json
from pathlib import Path

from rdkit import Chem


def load_predictions(path, modality):
    run = json.loads(path.with_suffix('.run.json').read_text(encoding='utf-8'))
    if run['split'] != 'test':
        raise ValueError('Expected TEST predictions after weights were fixed on VALID')
    rows = {}
    with path.open(encoding='utf-8') as stream:
        for line in stream:
            row = json.loads(line)
            if row['id'] in rows or row['split'] != 'test' or row['evaluation_metadata']['effective_input_mode'] != modality:
                raise ValueError(f"Duplicate ID, wrong split or wrong modality: {row['id']}")
            rows[row['id']] = row
    if not rows:
        raise ValueError(f'No predictions found: {path}')
    return run, rows


@lru_cache(maxsize=32768)
def skeleton(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f'Invalid SMILES in the strict candidate pool: {smiles}')
    Chem.RemoveStereochemistry(molecule)
    return Chem.MolToSmiles(molecule, isomericSmiles=True)


def outcomes(row):
    target = row['target_for_evaluation_only']
    target_skeleton = skeleton(target)
    diag = row['generation_diagnostics']
    ordered = sorted(diag['candidate_pool'], key=lambda c: (-c['generation_score'], c['canonical_smiles']))
    pool = [c['canonical_smiles'] for c in ordered]
    ranked = [c['canonical_smiles'] for c in row['beams']]
    if set(ranked) - set(pool):
        raise ValueError(f"Reranked candidates differ from the saved pool: {row['id']}")
    result = {'pool_exact': target in pool, 'pool_topology': any(skeleton(s) == target_skeleton for s in pool),
              'candidate_shortfall': bool(diag['shortfall'])}
    for name, order in [('generation', pool[:10]), ('reranked', ranked[:10])]:
        result[f'{name}_top10_topology'] = any(skeleton(s) == target_skeleton for s in order)
        for k in (1, 3, 10):
            result[f'{name}_top{k}_exact'] = target in order[:k]
    return result


def compare(hf_path, cf_path):
    hf_run, hf = load_predictions(hf_path, 'proton_only')
    cf_run, cf = load_predictions(cf_path, 'carbon_only')
    if hf_run['generation_policy'] != cf_run['generation_policy'] or hf.keys() != cf.keys():
        raise ValueError('HF and CF must use the same decoding policy and identical TEST IDs')
    bins = defaultdict(list)
    for mid in hf:
        h, c = hf[mid], cf[mid]
        truth = h['target_for_evaluation_only']
        hmeta, cmeta = h['evaluation_metadata'], c['evaluation_metadata']
        if truth != c['target_for_evaluation_only'] or any(hmeta[key] != cmeta[key] for key in ('input_molecular_formula', 'size_bin')):
            raise ValueError(f'Test target or metadata differs: {mid}')
        pair = (outcomes(h), outcomes(c))
        bins['overall'].append(pair)
        bins[hmeta['size_bin']].append(pair)
    report = {'hf_run': hf_run, 'cf_run': cf_run, 'groups': {}}
    for label, pairs in bins.items():
        group = {'samples': len(pairs)}
        for metric in pairs[0][0]:
            group[metric] = {'hf': sum(h[metric] for h, _ in pairs), 'cf': sum(c[metric] for _, c in pairs),
                             'hf_only': sum(h[metric] and not c[metric] for h, c in pairs),
                             'cf_only': sum(c[metric] and not h[metric] for h, c in pairs)}
        report['groups'][label] = group
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hf', required=True, type=Path)
    parser.add_argument('--cf', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Choose a fresh output path')
    result = compare(args.hf, args.cf)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2)
    print(json.dumps(result['groups']['overall'], indent=2))


if __name__ == '__main__':
    main()
