"""Create matched HF/CF data views from an existing frozen AtomPeak export."""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import ExitStack
import hashlib
import json
import pickle
from pathlib import Path
import shutil
import sys

import lmdb

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from atompeak_data.schema import BRANCHES, NUCLEI, SCHEMA_VERSION, branch_supervision


ARMS = {'hf': 'proton', 'cf': 'carbon', 'hf_shift': 'proton'}


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding='utf-8')


def view_record(blob, keep):
    record = pickle.loads(blob)
    hidden = next(branch for branch in BRANCHES if branch != keep)
    record['spectra'][hidden] = {'peaks': [], 'metadata': {}}
    record['h_nmr_peaks' if hidden == 'proton' else 'c_nmr_peaks'] = []
    record['branch_present'][hidden] = False
    record['conditions']['by_nucleus'][NUCLEI[hidden]] = {}
    record['provenance']['selected_spectrum_ids'][hidden] = None
    group_key = 'proton_environments' if hidden == 'proton' else 'carbon_groups'
    _, record['alignment_supervision'][hidden] = branch_supervision(
        None, record['candidate_environments'][group_key], None)
    return record


def index_row(record, keep):
    labels = record['alignment_supervision'][keep]
    supervised = any(map(any, labels['pair_known_mask'])) or any(labels['allowed_row_mask']) or bool(labels['set_peak_membership'])
    return {'id': record['id'], 'split': record['split'], 'connectivity_key': record['connectivity_key'],
            'mode': keep + '_only', 'has_supervision': supervised,
            'revision': record['provenance'].get('molecule_revision'),
            'selected_spectrum_ids': record['provenance']['selected_spectrum_ids']}


def make_views(source, output, arms):
    manifest = json.loads((source / 'manifest.json').read_text(encoding='utf-8'))
    if manifest.get('schema_version') != SCHEMA_VERSION or manifest.get('state') != 'complete':
        raise ValueError('Expected a completed AtomPeak partial export')
    if len(arms) != len(set(arms)):
        raise ValueError('Duplicate study arms')
    for arm in arms:
        if (output / arm).exists():
            raise ValueError(f'Output already exists: {output / arm}')
    for name, checksum in manifest['vocab_sha256'].items():
        if hashlib.sha256((source / 'dict/bart' / name).read_bytes()).hexdigest() != checksum:
            raise ValueError(f'Source vocabulary checksum differs: {name}')

    manifests, fingerprints = {}, {}
    for arm in arms:
        target = output / arm
        (target / 'lmdb_dataset').mkdir(parents=True, exist_ok=False)
        (target / 'dict/bart').mkdir(parents=True)
        for name in manifest['vocab_sha256']:
            shutil.copyfile(source / 'dict/bart' / name, target / 'dict/bart' / name)
        rich_h = arm != 'hf_shift'
        tokenization = {**manifest['tokenization'], 'use_formula': True, 'use_solvent': False,
                        'use_temperature': False, 'use_field_strength': False, 'use_range_shift': False,
                        'use_h_category': rich_h, 'use_h_nh': rich_h, 'use_h_jvalue': rich_h}
        value = {k: v for k, v in manifest.items()
                 if k not in ('classes', 'warnings', 'skipped_reasons', 'elapsed_seconds')}
        value.update(state='running', tokenization=tokenization,
                     single_nucleus_study={'arm': arm, 'visible_nucleus': NUCLEI[ARMS[arm]],
                                           'source_snapshot_fingerprint': manifest['source_snapshot_fingerprint'],
                                           'cohort': 'same_dual_records_and_original_splits'})
        manifests[arm] = value
        fingerprints[arm] = hashlib.sha256(json.dumps(tokenization, sort_keys=True).encode())
        dump(target / 'manifest.json', value)

    counts = Counter()
    labelled = {arm: Counter() for arm in arms}
    matched_ids = set()
    connectivity_splits = {}
    ordered_ids = hashlib.sha256()
    for split in ('train', 'valid', 'test'):
        input_file = source / 'lmdb_dataset' / f'{split}.lmdb'
        with ExitStack() as stack:
            env = stack.enter_context(lmdb.open(str(input_file), subdir=False, readonly=True, lock=False))
            source_tx = stack.enter_context(env.begin())
            writers, indexes = {}, {}
            for arm in arms:
                writer = stack.enter_context(lmdb.open(str(output / arm / 'lmdb_dataset' / f'{split}.lmdb'),
                    subdir=False, map_size=max(input_file.stat().st_size * 2, 16 * 1024**2)))
                writers[arm] = stack.enter_context(writer.begin(write=True))
                indexes[arm] = stack.enter_context((output / arm / 'records_index.jsonl').open(
                    'x' if split == 'train' else 'a', encoding='utf-8'))
            for _, blob in source_tx.cursor():
                row = pickle.loads(blob)
                if not all(row['branch_present'][b] and row['spectra'][b]['peaks'] for b in BRANCHES):
                    continue
                mid = row['id']
                if row['split'] != split or mid in matched_ids:
                    raise ValueError(f'Split mismatch or duplicate record: {mid}')
                if connectivity_splits.setdefault(row['connectivity_key'], split) != split:
                    raise ValueError(f'Connectivity occurs in multiple splits: {mid}')
                matched_ids.add(mid)
                ordered_ids.update(f'{split}\t{mid}\n'.encode())
                key = f'{counts[split]:012d}'.encode()
                for arm in arms:
                    view = view_record(blob, ARMS[arm])
                    payload = pickle.dumps(view, protocol=4)
                    writers[arm].put(key, payload)
                    entry = index_row(view, ARMS[arm])
                    indexes[arm].write(json.dumps(entry, ensure_ascii=False) + '\n')
                    labelled[arm][split] += entry['has_supervision']
                    fingerprints[arm].update(f'{split}\t{mid}\n'.encode() + payload)
                counts[split] += 1
        if not counts[split]:
            raise ValueError(f'No paired records in {split}')

    for arm in arms:
        value = manifests[arm]
        value.update(state='complete', records=sum(counts.values()), selected_molecules=sum(counts.values()),
                     split_counts=dict(counts), source_snapshot_fingerprint=fingerprints[arm].hexdigest())
        value['single_nucleus_study']['ordered_cohort_sha256'] = ordered_ids.hexdigest()
        value['single_nucleus_study']['supervised_records'] = dict(labelled[arm])
        value['single_nucleus_study']['natural_train_labelled_fraction'] = labelled[arm]['train'] / counts['train']
        dump(output / arm / 'manifest.json', value)
    return {'splits': dict(counts), 'ordered_cohort_sha256': ordered_ids.hexdigest(),
            'arms': {arm: str(output / arm) for arm in arms}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--output-root', required=True, type=Path)
    parser.add_argument('--arms', choices=ARMS, nargs='+', default=['hf', 'cf'])
    args = parser.parse_args()
    print(json.dumps(make_views(args.source, args.output_root, args.arms), indent=2))


if __name__ == '__main__':
    main()
