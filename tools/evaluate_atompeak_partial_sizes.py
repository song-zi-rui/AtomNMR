"""Add size/mode statistics to saved partial-v1 JSONL without model inference."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pickle
import sys
from pathlib import Path

os.environ['NMRPEAK_LIGHT_IMPORT'] = '1'
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rdkit import RDLogger
from nmrpeak.atompeak.generation_metrics_partial import GenerationMetrics, parse_structure


def evaluate_saved(predictions, data, split, output):
    import lmdb
    if output.exists():
        raise ValueError('Choose a new metrics file; existing results are never overwritten')
    metrics = GenerationMetrics()
    skipped_missing_input = 0
    with lmdb.open(str(data/'lmdb_dataset'/f'{split}.lmdb'), subdir=False, readonly=True, lock=False) as env:
        with env.begin() as tx:
            records = {r['id']: r for _, b in tx.cursor() if (r := pickle.loads(b))}
    expected = set()
    modes = set()
    with predictions.open(encoding='utf-8') as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            mid = row['id']
            if mid not in records:
                raise ValueError('Prediction ID not present in selected split: ' + mid)
            record = records[mid]
            mode = row.get('input_mode', 'as_stored')
            if mode not in ('as_stored', 'carbon_only', 'proton_only'):
                raise ValueError('Unsupported input mode: ' + mode)
            modes.add(mode)
            if mode != 'as_stored':
                keep = 'carbon' if mode == 'carbon_only' else 'proton'
                if not record['branch_present'][keep]:
                    raise ValueError('Prediction uses a missing input branch: ' + mid)
                record = {**record, 'branch_present': {b: b == keep for b in ('carbon', 'proton')}}
            target = parse_structure(row.get('target_for_evaluation_only'))
            truth = parse_structure(record['smiles_canonical'])
            if target is None or target['canonical'] != truth['canonical']:
                raise ValueError('Stored target differs from frozen split: ' + mid)
            metrics.update(record, row['beams'])
    if len(modes) > 1:
        raise ValueError('Use one input-mode experiment per file')
    mode = next(iter(modes), 'as_stored')
    for mid, record in records.items():
        keep = 'carbon' if mode == 'carbon_only' else 'proton'
        if mode == 'as_stored' or record['branch_present'][keep]:
            expected.add(mid)
        else:
            skipped_missing_input += 1
    result = {**metrics.result(), 'prediction_file': str(predictions.resolve()),
              'predictions_sha256': hashlib.sha256(predictions.read_bytes()).hexdigest(),
              'data': str(data.resolve()), 'split': split, 'input_mode': mode,
              'eligible_split_records': len(expected), 'unprocessed_eligible_records': len(expected-metrics.seen),
              'skipped_missing_input': skipped_missing_input}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--predictions', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--split', choices=('train', 'valid', 'test'), default='test')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.error')
    result = evaluate_saved(args.predictions, args.data, args.split, args.output)
    print(json.dumps({'processed': result['processed'], 'by_size': result['by_size'],
                      'output': str(args.output.resolve())}, ensure_ascii=False))


if __name__ == '__main__':
    main()
