#!/usr/bin/env python
import argparse
import collections
import json
import pickle

import lmdb


ASSIGNMENT_KEYS = {
    "assignment_matrix",
    "atom_peak_assignment",
    "atom_peak_assignments",
    "peak_atom_assignment",
    "peak_atom_assignments",
    "assignments",
}

PEAK_KEYS = {
    "c_nmr_peaks",
    "cnmr_peaks",
    "carbon_peaks",
    "h_nmr_peaks",
    "hnmr_peaks",
    "proton_peaks",
    "peaks",
    "spectrum_peaks",
}

ATOM_KEYS = {
    "atom_groups",
    "assignment_atom_groups",
    "assigned_atoms",
    "atoms",
    "atom_symbols",
    "smiles",
    "smiles_canonical",
    "canonical_smiles",
}


def load_lmdb_items(path, limit):
    env = lmdb.open(path, subdir=False, readonly=True, lock=False, readahead=False, meminit=False, max_readers=16)
    with env.begin() as txn:
        for idx, (_, value) in enumerate(txn.cursor()):
            if idx >= limit:
                break
            yield pickle.loads(value)


def summarise_value(value):
    if isinstance(value, list):
        return {"type": "list", "len": len(value)}
    if isinstance(value, tuple):
        return {"type": "tuple", "len": len(value)}
    if isinstance(value, dict):
        return {"type": "dict", "keys": sorted(list(value.keys()))[:20]}
    shape = getattr(value, "shape", None)
    if shape is not None:
        return {"type": type(value).__name__, "shape": list(shape)}
    return {"type": type(value).__name__}


def main():
    parser = argparse.ArgumentParser(description="Audit atom-peak assignment fields in an NMRPeak-style LMDB file.")
    parser.add_argument("lmdb_path")
    parser.add_argument("--limit", type=int, default=100)
    args = parser.parse_args()

    key_counts = collections.Counter()
    assignment_hits = collections.Counter()
    peak_hits = collections.Counter()
    atom_hits = collections.Counter()
    examples = {}
    total = 0

    for item in load_lmdb_items(args.lmdb_path, args.limit):
        total += 1
        if not isinstance(item, dict):
            continue
        for key in item:
            key_counts[key] += 1
            examples.setdefault(key, summarise_value(item[key]))
        for key in ASSIGNMENT_KEYS & set(item):
            assignment_hits[key] += 1
        for key in PEAK_KEYS & set(item):
            peak_hits[key] += 1
        for key in ATOM_KEYS & set(item):
            atom_hits[key] += 1

    report = {
        "sampled_records": total,
        "assignment_fields": dict(assignment_hits),
        "peak_fields": dict(peak_hits),
        "atom_or_smiles_fields": dict(atom_hits),
        "all_fields": dict(key_counts.most_common()),
        "field_examples": examples,
        "ready_for_assignment_training": bool(assignment_hits and peak_hits and atom_hits),
    }
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
