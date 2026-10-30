#!/usr/bin/env python
"""Evaluate fixed-beam Base/UOT inference output without target leakage."""

import argparse
from collections import Counter
import json
import math
import pickle

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors


FINGERPRINT_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048,
)


def molecule_record(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        mol = None
    if mol is None:
        return None
    return {
        "mol": mol,
        "canonical": Chem.MolToSmiles(
            mol, canonical=True, isomericSmiles=True
        ),
        "canonical_no_isomer": Chem.MolToSmiles(
            mol, canonical=True, isomericSmiles=False
        ),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "fingerprint": FINGERPRINT_GENERATOR.GetFingerprint(mol),
    }


def evaluate(samples, top_k):
    maximum_k = max(top_k)
    rejection_counts = Counter()
    raw_total = 0
    raw_valid = 0
    raw_unique = 0
    retained_total = 0
    metrics = {
        k: {
            "recall_hits": 0,
            "recall_no_isomer_hits": 0,
            "valid_slots": 0,
            "unique_candidates": 0,
            "formula_matches": 0,
            "similarity_sum": 0.0,
        }
        for k in top_k
    }

    evaluated_samples = 0
    for sample in samples:
        target = molecule_record(sample.get("target_smiles"))
        if target is None:
            continue
        evaluated_samples += 1

        raw_records = [
            molecule_record(smiles)
            for smiles in sample.get("beam_candidates", [])
        ]
        raw_total += len(raw_records)
        valid_raw = [record for record in raw_records if record is not None]
        raw_valid += len(valid_raw)
        raw_unique += len(
            {record["canonical"] for record in valid_raw}
        )

        rows = sample.get("reranked", [])
        for row in rows:
            if not row.get("kept", row.get("valid", False)):
                rejection_counts[
                    row.get("rejection_reason") or "rejected"
                ] += 1
        kept_rows = [
            row
            for row in rows
            if row.get("kept", row.get("valid", False))
            and math.isfinite(float(row.get("joint_score", float("-inf"))))
        ]
        retained_total += len(kept_rows)

        ranked = []
        for row in kept_rows[:maximum_k]:
            record = molecule_record(
                row.get("canonical_smiles") or row.get("smiles")
            )
            if record is not None:
                ranked.append(record)

        for k in top_k:
            candidates = ranked[:k]
            canonical = [record["canonical"] for record in candidates]
            canonical_no_isomer = [
                record["canonical_no_isomer"] for record in candidates
            ]
            values = metrics[k]
            values["recall_hits"] += target["canonical"] in canonical
            values["recall_no_isomer_hits"] += (
                target["canonical_no_isomer"] in canonical_no_isomer
            )
            values["valid_slots"] += len(candidates)
            values["unique_candidates"] += len(set(canonical))
            values["formula_matches"] += sum(
                record["formula"] == target["formula"]
                for record in candidates
            )
            if candidates:
                values["similarity_sum"] += max(
                    DataStructs.TanimotoSimilarity(
                        target["fingerprint"], record["fingerprint"]
                    )
                    for record in candidates
                )

    sample_count = max(evaluated_samples, 1)
    result = {
        "samples": evaluated_samples,
        "raw_beam": {
            "validity": raw_valid / max(raw_total, 1),
            "uniqueness": raw_unique / max(raw_valid, 1),
            "average_size": raw_total / sample_count,
        },
        "preprocessing": {
            "average_retained": retained_total / sample_count,
            "rejection_counts": dict(sorted(rejection_counts.items())),
        },
        "top_k": {},
    }
    for k in top_k:
        values = metrics[k]
        valid_slots = values["valid_slots"]
        result["top_k"][str(k)] = {
            "recall": values["recall_hits"] / sample_count,
            "recall_no_isomer": (
                values["recall_no_isomer_hits"] / sample_count
            ),
            "validity": valid_slots / (sample_count * k),
            "uniqueness": (
                values["unique_candidates"] / max(valid_slots, 1)
            ),
            "formula_match": (
                values["formula_matches"] / max(valid_slots, 1)
            ),
            "max_similarity": values["similarity_sum"] / sample_count,
        }
    return result


def print_result(label, result):
    print(f"\n{label}")
    print("=" * 92)
    print(f"samples={result['samples']}")
    raw = result["raw_beam"]
    print(
        "raw_beam: "
        f"size={raw['average_size']:.2f}, "
        f"validity={raw['validity']:.4f}, "
        f"uniqueness={raw['uniqueness']:.4f}"
    )
    preprocessing = result["preprocessing"]
    print(
        "preprocessing: "
        f"average_retained={preprocessing['average_retained']:.2f}, "
        f"rejections={preprocessing['rejection_counts']}"
    )
    print(
        f"{'Top-k':<8}{'Recall':>11}{'No isomer':>13}"
        f"{'Validity':>12}{'Unique':>11}{'Formula':>11}"
        f"{'Max sim.':>12}"
    )
    for k, values in result["top_k"].items():
        print(
            f"{k:<8}{values['recall']:>11.4f}"
            f"{values['recall_no_isomer']:>13.4f}"
            f"{values['validity']:>12.4f}"
            f"{values['uniqueness']:>11.4f}"
            f"{values['formula_match']:>11.4f}"
            f"{values['max_similarity']:>12.4f}"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--label", default="AtomPeak generation evaluation")
    parser.add_argument("--top-k", type=int, nargs="+", default=[1, 3, 5, 10])
    parser.add_argument("--output-json")
    args = parser.parse_args()

    with open(args.input, "rb") as handle:
        samples = pickle.load(handle)
    result = evaluate(samples, sorted(set(args.top_k)))
    print_result(args.label, result)

    output_json = args.output_json or f"{args.input}.metrics.json"
    with open(output_json, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False)
    print(f"\nmetrics_json={output_json}")


if __name__ == "__main__":
    main()
