#!/usr/bin/env python
"""Fit validation-only, per-beam calibration for C/H-UOT candidate reranking."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import pickle
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from rdkit import Chem


def load_calibration_utilities():
    """Load the shared scorer without importing the complete Uni-Core package."""

    utility_path = (
        Path(__file__).resolve().parents[1]
        / "nmrpeak"
        / "utils"
        / "rerank_calibration.py"
    )
    spec = importlib.util.spec_from_file_location("atompeak_rerank_calibration", utility_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load calibration utilities from {utility_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_calibration_utilities = load_calibration_utilities()
DEFAULT_STANDARDIZATION = _calibration_utilities.DEFAULT_STANDARDIZATION
apply_calibration = _calibration_utilities.apply_calibration


Signal = Tuple[float, float, float]
FeatureRow = Tuple[Signal, bool]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Fit generation/C-UOT/H-UOT weights on validation candidate beams. "
            "The test set is never read by this program."
        )
    )
    parser.add_argument("input", help="Raw dual-UOT validation .pkl from infer_atompeak.py")
    parser.add_argument("output", help="Calibration JSON to write")
    parser.add_argument(
        "--step",
        type=float,
        default=0.05,
        help="Simplex grid step for three non-negative weights (default: 0.05).",
    )
    parser.add_argument(
        "--standardization",
        choices=("robust_zscore", "zscore"),
        default="robust_zscore",
        help="Per-beam score standardization method.",
    )
    parser.add_argument(
        "--clip",
        type=float,
        default=3.0,
        help="Absolute cap after score standardization (0 disables clipping).",
    )
    parser.add_argument(
        "--target-match",
        choices=("exact", "connectivity"),
        default="exact",
        help="Use isomeric canonical SMILES or connectivity-only SMILES as the correct label.",
    )
    return parser.parse_args()


def canonical_smiles(smiles: str, target_match: str):
    if not smiles:
        return None
    try:
        molecule = Chem.MolFromSmiles(str(smiles))
    except Exception:
        molecule = None
    if molecule is None:
        return None
    return Chem.MolToSmiles(
        molecule,
        canonical=True,
        isomericSmiles=target_match == "exact",
    )


def row_key(row: Dict, target_match: str):
    canonical = row.get("canonical_smiles")
    if not canonical:
        return None
    return canonical_smiles(canonical, target_match)


def feature_beams(
    samples: Iterable[Dict],
    standardization: Dict,
    target_match: str,
) -> Tuple[List[List[FeatureRow]], Dict[str, int]]:
    """Convert raw rows to normalized features and validation correctness labels."""

    calibration = {
        "weights": {"generation": 1.0, "carbon": 1.0, "proton": 1.0},
        "standardization": standardization,
    }
    beams: List[List[FeatureRow]] = []
    counts = {
        "samples": 0,
        "target_invalid": 0,
        "rankable": 0,
        "oracle_beams": 0,
        "pairwise_beams": 0,
    }
    for sample in samples:
        counts["samples"] += 1
        target = canonical_smiles(sample.get("target_smiles"), target_match)
        if target is None:
            counts["target_invalid"] += 1
            continue

        rows = apply_calibration(sample.get("reranked", []), calibration)
        features: List[FeatureRow] = []
        for row in rows:
            if not row.get("kept", False) or not math.isfinite(
                float(row.get("joint_score", float("-inf")))
            ):
                continue
            candidate = row_key(row, target_match)
            if candidate is None:
                continue
            signal = (
                float(row["generation_normalized"]),
                float(row["carbon_uot_normalized"]),
                float(row["proton_uot_normalized"]),
            )
            features.append((signal, candidate == target))

        if not features:
            continue
        counts["rankable"] += 1
        correct = sum(is_correct for _, is_correct in features)
        if correct:
            counts["oracle_beams"] += 1
        if 0 < correct < len(features):
            counts["pairwise_beams"] += 1
            beams.append(features)
    return beams, counts


def simplex_grid(step: float):
    if not 0 < step <= 1:
        raise ValueError("--step must lie in (0, 1]")
    units = round(1.0 / step)
    if not math.isclose(units * step, 1.0, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError("--step must divide 1 exactly, e.g. 0.1, 0.05 or 0.025")
    for generation in range(units + 1):
        for carbon in range(units - generation + 1):
            proton = units - generation - carbon
            yield (
                generation / units,
                carbon / units,
                proton / units,
            )


def score(signal: Signal, weights: Signal) -> float:
    return sum(value * weight for value, weight in zip(signal, weights))


def evaluate_weights(beams: Sequence[Sequence[FeatureRow]], weights: Signal):
    wins = 0.0
    pair_count = 0
    top1_correct = 0
    for beam in beams:
        values = [(score(signal, weights), correct) for signal, correct in beam]
        correct_scores = [value for value, correct in values if correct]
        incorrect_scores = [value for value, correct in values if not correct]
        best_correct = max(correct_scores)
        for incorrect in incorrect_scores:
            if math.isclose(best_correct, incorrect, abs_tol=1e-12):
                wins += 0.5
            elif best_correct > incorrect:
                wins += 1.0
            pair_count += 1
        best_score = max(value for value, _ in values)
        if any(
            correct and math.isclose(value, best_score, abs_tol=1e-12)
            for value, correct in values
        ):
            top1_correct += 1

    return {
        "pairwise_win_rate": wins / pair_count if pair_count else float("nan"),
        "top1_oracle_rate": top1_correct / len(beams) if beams else float("nan"),
        "pair_count": pair_count,
        "beam_count": len(beams),
    }


def select_weights(beams: Sequence[Sequence[FeatureRow]], step: float):
    best_weights = None
    best_metrics = None
    for weights in simplex_grid(step):
        metrics = evaluate_weights(beams, weights)
        candidate_key = (
            metrics["pairwise_win_rate"],
            metrics["top1_oracle_rate"],
            weights[0],
            weights[1],
        )
        if best_metrics is None or candidate_key > best_metrics[0]:
            best_weights = weights
            best_metrics = (candidate_key, metrics)
    if best_weights is None:
        raise RuntimeError(
            "No validation beam contains both a correct and an incorrect retained candidate."
        )
    return best_weights, best_metrics[1]


def main():
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    with input_path.open("rb") as handle:
        samples = pickle.load(handle)
    if not isinstance(samples, list):
        raise TypeError("The inference pickle must contain a list of per-spectrum records")

    standardization = dict(DEFAULT_STANDARDIZATION)
    standardization.update(
        {
            "method": args.standardization,
            "clip": float(args.clip),
            "scope": "per_beam",
        }
    )
    beams, counts = feature_beams(samples, standardization, args.target_match)
    weights, metrics = select_weights(beams, args.step)
    ablation = {
        "generation_only": evaluate_weights(beams, (1.0, 0.0, 0.0)),
        "carbon_only": evaluate_weights(beams, (0.0, 1.0, 0.0)),
        "proton_only": evaluate_weights(beams, (0.0, 0.0, 1.0)),
        "carbon_proton_equal": evaluate_weights(beams, (0.0, 0.5, 0.5)),
        "equal_three_way": evaluate_weights(beams, (1.0 / 3.0,) * 3),
    }
    result = {
        "version": 1,
        "weights": {
            "generation": weights[0],
            "carbon": weights[1],
            "proton": weights[2],
        },
        "standardization": standardization,
        "selection": {
            "split": "validation",
            "objective": "pairwise_win_rate",
            "target_match": args.target_match,
            "grid_step": args.step,
            **counts,
            **metrics,
        },
        "validation_ablation": ablation,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")

    print("Calibration fitted on validation candidates only")
    print(f"samples={counts['samples']}")
    print(f"target_invalid={counts['target_invalid']}")
    print(
        "valid_target_without_rankable_beam="
        + str(counts['samples'] - counts['target_invalid'] - counts['rankable'])
    )
    print(f"rankable_beams={counts['rankable']}")
    print(f"oracle_beams={counts['oracle_beams']}")
    print(f"pairwise_beams={metrics['beam_count']}")
    print(f"pairwise_win_rate={metrics['pairwise_win_rate']:.6f}")
    print(f"top1_oracle_rate={metrics['top1_oracle_rate']:.6f}")
    print("weights=" + json.dumps(result["weights"], sort_keys=True))
    for name, value in ablation.items():
        print(
            f"{name}: pairwise={value['pairwise_win_rate']:.6f}, "
            f"top1_oracle={value['top1_oracle_rate']:.6f}"
        )
    print(f"output={output_path}")


if __name__ == "__main__":
    main()
