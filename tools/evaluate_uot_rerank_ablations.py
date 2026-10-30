#!/usr/bin/env python
"""Re-rank a fixed Dual-UOT test beam with pre-specified score ablations.

This script never generates new candidates and never fits weights. It reads the
raw generation/C-UOT/H-UOT scores saved in an inference pickle, applies the
validation-fitted per-beam standardization, and writes one result pickle per
ablation. Use it before evaluating the pickles with
tools/evaluate_atompeak_generation.py.
"""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import pickle
from pathlib import Path
from typing import Dict, Iterable, Mapping


MODE_WEIGHTS = {
    "generation": {"generation": 1.0, "carbon": 0.0, "proton": 0.0},
    "carbon": {"generation": 0.0, "carbon": 1.0, "proton": 0.0},
    "proton": {"generation": 0.0, "carbon": 0.0, "proton": 1.0},
    "carbon_proton": {"generation": 0.0, "carbon": 0.5, "proton": 0.5},
}


def load_calibration_utilities():
    utility_path = (
        Path(__file__).resolve().parents[1]
        / "nmrpeak"
        / "utils"
        / "rerank_calibration.py"
    )
    spec = importlib.util.spec_from_file_location(
        "atompeak_rerank_calibration", utility_path
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load calibration utilities from {utility_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CALIBRATION_UTILITIES = load_calibration_utilities()
apply_calibration = CALIBRATION_UTILITIES.apply_calibration
load_calibration = CALIBRATION_UTILITIES.load_calibration


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Apply fixed generation/C-UOT/H-UOT reranking ablations to an "
            "already generated test beam."
        )
    )
    parser.add_argument(
        "input",
        help="Dual-UOT test inference pickle containing raw score fields.",
    )
    parser.add_argument(
        "calibration",
        help="Validation-fitted calibration.json; its standardization is reused.",
    )
    parser.add_argument(
        "output_dir",
        help="Directory for uot_<mode>_test.pkl files.",
    )
    parser.add_argument(
        "--modes",
        nargs="+",
        choices=(*MODE_WEIGHTS, "joint"),
        default=(*MODE_WEIGHTS, "joint"),
        help="Pre-specified ablations to write.",
    )
    return parser.parse_args()


def fixed_calibration(
    calibration: Mapping[str, object], mode: str
) -> Dict[str, object]:
    config = copy.deepcopy(dict(calibration))
    if mode != "joint":
        config["weights"] = dict(MODE_WEIGHTS[mode])
    return config


def rerank_samples(
    samples: Iterable[Mapping[str, object]],
    calibration: Mapping[str, object],
    mode: str,
):
    result = []
    for sample in samples:
        row = copy.deepcopy(dict(sample))
        row["reranked"] = apply_calibration(
            row.get("reranked", []), calibration
        )
        row["rerank_mode"] = f"uot_ablation_{mode}"
        row["rerank_calibration"] = calibration
        result.append(row)
    return result


def main():
    args = parse_args()
    with Path(args.input).open("rb") as handle:
        samples = pickle.load(handle)
    if not isinstance(samples, list):
        raise TypeError("Input pickle must contain a list of per-spectrum records")

    calibration = load_calibration(args.calibration)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "input": str(Path(args.input).resolve()),
        "calibration": str(Path(args.calibration).resolve()),
        "modes": {},
    }

    for mode in args.modes:
        mode_calibration = fixed_calibration(calibration, mode)
        output_path = output_dir / f"uot_{mode}_test.pkl"
        with output_path.open("wb") as handle:
            pickle.dump(
                rerank_samples(samples, mode_calibration, mode), handle
            )
        manifest["modes"][mode] = {
            "path": str(output_path.resolve()),
            "weights": mode_calibration["weights"],
            "standardization": mode_calibration["standardization"],
        }
        print(
            f"{mode}: weights={json.dumps(mode_calibration['weights'], sort_keys=True)}"
        )
        print(f"{mode}: output={output_path}")

    (output_dir / "ablation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()

