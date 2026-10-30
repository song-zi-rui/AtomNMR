#!/usr/bin/env python
import argparse
import json
import pickle
from pathlib import Path

from nmrpeak.utils.smiles_hit_metrics import (
    evaluate_smiles_hit,
    print_evaluation_results,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate the pickle produced by nmrpeak.infer."
    )
    parser.add_argument("input", help="Inference .out.pkl file")
    parser.add_argument("--beam-size", type=int, required=True)
    parser.add_argument("--k-list", type=int, nargs="+", default=[1, 5, 10])
    parser.add_argument("--output", help="Output JSON path")
    return parser.parse_args()


def join_tokens(tokens):
    if isinstance(tokens, str):
        return tokens
    return "".join(str(token) for token in tokens)


def flatten_inference(raw_batches, beam_size):
    targets = []
    predictions = []
    for batch in raw_batches:
        if not isinstance(batch, (tuple, list)) or len(batch) < 3:
            raise ValueError("Unexpected inference batch format")
        batch_predictions, _, batch_targets = batch
        if len(batch_predictions) != len(batch_targets):
            raise ValueError("Prediction/target batch sizes do not match")
        for beams, target in zip(batch_predictions, batch_targets):
            target_smiles = join_tokens(target)
            beam_smiles = [join_tokens(tokens) for tokens in beams[:beam_size]]
            if len(beam_smiles) < beam_size:
                beam_smiles.extend([""] * (beam_size - len(beam_smiles)))
            targets.append(target_smiles)
            predictions.extend(beam_smiles)
    return targets, predictions


def main():
    args = parse_args()
    k_list = sorted(set(k for k in args.k_list if 0 < k <= args.beam_size))
    if not k_list:
        raise ValueError("--k-list must contain a value between 1 and beam size")

    input_path = Path(args.input)
    with input_path.open("rb") as handle:
        raw_batches = pickle.load(handle)
    targets, predictions = flatten_inference(raw_batches, args.beam_size)
    metrics = evaluate_smiles_hit(
        target_smiles=targets,
        predict_smiles=predictions,
        k_list=k_list,
        beam_size=args.beam_size,
    )
    print(f"samples={len(targets)}")
    print_evaluation_results(metrics, k_list)

    output_path = (
        Path(args.output)
        if args.output
        else input_path.with_suffix(".metrics.json")
    )
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
    print(f"metrics_saved={output_path}")


if __name__ == "__main__":
    main()
