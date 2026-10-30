#!/usr/bin/env python
import argparse
import os

from unicore.data import Dictionary, TokenizeDataset

from nmrpeak.data import (
    KeyDataset,
    RandomSmilesDataset,
    SmilesTokenDataset,
    TokenizedNmrDataset,
)
from nmrpeak.data.composite_lmdb_dataset import load_composite_lmdb_dataset
from nmrpeak.data.spectrum_filter_dataset import SpectrumFilterDataset


def parse_args():
    parser = argparse.ArgumentParser(
        description="Scan sequence lengths and token IDs before GPU training."
    )
    parser.add_argument("data")
    parser.add_argument("--dict-path", required=True)
    parser.add_argument("--split", default="train")
    parser.add_argument("--data-sources", default="nmrexp,nmrshiftdb2")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--max-seq-len", type=int, default=512)
    parser.add_argument("--model-max-positions", type=int, default=1024)
    parser.add_argument("--shift-length", type=int, default=2795)
    parser.add_argument("--allow-single-nucleus", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    dictionary = Dictionary.load(
        os.path.join(args.dict_path, "bart/total_dict.txt")
    )
    raw = load_composite_lmdb_dataset(
        args.data,
        args.split,
        data_sources=args.data_sources,
    )
    dataset = SpectrumFilterDataset(
        raw,
        use_cnmr=True,
        use_hnmr=True,
        allow_single_nucleus=args.allow_single_nucleus,
    )

    spectrum_text = TokenizedNmrDataset(
        dataset,
        use_cnmr=True,
        use_hnmr=True,
        use_h_nh=False,
        use_h_jvalue=False,
        use_h_category=False,
        use_range_shift=False,
        use_formula=True,
        use_spec_end_token=True,
    )
    spectrum_tokens = TokenizeDataset(
        spectrum_text,
        dictionary,
        max_seq_len=args.max_seq_len,
    )
    smiles = KeyDataset(dataset, "smiles_canonical")
    smiles = RandomSmilesDataset(smiles, seed=1, prob=0)
    smiles_tokens = SmilesTokenDataset(
        smiles, shift_length=args.shift_length
    )

    total = len(dataset) if args.limit <= 0 else min(args.limit, len(dataset))
    maximum = {
        "spec_length": (0, -1),
        "smiles_length": (0, -1),
        "spec_token": (-1, -1),
        "smiles_token": (-1, -1),
    }
    failures = []
    for index in range(total):
        spec = spectrum_tokens[index]
        mol = smiles_tokens[index]
        spec_length = len(spec) + 2
        smiles_length = len(mol) + 2
        spec_max = int(spec.max()) if len(spec) else -1
        smiles_max = int(mol.max()) if len(mol) else -1

        for name, value in (
            ("spec_length", spec_length),
            ("smiles_length", smiles_length),
            ("spec_token", spec_max),
            ("smiles_token", smiles_max),
        ):
            if value > maximum[name][0]:
                maximum[name] = (value, index)

        reasons = []
        if spec_length > args.model_max_positions:
            reasons.append(f"spec_length={spec_length}")
        if smiles_length > args.model_max_positions:
            reasons.append(f"smiles_length={smiles_length}")
        if spec_max >= len(dictionary):
            reasons.append(f"spec_token={spec_max}")
        if smiles_max >= len(dictionary):
            reasons.append(f"smiles_token={smiles_max}")
        if reasons and len(failures) < 20:
            item = dataset[index]
            failures.append(
                {
                    "index": index,
                    "smiles": item.get("smiles_canonical"),
                    "reasons": reasons,
                }
            )
        if (index + 1) % 10000 == 0:
            print(f"scanned={index + 1}/{total}", flush=True)

    print(f"dictionary_size={len(dictionary)}")
    print(f"records_scanned={total}")
    for name, (value, index) in maximum.items():
        print(f"max_{name}={value} at_index={index}")
    if failures:
        print("status=FAILED")
        for failure in failures:
            print(failure)
        raise SystemExit(1)
    print("status=OK")


if __name__ == "__main__":
    main()
