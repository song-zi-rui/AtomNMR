#!/usr/bin/env python
import argparse

import numpy as np

from nmrpeak.data.atompeak_gnn_dataset_v2 import AtomPeakGnnDatasetV2
from nmrpeak.data.composite_lmdb_dataset import load_composite_lmdb_dataset
from nmrpeak.data.spectrum_filter_dataset import SpectrumFilterDataset


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate the combined LMDB contract for dual UOT training."
    )
    parser.add_argument("data")
    parser.add_argument("--split", default="train")
    parser.add_argument("--data-sources", default="nmrexp,nmrshiftdb2")
    parser.add_argument("--samples", type=int, default=500)
    parser.add_argument("--allow-single-nucleus", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
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
    features = AtomPeakGnnDatasetV2(dataset)
    checked = min(args.samples, len(features))
    supervised = {"carbon": 0, "proton": 0}

    for index in range(checked):
        item = features[index]
        atom_count = item["graph"]["node_features"].shape[0]
        for branch in ("carbon", "proton"):
            data = item[branch]
            peak_count = data["peak_shifts"].shape[0]
            group_count = data["group_membership"].shape[0]
            assert data["group_membership"].shape == (group_count, atom_count)
            assert data["target"].shape == (peak_count, group_count)
            assert data["target_fraction"].shape == (peak_count, group_count)
            assert data["hard_mask"].shape == (peak_count, group_count)
            assert np.all(data["peak_shifts"][:-1] >= data["peak_shifts"][1:])
            assert np.all(data["target"] <= data["hard_mask"])
            supervised[branch] += int(data["assignment_valid"] > 0)

    print(f"usable_records={len(dataset)}")
    print(f"checked_records={checked}")
    print(f"carbon_supervised={supervised['carbon']}")
    print(f"proton_supervised={supervised['proton']}")
    print("dataset_contract=OK")


if __name__ == "__main__":
    main()
