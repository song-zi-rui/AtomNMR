#!/usr/bin/env python
"""Experimental-only NMRShiftDB2 guard for the existing AtomPeak dataset pipeline.

All existing stages/arguments remain available. NMRexp is not reclassified by
this wrapper. Use a NEW output directory; old mixed-source datasets are untouched.
"""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.filter_nmrshiftdb2_experimental import POLICY, assert_experimental_ch, sha256
from tools.shiftdb_atom_identity import POLICY as IDENTITY_POLICY


def main(argv=None):
    from tools import prepare_atompeak_dataset as pipeline

    original_validate = pipeline.validate_args
    original_record = pipeline.clean_converter.molecule_record
    original_schema = pipeline.SCHEMA_VERSION
    source = {}

    def validate(args):
        original_validate(args)
        counts = assert_experimental_ch(args.sdf_path)
        signature = {"policy": POLICY, "assignment_identity_policy": IDENTITY_POLICY, "sdf_sha256": sha256(args.sdf_path), "spectra": counts}
        output = Path(args.output_dir).resolve()
        marker = output / "experimental_input.json"
        if output.exists() and any(output.iterdir()):
            if not marker.exists() or json.loads(marker.read_text(encoding="utf-8")) != signature:
                raise ValueError("Use a new output directory: this directory is not the same experimental-only run.")
        output.mkdir(parents=True, exist_ok=True)
        marker.write_text(json.dumps(signature, indent=2) + "\n", encoding="utf-8")
        source.update(signature, sdf_path=str(Path(args.sdf_path).resolve()))
        pipeline.log(f"Experimental-only NMRShiftDB2 input validated: {counts}")

    def molecule_record(mol, args):
        record, reason = original_record(mol, args)
        if record is not None:
            record["provenance"].update(source)
            record["provenance"]["source"] = source["sdf_path"]
            if mol.HasProp("Original SDF record index"):
                record["provenance"]["original_sdf_record_index"] = int(mol.GetProp("Original SDF record index"))
            if mol.HasProp("Original SDF SHA256"):
                record["provenance"]["original_sdf_sha256"] = mol.GetProp("Original SDF SHA256")
        return record, reason

    pipeline.validate_args = validate
    pipeline.clean_converter.molecule_record = molecule_record
    pipeline.SCHEMA_VERSION = "atompeak_pipeline_v2_experimental_atom_identity_v1"
    try:
        return pipeline.main(argv)
    finally:
        pipeline.validate_args = original_validate
        pipeline.clean_converter.molecule_record = original_record
        pipeline.SCHEMA_VERSION = original_schema


if __name__ == "__main__":
    raise SystemExit(main())
