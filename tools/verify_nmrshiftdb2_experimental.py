#!/usr/bin/env python
"""Read the final C/H exports through RDKit and independently verify their counts."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys

from rdkit import Chem, RDLogger, rdBase

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from tools.filter_nmrshiftdb2_experimental import SPECTRUM, classify, sha256


def verify(output):
    output = Path(output).resolve()
    target = output / "validation_report.json"
    if target.exists():
        raise FileExistsError(target)
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    report = json.loads((output / "report.json").read_text(encoding="utf-8"))
    if sha256(manifest["source"]["path"]) != manifest["source"]["sha256"]:
        raise ValueError("Original SDF hash mismatch")
    for name, expected in manifest["artifacts"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"Artifact hash mismatch: {name}")
    result = {"passed": True, "rdkit_version": rdBase.rdkitVersion,
              "source_unchanged": True, "manifest_artifact_hashes_valid": True,
              "canonical_structure_definition": "RDKit sanitized canonical isomeric SMILES; no salt stripping or tautomer merging",
              "pools": {}}
    RDLogger.DisableLog("rdApp.*")
    for pool, filename in (("experimental_ch_usable", "experimental_ch.sd"),
                           ("experimental_ch_both_usable", "experimental_ch_both.sd")):
        counts, spectra = Counter(), Counter()
        structures = set()
        mols = Chem.SDMolSupplier(str(output / filename), sanitize=True, removeHs=False)
        for index, mol in enumerate(mols, 1):
            if mol is None:
                raise ValueError(f"RDKit cannot parse exported molecule: {filename}:{index}")
            counts["molecule_records"] += 1
            props = {name: mol.GetProp(name) for name in mol.GetPropNames()}
            local_types = set()
            for name in props:
                match = SPECTRUM.fullmatch(name)
                if not match:
                    continue
                nucleus, sid = match.groups()
                if nucleus not in {"13C", "1H"} or classify(props, sid, nucleus)[0] != "experimental":
                    raise ValueError(f"Non-experimental C/H field: {filename}:{index}:{name}")
                if not props[name].strip():
                    raise ValueError(f"Empty training field: {filename}:{index}:{name}")
                spectra[nucleus] += 1
                local_types.add(nucleus)
            if not local_types:
                raise ValueError(f"RDKit lost spectrum fields: {filename}:{index}")
            if local_types == {"13C", "1H"}:
                counts["molecules_with_both_ch"] += 1
            elif pool == "experimental_ch_both_usable":
                raise ValueError(f"Missing one nucleus in both-nuclei pool: {index}")
            structures.add(Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True))
        expected = report["pools"][pool]
        if dict(spectra) != expected["by_spectrum_type"]:
            raise ValueError(f"RDKit spectrum counts disagree with report: {pool}: {spectra}")
        for key in ("molecule_records", "molecules_with_both_ch"):
            if counts[key] != expected[key]:
                raise ValueError(f"RDKit molecule counts disagree with report: {pool}:{key}")
        result["pools"][pool] = {**counts, "by_spectrum_type": dict(spectra),
                                 "unique_canonical_isomeric_smiles": len(structures)}
        print(f"[verify] {filename}: counts and experimental sources confirmed", flush=True)
    temporary = target.with_suffix(".json.partial")
    temporary.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.rename(target)
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/nmrshiftdb2_experimental_v1")
    verify(parser.parse_args().output_dir)
