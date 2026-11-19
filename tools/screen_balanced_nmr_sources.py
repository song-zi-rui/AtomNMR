"""Screen downloaded computational shifts and an independent experimental carbon panel."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sqlite3
import zipfile

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/nmr_atom_assignment_sources"
OUT = ROOT / "outputs/nmr_balance_selection_20261119"
DB = ROOT / "data/assignment_database_v1/databases/assignments.sqlite"
SPLIT = ROOT / "data/assignment_database_v1/exports/atompeak_complete_ch_10003_20261104_training/records_index.jsonl"
GYM = ROOT / "outputs/experimental_2d_nmr_overlap_20261112/2dnmrgym_overlap.csv"
CASPRE = DATA / "caspre_jeol/caspre_train_CDCl3_JEOL.zip"
DFT = DATA / "cascade_dft8k"
SPEC = {
    "sulfonamide": "[SX4](=[OX1])(=[OX1])[NX3]",
    "aromatic_sulfur": "[s]",
    "thioether": "[#6]-[SX2]-[#6]",
    "amide": "[CX3](=[OX1])[NX3]",
}
PATTERNS = {name: Chem.MolFromSmarts(smarts) for name, smarts in SPEC.items()}
DFT_ELEMENTS = {1, 6, 7, 8, 9, 16, 17, 35}


def csv_rows(data):
    return csv.DictReader(io.StringIO(data.decode("utf-8-sig"), newline=""))


def size_band(n):
    if n < 16:
        return "<=15"
    if n < 26:
        return "16-25"
    if n < 36:
        return "26-35"
    return "36-45" if n < 46 else "46+"


def feature_record(mol):
    heavy = Chem.RemoveHs(mol)
    key = Chem.MolToInchiKey(heavy)
    if len(key) != 27:
        raise ValueError("InChIKey could not be computed")
    # Strip stereo from the scaffold input only. Some scaffold bond stereo
    # references the side chains that were removed.
    scaffold_input = Chem.Mol(heavy)
    Chem.RemoveStereochemistry(scaffold_input)
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=scaffold_input)
    elements = {a.GetAtomicNum() for a in mol.GetAtoms()}
    rings = heavy.GetRingInfo().AtomRings()
    ar = [r for r in rings if all(heavy.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    het_ar = sum(any(heavy.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in r) for r in ar)
    c0 = sum(a.GetAtomicNum() == 6 and a.GetTotalNumHs() == 0 for a in heavy.GetAtoms())
    flags = {n for n, q in PATTERNS.items() if heavy.HasSubstructMatch(q)}
    if ar and "amide" in flags:
        flags.add("aryl_amide")
    if 16 in elements:
        flags.add("sulfur")
    if 9 in elements:
        flags.add("fluorinated")
    if 17 in elements:
        flags.add("chlorinated")
    if het_ar:
        flags.add("heteroaromatic")
    if het_ar and elements & {9, 17, 35}:
        flags.add("halogen_heteroaromatic")
    if len(rings) >= 3:
        flags.add("polycyclic")
    if c0 >= 12:
        flags.add("nonprotonated_c")
    stereo_known = all(s.specified == Chem.StereoSpecified.Specified for s in Chem.FindPotentialStereo(heavy))
    organic = (
        elements <= DFT_ELEMENTS
        and len(Chem.GetMolFrags(heavy)) == 1
        and all(not a.GetFormalCharge() and not a.GetNumRadicalElectrons() for a in heavy.GetAtoms())
    )
    n = heavy.GetNumHeavyAtoms()
    return {
        "canonical_smiles": Chem.MolToSmiles(heavy),
        "inchikey": key,
        "connectivity": key[:14],
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "heavy_atoms": n,
        "band": size_band(n),
        "scaffold": scaffold or Chem.MolToSmiles(scaffold_input),
        "features": sorted(flags),
        "nonprotonated_carbon_count": c0,
        "heteroaromatic_rings": het_ar,
        "structure_has_defined_stereo": stereo_known,
        "neutral_single_component_supported_elements": organic,
    }


def shifts_by_index(shift_rows, field):
    result = {}
    for row in shift_rows:
        idx, value = int(row["atom_index"]), float(row[field])
        if idx in result or not math.isfinite(value):
            raise ValueError("duplicate atom index or nonfinite shift")
        result[idx] = value
    return result


def write_jsonl(path, rows):
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        for r in rows:
            stream.write(json.dumps(r, ensure_ascii=False, allow_nan=False) + "\n")


def annotate_overlap(rows, connection):
    by_layer = defaultdict(list)
    for row in rows:
        by_layer[row["connectivity"]].append(row)
    for layer, members in by_layer.items():
        # Complete keys start with the connectivity block followed by a hyphen.
        hits = list(connection.execute(
            "SELECT id,inchikey,queue FROM molecules WHERE inchikey>=? AND inchikey<?",
            (layer + "-", layer + "."),
        ))
        for row in members:
            exact = [{"molecule_id": mid, "queue": q} for mid, key, q in hits if key == row["inchikey"]]
            row["main_database_match"] = "exact" if exact else "connectivity_only" if hits else "new"
            row["main_database_exact_hits"] = exact


def select(rows, quotas, preferred):
    picked, used = [], set()
    scaffolds = Counter()
    for band, count in quotas.items():
        pool = [r for r in rows if r["band"] == band]
        order = preferred[band]
        feature_counts = Counter()
        for _ in range(count):
            best = None
            for r in pool:
                if r["connectivity"] in used:
                    continue
                relevant = set(r["features"]) & set(order)
                if not relevant:
                    continue
                feature = min(relevant, key=lambda key: (feature_counts[key], order.index(key)))
                if scaffolds[r["scaffold"]] >= 2:
                    continue
                score = (
                    feature_counts[feature],
                    scaffolds[r["scaffold"]],
                    -len(relevant),
                    hashlib.sha256(r["connectivity"].encode()).hexdigest(),
                )
                if best is None or score < best[0]:
                    best = score, r, feature
            if best is None:
                break
            _, row, feature = best
            picked.append({**row, "selection_group": feature})
            used.add(row["connectivity"])
            scaffolds[row["scaffold"]] += 1
            feature_counts[feature] += 1
    return picked


def read_caspre(*, include_source_structure=False):
    rejected = Counter()
    rows = []
    with zipfile.ZipFile(CASPRE) as archive:
        stem = "caspre_train_CDCl3_JEOL"
        metadata = {r["record"]: r for r in csv_rows(archive.read(stem + "_molecules.csv"))}
        assignments = defaultdict(list)
        for r in csv_rows(archive.read(stem + "_shifts.csv")):
            assignments[r["record"]].append(r)
        source = gzip.decompress(archive.read(stem + ".sdf.gz"))
    if len(metadata) != len(source.split(b"$$$$")) - 1:
        raise ValueError("CASPRE structure and index lengths differ")
    for record, mol in enumerate(Chem.ForwardSDMolSupplier(io.BytesIO(source), removeHs=False), 1):
        base = metadata[str(record)]
        if mol is None or mol.GetProp("_Name").strip() != base["molecule_id"]:
            rejected["structure_record_mismatch"] += 1
            continue
        try:
            r = feature_record(mol)
            raw = assignments[str(record)]
            shifts = shifts_by_index(raw, "shift_ppm")
        except (ValueError, RuntimeError):
            rejected["invalid_atom_or_structure_data"] += 1
            continue
        if r["inchikey"] != base["inchikey"]:
            rejected["structure_key_mismatch"] += 1
            continue
        if base["solvent"] != "CDCl3" or any(x["solvent"] != base["solvent"] or x["molecule_id"] != base["molecule_id"] for x in raw):
            rejected["source_metadata_mismatch"] += 1
            continue
        carbons = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6}
        if not carbons or set(shifts) != carbons or len(shifts) != int(base["n_shifts"]):
            rejected["missing_carbon_or_inconsistent_atom_index"] += 1
            continue
        rows.append({
            **r,
            "source": "CASPRE-JEOL",
            "source_id": base["molecule_id"],
            "source_sdf_record_1based": record,
            "compound_name": base["compound_name"],
            "solvent": base["solvent"],
            "label_origin": "experimental_source_assignment",
            "atom_index_basis": "source_sdf_explicit_h_zero_based",
            "c_shifts": [{"atom_index": i, "ppm": value} for i, value in sorted(shifts.items())],
            **({"source_mol": mol} if include_source_structure else {}),
        })
    return rows, {"source_records": len(metadata), "invalid": dict(rejected)}


def read_dft(*, include_source_structure=False):
    assignments = defaultdict(list)
    with gzip.open(DFT / "DFT8K.csv.gz", "rb") as stream:
        for r in csv_rows(stream.read()):
            assignments[r["mol_id"]].append(r)
    rejected, rows, seen = Counter(), [], set()
    with gzip.open(DFT / "DFT.sdf.gz", "rb") as stream:
        for record, mol in enumerate(Chem.ForwardSDMolSupplier(stream, removeHs=False), 1):
            if mol is None:
                rejected["invalid_structure"] += 1
                continue
            mid = mol.GetProp("_Name").strip()
            if mid not in assignments:
                continue
            seen.add(mid)
            try:
                r = feature_record(mol)
                raw = assignments[mid]
                shifts = shifts_by_index(raw, "Shift")
                actual = {a.GetIdx(): a.GetAtomicNum() for a in mol.GetAtoms() if a.GetAtomicNum() in {1, 6}}
                if actual != {int(x["atom_index"]): int(x["atom_type"]) for x in raw}:
                    raise ValueError("incorrect element or missing H/C label")
            except (ValueError, RuntimeError):
                rejected["invalid_or_incomplete_atom_mapping"] += 1
                continue
            if set(actual.values()) != {1, 6}:
                rejected["missing_nucleus"] += 1
                continue
            h, c = [], []
            for idx, value in sorted(shifts.items()):
                rec = {"atom_index": idx, "ppm": value}
                if actual[idx] == 6:
                    c.append(rec)
                else:
                    neighbors = mol.GetAtomWithIdx(idx).GetNeighbors()
                    rec["attached_atom_index"] = neighbors[0].GetIdx() if len(neighbors) == 1 else None
                    h.append(rec)
            rows.append({
                **r,
                "source": "CASCADE-DFT8K",
                "source_id": mid,
                "source_sdf_record_1based": record,
                "label_origin": "dft_calculated",
                "atom_index_basis": "source_sdf_explicit_h_zero_based",
                "h_shifts": h,
                "c_shifts": c,
                **({"source_mol": mol} if include_source_structure else {}),
            })
    rejected["label_id_not_found_in_sdf"] = len(set(assignments) - seen)
    return rows, {"source_records": len(assignments), "invalid": dict(rejected)}


def group_counts(rows):
    return {
        "records": len(rows),
        "distinct_connectivity": len({r["connectivity"] for r in rows}),
        "by_band": dict(Counter(r["band"] for r in rows)),
        "by_match": dict(Counter(r["main_database_match"] for r in rows)),
        "by_group": dict(Counter(r["selection_group"] for r in rows if "selection_group" in r)),
        "features_by_band": {
            band: dict(Counter(f for r in rows if r["band"] == band for f in r["features"]))
            for band in ("<=15", "16-25", "26-35", "36-45", "46+")
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=("all", "caspre", "dft8k"), default="all")
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")
    OUT.mkdir(parents=True, exist_ok=True)
    split = {json.loads(line)["connectivity_key"] for line in SPLIT.read_text(encoding="utf-8").splitlines()}
    with GYM.open(encoding="utf-8-sig", newline="") as f:
        gym = {r["connectivity"] for r in csv.DictReader(f)}
    report_path = OUT / "screening_summary.json"
    report = json.loads(report_path.read_text(encoding="utf-8")) if report_path.is_file() else {}
    base_groups = ["sulfur", "aryl_amide", "halogen_heteroaromatic", "heteroaromatic", "polycyclic", "nonprotonated_c"]
    groups = {band: base_groups for band in ("16-25", "26-35")}
    big = ["sulfur", "fluorinated", "chlorinated", "nonprotonated_c", "heteroaromatic", "polycyclic"]
    groups.update({"36-45": big, "46+": big})
    with sqlite3.connect(DB.as_uri() + "?mode=ro", uri=True) as db:
        db.execute("PRAGMA query_only=ON")
        db.execute("BEGIN")
        holdout = None
        if args.source in {"all", "caspre"} and CASPRE.exists():
            rows, summary = read_caspre()
            annotate_overlap(rows, db)
            ready = [
                r for r in rows
                if r["connectivity"] not in split | gym
                and r["structure_has_defined_stereo"]
                and r["neutral_single_component_supported_elements"]
                and r["heavy_atoms"] >= 16
            ]
            picked = select(ready, {"16-25": 20, "26-35": 40, "36-45": 25, "46+": 15}, groups)
            holdout = {r["connectivity"] for r in picked}
            write_jsonl(OUT / "caspre_experimental_13c_panel.jsonl", picked)
            report["caspre"] = {**summary, "validated": group_counts(rows), "eligible": group_counts(ready), "selected": group_counts(picked)}
            print("CASPRE", len(rows), len(ready), len(picked), flush=True)

        if args.source in {"all", "dft8k"} and all((DFT / f).exists() for f in ("DFT8K.csv.gz", "DFT.sdf.gz")):
            if holdout is None:
                panel_path = OUT / "caspre_experimental_13c_panel.jsonl"
                holdout = {json.loads(line)["connectivity"] for line in panel_path.read_text(encoding="utf-8").splitlines()} if panel_path.is_file() else set()
            rows, summary = read_dft()
            annotate_overlap(rows, db)
            excluded = gym | split | holdout
            ready = [
                r for r in rows if r["connectivity"] not in excluded
                and r["structure_has_defined_stereo"] and r["neutral_single_component_supported_elements"]
                and 16 <= r["heavy_atoms"] <= 45
            ]
            picked = select(ready, {"16-25": 45, "26-35": 40, "36-45": 15}, groups)
            write_jsonl(OUT / "cascade_dft8k_training_pilot.jsonl", picked)
            report["dft8k"] = {**summary, "validated": group_counts(rows), "eligible": group_counts(ready), "selected": group_counts(picked)}
            assert not ({r["connectivity"] for r in picked} & excluded)
            print("DFT8K", len(rows), len(ready), len(picked), flush=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
