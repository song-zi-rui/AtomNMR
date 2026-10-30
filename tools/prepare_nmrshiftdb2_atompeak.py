#!/usr/bin/env python
"""Prepare NMRShiftDB2 records for carbon/proton peak-environment alignment.

The converter keeps three kinds of information separate:

1. model inputs: molecular graph candidates and experimental peak attributes;
2. supervision only: NMRShiftDB2 atom assignments and target matrices;
3. provenance: source fields, spectrum conditions, and quality flags.

Unlike the first prototype converter, signals are never merged after rounding
their shifts. Repeated records are merged only when they belong to the same SDF
spectrum field and have the same original shift text and extra field.
"""

import argparse
import json
import pickle
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import lmdb
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from tools.shiftdb_atom_identity import canonical_identity, identity_for_current, source_peaks, POLICY
from tools.filter_nmrshiftdb2_experimental import classify


SPECTRUM_RE = re.compile(r"^Spectrum\s+(.+?)\s+(\d+)$")
EXTRA_RE = re.compile(r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+))?\s*([A-Za-z]+)?\s*$")
CARBON_CLASS = {"S": "C", "D": "CH", "T": "CH2", "Q": "CH3"}


def parse_indexed_property(value):
    """Parse strings such as '0:CDCl3 1:DMSO-D6' without splitting values."""
    matches = list(re.finditer(r"(?:^|\s)(\d+):", value or ""))
    parsed = {}
    for pos, match in enumerate(matches):
        end = matches[pos + 1].start() if pos + 1 < len(matches) else len(value)
        parsed[match.group(1)] = value[match.end() : end].strip()
    return parsed


def indexed_metadata(mol, property_name, spectrum_id):
    if not mol.HasProp(property_name):
        return None
    return parse_indexed_property(mol.GetProp(property_name)).get(str(spectrum_id))


def spectrum_fields(mol, nucleus):
    fields = []
    for name in mol.GetPropNames():
        match = SPECTRUM_RE.match(name)
        if match and match.group(1) == nucleus:
            fields.append((int(match.group(2)), name))
    return sorted(fields)


def signal_count(mol, field):
    return sum(bool(part.strip()) for part in mol.GetProp(field).split("|"))


def spectrum_meta(mol, field):
    match = SPECTRUM_RE.match(field)
    sid = match.group(2)
    return {
        "spectrum_id": sid,
        "source_field": field,
        "solvent": indexed_metadata(mol, "Solvent", sid),
        "temperature_k": indexed_metadata(mol, "Temperature [K]", sid),
        "field_strength_mhz": indexed_metadata(mol, "Field Strength [MHz]", sid),
        "assignment_method": indexed_metadata(mol, "Assignment Method", sid),
        "program": indexed_metadata(mol, "Program", sid),
        "raw_signal_records": signal_count(mol, field),
    }


def normalise_condition(value):
    if value is None:
        return None
    value = re.sub(r"\s+", " ", value.strip().lower())
    return None if value in {"", "unreported", "unknown"} else value


def pair_score(c_meta, h_meta):
    score = 0.0
    for key, weight in (("solvent", 6.0), ("temperature_k", 2.0)):
        left = normalise_condition(c_meta.get(key))
        right = normalise_condition(h_meta.get(key))
        if left and right:
            score += weight if left == right else -weight
    score += min(c_meta["raw_signal_records"], 20) * 0.01
    score += min(h_meta["raw_signal_records"], 20) * 0.01
    return score


def choose_spectra(mol, use_cnmr=True, use_hnmr=True):
    carbon = [(field, spectrum_meta(mol, field)) for _, field in spectrum_fields(mol, "13C")]
    proton = [(field, spectrum_meta(mol, field)) for _, field in spectrum_fields(mol, "1H")]
    if use_cnmr and use_hnmr and carbon and proton:
        candidates = [(pair_score(cm, hm), cf, cm, hf, hm) for cf, cm in carbon for hf, hm in proton]
        _, c_field, c_meta, h_field, h_meta = max(candidates, key=lambda row: row[0])
        matched = all(
            normalise_condition(c_meta.get(key)) == normalise_condition(h_meta.get(key))
            for key in ("solvent",)
            if normalise_condition(c_meta.get(key)) and normalise_condition(h_meta.get(key))
        )
        return (c_field, c_meta), (h_field, h_meta), matched
    c_selected = carbon[0] if use_cnmr and carbon else (None, None)
    h_selected = proton[0] if use_hnmr and proton else (None, None)
    return c_selected, h_selected, None


def parse_extra(raw, nucleus):
    match = EXTRA_RE.match(raw or "")
    numeric = float(match.group(1)) if match and match.group(1) else None
    code = match.group(2).upper() if match and match.group(2) else None
    parsed = {"raw": raw, "reported_uncertainty": numeric, "database_code": code}
    if nucleus == "13C":
        parsed["carbon_hydrogen_class"] = CARBON_CLASS.get(code)
    else:
        parsed["reported_multiplicity"] = code.lower() if code else None
    return parsed


def parse_signal(raw):
    parts = raw.split(";")
    if len(parts) < 3:
        return None
    shift_text = parts[0].strip()
    try:
        shift = float(shift_text)
    except ValueError:
        return None
    atom_indices = []
    for token in re.split(r"[, _]+", parts[-1].strip()):
        try:
            atom_indices.append(int(token))
        except ValueError:
            continue
    return shift, shift_text, parts[1].strip(), atom_indices


def parse_peak_field(mol, nucleus, field, *, identity=None, record_conflict=False):
    if field is None:
        return [], 0
    raw_peaks, labels, validation = source_peaks(mol.GetProp(field), nucleus,
        identity or identity_for_current(mol), record_conflict=record_conflict)
    skipped = sum(bool(f['issues']) for p in raw_peaks for f in p['source_assignment_fragments'])
    peaks = []
    prefix = "C" if nucleus == "13C" else "H"
    for index, original in enumerate(sorted(raw_peaks, key=lambda x: x["shift_ppm"], reverse=True)):
        item = {'shift_ppm': original['shift_ppm'],
                'shift_text': original['source_assignment_fragments'][0]['shift_text'],
                'extra': parse_extra(original['source_signal_code'], nucleus),
                'assigned_atom_indices': labels.get(original['peak_id'], {}).get('atom_indices', []),
                'assignment_record_count': len(original['source_assignment_fragments']),
                'source_assignment_fragments': original['source_assignment_fragments'],
                'source_assignment_validation': validation}
        item["peak_id"] = f"{prefix}p{index}"
        item["nucleus"] = nucleus
        item["assigned_atom_indices"] = sorted(item["assigned_atom_indices"])
        item["experimental_mass"] = 1.0
        item["mass_source"] = "unit_peak_mass_no_reliable_integral"
        item["input_confidence"] = 1.0 if nucleus == "13C" else 0.7
        peaks.append(item)
    return peaks, skipped


def fragment_smiles(mol, center_idx, radius):
    bonds = set()
    frontier = {center_idx}
    atoms = {center_idx}
    for _ in range(radius):
        next_frontier = set()
        for idx in frontier:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                bonds.add(bond.GetIdx())
                other = bond.GetOtherAtomIdx(idx)
                atoms.add(other)
                next_frontier.add(other)
        frontier = next_frontier
    try:
        return Chem.MolFragmentToSmiles(
            mol,
            atomsToUse=sorted(atoms),
            bondsToUse=sorted(bonds),
            rootedAtAtom=center_idx,
            canonical=True,
            isomericSmiles=True,
        )
    except Exception:
        return ""


def atom_feature(mol, atom, symmetry_class):
    neighbors = []
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        neighbors.append(
            {
                "atom_index": other.GetIdx(),
                "element": other.GetSymbol(),
                "bond_type": str(bond.GetBondType()),
                "conjugated": bool(bond.GetIsConjugated()),
            }
        )
    return {
        "atom_index": atom.GetIdx(),
        "element": atom.GetSymbol(),
        "atomic_number": atom.GetAtomicNum(),
        "formal_charge": atom.GetFormalCharge(),
        "degree": atom.GetDegree(),
        "total_h_count": int(atom.GetTotalNumHs(includeNeighbors=True)),
        "hybridization": str(atom.GetHybridization()),
        "aromatic": bool(atom.GetIsAromatic()),
        "in_ring": bool(atom.IsInRing()),
        "chirality": str(atom.GetChiralTag()),
        "symmetry_class": int(symmetry_class),
        "neighbors": sorted(neighbors, key=lambda row: row["atom_index"]),
        "radius1_smiles": fragment_smiles(mol, atom.GetIdx(), 1),
        "radius2_smiles": fragment_smiles(mol, atom.GetIdx(), 2),
    }


def molecular_graph(mol):
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=True))
    nodes = [atom_feature(mol, atom, ranks[atom.GetIdx()]) for atom in mol.GetAtoms()]
    edges = []
    for bond in mol.GetBonds():
        edges.append(
            {
                "source": bond.GetBeginAtomIdx(),
                "target": bond.GetEndAtomIdx(),
                "bond_type": str(bond.GetBondType()),
                "aromatic": bool(bond.GetIsAromatic()),
                "conjugated": bool(bond.GetIsConjugated()),
                "in_ring": bool(bond.IsInRing()),
            }
        )
    return {"nodes": nodes, "edges": edges, "node_count": len(nodes), "edge_count": len(edges)}, ranks


def build_candidate_groups(mol, ranks, nucleus):
    by_rank = defaultdict(list)
    for atom in mol.GetAtoms():
        if nucleus == "13C" and atom.GetSymbol() == "C":
            by_rank[ranks[atom.GetIdx()]].append(atom.GetIdx())
        elif nucleus == "1H" and atom.GetSymbol() != "H" and atom.GetTotalNumHs(includeNeighbors=True) > 0:
            by_rank[ranks[atom.GetIdx()]].append(atom.GetIdx())

    groups = []
    prefix = "Cg" if nucleus == "13C" else "Hg"
    for index, (_, atom_indices) in enumerate(sorted(by_rank.items(), key=lambda row: min(row[1]))):
        atoms = [mol.GetAtomWithIdx(idx) for idx in atom_indices]
        group = {
            "group_id": f"{prefix}{index}",
            "nucleus": nucleus,
            "center_atom_indices": atom_indices,
            "construction": "rdkit_graph_symmetry_class",
            "graph_equivalent_count": len(atom_indices),
            "representative_atom_index": atom_indices[0],
            "representative_environment": {
                "element": atoms[0].GetSymbol(),
                "hybridization": str(atoms[0].GetHybridization()),
                "aromatic": bool(atoms[0].GetIsAromatic()),
                "formal_charge": atoms[0].GetFormalCharge(),
                "radius1_smiles": fragment_smiles(mol, atom_indices[0], 1),
                "radius2_smiles": fragment_smiles(mol, atom_indices[0], 2),
            },
        }
        if nucleus == "13C":
            group["transport_mass"] = float(len(atom_indices))
            group["mass_source"] = "graph_equivalent_carbon_count"
        else:
            h_count = sum(atom.GetTotalNumHs(includeNeighbors=True) for atom in atoms)
            group["observable_h_count"] = int(h_count)
            group["transport_mass"] = float(max(h_count, 1))
            group["mass_source"] = "rdkit_attached_h_count"
        groups.append(group)
    return groups


def assignment_targets(peaks, groups):
    binary = np.zeros((len(peaks), len(groups)), dtype=np.float32)
    fractional = np.zeros_like(binary)
    for p_idx, peak in enumerate(peaks):
        assigned = set(peak["assigned_atom_indices"])
        for g_idx, group in enumerate(groups):
            members = set(group["center_atom_indices"])
            overlap = len(assigned & members)
            if overlap:
                binary[p_idx, g_idx] = 1.0
                fractional[p_idx, g_idx] = overlap / max(len(assigned), 1)
    return binary, fractional


def alignment_branch(peaks, groups):
    target, target_fraction = assignment_targets(peaks, groups)
    return {
        "peak_ids": [peak["peak_id"] for peak in peaks],
        "group_ids": [group["group_id"] for group in groups],
        "peak_mass": [peak["experimental_mass"] for peak in peaks],
        "group_mass": [group["transport_mass"] for group in groups],
        "graph_existence_mask": np.ones((len(peaks), len(groups)), dtype=np.float32),
        "target_matrix": target,
        "target_fraction_matrix": target_fraction,
        "target_source": "nmrshiftdb2_atom_indices_supervision_only",
        "assignment_valid_by_peak": [bool(p['assigned_atom_indices']) for p in peaks],
    }


def legacy_fields(carbon_peaks, proton_peaks, c_groups, h_groups, c_target, h_target):
    c_legacy = [
        {"delta (ppm)": p["shift_ppm"], "atom_indices": p["assigned_atom_indices"], "spectrum_id": p["spectrum_id"]}
        for p in carbon_peaks
    ]
    h_legacy = [
        {
            "centroid": p["shift_ppm"],
            "delta": p["shift_ppm"],
            "rangeMax": p["shift_ppm"],
            "rangeMin": p["shift_ppm"],
            "nH": 1,
            "nH_source": "unit_placeholder_not_experimental_integral",
            "category": p["extra"].get("reported_multiplicity") or "unk",
            "j_values": "_",
            "atom_indices": p["assigned_atom_indices"],
            "spectrum_id": p["spectrum_id"],
        }
        for p in proton_peaks
    ]
    groups = []
    for group in c_groups:
        groups.append({"element": "C", "nucleus": "13C", "atom_indices": group["center_atom_indices"], "count": group["transport_mass"]})
    for group in h_groups:
        groups.append({"element": "H", "nucleus": "1H", "atom_indices": group["center_atom_indices"], "count": group["transport_mass"]})
    combined = np.zeros((len(c_legacy) + len(h_legacy), len(groups)), dtype=np.float32)
    combined[: len(c_legacy), : len(c_groups)] = c_target
    combined[len(c_legacy) :, len(c_groups) :] = h_target
    return c_legacy or None, h_legacy or None, groups, combined


def molecule_record(mol, args):
    try:
        source = mol
        identity = canonical_identity(source, 'shiftdb:' + (source.GetProp('nmrshiftdb2 ID') if source.HasProp('nmrshiftdb2 ID') else 'record'), sanitize=True)
        mol = identity.processed
        # SDF properties contain references in SOURCE order, despite being
        # copied onto the processed molecule. Always parse through identity.
        for name in source.GetPropNames(): mol.SetProp(name, source.GetProp(name))
    except Exception as exc:
        return None, "atom_identity_mapping_failed" if 'parser_error:' in str(exc) else "sanitize_failed"

    c_selected, h_selected, condition_match = choose_spectra(mol, args.use_cnmr, args.use_hnmr)
    c_field, c_meta = c_selected
    h_field, h_meta = h_selected
    if args.use_cnmr and c_field is None:
        return None, "missing_13c"
    if args.use_hnmr and h_field is None:
        return None, "missing_1h"

    props = {n: source.GetProp(n) for n in source.GetPropNames()}
    record_conflict = any(source_peaks(source.GetProp(f), '13C', identity)[2]['status'] == 'source_assignment_inconsistency'
                          for sid,f in spectrum_fields(source, '13C') if classify(props, str(sid), '13C')[0] == 'experimental')
    carbon_peaks, skipped_c = parse_peak_field(mol, "13C", c_field, identity=identity)
    proton_peaks, skipped_h = parse_peak_field(mol, "1H", h_field, identity=identity, record_conflict=record_conflict)
    if any(p['source_assignment_validation']['status'] != 'identity_validated' for p in carbon_peaks + proton_peaks):
        # Legacy alignment readers do not all implement per-peak masks.
        # Retain the complete spectrum in the audit/annotation DB, but never
        # serialize corrupt or uncertain supervision into their training LMDB.
        return None, "source_assignment_identity_requires_check"
    if args.use_cnmr and not carbon_peaks:
        return None, "no_valid_13c_assignments"
    if args.use_hnmr and not proton_peaks:
        return None, "no_valid_1h_assignments"

    graph, ranks = molecular_graph(mol)
    c_groups = build_candidate_groups(mol, ranks, "13C") if carbon_peaks else []
    h_groups = build_candidate_groups(mol, ranks, "1H") if proton_peaks else []
    c_alignment = alignment_branch(carbon_peaks, c_groups)
    h_alignment = alignment_branch(proton_peaks, h_groups)
    for peak in carbon_peaks:
        peak["spectrum_id"] = c_meta["spectrum_id"]
    for peak in proton_peaks:
        peak["spectrum_id"] = h_meta["spectrum_id"]

    c_legacy, h_legacy, legacy_groups, legacy_target = legacy_fields(
        carbon_peaks,
        proton_peaks,
        c_groups,
        h_groups,
        c_alignment["target_matrix"],
        h_alignment["target_matrix"],
    )
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    record = {
        "schema_version": "nmrpeak_atompeak_v2",
        "identifiers": {
            "nmrshiftdb2_id": mol.GetProp("nmrshiftdb2 ID") if mol.HasProp("nmrshiftdb2 ID") else "",
            "inchi": mol.GetProp("INChI") if mol.HasProp("INChI") else "",
            "inchi_key": mol.GetProp("INChI key") if mol.HasProp("INChI key") else "",
        },
        "molecule": {
            "smiles": smiles,
            "canonical_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
            "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
            "murcko_scaffold": MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=True),
        },
        "molecular_graph": graph,
        "spectra": {
            "carbon": {"metadata": c_meta, "peaks": carbon_peaks} if c_meta else None,
            "proton": {"metadata": h_meta, "peaks": proton_peaks} if h_meta else None,
            "condition_pair_match": condition_match,
        },
        "candidate_environments": {"carbon_groups": c_groups, "proton_environments": h_groups},
        "alignment_supervision": {"carbon": c_alignment, "proton": h_alignment},
        "quality": {
            "skipped_13c_signal_records": skipped_c,
            "skipped_1h_signal_records": skipped_h,
            "has_carbon_supervision": bool(c_alignment["target_matrix"].sum()),
            "has_proton_supervision": bool(h_alignment["target_matrix"].sum()),
            "experimental_1h_integral_available": False,
            "experimental_1h_j_available": False,
            "notes": [
                "graph symmetry groups are computed candidates, not database equivalence labels",
                "NMRShiftDB2 atom indices are stored as supervision and are not peak input features",
                "proton peak mass is unit mass because reliable experimental integrals are unavailable",
            ],
        },
        "provenance": {
            "atom_identity": identity.to_dict(),
            "assignment_identity_policy": POLICY,
            "original_molblock": Chem.MolToMolBlock(source),
            "source": "nmrshiftdb2withsignals.sd",
            "signal_merge_rule": "same_spectrum_field_and_exact_shift_text_and_exact_extra",
            "shift_rounding_used": False,
            "candidate_group_rule": "RDKit CanonicalRankAtoms breakTies=False includeChirality=True",
        },
        # Backward-compatible fields for the current NMRPeak prototype.
        "smiles": smiles,
        "smiles_canonical": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
        "atoms": [atom.GetSymbol() for atom in mol.GetAtoms()],
        "coordinates": np.zeros((mol.GetNumAtoms(), 3), dtype=np.float32),
        "c_nmr_peaks": c_legacy,
        "h_nmr_peaks": h_legacy,
        "atom_groups": legacy_groups,
        "assignment_matrix": legacy_target,
    }
    return record, "ok"


def json_ready(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_ready(item) for item in value]
    return value


def split_records(records, train_ratio, valid_ratio, seed, split_mode):
    rng = random.Random(seed)
    if split_mode == "random":
        units = [[record] for record in records]
    else:
        grouped = defaultdict(list)
        for record in records:
            key = record["molecule"]["murcko_scaffold"] or record["molecule"]["canonical_smiles"]
            grouped[key].append(record)
        units = list(grouped.values())
    rng.shuffle(units)
    target_train = len(records) * train_ratio
    target_valid = len(records) * valid_ratio
    splits = {"train": [], "valid": [], "test": []}
    for unit in units:
        if len(splits["train"]) < target_train:
            split = "train"
        elif len(splits["valid"]) < target_valid:
            split = "valid"
        else:
            split = "test"
        splits[split].extend(unit)
    return splits


def write_lmdb(path, records, map_size):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    env = lmdb.open(str(path), subdir=False, map_size=map_size, lock=True)
    with env.begin(write=True) as txn:
        for index, record in enumerate(records):
            txn.put(str(index).encode("ascii"), pickle.dumps(record, protocol=pickle.HIGHEST_PROTOCOL))
    env.sync()
    env.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdf", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--max-records", type=int)
    parser.add_argument("--sample-inchi-key", help="Also write the requested record as sample_record.json")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--train-ratio", type=float, default=0.9)
    parser.add_argument("--valid-ratio", type=float, default=0.05)
    parser.add_argument("--split-mode", choices=("scaffold", "random"), default="scaffold")
    parser.add_argument("--map-size-gb", type=float, default=16.0)
    parser.add_argument("--use-cnmr", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--use-hnmr", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()

    RDLogger.DisableLog("rdApp.*")
    supplier = Chem.SDMolSupplier(args.sdf, removeHs=False, sanitize=False)
    records = []
    reasons = defaultdict(int)
    sample = None
    for mol_index, mol in enumerate(supplier):
        if args.limit is not None and mol_index >= args.limit:
            break
        if mol is None:
            reasons["mol_none"] += 1
            continue
        record, reason = molecule_record(mol, args)
        reasons[reason] += 1
        if record is None:
            continue
        records.append(record)
        if args.sample_inchi_key and record["identifiers"]["inchi_key"] == args.sample_inchi_key:
            sample = record
        if args.max_records is not None and len(records) >= args.max_records and not args.sample_inchi_key:
            break

    if not records:
        raise RuntimeError(f"No records converted: {dict(reasons)}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    splits = split_records(records, args.train_ratio, args.valid_ratio, args.seed, args.split_mode)
    map_size = int(args.map_size_gb * 1024**3)
    for split, items in splits.items():
        write_lmdb(out_dir / "lmdb_dataset" / f"{split}.lmdb", items, map_size)

    if sample is None:
        sample = records[0]
    with open(out_dir / "sample_record.json", "w", encoding="utf-8") as handle:
        json.dump(json_ready(sample), handle, ensure_ascii=False, indent=2)

    summary = {
        "schema_version": "nmrpeak_atompeak_v2",
        "total_converted": len(records),
        "split_mode": args.split_mode,
        "split_sizes": {key: len(value) for key, value in splits.items()},
        "conversion_reasons": dict(reasons),
        "sample_inchi_key": sample["identifiers"]["inchi_key"],
        "sample_carbon_shape": list(sample["alignment_supervision"]["carbon"]["target_matrix"].shape),
        "sample_proton_shape": list(sample["alignment_supervision"]["proton"]["target_matrix"].shape),
    }
    with open(out_dir / "conversion_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
