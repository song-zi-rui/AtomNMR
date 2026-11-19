"""Audit the published IBM computational NMR Parquet against AtomNMR."""

from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import sqlite3

import numpy as np
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdDetermineBonds

from screen_balanced_nmr_sources import (
    DB, GYM, ROOT, SPLIT, annotate_overlap, feature_record, write_jsonl,
)

SOURCE = ROOT / "data/nmr_atom_assignment_sources/ibm_computational_nmr/NMR_data.parquet"
OUT = ROOT / "outputs/ibm_computational_nmr_audit_20261119"
MD5 = "5ba7fa1cfa1ea03a8c46dab8ed3c3501"
SHIELD = re.compile(r"^\s*([CH])\s+\d+\s*/\s*(\d+)\s+SHIELD\s*:\s*(.*)$", re.M)
NUCLEI = {"h": "1H", "c": "13C"}
SIZE_BANDS = ("<=15", "16-25", "26-35", "36-45", "46+")


def finite_array(values, shape):
    array = np.asarray(values, dtype=float)
    return array.shape == shape and bool(np.isfinite(array).all())


def frame_labels(row):
    """Check source indices/values against the original CPMD output."""
    atoms = row["atoms"]
    ids = {n: [i for i, el in enumerate(atoms) if el == n.upper()] for n in NUCLEI}
    if not all(ids.values()):
        raise ValueError("missing_nucleus")
    frames = [(k, frame) for k, frame in row["frames"].items() if frame is not None]
    if len(frames) != row["number_of_frames"] or len(frames) != len(set(row["frame_ids"])) or len(frames) < 1:
        raise ValueError("frame_count_inconsistent")
    per_nucleus = defaultdict(list)
    expected = {i: n.upper() for n in ids for i in ids[n]}
    for _, frame in frames:
        shielding = {}
        for atom, position, values in SHIELD.findall(frame["nmr_cpmd_text"]):
            index = int(position) - 1
            if index in shielding:
                raise ValueError("duplicate_index_in_raw_output")
            shielding[index] = (atom, float(values.split()[-1]))
        if {i: x[0] for i, x in shielding.items()} != expected:
            raise ValueError("raw_atom_index_or_element_mismatch")
        for n, indices in ids.items():
            values = frame[f"{n}_nmr_peaks_unsorted"]
            if not finite_array(values, (len(indices),)):
                raise ValueError("missing_or_nonfinite_frame_shifts")
            if not np.allclose(values, [shielding[i][1] for i in indices], rtol=0, atol=1e-8):
                raise ValueError("raw_shift_order_or_value_mismatch")
            per_nucleus[n].append(values)
    labels = {}
    grouping_ok = {}
    published_std_ok = {}
    group_indices = {}
    for n, values in per_nucleus.items():
        array = np.asarray(values)
        means = array.mean(axis=0)
        std = array.std(axis=0)
        reported = row["averaged_frames"][f"{n}_nmr_peaks_frame_ave"]
        if not finite_array(reported, (len(ids[n]),)) or not np.allclose(means, reported, atol=1e-8, rtol=0):
            raise ValueError("frame_mean_mismatch")
        labels[n] = [
            {"atom_index": i, "ppm": float(mean), "frame_std_ppm": float(s)}
            for i, mean, s in zip(ids[n], means, std)
        ]
        published_std_ok[n] = bool(np.allclose(std, row[f"{n}_nmr_std"], rtol=0, atol=1e-8))
        groups = row[f"{n}_atoms_group_index"]
        indices = [i for g in groups for i in g]
        mean_by_atom = dict(zip(ids[n], means))
        grouped = row["averaged_frames"][f"{n}_nmr_peaks_grouped_frame_ave"]
        grouping_ok[n] = (
            all(groups) and sorted(indices) == ids[n]
            and finite_array(grouped, (len(groups),))
            and bool(np.allclose(grouped, [np.mean([mean_by_atom[i] for i in g]) for g in groups], rtol=0, atol=1e-8))
        )
        group_indices[n] = groups
    return labels, grouping_ok, published_std_ok, group_indices


def geometry_atom_order(row, template, *, allow_unspecified_stereo=False):
    """Recover the IBM XYZ order via a complete, element-matched graph isomorphism."""
    atoms = row["atoms"]
    xyz = json.loads(row["xyz"])
    explicit = Chem.AddHs(template)
    if Counter(atoms) != Counter(a.GetSymbol() for a in explicit.GetAtoms()):
        return "element_count_conflict", None
    if not finite_array(xyz, (len(atoms), 3)):
        return "invalid_geometry", None
    xyz_block = f"{len(atoms)}\nIBM atom order\n" + "".join(
        el + " " + " ".join(map(str, pos)) + "\n" for el, pos in zip(atoms, xyz)
    )
    params = Chem.AdjustQueryParameters()
    params.makeBondsGeneric = True
    query = Chem.AdjustQueryProperties(explicit, params)
    key = Chem.MolToInchiKey(template)
    found_graph = False
    for factor in (None, 1.15, 1.20, 1.25, 1.30):
        geom = Chem.MolFromXYZBlock(xyz_block)
        if geom is None:
            return "invalid_geometry", None
        if factor is None:
            rdDetermineBonds.DetermineConnectivity(geom)
        else:
            rdDetermineBonds.DetermineConnectivity(geom, covFactor=factor, useVdw=True)
        if geom.GetNumBonds() != explicit.GetNumBonds():
            continue
        geom.UpdatePropertyCache(strict=False)
        matches = geom.GetSubstructMatches(query, uniquify=False, maxMatches=16)
        found_graph |= bool(matches)
        for mapping in matches:
            inverse = [0] * len(mapping)
            for i, j in enumerate(mapping):
                inverse[j] = i
            ordered = Chem.RenumberAtoms(explicit, inverse)
            ordered.AddConformer(Chem.Conformer(geom.GetConformer()), assignId=True)
            inferred = Chem.Mol(ordered)
            Chem.AssignStereochemistryFrom3D(inferred, replaceExistingTags=True)
            inferred_heavy = Chem.RemoveHs(inferred)
            if Chem.MolToInchiKey(inferred_heavy) == key:
                return "validated_default" if factor is None else "validated_alternate_radius", ordered
            if (
                allow_unspecified_stereo
                and inferred_heavy.HasSubstructMatch(template, useChirality=True)
                and Chem.MolToSmiles(inferred_heavy, isomericSmiles=False)
                == Chem.MolToSmiles(template, isomericSmiles=False)
            ):
                return "validated_declared_stereo_only", ordered
    return "geometry_stereo_difference" if found_graph else "geometry_graph_mismatch", None


def unkeyed_main_rows(connection):
    rowids = connection.execute(
        "SELECT rowid FROM molecules INDEXED BY molecule_inchikey "
        "WHERE inchikey IS NULL OR length(inchikey)!=27"
    ).fetchall()
    by_smiles = defaultdict(list)
    for rowid, in rowids:
        mid, smiles, queue = connection.execute(
            "SELECT id,smiles,queue FROM molecules WHERE rowid=?", (rowid,)
        ).fetchone()
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            by_smiles[Chem.MolToSmiles(mol, isomericSmiles=False)].append((mid, queue))
    return by_smiles, len(rowids)


def distribution(rows):
    return {
        "count": len(rows),
        "by_band": dict(Counter(r["band"] for r in rows)),
        "by_size": dict(sorted(Counter(r["heavy_atoms"] for r in rows).items())),
        "features": dict(sorted(Counter(k for r in rows for k in r["features"]).items())),
        "features_26_or_more": dict(sorted(Counter(
            k for r in rows if r["heavy_atoms"] >= 26 for k in r["features"]
        ).items())),
    }


def quantiles(values):
    if not values:
        return {}
    return {str(p): round(float(np.percentile(values, p)), 4) for p in (0, 1, 10, 50, 90, 99, 100)}


def main():
    RDLogger.DisableLog("rdApp.*")
    digest = hashlib.md5(SOURCE.read_bytes()).hexdigest()
    if digest != MD5:
        raise ValueError("Parquet checksum does not match the published file")
    OUT.mkdir(parents=True, exist_ok=True)
    rows, payloads, errors = [], {}, {}
    metrics = defaultdict(list)
    grouping = Counter()
    std_agree = Counter()
    seen_ids = set()
    for batch in pq.ParquetFile(SOURCE).iter_batches(batch_size=48):
        for row in batch.to_pylist():
            source_id = row["id"]
            if source_id in seen_ids:
                raise ValueError("Repeated numeric source identifier")
            seen_ids.add(source_id)
            mol = Chem.MolFromSmiles(row["smiles"])
            if mol is None:
                raise ValueError(f"Could not parse source structure {source_id}")
            info = {
                **feature_record(mol), "source_id": source_id, "source_smiles": row["smiles"],
                "molecular_weight": round(Descriptors.MolWt(mol), 3),
                "frame_count": row["number_of_frames"],
                "element_counts": dict(Counter(row["atoms"])),
            }
            rows.append(info)
            try:
                shifts, group_checks, std_checks, groups = frame_labels(row)
            except (ValueError, TypeError, KeyError, IndexError) as exc:
                errors[source_id] = str(exc)
                info["raw_shift_audit"] = "failed"
                continue
            info["raw_shift_audit"] = "passed"
            info["grouping_audit"] = group_checks
            grouping.update(key for key, ok in group_checks.items() if ok)
            std_agree.update(key for key, ok in std_checks.items() if ok)
            for n, values in shifts.items():
                info[f"{n}_atoms_with_shifts"] = len(values)
                metrics[f"{n}_ppm"].extend(v["ppm"] for v in values)
                metrics[f"{n}_frame_std_ppm"].extend(v["frame_std_ppm"] for v in values)
            try:
                status, ordered = geometry_atom_order(row, mol)
            except (RuntimeError, ValueError):
                status, ordered = "geometry_mapping_error", None
            info["geometry_mapping"] = status
            if ordered is None:
                continue
            for shift in shifts["h"]:
                h = ordered.GetAtomWithIdx(shift["atom_index"])
                neighbors = h.GetNeighbors()
                if len(neighbors) != 1:
                    raise ValueError(f"Unexpected hydrogen degree on source {source_id}")
                shift["attached_atom_index"] = neighbors[0].GetIdx()
            block = Chem.MolToMolBlock(ordered)
            restored = Chem.MolFromMolBlock(block, removeHs=False)
            if (
                restored is None or Chem.MolToInchiKey(Chem.RemoveHs(restored)) != info["inchikey"]
                or [a.GetSymbol() for a in restored.GetAtoms()] != row["atoms"]
            ):
                info["geometry_mapping"] = "molblock_round_trip_failed"
                continue
            payloads[source_id] = {
                "atom_index_basis": "zero_based_IBM_XYZ_atom_order",
                "source_order_molblock": block,
                "frame_ids": row["frame_ids"],
                "chemical_shift_reference": "TMS",
                "label_origin": "isolated_molecule_computational_DFT",
                "h_shifts": shifts["h"], "c_shifts": shifts["c"],
                "c_source_groups": groups["c"] if group_checks["c"] else [],
                "h_source_groups": groups["h"] if group_checks["h"] else [],
            }
    with sqlite3.connect(DB.as_uri() + "?mode=ro", uri=True) as connection:
        connection.execute("PRAGMA query_only=ON")
        connection.execute("BEGIN")
        annotate_overlap(rows, connection)
        unkeyed, unkeyed_count = unkeyed_main_rows(connection)
        main_count = connection.execute("SELECT count(*) FROM molecules").fetchone()[0]
        for row in rows:
            if row["main_database_match"] != "new" or not unkeyed:
                continue
            mol = Chem.MolFromSmiles(row["canonical_smiles"])
            hits = unkeyed.get(Chem.MolToSmiles(mol, isomericSmiles=False), [])
            if hits:
                row["main_database_match"] = "unkeyed_graph_hit"
                row["unkeyed_graph_hits"] = [{"molecule_id": mid, "queue": q} for mid, q in hits]
    export_layers = {json.loads(l)["connectivity_key"]: json.loads(l)["split"]
                     for l in SPLIT.read_text(encoding="utf-8").splitlines()}
    with GYM.open(encoding="utf-8-sig", newline="") as f:
        gym_layers = {r["connectivity"] for r in csv.DictReader(f)}
    dft_file = ROOT / "outputs/nmr_balance_selection_20261119/dft8k_dual_novel_atom_nmr.jsonl"
    dft_layers = {json.loads(l)["connectivity"] for l in dft_file.read_text(encoding="utf-8").splitlines()}
    chosen = []
    seen = set()
    reasons = Counter()
    for row in sorted(rows, key=lambda r: (-r["heavy_atoms"], r["source_id"])):
        row["previous_atompeak_split"] = export_layers.get(row["connectivity"])
        row["in_2dnmrgym_eval"] = row["connectivity"] in gym_layers
        row["in_prior_novel_dft8k"] = row["connectivity"] in dft_layers
        why = []
        if row["main_database_match"] != "new":
            why.append("already_in_main_by_connectivity")
        if row["heavy_atoms"] < 16:
            why.append("below_current_size_priority")
        if row["source_id"] not in payloads:
            why.append("atom_mapping_not_verified")
        if not row["structure_has_defined_stereo"]:
            why.append("unspecified_source_stereochemistry")
        if not row["neutral_single_component_supported_elements"]:
            why.append("outside_previous_first_batch_element_or_charge_scope")
        if row["previous_atompeak_split"] or row["in_2dnmrgym_eval"] or row["in_prior_novel_dft8k"]:
            why.append("overlap_with_previous_exports_or_eval")
        if not why:
            if row["connectivity"] in seen:
                why.append("source_connectivity_duplicate")
            else:
                seen.add(row["connectivity"])
                chosen.append({**row, **payloads[row["source_id"]]})
        row["screening_reasons"] = why
        reasons.update(why)
    write_jsonl(OUT / "ibm_all_records_audit.jsonl", rows)
    write_jsonl(OUT / "ibm_novel_dual_candidates.jsonl", chosen)
    result = {
        "source_bytes": SOURCE.stat().st_size, "source_md5": digest,
        "main_database_records": main_count, "main_missing_valid_key": unkeyed_count,
        "all": distribution(rows),
        "distinct_inchikeys": len({r["inchikey"] for r in rows}),
        "distinct_connectivities": len({r["connectivity"] for r in rows}),
        "molecular_weight_quantiles": quantiles([r["molecular_weight"] for r in rows]),
        "at_least_500_da": sum(r["molecular_weight"] >= 500 for r in rows),
        "frame_count_histogram": dict(sorted(Counter(r["frame_count"] for r in rows).items())),
        "frame_counts_total": sum(r["frame_count"] for r in rows),
        "raw_shift_errors": errors,
        "grouped_frame_mean_verified": dict(grouping),
        "published_std_equals_population_std": dict(std_agree),
        "geometry_mapping": dict(Counter(r.get("geometry_mapping") for r in rows)),
        "total_shifts_by_nucleus": {n: len(metrics[f"{n}_ppm"]) for n in NUCLEI},
        "atom_mean_ppm_quantiles": {n: quantiles(metrics[f"{n}_ppm"]) for n in NUCLEI},
        "atom_frame_sd_ppm_quantiles": {n: quantiles(metrics[f"{n}_frame_std_ppm"]) for n in NUCLEI},
        "match_counts_by_size": {
            band: dict(Counter(r["main_database_match"] for r in rows if r["band"] == band))
            for band in SIZE_BANDS
        },
        "prior_export_split_overlap": dict(Counter(r["previous_atompeak_split"] for r in rows if r["previous_atompeak_split"])),
        "overlaps_2dnmrgym": sum(r["in_2dnmrgym_eval"] for r in rows),
        "overlaps_novel_dft8k": sum(r["in_prior_novel_dft8k"] for r in rows),
        "db_queues_for_exact": dict(Counter(
            h["queue"] for r in rows for h in r["main_database_exact_hits"]
        )),
        "novel_before_quality": distribution([r for r in rows if r["main_database_match"] == "new"]),
        "novel_with_mapped_atoms_and_size_priority": distribution([
            r for r in rows if r["main_database_match"] == "new"
            and r["source_id"] in payloads and r["heavy_atoms"] >= 16
        ]),
        "selected": distribution(chosen),
        "screening_reason_counts_nonexclusive": dict(reasons),
    }
    (OUT / "summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print("Rows", len(rows), "with source labels", len(rows)-len(errors),
          "mapped", len(payloads), "novel screened", len(chosen), flush=True)


if __name__ == "__main__":
    main()
