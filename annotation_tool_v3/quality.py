"""Conservative admission checks; passing is not proof of chemical correctness."""
from __future__ import annotations

import json
import math
import re
from functools import lru_cache

import numpy as np
from rdkit import Chem
from scipy.optimize import linear_sum_assignment

POLICY = "nmrexp-link-v2"
THRESHOLDS = {
    "13C": {"tolerance": 1.0, "coverage": .90, "median": .30, "p95": 1.0},
    "1H": {"tolerance": .10, "coverage": .90, "median": .03, "p95": .10},
}


def peak_match(local, remote, tolerance):
    """Maximize in-tolerance cardinality, then minimize error (no peak reuse)."""
    if not local or not remote:
        return {"f1": 0., "local_coverage": 0., "remote_coverage": 0.,
                "median": None, "p95": None, "mapping": []}
    a, b = np.asarray(local, dtype=float), np.asarray(remote, dtype=float)
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("nonfinite shift")
    ai, bi = np.argsort(a, kind="stable"), np.argsort(b, kind="stable")
    if len(a) == len(b) and np.all(abs(a[ai] - b[bi]) <= tolerance + 1e-8):
        left, right = ai, bi
    else:
        distance = abs(a[:, None] - b[None, :])
        valid = distance <= tolerance + 1e-8
        cost = np.where(valid, distance / max(tolerance, 1e-12), min(len(a), len(b)) + 1.)
        left, right = linear_sum_assignment(cost)
        keep = valid[left, right]
        left, right = left[keep], right[keep]
    errors = abs(a[left] - b[right])
    n = len(errors)
    return {"f1": 2*n/(len(a)+len(b)), "local_coverage": n/len(a), "remote_coverage": n/len(b),
            "median": float(np.median(errors)) if n else None,
            "p95": float(np.percentile(errors, 95)) if n else None,
            "mapping": [[int(i), int(j), float(e)] for i, j, e in zip(left, right, errors)]}


@lru_cache(maxsize=8192)
def structure_info(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("invalid_structure")
    carbons = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    return {"canonical": Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True),
            "components": len(Chem.GetMolFrags(mol)), "carbons": len(carbons),
            "carbon_h": sum(a.GetTotalNumHs(includeNeighbors=True) for a in carbons),
            "total_h": sum(a.GetTotalNumHs(includeNeighbors=True) for a in mol.GetAtoms() if a.GetAtomicNum() != 1),
            "ch2": sum(a.GetTotalNumHs(includeNeighbors=True) == 2 for a in carbons),
            "ch3": sum(a.GetTotalNumHs(includeNeighbors=True) == 3 for a in carbons)}


def bbox(value):
    vals = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(value or ""))
    if len(vals) != 4:
        return None
    vals = tuple(float(x) for x in vals)
    return vals if all(math.isfinite(x) and 0 <= x <= 1 for x in vals) and vals[0] < vals[2] and vals[1] < vals[3] else None


def source_anchor(h, c):
    reasons = []
    if not h.get("Filename") or h.get("Filename") != c.get("Filename"):
        reasons.append("source_doi_conflict_or_missing")
    for field in ("Page_in_file_mol", "Page_in_file_para"):
        if h.get(field) is None or c.get(field) is None or h.get(field) != c.get(field):
            reasons.append("source_" + field + "_conflict_or_missing")
    if h.get("Page_in_file_mol") != h.get("Page_in_file_para"):
        reasons.append("structure_text_on_different_pages")
    for field in ("Location_in_page_mol", "Location_in_page_para"):
        a, b = bbox(h.get(field)), bbox(c.get(field))
        if a is None or b is None:
            reasons.append("source_" + field + "_missing")
        elif max(abs(x-y) for x,y in zip(a,b)) > .002:
            reasons.append("source_" + field + "_conflict")
    return reasons


def evaluate(smiles, h, c, solvent, references=()):
    """h/c are (source row, raw metadata, parsed peaks, retrieval F1, error)."""
    reasons, warnings, metrics = [], [], {}
    info = structure_info(smiles)
    for item in (h, c):
        if structure_info(item[1]["SMILES"])["canonical"] != info["canonical"]:
            reasons.append("source_structure_mismatch")
        if not solvent or solvent.lower() == "mixed" or item[0]["solvent"] != solvent:
            reasons.append("source_solvent_conflict_or_missing")
    if info["components"] != 1:
        reasons.append("multicomponent_structure")
    reasons.extend(source_anchor(h[1], c[1]))
    hf = re.search(r"\d+(?:\.\d+)?", str(h[1].get("NMR_frequency") or ""))
    cf = re.search(r"\d+(?:\.\d+)?", str(c[1].get("NMR_frequency") or ""))
    if hf and cf and float(cf[0]) > 0 and not 3.75 <= float(hf[0])/float(cf[0]) <= 4.10:
        reasons.append("source_field_ratio_conflict")
    hp, cp = h[2], c[2]
    if not hp or not cp:
        reasons.append("empty_spectrum")
    if len(cp) > info["carbons"]:
        reasons.append("carbon_peak_count_exceeds_carbon_atoms")
    integrals = [p.get("reported_nH") for p in hp]
    if integrals and all(x is not None and x > 0 and abs(x-round(x)) < 1e-6 for x in integrals):
        total = sum(integrals)
        metrics["hydrogen_accounting"] = {"reported": total, "carbon_bound": info["carbon_h"], "formula_total": info["total_h"]}
        if total > info["total_h"] + .01 or total < info["carbon_h"] - .01:
            reasons.append("hydrogen_integral_structure_conflict")
    else:
        warnings.append("hydrogen_integrals_incomplete_or_noninteger")
    # A targeted, explainable warning, not a universal multiplicity predictor.
    triplets = [p for p in hp if str(p.get("multiplicity", "")).lower() == "t"
                and p.get("reported_nH") == 2 and 1.5 <= p["shift_ppm"] <= 4.5]
    if info["ch2"] == 0 and info["ch3"] and len(triplets) >= 2:
        reasons.append("two_2H_triplets_but_no_CH2_in_structure")
    by_nucleus = {}
    for ref in references:
        by_nucleus.setdefault(ref["nucleus"], []).append(ref)
    for nuc, refs in by_nucleus.items():
        policy = THRESHOLDS[nuc]
        remote = hp if nuc == "1H" else cp
        checks = []
        for ref in refs:
            if not ref["solvent"] or ref["solvent"] != solvent:
                continue
            peaks = json.loads(ref["peaks_json"])
            match = peak_match([p["shift_ppm"] for p in peaks], [p["shift_ppm"] for p in remote], policy["tolerance"])
            match["spectrum_id"] = ref["id"]
            match["passed"] = (min(match["local_coverage"], match["remote_coverage"]) >= policy["coverage"]
                               and match["median"] is not None and match["median"] <= policy["median"] + 1e-8
                               and match["p95"] <= policy["p95"] + 1e-8)
            checks.append(match)
        metrics[nuc] = checks
        if not any(x["passed"] for x in checks):
            reasons.append("shared_" + nuc + "_spectrum_conflict" if checks else "shared_" + nuc + "_no_same_known_solvent")
    return {"policy": POLICY, "decision": "quarantined" if reasons else "passed",
            "reasons": sorted(set(reasons)), "warnings": warnings, "metrics": metrics,
            "source_ids": {"1H": h[0]["id"], "13C": c[0]["id"]},
            "source_document_verified": False, "source_anchor_check": "extracted_metadata_only",
            "temperature_match": "unknown", "same_physical_sample": "not_established"}


def references(c, mid):
    return c.execute("SELECT id,nucleus,solvent,temperature_k,peaks_json FROM spectra WHERE molecule_id=? AND source='nmrshiftdb2'", (mid,)).fetchall() if mid else []


def write_result(c, mid, candidate, result, spectrum_ids=()):
    from . import core as db
    c.execute("INSERT INTO candidate_quality VALUES(?,?,?,?,?) ON CONFLICT(smiles) DO UPDATE SET policy=excluded.policy,decision=excluded.decision,details_json=excluded.details_json,updated_at=excluded.updated_at",
              (candidate["smiles"], POLICY, result["decision"], db.dumps(result), db.now()))
    changed = False
    for sid in spectrum_ids:
        old = c.execute("SELECT decision FROM spectrum_quality WHERE spectrum_id=?", (sid,)).fetchone()
        was_active = not old or old[0] == "passed"
        c.execute("INSERT INTO spectrum_quality VALUES(?,?,?,?,?) ON CONFLICT(spectrum_id) DO UPDATE SET policy=excluded.policy,decision=excluded.decision,details_json=excluded.details_json,updated_at=excluded.updated_at",
                  (sid, POLICY, result["decision"], db.dumps(result), db.now()))
        changed |= was_active != (result["decision"] == "passed")
    if mid and changed:
        c.execute("UPDATE molecules SET revision=revision+1,task_done=0 WHERE id=?", (mid,))
        db.refresh_molecule(c, mid)
