#!/usr/bin/env python
"""Leakage-safe entry point for the enriched NMRShiftDB2 atom-peak schema.

The implementation lives in ``prepare_nmrshiftdb2_atompeak.py``. This entry
point enforces the final public schema: database atom indices are moved out of
the spectrum input objects and kept exclusively under alignment supervision.
"""

import copy
import importlib.util
import re
from pathlib import Path


def load_base():
    path = Path(__file__).with_name("prepare_nmrshiftdb2_atompeak.py")
    spec = importlib.util.spec_from_file_location("atompeak_base", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


base = load_base()
_base_molecule_record = base.molecule_record


def _unique_signal_count(mol, field):
    identities = set()
    for raw in filter(str.strip, mol.GetProp(field).split("|")):
        parsed = base.parse_signal(raw)
        if parsed is not None:
            _, shift_text, extra, _ = parsed
            identities.add((shift_text, extra))
    return len(identities)


def choose_spectra(mol, use_cnmr=True, use_hnmr=True):
    carbon = [(field, base.spectrum_meta(mol, field)) for _, field in base.spectrum_fields(mol, "13C")]
    proton = [(field, base.spectrum_meta(mol, field)) for _, field in base.spectrum_fields(mol, "1H")]
    if use_cnmr and use_hnmr and carbon and proton:
        candidates = []
        for c_field, c_meta in carbon:
            for h_field, h_meta in proton:
                score = base.pair_score(c_meta, h_meta)
                score += 0.05 * _unique_signal_count(mol, c_field)
                score += 0.05 * _unique_signal_count(mol, h_field)
                candidates.append((score, c_field, c_meta, h_field, h_meta))
        _, c_field, c_meta, h_field, h_meta = max(candidates, key=lambda row: row[0])
        c_solvent = base.normalise_condition(c_meta.get("solvent"))
        h_solvent = base.normalise_condition(h_meta.get("solvent"))
        matched = None if not c_solvent or not h_solvent else c_solvent == h_solvent
        return (c_field, c_meta), (h_field, h_meta), matched
    c_selected = carbon[0] if use_cnmr and carbon else (None, None)
    h_selected = proton[0] if use_hnmr and proton else (None, None)
    return c_selected, h_selected, None


def _clean_branch(record, branch_name):
    spectrum = record["spectra"].get(branch_name)
    supervision = record["alignment_supervision"][branch_name]
    if not spectrum:
        supervision["assigned_atom_indices_by_peak"] = []
        supervision["assignment_record_count_by_peak"] = []
        return

    clean_peaks = []
    assigned = []
    record_counts = []
    source_fragments = []
    for peak in spectrum["peaks"]:
        assigned.append(list(peak.get("assigned_atom_indices", [])))
        record_counts.append(int(peak.get("assignment_record_count", 0)))
        source_fragments.append(peak.get('source_assignment_fragments', []))
        clean = copy.deepcopy(peak)
        clean.pop("assigned_atom_indices", None)
        clean.pop("assignment_record_count", None)
        clean.pop('source_assignment_fragments', None)
        clean.pop('source_assignment_validation', None)
        clean_peaks.append(clean)
    spectrum["peaks"] = clean_peaks
    supervision["assigned_atom_indices_by_peak"] = assigned
    supervision["assignment_record_count_by_peak"] = record_counts
    supervision['source_assignment_fragments_by_peak'] = source_fragments


def molecule_record(mol, args):
    record, reason = _base_molecule_record(mol, args)
    if record is None:
        return None, reason
    _clean_branch(record, "carbon")
    _clean_branch(record, "proton")
    for field in ("c_nmr_peaks", "h_nmr_peaks"):
        for peak in record.get(field) or []:
            peak.pop("atom_indices", None)
    record["schema_version"] = "nmrpeak_atompeak_v2.1_leakage_safe"
    record["quality"]["notes"].append(
        "assigned atom indices are stored only in alignment_supervision and never in spectrum inputs"
    )
    return record, reason


base.choose_spectra = choose_spectra
base.molecule_record = molecule_record


if __name__ == "__main__":
    base.main()
