"""Portable exports from one consistent SQLite read snapshot."""
from __future__ import annotations

import hashlib
import json
import pickle
import sqlite3
import uuid
from pathlib import Path

from . import core as db
from .export_policy import Policy


def training_view(record):
    value = dict(record)
    current_spectra = [s for s in record['spectra'] if s.get('manual_state', 'active') == 'active'
                       and s.get('quality', {}).get('decision') != 'quarantined']
    value['spectra'] = current_spectra
    spectra = {s["id"]: s for s in current_spectra}
    pairs = sorted([p for p in record['pairs'] if p['c_id'] in spectra and p['h_id'] in spectra], key=lambda p: (
        -sum(bool(spectra[sid].get('version_info')) for sid in (p['c_id'], p['h_id'])),
        p["solvent_match"] != "same_known", "cross_database" in p["basis"], p["id"]))
    pair = pairs[0] if pairs else None
    value['pairs'] = pairs
    chosen = {}
    for nucleus in ("13C", "1H"):
        sid = pair["c_id" if nucleus == "13C" else "h_id"] if pair else None
        candidates = sorted([s for s in current_spectra if s['nucleus'] == nucleus],
                            key=lambda s: (not bool(s.get('version_info')), -(s.get('version_info') or {}).get('number', 0), s['id']))
        chosen[nucleus] = spectra.get(sid) or next(iter(candidates), None)
    value.update(smiles_canonical=record["smiles"], molecular_formula=record["formula"], selected_pair_id=pair["id"] if pair else None)
    value["selected_spectrum_ids"] = {nuc: s["id"] if s else None for nuc,s in chosen.items()}
    value["paired_solvent_known_match"] = bool(pair and pair["solvent_match"] == "same_known")
    value["pairing_basis"] = pair["basis"] if pair else "unpaired_independent_spectra"
    value["conditions"] = {k: v["conditions"] if v else None for k, v in chosen.items()}
    value["alignment_supervision"] = {}
    for nuc, name, branch in (("13C", "c_nmr_peaks", "carbon"), ("1H", "h_nmr_peaks", "proton")):
        peaks = chosen[nuc]["peaks"] if chosen[nuc] else []
        value[name] = []
        for p in peaks:
            # Model-facing peak features must not include source atom IDs,
            # raw third-field assignments, or the annotation answer itself.
            q = {k:p[k] for k in ('peak_id','shift_ppm','multiplicity','reported_nH','j_values',
                                 'line_width_hz','integral','range_min','range_max','source_signal_code') if k in p}
            # An accepted atom map does not repair contradictory S/D metadata.
            # Raw codes remain in record/audit, not as a trustworthy model feature.
            code_ok = not any(i.get('code')=='carbon_hcount_code_conflict'
                for f in p.get('source_assignment_fragments',[]) for i in f.get('issues',[]))
            if p.get('signal_code_revision'):
                code_ok = p.get('effective_signal_code_eligible', False)
                q['source_signal_code'] = p['effective_signal_code']
                q['signal_code_origin'] = 'manual_revision'
                if code_ok:q['carbon_hydrogen_class'] = p['signal_code_revision']['code']
            q['source_signal_code_eligible'] = code_ok
            if not code_ok:q.pop('source_signal_code',None)
            if nuc == "13C":
                q["delta (ppm)"] = p["shift_ppm"]
            else:
                q.update(delta=p["shift_ppm"], centroid=p["shift_ppm"], category=p.get("multiplicity"), nH=p.get("reported_nH"),
                         rangeMin=p.get("range_min", p["shift_ppm"]), rangeMax=p.get("range_max", p["shift_ppm"]))
            value[name].append(q)
        value["alignment_supervision"][branch] = {
            "assigned_atom_indices_by_peak": [p["assignment"]["atom_indices"] for p in peaks],
            "assignment_valid_by_peak": [p["assignment"]["status"] in db.ASSIGNED - {"source_group"} and p['assignment'].get('source_identity_eligible', True) for p in peaks],
            "assignment_present_by_peak": [p["assignment"]["status"] in db.ASSIGNED and p['assignment'].get('source_identity_eligible', True) for p in peaks],
            "assignment_status_by_peak": [p["assignment"]["status"] for p in peaks],
            "assignment_origin_by_peak": [p['assignment'].get('origin') for p in peaks],
            "ai_provenance_by_peak": [p['assignment'].get('ai') for p in peaks],
            "source_assignment_fragments_by_peak": [p.get('source_assignment_fragments', []) for p in peaks],
            "assignment_groups": [{**g, 'peak_indices': [i for i,p in enumerate(peaks) if p['peak_id'] in g['peak_ids']],
                                   'supervision_scope': 'set_only_no_within_group_order'}
                                  for g in (chosen[nuc].get('assignment_groups', []) if chosen[nuc] else [])],
            "target_resolution": "carbon_atom" if nuc == "13C" else "attached_heavy_atom"}
    # Canonical SMILES atom ordering is explicit; no fictional 3-D conformer is generated.
    value["export_note"] = "Portable annotation LMDB; use per-peak masks. No coordinates or automatic train/valid/test split."
    return value


def parquet_schema():
    import pyarrow as pa
    peak = pa.struct([("peak_id", pa.string()), ("shift_ppm", pa.float64()), ("multiplicity", pa.string()),
                      ("reported_nH", pa.float64()), ("j_values", pa.list_(pa.float64())),
                      ("status", pa.string()), ("atom_indices", pa.list_(pa.int32())),
                      ("annotation_origin", pa.string()), ("ai_provenance_json", pa.string())])
    spectrum = pa.struct([("id", pa.string()), ("nucleus", pa.string()), ("source", pa.string()),
                          ("solvent", pa.string()), ("frequency_mhz", pa.float64()), ("temperature_k", pa.float64()),
                          ("peaks", pa.list_(peak))])
    return pa.schema([("id", pa.string()), ("smiles", pa.string()), ("formula", pa.string()), ("origin", pa.string()),
                      ("queue", pa.string()), ("revision", pa.int64()), ("spectra", pa.list_(spectrum)),
                      ("record_json", pa.string())])


def parquet_row(r):
    v = {k: r[k] for k in ("id", "smiles", "formula", "origin", "queue", "revision")}
    v["record_json"] = db.dumps(r)
    v["spectra"] = []
    for s in r["spectra"]:
        item = {k: s[k] for k in ("id", "nucleus", "source", "solvent", "frequency_mhz", "temperature_k")}
        item["peaks"] = [{**{k: p.get(k) for k in ("peak_id", "shift_ppm", "multiplicity", "reported_nH", "j_values")},
                          "status": p["assignment"]["status"], "atom_indices": p["assignment"]["atom_indices"],
                          "annotation_origin": p['assignment'].get('origin'),
                          "ai_provenance_json": db.dumps(p['assignment'].get('ai'))} for p in s["peaks"]]
        v["spectra"].append(item)
    return v


def export_database(path, output_root, fmt, filters=None, molecule_id=None, training_options=None, progress=None, export_options=None):
    policy=Policy(export_options)
    if fmt not in ("jsonl", "parquet", "lmdb", "atompeak"):
        raise ValueError("Unsupported format")
    folder = Path(output_root) / (db.now().replace(":", "").replace("+", "_") + "_" + uuid.uuid4().hex[:8])
    if fmt == 'atompeak':
        from atompeak_data.export import export_training_dataset
        return export_training_dataset(path, folder, filters=filters, molecule_id=molecule_id,
            archive=True, progress=progress, export_options=export_options, **(training_options or {}))
    folder.mkdir(parents=True)
    target = folder / ("records." + fmt)
    c = sqlite3.connect(f"file:{Path(path).as_posix()}?mode=ro", uri=True, timeout=60)
    c.row_factory = sqlite3.Row
    c.execute("BEGIN")
    where, args = db.selection_where(filters or {})
    if molecule_id:
        where, args = " WHERE id=?", [molecule_id]
    ids = c.execute("SELECT id FROM molecules" + where + " ORDER BY id", args)
    total=c.execute('SELECT count(*) FROM molecules'+where,args).fetchone()[0]
    count = 0
    enforcement = c.execute("SELECT value FROM meta WHERE key='quality_export_enforced'").fetchone()
    strict = bool(enforcement and json.loads(enforcement[0]))
    skipped = 0
    def eligible_records():
        nonlocal skipped
        for row in ids:
            value = db.record(c, row[0], require_quality=strict)
            accepted=policy.accept(value)
            if progress and policy.scanned%100==0:progress({'processed':policy.scanned,'total':total,'exported':policy.kept})
            if not accepted:
                skipped += 1
                continue
            yield value
    try:
        if fmt == "jsonl":
            with target.open("w", encoding="utf-8", newline="\n") as f:
                for value in eligible_records():
                    f.write(db.dumps(value) + "\n")
                    count += 1
        elif fmt == "parquet":
            import pyarrow as pa
            import pyarrow.parquet as pq
            schema = parquet_schema()
            with pq.ParquetWriter(target, schema, compression="zstd") as writer:
                batch = []
                for value in eligible_records():
                    batch.append(parquet_row(value))
                    count += 1
                    if len(batch) >= 250:
                        writer.write_table(pa.Table.from_pylist(batch, schema=schema))
                        batch.clear()
                if batch:
                    writer.write_table(pa.Table.from_pylist(batch, schema=schema))
        else:
            import lmdb
            env = lmdb.open(str(target), subdir=False, map_size=32 * 1024 * 1024)
            try:
                batch = []
                def write_batch(batch):
                    while True:
                        try:
                            with env.begin(write=True) as txn:
                                for key, value in batch:
                                    txn.put(key, value)
                            return
                        except lmdb.MapFullError:
                            env.set_mapsize(env.info()["map_size"] * 2)
                for value in eligible_records():
                    batch.append((str(count).encode(), pickle.dumps(training_view(value), protocol=4)))
                    count += 1
                    if count % 250 == 0:
                        write_batch(batch)
                        batch.clear()
                if batch:
                    write_batch(batch)
                env.sync()
            finally:
                env.close()
        digest = hashlib.sha256()
        with target.open("rb") as f:
            for chunk in iter(lambda: f.read(1024*1024), b""):
                digest.update(chunk)
        manifest = {"schema_version": db.SCHEMA, "format": fmt, "records": count, "created_at": db.now(),
                    "source_assignment_identity_policy": "shiftdb-visible-atom-identity-v1",
                    "source_assignment_flags_are_not_chemical_correctness_certificates": True,
                    "quality_export_enforced": strict, "skipped_no_eligible_spectra": policy.reasons.get('no_eligible_spectra',0),
                    'skipped_by_export_selection':skipped,
                    "quarantined_spectra_excluded": True,
                    "manual_holds_exclusions_references_superseded_excluded": True,
                    "filters": filters or {}, "molecule_id": molecule_id, "sha256": digest.hexdigest(),
                    "export_selection":policy.report(),
                    "atom_index_basis": "zero_based_canonical_smiles", "h_assignment": "attached_heavy_atom",
                    "review_required": False, "has_raw_fid": False, "experimental_only": True,
                    "application_compatibility": "annotation-tool-v3; not a drop-in replacement for every NMRPeak training reader",
                    "scope": "curated molecules only; unresolved source candidates are not fabricated",
                    "spectrum_version_policy": "current revised spectra only; originals and historical versions are reference-only",
                    "split_policy": "not split; group identical connectivity when generating training splits"}
        (folder / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
        return {"file": str(target), "manifest": str(folder / "manifest.json"), "count": count,'selection':policy.report()}
    finally:
        c.close()


def preview_database(path,filters=None,molecule_id=None,export_options=None,progress=None,require_quality=False):
    policy=Policy(export_options)
    c=sqlite3.connect(Path(path).resolve().as_uri()+'?mode=ro',uri=True,timeout=60);c.row_factory=sqlite3.Row
    try:
        c.execute('BEGIN');where,args=db.selection_where(filters or {})
        if molecule_id:where,args=' WHERE id=?',[molecule_id]
        strict=c.execute("SELECT value FROM meta WHERE key='quality_export_enforced'").fetchone()
        strict=require_quality or bool(strict and json.loads(strict[0]));total=c.execute('SELECT count(*) FROM molecules'+where,args).fetchone()[0]
        examples=[]
        for row in c.execute('SELECT id FROM molecules'+where+' ORDER BY id',args):
            record=db.record(c,row[0],require_quality=strict)
            if policy.accept(record) and len(examples)<20:examples.append({'id':record['id'],'formula':record['formula']})
            if progress and policy.scanned%100==0:progress({'processed':policy.scanned,'total':total,'exported':policy.kept})
        return {'preview':True,'count':policy.kept,'selection':policy.report(),'examples':examples,
                'note':'预览和正式导出各自读取一致性快照；期间数据库修改可能使数量变化。训练包还会检查序列长度及可用监督，最终数量以导出清单为准。'}
    finally:c.close()
