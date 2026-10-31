"""Resumable, read-only-source migration. Run from the project root with -m."""
from __future__ import annotations

import argparse
import ast
import base64
import json
import math
import pickle
import re
import sqlite3
import sys
import time
import warnings
from collections import defaultdict
from pathlib import Path

import lmdb
import numpy as np
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger
from scipy.optimize import linear_sum_assignment

from . import core as db
from tools.filter_nmrshiftdb2_experimental import iter_records, parse_record, classify, indexed, SPECTRUM, sha256
from tools.shiftdb_atom_identity import canonical_identity, identity_for_current, source_peaks as shiftdb_source_peaks, parse_fragments, validate_fragments
from tools.shiftdb_numbering import assess_numbering

SOURCE_SDF = db.ROOT / "data/nmrshiftdb2_experimental_v1/experimental_ch.sd"
LMDB_DIR = db.ROOT / "data/NMRexp/lmdb_dataset"
PARQUET_DIR = db.ROOT / "data/NMRexp/source_parquet/ml_ready"
INDEX = db.ROOT / "data/atompeak_strict/pipeline.sqlite"
warnings.filterwarnings("ignore", category=DeprecationWarning)
RDLogger.DisableLog("rdApp.warning")


def log(value):
    print(f"[{db.now()}] {value}", flush=True)


def clean(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(x) for x in value]
    return value


def number(value):
    if value is None:
        return None
    match = re.search(r"[-+]?\d+(?:\.\d+)?", str(value))
    result = float(match[0]) if match else None
    return result if result is not None and math.isfinite(result) else None


def solvent(raw):
    if raw is None:
        return None
    compact = re.sub(r"[\s_()\-]", "", str(raw)).lower()
    names = {"CDCl3": ["cdcl3", "chloroformd1", "chloroformd"],
             "DMSO-d6": ["dmsod6", "dimethylsulfoxided6", "dimethylsulphoxided6"],
             "CD2Cl2": ["cd2cl2", "dichloromethaned2"], "CD3OD": ["cd3od", "methanold4"],
             "acetone-d6": ["acetoned6", "cd3cocd3"], "C6D6": ["c6d6", "benzened6"],
             "CD3CN": ["cd3cn", "acetonitriled3"], "D2O": ["d2o", "deuteriumoxide"],
             "THF-d8": ["thfd8", "tetrahydrofurand8"], "pyridine-d5": ["pyridined5"]}
    found = [name for name, terms in names.items() if any(term in compact for term in terms)]
    if len(found) == 1 and not any(x in compact for x in ("mixed", "+", "/", ":")):
        return found[0]
    if compact in {"", "none", "unknown", "na", "n/a", "notknown", "nan", "unreported", "notspecified"}:
        return None
    return str(raw).strip()


def shiftdb_peaks(value, mol, mapping, nucleus, *, identity=None, record_conflict=False):
    identity = identity or identity_for_current(mol, mapping)
    peaks, labels, validation = shiftdb_source_peaks(value, nucleus, identity, record_conflict=record_conflict)
    if validation['issues']:
        warnings.warn(f"ShiftDB2 {nucleus} {validation['status']}: {validation['issue_counts']}", RuntimeWarning)
    return peaks, labels


def normalize_shiftdb(c):
    if db.cursor(c, "shiftdb_signal_groups_v2") == "done":
        return
    changed = set()
    count = 0
    rows = c.execute("SELECT s.*,m.smiles FROM spectra s JOIN molecules m ON m.id=s.molecule_id WHERE s.source='nmrshiftdb2'").fetchall()
    for row in rows:
        if c.execute("SELECT 1 FROM annotations WHERE spectrum_id=? AND origin!='nmrshiftdb2' LIMIT 1", (row["id"],)).fetchone():
            db.issue(c,"normalize",row["id"],"manual_annotation_preserved_not_renumbered")
            continue
        source = json.loads(row["source_json"])
        value = source["original_properties"][source["spectrum_field"]]
        if not source.get('original_molblock'):
            # Old integer-only maps cannot tell removed explicit H from its
            # parent. Leave old data intact until the raw-SDF audit repairs it.
            db.issue(c, 'normalize', row['id'], 'legacy_source_identity_requires_raw_sdf_audit')
            continue
        original = Chem.MolFromMolBlock(source['original_molblock'], sanitize=False, removeHs=False)
        identity = canonical_identity(original, source.get('atom_identity', {}).get('source_key', row['source_key']), sanitize=True)
        if identity.smiles != row['smiles']:
            db.issue(c, 'normalize', row['id'], 'source_structure_changed_not_renumbered')
            continue
        conflict = any(validate_fragments(parse_fragments(v, '13C', identity), '13C')['status']=='source_assignment_inconsistency'
                       for k,v in source['original_properties'].items() if (m:=SPECTRUM.fullmatch(k)) and m[1]=='13C'
                       and classify(source['original_properties'],m[2],m[1])[0]=='experimental')
        peaks, labels = shiftdb_peaks(value, identity.processed, identity.legacy_mapping, row['nucleus'],
                                    identity=identity, record_conflict=conflict and row['nucleus']=='1H')
        conditions = json.loads(row["conditions_json"])
        conditions["solvent"] = solvent(conditions.get("solvent_raw"))
        conditions["frequency_scope"] = "source_reported_unspecified"
        c.execute("UPDATE spectra SET peaks_json=?,conditions_json=?,solvent=? WHERE id=?",(db.dumps(peaks),db.dumps(conditions),conditions["solvent"],row["id"]))
        c.execute("DELETE FROM annotations WHERE spectrum_id=? AND origin='nmrshiftdb2'",(row["id"],))
        for pid,a in labels.items():
            c.execute("INSERT INTO annotations VALUES(?,?,?,?,?,?,?)",(row["id"],pid,a["status"],db.dumps(a["atom_indices"]),a.get("note",""),"nmrshiftdb2",db.now()))
        changed.add(row["molecule_id"])
        count+=1
        if count%2000==0:
            c.commit()
            log(f"Normalize ShiftDB source signal groups {count:,}")
    for mid in changed:
        db.refresh_molecule(c,mid)
    for pair in c.execute("SELECT * FROM pairs").fetchall():
        db.add_pair(c,pair["molecule_id"],pair["c_id"],pair["h_id"],pair["basis"])
    c.execute("UPDATE candidates SET state='duplicate_complete' WHERE state='awaiting_source' AND molecule_id IN (SELECT id FROM molecules WHERE queue='complete')")
    db.progress(c,"shiftdb_signal_groups_v2","done")


def migrate_shiftdb(c, path=SOURCE_SDF):
    stage = "shiftdb"
    if db.cursor(c, stage) == "done":
        return
    start = int(db.cursor(c, stage))
    digest = sha256(path)
    old = c.execute("SELECT value FROM meta WHERE key='shiftdb_sha256'").fetchone()
    if old and json.loads(old[0]) != digest:
        raise ValueError("Source SDF changed; use a new database")
    db.set_meta(c, "shiftdb_sha256", digest)
    db.set_meta(c, "shiftdb_path", str(path))
    for i, _, raw in iter_records(path):
        if i <= start:
            continue
        try:
            block, props, _ = parse_record(raw)
            original = Chem.MolFromMolBlock(block.decode(), removeHs=False, sanitize=False)
            if original is None:
                raise ValueError("invalid_structure")
            identity = canonical_identity(original, f"shiftdb:{props.get('Original SDF record index', i)}", sanitize=True)
            smi, mol, mapping = identity.smiles, identity.processed, identity.legacy_mapping
            # A malformed experimental carbon identity can taint otherwise
            # plausible proton indices; don't silently trust that record's H.
            record_conflict = any(
                validate_fragments(parse_fragments(v, '13C', identity), '13C')['status'] == 'source_assignment_inconsistency'
                for k,v in props.items() if (m:=SPECTRUM.fullmatch(k)) and m[1]=='13C'
                and classify(props,m[2],m[1])[0]=='experimental')
            items = []
            for key, value in props.items():
                match = SPECTRUM.fullmatch(key)
                if not match or match[1] not in ("13C", "1H"):
                    continue
                if classify(props, match[2], match[1])[0] != "experimental":
                    db.issue(c, stage, f"{i}:{key}", "not_experimental")
                    continue
                peaks, labels = shiftdb_peaks(value, mol, mapping, match[1], identity=identity,
                    record_conflict=record_conflict and match[1]=='1H')
                if not peaks:
                    continue
                conditions = {}
                for source_name, out in [("Solvent", "solvent_raw"), ("Temperature [K]", "temperature_raw"), ("Field Strength [MHz]", "frequency_raw")]:
                    conditions[out] = indexed(props.get(source_name, "")).get(match[2])
                conditions.update(solvent=solvent(conditions["solvent_raw"]), temperature_k=number(conditions["temperature_raw"]), frequency_mhz=number(conditions["frequency_raw"]))
                validation=shiftdb_source_peaks(value,match[1],identity,
                    record_conflict=record_conflict and match[1]=='1H')[2]
                if validation['status']!='identity_validated':
                    validation['numbering_workflow']=assess_numbering(value,match[1],identity,
                        source_inchi=props.get('INChI'),record_conflict=record_conflict and match[1]=='1H')
                items.append((key, match[1], conditions, peaks, labels, validation))
            if items:
                mid = db.add_molecule(c, smi, "nmrshiftdb2")
                imported = []
                for key, nucleus, conditions, peaks, labels, validation in items:
                    sid = db.add_spectrum(c, mid, nucleus, "nmrshiftdb2", f"shiftdb:{i}:{key}", conditions, peaks,
                        {"source_file": str(path.relative_to(db.ROOT)), "source_sha256": digest,
                         "source_record_index": props.get("Original SDF record index", str(i)), "nmrshiftdb2_id": props.get("nmrshiftdb2 ID"),
                         "spectrum_field": key, "experimental_filter": "nmrshiftdb2-experimental-v1",
                         "original_atom_mapping": mapping, "original_properties": props,
                         "original_molblock": block.decode(), "atom_identity": identity.to_dict(),
                         "source_assignment_validation": validation}, labels)
                    imported.append((sid, nucleus, conditions))
                for cs, cn, cc in imported:
                    for hs, hn, hc in imported:
                        if cn == "13C" and hn == "1H" and cc["solvent"] == hc["solvent"] and cc["solvent"]:
                            db.add_pair(c, mid, cs, hs, "same_shiftdb_record_sample_unknown")
                db.refresh_molecule(c, mid)
        except sqlite3.Error:
            raise
        except Exception as exc:
            db.issue(c, stage, str(i), str(exc))
        if i % 500 == 0:
            db.progress(c, stage, i)
            log(f"ShiftDB2 {i:,} records")
    if sha256(path) != digest:
        raise ValueError("Input changed during migration")
    c.execute("UPDATE candidates SET state='duplicate_complete' WHERE state='awaiting_source' AND molecule_id IN (SELECT id FROM molecules WHERE queue='complete')")
    db.progress(c, stage, "done")


def build_candidates(c, index=INDEX):
    if db.cursor(c, "candidate_index") == "done":
        return
    source = sqlite3.connect(f"file:{index.as_posix()}?mode=ro", uri=True)
    source.row_factory = sqlite3.Row
    # Old keep_record exclusions were built from mixed ShiftDB spectra: do not reuse them.
    rows = source.execute("SELECT canonical_smiles,exact_structure_key,source_split,source_locator FROM source_records WHERE source='nmrexp' AND has_carbon=1 AND has_proton=1 AND canonical_smiles IS NOT NULL ORDER BY id")
    count = 0
    for r in rows:
        key = base64.b64decode(r["source_locator"]).decode()
        c.execute("INSERT OR IGNORE INTO candidates(smiles,exact_key,lmdb_split,lmdb_key) VALUES(?,?,?,?)", (r["canonical_smiles"], r["exact_structure_key"], r["source_split"], key))
        count += 1
        if count % 50000 == 0:
            c.commit()
            log(f"NMRPeak index {count:,}")
    source.close()
    db.progress(c, "candidate_index", "done")


def shifts(record, nucleus):
    result = []
    for p in record.get("c_nmr_peaks" if nucleus == "13C" else "h_nmr_peaks", []) or []:
        raw = next((p[k] for k in (("delta (ppm)", "shift_ppm", "shift") if nucleus == "13C" else ("centroid", "delta", "shift_ppm")) if p.get(k) is not None), None) if isinstance(p, dict) else p
        try:
            val = float(raw)
            if math.isfinite(val):
                result.append(val)
        except (ValueError, TypeError):
            pass
    return result


def hydrate_candidates(c):
    if db.cursor(c, "candidate_details") == "done":
        return
    envs = {name: lmdb.open(str(LMDB_DIR / f"{name}.lmdb"), subdir=False, readonly=True, lock=False, readahead=False) for name in ("train", "valid", "test")}
    transactions = {name: env.begin() for name, env in envs.items()}
    count = 0
    while True:
        rows = c.execute("SELECT * FROM candidates WHERE state='indexed' ORDER BY lmdb_split,lmdb_key LIMIT 1000").fetchall()
        if not rows:
            break
        for row in rows:
            try:
                raw = transactions[row["lmdb_split"]].get(row["lmdb_key"].encode())
                if raw is None:
                    raise ValueError("LMDB key missing")
                value = pickle.loads(raw)  # Only trusted, user-provided local dataset files.
                cs, hs = shifts(value, "13C"), shifts(value, "1H")
                if not cs or not hs:
                    raise ValueError("Indexed entry no longer contains both spectra")
                old = c.execute("SELECT id,queue FROM molecules WHERE smiles=?", (row["smiles"],)).fetchone()
                state = "duplicate_complete" if old and old["queue"] == "complete" else "awaiting_source"
                c.execute("UPDATE candidates SET raw_smiles=?,solvent=?,c_shifts=?,h_shifts=?,state=?,molecule_id=? WHERE smiles=?",
                    (value.get("smiles") or value.get("smiles_canonical"), solvent(value.get("nmr_solvent")), db.dumps(cs), db.dumps(hs), state, old["id"] if old else None, row["smiles"]))
            except Exception as exc:
                c.execute("UPDATE candidates SET state='unavailable',reason=? WHERE smiles=?", (str(exc), row["smiles"]))
        count += len(rows)
        c.commit()
        if count % 10000 == 0:
            log(f"NMRPeak details {count:,}")
    for env in envs.values():
        env.close()
    db.progress(c, "candidate_details", "done")


def index_sources(c):
    # Build source mapping for all indexed structures, without pretending LMDB index == Parquet row.
    wanted = {x[0] for x in c.execute("SELECT raw_smiles FROM candidates WHERE state='awaiting_source'")}
    wanted.update(x[0] for x in c.execute("SELECT smiles FROM candidates WHERE state='awaiting_source'"))
    for path in sorted(PARQUET_DIR.glob("*.parquet")):
        stage = "parquet:" + path.name
        if db.cursor(c, stage) == "done":
            continue
        start = int(db.cursor(c, stage))
        offset = 0
        for batch in pq.ParquetFile(path).iter_batches(batch_size=8192):
            if offset + batch.num_rows <= start:
                offset += batch.num_rows
                continue
            values = batch.to_pylist()
            inserts = []
            for j, row in enumerate(values):
                if offset + j < start or row.get("SMILES") not in wanted:
                    continue
                nuc = str(row.get("NMR_type") or "").split()[0]
                if nuc not in ("13C", "1H"):
                    continue
                row = clean(row)
                inserts.append((f"{path.name}:{offset+j}", row["SMILES"], nuc, row.get("Filename"), solvent(row.get("NMR_solvent")), row.get("Page_in_file_para"), path.name, offset+j, db.dumps(row)))
            c.executemany("INSERT OR IGNORE INTO source_rows VALUES(?,?,?,?,?,?,?,?,?)", inserts)
            offset += batch.num_rows
            db.progress(c, stage, offset)
            if offset % 131072 == 0:
                log(f"Source {path.name}: {offset:,}")
        db.set_meta(c, "parquet_sha256:" + path.name, sha256(path))
        db.progress(c, stage, "done")


def source_peaks(row, nucleus):
    result = []
    values = ast.literal_eval(row["NMR_processed"])
    for i, p in enumerate(values):
        pid = ("C" if nucleus == "13C" else "H") + str(i+1).zfill(3)
        if nucleus == "1H":
            category, js, nh, high, low = p[:5]
            high, low = max(float(high), float(low)), min(float(high), float(low))
            peak = {"peak_id": pid, "shift_ppm": (high+low)/2, "range_min": low, "range_max": high,
                    "multiplicity": category, "reported_nH": number(nh), "raw_value": p}
        else:
            val, category, js = p[:3]
            try:
                shift = float(val)
                peak = {"peak_id": pid, "shift_ppm": shift}
            except (ValueError, TypeError):
                matches = re.findall(r"-?\d+(?:\.\d+)?", str(val).replace("–", " ").replace("-", " "))
                if len(matches) != 2:
                    raise ValueError("unparsed carbon range")
                low, high = sorted(map(float, matches))
                peak = {"peak_id": pid, "shift_ppm": (low+high)/2, "range_min": low, "range_max": high}
            peak.update(multiplicity=category, reported_nH=None, raw_value=p)
        if isinstance(js, str):
            try:
                js = ast.literal_eval(js)
            except (ValueError, SyntaxError):
                js = [js]
        peak["j_values"] = [number(x) for x in (js or []) if number(x) is not None]
        peak["j_unit"] = "Hz" if peak["j_values"] and all("hz" in str(x).lower() for x in js) else None
        if peak["j_values"] and not peak["j_unit"]:
            peak["j_values_unconfirmed"] = peak["j_values"]
            peak["j_values"] = None
        if not math.isfinite(peak["shift_ppm"]):
            raise ValueError("nonfinite shift")
        result.append(peak)
    return result


def fingerprint(local, remote, tolerance):
    from .quality import peak_match
    match = peak_match(local, remote, tolerance)
    return match["f1"], match["median"] if match["median"] is not None else 1e9


def recover_one(c, candidate, rows=None, rematch=False):
    from . import quality
    if rows is None:
        rows = c.execute("SELECT * FROM source_rows WHERE raw_smiles IN (?,?)", (candidate["raw_smiles"], candidate["smiles"])).fetchall()
    mid = candidate["molecule_id"] if "molecule_id" in candidate.keys() else None
    if not mid:
        existing = c.execute("SELECT id FROM molecules WHERE smiles=?", (candidate["smiles"],)).fetchone()
        mid = existing[0] if existing else None
    refs = quality.references(c, mid)
    groups = defaultdict(lambda: {"13C": [], "1H": []})
    locals_ = {"13C": json.loads(candidate["c_shifts"]), "1H": json.loads(candidate["h_shifts"])}
    for r in rows:
        if not r["solvent"] or r["solvent"].lower() == "mixed":
            continue
        if not rematch and candidate["solvent"] and r["solvent"] != candidate["solvent"]:
            continue
        if rematch and refs and not any(x["solvent"] == r["solvent"] for x in refs):
            continue
        raw = json.loads(r["raw_json"])
        try:
            peaks = source_peaks(raw, r["nucleus"])
            f1, err = fingerprint(locals_[r["nucleus"]], [p["shift_ppm"] for p in peaks], .035 if r["nucleus"] == "1H" else .30)
        except Exception:
            continue
        if rematch or f1 >= .95:
            groups[(r["doi_encoded"], r["solvent"], r["raw_smiles"])][r["nucleus"]].append((r, raw, peaks, f1, err))
    pairs = []
    failed_checks = defaultdict(int)
    for (doi, sol, _), g in groups.items():
        for h in g["1H"]:
            for cs in g["13C"]:
                if not doi or h[0]["page"] is None or cs[0]["page"] is None or abs(h[0]["page"]-cs[0]["page"]) > 3:
                    continue
                hf, cf = number(h[1].get("NMR_frequency")), number(cs[1].get("NMR_frequency"))
                if hf and cf and not 3.75 <= hf/cf <= 4.10:
                    continue
                score = h[3]+cs[3] - .1*h[4] - .01*cs[4] - .002*abs(h[0]["page"]-cs[0]["page"])
                try:
                    check = quality.evaluate(candidate["smiles"], h, cs, sol, refs)
                except (ValueError, TypeError, KeyError) as exc:
                    failed_checks["source_parse_error:" + str(exc)] += 1
                    continue
                if check["decision"] != "passed":
                    for reason in check["reasons"]:
                        failed_checks[reason] += 1
                    continue
                if rematch:
                    # A deterministic ranking score, NOT a calibrated probability.
                    shared = [max(m["f1"] - m["median"]/10
                                  for m in check["metrics"][nuc] if m["passed"])
                              for nuc in ("13C", "1H") if check["metrics"].get(nuc)]
                    if not shared:
                        continue
                    score = min(shared)
                    check["selection_basis"] = "shared_experimental_spectrum_not_lmdb_fingerprint"
                pairs.append((score, h, cs, sol, check))
    pairs.sort(key=lambda x: (-x[0], x[1][0]["id"], x[2][0]["id"]))
    if not pairs:
        return None, "no_source_pair_passing_quality_checks:" + db.dumps(dict(sorted(failed_checks.items())))
    best = pairs[0]
    def identity(p):
        return (p[1][0]["doi_encoded"], p[1][0]["page"], p[2][0]["page"], p[1][0]["raw_json"], p[2][0]["raw_json"])
    if any(best[0]-p[0] < (.02 if rematch else .005) and identity(best) != identity(p) for p in pairs[1:]):
        return None, "multiple_source_pairs"
    return best, None


def import_recovered(c, candidate, found):
    from . import quality
    _, h, cs, sol = found[:4]
    # Validate exact source structure at the boundary; mapping is never based on LMDB ordinal.
    smi, mol, _ = db.canonicalize(Chem.MolFromSmiles(h[1]["SMILES"]))
    if smi != candidate["smiles"]:
        raise ValueError("source_structure_mismatch")
    old = c.execute("SELECT * FROM molecules WHERE smiles=?", (smi,)).fetchone()
    check = quality.evaluate(smi, h, cs, sol, quality.references(c, old["id"] if old else None))
    if check["decision"] != "passed":
        raise ValueError("quality_gate:" + ",".join(check["reasons"]))
    if old and old["queue"] == "complete":
        return old["id"], "duplicate_complete"
    mid = old["id"] if old else db.add_molecule(c, smi, "NMRexp")
    added = {}
    for nucleus, item in (("13C", cs), ("1H", h)):
        # Existing experimental spectra are preserved, never copied as labels onto a new measurement.
        existing = c.execute("SELECT id,solvent,nucleus FROM spectra WHERE molecule_id=? AND nucleus=? AND" + db.ACTIVE, (mid, nucleus)).fetchall()
        if old and existing:
            continue
        if old:
            opposite = c.execute("SELECT id,solvent FROM spectra WHERE molecule_id=? AND nucleus!=? AND" + db.ACTIVE, (mid, nucleus)).fetchall()
            if not any(x["solvent"] == sol for x in opposite):
                continue
        r, raw, peaks, f1, err = item
        conditions = {"solvent": sol, "solvent_raw": raw.get("NMR_solvent"), "frequency_mhz": number(raw.get("NMR_frequency")),
                      "frequency_raw": raw.get("NMR_frequency"), "temperature_k": None, "ph": None}
        provenance = {"source_shard": r["shard"], "source_row_index": r["row_index"], "source_row_index_base": 0,
                      "doi": (r["doi_encoded"] or "").replace("_", "/", 1), "page_index_scope": "filtered_nmr_pages",
                      "original_row": raw, "lmdb_split": candidate["lmdb_split"], "lmdb_key": candidate["lmdb_key"],
                      "recovery_f1": f1, "recovery_median_error_ppm": err,
                      "retrieval_fingerprint_is_not_chemical_validation": True,
                      "selection_basis": found[4].get("selection_basis", "lmdb_retrieval_plus_quality_checks") if len(found) > 4 else "quality_checked_import"}
        added[nucleus] = db.add_spectrum(c, mid, nucleus, "NMRexp", "nmrexp:"+r["id"], conditions, peaks, provenance)
    if old and not added:
        return mid, "overlap_not_complementable"
    if not old:
        db.add_pair(c, mid, added["13C"], added["1H"], "same_doi_structure_solvent_sample_unknown")
    else:
        for nucleus, sid in added.items():
            for other in c.execute("SELECT id,nucleus FROM spectra WHERE molecule_id=? AND nucleus!=? AND solvent=? AND" + db.ACTIVE, (mid, nucleus, sol)):
                db.add_pair(c, mid, sid if nucleus == "13C" else other["id"], sid if nucleus == "1H" else other["id"], "cross_database_same_structure_solvent_sample_unknown")
    quality.write_result(c, mid, candidate, check, added.values())
    if old and added:
        c.execute("UPDATE molecules SET revision=revision+1,task_done=0 WHERE id=?", (mid,))
    db.refresh_molecule(c, mid)
    return mid, "complemented" if old else "imported"


def recovery_write(c, operation, attempts=8):
    """Retry an entire short transaction; never upgrade an old read snapshot."""
    for attempt in range(attempts):
        try:
            c.execute("BEGIN IMMEDIATE")
            result = operation()
            c.commit()
            return result
        except sqlite3.Error as exc:
            c.rollback()
            if not db.is_busy_error(exc) or attempt + 1 == attempts:
                raise
            time.sleep(min(.05 * 2 ** attempt, .5))
        except Exception:
            c.rollback()
            raise


def recover(c, limit=0):
    if c.in_transaction:
        raise RuntimeError("Recovery must start outside a write transaction")
    c.execute("PRAGMA cache_size=-262144")
    old_timeout = c.execute("PRAGMA busy_timeout").fetchone()[0]
    c.execute("PRAGMA busy_timeout=250")
    count = 0
    previous = db.cursor(c, "recover")
    base = int(previous) if previous.isdigit() else 0
    try:
        # Read / fingerprint without the writer lock. Each candidate is its own
        # atomic commit, leaving room for interactive saves between candidates.
        while not limit or count < limit:
            rows = c.execute("SELECT * FROM candidates WHERE state='awaiting_source' ORDER BY molecule_id IS NULL,smiles LIMIT ?", (min(100, limit-count) if limit else 100,)).fetchall()
            if not rows:
                break
            for cand in rows:
                error = None
                try:
                    found, reason = recover_one(c, cand)
                except sqlite3.Error:
                    raise
                except Exception as exc:
                    found, reason, error = None, str(exc), str(exc)

                def write_candidate():
                    current = c.execute("SELECT state FROM candidates WHERE smiles=?", (cand["smiles"],)).fetchone()
                    if not current or current[0] != "awaiting_source":
                        return
                    c.execute("SAVEPOINT candidate_import")
                    try:
                        if error:
                            raise ValueError(error)
                        if found is None:
                            state, mid, pair = "source_unresolved", cand["molecule_id"], None
                        else:
                            mid, state = import_recovered(c, cand, found)
                            pair = db.dumps({"h_source_id": found[1][0]["id"], "c_source_id": found[2][0]["id"]})
                        c.execute("UPDATE candidates SET state=?,reason=?,molecule_id=?,source_pair_json=? WHERE smiles=?", (state, reason, mid, pair, cand["smiles"]))
                    except sqlite3.Error:
                        raise
                    except Exception as exc:
                        c.execute("ROLLBACK TO candidate_import")
                        db.issue(c, "recover", cand["smiles"], str(exc))
                        c.execute("UPDATE candidates SET state='source_unresolved',reason=? WHERE smiles=?", (str(exc), cand["smiles"]))
                    finally:
                        c.execute("RELEASE candidate_import")

                recovery_write(c, write_candidate)
                count += 1
            recovery_write(c, lambda: c.execute("INSERT INTO progress VALUES('recover',?,?) ON CONFLICT(stage) DO UPDATE SET cursor=excluded.cursor,updated_at=excluded.updated_at", (str(base + count), db.now())))
            if count % 2000 == 0:
                log(f"Recovered candidates {base + count:,}")
        if not c.execute("SELECT 1 FROM candidates WHERE state='awaiting_source' LIMIT 1").fetchone():
            recovery_write(c, lambda: c.execute("INSERT INTO progress VALUES('recover','done',?) ON CONFLICT(stage) DO UPDATE SET cursor=excluded.cursor,updated_at=excluded.updated_at", (db.now(),)))
    finally:
        c.execute(f"PRAGMA busy_timeout={old_timeout}")


def migrate_legacy(c):
    folder = db.ROOT / "data/atom_assignment_annotation_v2"
    targets = {x["sample_id"]: x for x in (json.loads(line) for line in (folder / "new_targets_200.provenance_checked.jsonl").read_text(encoding="utf-8").splitlines() if line)}
    sessions = json.loads((folder / "annotation_sessions.json").read_text(encoding="utf-8"))
    for index, session in enumerate(sessions, 1):
        task_id = f"TASK-{index:04d}"
        file = folder / "annotations_v2" / f"{task_id}.json"
        legacy_id = session["sample_id"]
        if not file.exists() or legacy_id not in targets:
            continue
        previous = c.execute("SELECT state FROM legacy_annotations WHERE task_id=?", (task_id,)).fetchone()
        if previous and previous[0] not in {"waiting_target", "partially_imported_peak_ambiguity"}:
            continue
        payload = json.loads(file.read_text(encoding="utf-8"))
        # These files can contain generated suggestions, not human annotations.
        # Preserve their original bytes semantically, without promoting them to labels.
        if str(payload.get("reviewer", "")).lower() in {"gpt", "chatgpt", "ai"}:
            c.execute("INSERT INTO legacy_annotations VALUES(?,?,?,?)", (task_id,legacy_id,db.dumps(payload),"reference_only_not_human_labels"))
            continue
        actionable = [a for nuc in payload.get("assignments", {}).values() for a in nuc.values() if a.get("status") != "unreviewed"]
        if actionable and all(a.get("status") in db.ASSIGNED and not a.get("atom_indices") for a in actionable):
            c.execute("INSERT INTO legacy_annotations VALUES(?,?,?,?) ON CONFLICT(task_id) DO UPDATE SET payload_json=excluded.payload_json,state=excluded.state", (task_id,legacy_id,db.dumps(payload),"invalid_old_assignment_preserved_not_imported"))
            continue
        old = targets[legacy_id]
        smi, mol, mapping = db.canonicalize(Chem.MolFromSmiles(old["structure"]["smiles_canonical"]))
        target = c.execute("SELECT id FROM molecules WHERE smiles=?", (smi,)).fetchone()
        state = "waiting_target"
        if target:
            state = "imported"
            for nuc, rows in payload.get("assignments", {}).items():
                local = {p["peak_id"]: p for p in old.get("spectra", {}).get(nuc, {}).get("peaks", [])}
                for pid, a in rows.items():
                    if a.get("status") == "unreviewed" or pid not in local:
                        continue
                    origin = str(a.get("label_origin", "")).lower()
                    if any(term in origin for term in ("gpt", "pseudo", "model", "predict", "dft", "generated")):
                        state = "reference_labels_not_imported"
                        continue
                    matches = []
                    old_source = old.get("provenance", {}).get("c_source" if nuc == "13C" else "h_source", {})
                    for spectrum in c.execute("SELECT id,peaks_json,source_json FROM spectra WHERE molecule_id=? AND nucleus=? AND source='NMRexp'", (target[0], nuc)):
                        actual_source = json.loads(spectrum["source_json"])
                        if not old_source or old_source.get("source_shard") != actual_source.get("source_shard") or old_source.get("source_row") != actual_source.get("source_row_index"):
                            continue
                        for peak in json.loads(spectrum["peaks_json"]):
                            if abs(peak["shift_ppm"]-local[pid]["shift_ppm"]) < 1e-5:
                                matches.append((spectrum["id"], peak["peak_id"]))
                    if len(matches) != 1:
                        state = "partially_imported_peak_ambiguity"
                        continue
                    status = {"unresolved": "ambiguous", "unobserved": "excluded", "impurity": "impurity"}.get(a.get("status"), a.get("status"))
                    if status not in db.STATES:
                        continue
                    old_atoms = a.get("atom_indices", [])
                    if (status in db.ASSIGNED and not old_atoms) or any(i not in mapping for i in old_atoms):
                        state = "invalid_old_assignment_preserved_not_imported"
                        continue
                    sid, new_pid = matches[0]
                    atoms = [mapping[i] for i in old_atoms]
                    existing = c.execute("SELECT * FROM annotations WHERE spectrum_id=? AND peak_id=?", (sid,new_pid)).fetchone()
                    if existing:
                        if existing["status"] != status or json.loads(existing["atoms_json"]) != atoms:
                            state = "legacy_duplicate_conflict_preserved"
                        continue
                    c.execute("INSERT INTO annotations VALUES(?,?,?,?,?,?,?)", (sid,new_pid,status,db.dumps(atoms),a.get("note", ""),"manual_v2",db.now()))
            db.refresh_molecule(c, target[0])
        c.execute("INSERT INTO legacy_annotations VALUES(?,?,?,?) ON CONFLICT(task_id) DO UPDATE SET payload_json=excluded.payload_json,state=excluded.state", (task_id,legacy_id,db.dumps(payload),state))
    c.commit()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["all", "shiftdb", "candidates", "sources", "recover", "legacy", "report"])
    p.add_argument("--db", type=Path, default=db.DEFAULT_DB)
    p.add_argument("--limit", type=int, default=0, help="Recovery limit; 0 means all. Safe to resume.")
    args = p.parse_args()
    c = db.connect(args.db)
    db.initialize(c)
    if args.stage in ("all", "shiftdb"):
        migrate_shiftdb(c)
        normalize_shiftdb(c)
    if args.stage in ("all", "candidates"):
        build_candidates(c)
        hydrate_candidates(c)
    if args.stage in ("all", "sources"):
        index_sources(c)
    if args.stage in ("all", "recover"):
        recover(c, args.limit)
    if args.stage in ("all", "legacy"):
        migrate_legacy(c)
    report = db.stats(c)
    from .paths import artifacts_root
    (artifacts_root(args.db) / "migration_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    log(db.dumps(report))
    c.close()


if __name__ == "__main__":
    main()
