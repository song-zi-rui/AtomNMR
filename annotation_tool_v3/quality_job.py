"""Resumable QC, conservative re-supplementation and batched source recovery."""
from __future__ import annotations

import argparse
import json
import os
import threading
import time
from collections import Counter, defaultdict
from pathlib import Path

from . import core as db, migrate, quality
from .paths import artifacts_root


def source_item(row):
    raw = json.loads(row["raw_json"])
    return row, raw, migrate.source_peaks(raw, row["nucleus"]), None, None


def audit_one(c, candidate):
    pair = json.loads(candidate["source_pair_json"] or "{}")
    items = {nuc: c.execute("SELECT * FROM source_rows WHERE id=?", (pair.get(key),)).fetchone()
             for nuc, key in (("1H", "h_source_id"), ("13C", "c_source_id"))}
    try:
        if not all(items.values()):
            raise ValueError("missing_original_source_pair")
        h, carbon = source_item(items["1H"]), source_item(items["13C"])
        return quality.evaluate(candidate["smiles"], h, carbon, items["1H"]["solvent"], quality.references(c, candidate["molecule_id"]))
    except (ValueError, TypeError, KeyError, SyntaxError) as exc:
        return {"policy": quality.POLICY, "decision": "quarantined", "reasons": [str(exc)], "source_document_verified": False}


def apply_audit(c, candidate, result):
    # All original NMRexp spectra on this molecule stay in SQLite, including labels.
    pair = json.loads(candidate["source_pair_json"] or "{}")
    keys = ["nmrexp:" + str(pair.get(k)) for k in ("h_source_id", "c_source_id")]
    ids = [x[0] for x in c.execute("SELECT id FROM spectra WHERE molecule_id=? AND source='NMRexp' AND source_key IN (?,?)", (candidate["molecule_id"], *keys))]
    quality.write_result(c, candidate["molecule_id"], candidate, result, ids)
    if result["decision"] != "passed":
        c.execute("UPDATE candidates SET state='quality_held',reason=? WHERE smiles=?", (db.dumps(result["reasons"]), candidate["smiles"]))
    c.execute("INSERT INTO quality_history(molecule_id,smiles,old_pair_json,new_pair_json,details_json,created_at) VALUES(?,?,?,?,?,?)",
              (candidate["molecule_id"], candidate["smiles"], candidate["source_pair_json"], candidate["source_pair_json"], db.dumps(result), db.now()))


def commit_prepared(c, values, apply, max_items=16, max_seconds=.035):
    """One writer, small durable transactions, with per-item rollback boundaries."""
    pos = 0
    while pos < len(values):
        def transaction():
            start = time.perf_counter()
            consumed = 0
            for value in values[pos:pos + max_items]:
                apply(*value)
                consumed += 1
                if time.perf_counter() - start >= max_seconds:
                    break
            return consumed
        pos += migrate.recovery_write(c, transaction)


def audit(c, where, progress=None, limit=0):
    count = 0
    last_smiles = ""
    decisions = Counter()
    reasons = Counter()
    while not limit or count < limit:
        rows = c.execute("SELECT candidates.* FROM candidates LEFT JOIN candidate_quality q USING(smiles) WHERE " + where +
                         " AND candidates.smiles>? AND (q.policy IS NULL OR q.policy!=?) ORDER BY candidates.smiles LIMIT ?",
                         (last_smiles, quality.POLICY, min(128, limit-count) if limit else 128)).fetchall()
        if not rows:
            break
        prepared = [(r, audit_one(c, r)) for r in rows]
        commit_prepared(c, prepared, lambda r, check: apply_audit(c, r, check))
        for r, check in prepared:
            decisions[check["decision"]] += 1
            reasons.update(check.get("reasons", []))
        count += len(rows)
        last_smiles = rows[-1]["smiles"]
        if progress:
            progress(count, dict(decisions))
    return {"processed": count, "decisions": dict(decisions), "reasons": dict(reasons)}


def rematch(c, progress=None):
    # Include previously failed overlap attempts, not only the original 148.
    rows = c.execute("SELECT candidates.* FROM molecules JOIN candidates INDEXED BY candidate_molecule ON candidates.molecule_id=molecules.id "
                     "WHERE molecules.origin='nmrshiftdb2' AND molecules.queue='waiting_spectrum' ORDER BY candidates.smiles").fetchall()
    counts = Counter()
    examples = []
    for index, candidate in enumerate(rows, 1):
        found, reason = migrate.recover_one(c, candidate, rematch=True)
        if found:
            def write():
                mid, status = migrate.import_recovered(c, candidate, found)
                pair = db.dumps({"h_source_id": found[1][0]["id"], "c_source_id": found[2][0]["id"]})
                c.execute("INSERT INTO quality_history(molecule_id,smiles,old_pair_json,new_pair_json,details_json,created_at) VALUES(?,?,?,?,?,?)",
                          (mid, candidate["smiles"], candidate["source_pair_json"], pair, db.dumps(found[4]), db.now()))
                c.execute("UPDATE candidates SET state=?,reason=NULL,source_pair_json=? WHERE smiles=?", (status,pair,candidate["smiles"]))
                return status
            status = migrate.recovery_write(c, write)
            counts[status] += 1
            examples.append({"molecule_id": candidate["molecule_id"], "h_source_id": found[1][0]["id"], "c_source_id": found[2][0]["id"]})
        else:
            counts["no_unique_safe_replacement"] += 1
            migrate.recovery_write(c, lambda: c.execute("UPDATE candidates SET reason=? WHERE smiles=?", (reason, candidate["smiles"])))
        if progress:
            progress(index, dict(counts))
    return {"searched": len(rows), "results": dict(counts), "new_pairs": examples}


def recover_batched(c, progress=None, limit=0):
    count, decisions = 0, Counter()
    while not limit or count < limit:
        rows = c.execute("SELECT * FROM candidates WHERE state='awaiting_source' ORDER BY molecule_id IS NULL,smiles LIMIT ?", (min(128, limit-count) if limit else 128,)).fetchall()
        if not rows:
            break
        wanted = sorted({s for r in rows for s in (r["smiles"], r["raw_smiles"]) if s})
        grouped = defaultdict(list)
        for r in c.execute("SELECT * FROM source_rows WHERE raw_smiles IN (" + ",".join("?" for _ in wanted) + ")", wanted):
            grouped[r["raw_smiles"]].append(r)
        prepared = []
        for candidate in rows:
            sources = [r for smi in {candidate["smiles"],candidate["raw_smiles"]} for r in grouped.get(smi, [])]
            try:
                found, reason = migrate.recover_one(c, candidate, rows=sources)
            except (ValueError, TypeError, KeyError, SyntaxError) as exc:
                found, reason = None, "parse_error:" + str(exc)
            prepared.append((candidate, found, reason))
        outcomes = {}
        def apply(candidate, found, reason):
            current = c.execute("SELECT state FROM candidates WHERE smiles=?", (candidate["smiles"],)).fetchone()
            if current[0] != "awaiting_source":
                return
            c.execute("SAVEPOINT quality_candidate")
            try:
                mid, state = migrate.import_recovered(c, candidate, found) if found else (candidate["molecule_id"], "source_unresolved")
                pair = db.dumps({"h_source_id": found[1][0]["id"], "c_source_id": found[2][0]["id"]}) if found else None
                c.execute("UPDATE candidates SET state=?,reason=?,molecule_id=?,source_pair_json=? WHERE smiles=?", (state,reason,mid,pair,candidate["smiles"]))
            except (ValueError, TypeError, KeyError) as exc:
                c.execute("ROLLBACK TO quality_candidate")
                state = "source_unresolved"
                c.execute("UPDATE candidates SET state=?,reason=? WHERE smiles=?", (state,str(exc),candidate["smiles"]))
            finally:
                c.execute("RELEASE quality_candidate")
            outcomes[candidate["smiles"]] = state
        commit_prepared(c, prepared, apply)
        decisions.update(outcomes.values())
        count += len(rows)
        if progress:
            progress(count, dict(decisions))
    return {"processed": count, "decisions": dict(decisions)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("priority", "full", "audit", "recover"))
    parser.add_argument("--db", type=Path, default=db.DEFAULT_DB)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()
    c = db.connect(args.db, timeout=.25)
    db.initialize(c)
    db.set_meta(c, "quality_export_enforced", True)
    c.commit()
    c.execute("PRAGMA cache_size=-131072")
    folder = artifacts_root(args.db)
    state = {"pid": os.getpid(), "state": "running", "started_at": db.now(), "policy": quality.POLICY}
    state_file = folder / "migration_job.json"
    lock_file = folder / "quality_job.lock"
    # OS-held lock is released automatically after crashes; no stale PID guessing.
    import msvcrt
    lock_handle = lock_file.open("a+b")
    lock_handle.seek(0)
    msvcrt.locking(lock_handle.fileno(), msvcrt.LK_NBLCK, 1)
    state_lock, stopped = threading.Lock(), threading.Event()
    def publish():
        with state_lock:
            state["updated_at"] = db.now()
            temporary = state_file.with_suffix(".json.tmp")
            temporary.write_text(json.dumps(state, ensure_ascii=False, indent=2), encoding="utf-8")
            temporary.replace(state_file)
    def heartbeat():
        while not stopped.wait(5):
            publish()
    worker = threading.Thread(target=heartbeat, daemon=True)
    results = {"policy": quality.POLICY, "thresholds": quality.THRESHOLDS, "started_at": db.now()}
    def phase(name, total, operation):
        start = time.perf_counter()
        state.update(phase=name, total=total, processed=0, rate_per_second=None, phase_eta_seconds=None)
        publish()
        def update(count, counts):
            elapsed = time.perf_counter()-start
            state.update(processed=count, outcomes=counts, elapsed_seconds=round(elapsed,2),
                         rate_per_second=round(count/elapsed,2), phase_eta_seconds=round(max(0,total-count)*elapsed/count) if count else None)
        result = operation(update)
        result["elapsed_seconds"] = round(time.perf_counter()-start, 3)
        results[name] = result
        (folder / ("quality_" + args.stage + "_report.json")).write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding="utf-8")
        print(name, db.dumps(result), flush=True)
    try:
        publish()
        worker.start()
        if args.stage == "priority" or (args.stage == "full" and db.cursor(c, "quality_priority:" + quality.POLICY) != "done"):
            n = c.execute("SELECT count(*) FROM candidates WHERE state='complemented'").fetchone()[0]
            phase("audit_complements", n, lambda cb: audit(c,"candidates.state='complemented'",cb))
            n = c.execute("SELECT count(*) FROM molecules JOIN candidates INDEXED BY candidate_molecule ON candidates.molecule_id=molecules.id WHERE molecules.origin='nmrshiftdb2' AND molecules.queue='waiting_spectrum'").fetchone()[0]
            phase("rematch_missing",n,lambda cb: rematch(c,cb))
            migrate.recovery_write(c, lambda: c.execute("INSERT INTO progress VALUES(?, 'done', ?) ON CONFLICT(stage) DO UPDATE SET cursor='done',updated_at=excluded.updated_at", ("quality_priority:" + quality.POLICY, db.now())))
        if args.stage in ("audit", "full"):
            n = c.execute("SELECT count(*) FROM candidates LEFT JOIN candidate_quality q USING(smiles) WHERE candidates.state='imported' AND (q.policy IS NULL OR q.policy!=?)", (quality.POLICY,)).fetchone()[0]
            phase("audit_imported",min(n,args.limit) if args.limit else n,lambda cb: audit(c,"candidates.state='imported'",cb,args.limit))
        if args.stage in ("recover", "full"):
            n = c.execute("SELECT count(*) FROM candidates WHERE state='awaiting_source'").fetchone()[0]
            phase("recover_remaining",min(n,args.limit) if args.limit else n,lambda cb: recover_batched(c,cb,args.limit))
        if args.stage == "full" and not args.limit:
            def final_report(callback):
                integrity = [r[0] for r in c.execute("PRAGMA quick_check")]
                foreign_keys = [tuple(r) for r in c.execute("PRAGMA foreign_key_check")]
                if integrity != ["ok"] or foreign_keys:
                    raise RuntimeError("Database integrity check failed")
                counts = db.stats(c)
                counts.update(integrity=integrity, foreign_key_violations=foreign_keys,
                              source_document_verified=False, automatic_pass_is_not_manual_validation=True)
                (folder / "quality_final_database_report.json").write_text(json.dumps(counts,ensure_ascii=False,indent=2),encoding="utf-8")
                callback(1,{"integrity":"ok"})
                return counts
            phase("integrity_report",1,final_report)
        state.update(state="finished",finished_at=db.now(),scope=args.stage,limited=bool(args.limit))
    except Exception as exc:
        state.update(state="failed",error=str(exc),failed_at=db.now())
        raise
    finally:
        stopped.set()
        worker.join(timeout=2)
        publish()
        c.close()
        lock_handle.close()


if __name__ == "__main__":
    main()
