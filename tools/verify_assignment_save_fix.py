"""Apply query indexes and verify label preservation; production save probes roll back."""
import argparse
import hashlib
import json
import statistics
import tempfile
import time
from pathlib import Path
from contextlib import closing

from annotation_tool_v3 import core as db
from annotation_tool_v3.maintenance import read_connection


def annotation_digest(c):
    digest = hashlib.sha256()
    count = 0
    for row in c.execute("SELECT * FROM annotations ORDER BY spectrum_id,peak_id"):
        digest.update(db.dumps(list(row)).encode())
        digest.update(b"\n")
        count += 1
    return {"rows": count, "sha256": digest.hexdigest()}


def timed(operation, repeats=5):
    values = []
    for _ in range(repeats):
        started = time.perf_counter()
        operation()
        values.append(1000 * (time.perf_counter()-started))
    return {"milliseconds": values, "median_ms": statistics.median(values)}


class RollbackCommit:
    """Exercise all production save SQL without persisting changes to labels/revision."""
    def __init__(self, connection):
        self.connection = connection

    def __getattr__(self, name):
        return getattr(self.connection, name)

    def commit(self):
        self.connection.rollback()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=db.DEFAULT_DB)
    parser.add_argument("--backup", type=Path, required=True)
    parser.add_argument("--molecule", default="M0f34e064bf2200828e6bdf1e")
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if not args.backup.is_file() or args.backup.resolve() == args.db.resolve():
        raise ValueError("A separate existing backup is required")
    result = {"created_at": db.now(), "backup": str(args.backup.resolve()), "molecule": args.molecule}
    with closing(read_connection(args.backup)) as backup:
        result["backup_annotations"] = annotation_digest(backup)
    with closing(db.connect(args.db)) as c:
        result["before_annotations"] = annotation_digest(c)
        before = db.record(c, args.molecule)
        result["before_record_read"] = timed(lambda: db.record(c, args.molecule), 3)
        print("Applying query indexes; no annotation updates", flush=True)
        db.initialize(c)
        result["pair_query_plan"] = [list(r) for r in c.execute("EXPLAIN QUERY PLAN SELECT * FROM pairs WHERE molecule_id=? ORDER BY id", (args.molecule,))]
        result["after_record_read"] = timed(lambda: db.record(c, args.molecule))
        payload = {"revision": before["revision"], "task_done": bool(before["task_done"]), "note": before["note"],
                   "assignments": [{"spectrum_id": s["id"], "peak_id": p["peak_id"], "status": p["assignment"]["status"],
                                    "atom_indices": p["assignment"]["atom_indices"], "note": p["assignment"]["note"]}
                                   for s in before["spectra"] for p in s["peaks"] if p["assignment"].get("origin") == "manual"]}
        result["production_save_rolled_back"] = timed(lambda: db.save_annotations(RollbackCommit(c), args.molecule, payload))
        result["record_unchanged"] = before == db.record(c, args.molecule)
        result["after_annotations"] = annotation_digest(c)
        result["annotations_unchanged"] = result["before_annotations"] == result["after_annotations"]
        result["backup_labels_match"] = result["backup_annotations"] == result["after_annotations"]
        result["revision_after"] = db.record(c, args.molecule)["revision"]
        if not result["record_unchanged"] or not result["annotations_unchanged"]:
            raise RuntimeError("Database changed during verification; inspect concurrent user activity, never auto-restore")
        # Verify real durable commits too, but only in a small disposable copy of
        # this molecule. Its timings are not full-production commit benchmarks.
        with tempfile.TemporaryDirectory(prefix="nmr-save-benchmark-") as folder:
            with closing(db.connect(Path(folder)/"test.sqlite")) as target:
                db.initialize(target)
                queries = {
                    "molecules": ("SELECT * FROM molecules WHERE id=?", (args.molecule,)),
                    "spectra": ("SELECT * FROM spectra WHERE molecule_id=?", (args.molecule,)),
                    "annotations": ("SELECT a.* FROM annotations a JOIN spectra s ON a.spectrum_id=s.id WHERE s.molecule_id=?", (args.molecule,)),
                    "pairs": ("SELECT * FROM pairs WHERE molecule_id=?", (args.molecule,)),
                }
                for table, (query, values) in queries.items():
                    for row in c.execute(query, values):
                        target.execute(f"INSERT INTO {table} VALUES({','.join('?' for _ in row)})", tuple(row))
                target.commit()
                def real_save():
                    result_row = db.save_annotations(target, args.molecule, payload)
                    payload["revision"] = result_row["revision"]
                result["small_copy_durable_save"] = timed(real_save)
    args.report.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
