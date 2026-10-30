#!/usr/bin/env python
"""Audit NMRShiftDB2–BMRB outputs against the strict acceptance invariants."""

from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import sqlite3
from collections import defaultdict
from pathlib import Path

import lmdb


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data" / "nmrshiftdb2_bmrb_enriched"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()
    output = Path(args.output_dir).expanduser().resolve()
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    report = json.loads((output / "enrichment_report.json").read_text(encoding="utf-8"))
    conn = sqlite3.connect(output / "bmrb_enrichment.sqlite")
    conn.row_factory = sqlite3.Row
    checks: dict[str, object] = {}
    failures: list[str] = []

    checks["sqlite_integrity"] = conn.execute("PRAGMA integrity_check").fetchone()[0]
    if checks["sqlite_integrity"] != "ok":
        failures.append("SQLite integrity check failed")

    source = Path(manifest["source_sdf"]["path"])
    checks["source_sha256"] = sha256_file(source)
    checks["source_unchanged"] = checks["source_sha256"] == manifest["source_sdf"]["sha256"]
    if not checks["source_unchanged"]:
        failures.append("Source SDF SHA256 changed")

    snapshot_failures = []
    for item in manifest["bmrb_files"]:
        path = output / "bmrb_snapshot" / item["name"]
        if not path.exists() or sha256_file(path) != item["sha256"]:
            snapshot_failures.append(item["name"])
    checks["bmrb_snapshot_hash_failures"] = snapshot_failures
    if snapshot_failures:
        failures.append("BMRB snapshot checksum mismatch")

    counts = dict(conn.execute("SELECT status,COUNT(*) FROM local_spectra GROUP BY status"))
    checks["status_counts"] = counts
    checks["record_counts"] = {
        "total": conn.execute("SELECT COUNT(*) FROM local_records").fetchone()[0],
        "valid": conn.execute("SELECT COUNT(*) FROM local_records WHERE valid=1").fetchone()[0],
        "invalid": conn.execute("SELECT COUNT(*) FROM local_records WHERE valid=0").fetchone()[0],
    }
    if checks["record_counts"] != {"total": 58187, "valid": 58029, "invalid": 158}:
        failures.append("Source record accounting differs from the accepted snapshot")

    jsonl_lines = 0
    provenance_errors = 0
    non_auto_transfer_errors = 0
    with (output / "records.jsonl").open("r", encoding="utf-8") as handle:
        for line in handle:
            jsonl_lines += 1
            record = json.loads(line)
            for spectrum in record["spectra"]:
                auto = spectrum["match"]["status"] == "A_auto"
                for value in spectrum["conditions"].values():
                    if not isinstance(value, dict) or value.get("transfer_status") != "filled_missing":
                        continue
                    if not auto:
                        non_auto_transfer_errors += 1
                    if not value.get("source_table") or not value.get("source_row_id"):
                        provenance_errors += 1
                for peak in spectrum["peaks"]:
                    for value in peak.get("provenance", {}).values():
                        if value.get("transfer_status") not in {"filled_missing", "audit_only_unit_unknown"}:
                            continue
                        if not auto:
                            non_auto_transfer_errors += 1
                        if not value.get("source_table") or not value.get("source_row_id"):
                            provenance_errors += 1
    checks["jsonl_lines"] = jsonl_lines
    checks["provenance_errors"] = provenance_errors
    checks["non_auto_transfer_errors"] = non_auto_transfer_errors
    if jsonl_lines != checks["record_counts"]["valid"]:
        failures.append("records.jsonl does not contain every valid source record")
    if provenance_errors or non_auto_transfer_errors:
        failures.append("Field provenance or A_auto transfer invariant failed")

    split_by_spectrum: dict[int, str] = {}
    lmdb_counts = {}
    for split in ("train", "valid", "test"):
        path = output / "lmdb_dataset" / f"{split}.lmdb"
        env = lmdb.open(str(path), subdir=False, readonly=True, lock=False)
        with env.begin() as txn:
            length_raw = txn.get(b"__len__")
            if length_raw is None:
                failures.append(f"{split}.lmdb has no __len__ key")
                length = 0
            else:
                length = int(length_raw)
            for index in range(length):
                payload = txn.get(str(index).encode("ascii"))
                if payload is None:
                    failures.append(f"{split}.lmdb missing item {index}")
                    break
                item = pickle.loads(payload)
                spectrum_id = int(item["local_spectrum_id"])
                if spectrum_id in split_by_spectrum:
                    failures.append(f"Spectrum {spectrum_id} appears in multiple splits")
                split_by_spectrum[spectrum_id] = split
        env.close()
        lmdb_counts[split] = length
    checks["lmdb_counts"] = lmdb_counts

    connectivity_splits: dict[str, set[str]] = defaultdict(set)
    for spectrum_id, split in split_by_spectrum.items():
        row = conn.execute(
            """SELECT r.connectivity_key FROM local_spectra s JOIN local_records r
               ON r.record_index=s.record_index WHERE s.id=?""", (spectrum_id,)
        ).fetchone()
        if row and row["connectivity_key"]:
            connectivity_splits[row["connectivity_key"]].add(split)
    leakage = {key: sorted(values) for key, values in connectivity_splits.items() if len(values) > 1}
    checks["connectivity_split_leakage"] = leakage
    if leakage:
        failures.append("Connectivity groups cross train/valid/test splits")

    checks["review_rows"] = max(sum(1 for _ in (output / "review_queue.csv").open(encoding="utf-8-sig")) - 1, 0)
    if checks["review_rows"] != counts.get("B_review", 0):
        failures.append("Review queue does not contain every B_review spectrum")
    checks["negative_control"] = report["negative_control"]
    if report["negative_control"]["shift_only_strict_matches"] != 0:
        failures.append("Shuffled negative control produced strict matches")

    partials = [str(path) for path in output.rglob("*.partial")]
    checks["partial_files"] = partials
    if partials:
        failures.append("Partial output files remain")

    result = {"passed": not failures, "failures": failures, "checks": checks}
    (output / "audit_report.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    conn.close()
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
