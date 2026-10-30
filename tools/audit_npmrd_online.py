"""Recheck NP-MRD's submitted Actual columns against locally mapped assignments.

This is a source-level audit, not a re-assignment of spectra or paper-by-paper
verification. All downloads and derived results go into a separate dated folder.
"""
from __future__ import annotations

import argparse
import collections
import concurrent.futures
import csv
import hashlib
import io
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.prepare_npmrd_experimental import page_info, parse_user_table, write_json, write_jsonl


def load_rows(path):
    return [json.loads(line) for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]


def signature(rows):
    return sorted((r["symbol"], r["source_atom_id"], r["shift_ppm"]) for r in rows)


def parse_metadata(text):
    """Read explicit Actual/Predicted headings, not positional headerless guesses."""
    lines = text.splitlines()
    head = next((i for i, line in enumerate(lines) if line.startswith("Atom Type,Atom No,")), None)
    if head is None:
        raise ValueError("metadata_assignment_header_missing")
    metadata = dict(line.split(": ", 1) for line in lines[:head] if ": " in line)
    reader = csv.DictReader(io.StringIO("\n".join(lines[head:])))
    required = {"Atom Type", "Atom No", "Actual Shift(ppm)", "Predicted Shift(ppm)"}
    if not required.issubset(reader.fieldnames or []):
        raise ValueError("actual_predicted_columns_missing")
    actual, predicted = [], []
    for row in reader:
        if row["Atom Type"] not in {"H", "C"}:
            continue
        for category, result in (("Actual", actual), ("Predicted", predicted)):
            buf = io.StringIO()
            csv.writer(buf).writerow([row["Atom Type"], row["Atom No"], row[f"{category} Shift(ppm)"],
                                     row.get(f"{category} Multiplet Type"), row.get(f"{category} J coupling(Hz)")])
            parsed, missing, errors = parse_user_table(buf.getvalue())
            if errors:
                raise ValueError(f"{category.lower()}_metadata_rows_invalid: {errors}")
            result.extend(parsed)
    return {"metadata": metadata, "columns": reader.fieldnames, "actual": actual, "predicted": predicted,
            "spectral_files": lines[next((i for i, line in enumerate(lines) if line.startswith("Spectral File")), len(lines)):],
            "actual_j_unit_hz_explicit": "Actual J coupling(Hz)" in reader.fieldnames}


def download(url, output, name):
    raw = output / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    target = raw / name
    stamp = raw / (name + ".download.json")
    if target.exists() and stamp.exists():
        meta = json.loads(stamp.read_text(encoding="utf-8"))
        if meta["url"] == url and meta["sha256"] == hashlib.sha256(target.read_bytes()).hexdigest():
            return target.read_text(encoding="utf-8-sig", errors="replace"), meta
    for attempt in range(2):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (NMR source audit)", "Accept": "*/*"})
            with urllib.request.urlopen(request, timeout=25) as response:
                data = response.read(3 * 1024 * 1024 + 1)
                if len(data) > 3 * 1024 * 1024:
                    raise ValueError("unexpectedly_large_source_document")
                meta = {"url": url, "resolved_url": response.url, "http_status": response.status,
                        "checked_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                        "content_type": response.headers.get("Content-Type"), "bytes": len(data),
                        "sha256": hashlib.sha256(data).hexdigest()}
            target.write_bytes(data)
            write_json(stamp, meta)
            return data.decode("utf-8-sig", errors="replace"), meta
        except (OSError, urllib.error.URLError):
            if attempt:
                raise
            time.sleep(0.5)


def inspect(record, output, existing_mapped):
    sid = record["submission_id"]
    member = record["record_id"] + ".txt"
    url = f"https://np-mrd.org/natural_products/{sid}/submission_spectrum_view?spectrum_type=simulation"
    result = {"record_id": record["record_id"], "np_mrd_id": record["np_mrd_id"],
              "existing_mapped": existing_mapped, "dedup_status": record["dedup"]["status"],
              "previous_reason": record.get("reason"), "page_url": url, "checks": {}, "errors": []}
    local = [r for group in record["assignments"].values() for r in group] if existing_mapped else record["assignment_rows"]
    checks = result["checks"]
    try:
        text, stamp = download(url, output, f"submission_{sid}.html")
        info = page_info(text)
        result["page"] = {"title": info["title"], "conditions": info["conditions_raw"], "references": info["references"], "download": stamp}
        checks["same_compound_page"] = info["conditions_raw"].get("NP-MRD ID") == record["np_mrd_id"]
        checks["official_experimental_or_simulated_type"] = bool(re.search(r",\s*(simulated|experimental)\)", info["title"] or "", re.I))
        link = next((link for link in info["links"] if Path(urllib.parse.urlparse(link).path).name == member), None)
        checks["exact_archive_member_linked"] = link is not None
        if link:
            table, table_stamp = download(urllib.parse.urljoin(url, link), output, f"{record['record_id']}.txt")
            live, missing, errors = parse_user_table(table)
            result["user_table_download"] = table_stamp
            checks["live_user_table_matches_local"] = not errors and signature(live) == signature(local)
    except Exception as exc:
        result["errors"].append(f"page_or_table: {type(exc).__name__}: {exc}")
    try:
        meta_url = f"https://np-mrd.org/submissions/{sid}/download_submitted_chemical_shift_data"
        text, stamp = download(meta_url, output, f"submission_{sid}.metadata.txt")
        parsed = parse_metadata(text)
        result["submitted_metadata"] = parsed
        result["metadata_download"] = stamp
        checks["same_compound_metadata"] = parsed["metadata"].get("NP-MRD ID") == record["np_mrd_id"]
        checks["actual_shifts_match_local"] = signature(parsed["actual"]) == signature(local)
        checks["actual_predicted_differ"] = signature(parsed["actual"]) != signature(parsed["predicted"])
        checks["explicit_actual_shift_column"] = True
    except Exception as exc:
        result["errors"].append(f"metadata: {type(exc).__name__}: {exc}")
    required = ("same_compound_page", "official_experimental_or_simulated_type", "exact_archive_member_linked",
                "live_user_table_matches_local", "same_compound_metadata", "actual_shifts_match_local", "explicit_actual_shift_column")
    result["passed_source_level_audit"] = existing_mapped and all(checks.get(key, False) for key in required)
    result["paper_values_independently_verified"] = False
    return result


def run(source, output):
    output.mkdir(parents=True, exist_ok=True)
    candidates = load_rows(source / "novel_exact.jsonl")
    review = load_rows(source / "review_queue.jsonl")
    examples = {}
    for row in review:
        examples.setdefault(row["reason"], row)
    tasks = [(r, True) for r in candidates] + [(r, False) for r in examples.values()]
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        jobs = [pool.submit(inspect, row, output, mapped) for row, mapped in tasks]
        for job in concurrent.futures.as_completed(jobs):
            results.append(job.result())
            if len(results) % 10 == 0:
                print(f"[online-audit] {len(results)}/{len(jobs)} passed={sum(r['passed_source_level_audit'] for r in results)}", flush=True)
                write_jsonl(output / "audit.partial.jsonl", sorted(results, key=lambda r: r["record_id"]))
    results.sort(key=lambda r: r["record_id"])
    write_jsonl(output / "audit.jsonl", results)
    passed = {r["record_id"]: r for r in results if r["passed_source_level_audit"]}
    checked = []
    for record in candidates:
        if record["record_id"] not in passed:
            continue
        row = json.loads(json.dumps(record))
        audit = passed[row["record_id"]]
        row["data_origin"] = "npmrd_actual_assignment_column_experimental_origin_per_official_classification"
        row["source_verification"] = {"level": "official_database_actual_column_verified",
            "checked_at": audit["metadata_download"]["checked_at"],
            "official_definition": "https://doi.org/10.1093/nar/gkae1067",
            "metadata_url": audit["metadata_download"]["url"],
            "metadata_sha256": audit["metadata_download"]["sha256"],
            "paper_values_independently_verified": False,
            "raw_spectrum_reassigned": False}
        checked.append(row)
    strict = [r for r in checked if r["dedup"]["status"] == "novel_connectivity"]
    write_jsonl(output / "novel_exact.source_checked.jsonl", checked)
    write_jsonl(output / "novel_connectivity.source_checked.jsonl", strict)

    def stats(rows):
        return {"records": len(rows), "molecules": len({r["structure"]["inchikey"] for r in rows}),
                "h_records": sum(bool(r["assignments"]["1H"]) for r in rows),
                "c_records": sum(bool(r["assignments"]["13C"]) for r in rows),
                "h_atom_assignments": sum(len(r["assignments"]["1H"]) for r in rows),
                "c_atom_assignments": sum(len(r["assignments"]["13C"]) for r in rows)}
    report = {"created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
              "scope": "All 82 previously mapped exact-novel candidates; one example for each prior review reason. No new review record is promoted.",
              "candidate_records_checked": len(candidates), "review_examples_checked": len(examples),
              "novel_exact": stats(checked), "novel_connectivity": stats(strict),
              "candidate_failures": [r for r in results if r["existing_mapped"] and not r["passed_source_level_audit"]],
              "review_examples": [r for r in results if not r["existing_mapped"]],
              "limitations": ["Official source-level experimental-origin evidence; not paper-by-paper numeric validation.",
                              "Existing atom mappings and InChIKey dedup results are reused.",
                              "Only representative review records were re-requested; 440 prior review records are not all reprocessed."]}
    write_json(output / "report.json", report)
    print(json.dumps({k: v for k, v in report.items() if k not in {"candidate_failures", "review_examples"}}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=Path("data/np_mrd_experimental"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/np_mrd_experimental/online_audit_20260916"))
    args = parser.parse_args()
    run(args.source_dir, args.output_dir)
