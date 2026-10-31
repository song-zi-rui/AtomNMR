"""Read-only audit of recorded DeepSeek usage. Never imports the API client or keys.

The compact representations below are OFFLINE size experiments, not production
request formats. Their character counts are not provider token estimates.
"""
from __future__ import annotations

import argparse
import collections
import copy
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import statistics


ROOT = Path(__file__).resolve().parents[1]
HARD = {"assigned", "equivalent_group", "overlap"}


def dumps(value):
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), allow_nan=False)


def parsed(value, fallback=None):
    return json.loads(value) if value else fallback


def distribution(values):
    if not values:
        return {"n": 0}
    ordered = sorted(values)
    return {"n": len(values), "sum": sum(values), "mean": statistics.mean(values),
            "median": statistics.median(values),
            "p90": ordered[max(0, math.ceil(len(ordered) * .9) - 1)], "max": max(values)}


def usage_parts(usage):
    if not usage or not isinstance(usage.get("total_tokens"), int):
        return None
    reasoning = (usage.get("completion_tokens_details") or {}).get("reasoning_tokens")
    assert isinstance(reasoning, int), "Cannot infer reasoning token count for this provider response"
    result = {k: usage[k] for k in ("prompt_tokens", "completion_tokens", "total_tokens")}
    result.update(reasoning_tokens=reasoning,
                  final_tokens=usage["completion_tokens"] - reasoning,
                  cache_hit_tokens=usage.get("prompt_cache_hit_tokens", 0),
                  cache_miss_tokens=usage.get("prompt_cache_miss_tokens", 0))
    assert result["total_tokens"] == result["prompt_tokens"] + result["completion_tokens"]
    assert result["completion_tokens"] >= reasoning
    return result


def pack_rows(rows):
    """Common columns once, rare fields in sparse extras; missing != explicit null."""
    if not rows:
        return {"columns": [], "rows": []}
    common = [k for k in rows[0] if all(k in row for row in rows)]
    result = {"columns": common, "rows": [[r[k] for k in common] for r in rows]}
    extras = {str(i): {k: v for k, v in r.items() if k not in common}
              for i, r in enumerate(rows) if any(k not in common for k in r)}
    if extras:
        result["extras"] = extras
    return result


def unpack_rows(table):
    return [dict(zip(table["columns"], row)) | table.get("extras", {}).get(str(i), {})
            for i, row in enumerate(table["rows"])]


def compact_context(context):
    """Retain all chemical values and atom indices; spectrum/peak aliases are local."""
    result = copy.deepcopy(context)
    result["atoms"] = pack_rows(result["atoms"])
    result["bonds"] = pack_rows(result["bonds"])
    spectra_map, peak_map = {}, {}
    counts = collections.Counter()
    for spectrum in result["spectra"]:
        nucleus = "C" if spectrum["nucleus"] == "13C" else "H"
        counts[nucleus] += 1
        original = spectrum["spectrum_id"]
        alias = nucleus + str(counts[nucleus])
        spectra_map[original] = alias
        spectrum["spectrum_id"] = alias
        peak_map[original] = {}
        for i, peak in enumerate(spectrum["peaks"]):
            peak_map[original][peak["peak_id"]] = "p" + str(i)
            peak["peak_id"] = "p" + str(i)
        spectrum["peaks"] = pack_rows(spectrum["peaks"])
    for pair in result["pairs"]:
        pair["c_id"] = spectra_map[pair["c_id"]]
        pair["h_id"] = spectra_map[pair["h_id"]]
    # Check reversibility for every actual input, including rare optional peak fields.
    restored = copy.deepcopy(result)
    restored["atoms"], restored["bonds"] = unpack_rows(restored["atoms"]), unpack_rows(restored["bonds"])
    inverse = {v: k for k, v in spectra_map.items()}
    for spectrum in restored["spectra"]:
        spectrum["spectrum_id"] = inverse[spectrum["spectrum_id"]]
        spectrum["peaks"] = unpack_rows(spectrum["peaks"])
        peaks_inverse = {v: k for k, v in peak_map[spectrum["spectrum_id"]].items()}
        for peak in spectrum["peaks"]:
            peak["peak_id"] = peaks_inverse[peak["peak_id"]]
    for pair in restored["pairs"]:
        pair["c_id"], pair["h_id"] = inverse[pair["c_id"]], inverse[pair["h_id"]]
    assert restored == context
    return result, spectra_map, peak_map


def compact_output(output, spectra_map, peak_map, omit_definite_reasons=False):
    grouped = collections.defaultdict(list)
    for assignment in output["assignments"]:
        row = copy.deepcopy(assignment)
        sid = row.pop("spectrum_id")
        row["peak_id"] = peak_map[sid][row["peak_id"]]
        if omit_definite_reasons and row["status"] in HARD:
            row.pop("brief_reason", None)
        grouped[spectra_map[sid]].append(row)
    return {"summary": output["summary"], "warnings": output["warnings"],
            "assignments": {sid: pack_rows(rows) for sid, rows in grouped.items()}}


def audit(path, job_id=None):
    with sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True) as c:
        c.row_factory = sqlite3.Row
        c.execute("PRAGMA query_only=ON")
        c.execute("BEGIN")
        job = c.execute("SELECT * FROM auto_jobs WHERE id=?", (job_id,)).fetchone() if job_id else c.execute(
            "SELECT * FROM auto_jobs ORDER BY created_at DESC,id DESC LIMIT 1").fetchone()
        if not job:
            raise ValueError("No matching continuous annotation job")
        rows = [dict(r) for r in c.execute("SELECT * FROM auto_items WHERE job_id=? ORDER BY id", (job["id"],))]
        legacy = [dict(r) for r in c.execute("SELECT state,usage_json FROM items")]
    config = json.loads(job["config_json"])
    states, failures, refs = collections.Counter(), collections.Counter(), collections.Counter()
    parts, peak_fields, statuses = collections.Counter(), collections.Counter(), collections.Counter()
    items, measured, sizes, output_sizes = [], [], [], []
    complexity = collections.defaultdict(list)
    for row in rows:
        context = json.loads(row["input_json"])
        usage = usage_parts(parsed(row["usage_json"], {}))
        diagnostics = parsed(row["diagnostics_json"], {}) or {}
        output = parsed(row["output_json"])
        states[row["state"]] += 1
        for reference in json.loads(row["reference_json"]):
            refs[reference["source"]] += 1
        for key, value in context.items():
            parts[key] += len(dumps({key: value})) - 2
        for spectrum in context["spectra"]:
            for peak in spectrum["peaks"]:
                peak_fields.update(peak.keys())
        compact, spectra_map, peak_map = compact_context(context)
        sizes.append({"original": len(dumps(context)), "compact": len(dumps(compact))})
        if output:
            wire = {k: output[k] for k in ("summary", "warnings", "assignments")}
            output_sizes.append({"original": len(dumps(wire)),
                                 "compact": len(dumps(compact_output(output, spectra_map, peak_map))),
                                 "compact_without_definite_reasons": len(dumps(compact_output(output, spectra_map, peak_map, True)))})
            statuses.update(a["status"] for a in output["assignments"])
        error_kind = None
        if row["state"] == "error":
            error_kind = "output_length" if diagnostics.get("finish_reason") == "length" else (
                "validation_after_stop" if diagnostics.get("finish_reason") == "stop" else "no_usage_returned")
            failures[error_kind] += 1
        n = len(context["atoms"])
        item = {"id": row["id"], "molecule_id": row["molecule_id"], "formula": row["formula"],
                "state": row["state"], "atoms": n, "spectra": len(context["spectra"]),
                "peaks": sum(len(s["peaks"]) for s in context["spectra"]),
                "error_kind": error_kind, "finish_reason": diagnostics.get("finish_reason"),
                "error": row["error"], "usage": usage, "elapsed": row["elapsed"],
                "input_characters": sizes[-1]}
        items.append(item)
        if usage:
            measured.append(item)
            complexity["<=15" if n <= 15 else "16-25" if n <= 25 else "26-35" if n <= 35 else ">=36"].append(item)
    metrics = {key: distribution([r["usage"][key] for r in measured]) for key in (
        "prompt_tokens", "completion_tokens", "reasoning_tokens", "final_tokens", "total_tokens", "cache_hit_tokens", "cache_miss_tokens")}
    totals = {k: v["sum"] for k, v in metrics.items()}
    grouped = {}
    for name, group in {"done": [r for r in measured if r["state"] == "done"],
                        **{name: [r for r in measured if r["error_kind"] == name] for name in failures}}.items():
        grouped[name] = {"n": len(group), **{key: sum(r["usage"][key] for r in group) for key in totals}}
    size_totals = {k: sum(v[k] for v in sizes) for k in ("original", "compact")}
    output_totals = {k: sum(v[k] for v in output_sizes) for k in ("original", "compact", "compact_without_definite_reasons")}
    legacy_usage = [usage_parts(parsed(r["usage_json"], {})) for r in legacy]
    legacy_total = sum(u["total_tokens"] for u in legacy_usage if u)
    return {"generated_at": dt.datetime.now(dt.timezone.utc).isoformat(), "source": str(path.resolve()),
            "job": {k: job[k] for k in ("id", "state", "scope", "processed", "success", "failures", "applied_peaks")},
            "config": {k: config.get(k) for k in ("model", "thinking", "max_tokens", "timeout", "reasoning_effort")},
            "system_prompt_characters": len(config.get("system_prompt", "")),
            "row_count": len(rows), "usage_count": len(measured), "unreported_usage_count": len(rows) - len(measured),
            "states": dict(states), "failures": dict(failures), "reference_labels": dict(refs),
            "totals": totals, "distributions": metrics, "groups": grouped,
            "input_character_components": dict(parts), "input_characters": size_totals,
            "output_characters": output_totals, "output_status_counts": dict(statuses),
            "peak_field_counts": dict(peak_fields),
            "complexity": {name: {"n": len(group), "mean_reasoning": statistics.mean(r["usage"]["reasoning_tokens"] for r in group),
                                  "mean_total": statistics.mean(r["usage"]["total_tokens"] for r in group),
                                  "length_failures": sum(r["error_kind"] == "output_length" for r in group)} for name, group in complexity.items()},
            "legacy": {"rows": len(legacy), "usage_count": sum(u is not None for u in legacy_usage), "recorded_total": legacy_total},
            "snapshot_sha256": hashlib.sha256(dumps(rows).encode()).hexdigest(),
            "notes": ["API usage is authoritative; absent usage is unknown, not zero.",
                      "final_tokens = completion_tokens - reasoning_tokens; truncated requests include incomplete final output.",
                      "Character compression is an offline representation experiment, not a measured token saving.",
                      "The compact input retains all values and is round-trip checked; protocol instructions and local ID maps are excluded from its size.",
                      "Output character analysis uses normalized saved JSON; original raw provider content is not retained.",
                      "No API calls, key access, database writes, or production behavior changes."], "items": items}


def render_svg(report):
    totals = report["totals"]
    segments = [("输入", totals["prompt_tokens"], "#417da0"),
                ("内部思考", totals["reasoning_tokens"], "#c68135"),
                ("最终输出", totals["final_tokens"], "#258777")]
    elements = ['<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="300" viewBox="0 0 1120 300">',
                '<rect width="1120" height="300" rx="16" fill="#f6f8fa"/>',
                '<g font-family="Microsoft YaHei,Arial,sans-serif" fill="#183442">',
                '<text x="36" y="48" font-size="24" font-weight="bold">当前自动标注：token 消耗构成</text>',
                f'<text x="36" y="80" font-size="15">{report["row_count"]} 次请求，{report["usage_count"]} 次有用量；合计 {totals["total_tokens"]:,} token（未回报部分未知）</text>']
    x = 36
    for label, value, color in segments:
        width = 1048 * value / totals["total_tokens"]
        elements += [f'<rect x="{x:.2f}" y="110" width="{width:.2f}" height="56" fill="{color}"/>',
                     f'<text x="{x + width / 2:.2f}" y="146" font-size="19" text-anchor="middle" fill="white">{value / totals["total_tokens"]:.1%}</text>']
        x += width
    for i, (label, value, color) in enumerate(segments):
        x = 36 + i * 350
        elements += [f'<rect x="{x}" y="194" width="14" height="14" fill="{color}"/>',
                     f'<text x="{x + 24}" y="207" font-size="17">{label}：{value:,}</text>']
    elements += ['<text x="36" y="254" font-size="15">数据库查询、RDKit 处理、校验、保存、统计与页面刷新：本地操作，API token 为 0。</text>', '</g></svg>']
    return "\n".join(elements)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=ROOT / "data/assignment_database_v1/databases/ai_annotations.sqlite")
    parser.add_argument("--job-id")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "docs/ai_token_audit_20261028")
    args = parser.parse_args()
    report = audit(args.db, args.job_id)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "token_usage.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (args.output_dir / "token_flow.svg").write_text(render_svg(report), encoding="utf-8")
    print(json.dumps({k: report[k] for k in ("job", "totals", "failures", "input_characters", "output_characters", "complexity")}, ensure_ascii=True, indent=2))


if __name__ == "__main__":
    main()
