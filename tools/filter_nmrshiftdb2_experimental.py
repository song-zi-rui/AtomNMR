#!/usr/bin/env python
"""Audit every SDF spectrum; export experimental-only copies without renumbering atoms.

Source status is inferred from the indexed metadata in the NMRShiftDB2 export,
not from a spectrum's ordinal, chemical shifts, or assignment-method wording.
This is a source filter, not independent verification of the original experiment.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem, RDLogger, rdBase


POLICY = "nmrshiftdb2-experimental-v1"
ROOT = Path(__file__).resolve().parents[1]
SPECTRUM = re.compile(r"^Spectrum (.+) (\d+)$")
RAW_DATA = re.compile(r"^rawdata (.+) (\d+)$")
HEADER = re.compile(rb"^>[^\r\n]*<([^>]+)>[^\r\n]*\r?\n", re.M)
INDEX = re.compile(r"(?:^|\s)(\d+):")
EXPERIMENT_FIELDS = ("Solvent", "Temperature [K]", "Field Strength [MHz]", "Assignment Method")
COMPUTATION_FIELDS = (
    "Program", "NMRProgram", "NMRMethod", "NMRBasisSet", "NMRModel", "NMRLocalis",
    "NMRStandard", "GeomProgram", "GeomMethod", "GeomModel", "GeomBasisSet",
)
INDEXED_FIELDS = set(EXPERIMENT_FIELDS + COMPUTATION_FIELDS + ("Machine",))
NON_NMR_TYPES = {"Cyclic voltammetry"}
OUTPUTS = {
    "experimental_all_nmr": "experimental_all_nmr.sd",
    "experimental_ch_usable": "experimental_ch.sd",
    "experimental_ch_both_usable": "experimental_ch_both.sd",
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def indexed_items(value):
    matches = list(INDEX.finditer(value))
    result = []
    for i, match in enumerate(matches):
        end = matches[i + 1].start() if i + 1 < len(matches) else len(value)
        result.append((match[1], value[match.end():end].strip()))
    return result


def indexed(value):
    result = {}
    for sid, item in indexed_items(value):
        if sid in result:
            raise ValueError(f"Duplicate indexed metadata ID: {sid}")
        result[sid] = item
    return result


def iter_records(path):
    """Yield original bytes, one-based record index and source start line."""
    lines = []
    start = 1
    index = 0
    with Path(path).open("rb") as handle:
        for line_number, line in enumerate(handle, 1):
            lines.append(line)
            if line.strip() == b"$$$$":
                index += 1
                yield index, start, b"".join(lines)
                lines = []
                start = line_number + 1
    if b"".join(lines).strip():
        raise ValueError(f"Unterminated SDF record after line {start}")


def parse_record(raw):
    matches = list(HEADER.finditer(raw))
    if not matches:
        raise ValueError("Record has no properties; cannot audit safely")
    props, sections = {}, []
    for i, match in enumerate(matches):
        end = matches[i + 1].start() if i + 1 < len(matches) else raw.rfind(b"$$$$")
        name = match[1].decode("utf-8")
        if name in props:
            raise ValueError(f"Duplicate SDF property: {name}")
        props[name] = raw[match.end():end].decode("utf-8").strip()
        sections.append((name, raw[match.start():end]))
    return raw[:matches[0].start()], props, sections


def classify(props, sid, spectrum_type):
    if spectrum_type in NON_NMR_TYPES:
        return "non_nmr", {"reason": "not_an_nmr_experiment"}
    evidence = {}
    malformed = []
    for name in INDEXED_FIELDS:
        if name not in props:
            continue
        values = indexed_items(props[name])
        if not values and props[name].strip():
            malformed.append(name)
        local = [value for key, value in values if key == sid]
        if local:
            evidence[name] = local[0] if len(local) == 1 else local
    if malformed:
        return "unknown", {"reason": "unindexed_metadata", "fields": sorted(malformed)}
    calculation = any(name in evidence for name in COMPUTATION_FIELDS)
    experiment = any(name in evidence for name in EXPERIMENT_FIELDS)
    if calculation and experiment:
        return "conflict", evidence
    if calculation:
        return "computed", evidence
    if experiment:
        return "experimental", evidence
    return "unknown", evidence


def filter_record(raw, keep_fields, source_index, source_sha):
    molblock, props, sections = parse_record(raw)
    keep_ids = {SPECTRUM.fullmatch(name)[2] for name in keep_fields}
    pieces = [molblock]
    def property_bytes(name, value):
        # Some source properties start with a spurious blank line/CR. RDKit
        # otherwise silently stops reading the remaining properties. Normalize
        # only property framing, never the molecule block or stripped values.
        return f"> <{name}>\n{value}\n\n".encode("utf-8")

    for name, section in sections:
        match = SPECTRUM.fullmatch(name)
        raw_match = RAW_DATA.fullmatch(name)
        if match:
            if name in keep_fields:
                pieces.append(property_bytes(name, props[name]))
        elif raw_match:
            if f"Spectrum {raw_match[1]} {raw_match[2]}" in keep_fields:
                pieces.append(property_bytes(name, props[name]))
        elif name in INDEXED_FIELDS:
            values = [(k, v) for k, v in indexed_items(props[name]) if k in keep_ids]
            if values:
                value = " ".join(f"{k}:{v}" for k, v in values)
                pieces.append(property_bytes(name, value))
        elif name not in {"Experimental filter policy", "Original SDF record index", "Original SDF SHA256"}:
            pieces.append(property_bytes(name, props[name]))
    for name, value in (
        ("Experimental filter policy", POLICY),
        ("Original SDF record index", str(source_index)),
        ("Original SDF SHA256", source_sha),
    ):
        pieces.append(property_bytes(name, value))
    pieces.append(b"$$$$\n")
    result = b"".join(pieces)
    new_block, new_props, _ = parse_record(result)
    assert molblock == new_block, "Molecular block changed"
    assert {name for name in new_props if SPECTRUM.fullmatch(name)} == set(keep_fields)
    for name in keep_fields:
        assert props[name] == new_props[name], "Spectrum values or atom indices changed"
        match = SPECTRUM.fullmatch(name)
        assert classify(new_props, match[2], match[1])[0] == "experimental"
    return result


def validate_molecule(molblock):
    try:
        mol = Chem.MolFromMolBlock(molblock.decode("utf-8"), sanitize=False, removeHs=False)
        if mol is None:
            return None, "mol_parse_failed"
        Chem.SanitizeMol(mol)
        return mol, None
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"


def signal_quality(value, mol, nucleus):
    """Keep source occurrences and reject spectrum-level identity conflicts.

    This filter assesses assignment usability, never rewrites the source SDF.
    """
    from tools.shiftdb_atom_identity import identity_for_current, parse_fragments, validate_fragments
    if mol is not None:
        fragments = parse_fragments(value, nucleus, identity_for_current(mol))
        validation = validate_fragments(fragments, nucleus)
        blocked = validation['status'] != 'identity_validated'
        usable = [f for f in fragments if f['shift_ppm'] is not None and f['assignments'] and not blocked]
        return {'raw_assignment_records':len(fragments),
                'nonfinite_assignment_records':sum(f['shift_ppm'] is None and f['raw_signal'].split(';')[0].strip().lower() in {'nan','inf','-inf','+inf','infinity'} for f in fragments),
                'finite_shift_records':sum(f['shift_ppm'] is not None for f in fragments),
                'usable_assignment_records':len(usable),
                'unusable_assignment_records':len(fragments)-len(usable),
                'unique_usable_peaks':len({(f['shift_text'],f['source_signal_code']) for f in usable}),
                'assignment_identity_status':validation['status'],
                'assignment_identity_issue_counts':validation['issue_counts']}
    count = Counter()
    unique = set()
    for raw in filter(str.strip, value.split("|")):
        count["raw_assignment_records"] += 1
        parts = raw.split(";")
        if len(parts) < 3:
            count["unusable_assignment_records"] += 1
            continue
        try:
            shift = float(parts[0].strip())
        except ValueError:
            count["unusable_assignment_records"] += 1
            continue
        if not math.isfinite(shift):
            count["nonfinite_assignment_records"] += 1
            continue
        count["finite_shift_records"] += 1
        # No structure was parsed: record the value, never assert an assignment.
    count["unique_usable_peaks"] = len(unique)
    return dict(count)


class Pool:
    def __init__(self):
        self.counts = Counter()
        self.types = Counter()
        self.keys = set()
        self.connectivity = set()
        self.ids = set()

    def add(self, props, spectra, include_empty=False):
        if not spectra and not include_empty:
            return
        self.counts["molecule_records"] += 1
        self.counts["spectra"] += len(spectra)
        types = {row["spectrum_type"] for row in spectra}
        self.types.update(row["spectrum_type"] for row in spectra)
        self.counts["molecules_with_13c"] += int("13C" in types)
        self.counts["molecules_with_1h"] += int("1H" in types)
        self.counts["molecules_with_both_ch"] += int({"13C", "1H"} <= types)
        self.counts["molecules_with_either_ch"] += int(bool({"13C", "1H"} & types))
        key = props.get("INChI key", "").strip()
        if key:
            self.keys.add(key)
            self.connectivity.add(key.split("-")[0])
            if not re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", key):
                self.counts["malformed_source_inchi_key_records"] += 1
        else:
            self.counts["missing_inchi_key"] += 1
        self.ids.add(props.get("nmrshiftdb2 ID", ""))

    def summary(self):
        return {**dict(self.counts), "by_spectrum_type": dict(sorted(self.types.items())),
                "unique_source_inchi_keys": len(self.keys),
                "unique_source_connectivity_keys": len(self.connectivity),
                "unique_source_record_ids": len(self.ids - {""})}


def json_line(handle, value):
    handle.write(json.dumps(value, ensure_ascii=False, sort_keys=True) + "\n")


def markdown_report(report):
    pools = report["pools"]
    lines = ["# NMRShiftDB2 实验谱筛选报告", "", f"策略：`{POLICY}`", "",
             "按每张谱的索引关联来源字段，不以 Spectrum 0 或整个分子判断。",
             "experimental 表示按本次 SDF 的来源/条件字段归入实验类；不是逐篇文献或原始 FID 的独立复核。",
             "实验条件字段存在但内容为 Unreported，不代表溶剂/温度实际已知。", "",
             "|范围|分子记录|不同源 InChIKey 字段值|C 谱|H 谱|同时有 C/H 的记录|", "|---|---:|---:|---:|---:|---:|"]
    labels = {"source_all": "原始文件（所有来源）", "experimental_all_nmr": "全部实验 NMR 类型",
              "experimental_ch_raw": "实验 C/H（未经可用性检查）",
              "experimental_ch_usable": "实验 C/H 可用池", "experimental_ch_both_usable": "实验 C/H 双谱可用池"}
    for name, label in labels.items():
        p = pools[name]
        t = p["by_spectrum_type"]
        lines.append(f"|{label}|{p.get('molecule_records', 0):,}|{p['unique_source_inchi_keys']:,}|{t.get('13C', 0):,}|{t.get('1H', 0):,}|{p.get('molecules_with_both_ch', 0):,}|")
    lines += ["", "## 谱来源（计数包含空峰表和无效结构，确保源文件全部入账）", "",
              "|谱类型|experimental|computed|unknown|conflict|non_nmr|", "|---|---:|---:|---:|---:|---:|"]
    for nucleus, counts in report["classification_by_type"].items():
        lines.append(f"|{nucleus}|" + "|".join(str(counts.get(s, 0)) for s in ("experimental", "computed", "unknown", "conflict", "non_nmr")) + "|")
    lines += ["", "## 规则与边界", "",
              "- Program、NMRProgram 或其他计算条件条目：computed；计算与实验元数据同索引并存：conflict，排除。",
              "- 无计算条件，且有该谱索引的 Solvent/Temperature/Field Strength/Assignment Method 条目：experimental。",
              "- 没有足够来源依据或条件格式无法定位到谱：unknown，排除。Assignment Method 中提到预测辅助不自动判为计算谱。",
              "- experimental_all_nmr.sd 保留全部实验类型及原始问题记录，仅适合归档，不直接作为训练入口。",
              "- experimental_ch.sd 只保留实验 13C/1H、RDKit 可解析并通过 sanitize 的结构，以及至少一个当前转换器可读取的有限位移归属峰；含非有限位移的谱整张排除。",
              "- 同一谱索引的条件字段有重复值时原文全部保留在实验归档及审计文件，但该谱暂不进入可用池，避免转换器静默选取最后一个值。",
              "- experimental_ch_both.sd 进一步要求同条分子记录同时有可用实验 C/H；不保证溶剂、温度相同。",
              "- 原始分子块、原子顺序、保留的 Spectrum 内容不变；删除被排除谱及对应索引元数据，不重编号。",
              "- 派生文件统一属性头及分隔空行的格式，修复源文件少数异常换行导致 RDKit 漏读谱字段的问题；数值不变。",
              "- 部分峰记录仍可能无归属/归属无效，原文不改，由转换器跳过；详见 spectra_audit.jsonl 的 peak_quality。",
              "- 分子记录数不等于去重结构数；InChIKey 及连接层统计使用源文件字段，不重新计算。少量源标识格式不完整，malformed_source_inchi_key_records 单独计数；这些不同字段值不能全部视为已验证的不同结构。本步骤未做分子合并。",
              "- 可用池不是最终 LMDB 样本数：后续仍有结构去重、划分、谱对选择等步骤。",
              "- 未修改已有 LMDB、GNN 数据库或模型；此前训练产物需要用实验池重建后再训练。", "",
              "## 排除原因及计数", "", "```json", json.dumps(report["counters"], ensure_ascii=False, indent=2), "```", "",
              "## 来源与完整性", "", f"原文件 SHA256：`{report['source_sha256']}`", "",
              f"源文件运行前后 SHA256 一致：`{report['source_unchanged']}`；RDKit：`{report['rdkit_version']}`。", "",
              "[NMRShiftDB2 官方说明 §2.3 / §2.4](https://nmrshiftdb.nmr.uni-koeln.de/nmrshiftdbhtml/using.html)：库同时收录实验和计算谱，两类谱有不同的来源条件字段。", ""]
    return "\n".join(lines)


def run(source, output, max_records=None):
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Output directory must be new or empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    source_sha = sha256(source)
    counters, by_type, programs = Counter(), defaultdict(Counter), Counter()
    nonempty_by_type = defaultdict(Counter)
    pools = {name: Pool() for name in ["source_all", "experimental_ch_raw", *OUTPUTS]}
    artifact_names = [*OUTPUTS.values(), "spectra_audit.jsonl", "records_audit.jsonl"]
    RDLogger.DisableLog("rdApp.*")
    with ExitStack() as stack:
        sdfs = {name: stack.enter_context((output / (filename + ".partial")).open("wb")) for name, filename in OUTPUTS.items()}
        spectra_file = stack.enter_context((output / "spectra_audit.jsonl.partial").open("w", encoding="utf-8", newline="\n"))
        records_file = stack.enter_context((output / "records_audit.jsonl.partial").open("w", encoding="utf-8", newline="\n"))
        for index, line, raw in iter_records(source):
            if max_records is not None and index > max_records:
                break
            block, props, _ = parse_record(raw)
            counters["source_records"] += 1
            formatting_repair = any(
                raw[m.end():m.end() + 1] in (b"\r", b"\n") and bool(props[m[1].decode("utf-8")])
                for m in HEADER.finditer(raw)
            )
            if formatting_repair:
                counters["records_with_leading_blank_property_value_repaired"] += 1
            mol, error = validate_molecule(block)
            counters["invalid_structure_records" if error else "valid_structure_records"] += 1
            if error == "mol_parse_failed":
                counters["structure_parse_failed"] += 1
            elif error:
                counters["structure_sanitize_failed"] += 1
            spectra, exp_ch, exp_all, usable = [], [], [], []
            for field, value in props.items():
                match = SPECTRUM.fullmatch(field)
                if not match:
                    continue
                nucleus, sid = match.groups()
                status, evidence = classify(props, sid, nucleus)
                counters["source_spectra_fields"] += 1
                by_type[nucleus][status] += 1
                if value.strip():
                    nonempty_by_type[nucleus][status] += 1
                quality = signal_quality(value, mol, nucleus) if nucleus in {"13C", "1H"} else None
                row = {"source_record_index": index, "source_start_line": line,
                       "record_id": props.get("nmrshiftdb2 ID"), "inchi_key": props.get("INChI key"),
                       "source_field": field, "spectrum_id": sid, "spectrum_type": nucleus,
                       "source_class": status, "source_evidence": evidence, "structure_valid": not error,
                       "peak_quality": quality, "training_pool_eligible": False}
                spectra.append(row)
                if status == "computed":
                    for name in ("Program", "NMRProgram"):
                        if name in evidence:
                            program = evidence[name]
                            programs[(nucleus, name, program if isinstance(program, str) else json.dumps(program))] += 1
                if status == "experimental":
                    exp_all.append(row)
                    if nucleus in {"13C", "1H"}:
                        exp_ch.append(row)
                        reason = None
                        if error:
                            reason = "invalid_structure"
                        elif any(isinstance(value, list) for value in evidence.values()):
                            reason = "duplicate_condition_metadata"
                        elif quality.get("nonfinite_assignment_records", 0):
                            reason = "nonfinite_shifts"
                        elif not quality.get("raw_assignment_records", 0):
                            reason = "empty_spectrum"
                        elif not quality.get("unique_usable_peaks", 0):
                            reason = "no_usable_assigned_peak"
                        if reason:
                            counters[f"experimental_{nucleus}_excluded_{reason}"] += 1
                            row["training_exclusion_reason"] = reason
                        else:
                            row["training_pool_eligible"] = True
                            usable.append(row)
                            if quality.get("unusable_assignment_records", 0):
                                counters[f"usable_{nucleus}_spectra_with_skipped_signal_records"] += 1
                json_line(spectra_file, row)
            pools["source_all"].add(props, spectra, include_empty=True)
            if not spectra:
                counters["source_records_without_spectrum"] += 1
            pools["experimental_ch_raw"].add(props, exp_ch)
            both = usable if {row["spectrum_type"] for row in usable} == {"13C", "1H"} else []
            selected = {"experimental_all_nmr": exp_all, "experimental_ch_usable": usable,
                        "experimental_ch_both_usable": both}
            for name, rows in selected.items():
                pools[name].add(props, rows)
                if rows:
                    sdfs[name].write(filter_record(raw, {row["source_field"] for row in rows}, index, source_sha))
            json_line(records_file, {"source_record_index": index, "source_start_line": line,
                                     "record_id": props.get("nmrshiftdb2 ID"), "inchi_key": props.get("INChI key"),
                                     "structure_error": error,
                                     "property_leading_blank_repaired": bool(formatting_repair),
                                     "source_spectrum_count": len(spectra),
                                     "exported_fields": {name: [row["source_field"] for row in rows] for name, rows in selected.items()}})
            if index % 5000 == 0:
                print(f"[filter] processed {index:,} source records", flush=True)
    source_unchanged = sha256(source) == source_sha
    if not source_unchanged:
        raise RuntimeError("Source changed during filtering; outputs remain .partial")
    summaries = {name: pool.summary() for name, pool in pools.items()}
    assert sum(sum(counts.values()) for counts in by_type.values()) == counters["source_spectra_fields"]
    assert counters["source_records"] == counters["invalid_structure_records"] + counters["valid_structure_records"]
    for filename in artifact_names:
        (output / (filename + ".partial")).rename(output / filename)
    report = {"policy": POLICY, "source_sdf": str(source), "source_sha256": source_sha,
              "source_unchanged": source_unchanged, "rdkit_version": rdBase.rdkitVersion,
              "partial_scan": max_records is not None, "counters": dict(sorted(counters.items())),
              "pools": summaries, "classification_by_type": {k: dict(v) for k, v in sorted(by_type.items())},
              "nonempty_spectrum_fields_by_type": {k: dict(v) for k, v in sorted(nonempty_by_type.items())},
              "computed_programs": [{"spectrum_type": k[0], "field": k[1], "value": k[2], "count": v} for k, v in sorted(programs.items())],
              "outputs": {name: {"filename": filename, "sha256": sha256(output / filename)} for name, filename in OUTPUTS.items()}}
    (output / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (output / "report.md").write_text(markdown_report(report), encoding="utf-8")
    manifest = {"policy": POLICY, "generated_at": datetime.now(timezone.utc).isoformat(),
                "script_sha256": sha256(__file__), "source": {"path": str(source), "sha256": source_sha},
                "parameters": {"max_records": max_records},
                "artifacts": {name: sha256(output / name) for name in artifact_names + ["report.json", "report.md"]}}
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(output), "counters": report["counters"], "pools": summaries}, ensure_ascii=False, indent=2))
    return report


def assert_experimental_ch(path):
    """Fail closed at a training entry point if any C/H source is not experimental."""
    counts = Counter()
    for index, _, raw in iter_records(path):
        _, props, _ = parse_record(raw)
        for name in props:
            match = SPECTRUM.fullmatch(name)
            if not match or match[1] not in {"13C", "1H"}:
                continue
            status, _ = classify(props, match[2], match[1])
            if status != "experimental":
                raise ValueError(f"Experimental-only training refused: record {index}, {name}, status={status}. Run filter_nmrshiftdb2_experimental.py first.")
            counts[match[1]] += 1
    if not counts:
        raise ValueError("No experimental C/H spectra in training source")
    return dict(counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sdf", type=Path, default=ROOT / "data/nmrshiftdb2/nmrshiftdb2withsignals.sd")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/nmrshiftdb2_experimental_v1")
    parser.add_argument("--max-records", type=int, help="Smoke test only; output report is marked partial")
    args = parser.parse_args()
    if args.max_records is not None and args.max_records < 1:
        parser.error("--max-records must be positive")
    run(args.source_sdf, args.output_dir, args.max_records)


if __name__ == "__main__":
    main()
