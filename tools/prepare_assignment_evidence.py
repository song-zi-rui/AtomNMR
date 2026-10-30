#!/usr/bin/env python
"""Add auditable spectral evidence to existing suggestions; never read gold labels."""
import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import types

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
# Same lightweight import boundary as tools/nmr_assignment_assistant.py.
if "nmrpeak" not in sys.modules:
    package = types.ModuleType("nmrpeak")
    package.__path__ = [str(Path(__file__).resolve().parents[1] / "nmrpeak")]
    package.__package__ = "nmrpeak"
    sys.modules["nmrpeak"] = package
from annotation_tool_v2.spectral_evidence import build_evidence, present, numbers, rerank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--training-report", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.report.exists():
        parser.error("output/report already exists; choose new names (no silent overwrite)")
    vocab = json.loads(args.training_report.read_text(encoding="utf-8"))["solvent_vocabulary"]
    counts, solvents, scopes = Counter(), Counter(), Counter()
    details = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_name(args.output.name + ".partial")
    with args.input.open(encoding="utf-8") as source, partial.open("x", encoding="utf-8", newline="\n") as target:
        for line in source:
            if not line.strip():
                continue
            sample = json.loads(line)
            if "spectral_evidence" in sample:
                raise ValueError("input already enriched; use the original v2 file")
            sample["spectral_evidence"] = evidence = build_evidence(sample, vocab)
            rerank(sample)
            # Keep nucleus-specific conditions; never fabricate mixture composition.
            sample.setdefault("conditions", {})["by_nucleus"] = evidence["conditions_by_nucleus"]
            target.write(json.dumps(sample, ensure_ascii=False, allow_nan=False) + "\n")
            target.flush()
            counts["records"] += 1
            counts["source_nmrexp"] += sample.get("source", {}).get("database") == "NMRexp"
            counts["verified_provenance"] += sample.get("provenance", {}).get("status") == "verified"
            cond = evidence["conditions_by_nucleus"]
            counts["same_reported_h_c_solvent"] += present(cond["1H"]["solvent_raw"]) and cond["1H"]["solvent_raw"] == cond["13C"]["solvent_raw"]
            counts["both_frequencies_reported"] += all(present(c["frequency_raw"]) for c in cond.values())
            counts["meaningful_notes"] += any(present(c["note"]) for c in cond.values())
            counts["other_solvent_embedding"] += evidence["solvent_model_status"] == "other_embedding"
            original_c = ast.literal_eval(sample.get("provenance", {}).get("c_source", {}).get("NMR_processed") or "[]")
            original_rich = [r for r in original_c if any(v is not None for v in r[1:])]
            counts["source_c_multiplicity_or_j_entries"] += len(original_rich)
            counts["source_c_multiplicity_or_j_records"] += bool(original_rich)
            not_unique = sum(v["alignment"] != "unique" for v in evidence["carbon_peaks"].values())
            counts["c_metadata_alignment_unresolved_peaks"] += not_unique
            counts["c_metadata_alignment_unresolved_records"] += bool(not_unique)
            scopes.update(evidence["reference_scope_atom_counts"])
            solvents[sample["conditions"].get("solvent_normalized", "unknown")] += 1
            for nucleus in ("13C", "1H"):
                peaks = sample["spectra"][nucleus]["peaks"]
                counts[nucleus + "_peaks"] += len(peaks)
                fields = Counter()
                for p in peaks:
                    meta = evidence["carbon_peaks"][p["peak_id"]] if nucleus == "13C" else p
                    fields["multiplicity_peaks"] += present(meta.get("multiplicity"))
                    fields["j_peaks"] += bool(numbers(meta.get("j_values")))
                    fields["integral_peaks"] += present(meta.get("nH"))
                    fields["range_peaks"] += p.get("range_min") is not None and p.get("range_max") is not None and p["range_min"] != p["range_max"]
                    if nucleus == "13C":
                        fields["raw_confirmed_j_peaks"] += meta["model_eligible"]
                for key, count in fields.items():
                    counts[nucleus + "_" + key] += count
                    counts[nucleus + "_" + key.replace("_peaks", "_records")] += bool(count)
            counts["reranked_records"] += bool(evidence.get("ranking_rules"))
            counts["changed_records"] += bool(evidence.get("changed_peak_ids"))
            counts["changed_peaks"] += len(evidence.get("changed_peak_ids", []))
            if evidence.get("ranking_rules"):
                details.append({"sample_id": sample["sample_id"], "rules": evidence["ranking_rules"],
                                "changed_peak_ids": evidence.get("changed_peak_ids", [])})
    partial.rename(args.output)
    report = {"schema_version": "spectral-evidence-audit-v1", "input": str(args.input),
              "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
              "training_report_sha256": hashlib.sha256(args.training_report.read_bytes()).hexdigest(),
              "output_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
              "counts": dict(sorted(counts.items())), "solvents": dict(solvents),
              "reference_scope_atoms": dict(scopes), "reranking_details": details,
              "performance_claim": "No gold evaluation or retraining; changes are review hypotheses, not accuracy gains."}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
    print(json.dumps({k: v for k, v in report.items() if k != "reranking_details"}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
