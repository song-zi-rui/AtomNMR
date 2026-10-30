#!/usr/bin/env python
"""Audit ^13C decoupling evidence without using a proposed molecular structure.

The tool deliberately separates two questions which must not be allowed to
confirm one another:

1. What does the source text say about the acquisition/decoupling protocol?
2. Given that protocol, may a reported splitting be used as a *soft*
   structure-consistency factor later?

For an ordinary unannotated ``13C NMR`` entry the conservative default is
``13C{1H}``, with every other decoupling state recorded as ``unknown``.  It is
an assumption, not a discovered fact.  In particular, absence of a ``31P`` or
``19F`` annotation never means that nucleus was confirmed not to be decoupled.

Supported inputs are the raw NMRexp parquet shards and NMRShiftDB2 SDF.  The
output contains only acquisition-text evidence and is intentionally independent
of SMILES, atom assignments, AtomPeak candidates, and predicted shifts.  This
keeps protocol assessment from becoming a self-confirming structure rule.

Examples
--------
Smoke-test NMRexp (writes a resumable, small audit):

    python tools/assess_c13_decoupling.py nmrexp \
      --parquet-dir data/NMRexp/source_parquet/ml_ready \
      --output-dir data/c13_decoupling_audit/smoke --max-records 1000

Run the full NMRexp audit, retaining a compact record for every spectrum:

    python tools/assess_c13_decoupling.py nmrexp \
      --parquet-dir data/NMRexp/source_parquet/ml_ready \
      --output-dir data/c13_decoupling_audit/nmrexp_full --write-all

Audit NMRShiftDB2 metadata (which normally has no acquisition protocol text):

    python tools/assess_c13_decoupling.py nmrshiftdb2 \
      --sdf data/nmrshiftdb2/nmrshiftdb2withsignals.sd \
      --output-dir data/c13_decoupling_audit/nmrshiftdb2_full
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Iterator


SCHEMA_VERSION = "c13-decoupling-assessment-v1"
DEFAULT_HYPOTHESIS = "C13_H_decoupled_other_unknown"
MISSING = {"", "_", "-", "none", "null", "unknown", "nan", "n/a"}
SPECTRUM_13C_RE = re.compile(r"^Spectrum\s+13C\s+(\d+)$")
SPECTRUM_RE = re.compile(r"^Spectrum\s+(.+?)\s+(\d+)$")


NUCLEUS_PATTERNS = {
    "1H": r"(?:\^?\s*1\s*H|protons?|\bH\b)",
    "19F": r"(?:\^?\s*19\s*F|fluorine|\bF\b)",
    "31P": r"(?:\^?\s*31\s*P|phosphorus|\bP\b)",
}
NORM_NUCLEUS = {
    "1H": "1H", "H": "1H", "PROTON": "1H", "PROTONS": "1H",
    "19F": "19F", "F": "19F", "FLUORINE": "19F",
    "31P": "31P", "P": "31P", "PHOSPHORUS": "31P",
}
EDITED_RE = re.compile(r"\b(?:DEPT(?:[- ]?\d+)?|APT|JMOD|attached[- ]proton)\b", re.I)
CALCULATED_RE = re.compile(
    r"\b(?:HOSE|ACD(?:/LABS)?|CSEARCH|PREDICT(?:ED|ION)?|CALCULATED|SIMULATED)\b", re.I
)
GENERIC_DECOUPLING_RE = re.compile(r"\bdecoupl(?:ed|ing)\b", re.I)
GENERIC_J_RE = re.compile(r"(?:\b\d+\s*)?J\s*(?:[A-Za-z]{1,2})?\s*=\s*\d+(?:\.\d+)?", re.I)
MULTIPLICITY_RE = re.compile(
    r"\(\s*(?:br\.?\s*)?(s|d|t|q|m|dd|dt|td|dq|qd|tt|ddd|ddt|dtd|tdd)\s*(?:,|\))",
    re.I,
)
J_ASSIGNMENT_RE = re.compile(
    r"(?:(?P<order>\d+)\s*)?J\s*(?P<label>[A-Za-z]{0,2})\s*=\s*"
    r"(?P<values>\d+(?:\.\d+)?(?:\s*(?:,|and)\s*\d+(?:\.\d+)?)*)\s*(?P<unit>Hz)?",
    re.I,
)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def clean_text(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return None if text.lower() in MISSING else text


def compact_text(value: Any, limit: int = 800) -> str | None:
    text = clean_text(value)
    if text is None:
        return None
    text = re.sub(r"\s+", " ", text)
    return text if len(text) <= limit else text[: limit - 1] + "…"


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, ensure_ascii=True, separators=(",", ":"), default=str)


def stable_id(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()[:20]


def path_signature(path: Path) -> dict[str, Any]:
    path = path.expanduser().resolve()
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def parse_indexed_property(value: str | None) -> dict[str, str]:
    """Read NMRShiftDB2's ``0:value 1:value`` metadata convention."""
    value = value or ""
    matches = list(re.finditer(r"(?:^|\s)(\d+):", value))
    result: dict[str, str] = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(value)
        result[match.group(1)] = value[match.end():end].strip()
    return result


def indexed_sdf_value(properties: dict[str, str], name: str, spectrum_id: str) -> str | None:
    return parse_indexed_property(properties.get(name)).get(str(spectrum_id))


def normalize_nuclei(text: str) -> set[str]:
    """Normalize only an explicitly delimited decoupling-nucleus expression."""
    normalized: set[str] = set()
    for raw in re.findall(r"\^?\s*\d*\s*[A-Za-z]+", text.upper()):
        token = re.sub(r"\s+", "", raw).lstrip("^")
        if token in NORM_NUCLEUS:
            normalized.add(NORM_NUCLEUS[token])
    return normalized


def named_decoupling_evidence(text: str) -> tuple[set[str], set[str], list[str]]:
    """Return declared decoupled nuclei, declared coupled nuclei, and excerpts."""
    decoupled: set[str] = set()
    coupled: set[str] = set()
    excerpts: list[str] = []
    for nucleus, expression in NUCLEUS_PATTERNS.items():
        direct = re.compile(
            rf"{expression}\s*(?:-|\s)?(?:broadband\s*)?decoupl(?:ed|ing)\b", re.I
        )
        reverse = re.compile(
            rf"\bdecoupl(?:ed|ing)\s*(?:of\s*)?{expression}\b", re.I
        )
        for match in list(direct.finditer(text)) + list(reverse.finditer(text)):
            decoupled.add(nucleus)
            excerpts.append(match.group(0))
        coupled_pattern = re.compile(
            rf"{expression}\s*(?:-|\s)?coupled\b", re.I
        )
        for match in coupled_pattern.finditer(text):
            # "H decoupled" does not match this expression because the word
            # starts with "de"; the boundary after "coupled" is intentional.
            coupled.add(nucleus)
            excerpts.append(match.group(0))
    if re.search(r"\b(?:off[- ]resonance|gated[- ]decoupl(?:ed|ing)|proton[- ]coupled)\b", text, re.I):
        coupled.add("1H")
        excerpts.append("off-resonance/gated/proton-coupled")
    return decoupled, coupled, sorted(set(excerpts), key=str.lower)


def curly_decoupling_evidence(text: str) -> tuple[set[str], list[str]]:
    decoupled: set[str] = set()
    excerpts: list[str] = []
    for match in re.finditer(r"(?:\^?\s*13\s*C|\bC)\s*\{([^{}]+)\}", text, re.I):
        nuclei = normalize_nuclei(match.group(1))
        if nuclei:
            decoupled.update(nuclei)
            excerpts.append(match.group(0))
    return decoupled, excerpts


def decode_j_partner(label: str | None) -> str | None:
    compact = re.sub(r"[^A-Za-z]", "", label or "").upper()
    if compact in {"CF", "FC", "HF", "FH"}:
        return "19F"
    if compact in {"CP", "PC", "HP", "PH"}:
        return "31P"
    if compact in {"CH", "HC", "HH"}:
        return "1H"
    return None


def j_evidence(text: str) -> dict[str, Any]:
    values: list[float] = []
    partners: set[str] = set()
    units_explicit = False
    labels: list[str] = []
    for match in J_ASSIGNMENT_RE.finditer(text):
        raw_values = re.findall(r"\d+(?:\.\d+)?", match.group("values"))
        values.extend(float(value) for value in raw_values)
        partner = decode_j_partner(match.group("label"))
        if partner:
            partners.add(partner)
        if match.group("label"):
            labels.append(match.group("label").upper())
        units_explicit = units_explicit or bool(match.group("unit"))
    return {
        "reported": bool(values) or bool(GENERIC_J_RE.search(text)),
        "values_hz": sorted(set(values)),
        "unit_explicit": units_explicit,
        "partner_nuclei": sorted(partners),
        "raw_labels": sorted(set(labels)),
        "multiplicity_labels": sorted(set(match.group(1).lower() for match in MULTIPLICITY_RE.finditer(text))),
    }


def evidence_text(record: dict[str, Any]) -> str:
    return "\n".join(
        value for value in (
            clean_text(record.get("experiment_type_raw")),
            clean_text(record.get("note_raw")),
            clean_text(record.get("shift_text_raw")),
            clean_text(record.get("program_raw")),
        ) if value
    )


def hypothesis(identifier: str, role: str, description: str) -> dict[str, str]:
    return {"id": identifier, "role": role, "description": description}


def assess_record(record: dict[str, Any]) -> dict[str, Any]:
    """Assess an acquisition record without reading structure/candidate fields.

    This function is deliberately pure so it can be unit-tested and reused by
    later data materialization code.  A caller may supply extra keys (including
    SMILES), but they are not read here.
    """
    if str(record.get("observed_nucleus") or "13C") != "13C":
        raise ValueError("assess_record accepts only observed_nucleus='13C'")

    text = evidence_text(record)
    brace_decoupled, brace_excerpts = curly_decoupling_evidence(text)
    named_decoupled, named_coupled, named_excerpts = named_decoupling_evidence(text)
    declared_decoupled = brace_decoupled | named_decoupled
    declared_coupled = named_coupled - declared_decoupled
    j = j_evidence(clean_text(record.get("shift_text_raw")) or "")
    edited_match = EDITED_RE.search(text)
    calculated_match = CALCULATED_RE.search(text)
    generic_decoupling = bool(GENERIC_DECOUPLING_RE.search(text))
    protocol_excerpts = sorted(set(brace_excerpts + named_excerpts), key=str.lower)
    review_reasons: list[str] = []

    # The default is a prior needed to choose downstream rules, not a claim
    # recovered from the source.  Other nuclei intentionally stay unknown.
    assumed_decoupled = set() if "1H" in declared_coupled else {"1H"}
    other_state = {
        nucleus: "declared_decoupled" if nucleus in declared_decoupled else "unknown"
        for nucleus in ("19F", "31P")
    }

    conflict = bool(declared_decoupled & declared_coupled)
    if calculated_match:
        primary = "C13_nonexperimental_prediction_or_calculation"
        decision = "exclude_nonexperimental_prediction_or_calculation"
        protocol_confidence = "not_applicable"
        carbon_type_policy = "disabled"
        splitting_policy = "disabled"
        review_reasons.append("source metadata identifies a prediction/calculation program, not an experimental acquisition")
        hypotheses = [hypothesis(
            primary,
            "excluded",
            "Do not use this record for experimental splitting, decoupling, or spin-constraint learning.",
        )]
    elif edited_match:
        primary = "C13_edited_CH_type_experiment"
        decision = "edited_experiment_declared"
        protocol_confidence = "declared"
        carbon_type_policy = "eligible_only_if_explicit_edited_type_is_extracted"
        splitting_policy = "disabled"
        hypotheses = [hypothesis(primary, "declared", "DEPT/APT/JMOD-like experiment; do not treat its labels as line splitting.")]
    elif conflict:
        primary = "C13_protocol_conflict"
        decision = "review_conflicting_protocol_metadata"
        protocol_confidence = "conflicting"
        carbon_type_policy = "disabled"
        splitting_policy = "disabled"
        review_reasons.append("same source text declares a nucleus both coupled and decoupled")
        hypotheses = [
            hypothesis(DEFAULT_HYPOTHESIS, "parallel", "Default only; source metadata is internally inconsistent."),
            hypothesis("C13_proton_coupled", "parallel", "Alternative required because protocol text conflicts."),
        ]
    elif "1H" in declared_coupled:
        primary = "C13_proton_coupled"
        decision = "proton_coupled_declared"
        protocol_confidence = "declared"
        carbon_type_policy = "disabled"
        splitting_policy = "reported_CH_only_soft"
        hypotheses = [hypothesis(primary, "declared", "Source explicitly indicates a proton-coupled/off-resonance/gated carbon experiment.")]
    elif "1H" in declared_decoupled:
        primary = DEFAULT_HYPOTHESIS
        decision = "proton_decoupling_declared"
        protocol_confidence = "declared"
        carbon_type_policy = "disabled"
        splitting_policy = "reported_heteronuclear_only_soft"
        hypotheses = [hypothesis(primary, "declared", "Source explicitly declares 1H decoupling; other nuclei remain source-specific.")]
    else:
        primary = DEFAULT_HYPOTHESIS
        decision = "default_H_decoupled_assumption"
        protocol_confidence = "assumed"
        carbon_type_policy = "disabled"
        splitting_policy = "reported_heteronuclear_only_soft"
        hypotheses = [hypothesis(primary, "default", "Conservative prior for bare '13C NMR'; it is not verified metadata.")]
        review_reasons.append("no explicit 1H decoupling or pulse-program metadata")

    if generic_decoupling and not declared_decoupled and not declared_coupled:
        review_reasons.append("source says 'decoupled' but does not name the nucleus")

    partner_set = set(j["partner_nuclei"])
    if j["reported"] and not calculated_match:
        if not j["unit_explicit"]:
            review_reasons.append("reported J-like value has no explicit Hz unit")
        if "1H" in partner_set:
            if "1H" in declared_decoupled:
                decision = "review_reported_CH_coupling_conflicts_with_H_decoupling"
                protocol_confidence = "conflicting"
                review_reasons.append("reported C-H coupling conflicts with declared 1H decoupling")
            elif "1H" not in declared_coupled:
                decision = "parallel_reanalysis_required"
                review_reasons.append("reported C-H coupling is incompatible with an unverified default 13C{1H} assumption")
                hypotheses.append(hypothesis("C13_proton_coupled", "parallel", "Reported C-H J needs a proton-coupled alternative."))
            splitting_policy = "reported_CH_only_soft"
        elif partner_set & {"19F", "31P"}:
            # Observing a labelled J is positive evidence that this particular
            # coupling survived into the source description.  It is still not
            # licence to infer anything from *unreported* P/F splittings.
            observed = sorted(partner_set & {"19F", "31P"})
            explicit_conflict = [nucleus for nucleus in observed if nucleus in declared_decoupled]
            if explicit_conflict:
                decision = "review_reported_heteronuclear_coupling_conflicts_with_decoupling"
                protocol_confidence = "conflicting"
                review_reasons.append("reported heteronuclear J conflicts with a named decoupling declaration: " + ", ".join(explicit_conflict))
            else:
                if decision == "default_H_decoupled_assumption":
                    decision = "default_H_decoupled_with_reported_heteronuclear_evidence"
                hypotheses.append(hypothesis(
                    "C13_H_decoupled_heteronuclear_visible",
                    "observed_evidence",
                    "Use only explicitly reported P/F coupling as a soft candidate-specific factor.",
                ))
            splitting_policy = "reported_heteronuclear_soft"
        elif not partner_set:
            # A bare 'd, J=...' can be an H, F, P, or other coupling.  The
            # correct action is sensitivity analysis, not an invented partner.
            decision = "parallel_reanalysis_required"
            review_reasons.append("reported J has no coupled-nucleus identity")
            hypotheses.extend([
                hypothesis("C13_proton_coupled", "parallel", "One possible origin of an unlabelled carbon J."),
                hypothesis("C13_H_decoupled_heteronuclear_visible", "parallel", "Alternative heteronuclear origin; partner remains unknown."),
            ])
            splitting_policy = "unresolved_J_soft_only"

    if not j["reported"] and not declared_decoupled and not declared_coupled and not edited_match:
        # A centre-shift list cannot establish whether any coupling was absent;
        # do not activate a spin constraint merely because the default is common.
        splitting_policy = "disabled"

    requires_parallel = decision.startswith("parallel_") or decision.startswith("review_")
    input_capability = record.get("input_capability") or (
        "reported_peak_text" if clean_text(record.get("shift_text_raw")) else "shift_list_only"
    )
    result = {
        "schema_version": SCHEMA_VERSION,
        "record_id": record["record_id"],
        "source": record.get("source"),
        "source_locator": record.get("source_locator"),
        "observed_nucleus": "13C",
        "input_capability": input_capability,
        "source_metadata": {
            "experiment_type_raw": compact_text(record.get("experiment_type_raw")),
            "note_raw": compact_text(record.get("note_raw")),
            "program_raw": compact_text(record.get("program_raw")),
            "shift_text_raw": compact_text(record.get("shift_text_raw")),
            # This provenance may help a later human review locate companion
            # 19F/31P spectra.  It is deliberately not used to infer a C13
            # decoupling mode in this acquisition-only classifier.
            "related_source_spectrum_types": record.get("related_source_spectrum_types") or [],
        },
        "protocol": {
            "declared_decoupled_nuclei": sorted(declared_decoupled),
            "declared_coupled_nuclei": sorted(declared_coupled),
            "assumed_decoupled_nuclei": sorted(assumed_decoupled),
            "other_decoupling_state": other_state,
            "protocol_excerpts": protocol_excerpts,
            "primary_hypothesis": primary,
            "hypotheses": {item["id"]: item for item in {item["id"]: item for item in hypotheses}.values()},
            "confidence": protocol_confidence,
        },
        "splitting_evidence": {
            "reported_J": j["reported"],
            "j_values_hz": j["values_hz"],
            "j_unit_explicit": j["unit_explicit"],
            "partner_nuclei": j["partner_nuclei"],
            "raw_labels": j["raw_labels"],
            "multiplicity_labels": j["multiplicity_labels"],
        },
        "policy": {
            "decision": decision,
            "requires_parallel_reanalysis": requires_parallel,
            "carbon_hydrogen_class_constraint": carbon_type_policy,
            "splitting_constraint": splitting_policy,
            "requires_unique_raw_peak_alignment": bool(j["reported"]),
            "eligible_for_direct_training_constraint": False,
            "missing_other_decoupling_means": "unknown_not_not_decoupled",
            "candidate_structure_used_in_assessment": False,
        },
        "review_reasons": sorted(set(review_reasons)),
    }
    return result


def iter_nmrexp(parquet_dir: Path, max_records: int | None) -> Iterator[dict[str, Any]]:
    try:
        import pyarrow.parquet as pq
    except ImportError as exc:  # pragma: no cover - environment specific
        raise RuntimeError("NMRexp scan requires pyarrow. Install the project runtime with pyarrow available.") from exc
    columns = ["Filename", "NMR_type", "NMR_frequency", "NMR_solvent", "NMR_shift_text", "NMR_note", "NMR_processed"]
    count = 0
    for shard in sorted(parquet_dir.glob("*.parquet")):
        row_offset = 0
        reader = pq.ParquetFile(shard)
        for batch in reader.iter_batches(batch_size=65536, columns=columns):
            for local_index, row in enumerate(batch.to_pylist()):
                if not str(row.get("NMR_type") or "").strip().lower().startswith("13c"):
                    continue
                source_row = row_offset + local_index
                locator = f"{shard.name}:{source_row}"
                yield {
                    "record_id": "nmrexp:" + stable_id({"shard": shard.name, "row": source_row}),
                    "source": "NMRexp",
                    "source_locator": locator,
                    "observed_nucleus": "13C",
                    "experiment_type_raw": row.get("NMR_type"),
                    "note_raw": row.get("NMR_note"),
                    "program_raw": None,
                    "shift_text_raw": row.get("NMR_shift_text"),
                    "input_capability": "reported_peak_text" if clean_text(row.get("NMR_shift_text")) else "shift_list_only",
                }
                count += 1
                if max_records is not None and count >= max_records:
                    return
            row_offset += batch.num_rows


def iter_nmrshiftdb2(sdf_path: Path, max_records: int | None) -> Iterator[dict[str, Any]]:
    try:
        from rdkit import Chem, RDLogger
    except ImportError as exc:  # pragma: no cover - environment specific
        raise RuntimeError("NMRShiftDB2 scan requires RDKit.") from exc
    RDLogger.DisableLog("rdApp.*")
    count = 0
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    for mol_index, mol in enumerate(supplier):
        if mol is None:
            continue
        properties = {name: mol.GetProp(name) for name in mol.GetPropNames()}
        db_id = properties.get("nmrshiftdb2 ID") or str(mol_index)
        related_types = sorted({match.group(1) for name in properties if (match := SPECTRUM_RE.match(name))})
        for field in sorted(properties):
            match = SPECTRUM_13C_RE.match(field)
            if not match:
                continue
            spectrum_id = match.group(1)
            yield {
                "record_id": "nmrshiftdb2:" + stable_id({"id": db_id, "spectrum_id": spectrum_id}),
                "source": "NMRShiftDB2",
                "source_locator": f"nmrshiftdb2_id={db_id};spectrum_id={spectrum_id}",
                "observed_nucleus": "13C",
                "experiment_type_raw": "13C NMR",
                "note_raw": indexed_sdf_value(properties, "NMR Note", spectrum_id),
                "program_raw": " | ".join(
                    value for value in (
                        indexed_sdf_value(properties, "Program", spectrum_id),
                        indexed_sdf_value(properties, "NMRMethod", spectrum_id),
                        indexed_sdf_value(properties, "NMRProgram", spectrum_id),
                        indexed_sdf_value(properties, "BasisSet", spectrum_id),
                    ) if value
                ) or None,
                "shift_text_raw": properties.get(field),
                "input_capability": "assigned_shift_list",
                "related_source_spectrum_types": related_types,
            }
            count += 1
            if max_records is not None and count >= max_records:
                return


def source_signature(args: argparse.Namespace) -> dict[str, Any]:
    if args.source == "nmrexp":
        directory = args.parquet_dir.expanduser().resolve()
        files = sorted(directory.glob("*.parquet"))
        if not files:
            raise FileNotFoundError(f"No parquet shards found under {directory}")
        return {"source": "NMRexp", "files": [path_signature(path) for path in files]}
    sdf = args.sdf.expanduser().resolve()
    if not sdf.is_file():
        raise FileNotFoundError(sdf)
    return {"source": "NMRShiftDB2", "sdf": path_signature(sdf)}


def assessment_iter(args: argparse.Namespace) -> Iterable[dict[str, Any]]:
    if args.source == "nmrexp":
        return iter_nmrexp(args.parquet_dir.expanduser().resolve(), args.max_records)
    return iter_nmrshiftdb2(args.sdf.expanduser().resolve(), args.max_records)


def interesting(result: dict[str, Any]) -> bool:
    return (
        result["policy"]["decision"] != "default_H_decoupled_assumption"
        or result["splitting_evidence"]["reported_J"]
        or bool(result["protocol"]["declared_decoupled_nuclei"])
        or bool(result["protocol"]["declared_coupled_nuclei"])
    )


def review_row(result: dict[str, Any]) -> dict[str, Any]:
    metadata = result["source_metadata"]
    return {
        "record_id": result["record_id"],
        "source": result.get("source"),
        "source_locator": result.get("source_locator"),
        "decision": result["policy"]["decision"],
        "primary_hypothesis": result["protocol"]["primary_hypothesis"],
        "reported_J": result["splitting_evidence"]["reported_J"],
        "partner_nuclei": ";".join(result["splitting_evidence"]["partner_nuclei"]),
        "reasons": " | ".join(result["review_reasons"]),
        "nmr_type": metadata["experiment_type_raw"] or "",
        "nmr_note": metadata["note_raw"] or "",
        "nmr_shift_text": metadata["shift_text_raw"] or "",
    }


REVIEW_FIELDS = [
    "record_id", "source", "source_locator", "decision", "primary_hypothesis", "reported_J",
    "partner_nuclei", "reasons", "nmr_type", "nmr_note", "nmr_shift_text",
]


def ensure_run_state(output_dir: Path, args: argparse.Namespace, signature: dict[str, Any]) -> tuple[dict[str, Any], Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    state_path = output_dir / "run_state.json"
    config = {
        "schema_version": SCHEMA_VERSION,
        "source": args.source,
        "max_records": args.max_records,
        "write_all": bool(args.write_all),
        "source_signature": signature,
    }
    if state_path.exists():
        if not args.resume:
            raise RuntimeError(f"{state_path} already exists; use --resume or choose a new --output-dir")
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state.get("config") != config:
            raise RuntimeError("Saved run configuration or input signature differs. Use a fresh output directory; do not mix audits.")
        if state.get("completed"):
            raise RuntimeError("This audit is already complete. Use its outputs or choose a fresh output directory.")
        return state, state_path
    if args.resume:
        raise RuntimeError(f"--resume was requested but {state_path} does not exist")
    occupied = [name for name in ("assessments.jsonl", "overrides.jsonl", "review_queue.csv") if (output_dir / name).exists()]
    if occupied:
        raise RuntimeError("Output directory contains audit files but no compatible run state: " + ", ".join(occupied))
    return {
        "schema_version": SCHEMA_VERSION,
        "config": config,
        "started_at": utc_now(),
        "updated_at": utc_now(),
        "completed": False,
        "processed_records": 0,
        "counters": {},
    }, state_path


def counter_add(counters: Counter, result: dict[str, Any]) -> None:
    counters["records"] += 1
    counters["decision:" + result["policy"]["decision"]] += 1
    counters["protocol_confidence:" + result["protocol"]["confidence"]] += 1
    counters["reported_J"] += int(result["splitting_evidence"]["reported_J"])
    counters["requires_parallel_reanalysis"] += int(result["policy"]["requires_parallel_reanalysis"])
    for nucleus in result["protocol"]["declared_decoupled_nuclei"]:
        counters["declared_decoupled:" + nucleus] += 1
    for nucleus in result["splitting_evidence"]["partner_nuclei"]:
        counters["reported_J_partner:" + nucleus] += 1


def write_report(output_dir: Path, state: dict[str, Any]) -> None:
    counters = state.get("counters", {})
    report = {
        "schema_version": SCHEMA_VERSION,
        "source": state["config"]["source"],
        "generated_at": utc_now(),
        "completed": bool(state.get("completed")),
        "processed_records": state.get("processed_records", 0),
        "counters": counters,
        "default_policy": {
            "primary_hypothesis": DEFAULT_HYPOTHESIS,
            "meaning": "Assume 1H decoupling only when no direct protocol metadata exists; all other decoupling states remain unknown.",
            "not_a_claim": "Missing 31P/19F text is not evidence that those nuclei were not decoupled.",
        },
        "use_in_atompeak": {
            "safe_default": "Use shifts only; disable carbon splitting factors for default-assumption records with no reported J.",
            "reported_heteronuclear_J": "May be used only as a soft candidate-specific factor after unique raw-to-local peak alignment, explicit Hz-unit validation, and source-provenance validation.",
            "parallel_reanalysis": "Score all declared/plausible protocol hypotheses separately, then report candidate stability. Do not select a protocol because a candidate SMILES benefits from it.",
        },
    }
    write_json(output_dir / "decoupling_report.json", report)
    lines = [
        "# ^13C 去耦信息审计报告", "",
        f"- 数据源：`{report['source']}`", f"- 已处理 ^13C 谱：{report['processed_records']:,}",
        "- 默认先验：未显式注明时按 `^13C{^1H}` 处理；^19F/^31P 等其余去耦状态为 `unknown`，不是 `not-decoupled`。", "",
        "## 统计", "",
    ]
    for key, value in sorted(counters.items()):
        lines.append(f"- `{key}`：{value:,}")
    lines.extend([
        "", "## 下游使用规则", "",
        "- 只有中心位移列表且无明确 J/实验模式：只用于位移匹配，不启用碳裂分约束。",
        "- 原文明确报告 C-F/C-P 等 J：先完成唯一的原文峰—本地峰对齐、Hz 单位和来源核验，才可作为软的候选结构一致性证据。",
        "- `parallel_reanalysis_required`：在多个实验假设下并行评分，检查候选结构和归属是否稳定；不能让候选结构反过来决定实验模式。",
    ])
    (output_dir / "decoupling_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> int:
    output_dir = args.output_dir.expanduser().resolve()
    signature = source_signature(args)
    state, state_path = ensure_run_state(output_dir, args, signature)
    completed_before = int(state.get("processed_records", 0))
    counters = Counter({key: int(value) for key, value in state.get("counters", {}).items()})
    state["updated_at"] = utc_now()
    write_json(state_path, state)

    all_path = output_dir / "assessments.jsonl"
    overrides_path = output_dir / "overrides.jsonl"
    review_path = output_dir / "review_queue.csv"
    all_handle = all_path.open("a", encoding="utf-8", newline="\n") if args.write_all else None
    overrides_handle = overrides_path.open("a", encoding="utf-8", newline="\n")
    review_exists = review_path.exists() and review_path.stat().st_size > 0
    review_handle = review_path.open("a", encoding="utf-8-sig" if not review_exists else "utf-8", newline="")
    writer = csv.DictWriter(review_handle, fieldnames=REVIEW_FIELDS)
    if not review_exists:
        writer.writeheader()

    try:
        for source_index, record in enumerate(assessment_iter(args), start=1):
            if source_index <= completed_before:
                continue
            result = assess_record(record)
            payload = json.dumps(result, ensure_ascii=False, separators=(",", ":")) + "\n"
            if all_handle is not None:
                all_handle.write(payload)
            if interesting(result):
                overrides_handle.write(payload)
            if result["policy"]["requires_parallel_reanalysis"]:
                writer.writerow(review_row(result))
            counter_add(counters, result)
            state["processed_records"] = source_index
            if source_index % args.checkpoint_every == 0:
                all_handle and all_handle.flush()
                overrides_handle.flush()
                review_handle.flush()
                state["counters"] = dict(sorted(counters.items()))
                state["updated_at"] = utc_now()
                write_json(state_path, state)
                print(f"[audit] processed {source_index:,} ^13C records", flush=True)
    except KeyboardInterrupt:
        print("[audit] interrupted; progress was saved and can be continued with --resume", file=sys.stderr)
        state["completed"] = False
        raise
    finally:
        if all_handle is not None:
            all_handle.close()
        overrides_handle.close()
        review_handle.close()

    state["counters"] = dict(sorted(counters.items()))
    state["updated_at"] = utc_now()
    state["completed"] = True
    state["finished_at"] = utc_now()
    write_json(state_path, state)
    write_report(output_dir, state)
    print(json.dumps({
        "output_dir": str(output_dir), "processed_records": state["processed_records"],
        "reported_J": counters["reported_J"],
        "parallel_reanalysis": counters["requires_parallel_reanalysis"],
        "report": str(output_dir / "decoupling_report.md"),
    }, ensure_ascii=False, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest="source", required=True)
    for source in ("nmrexp", "nmrshiftdb2"):
        child = subparsers.add_parser(source)
        child.add_argument("--output-dir", type=Path, required=True)
        child.add_argument("--max-records", type=int, help="Maximum ^13C records; useful for a smoke test.")
        child.add_argument("--write-all", action="store_true", help="Write every assessed spectrum; otherwise retain only protocol/J exceptions plus summary counts.")
        child.add_argument("--resume", action="store_true", help="Continue an interrupted run with identical input and options.")
        child.add_argument("--checkpoint-every", type=int, default=5000)
    subparsers.choices["nmrexp"].add_argument(
        "--parquet-dir", type=Path, default=Path("data/NMRexp/source_parquet/ml_ready"),
    )
    subparsers.choices["nmrshiftdb2"].add_argument(
        "--sdf", type=Path, default=Path("data/nmrshiftdb2/nmrshiftdb2withsignals.sd"),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.max_records is not None and args.max_records <= 0:
        raise SystemExit("--max-records must be positive")
    if args.checkpoint_every <= 0:
        raise SystemExit("--checkpoint-every must be positive")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
