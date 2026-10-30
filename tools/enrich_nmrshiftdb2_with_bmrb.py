#!/usr/bin/env python
"""Enrich every NMRShiftDB2 1D 1H/13C spectrum from BMRB Metabolomics.

The pipeline is deliberately conservative and restartable.  The source SDF is
never modified.  All candidate decisions and transferred fields are written to
an SQLite audit database before model-ready JSONL/LMDB files are materialised.

Typical usage::

    python tools/enrich_nmrshiftdb2_with_bmrb.py fetch
    python tools/enrich_nmrshiftdb2_with_bmrb.py index --resume
    python tools/enrich_nmrshiftdb2_with_bmrb.py match --resume
    python tools/enrich_nmrshiftdb2_with_bmrb.py export
    python tools/enrich_nmrshiftdb2_with_bmrb.py report

Use ``all`` to run every stage.  Network access is needed only by ``fetch``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import pickle
import re
import shutil
import sqlite3
import statistics
import sys
import time
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import lmdb
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE_SDF = ROOT / "data" / "nmrshiftdb2" / "nmrshiftdb2withsignals.sd"
DEFAULT_OUTPUT_DIR = ROOT / "data" / "nmrshiftdb2_bmrb_enriched"
DEFAULT_BMRB_DIR = DEFAULT_OUTPUT_DIR / "bmrb_snapshot"
BMRB_BASE_URL = "https://bmrb.io/ftp/pub/bmrb/relational_tables/metabolomics"
REQUIRED_TABLES = (
    "Chem_comp.csv",
    "Chem_comp_SMILES.csv",
    "Chem_comp_descriptor.csv",
    "Chem_comp_db_link.csv",
    "Experiment.csv",
    "NMR_spectrometer.csv",
    "Peak.csv",
    "Peak_char.csv",
    "Peak_coupling.csv",
    "Sample.csv",
    "Sample_component.csv",
    "Sample_condition_variable.csv",
    "Spectral_peak_list.csv",
)
SPECTRUM_RE = re.compile(r"^Spectrum\s+(.+?)\s+(\d+)$")
INDEX_RE = re.compile(r"(?:^|\s)(\d+):")
NUMBER_RE = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)")
SCHEMA_VERSION = "nmrshiftdb2_bmrb_enriched_v1_atom_identity_v1"


@dataclass(frozen=True)
class MatchPolicy:
    name: str = "strict-v1"
    carbon_tolerance: float = 0.25
    proton_tolerance: float = 0.03
    min_coverage: float = 0.95
    min_count_ratio: float = 0.90
    carbon_median_max: float = 0.05
    carbon_p95_max: float = 0.15
    proton_median_max: float = 0.01
    proton_p95_max: float = 0.03
    min_margin: float = 0.05
    review_coverage: float = 0.90
    review_count_ratio: float = 0.80
    temperature_tolerance_k: float = 3.0
    field_tolerance_mhz: float = 5.0


POLICIES = {"strict-v1": MatchPolicy()}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def log(message: str) -> None:
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def json_dumps(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def clean_text(value: Any) -> str | None:
    if value is None:
        return None
    text = re.sub(r"\s+", " ", str(value).strip())
    return None if not text or text.lower() in {".", "?", "na", "n/a", "none", "null", "unknown", "unreported"} else text


def safe_float(value: Any) -> float | None:
    text = clean_text(value)
    if text is None:
        return None
    match = NUMBER_RE.search(text.replace(",", ""))
    if not match:
        return None
    try:
        result = float(match.group())
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def normalise_formula(value: Any) -> str | None:
    text = clean_text(value)
    if text is None:
        return None
    pieces = re.findall(r"([A-Z][a-z]?)(\d*)", text.replace(" ", ""))
    if not pieces:
        return re.sub(r"\s+", "", text)
    counts: dict[str, int] = defaultdict(int)
    for element, count in pieces:
        counts[element] += int(count or 1)
    order = []
    if "C" in counts:
        order.append("C")
    if "H" in counts:
        order.append("H")
    order.extend(sorted(key for key in counts if key not in {"C", "H"}))
    return "".join(element + (str(counts[element]) if counts[element] != 1 else "") for element in order)


SOLVENT_ALIASES = {
    "cdcl3": "chloroform_d",
    "chloroform-d": "chloroform_d",
    "chloroform-d1": "chloroform_d",
    "chloroform d": "chloroform_d",
    "dmso-d6": "dmso_d6",
    "dmsod6": "dmso_d6",
    "dimethyl sulfoxide-d6": "dmso_d6",
    "d2o": "water_d2",
    "water-d2": "water_d2",
    "water d2": "water_d2",
    "deuterium oxide": "water_d2",
    "cd3od": "methanol_d4",
    "methanol-d4": "methanol_d4",
    "acetone-d6": "acetone_d6",
    "c6d6": "benzene_d6",
    "benzene-d6": "benzene_d6",
    "cd3cn": "acetonitrile_d3",
    "acetonitrile-d3": "acetonitrile_d3",
    "pyridine-d5": "pyridine_d5",
    "thf-d8": "tetrahydrofuran_d8",
    "tetrahydrofuran-d8": "tetrahydrofuran_d8",
}


def normalise_solvent(value: Any) -> str | None:
    text = clean_text(value)
    if text is None:
        return None
    normal = text.lower().replace("₂", "2").replace("₃", "3").replace("₄", "4").replace("₅", "5").replace("₆", "6").replace("₈", "8")
    compact = re.sub(r"[^a-z0-9]+", "", normal)
    detected = []
    signatures = (
        ("chloroform_d", ("cdcl3", "chloroformd1", "deuterochloroform")),
        ("dmso_d6", ("dmsod6", "dimethylsulfoxided6")),
        ("water_d2", ("d2o", "deuteriumoxide", "heavywater")),
        ("methanol_d4", ("cd3od", "methanold4")),
        ("acetone_d6", ("acetoned6",)),
        ("benzene_d6", ("c6d6", "benzened6")),
        ("acetonitrile_d3", ("cd3cn", "acetonitriled3")),
        ("pyridine_d5", ("pyridined5",)),
        ("tetrahydrofuran_d8", ("thfd8", "tetrahydrofurand8")),
    )
    for canonical, aliases in signatures:
        if any(alias in compact for alias in aliases):
            detected.append(canonical)
    if detected:
        return "+".join(sorted(dict.fromkeys(detected)))
    normal = re.sub(r"[\[\](){}]", "", normal)
    normal = re.sub(r"\s+", " ", normal).strip()
    if normal in SOLVENT_ALIASES:
        return SOLVENT_ALIASES[normal]
    # BMRB solvent systems are sometimes percentages or mixtures.  Preserve a
    # deterministic canonical representation without pretending they are pure.
    parts = [part.strip() for part in re.split(r"[+/;,]", normal) if part.strip()]
    mapped = [SOLVENT_ALIASES.get(part, re.sub(r"[^a-z0-9]+", "_", part).strip("_")) for part in parts]
    return "+".join(sorted(dict.fromkeys(filter(None, mapped)))) or None


def parse_indexed_property(value: str | None) -> dict[str, str]:
    matches = list(INDEX_RE.finditer(value or ""))
    result = {}
    for position, match in enumerate(matches):
        end = matches[position + 1].start() if position + 1 < len(matches) else len(value or "")
        result[match.group(1)] = (value or "")[match.end() : end].strip()
    return result


def indexed_metadata(mol: Chem.Mol, name: str, spectrum_id: str) -> str | None:
    if not mol.HasProp(name):
        return None
    return clean_text(parse_indexed_property(mol.GetProp(name)).get(str(spectrum_id)))


def parse_local_signal(raw: str, nucleus: str) -> dict[str, Any] | None:
    parts = raw.split(";")
    if not parts:
        return None
    shift = safe_float(parts[0])
    if shift is None:
        return None
    extra = parts[1].strip() if len(parts) > 1 else ""
    atoms: list[int] = []
    if len(parts) > 2:
        for token in re.split(r"[, _]+", parts[-1].strip()):
            try:
                atoms.append(int(token))
            except ValueError:
                pass
    multiplicity = None
    if nucleus == "1H":
        codes = re.findall(r"[A-Za-z]+", extra)
        multiplicity = codes[-1].lower() if codes else None
    return {
        "shift_ppm": shift,
        "shift_text": parts[0].strip(),
        "multiplicity_raw": multiplicity,
        "atom_indices": sorted(set(atoms)),
        "raw_signal": raw.strip(),
    }


def local_spectrum_rows(mol: Chem.Mol) -> Iterator[dict[str, Any]]:
    from tools.shiftdb_atom_identity import canonical_identity, identity_for_current, source_peaks, parse_fragments
    from tools.filter_nmrshiftdb2_experimental import classify
    mapping_failure = None
    try:
        identity = canonical_identity(mol, 'shiftdb:' + (mol.GetProp('nmrshiftdb2 ID') if mol.HasProp('nmrshiftdb2 ID') else 'record'), sanitize=True)
    except Exception as exc:
        # Exact-structure searching may retain legacy valence structures, but
        # none of their unresolved assignment integers become canonical labels.
        mol.UpdatePropertyCache(strict=False)
        identity = identity_for_current(mol)
        mapping_failure = str(exc)
    props = {n:mol.GetProp(n) for n in mol.GetPropNames()}
    fields = []
    for name in mol.GetPropNames():
        match = SPECTRUM_RE.match(name)
        if match and match.group(1) in {"1H", "13C"}:
            fields.append((match.group(1), int(match.group(2)), name))
    record_conflict = any(source_peaks(mol.GetProp(f), '13C', identity)[2]['status'] == 'source_assignment_inconsistency'
                          for n,s,f in fields if n=='13C' and classify(props,str(s),n)[0]=='experimental')
    for nucleus, numeric_id, field in sorted(fields, key=lambda row: (row[0], row[1])):
        raw_peaks, labels, validation = source_peaks(mol.GetProp(field), nucleus, identity,
            record_conflict=record_conflict and nucleus=='1H')
        if mapping_failure:
            validation = {**validation, 'status':'atom_identity_mapping_failed', 'mapping_error':mapping_failure}
            labels = {}
        peaks = []
        for p in sorted(raw_peaks, key=lambda item:item['shift_ppm'], reverse=True):
            peaks.append({**p, 'shift_text':p['source_assignment_fragments'][0]['shift_text'],
                          'multiplicity_raw': parse_local_signal(p['raw_signal'].split('|')[0], nucleus)['multiplicity_raw'],
                          'atom_indices':labels.get(p['peak_id'],{}).get('atom_indices',[]),
                          'atom_index_basis':'zero_based_canonical_smiles',
                          'atom_identity':identity.to_dict(), 'source_assignment_validation':validation})
        invalid = sum(f['shift_ppm'] is None for f in parse_fragments(mol.GetProp(field),nucleus,identity))
        spectrum_id = str(numeric_id)
        yield {
            "nucleus": nucleus,
            "spectrum_id": spectrum_id,
            "source_field": field,
            "peaks": peaks,
            "invalid_signal_count": invalid,
            "solvent_raw": indexed_metadata(mol, "Solvent", spectrum_id),
            "temperature_k": safe_float(indexed_metadata(mol, "Temperature [K]", spectrum_id)),
            "field_strength_mhz": safe_float(indexed_metadata(mol, "Field Strength [MHz]", spectrum_id)),
        }


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def classify_bmrb_nucleus(row: Mapping[str, str]) -> str | None:
    if clean_text(row.get("Number_of_spectral_dimensions")) != "1":
        return None
    text = " ".join(
        str(row.get(key, "")) for key in ("Sf_framecode", "Experiment_name", "Name", "Details")
    ).lower()
    if any(tag in text for tag in ("hsqc", "hmbc", "cosy", "tocsy", "dept")):
        return None
    if re.search(r"(?:^|[^0-9])1d\s*1h(?:$|[^0-9])", text) or "spectral_peak_1h" in text or "spectral_peaks_1d_1h" in text:
        return "1H"
    if re.search(r"(?:^|[^0-9])1d\s*13c(?:$|[^0-9])", text) or "spectral_peak_13c" in text or "spectral_peaks_1d_13c" in text:
        return "13C"
    return None


def concentration_to_molar(value: Any, unit: Any) -> float | None:
    amount = safe_float(value)
    units = (clean_text(unit) or "").replace("μ", "u").replace("µ", "u").lower()
    factors = {"m": 1.0, "mol/l": 1.0, "mm": 1e-3, "mmol/l": 1e-3, "um": 1e-6, "umol/l": 1e-6, "nm": 1e-9, "nmol/l": 1e-9}
    if amount is None or units not in factors:
        return None
    return amount * factors[units]


def condition_value(conditions: Mapping[str, Any], key: str) -> Any:
    value = conditions.get(key)
    if isinstance(value, Mapping):
        return value.get("normalized_value")
    return value


def connect_db(output_dir: Path) -> sqlite3.Connection:
    output_dir.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(output_dir / "bmrb_enrichment.sqlite", timeout=120)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS progress(
            stage TEXT PRIMARY KEY, status TEXT NOT NULL, processed INTEGER NOT NULL DEFAULT 0,
            total INTEGER, updated_at TEXT NOT NULL, details TEXT
        );
        CREATE TABLE IF NOT EXISTS local_records(
            record_index INTEGER PRIMARY KEY, valid INTEGER NOT NULL, invalid_reason TEXT,
            nmrshiftdb2_id TEXT, inchi TEXT, inchi_key TEXT, connectivity_key TEXT,
            smiles TEXT, formula TEXT, molecular_weight REAL, ring_count INTEGER,
            spectrum_count INTEGER NOT NULL DEFAULT 0
        );
        CREATE INDEX IF NOT EXISTS idx_local_key ON local_records(inchi_key);
        CREATE TABLE IF NOT EXISTS local_spectra(
            id INTEGER PRIMARY KEY AUTOINCREMENT, record_index INTEGER NOT NULL,
            nucleus TEXT NOT NULL, spectrum_id TEXT NOT NULL, source_field TEXT NOT NULL,
            peaks_json TEXT NOT NULL, peak_count INTEGER NOT NULL,
            invalid_signal_count INTEGER NOT NULL DEFAULT 0,
            solvent_raw TEXT, solvent_norm TEXT, temperature_k REAL, field_strength_mhz REAL,
            status TEXT NOT NULL DEFAULT 'unmatched', best_bmrb_spectrum_id INTEGER,
            best_candidate_id INTEGER, decision_reason TEXT, UNIQUE(record_index,nucleus,spectrum_id),
            FOREIGN KEY(record_index) REFERENCES local_records(record_index)
        );
        CREATE INDEX IF NOT EXISTS idx_local_spectrum_nucleus ON local_spectra(nucleus);
        CREATE TABLE IF NOT EXISTS bmrb_compounds(
            entry_id TEXT PRIMARY KEY, comp_id TEXT, name TEXT, inchi TEXT, inchi_key TEXT,
            connectivity_key TEXT, smiles TEXT, formula TEXT, formal_charge TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_bmrb_key ON bmrb_compounds(inchi_key);
        CREATE TABLE IF NOT EXISTS bmrb_spectra(
            id INTEGER PRIMARY KEY AUTOINCREMENT, entry_id TEXT NOT NULL, peak_list_id TEXT NOT NULL,
            sf_id TEXT, experiment_id TEXT, sample_id TEXT, condition_list_id TEXT,
            spectrometer_id TEXT, nucleus TEXT NOT NULL, experiment_name TEXT,
            peaks_json TEXT NOT NULL, peak_count INTEGER NOT NULL, conditions_json TEXT NOT NULL,
            UNIQUE(entry_id,peak_list_id), FOREIGN KEY(entry_id) REFERENCES bmrb_compounds(entry_id)
        );
        CREATE INDEX IF NOT EXISTS idx_bmrb_spectrum_entry_nucleus ON bmrb_spectra(entry_id,nucleus);
        CREATE TABLE IF NOT EXISTS candidates(
            id INTEGER PRIMARY KEY AUTOINCREMENT, local_spectrum_id INTEGER NOT NULL,
            bmrb_spectrum_id INTEGER NOT NULL, rank INTEGER, matched_count INTEGER NOT NULL,
            local_coverage REAL NOT NULL, bmrb_coverage REAL NOT NULL, shift_f1 REAL NOT NULL,
            count_ratio REAL NOT NULL, median_error REAL, p95_error REAL,
            solvent_conflict INTEGER NOT NULL, temperature_conflict INTEGER NOT NULL,
            field_conflict INTEGER NOT NULL, structure_conflict INTEGER NOT NULL,
            decision TEXT, reason TEXT, mappings_json TEXT NOT NULL,
            UNIQUE(local_spectrum_id,bmrb_spectrum_id),
            FOREIGN KEY(local_spectrum_id) REFERENCES local_spectra(id),
            FOREIGN KEY(bmrb_spectrum_id) REFERENCES bmrb_spectra(id)
        );
        CREATE INDEX IF NOT EXISTS idx_candidates_local ON candidates(local_spectrum_id,rank);
        """
    )
    return conn


def set_meta(conn: sqlite3.Connection, key: str, value: Any) -> None:
    conn.execute("INSERT OR REPLACE INTO meta(key,value) VALUES(?,?)", (key, json_dumps(value) if not isinstance(value, str) else value))


def require_current_assignment_index(conn):
    from tools.shiftdb_atom_identity import POLICY as IDENTITY_POLICY
    row = conn.execute("SELECT value FROM meta WHERE key='assignment_identity_policy'").fetchone()
    if not row or row[0] != IDENTITY_POLICY:
        raise RuntimeError('Assignment identity policy changed. Rebuild the local index in a NEW output directory before matching/exporting.')


def set_progress(conn: sqlite3.Connection, stage: str, status: str, processed: int, total: int | None = None, details: str | None = None) -> None:
    conn.execute(
        "INSERT OR REPLACE INTO progress(stage,status,processed,total,updated_at,details) VALUES(?,?,?,?,?,?)",
        (stage, status, processed, total, utc_now(), details),
    )


def build_manifest(source_sdf: Path, bmrb_dir: Path, output_dir: Path, policy: MatchPolicy) -> dict[str, Any]:
    http_metadata_path = bmrb_dir / "http_metadata.json"
    http_metadata = json.loads(http_metadata_path.read_text(encoding="utf-8")) if http_metadata_path.exists() else {}
    files = []
    for name in REQUIRED_TABLES:
        path = bmrb_dir / name
        if path.exists() and path.stat().st_size:
            files.append({
                "name": name, "url": f"{BMRB_BASE_URL}/{name}", "size": path.stat().st_size,
                "sha256": sha256_file(path), "last_modified": http_metadata.get(name, {}).get("last_modified"),
                "etag": http_metadata.get(name, {}).get("etag"),
                "snapshot_file_mtime_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(),
            })
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at": utc_now(),
        "source_sdf": {"path": str(source_sdf.resolve()), "size": source_sdf.stat().st_size, "sha256": sha256_file(source_sdf)},
        "bmrb_base_url": BMRB_BASE_URL,
        "bmrb_files": files,
        "match_policy": policy.__dict__,
        "tool": str(Path(__file__).resolve()),
        "python": sys.version,
        "rdkit": getattr(sys.modules.get("rdkit"), "__version__", "unknown"),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return manifest


def run_fetch(args: argparse.Namespace) -> None:
    bmrb_dir = Path(args.bmrb_dir).expanduser().resolve()
    bmrb_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = bmrb_dir / "http_metadata.json"
    http_metadata = json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.exists() else {}
    for name in REQUIRED_TABLES:
        target = bmrb_dir / name
        if target.exists() and target.stat().st_size > 0 and not args.force:
            log(f"fetch: keep existing {name} ({target.stat().st_size:,} bytes)")
            if name not in http_metadata or http_metadata[name].get("metadata_error"):
                try:
                    head = urllib.request.Request(f"{BMRB_BASE_URL}/{name}", method="HEAD", headers={"User-Agent": "NMRPeak-BMRB-enrichment/1.0"})
                    with urllib.request.urlopen(head, timeout=60) as response:
                        http_metadata[name] = {"last_modified": response.headers.get("Last-Modified"), "etag": response.headers.get("ETag"), "checked_at": utc_now()}
                except Exception as exc:
                    http_metadata[name] = {"last_modified": None, "etag": None, "checked_at": utc_now(), "metadata_error": type(exc).__name__}
        else:
            temporary = target.with_suffix(target.suffix + ".partial")
            request = urllib.request.Request(f"{BMRB_BASE_URL}/{name}", headers={"User-Agent": "NMRPeak-BMRB-enrichment/1.0"})
            log(f"fetch: downloading {name}")
            with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as handle:
                shutil.copyfileobj(response, handle)
                http_metadata[name] = {"last_modified": response.headers.get("Last-Modified"), "etag": response.headers.get("ETag"), "checked_at": utc_now()}
            if temporary.stat().st_size == 0:
                raise RuntimeError(f"Downloaded empty BMRB table: {name}")
            temporary.replace(target)
    metadata_path.write_text(json.dumps(http_metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest = build_manifest(Path(args.source_sdf), bmrb_dir, Path(args.output_dir), POLICIES[args.match_policy])
    log(f"fetch: {len(manifest['bmrb_files'])} files recorded in manifest.json")


def prefer_smiles(rows: Sequence[Mapping[str, str]]) -> str | None:
    ranked = sorted(rows, key=lambda row: 0 if str(row.get("Type", "")).lower() == "isomeric" else 1)
    return clean_text(ranked[0].get("String")) if ranked else None


def run_index_bmrb(conn: sqlite3.Connection, bmrb_dir: Path) -> None:
    log("index: loading BMRB relational tables")
    tables = {name: read_csv(bmrb_dir / name) for name in REQUIRED_TABLES}
    smiles_by_entry: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tables["Chem_comp_SMILES.csv"]:
        smiles_by_entry[row["Entry_ID"]].append(row)

    with conn:
        conn.execute("DELETE FROM candidates")
        conn.execute("DELETE FROM bmrb_spectra")
        conn.execute("DELETE FROM bmrb_compounds")
        for row in tables["Chem_comp.csv"]:
            entry = row["Entry_ID"]
            inchi = clean_text(row.get("InChI_code"))
            key = None
            if inchi:
                try:
                    key = Chem.InchiToInchiKey(inchi)
                except Exception:
                    key = None
            smiles = prefer_smiles(smiles_by_entry.get(entry, []))
            conn.execute(
                "INSERT OR REPLACE INTO bmrb_compounds VALUES(?,?,?,?,?,?,?,?,?)",
                (entry, row.get("ID"), clean_text(row.get("Name")), inchi, key, key[:14] if key else None, smiles, normalise_formula(row.get("Formula")), clean_text(row.get("Formal_charge"))),
            )

    experiments = {(row["Entry_ID"], row["ID"]): row for row in tables["Experiment.csv"]}
    spectrometers = {(row["Entry_ID"], row["ID"]): row for row in tables["NMR_spectrometer.csv"]}
    samples = {(row["Entry_ID"], row["ID"]): row for row in tables["Sample.csv"]}
    components: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in tables["Sample_component.csv"]:
        components[(row["Entry_ID"], row["Sample_ID"])].append(row)
    variables: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in tables["Sample_condition_variable.csv"]:
        variables[(row["Entry_ID"], row["Sample_condition_list_ID"])].append(row)
    peak_chars: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in tables["Peak_char.csv"]:
        if clean_text(row.get("Spectral_dim_ID")) in {None, "1"} and safe_float(row.get("Chem_shift_val")) is not None:
            peak_chars[(row["Entry_ID"], row["Spectral_peak_list_ID"])].append(row)
    couplings: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in tables["Peak_coupling.csv"]:
        couplings[(row["Entry_ID"], row["Spectral_peak_list_ID"], row["Peak_ID"])].append(row)

    inserted = 0
    with conn:
        for row in tables["Spectral_peak_list.csv"]:
            nucleus = classify_bmrb_nucleus(row)
            if nucleus is None:
                continue
            entry, peak_list_id = row["Entry_ID"], row["ID"]
            chars = peak_chars.get((entry, peak_list_id), [])
            if not chars:
                continue
            exp = experiments.get((entry, row.get("Experiment_ID", "")), {})
            sample_id = clean_text(row.get("Sample_ID")) or clean_text(exp.get("Sample_ID"))
            condition_id = clean_text(row.get("Sample_condition_list_ID")) or clean_text(exp.get("Sample_condition_list_ID"))
            spectrometer_id = clean_text(exp.get("NMR_spectrometer_ID"))
            sample = samples.get((entry, sample_id or ""), {})
            component_rows = components.get((entry, sample_id or ""), [])
            solvent_names = [clean_text(c.get("Mol_common_name")) for c in component_rows if (clean_text(c.get("Type")) or "").lower() == "solvent"]
            solvent_raw = clean_text(sample.get("Solvent_system")) or "+".join(filter(None, solvent_names)) or None
            condition_rows = variables.get((entry, condition_id or ""), [])
            ph = temperature = None
            for variable in condition_rows:
                kind = (clean_text(variable.get("Type")) or "").lower()
                if kind == "ph":
                    ph = safe_float(variable.get("Val"))
                elif "temp" in kind:
                    temperature = safe_float(variable.get("Val"))
                    units = (clean_text(variable.get("Val_units")) or "K").lower()
                    if temperature is not None and units in {"c", "celsius", "deg c", "°c"}:
                        temperature += 273.15
            spectrometer = spectrometers.get((entry, spectrometer_id or ""), {})
            field = safe_float(spectrometer.get("Field_strength"))
            solutes = [c for c in component_rows if (clean_text(c.get("Type")) or "").lower() == "solute"]
            concentration = None
            concentration_ambiguous = len(solutes) != 1
            if len(solutes) == 1:
                solute = solutes[0]
                concentration = {
                    "raw_value": clean_text(solute.get("Concentration_val")),
                    "raw_unit": clean_text(solute.get("Concentration_val_units")),
                    "normalized_value": concentration_to_molar(solute.get("Concentration_val"), solute.get("Concentration_val_units")),
                    "normalized_unit": "M" if concentration_to_molar(solute.get("Concentration_val"), solute.get("Concentration_val_units")) is not None else None,
                    "component_id": clean_text(solute.get("ID")),
                    "component_name": clean_text(solute.get("Mol_common_name")),
                }
            peaks = []
            for char in sorted(chars, key=lambda item: safe_float(item.get("Chem_shift_val")) or 0.0, reverse=True):
                peak_id = char["Peak_ID"]
                j_values = sorted({value for c in couplings.get((entry, peak_list_id, peak_id), []) if (value := safe_float(c.get("Coupling_val"))) is not None and value > 0})
                peaks.append(
                    {
                        "peak_id": peak_id,
                        "shift_ppm": safe_float(char.get("Chem_shift_val")),
                        "shift_error": safe_float(char.get("Chem_shift_val_err")),
                        "multiplicity_raw": clean_text(char.get("Coupling_pattern")),
                        "j_values_hz": j_values,
                        "line_width_raw": clean_text(char.get("Line_width_val")),
                        "line_width_hz": None,
                        "line_width_unit_status": "unit_not_exposed_by_relational_table",
                        "intensity_raw": [clean_text(c.get("Intensity_val")) for c in couplings.get((entry, peak_list_id, peak_id), []) if clean_text(c.get("Intensity_val"))],
                    }
                )
            conditions = {
                "solvent": {"raw_value": solvent_raw, "normalized_value": normalise_solvent(solvent_raw)},
                "temperature_k": {"raw_value": temperature, "normalized_value": temperature},
                "field_strength_mhz": {"raw_value": clean_text(spectrometer.get("Field_strength")), "normalized_value": field},
                "ph": {"raw_value": ph, "normalized_value": ph},
                "analyte_concentration": concentration,
                "concentration_ambiguous": concentration_ambiguous,
            }
            conn.execute(
                """INSERT OR REPLACE INTO bmrb_spectra(
                    entry_id,peak_list_id,sf_id,experiment_id,sample_id,condition_list_id,
                    spectrometer_id,nucleus,experiment_name,peaks_json,peak_count,conditions_json
                ) VALUES(?,?,?,?,?,?,?,?,?,?,?,?)""",
                (entry, peak_list_id, row.get("Sf_ID"), row.get("Experiment_ID"), sample_id, condition_id, spectrometer_id, nucleus, clean_text(row.get("Experiment_name")), json_dumps(peaks), len(peaks), json_dumps(conditions)),
            )
            inserted += 1
        set_meta(conn, "bmrb_indexed_at", utc_now())
        set_progress(conn, "index_bmrb", "complete", inserted, inserted)
    actual_spectra = conn.execute("SELECT COUNT(*) FROM bmrb_spectra").fetchone()[0]
    log(f"index: BMRB compounds={len(tables['Chem_comp.csv']):,}; eligible 1D spectra={actual_spectra:,}")


def mol_identifiers(mol: Chem.Mol) -> tuple[str | None, str | None, str | None, str | None]:
    inchi = clean_text(mol.GetProp("INChI")) if mol.HasProp("INChI") else None
    key = clean_text(mol.GetProp("INChI key")) if mol.HasProp("INChI key") else None
    smiles = None
    formula = normalise_formula(inchi.split("/")[1]) if inchi and "/" in inchi else None
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        inchi = inchi or Chem.MolToInchi(mol)
        key = key or (Chem.InchiToInchiKey(inchi) if inchi else None)
        formula = formula or rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        pass
    return inchi, key, smiles, formula


def run_index_local(conn: sqlite3.Connection, source_sdf: Path, resume: bool, commit_interval: int) -> None:
    from tools.shiftdb_atom_identity import POLICY as IDENTITY_POLICY
    prior_identity = conn.execute("SELECT value FROM meta WHERE key='assignment_identity_policy'").fetchone()
    if resume and conn.execute('SELECT 1 FROM local_records LIMIT 1').fetchone() and (not prior_identity or prior_identity[0] != IDENTITY_POLICY):
        raise RuntimeError('Old local index lacks verified atom identities; use a NEW output directory and rebuild index, not --resume.')
    set_meta(conn, 'assignment_identity_policy', IDENTITY_POLICY)
    source_hash = sha256_file(source_sdf)
    previous_hash = conn.execute("SELECT value FROM meta WHERE key='source_sdf_sha256'").fetchone()
    if previous_hash and previous_hash[0] != source_hash:
        raise RuntimeError("Source SDF hash changed since the previous index run")
    if not resume:
        with conn:
            conn.execute("DELETE FROM candidates")
            conn.execute("DELETE FROM local_spectra")
            conn.execute("DELETE FROM local_records")
            conn.execute("DELETE FROM progress WHERE stage IN ('index_local','match','export','report')")
    last_row = conn.execute("SELECT MAX(record_index) FROM local_records").fetchone()[0]
    last_index = int(last_row) if resume and last_row is not None else -1
    # Parse first and treat sanitisation quality separately.  Hundreds of
    # NMRShiftDB2 structures contain legacy valence representations that RDKit
    # cannot fully sanitise, but their database InChIKey and spectra remain
    # usable for exact-key lookup.  Only records that cannot be parsed at all
    # are classified as invalid.
    with source_sdf.open("rb") as count_handle:
        total = sum(1 for line in count_handle if line.strip() == b"$$$$")
    set_meta(conn, "source_sdf_sha256", source_hash)
    set_meta(conn, "source_sdf_path", str(source_sdf.resolve()))
    set_progress(conn, "index_local", "running", max(last_index + 1, 0), total)
    conn.commit()
    valid = invalid = spectra_count = 0
    with source_sdf.open("rb") as sdf_handle:
        supplier = Chem.ForwardSDMolSupplier(sdf_handle, sanitize=False, removeHs=False, strictParsing=False)
        for index, mol in enumerate(supplier):
            if index <= last_index:
                continue
            if mol is None:
                conn.execute("INSERT OR REPLACE INTO local_records(record_index,valid,invalid_reason) VALUES(?,0,?)", (index, "rdkit_parse_failed"))
                invalid += 1
            else:
                try:
                    inchi, key, smiles, formula = mol_identifiers(mol)
                    spectra = list(local_spectrum_rows(mol))
                    nmr_id = clean_text(mol.GetProp("nmrshiftdb2 ID")) if mol.HasProp("nmrshiftdb2 ID") else None
                    molecular_weight = ring_count = None
                    try:
                        molecular_weight = Descriptors.MolWt(mol)
                        ring_count = rdMolDescriptors.CalcNumRings(mol)
                    except Exception:
                        pass
                    conn.execute(
                        """INSERT OR REPLACE INTO local_records(
                            record_index,valid,nmrshiftdb2_id,inchi,inchi_key,connectivity_key,smiles,
                            formula,molecular_weight,ring_count,spectrum_count
                        ) VALUES(?,1,?,?,?,?,?,?,?,?,?)""",
                        (index, nmr_id, inchi, key, key[:14] if key else None, smiles, normalise_formula(formula), molecular_weight, ring_count, len(spectra)),
                    )
                    for spectrum in spectra:
                        conn.execute(
                            """INSERT OR REPLACE INTO local_spectra(
                                record_index,nucleus,spectrum_id,source_field,peaks_json,peak_count,
                                invalid_signal_count,solvent_raw,solvent_norm,temperature_k,field_strength_mhz,
                                status,decision_reason
                            ) VALUES(?,?,?,?,?,?,?,?,?,?,?,'unmatched',NULL)""",
                            (index, spectrum["nucleus"], spectrum["spectrum_id"], spectrum["source_field"], json_dumps(spectrum["peaks"]), len(spectrum["peaks"]), spectrum["invalid_signal_count"], spectrum["solvent_raw"], normalise_solvent(spectrum["solvent_raw"]), spectrum["temperature_k"], spectrum["field_strength_mhz"]),
                        )
                        spectra_count += 1
                    valid += 1
                except Exception as exc:
                    conn.execute("INSERT OR REPLACE INTO local_records(record_index,valid,invalid_reason) VALUES(?,0,?)", (index, f"sanitize_or_index_failed:{type(exc).__name__}"))
                    invalid += 1
            if (index + 1) % commit_interval == 0:
                set_progress(conn, "index_local", "running", index + 1, total)
                conn.commit()
                log(f"index: local SDF {index + 1:,}/{total:,}")
    with conn:
        set_progress(conn, "index_local", "complete", total, total)
        set_meta(conn, "local_indexed_at", utc_now())
    counts = conn.execute("SELECT SUM(valid),SUM(CASE WHEN valid=0 THEN 1 ELSE 0 END) FROM local_records").fetchone()
    spectra_total = conn.execute("SELECT COUNT(*) FROM local_spectra").fetchone()[0]
    log(f"index: local valid={counts[0]:,}; invalid={counts[1]:,}; spectra={spectra_total:,}")


def percentile95(values: Sequence[float]) -> float | None:
    if not values:
        return None
    return float(np.percentile(np.asarray(values, dtype=float), 95))


def align_peaks(local_peaks: Sequence[Mapping[str, Any]], bmrb_peaks: Sequence[Mapping[str, Any]], tolerance: float, nucleus: str) -> dict[str, Any]:
    local = [float(p["shift_ppm"]) for p in local_peaks if safe_float(p.get("shift_ppm")) is not None]
    remote = [float(p["shift_ppm"]) for p in bmrb_peaks if safe_float(p.get("shift_ppm")) is not None]
    if not local or not remote:
        return {"matched_count": 0, "local_coverage": 0.0, "bmrb_coverage": 0.0, "shift_f1": 0.0, "count_ratio": 0.0, "median_error": None, "p95_error": None, "mappings": []}
    n, m = len(local), len(remote)
    size = n + m
    cost = np.zeros((size, size), dtype=float)
    cost[:n, :m] = 4.0
    for i, left in enumerate(local):
        for j, right in enumerate(remote):
            diff = abs(left - right)
            if diff <= tolerance:
                cost[i, j] = diff / tolerance
    cost[:n, m:] = 1.0
    cost[n:, :m] = 1.0
    rows, cols = linear_sum_assignment(cost)
    mappings = []
    uniqueness_gap = 0.002 if nucleus == "1H" else 0.01
    for i, j in zip(rows.tolist(), cols.tolist()):
        if i >= n or j >= m:
            continue
        diff = abs(local[i] - remote[j])
        if diff > tolerance:
            continue
        alternatives = sorted(abs(local[i] - value) for k, value in enumerate(remote) if k != j and abs(local[i] - value) <= tolerance)
        unique = not alternatives or alternatives[0] - diff >= uniqueness_gap
        mappings.append({"local_peak_index": i, "bmrb_peak_index": j, "error_ppm": diff, "annotation_unique": unique})
    errors = [item["error_ppm"] for item in mappings]
    local_cov = len(mappings) / n
    bmrb_cov = len(mappings) / m
    f1 = 2 * local_cov * bmrb_cov / (local_cov + bmrb_cov) if local_cov + bmrb_cov else 0.0
    return {
        "matched_count": len(mappings),
        "local_coverage": local_cov,
        "bmrb_coverage": bmrb_cov,
        "shift_f1": f1,
        "count_ratio": min(n, m) / max(n, m),
        "median_error": statistics.median(errors) if errors else None,
        "p95_error": percentile95(errors),
        "mappings": mappings,
    }


def condition_conflicts(local: sqlite3.Row | Mapping[str, Any], remote: Mapping[str, Any], policy: MatchPolicy) -> tuple[bool, bool, bool]:
    remote_solvent = condition_value(remote, "solvent")
    solvent_conflict = bool(local["solvent_norm"] and remote_solvent and normalise_solvent(remote_solvent) != local["solvent_norm"])
    remote_temp = safe_float(condition_value(remote, "temperature_k"))
    remote_field = safe_float(condition_value(remote, "field_strength_mhz"))
    temperature_conflict = bool(local["temperature_k"] is not None and remote_temp is not None and abs(local["temperature_k"] - remote_temp) > policy.temperature_tolerance_k)
    field_conflict = bool(local["field_strength_mhz"] is not None and remote_field is not None and abs(local["field_strength_mhz"] - remote_field) > policy.field_tolerance_mhz)
    return solvent_conflict, temperature_conflict, field_conflict


def candidate_is_strict(metrics: Mapping[str, Any], nucleus: str, policy: MatchPolicy) -> bool:
    median_max = policy.proton_median_max if nucleus == "1H" else policy.carbon_median_max
    p95_max = policy.proton_p95_max if nucleus == "1H" else policy.carbon_p95_max
    return (
        metrics["local_coverage"] >= policy.min_coverage
        and metrics["bmrb_coverage"] >= policy.min_coverage
        and metrics["count_ratio"] >= policy.min_count_ratio
        and metrics["median_error"] is not None and metrics["median_error"] <= median_max
        and metrics["p95_error"] is not None and metrics["p95_error"] <= p95_max
    )


def candidate_is_review(metrics: Mapping[str, Any], nucleus: str, policy: MatchPolicy) -> bool:
    median_max = 0.015 if nucleus == "1H" else 0.125
    return (
        metrics["local_coverage"] >= policy.review_coverage
        and metrics["bmrb_coverage"] >= policy.review_coverage
        and metrics["count_ratio"] >= policy.review_count_ratio
        and metrics["median_error"] is not None and metrics["median_error"] <= median_max
    )


def run_match(conn: sqlite3.Connection, policy: MatchPolicy, resume: bool, commit_interval: int) -> None:
    require_current_assignment_index(conn)
    # Recompute canonical solvents so policy/alias improvements do not require
    # a costly SDF re-index.
    with conn:
        for row in conn.execute("SELECT id,solvent_raw,solvent_norm FROM local_spectra WHERE solvent_raw IS NOT NULL"):
            canonical = normalise_solvent(row["solvent_raw"])
            if canonical != row["solvent_norm"]:
                conn.execute("UPDATE local_spectra SET solvent_norm=? WHERE id=?", (canonical, row["id"]))
    if not resume:
        with conn:
            conn.execute("DELETE FROM candidates")
            conn.execute("UPDATE local_spectra SET status='unmatched',best_bmrb_spectrum_id=NULL,best_candidate_id=NULL,decision_reason=NULL")
            conn.execute("DELETE FROM progress WHERE stage='match'")
    last = conn.execute("SELECT processed FROM progress WHERE stage='match' AND status='running'").fetchone()
    last_id = int(last[0]) if resume and last else 0
    total = conn.execute("SELECT COUNT(*) FROM local_spectra WHERE peak_count>0").fetchone()[0]
    rows = conn.execute(
        """SELECT s.*,r.inchi_key,r.formula FROM local_spectra s
           JOIN local_records r ON r.record_index=s.record_index
           WHERE s.peak_count>0 AND s.id>? ORDER BY s.id""",
        (last_id,),
    )
    processed = 0
    for local in rows:
        local_peaks = json.loads(local["peaks_json"])
        if not local["inchi_key"]:
            conn.execute("UPDATE local_spectra SET status='no-structure-candidate',decision_reason='missing_local_inchikey' WHERE id=?", (local["id"],))
        else:
            exact_compounds = conn.execute("SELECT entry_id,formula FROM bmrb_compounds WHERE inchi_key=?", (local["inchi_key"],)).fetchall()
            entry_ids = [row["entry_id"] for row in exact_compounds if not row["formula"] or not local["formula"] or row["formula"] == local["formula"]]
            if not entry_ids:
                conn.execute("UPDATE local_spectra SET status='no-structure-candidate',decision_reason='no_exact_full_inchikey' WHERE id=?", (local["id"],))
            else:
                placeholders = ",".join("?" for _ in entry_ids)
                remote_rows = conn.execute(
                    f"SELECT * FROM bmrb_spectra WHERE nucleus=? AND entry_id IN ({placeholders})",
                    (local["nucleus"], *entry_ids),
                ).fetchall()
                if not remote_rows:
                    conn.execute("UPDATE local_spectra SET status='C_reject',decision_reason='no_same_nucleus_1d_spectrum' WHERE id=?", (local["id"],))
                else:
                    candidates = []
                    tolerance = policy.proton_tolerance if local["nucleus"] == "1H" else policy.carbon_tolerance
                    for remote in remote_rows:
                        remote_peaks = json.loads(remote["peaks_json"])
                        metrics = align_peaks(local_peaks, remote_peaks, tolerance, local["nucleus"])
                        remote_conditions = json.loads(remote["conditions_json"])
                        solvent_conflict, temp_conflict, field_conflict = condition_conflicts(local, remote_conditions, policy)
                        metrics.update({"remote": remote, "solvent_conflict": solvent_conflict, "temperature_conflict": temp_conflict, "field_conflict": field_conflict})
                        candidates.append(metrics)
                    candidates.sort(key=lambda item: (-item["shift_f1"], item["median_error"] if item["median_error"] is not None else float("inf"), item["p95_error"] if item["p95_error"] is not None else float("inf"), item["remote"]["id"]))
                    for rank, candidate in enumerate(candidates, start=1):
                        conn.execute(
                            """INSERT OR REPLACE INTO candidates(
                                local_spectrum_id,bmrb_spectrum_id,rank,matched_count,local_coverage,bmrb_coverage,
                                shift_f1,count_ratio,median_error,p95_error,solvent_conflict,temperature_conflict,
                                field_conflict,structure_conflict,decision,reason,mappings_json
                            ) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,0,NULL,NULL,?)""",
                            (local["id"], candidate["remote"]["id"], rank, candidate["matched_count"], candidate["local_coverage"], candidate["bmrb_coverage"], candidate["shift_f1"], candidate["count_ratio"], candidate["median_error"], candidate["p95_error"], int(candidate["solvent_conflict"]), int(candidate["temperature_conflict"]), int(candidate["field_conflict"]), json_dumps(candidate["mappings"])),
                        )
                    best = candidates[0]
                    second_f1 = candidates[1]["shift_f1"] if len(candidates) > 1 else None
                    margin = best["shift_f1"] - second_f1 if second_f1 is not None else 1.0
                    conflicts = [name for name, value in (("solvent_conflict", best["solvent_conflict"]), ("temperature_conflict", best["temperature_conflict"]), ("field_conflict", best["field_conflict"])) if value]
                    strict = candidate_is_strict(best, local["nucleus"], policy)
                    if strict and not conflicts and margin >= policy.min_margin:
                        status, reason = "A_auto", "strict_unique_match"
                    elif candidate_is_review(best, local["nucleus"], policy) or strict:
                        status = "B_review"
                        reason = ",".join(conflicts or (["ambiguous_best_candidate"] if margin < policy.min_margin else ["near_strict_threshold"]))
                    else:
                        status, reason = "C_reject", "insufficient_peak_agreement"
                    best_candidate = conn.execute("SELECT id FROM candidates WHERE local_spectrum_id=? AND rank=1", (local["id"],)).fetchone()[0]
                    conn.execute("UPDATE candidates SET decision=?,reason=? WHERE id=?", (status, reason, best_candidate))
                    conn.execute(
                        "UPDATE local_spectra SET status=?,best_bmrb_spectrum_id=?,best_candidate_id=?,decision_reason=? WHERE id=?",
                        (status, best["remote"]["id"], best_candidate, reason, local["id"]),
                    )
        processed += 1
        if processed % commit_interval == 0:
            set_progress(conn, "match", "running", local["id"], total)
            conn.commit()
            log(f"match: processed {processed:,} spectra in this run")
    with conn:
        max_id = conn.execute("SELECT COALESCE(MAX(id),0) FROM local_spectra").fetchone()[0]
        set_progress(conn, "match", "complete", max_id, total)
        set_meta(conn, "matched_at", utc_now())
    summary = conn.execute("SELECT status,COUNT(*) n FROM local_spectra GROUP BY status ORDER BY status").fetchall()
    log("match: " + ", ".join(f"{row['status']}={row['n']:,}" for row in summary))


def provenance_value(raw: Any, normalized: Any, source_table: str, source_row_id: str, status: str) -> dict[str, Any] | None:
    if raw is None and normalized is None:
        return None
    return {"raw_value": raw, "normalized_value": normalized, "source_table": source_table, "source_row_id": source_row_id, "transfer_status": status}


def source_condition(local: sqlite3.Row, key: str) -> dict[str, Any] | None:
    mapping = {
        "solvent": (local["solvent_raw"], local["solvent_norm"]),
        "temperature_k": (local["temperature_k"], local["temperature_k"]),
        "field_strength_mhz": (local["field_strength_mhz"], local["field_strength_mhz"]),
    }
    if key not in mapping:
        return None
    raw, normalized = mapping[key]
    return provenance_value(raw, normalized, "NMRShiftDB2", f"record:{local['record_index']}:{local['source_field']}:{key}", "original")


def bmrb_condition(remote: sqlite3.Row, conditions: Mapping[str, Any], key: str) -> dict[str, Any] | None:
    value = conditions.get(key)
    if not isinstance(value, Mapping) or (value.get("raw_value") is None and value.get("normalized_value") is None):
        return None
    result = dict(value)
    result.update({"source_table": "BMRB relational tables", "source_row_id": f"{remote['entry_id']}:peak_list:{remote['peak_list_id']}:experiment:{remote['experiment_id']}:sample:{remote['sample_id']}", "transfer_status": "filled_missing"})
    return result


def build_enriched_spectrum(conn: sqlite3.Connection, local: sqlite3.Row) -> dict[str, Any]:
    source_id = f"record:{local['record_index']}:{local['source_field']}"
    peaks = json.loads(local["peaks_json"])
    conditions = {
        key: source_condition(local, key) for key in ("solvent", "temperature_k", "field_strength_mhz")
    }
    conditions.update({"ph": None, "analyte_concentration": None})
    match_info: dict[str, Any] = {"status": local["status"], "reason": local["decision_reason"]}
    if local["status"] == "A_auto" and local["best_bmrb_spectrum_id"]:
        remote = conn.execute("SELECT * FROM bmrb_spectra WHERE id=?", (local["best_bmrb_spectrum_id"],)).fetchone()
        candidate = conn.execute("SELECT * FROM candidates WHERE id=?", (local["best_candidate_id"],)).fetchone()
        remote_conditions = json.loads(remote["conditions_json"])
        for key in conditions:
            if conditions[key] is None:
                conditions[key] = bmrb_condition(remote, remote_conditions, key)
        remote_peaks = json.loads(remote["peaks_json"])
        mappings = json.loads(candidate["mappings_json"])
        by_local = {item["local_peak_index"]: item for item in mappings}
        for index, peak in enumerate(peaks):
            peak.setdefault("provenance", {})
            peak["provenance"]["shift_ppm"] = {"source_table": "NMRShiftDB2", "source_row_id": source_id, "transfer_status": "original"}
            mapping = by_local.get(index)
            if not mapping or not mapping["annotation_unique"]:
                peak.setdefault("j_values", "_")
                peak.setdefault("line_width_hz", None)
                continue
            remote_peak = remote_peaks[mapping["bmrb_peak_index"]]
            if local["nucleus"] == "1H" and not peak.get("multiplicity_raw") and remote_peak.get("multiplicity_raw"):
                peak["multiplicity_raw"] = remote_peak["multiplicity_raw"]
                peak["category"] = str(remote_peak["multiplicity_raw"]).lower().replace(" ", "_")
                peak["provenance"]["category"] = {"source_table": "Peak_char.csv", "source_row_id": f"{remote['entry_id']}:{remote['peak_list_id']}:{remote_peak['peak_id']}", "transfer_status": "filled_missing"}
            elif local["nucleus"] == "1H":
                peak["category"] = peak.get("multiplicity_raw") or "unk"
            j_values = remote_peak.get("j_values_hz") or []
            peak["j_values"] = "_".join(f"{value:g}" for value in j_values) + ("_" if j_values else "_")
            if j_values:
                peak["provenance"]["j_values"] = {"source_table": "Peak_coupling.csv", "source_row_id": f"{remote['entry_id']}:{remote['peak_list_id']}:{remote_peak['peak_id']}", "transfer_status": "filled_missing"}
            peak["line_width_hz"] = remote_peak.get("line_width_hz")
            if remote_peak.get("line_width_raw") is not None:
                peak["line_width_raw"] = remote_peak["line_width_raw"]
                peak["provenance"]["line_width_raw"] = {"source_table": "Peak_char.csv", "source_row_id": f"{remote['entry_id']}:{remote['peak_list_id']}:{remote_peak['peak_id']}", "transfer_status": "audit_only_unit_unknown"}
        match_info.update(
            {
                "bmrb_entry_id": remote["entry_id"], "peak_list_id": remote["peak_list_id"],
                "experiment_id": remote["experiment_id"], "sample_id": remote["sample_id"],
                "condition_list_id": remote["condition_list_id"], "spectrometer_id": remote["spectrometer_id"],
                "metrics": {key: candidate[key] for key in ("matched_count", "local_coverage", "bmrb_coverage", "shift_f1", "count_ratio", "median_error", "p95_error")},
            }
        )
    for peak in peaks:
        if local["nucleus"] == "1H":
            peak.setdefault("category", peak.get("multiplicity_raw") or "unk")
        peak.setdefault("j_values", "_")
        peak.setdefault("line_width_hz", None)
    return {
        "local_spectrum_id": local["id"], "nucleus": local["nucleus"], "spectrum_id": local["spectrum_id"],
        "source_field": local["source_field"], "conditions": conditions, "peaks": peaks, "match": match_info,
    }


def stable_split(connectivity_key: str, seed: int = 2026) -> str:
    number = int(hashlib.sha256(f"{seed}:{connectivity_key}".encode()).hexdigest()[:16], 16) / float(0xFFFFFFFFFFFFFFFF)
    return "train" if number < 0.8 else "valid" if number < 0.9 else "test"


def lmdb_record(record: sqlite3.Row, spectrum: Mapping[str, Any]) -> dict[str, Any]:
    peaks = spectrum["peaks"]
    item: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "source_record_index": record["record_index"],
        "local_spectrum_id": spectrum["local_spectrum_id"],
        "smiles_canonical": record["smiles"],
        "molecular_formula": record["formula"],
        "conditions": spectrum["conditions"],
        "match": spectrum["match"],
        "c_nmr_peaks": None,
        "h_nmr_peaks": None,
    }
    if spectrum["nucleus"] == "13C":
        item["c_nmr_peaks"] = [{"delta (ppm)": peak["shift_ppm"], "line_width_hz": peak.get("line_width_hz")} for peak in peaks]
    else:
        item["h_nmr_peaks"] = [
            {"centroid": peak["shift_ppm"], "delta": peak["shift_ppm"], "rangeMax": peak["shift_ppm"], "rangeMin": peak["shift_ppm"], "nH": 1, "category": peak.get("category") or "unk", "j_values": peak.get("j_values") or "_", "line_width_hz": peak.get("line_width_hz")}
            for peak in peaks
        ]
    return item


def write_lmdb(path: Path, records: Sequence[dict[str, Any]], map_size_gb: float) -> int:
    if path.exists():
        path.unlink()
    lock = Path(str(path) + "-lock")
    if lock.exists():
        lock.unlink()
    payloads = [pickle.dumps(record, protocol=pickle.HIGHEST_PROTOCOL) for record in records]
    estimated = sum(len(payload) + 128 for payload in payloads)
    required_map_size = max(128 * 1024**2, int(estimated * 2.0))
    maximum_map_size = int(map_size_gb * 1024**3)
    if required_map_size > maximum_map_size:
        raise RuntimeError(
            f"Estimated LMDB map size {required_map_size / 1024**3:.2f} GiB exceeds "
            f"--map-size-gb={map_size_gb:.2f}"
        )
    env = lmdb.open(str(path), subdir=False, map_size=required_map_size, lock=True)
    count = 0
    with env.begin(write=True) as txn:
        for count, payload in enumerate(payloads, start=1):
            txn.put(str(count - 1).encode("ascii"), payload)
        txn.put(b"__len__", str(count).encode("ascii"))
    env.sync()
    env.close()
    return count


def run_export(conn: sqlite3.Connection, output_dir: Path, map_size_gb: float) -> None:
    require_current_assignment_index(conn)
    source_hash = conn.execute("SELECT value FROM meta WHERE key='source_sdf_sha256'").fetchone()
    source_path = conn.execute("SELECT value FROM meta WHERE key='source_sdf_path'").fetchone()
    if not source_hash or not source_path or sha256_file(Path(source_path[0])) != source_hash[0]:
        raise RuntimeError("Source SDF hash verification failed before export")
    records_path = output_dir / "records.jsonl"
    partial = output_dir / "records.jsonl.partial"
    split_records: dict[str, list[dict[str, Any]]] = {"train": [], "valid": [], "test": []}
    valid_records = conn.execute("SELECT * FROM local_records WHERE valid=1 ORDER BY record_index")
    written = 0
    with partial.open("w", encoding="utf-8", newline="\n") as handle:
        for record in valid_records:
            local_spectra = conn.execute("SELECT * FROM local_spectra WHERE record_index=? ORDER BY nucleus,spectrum_id", (record["record_index"],)).fetchall()
            spectra = [build_enriched_spectrum(conn, row) for row in local_spectra]
            output = {
                "schema_version": SCHEMA_VERSION,
                "source_record_index": record["record_index"],
                "identifiers": {"nmrshiftdb2_id": record["nmrshiftdb2_id"], "inchi": record["inchi"], "inchi_key": record["inchi_key"]},
                "molecule": {"canonical_smiles": record["smiles"], "molecular_formula": record["formula"], "molecular_weight": record["molecular_weight"], "ring_count": record["ring_count"]},
                "spectra": spectra,
            }
            handle.write(json.dumps(output, ensure_ascii=False, separators=(",", ":")) + "\n")
            split = stable_split(record["connectivity_key"] or f"record-{record['record_index']}")
            for spectrum in spectra:
                if spectrum["peaks"]:
                    split_records[split].append(lmdb_record(record, spectrum))
            written += 1
    partial.replace(records_path)
    lmdb_dir = output_dir / "lmdb_dataset"
    lmdb_dir.mkdir(parents=True, exist_ok=True)
    lmdb_counts = {split: write_lmdb(lmdb_dir / f"{split}.lmdb", rows, map_size_gb) for split, rows in split_records.items()}
    review_path = output_dir / "review_queue.csv"
    review_rows = conn.execute(
        """SELECT l.record_index,l.nucleus,l.spectrum_id,l.status,l.decision_reason,
                  b.entry_id,b.peak_list_id,c.shift_f1,c.local_coverage,c.bmrb_coverage,
                  c.count_ratio,c.median_error,c.p95_error
           FROM local_spectra l LEFT JOIN candidates c ON c.id=l.best_candidate_id
           LEFT JOIN bmrb_spectra b ON b.id=l.best_bmrb_spectrum_id
           WHERE l.status='B_review' ORDER BY l.nucleus,c.shift_f1 DESC"""
    ).fetchall()
    columns = [description[0] for description in review_rows.description] if hasattr(review_rows, "description") else ["record_index", "nucleus", "spectrum_id", "status", "decision_reason", "entry_id", "peak_list_id", "shift_f1", "local_coverage", "bmrb_coverage", "count_ratio", "median_error", "p95_error"]
    # sqlite result lists do not expose .description; use the explicit stable schema above.
    columns = ["record_index", "nucleus", "spectrum_id", "status", "decision_reason", "entry_id", "peak_list_id", "shift_f1", "local_coverage", "bmrb_coverage", "count_ratio", "median_error", "p95_error"]
    with review_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        writer.writerows([tuple(row) for row in review_rows])
    with conn:
        set_meta(conn, "exported_at", utc_now())
        set_meta(conn, "lmdb_counts", lmdb_counts)
        set_progress(conn, "export", "complete", written, written, json_dumps(lmdb_counts))
    if sha256_file(Path(source_path[0])) != source_hash[0]:
        raise RuntimeError("Source SDF changed during export")
    log(f"export: JSONL records={written:,}; LMDB={lmdb_counts}; review={len(review_rows):,}")


def wilson_lower(successes: int, total: int, z: float = 1.959963984540054) -> float | None:
    if total <= 0:
        return None
    p = successes / total
    denominator = 1 + z * z / total
    centre = p + z * z / (2 * total)
    adjustment = z * math.sqrt((p * (1 - p) + z * z / (4 * total)) / total)
    return (centre - adjustment) / denominator


def run_negative_control(conn: sqlite3.Connection, policy: MatchPolicy) -> dict[str, Any]:
    local = conn.execute("SELECT id,nucleus,peaks_json FROM local_spectra WHERE status IN ('A_auto','B_review') ORDER BY nucleus,id").fetchall()
    remote_by_nucleus = {
        nucleus: conn.execute("SELECT id,entry_id,peaks_json FROM bmrb_spectra WHERE nucleus=? ORDER BY entry_id,id", (nucleus,)).fetchall()
        for nucleus in ("1H", "13C")
    }
    tested = passed = 0
    for row in local:
        remotes = remote_by_nucleus[row["nucleus"]]
        if len(remotes) < 2:
            continue
        remote = remotes[(row["id"] * 7919) % len(remotes)]
        tolerance = policy.proton_tolerance if row["nucleus"] == "1H" else policy.carbon_tolerance
        metrics = align_peaks(json.loads(row["peaks_json"]), json.loads(remote["peaks_json"]), tolerance, row["nucleus"])
        tested += 1
        passed += int(candidate_is_strict(metrics, row["nucleus"], policy))
    return {"tested": tested, "shift_only_strict_matches": passed, "rate": passed / tested if tested else None, "note": "Structure gate is intentionally omitted; production auto-matching still requires an exact full InChIKey."}


def report_statistics(conn: sqlite3.Connection, output_dir: Path, policy: MatchPolicy) -> dict[str, Any]:
    scalar = lambda sql, params=(): conn.execute(sql, params).fetchone()[0]
    status_by_nucleus = {}
    for row in conn.execute("SELECT nucleus,status,COUNT(*) n FROM local_spectra GROUP BY nucleus,status ORDER BY nucleus,status"):
        status_by_nucleus.setdefault(row["nucleus"], {})[row["status"]] = row["n"]
    coverage = {}
    for nucleus in ("1H", "13C"):
        coverage[nucleus] = {
            "spectra": scalar("SELECT COUNT(*) FROM local_spectra WHERE nucleus=?", (nucleus,)),
            "solvent_existing": scalar("SELECT COUNT(*) FROM local_spectra WHERE nucleus=? AND solvent_norm IS NOT NULL", (nucleus,)),
            "temperature_existing": scalar("SELECT COUNT(*) FROM local_spectra WHERE nucleus=? AND temperature_k IS NOT NULL", (nucleus,)),
            "field_existing": scalar("SELECT COUNT(*) FROM local_spectra WHERE nucleus=? AND field_strength_mhz IS NOT NULL", (nucleus,)),
        }
    additions = Counter({f"{nucleus}:{field}": 0 for nucleus in ("1H", "13C") for field in ("solvent", "temperature_k", "field_strength_mhz", "ph", "analyte_concentration")})
    peak_additions = Counter({f"1H:{field}": 0 for field in ("category", "j_values", "line_width_raw", "line_width_hz")})
    peak_additions.update({"13C:line_width_raw": 0, "13C:line_width_hz": 0})
    auto_rows = conn.execute("SELECT * FROM local_spectra WHERE status='A_auto'").fetchall()
    for local in auto_rows:
        enriched = build_enriched_spectrum(conn, local)
        for key, value in enriched["conditions"].items():
            if value and value.get("transfer_status") == "filled_missing":
                additions[f"{local['nucleus']}:{key}"] += 1
        for peak in enriched["peaks"]:
            provenance = peak.get("provenance", {})
            for field in ("category", "j_values", "line_width_raw"):
                if field in provenance:
                    peak_additions[f"{local['nucleus']}:{field}"] += 1
    for nucleus in ("1H", "13C"):
        coverage[nucleus]["solvent_after"] = coverage[nucleus]["solvent_existing"] + additions[f"{nucleus}:solvent"]
        coverage[nucleus]["temperature_after"] = coverage[nucleus]["temperature_existing"] + additions[f"{nucleus}:temperature_k"]
        coverage[nucleus]["field_after"] = coverage[nucleus]["field_existing"] + additions[f"{nucleus}:field_strength_mhz"]
        coverage[nucleus]["ph_after"] = additions[f"{nucleus}:ph"]
        coverage[nucleus]["analyte_concentration_after"] = additions[f"{nucleus}:analyte_concentration"]
    best_metrics = {}
    for nucleus in ("1H", "13C"):
        rows = conn.execute(
            """SELECT c.shift_f1,c.median_error,c.p95_error FROM candidates c
               JOIN local_spectra l ON l.id=c.local_spectrum_id WHERE c.rank=1 AND l.nucleus=?""", (nucleus,)
        ).fetchall()
        best_metrics[nucleus] = {
            "count": len(rows),
            "shift_f1_median": statistics.median([row["shift_f1"] for row in rows]) if rows else None,
            "median_error_median": statistics.median([row["median_error"] for row in rows if row["median_error"] is not None]) if rows else None,
        }
    all_mw = [row[0] for row in conn.execute("SELECT molecular_weight FROM local_records WHERE valid=1 AND molecular_weight IS NOT NULL")]
    matched_mw = [row[0] for row in conn.execute("SELECT DISTINCT r.molecular_weight FROM local_records r JOIN local_spectra s ON s.record_index=r.record_index WHERE s.status='A_auto' AND r.molecular_weight IS NOT NULL")]
    all_rings = [row[0] for row in conn.execute("SELECT ring_count FROM local_records WHERE valid=1 AND ring_count IS NOT NULL")]
    matched_rings = [row[0] for row in conn.execute("SELECT DISTINCT r.ring_count FROM local_records r JOIN local_spectra s ON s.record_index=r.record_index WHERE s.status='A_auto' AND r.ring_count IS NOT NULL")]
    element_all, element_matched = Counter(), Counter()
    for formula, matched in conn.execute("""SELECT r.formula,MAX(CASE WHEN s.status='A_auto' THEN 1 ELSE 0 END)
                                         FROM local_records r LEFT JOIN local_spectra s ON s.record_index=r.record_index
                                         WHERE r.valid=1 GROUP BY r.record_index,r.formula"""):
        elements = set(re.findall(r"[A-Z][a-z]?", formula or ""))
        element_all.update(elements)
        if matched:
            element_matched.update(elements)
    solvent_all = {row[0] or "missing": row[1] for row in conn.execute("SELECT solvent_norm,COUNT(*) FROM local_spectra GROUP BY solvent_norm ORDER BY COUNT(*) DESC LIMIT 20")}
    solvent_matched = {row[0] or "missing": row[1] for row in conn.execute("SELECT solvent_norm,COUNT(*) FROM local_spectra WHERE status='A_auto' GROUP BY solvent_norm ORDER BY COUNT(*) DESC LIMIT 20")}

    original_lengths, enriched_lengths = [], []
    for row in conn.execute("SELECT * FROM local_spectra WHERE peak_count>0"):
        enriched = build_enriched_spectrum(conn, row)
        peak_count = len(enriched["peaks"])
        original = peak_count + 1
        if row["nucleus"] == "1H":
            original += 2 * peak_count  # category and nH tokens used by NMRPeak-G H mode
            original += sum(len(re.findall(r"\d+(?:\.\d+)?", peak.get("j_values", "_"))) for peak in enriched["peaks"])
        # Five condition values plus start/end markers and one optional line-width token per peak.
        enriched_length = original + 7 + peak_count
        original_lengths.append(original)
        enriched_lengths.append(enriched_length)
    invalid_reasons = {row["invalid_reason"]: row["n"] for row in conn.execute("SELECT invalid_reason,COUNT(*) n FROM local_records WHERE valid=0 GROUP BY invalid_reason")}
    rejection_reasons = {row["decision_reason"]: row["n"] for row in conn.execute("SELECT decision_reason,COUNT(*) n FROM local_spectra WHERE status IN ('B_review','C_reject','no-structure-candidate') GROUP BY decision_reason")}
    negative = run_negative_control(conn, policy)
    audit_rows = conn.execute(
        """SELECT l.record_index,l.nucleus,l.spectrum_id,b.entry_id,b.peak_list_id,c.shift_f1,c.median_error,c.p95_error
           FROM local_spectra l JOIN candidates c ON c.id=l.best_candidate_id JOIN bmrb_spectra b ON b.id=l.best_bmrb_spectrum_id
           WHERE l.status='A_auto' ORDER BY l.nucleus,c.shift_f1,l.id"""
    ).fetchall()
    by_nucleus = defaultdict(list)
    for row in audit_rows:
        by_nucleus[row["nucleus"]].append(row)
    audit_sample = []
    for nucleus, rows in by_nucleus.items():
        if len(rows) <= 100:
            audit_sample.extend(rows)
        else:
            indices = np.linspace(0, len(rows) - 1, 100, dtype=int)
            audit_sample.extend(rows[index] for index in indices)
    audit_path = output_dir / "auto_audit_sample.csv"
    audit_columns = ["record_index", "nucleus", "spectrum_id", "entry_id", "peak_list_id", "shift_f1", "median_error", "p95_error", "manual_correct"]
    with audit_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle); writer.writerow(audit_columns)
        for row in audit_sample:
            writer.writerow([*tuple(row), ""])
    report = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": utc_now(),
        "source": {
            "records_total": scalar("SELECT COUNT(*) FROM local_records"),
            "records_valid": scalar("SELECT COUNT(*) FROM local_records WHERE valid=1"),
            "records_invalid": scalar("SELECT COUNT(*) FROM local_records WHERE valid=0"),
            "invalid_reasons": invalid_reasons,
            "spectra_total": scalar("SELECT COUNT(*) FROM local_spectra"),
        },
        "bmrb": {
            "compounds": scalar("SELECT COUNT(*) FROM bmrb_compounds"),
            "eligible_1d_spectra": scalar("SELECT COUNT(*) FROM bmrb_spectra"),
            "exact_structure_overlap": scalar("SELECT COUNT(DISTINCT l.inchi_key) FROM local_records l JOIN bmrb_compounds b ON b.inchi_key=l.inchi_key WHERE l.valid=1"),
        },
        "status_by_nucleus": status_by_nucleus,
        "existing_coverage": coverage,
        "condition_additions_by_spectrum": dict(additions),
        "peak_field_additions": dict(peak_additions),
        "best_candidate_metrics": best_metrics,
        "rejection_reasons": rejection_reasons,
        "distribution_bias": {
            "all_molecular_weight_mean": statistics.fmean(all_mw) if all_mw else None,
            "auto_matched_molecular_weight_mean": statistics.fmean(matched_mw) if matched_mw else None,
            "all_molecular_weight_median": statistics.median(all_mw) if all_mw else None,
            "auto_matched_molecular_weight_median": statistics.median(matched_mw) if matched_mw else None,
            "all_ring_count_mean": statistics.fmean(all_rings) if all_rings else None,
            "auto_matched_ring_count_mean": statistics.fmean(matched_rings) if matched_rings else None,
            "element_presence_all": dict(element_all.most_common()),
            "element_presence_auto_matched": dict(element_matched.most_common()),
            "solvent_top20_all": solvent_all,
            "solvent_top20_auto_matched": solvent_matched,
        },
        "token_length_impact_all_optional_fields": {
            "original_median": statistics.median(original_lengths) if original_lengths else None,
            "enriched_median": statistics.median(enriched_lengths) if enriched_lengths else None,
            "original_p95": percentile95(original_lengths),
            "enriched_p95": percentile95(enriched_lengths),
            "original_over_512": sum(value > 512 for value in original_lengths),
            "enriched_over_512": sum(value > 512 for value in enriched_lengths),
            "definition": "Per-spectrum estimate with H category/nH/J, all five condition families, end markers, and optional per-peak line-width tokens; formula/BOS/EOS excluded.",
        },
        "negative_control": negative,
        "manual_audit": {
            "sample_file": str(audit_path), "sample_size": len(audit_sample), "status": "pending_manual_review",
            "target_wilson_95_lower_bound": 0.98, "lower_bound_if_all_200_correct": wilson_lower(200, 200),
        },
        "model_value_assessment": {
            "training_performed": False,
            "conclusion": "No model-performance improvement is claimed without controlled retraining.",
            "likely_high_value_fields": ["1H multiplicity", "1H J coupling"],
            "likely_context_fields": ["solvent", "temperature", "field strength", "pH", "concentration"],
        },
    }
    return report


def markdown_report(report: Mapping[str, Any]) -> str:
    lines = ["# NMRShiftDB2–BMRB enrichment report", "", f"Generated: {report['generated_at']}", "", "## Dataset inventory", ""]
    source, bmrb = report["source"], report["bmrb"]
    lines.extend([
        f"- NMRShiftDB2 records: {source['records_total']:,} total; {source['records_valid']:,} valid; {source['records_invalid']:,} invalid.",
        f"- Indexed NMRShiftDB2 spectra: {source['spectra_total']:,}.",
        f"- BMRB compounds / eligible 1D spectra: {bmrb['compounds']:,} / {bmrb['eligible_1d_spectra']:,}.",
        f"- Exact full-InChIKey overlap: {bmrb['exact_structure_overlap']:,} unique structures.", "",
        "## Match decisions", "", "| Nucleus | A_auto | B_review | C_reject | No structure candidate |", "|---|---:|---:|---:|---:|",
    ])
    for nucleus in ("1H", "13C"):
        values = report["status_by_nucleus"].get(nucleus, {})
        lines.append(f"| {nucleus} | {values.get('A_auto',0):,} | {values.get('B_review',0):,} | {values.get('C_reject',0):,} | {values.get('no-structure-candidate',0):,} |")
    lines.extend(["", "## Added fields", ""])
    for key, value in sorted(report["condition_additions_by_spectrum"].items()):
        lines.append(f"- `{key}`: {value:,} spectra")
    for key, value in sorted(report["peak_field_additions"].items()):
        lines.append(f"- `{key}`: {value:,} peaks")
    if not report["condition_additions_by_spectrum"] and not report["peak_field_additions"]:
        lines.append("- No fields were transferred under the strict policy.")
    negative = report["negative_control"]
    lines.extend([
        "", "## Confidence and limitations", "",
        f"- Shift-only shuffled negative control: {negative['shift_only_strict_matches']:,}/{negative['tested']:,} passed strict peak thresholds. Production matching additionally requires an exact full InChIKey.",
        f"- A stratified audit sample of {report['manual_audit']['sample_size']:,} spectra was written to `auto_audit_sample.csv`; manual correctness labels are pending.",
        "- Line widths are retained as audit-only raw values because the public relational table does not expose an explicit unit column.",
        "- No NMRPeak-G retraining was performed, so this report does not claim a model-performance gain.",
        "", "## Expected model value", "",
        "J coupling and multiplicity are the most direct topology-related additions. Experimental conditions can help explain systematic shift variation, but the exact-overlap subset is small relative to the full source dataset; any global gain must be established by a controlled ablation study.", "",
    ])
    return "\n".join(lines)


def run_report(conn: sqlite3.Connection, output_dir: Path, policy: MatchPolicy) -> None:
    report = report_statistics(conn, output_dir, policy)
    (output_dir / "enrichment_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    (output_dir / "enrichment_report.md").write_text(markdown_report(report), encoding="utf-8")
    with conn:
        set_meta(conn, "reported_at", utc_now())
        set_progress(conn, "report", "complete", report["source"]["spectra_total"], report["source"]["spectra_total"])
    log(f"report: wrote enrichment_report.md/json and {report['manual_audit']['sample_size']} audit rows")


def validate_paths(args: argparse.Namespace, require_bmrb: bool = True, require_source: bool = True) -> tuple[Path, Path, Path]:
    source = Path(args.source_sdf).expanduser().resolve()
    bmrb = Path(args.bmrb_dir).expanduser().resolve()
    output = Path(args.output_dir).expanduser().resolve()
    if require_source and not source.exists():
        raise FileNotFoundError(source)
    if require_bmrb:
        missing = [name for name in REQUIRED_TABLES if not (bmrb / name).exists() or (bmrb / name).stat().st_size == 0]
        if missing:
            raise FileNotFoundError(f"Missing BMRB tables: {', '.join(missing)}; run fetch first")
    output.mkdir(parents=True, exist_ok=True)
    return source, bmrb, output


def run_stage(args: argparse.Namespace, stage: str) -> None:
    source, bmrb, output = validate_paths(args, require_bmrb=stage != "fetch")
    policy = POLICIES[args.match_policy]
    if stage == "fetch":
        run_fetch(args); return
    conn = connect_db(output)
    try:
        if stage == "index":
            run_index_bmrb(conn, bmrb)
            run_index_local(conn, source, args.resume, args.commit_interval)
            build_manifest(source, bmrb, output, policy)
        elif stage == "match":
            run_match(conn, policy, args.resume, args.commit_interval)
        elif stage == "export":
            run_export(conn, output, args.map_size_gb)
        elif stage == "report":
            run_report(conn, output, policy)
        else:
            raise ValueError(stage)
    finally:
        conn.close()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("fetch", "index", "match", "export", "report", "all"))
    parser.add_argument("--source-sdf", default=str(DEFAULT_SOURCE_SDF))
    parser.add_argument("--bmrb-dir", default=str(DEFAULT_BMRB_DIR))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--match-policy", choices=tuple(POLICIES), default="strict-v1")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--force", action="store_true", help="Redownload BMRB files during fetch")
    parser.add_argument("--commit-interval", type=int, default=500)
    parser.add_argument("--map-size-gb", type=float, default=8.0)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.commit_interval <= 0 or args.map_size_gb <= 0:
        raise ValueError("--commit-interval and --map-size-gb must be positive")
    RDLogger.DisableLog("rdApp.*")
    stages = ("fetch", "index", "match", "export", "report") if args.stage == "all" else (args.stage,)
    for stage in stages:
        log(f"stage={stage}; policy={args.match_policy}")
        run_stage(args, stage)
        # Downstream all-mode stages should reuse completed state.
        args.resume = True
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
