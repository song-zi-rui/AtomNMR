#!/usr/bin/env python
"""Build leakage-safe NMRShiftDB2 and NMRexp datasets with manual resume.

The pipeline is intentionally split into independent stages:

    index       inspect source records and write a resumable SQLite index
    split       deduplicate NMRShiftDB2 and assign leakage-safe data splits
    materialize write deterministic, resumable LMDB files
    audit       verify split isolation, deduplication, and output counts
    status      show saved progress without opening the source datasets

No stage submits or requeues a Slurm job. Re-running the same command resumes
from the last committed SQLite/LMDB batch.
"""

import argparse
import base64
import copy
import datetime as dt
import hashlib
import importlib.util
import json
import os
import pickle
import signal
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import lmdb
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold


SCHEMA_VERSION = "atompeak_pipeline_v2_atom_identity_v1"
DEFAULT_SPLITS = ("train", "valid", "test")
STOP_REQUESTED = False
STOP_SIGNAL = None


def log(message):
    stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{stamp}] {message}", flush=True)


def request_stop(signum, _frame):
    global STOP_REQUESTED, STOP_SIGNAL
    STOP_REQUESTED = True
    STOP_SIGNAL = signal.Signals(signum).name
    log(f"Received {STOP_SIGNAL}; finishing the current record and committing the active batch")


def install_signal_handlers():
    signal.signal(signal.SIGTERM, request_stop)
    signal.signal(signal.SIGINT, request_stop)
    if hasattr(signal, "SIGUSR1"):
        signal.signal(signal.SIGUSR1, request_stop)


def load_clean_converter():
    path = Path(__file__).with_name("prepare_nmrshiftdb2_atompeak_clean.py")
    spec = importlib.util.spec_from_file_location("atompeak_clean_converter", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


clean_converter = load_clean_converter()


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def plain(value):
    if isinstance(value, dict):
        return {str(key): plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    if isinstance(value, set):
        return sorted(plain(item) for item in value)
    if isinstance(value, bytes):
        return {"bytes_b64": base64.b64encode(value).decode("ascii")}
    if hasattr(value, "tolist"):
        return plain(value.tolist())
    if hasattr(value, "item"):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def canonical_json(value):
    return json.dumps(plain(value), ensure_ascii=True, sort_keys=True, separators=(",", ":"))


def stable_hash(value):
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def path_signature(path):
    path = Path(path).expanduser().resolve()
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def find_lmdb_path(root, split):
    root = Path(root).expanduser().resolve()
    candidates = (
        root / f"{split}.lmdb",
        root / "lmdb_dataset" / f"{split}.lmdb",
        root / split / "lmdb_dataset" / f"{split}.lmdb",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    matches = sorted(root.glob(f"**/{split}.lmdb"))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise FileNotFoundError(f"Could not find {split}.lmdb below {root}")
    raise RuntimeError(f"Multiple {split}.lmdb files found below {root}: {matches}")


def semantic_config(args):
    sdf_path = Path(args.sdf_path).expanduser().resolve()
    nmrexp_paths = {split: path_signature(find_lmdb_path(args.nmrexp_dir, split)) for split in DEFAULT_SPLITS}
    return {
        "schema_version": SCHEMA_VERSION,
        "sdf": path_signature(sdf_path),
        "nmrexp": nmrexp_paths,
        "nuclei": args.nuclei,
        "train_ratio": args.train_ratio,
        "valid_ratio": args.valid_ratio,
        "seed": args.seed,
    }


def config_hash(args):
    return stable_hash(semantic_config(args))


def connect_db(output_dir):
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    db_path = output_dir / "pipeline.sqlite"
    conn = sqlite3.connect(str(db_path), timeout=120)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=FULL")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA busy_timeout=120000")
    create_schema(conn)
    return conn


def create_schema(conn):
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS source_records (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            source TEXT NOT NULL,
            source_split TEXT NOT NULL DEFAULT '',
            source_path TEXT NOT NULL,
            source_locator TEXT NOT NULL,
            source_record_id TEXT NOT NULL,
            exact_structure_key TEXT,
            connectivity_key TEXT,
            scaffold_key TEXT,
            canonical_smiles TEXT,
            spectrum_hash TEXT,
            carbon_spectrum_hash TEXT,
            proton_spectrum_hash TEXT,
            solvent TEXT,
            frequency TEXT,
            has_carbon INTEGER NOT NULL DEFAULT 0,
            has_proton INTEGER NOT NULL DEFAULT 0,
            has_assignment INTEGER NOT NULL DEFAULT 0,
            quality_score INTEGER NOT NULL DEFAULT 0,
            duplicate_of INTEGER REFERENCES source_records(id),
            keep_record INTEGER NOT NULL DEFAULT 1,
            final_split TEXT,
            exclusion_reason TEXT,
            metadata_json TEXT NOT NULL DEFAULT '{}',
            UNIQUE(source, source_path, source_locator)
        );

        CREATE INDEX IF NOT EXISTS idx_source_records_structure
            ON source_records(source, exact_structure_key, connectivity_key);
        CREATE INDEX IF NOT EXISTS idx_source_records_scaffold
            ON source_records(source, scaffold_key);
        CREATE INDEX IF NOT EXISTS idx_source_records_split
            ON source_records(source, keep_record, final_split, id);
        CREATE INDEX IF NOT EXISTS idx_source_records_dedup
            ON source_records(
                source, exact_structure_key, carbon_spectrum_hash,
                proton_spectrum_hash, solvent, frequency
            );

        CREATE TABLE IF NOT EXISTS pipeline_progress (
            stage TEXT NOT NULL,
            source TEXT NOT NULL,
            split TEXT NOT NULL,
            last_source_locator TEXT,
            last_output_index INTEGER NOT NULL DEFAULT -1,
            processed_count INTEGER NOT NULL DEFAULT 0,
            kept_count INTEGER NOT NULL DEFAULT 0,
            total_count INTEGER,
            config_hash TEXT NOT NULL,
            status TEXT NOT NULL,
            updated_at TEXT NOT NULL,
            PRIMARY KEY(stage, source, split),
            CHECK(status IN ('pending', 'running', 'paused', 'completed', 'failed'))
        );

        CREATE TABLE IF NOT EXISTS split_assignments (
            record_id INTEGER PRIMARY KEY REFERENCES source_records(id),
            final_split TEXT NOT NULL,
            exclusion_reason TEXT
        );

        CREATE TABLE IF NOT EXISTS pipeline_meta (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        """
    )
    conn.commit()


def get_progress(conn, stage, source, split):
    return conn.execute(
        "SELECT * FROM pipeline_progress WHERE stage=? AND source=? AND split=?",
        (stage, source, split),
    ).fetchone()


def check_progress_config(row, expected_hash):
    if row is not None and row["config_hash"] != expected_hash:
        raise RuntimeError(
            "Configuration differs from the saved checkpoint. Refusing to mix datasets. "
            f"Saved={row['config_hash'][:12]}, current={expected_hash[:12]}. "
            "Use --reset-stage explicitly after confirming the change."
        )


def save_progress(
    conn,
    stage,
    source,
    split,
    cfg_hash,
    status,
    last_source_locator=None,
    last_output_index=-1,
    processed_count=0,
    kept_count=0,
    total_count=None,
):
    conn.execute(
        """
        INSERT INTO pipeline_progress(
            stage, source, split, last_source_locator, last_output_index,
            processed_count, kept_count, total_count, config_hash, status, updated_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(stage, source, split) DO UPDATE SET
            last_source_locator=excluded.last_source_locator,
            last_output_index=excluded.last_output_index,
            processed_count=excluded.processed_count,
            kept_count=excluded.kept_count,
            total_count=excluded.total_count,
            config_hash=excluded.config_hash,
            status=excluded.status,
            updated_at=excluded.updated_at
        """,
        (
            stage,
            source,
            split,
            last_source_locator,
            last_output_index,
            processed_count,
            kept_count,
            total_count,
            cfg_hash,
            status,
            utc_now(),
        ),
    )


def mark_failed(conn, stage, source, split, cfg_hash, error):
    row = get_progress(conn, stage, source, split)
    with conn:
        save_progress(
            conn,
            stage,
            source,
            split,
            cfg_hash,
            "failed",
            row["last_source_locator"] if row else None,
            row["last_output_index"] if row else -1,
            row["processed_count"] if row else 0,
            row["kept_count"] if row else 0,
            row["total_count"] if row else None,
        )
        conn.execute(
            "INSERT OR REPLACE INTO pipeline_meta(key, value) VALUES (?, ?)",
            (f"last_error:{stage}:{source}:{split}", str(error)),
        )


def structure_keys_from_mol(mol):
    work = Chem.Mol(mol)
    Chem.SanitizeMol(work)
    heavy = Chem.RemoveHs(work)
    Chem.SanitizeMol(heavy)
    canonical = Chem.MolToSmiles(heavy, canonical=True, isomericSmiles=True)
    connectivity_smiles = Chem.MolToSmiles(heavy, canonical=True, isomericSmiles=False)
    try:
        inchi_key = Chem.MolToInchiKey(heavy)
    except Exception:
        inchi_key = ""
    exact_key = f"INCHIKEY:{inchi_key}" if inchi_key else f"SMILES:{canonical}"
    connectivity = inchi_key.split("-")[0] if inchi_key else connectivity_smiles
    connectivity_key = f"CONNECTIVITY:{connectivity}"
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=heavy, includeChirality=False)
    scaffold_key = f"SCAFFOLD:{scaffold}" if scaffold else connectivity_key
    return exact_key, connectivity_key, scaffold_key, canonical


def structure_keys_from_smiles(smiles):
    if isinstance(smiles, bytes):
        smiles = smiles.decode("utf-8", errors="replace")
    if not smiles:
        raise ValueError("missing_smiles")
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError("invalid_smiles")
    return structure_keys_from_mol(mol)


def branch_fingerprint(record, branch_name):
    spectrum = (record.get("spectra") or {}).get(branch_name)
    supervision = (record.get("alignment_supervision") or {}).get(branch_name) or {}
    if not spectrum:
        return None
    peaks = []
    assignments = supervision.get("assigned_atom_indices_by_peak") or []
    for index, peak in enumerate(spectrum.get("peaks") or []):
        peaks.append(
            {
                "shift_text": peak.get("shift_text"),
                "shift_ppm": peak.get("shift_ppm"),
                "extra": peak.get("extra"),
                "assigned_atom_indices": assignments[index] if index < len(assignments) else [],
            }
        )
    payload = {"peaks": peaks}
    return stable_hash(payload)


def condition_value(record, field):
    values = {}
    for branch_name, nucleus in (("carbon", "13C"), ("proton", "1H")):
        spectrum = (record.get("spectra") or {}).get(branch_name)
        if spectrum:
            raw = (spectrum.get("metadata") or {}).get(field)
            values[nucleus] = clean_converter.base.normalise_condition(str(raw)) if raw is not None else None
    if values and len(set(values.values())) == 1:
        return next(iter(values.values()))
    return canonical_json(values)


def spectrum_variant_id(field, prefix):
    if field is None:
        return f"{prefix}-"
    match = clean_converter.base.SPECTRUM_RE.match(field)
    return f"{prefix}{match.group(2)}" if match else f"{prefix}{stable_hash(field)[:8]}"


def spectrum_pair_score(mol, carbon_item, proton_item):
    c_field, c_meta = carbon_item
    h_field, h_meta = proton_item
    score = clean_converter.base.pair_score(c_meta, h_meta)
    score += 0.05 * clean_converter._unique_signal_count(mol, c_field)
    score += 0.05 * clean_converter._unique_signal_count(mol, h_field)
    return score


def selected_spectrum_variants(mol, nuclei):
    carbon = [
        (field, clean_converter.base.spectrum_meta(mol, field))
        for _, field in clean_converter.base.spectrum_fields(mol, "13C")
    ]
    proton = [
        (field, clean_converter.base.spectrum_meta(mol, field))
        for _, field in clean_converter.base.spectrum_fields(mol, "1H")
    ]

    if nuclei == "carbon":
        raw_variants = [(item, None) for item in carbon]
    elif nuclei == "proton":
        raw_variants = [(None, item) for item in proton]
    elif carbon and proton:
        raw_variants = []
        covered_proton = set()
        for c_item in carbon:
            h_item = max(
                proton,
                key=lambda item: (
                    spectrum_pair_score(mol, c_item, item),
                    item[0],
                ),
            )
            raw_variants.append((c_item, h_item))
            covered_proton.add(h_item[0])
        for h_item in proton:
            if h_item[0] in covered_proton:
                continue
            c_item = max(
                carbon,
                key=lambda item: (
                    spectrum_pair_score(mol, item, h_item),
                    item[0],
                ),
            )
            raw_variants.append((c_item, h_item))
    elif nuclei == "either":
        raw_variants = [(item, None) for item in carbon] or [(None, item) for item in proton]
    else:
        raw_variants = []

    variants = []
    seen = set()
    for c_item, h_item in raw_variants:
        c_field = c_item[0] if c_item else None
        h_field = h_item[0] if h_item else None
        key = (c_field, h_field)
        if key in seen:
            continue
        seen.add(key)
        variants.append(
            {
                "variant_id": f"{spectrum_variant_id(c_field, 'C')}__{spectrum_variant_id(h_field, 'H')}",
                "carbon_field": c_field,
                "proton_field": h_field,
            }
        )
    return variants


def convert_nmrshift_variant(mol, variant):
    work = Chem.Mol(mol)
    selected = {variant.get("carbon_field"), variant.get("proton_field")}
    selected.discard(None)
    for nucleus in ("13C", "1H"):
        for _, field in clean_converter.base.spectrum_fields(work, nucleus):
            if field not in selected:
                work.ClearProp(field)
    use_cnmr = variant.get("carbon_field") is not None
    use_hnmr = variant.get("proton_field") is not None
    if not use_cnmr and not use_hnmr:
        return None, "no_selected_spectrum"
    return clean_converter.molecule_record(
        work, SimpleNamespace(use_cnmr=use_cnmr, use_hnmr=use_hnmr)
    )


def invalid_nmrshift_index_row(mol_index, sdf_path, source_record_id, suffix, reason, metadata=None):
    payload = {"sdf_molecule_index": mol_index, "reason": reason}
    if metadata:
        payload.update(metadata)
    return {
        "source": "nmrshiftdb2",
        "source_split": "",
        "source_path": str(sdf_path),
        "source_locator": f"{mol_index}|{suffix}",
        "source_record_id": f"{source_record_id}:{suffix}",
        "keep_record": 0,
        "exclusion_reason": reason,
        "metadata_json": canonical_json(payload),
    }


def nmrshift_index_rows(mol, mol_index, sdf_path, nuclei):
    default_record_id = f"sdf:{mol_index}"
    if mol is None:
        return [
            invalid_nmrshift_index_row(
                mol_index, sdf_path, default_record_id, "invalid", "invalid_sdf_molecule"
            )
        ]

    source_record_id = default_record_id
    if mol.HasProp("nmrshiftdb2 ID"):
        source_record_id = mol.GetProp("nmrshiftdb2 ID") or default_record_id
    variants = selected_spectrum_variants(mol, nuclei)
    if not variants:
        return [
            invalid_nmrshift_index_row(
                mol_index,
                sdf_path,
                source_record_id,
                "no_variant",
                f"no_spectrum_variant_for_nuclei:{nuclei}",
            )
        ]

    rows = []
    for variant in variants:
        variant_id = variant["variant_id"]
        variant_metadata = {
            "sdf_molecule_index": mol_index,
            "variant_id": variant_id,
            "carbon_source_field": variant.get("carbon_field"),
            "proton_source_field": variant.get("proton_field"),
        }
        try:
            record, reason = convert_nmrshift_variant(mol, variant)
        except Exception as exc:
            rows.append(
                invalid_nmrshift_index_row(
                    mol_index,
                    sdf_path,
                    source_record_id,
                    variant_id,
                    f"conversion_exception:{type(exc).__name__}:{exc}",
                    variant_metadata,
                )
            )
            continue
        if record is None:
            rows.append(
                invalid_nmrshift_index_row(
                    mol_index,
                    sdf_path,
                    source_record_id,
                    variant_id,
                    f"conversion_failed:{reason}",
                    variant_metadata,
                )
            )
            continue

        try:
            exact, connectivity, scaffold, canonical = structure_keys_from_smiles(
                record["smiles_canonical"]
            )
        except Exception as exc:
            variant_metadata["invalid_canonical_smiles"] = record.get("smiles_canonical")
            rows.append(
                invalid_nmrshift_index_row(
                    mol_index,
                    sdf_path,
                    source_record_id,
                    variant_id,
                    f"structure_key_failed:{type(exc).__name__}:{exc}",
                    variant_metadata,
                )
            )
            continue
        c_hash = branch_fingerprint(record, "carbon")
        h_hash = branch_fingerprint(record, "proton")
        c_peaks = len((((record.get("spectra") or {}).get("carbon") or {}).get("peaks") or []))
        h_peaks = len((((record.get("spectra") or {}).get("proton") or {}).get("peaks") or []))
        solvent = condition_value(record, "solvent")
        frequency = condition_value(record, "field_strength_mhz")
        quality_score = c_peaks + h_peaks + 20 * int(bool(c_peaks and h_peaks))
        quality_score += 2 * int(bool(solvent)) + 2 * int(bool(frequency))
        variant_metadata.update(
            {
                "schema_version": record.get("schema_version"),
                "carbon_peak_count": c_peaks,
                "proton_peak_count": h_peaks,
                "pairing_rule": "cover_every_spectrum_field_with_best_metadata_compatible_cross_nucleus_match",
            }
        )
        rows.append(
            {
                "source": "nmrshiftdb2",
                "source_split": "",
                "source_path": str(sdf_path),
                "source_locator": f"{mol_index}|{variant_id}",
                "source_record_id": f"{source_record_id}:{variant_id}",
                "exact_structure_key": exact,
                "connectivity_key": connectivity,
                "scaffold_key": scaffold,
                "canonical_smiles": canonical,
                "spectrum_hash": stable_hash({"13C": c_hash, "1H": h_hash}),
                "carbon_spectrum_hash": c_hash,
                "proton_spectrum_hash": h_hash,
                "solvent": solvent,
                "frequency": frequency,
                "has_carbon": int(bool(c_peaks)),
                "has_proton": int(bool(h_peaks)),
                "has_assignment": 1,
                "quality_score": quality_score,
                "keep_record": 1,
                "metadata_json": canonical_json(variant_metadata),
            }
        )
    return rows

def first_present(record, names, default=None):
    for name in names:
        if name in record and record[name] is not None:
            return record[name]
    return default


def nmrexp_index_row(record, split, source_path, key):
    locator = base64.b64encode(key).decode("ascii")
    source_record_id = f"{split}:{key.decode('utf-8', errors='backslashreplace')}"
    base = {
        "source": "nmrexp",
        "source_split": split,
        "source_path": str(source_path),
        "source_locator": locator,
        "source_record_id": source_record_id,
        "metadata_json": canonical_json({"lmdb_key_b64": locator}),
    }
    try:
        smiles = first_present(record, ("smiles_canonical", "canonical_smiles", "smiles"))
        exact, connectivity, scaffold, canonical = structure_keys_from_smiles(smiles)
    except Exception as exc:
        base.update(keep_record=0, exclusion_reason=f"invalid_nmrexp_structure:{exc}")
        return base
    c_peaks = first_present(record, ("c_nmr_peaks", "cnmr_peaks", "carbon_peaks"), []) or []
    h_peaks = first_present(record, ("h_nmr_peaks", "hnmr_peaks", "proton_peaks"), []) or []
    c_hash = stable_hash(c_peaks) if c_peaks else None
    h_hash = stable_hash(h_peaks) if h_peaks else None
    solvent = first_present(record, ("solvent", "Solvent"))
    frequency = first_present(record, ("frequency", "field_strength_mhz", "Field Strength [MHz]"))
    base.update(
        exact_structure_key=exact,
        connectivity_key=connectivity,
        scaffold_key=scaffold,
        canonical_smiles=canonical,
        spectrum_hash=stable_hash({"13C": c_hash, "1H": h_hash}),
        carbon_spectrum_hash=c_hash,
        proton_spectrum_hash=h_hash,
        solvent=str(solvent) if solvent is not None else None,
        frequency=str(frequency) if frequency is not None else None,
        has_carbon=int(bool(c_peaks)),
        has_proton=int(bool(h_peaks)),
        has_assignment=int(bool(record.get("alignment_supervision") or record.get("assignment_matrix"))),
        quality_score=len(c_peaks) + len(h_peaks),
        keep_record=1,
    )
    return base


SOURCE_COLUMNS = (
    "source",
    "source_split",
    "source_path",
    "source_locator",
    "source_record_id",
    "exact_structure_key",
    "connectivity_key",
    "scaffold_key",
    "canonical_smiles",
    "spectrum_hash",
    "carbon_spectrum_hash",
    "proton_spectrum_hash",
    "solvent",
    "frequency",
    "has_carbon",
    "has_proton",
    "has_assignment",
    "quality_score",
    "duplicate_of",
    "keep_record",
    "final_split",
    "exclusion_reason",
    "metadata_json",
)


def insert_source_rows(conn, rows):
    placeholders = ",".join("?" for _ in SOURCE_COLUMNS)
    columns = ",".join(SOURCE_COLUMNS)
    updates = ",".join(
        f"{column}=excluded.{column}"
        for column in SOURCE_COLUMNS
        if column not in {"source", "source_path", "source_locator"}
    )
    sql = (
        f"INSERT INTO source_records({columns}) VALUES ({placeholders}) "
        f"ON CONFLICT(source, source_path, source_locator) DO UPDATE SET {updates}"
    )
    defaults = {
        "source_split": "",
        "has_carbon": 0,
        "has_proton": 0,
        "has_assignment": 0,
        "quality_score": 0,
        "keep_record": 1,
        "metadata_json": "{}",
    }
    values = []
    for row in rows:
        values.append(
            tuple(row[column] if column in row else defaults.get(column) for column in SOURCE_COLUMNS)
        )
    conn.executemany(sql, values)


def maybe_test_pause(args, run_processed):
    global STOP_REQUESTED, STOP_SIGNAL
    if args.max_records_this_run and run_processed >= args.max_records_this_run:
        STOP_REQUESTED = True
        STOP_SIGNAL = "TEST_LIMIT"
        log(f"Reached --max-records-this-run={args.max_records_this_run}; simulating a graceful pause")


def index_nmrshift(conn, args, cfg_hash):
    global STOP_REQUESTED
    stage, source, split = "index", "nmrshiftdb2", "all"
    progress = get_progress(conn, stage, source, split)
    check_progress_config(progress, cfg_hash)
    if progress and progress["status"] == "completed":
        log("NMRShiftDB2 index already completed; skipping")
        return True

    sdf_path = Path(args.sdf_path).expanduser().resolve()
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    total = len(supplier)
    start = int(progress["last_source_locator"]) + 1 if progress and progress["last_source_locator"] else 0
    processed = progress["processed_count"] if progress else 0
    kept = progress["kept_count"] if progress else 0
    if progress:
        log(
            f"Existing progress database found. Stage: index, status: {progress['status']}. "
            f"Last committed NMRShiftDB2 molecule: {start - 1}; resuming from molecule: {start}"
        )
    else:
        log(f"Starting NMRShiftDB2 index: {total} SDF molecules")
    with conn:
        save_progress(conn, stage, source, split, cfg_hash, "running", str(start - 1) if start else None, -1, processed, kept, total)

    rows = []
    batch_molecules = 0
    batch_kept = 0
    run_processed = 0
    last_index = start - 1
    try:
        for mol_index in range(start, total):
            variant_rows = nmrshift_index_rows(
                supplier[mol_index], mol_index, sdf_path, args.nuclei
            )
            rows.extend(variant_rows)
            batch_molecules += 1
            batch_kept += sum(int(bool(row.get("keep_record"))) for row in variant_rows)
            run_processed += 1
            last_index = mol_index
            maybe_test_pause(args, run_processed)
            if batch_molecules >= args.commit_interval or STOP_REQUESTED:
                with conn:
                    insert_source_rows(conn, rows)
                    processed += batch_molecules
                    kept += batch_kept
                    save_progress(
                        conn,
                        stage,
                        source,
                        split,
                        cfg_hash,
                        "paused" if STOP_REQUESTED else "running",
                        str(last_index),
                        -1,
                        processed,
                        kept,
                        total,
                    )
                log(
                    f"NMRShiftDB2 index committed: molecules={processed}/{total}, "
                    f"usable spectrum variants={kept}"
                )
                rows, batch_molecules, batch_kept = [], 0, 0
            if STOP_REQUESTED:
                return False
        if batch_molecules:
            with conn:
                insert_source_rows(conn, rows)
                processed += batch_molecules
                kept += batch_kept
        with conn:
            save_progress(
                conn,
                stage,
                source,
                split,
                cfg_hash,
                "completed",
                str(last_index),
                -1,
                processed,
                kept,
                total,
            )
        log(
            f"NMRShiftDB2 index completed: molecules={processed}, "
            f"usable spectrum variants={kept}"
        )
        return True
    except Exception as exc:
        mark_failed(conn, stage, source, split, cfg_hash, exc)
        raise

def index_nmrexp_split(conn, args, cfg_hash, split):
    global STOP_REQUESTED
    stage, source = "index", "nmrexp"
    progress = get_progress(conn, stage, source, split)
    check_progress_config(progress, cfg_hash)
    if progress and progress["status"] == "completed":
        log(f"NMRexp {split} index already completed; skipping")
        return True

    source_path = find_lmdb_path(args.nmrexp_dir, split)
    env = lmdb.open(str(source_path), subdir=False, readonly=True, lock=False, readahead=False, max_readers=256)
    total = env.stat()["entries"]
    processed = progress["processed_count"] if progress else 0
    kept = progress["kept_count"] if progress else 0
    last_key = base64.b64decode(progress["last_source_locator"]) if progress and progress["last_source_locator"] else None
    if progress:
        shown = last_key.decode("utf-8", errors="backslashreplace") if last_key else "<none>"
        log(
            f"Existing progress database found. Stage: index/nmrexp/{split}, status: {progress['status']}. "
            f"Last committed LMDB key: {shown}"
        )
    else:
        log(f"Starting NMRexp {split} index: {total} LMDB records")
    with conn:
        save_progress(
            conn,
            stage,
            source,
            split,
            cfg_hash,
            "running",
            progress["last_source_locator"] if progress else None,
            -1,
            processed,
            kept,
            total,
        )

    rows = []
    batch_kept = 0
    run_processed = 0
    current_key = last_key
    try:
        with env.begin(write=False) as txn:
            cursor = txn.cursor()
            if last_key is None:
                available = cursor.first()
            else:
                available = cursor.set_range(last_key)
                if available and cursor.key() == last_key:
                    available = cursor.next()
            while available:
                current_key = bytes(cursor.key())
                try:
                    record = pickle.loads(cursor.value())
                    row = nmrexp_index_row(record, split, source_path, current_key)
                except Exception as exc:
                    locator = base64.b64encode(current_key).decode("ascii")
                    row = {
                        "source": "nmrexp",
                        "source_split": split,
                        "source_path": str(source_path),
                        "source_locator": locator,
                        "source_record_id": f"{split}:{current_key!r}",
                        "keep_record": 0,
                        "exclusion_reason": f"unpickle_or_parse_failed:{exc}",
                        "metadata_json": canonical_json({"lmdb_key_b64": locator}),
                    }
                rows.append(row)
                batch_kept += int(bool(row.get("keep_record")))
                run_processed += 1
                maybe_test_pause(args, run_processed)
                if len(rows) >= args.commit_interval or STOP_REQUESTED:
                    locator = base64.b64encode(current_key).decode("ascii")
                    with conn:
                        insert_source_rows(conn, rows)
                        processed += len(rows)
                        kept += batch_kept
                        save_progress(
                            conn,
                            stage,
                            source,
                            split,
                            cfg_hash,
                            "paused" if STOP_REQUESTED else "running",
                            locator,
                            -1,
                            processed,
                            kept,
                            total,
                        )
                    log(f"NMRexp {split} index committed: {processed}/{total}, usable={kept}")
                    rows, batch_kept = [], 0
                if STOP_REQUESTED:
                    break
                available = cursor.next()
        if STOP_REQUESTED:
            env.close()
            return False
        if rows:
            with conn:
                insert_source_rows(conn, rows)
                processed += len(rows)
                kept += batch_kept
        locator = base64.b64encode(current_key).decode("ascii") if current_key is not None else None
        with conn:
            save_progress(conn, stage, source, split, cfg_hash, "completed", locator, -1, processed, kept, total)
        log(f"NMRexp {split} index completed: processed={processed}, usable={kept}")
        env.close()
        return True
    except Exception as exc:
        env.close()
        mark_failed(conn, stage, source, split, cfg_hash, exc)
        raise


def run_index(conn, args, cfg_hash):
    if not index_nmrshift(conn, args, cfg_hash):
        return False
    for split in DEFAULT_SPLITS:
        if STOP_REQUESTED:
            return False
        if not index_nmrexp_split(conn, args, cfg_hash, split):
            return False
    with conn:
        conn.execute("INSERT OR REPLACE INTO pipeline_meta(key, value) VALUES ('semantic_config', ?)", (canonical_json(semantic_config(args)),))
    log("Stage index completed for NMRShiftDB2 and all NMRexp splits")
    return True


def require_index_complete(conn, cfg_hash):
    expected = [("nmrshiftdb2", "all")] + [("nmrexp", split) for split in DEFAULT_SPLITS]
    incomplete = []
    for source, split in expected:
        row = get_progress(conn, "index", source, split)
        check_progress_config(row, cfg_hash)
        if row is None or row["status"] != "completed":
            incomplete.append(f"{source}/{split}")
    if incomplete:
        raise RuntimeError(f"Index stage is not complete for: {', '.join(incomplete)}")


def deduplicate_nmrshift(conn):
    rows = conn.execute(
        """
        SELECT id, exact_structure_key, carbon_spectrum_hash, proton_spectrum_hash,
               COALESCE(solvent, '') AS solvent_key,
               COALESCE(frequency, '') AS frequency_key,
               quality_score
        FROM source_records
        WHERE source='nmrshiftdb2' AND exact_structure_key IS NOT NULL
        ORDER BY exact_structure_key, COALESCE(carbon_spectrum_hash, ''),
                 COALESCE(proton_spectrum_hash, ''), COALESCE(solvent, ''),
                 COALESCE(frequency, ''), quality_score DESC, id ASC
        """
    )
    duplicate_map = {}
    current_key = None
    winner = None
    for row in rows:
        key = (
            row["exact_structure_key"],
            row["carbon_spectrum_hash"] or "",
            row["proton_spectrum_hash"] or "",
            row["solvent_key"],
            row["frequency_key"],
        )
        if key != current_key:
            current_key, winner = key, row["id"]
        else:
            duplicate_map[row["id"]] = winner
    return duplicate_map


def deterministic_tie(seed, key):
    return hashlib.sha256(f"{seed}:{key}".encode("utf-8")).hexdigest()


def build_leakage_components(rows):
    """Merge records connected by scaffold, connectivity, or exact structure."""
    parent = {row["id"]: row["id"] for row in rows}
    first_for_token = {}

    def find(record_id):
        while parent[record_id] != record_id:
            parent[record_id] = parent[parent[record_id]]
            record_id = parent[record_id]
        return record_id

    def union(left, right):
        left_root = find(left)
        right_root = find(right)
        if left_root == right_root:
            return
        if left_root > right_root:
            left_root, right_root = right_root, left_root
        parent[right_root] = left_root

    for row in rows:
        tokens = (
            f"scaffold:{row['scaffold_key']}" if row["scaffold_key"] else None,
            f"connectivity:{row['connectivity_key']}" if row["connectivity_key"] else None,
            f"exact:{row['exact_structure_key']}" if row["exact_structure_key"] else None,
        )
        for token in filter(None, tokens):
            if token in first_for_token:
                union(row["id"], first_for_token[token])
            else:
                first_for_token[token] = row["id"]

    grouped = defaultdict(
        lambda: {
            "ids": [],
            "tokens": set(),
            "records": 0,
            "carbon": 0,
            "proton": 0,
            "both": 0,
        }
    )
    for row in rows:
        item = grouped[find(row["id"])]
        item["ids"].append(row["id"])
        item["tokens"].update(
            token
            for token in (
                row["scaffold_key"],
                row["connectivity_key"],
                row["exact_structure_key"],
            )
            if token
        )
        item["records"] += 1
        item["carbon"] += row["has_carbon"]
        item["proton"] += row["has_proton"]
        item["both"] += int(bool(row["has_carbon"] and row["has_proton"]))

    components = []
    for item in grouped.values():
        tokens = sorted(item.pop("tokens"))
        components.append(dict(item, key=f"COMPONENT:{stable_hash(tokens)}"))
    return components

def assign_scaffold_groups(groups, train_ratio, valid_ratio, seed):
    ratios = {"train": train_ratio, "valid": valid_ratio, "test": 1.0 - train_ratio - valid_ratio}
    totals = {name: sum(group[name] for group in groups) for name in ("records", "carbon", "proton", "both")}
    targets = {
        split: {name: totals[name] * ratio for name in totals}
        for split, ratio in ratios.items()
    }
    current = {split: {name: 0 for name in totals} for split in DEFAULT_SPLITS}
    assignments = {}
    ordered = sorted(groups, key=lambda group: (-group["records"], deterministic_tie(seed, group["key"])))
    for group in ordered:
        best_split = None
        best_score = None
        for split in DEFAULT_SPLITS:
            score = 0.0
            for name, total in totals.items():
                weight = 1.0 if name == "records" else 0.15
                candidate = current[split][name] + group[name]
                score += weight * ((candidate - targets[split][name]) / max(total, 1)) ** 2
                for other in DEFAULT_SPLITS:
                    if other != split:
                        score += weight * ((current[other][name] - targets[other][name]) / max(total, 1)) ** 2
            tie = (DEFAULT_SPLITS.index(split),)
            candidate_score = (score, tie)
            if best_score is None or candidate_score < best_score:
                best_score, best_split = candidate_score, split
        assignments[group["key"]] = best_split
        for name in totals:
            current[best_split][name] += group[name]
    return assignments, current, targets


def run_split(conn, args, cfg_hash):
    require_index_complete(conn, cfg_hash)
    progress = get_progress(conn, "split", "all", "all")
    check_progress_config(progress, cfg_hash)
    if progress and progress["status"] == "completed":
        log("Stage split already completed; nothing to do")
        return True
    with conn:
        save_progress(conn, "split", "all", "all", cfg_hash, "running")

    try:
        duplicate_map = deduplicate_nmrshift(conn)
        winner_rows = conn.execute(
            """
            SELECT id, scaffold_key, connectivity_key, exact_structure_key,
                   has_carbon, has_proton
            FROM source_records
            WHERE source='nmrshiftdb2' AND exact_structure_key IS NOT NULL
            ORDER BY id
            """
        ).fetchall()
        winner_rows = [row for row in winner_rows if row["id"] not in duplicate_map]
        groups = build_leakage_components(winner_rows)
        largest_component = max((group["records"] for group in groups), default=0)
        log(
            f"Built {len(groups)} scaffold/connectivity components; "
            f"largest component={largest_component} records"
        )

        scaffold_splits, counts, targets = assign_scaffold_groups(
            groups, args.train_ratio, args.valid_ratio, args.seed
        )

        nmrshift_assignments = {}
        exact_to_splits = defaultdict(set)
        connectivity_to_splits = defaultdict(set)
        winner_split = {}
        for group in groups:
            final_split = scaffold_splits[group["key"]]
            for record_id in group["ids"]:
                nmrshift_assignments[record_id] = (final_split, None)
                winner_split[record_id] = final_split
        for row in winner_rows:
            final_split = winner_split[row["id"]]
            exact_to_splits[row["exact_structure_key"]].add(final_split)
            connectivity_to_splits[row["connectivity_key"]].add(final_split)
        conflicts = [key for key, values in connectivity_to_splits.items() if len(values) > 1]
        if conflicts:
            raise RuntimeError(f"Connectivity groups span multiple NMRShiftDB2 splits: {conflicts[:5]}")

        total_valid = conn.execute(
            "SELECT COUNT(*) FROM source_records WHERE exact_structure_key IS NOT NULL"
        ).fetchone()[0]
        processed = 0
        with conn:
            conn.execute("DROP TABLE IF EXISTS split_assignments_new")
            conn.execute(
                """
                CREATE TABLE split_assignments_new (
                    record_id INTEGER PRIMARY KEY REFERENCES source_records(id),
                    final_split TEXT NOT NULL,
                    exclusion_reason TEXT
                )
                """
            )
            conn.executemany(
                "INSERT INTO split_assignments_new(record_id, final_split, exclusion_reason) VALUES (?, ?, ?)",
                [(record_id, split, reason) for record_id, (split, reason) in nmrshift_assignments.items()],
            )
            processed += len(nmrshift_assignments)
            read_cursor = conn.execute(
                """
                SELECT id, source_split, exact_structure_key, connectivity_key
                FROM source_records
                WHERE source='nmrexp' AND exact_structure_key IS NOT NULL
                ORDER BY id
                """
            )
            write_cursor = conn.cursor()
            batch = []
            for row in read_cursor:
                exact_splits = exact_to_splits.get(row["exact_structure_key"], set())
                connectivity_splits = connectivity_to_splits.get(row["connectivity_key"], set())
                matched = exact_splits | connectivity_splits
                heldout = matched & {"valid", "test"}
                if len(heldout) > 1:
                    raise RuntimeError(
                        f"NMRexp record {row['id']} connects to both NMRShiftDB2 valid and test"
                    )
                if heldout:
                    final_split = next(iter(heldout))
                    match_level = "exact" if exact_splits & heldout else "connectivity"
                    reason = f"cross_source_{match_level}_heldout_{final_split}"
                else:
                    final_split = row["source_split"]
                    reason = "cross_source_match_train" if "train" in matched else None
                batch.append((row["id"], final_split, reason))
                processed += 1
                if len(batch) >= args.commit_interval:
                    write_cursor.executemany(
                        "INSERT INTO split_assignments_new(record_id, final_split, exclusion_reason) VALUES (?, ?, ?)",
                        batch,
                    )
                    batch = []
                if STOP_REQUESTED:
                    raise InterruptedError("split interrupted by signal")
            if batch:
                write_cursor.executemany(
                    "INSERT INTO split_assignments_new(record_id, final_split, exclusion_reason) VALUES (?, ?, ?)",
                    batch,
                )

            conn.execute("DROP TABLE IF EXISTS split_assignments")
            conn.execute("ALTER TABLE split_assignments_new RENAME TO split_assignments")
            conn.execute(
                """
                UPDATE source_records
                SET duplicate_of=NULL,
                    final_split=NULL,
                    keep_record=CASE WHEN exact_structure_key IS NULL THEN 0 ELSE 1 END,
                    exclusion_reason=CASE
                        WHEN exact_structure_key IS NULL THEN exclusion_reason
                        ELSE NULL
                    END
                """
            )
            conn.executemany(
                """
                UPDATE source_records
                SET duplicate_of=?, keep_record=0,
                    exclusion_reason='exact_duplicate_within_nmrshiftdb2'
                WHERE id=?
                """,
                [(winner, duplicate) for duplicate, winner in duplicate_map.items()],
            )
            conn.execute(
                """
                UPDATE source_records
                SET final_split=(SELECT final_split FROM split_assignments WHERE record_id=source_records.id),
                    exclusion_reason=COALESCE(
                        (SELECT exclusion_reason FROM split_assignments WHERE record_id=source_records.id),
                        source_records.exclusion_reason
                    )
                WHERE id IN (SELECT record_id FROM split_assignments)
                """
            )
            save_progress(
                conn,
                "split",
                "all",
                "all",
                cfg_hash,
                "completed",
                "atomic_split_table",
                -1,
                total_valid,
                processed,
                total_valid,
            )
            conn.execute(
                "INSERT OR REPLACE INTO pipeline_meta(key, value) VALUES ('split_summary', ?)",
                (
                    canonical_json(
                        {
                            "nmrshiftdb2_kept": len(nmrshift_assignments),
                            "nmrshiftdb2_exact_duplicates": len(duplicate_map),
                            "leakage_component_count": len(groups),
                            "largest_component_records": largest_component,
                            "actual_counts": counts,
                            "targets": targets,
                        }
                    ),
                ),
            )
        log(
            "Stage split completed atomically: "
            f"NMRShiftDB2 kept={len(nmrshift_assignments)}, exact duplicates={len(duplicate_map)}"
        )
        for split in DEFAULT_SPLITS:
            log(f"  NMRShiftDB2 {split}: {counts[split]['records']} records")
        return True
    except InterruptedError:
        with conn:
            save_progress(conn, "split", "all", "all", cfg_hash, "paused")
        log("Split calculation paused; the next run will recompute only the split stage")
        return False
    except Exception as exc:
        mark_failed(conn, "split", "all", "all", cfg_hash, exc)
        raise


def require_split_complete(conn, cfg_hash):
    row = get_progress(conn, "split", "all", "all")
    check_progress_config(row, cfg_hash)
    if row is None or row["status"] != "completed":
        raise RuntimeError("Split stage is not completed")


def output_paths(output_dir, source, split):
    root = Path(output_dir).expanduser().resolve() / source / "lmdb_dataset"
    root.mkdir(parents=True, exist_ok=True)
    final = root / f"{split}.lmdb"
    partial = root / f"{split}.lmdb.partial"
    return partial, final


def lmdb_entry_count(path):
    env = lmdb.open(str(path), subdir=False, readonly=True, lock=False, readahead=False)
    count = env.stat()["entries"]
    env.close()
    return count


def read_nmrexp_record(env_cache, source_path, locator):
    if source_path not in env_cache:
        env_cache[source_path] = lmdb.open(
            source_path, subdir=False, readonly=True, lock=False, readahead=False, max_readers=256
        )
    key = base64.b64decode(locator)
    with env_cache[source_path].begin(write=False) as txn:
        value = txn.get(key)
    if value is None:
        raise KeyError(f"NMRexp source key disappeared: {source_path}:{key!r}")
    return pickle.loads(value)


def attach_pipeline_metadata(record, row):
    output = copy.deepcopy(record)
    output["smiles_canonical"] = row["canonical_smiles"]
    output["atompeak_pipeline"] = {
        "schema_version": SCHEMA_VERSION,
        "source": row["source"],
        "source_record_id": row["source_record_id"],
        "source_locator": row["source_locator"],
        "final_split": row["final_split"],
        "exact_structure_key": row["exact_structure_key"],
        "connectivity_key": row["connectivity_key"],
        "scaffold_key": row["scaffold_key"],
        "spectrum_hash": row["spectrum_hash"],
        "cross_source_policy": row["exclusion_reason"],
    }
    return output


def recover_completed_materialization(conn, progress, cfg_hash, source, split, partial, final, expected):
    if not final.exists() or partial.exists():
        return False
    count = lmdb_entry_count(final)
    if count != expected:
        raise RuntimeError(
            f"Completed-looking file {final} has {count} entries, expected {expected}; "
            "use --reset-stage materialize --overwrite-completed after inspection"
        )
    with conn:
        save_progress(
            conn,
            "materialize",
            source,
            split,
            cfg_hash,
            "completed",
            progress["last_source_locator"] if progress else None,
            progress["last_output_index"] if progress else -1,
            expected,
            expected,
            expected,
        )
    log(f"Recovered completed materialization metadata from {final}")
    return True


def materialize_one(conn, args, cfg_hash, source, split):
    global STOP_REQUESTED
    progress = get_progress(conn, "materialize", source, split)
    check_progress_config(progress, cfg_hash)
    expected = conn.execute(
        "SELECT COUNT(*) FROM source_records WHERE source=? AND keep_record=1 AND final_split=?",
        (source, split),
    ).fetchone()[0]
    partial, final = output_paths(args.output_dir, source, split)
    if progress and progress["status"] == "completed" and final.exists():
        log(f"Materialize {source}/{split} already completed; skipping")
        return True
    if recover_completed_materialization(conn, progress, cfg_hash, source, split, partial, final, expected):
        return True
    if final.exists() and not args.overwrite_completed:
        raise RuntimeError(
            f"Refusing to overwrite completed LMDB {final}. Use --reset-stage materialize "
            "together with --overwrite-completed."
        )
    if progress and progress["processed_count"] and not partial.exists():
        raise RuntimeError(
            f"Saved materialize progress exists but {partial} is missing. "
            "Use --reset-stage materialize after inspection."
        )

    last_id = progress["last_output_index"] if progress else -1
    processed = progress["processed_count"] if progress else 0
    kept = progress["kept_count"] if progress else 0
    with conn:
        save_progress(
            conn,
            "materialize",
            source,
            split,
            cfg_hash,
            "running",
            progress["last_source_locator"] if progress else None,
            last_id,
            processed,
            kept,
            expected,
        )
    if progress:
        log(
            f"Resuming materialize {source}/{split}: committed={processed}/{expected}, "
            f"last record id={last_id}"
        )
    else:
        log(f"Starting materialize {source}/{split}: expected={expected}")

    env = lmdb.open(
        str(partial),
        subdir=False,
        map_size=int(args.map_size_gb * 1024**3),
        lock=True,
        sync=True,
        metasync=True,
        map_async=False,
    )
    sdf_supplier = None
    if source == "nmrshiftdb2":
        sdf_supplier = Chem.SDMolSupplier(str(Path(args.sdf_path).expanduser().resolve()), removeHs=False, sanitize=False)
    source_envs = {}
    txn = env.begin(write=True)
    batch_count = 0
    run_processed = 0
    current_locator = progress["last_source_locator"] if progress else None
    current_id = last_id

    try:
        rows = conn.execute(
            """
            SELECT * FROM source_records
            WHERE source=? AND keep_record=1 AND final_split=? AND id>?
            ORDER BY id
            """,
            (source, split, last_id),
        )
        for row in rows:
            current_id = row["id"]
            current_locator = row["source_locator"]
            if source == "nmrshiftdb2":
                molecule_index = int(current_locator.split("|", 1)[0])
                mol = sdf_supplier[molecule_index]
                if mol is None:
                    raise RuntimeError(
                        f"SDF molecule {molecule_index} became unreadable during materialize"
                    )
                metadata = json.loads(row["metadata_json"])
                variant = {
                    "carbon_field": metadata.get("carbon_source_field"),
                    "proton_field": metadata.get("proton_source_field"),
                }
                record, reason = convert_nmrshift_variant(mol, variant)
                if record is None:
                    raise RuntimeError(
                        f"SDF spectrum variant {current_locator} no longer converts: {reason}"
                    )
                current_hash = stable_hash(
                    {
                        "13C": branch_fingerprint(record, "carbon"),
                        "1H": branch_fingerprint(record, "proton"),
                    }
                )
                if current_hash != row["spectrum_hash"]:
                    raise RuntimeError(
                        f"Spectrum fingerprint changed for SDF variant {current_locator}"
                    )
            else:
                record = read_nmrexp_record(source_envs, row["source_path"], current_locator)
            output = attach_pipeline_metadata(record, row)
            key = f"{current_id:012d}".encode("ascii")
            txn.put(key, pickle.dumps(output, protocol=pickle.HIGHEST_PROTOCOL), overwrite=True)
            batch_count += 1
            run_processed += 1
            maybe_test_pause(args, run_processed)
            if batch_count >= args.commit_interval or STOP_REQUESTED:
                txn.commit()
                env.sync()
                processed += batch_count
                kept += batch_count
                with conn:
                    save_progress(
                        conn,
                        "materialize",
                        source,
                        split,
                        cfg_hash,
                        "paused" if STOP_REQUESTED else "running",
                        current_locator,
                        current_id,
                        processed,
                        kept,
                        expected,
                    )
                log(f"Materialize {source}/{split} committed: {processed}/{expected}")
                batch_count = 0
                if STOP_REQUESTED:
                    for source_env in source_envs.values():
                        source_env.close()
                    env.close()
                    return False
                txn = env.begin(write=True)
        if batch_count:
            txn.commit()
            env.sync()
            processed += batch_count
            kept += batch_count
        else:
            txn.abort()
        for source_env in source_envs.values():
            source_env.close()
        entries = env.stat()["entries"]
        env.sync()
        env.close()
        if entries != expected or processed != expected:
            raise RuntimeError(
                f"Materialize count mismatch for {source}/{split}: LMDB={entries}, "
                f"progress={processed}, expected={expected}"
            )
        if final.exists():
            if not args.overwrite_completed:
                raise RuntimeError(f"Refusing to replace existing {final}")
            final.unlink()
        os.replace(str(partial), str(final))
        lock_path = Path(str(partial) + "-lock")
        if lock_path.exists():
            lock_path.unlink()
        with conn:
            save_progress(
                conn,
                "materialize",
                source,
                split,
                cfg_hash,
                "completed",
                current_locator,
                current_id,
                processed,
                kept,
                expected,
            )
        log(f"Materialize {source}/{split} completed and finalized: {final}")
        return True
    except Exception as exc:
        try:
            txn.abort()
        except Exception:
            pass
        for source_env in source_envs.values():
            source_env.close()
        env.close()
        mark_failed(conn, "materialize", source, split, cfg_hash, exc)
        raise


def run_materialize(conn, args, cfg_hash):
    require_split_complete(conn, cfg_hash)
    for source in ("nmrshiftdb2", "nmrexp"):
        for split in DEFAULT_SPLITS:
            if STOP_REQUESTED:
                return False
            if not materialize_one(conn, args, cfg_hash, source, split):
                return False
    log("Stage materialize completed for all source/split outputs")
    return True


def materialized_counts(conn, output_dir):
    result = {}
    for source in ("nmrshiftdb2", "nmrexp"):
        result[source] = {}
        for split in DEFAULT_SPLITS:
            expected = conn.execute(
                "SELECT COUNT(*) FROM source_records WHERE source=? AND keep_record=1 AND final_split=?",
                (source, split),
            ).fetchone()[0]
            _, final = output_paths(output_dir, source, split)
            result[source][split] = {
                "expected": expected,
                "actual": lmdb_entry_count(final) if final.exists() else None,
                "path": str(final),
            }
    return result


def run_audit(conn, args, cfg_hash):
    require_split_complete(conn, cfg_hash)
    progress = get_progress(conn, "audit", "all", "all")
    check_progress_config(progress, cfg_hash)
    if progress and progress["status"] == "completed":
        log("Stage audit already completed; nothing to do")
        return True
    with conn:
        save_progress(conn, "audit", "all", "all", cfg_hash, "running")
    checks = {}
    checks["kept_nmrshiftdb2_duplicate_groups"] = conn.execute(
        """
        SELECT COUNT(*) FROM (
            SELECT exact_structure_key, carbon_spectrum_hash, proton_spectrum_hash,
                   COALESCE(solvent, ''), COALESCE(frequency, ''), COUNT(*) AS n
            FROM source_records
            WHERE source='nmrshiftdb2' AND keep_record=1
            GROUP BY exact_structure_key, carbon_spectrum_hash, proton_spectrum_hash,
                     COALESCE(solvent, ''), COALESCE(frequency, '')
            HAVING n > 1
        )
        """
    ).fetchone()[0]
    checks["nmrshiftdb2_connectivity_cross_split"] = conn.execute(
        """
        SELECT COUNT(*) FROM (
            SELECT connectivity_key, COUNT(DISTINCT final_split) AS n
            FROM source_records
            WHERE source='nmrshiftdb2' AND keep_record=1
            GROUP BY connectivity_key HAVING n > 1
        )
        """
    ).fetchone()[0]
    checks["nmrexp_train_exact_overlap_nmrshift_heldout"] = conn.execute(
        """
        SELECT COUNT(DISTINCT e.exact_structure_key)
        FROM source_records e JOIN source_records s
          ON e.exact_structure_key=s.exact_structure_key
        WHERE e.source='nmrexp' AND e.keep_record=1 AND e.final_split='train'
          AND s.source='nmrshiftdb2' AND s.keep_record=1 AND s.final_split IN ('valid','test')
        """
    ).fetchone()[0]
    checks["nmrexp_train_connectivity_overlap_nmrshift_heldout"] = conn.execute(
        """
        SELECT COUNT(DISTINCT e.connectivity_key)
        FROM source_records e JOIN source_records s
          ON e.connectivity_key=s.connectivity_key
        WHERE e.source='nmrexp' AND e.keep_record=1 AND e.final_split='train'
          AND s.source='nmrshiftdb2' AND s.keep_record=1 AND s.final_split IN ('valid','test')
        """
    ).fetchone()[0]
    counts = {
        row["source"]: {
            split: conn.execute(
                "SELECT COUNT(*) FROM source_records WHERE source=? AND keep_record=1 AND final_split=?",
                (row["source"], split),
            ).fetchone()[0]
            for split in DEFAULT_SPLITS
        }
        for row in conn.execute("SELECT DISTINCT source FROM source_records")
    }
    outputs = materialized_counts(conn, args.output_dir)
    for source, source_counts in outputs.items():
        for split, values in source_counts.items():
            if values["actual"] is not None:
                checks[f"lmdb_count_mismatch_{source}_{split}"] = int(values["actual"] != values["expected"])
    passed = all(value == 0 for value in checks.values())
    report = {
        "schema_version": SCHEMA_VERSION,
        "passed": passed,
        "checks": checks,
        "kept_counts": counts,
        "materialized": outputs,
        "created_at": utc_now(),
    }
    report_path = Path(args.output_dir).expanduser().resolve() / "audit_report.json"
    temp_path = report_path.with_suffix(".json.partial")
    temp_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    os.replace(str(temp_path), str(report_path))
    with conn:
        save_progress(
            conn,
            "audit",
            "all",
            "all",
            cfg_hash,
            "completed" if passed else "failed",
            str(report_path),
            -1,
            sum(sum(values.values()) for values in counts.values()),
            sum(sum(values.values()) for values in counts.values()),
            sum(sum(values.values()) for values in counts.values()),
        )
    log(json.dumps(report, ensure_ascii=False, indent=2))
    if not passed:
        raise RuntimeError(f"Dataset audit failed; inspect {report_path}")
    return True


def completed_output_files(output_dir):
    output_dir = Path(output_dir).expanduser().resolve()
    return [
        output_dir / source / "lmdb_dataset" / f"{split}.lmdb"
        for source in ("nmrshiftdb2", "nmrexp")
        for split in DEFAULT_SPLITS
    ]


def partial_output_files(output_dir):
    output_dir = Path(output_dir).expanduser().resolve()
    return list(output_dir.glob("*/lmdb_dataset/*.lmdb.partial*"))


def reset_stage(conn, output_dir, stage, overwrite_completed):
    finals = [path for path in completed_output_files(output_dir) if path.exists()]
    if stage in {"index", "split", "materialize"} and finals and not overwrite_completed:
        raise RuntimeError(
            "Completed LMDB files exist. Refusing to reset an upstream stage without "
            "--overwrite-completed. Raw source data are never touched."
        )
    if overwrite_completed:
        for path in finals:
            path.unlink()
            lock_path = Path(str(path) + "-lock")
            if lock_path.exists():
                lock_path.unlink()
    for path in partial_output_files(output_dir):
        path.unlink()
    with conn:
        if stage == "index":
            conn.execute("DELETE FROM split_assignments")
            conn.execute("UPDATE source_records SET duplicate_of=NULL")
            conn.execute("DELETE FROM source_records")
            conn.execute("DELETE FROM sqlite_sequence WHERE name='source_records'")
            conn.execute("DELETE FROM pipeline_progress")
            conn.execute("DELETE FROM pipeline_meta")
        elif stage == "split":
            conn.execute("DELETE FROM split_assignments")
            conn.execute(
                """
                UPDATE source_records SET
                    duplicate_of=NULL,
                    keep_record=CASE WHEN exact_structure_key IS NULL THEN 0 ELSE 1 END,
                    final_split=NULL,
                    exclusion_reason=CASE WHEN exact_structure_key IS NULL THEN exclusion_reason ELSE NULL END
                """
            )
            conn.execute("DELETE FROM pipeline_progress WHERE stage IN ('split','materialize','audit')")
        elif stage == "materialize":
            conn.execute("DELETE FROM pipeline_progress WHERE stage IN ('materialize','audit')")
        elif stage == "audit":
            conn.execute("DELETE FROM pipeline_progress WHERE stage='audit'")
        else:
            raise ValueError(f"Unknown reset stage: {stage}")
    if stage == "audit":
        report = Path(output_dir).expanduser().resolve() / "audit_report.json"
        if report.exists():
            report.unlink()
    log(f"Reset derived stage '{stage}'. Source SDF and NMRexp LMDB files were not modified")


def run_status(output_dir):
    output_dir = Path(output_dir).expanduser().resolve()
    db_path = output_dir / "pipeline.sqlite"
    if not db_path.exists():
        print(f"No progress database found: {db_path}")
        return 0
    conn = connect_db(output_dir)
    rows = conn.execute(
        "SELECT * FROM pipeline_progress ORDER BY CASE stage WHEN 'index' THEN 1 WHEN 'split' THEN 2 WHEN 'materialize' THEN 3 ELSE 4 END, source, split"
    ).fetchall()
    if not rows:
        print(f"Progress database exists but contains no stage records: {db_path}")
    else:
        print(f"Progress database: {db_path}")
        print("stage        source          split   status      processed/total   kept   last locator / output id   updated")
        print("-" * 132)
        for row in rows:
            total = "?" if row["total_count"] is None else str(row["total_count"])
            locator = row["last_source_locator"] or "-"
            if len(locator) > 24:
                locator = locator[:21] + "..."
            marker = f"{locator} / {row['last_output_index']}"
            print(
                f"{row['stage']:<12} {row['source']:<15} {row['split']:<7} {row['status']:<11} "
                f"{row['processed_count']:>9}/{total:<9} {row['kept_count']:>7}   {marker:<30} {row['updated_at']}"
            )
    print("\nIndexed records:")
    for row in conn.execute(
        """
        SELECT source, COALESCE(final_split, '-') AS final_split,
               COUNT(*) AS total, SUM(keep_record) AS kept
        FROM source_records GROUP BY source, COALESCE(final_split, '-')
        ORDER BY source, final_split
        """
    ):
        print(f"  {row['source']:<15} {row['final_split']:<7} total={row['total']:<9} kept={row['kept'] or 0}")
    errors = conn.execute("SELECT key, value FROM pipeline_meta WHERE key LIKE 'last_error:%' ORDER BY key").fetchall()
    if errors:
        print("\nLast recorded errors:")
        for row in errors:
            print(f"  {row['key']}: {row['value']}")
    conn.close()
    return 0


def add_common_arguments(parser):
    parser.add_argument("sdf_path", help="Original NMRShiftDB2 SDF file")
    parser.add_argument("nmrexp_dir", help="NMRexp root or lmdb_dataset directory")
    parser.add_argument("output_dir", help="Derived dataset output directory")
    parser.add_argument("--commit-interval", type=int, default=500)
    parser.add_argument("--map-size-gb", type=float, default=32.0)
    parser.add_argument("--nuclei", choices=("both", "either", "carbon", "proton"), default="both")
    parser.add_argument("--train-ratio", type=float, default=0.8)
    parser.add_argument("--valid-ratio", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--reset-stage", choices=("index", "split", "materialize", "audit"))
    parser.add_argument("--overwrite-completed", action="store_true")
    parser.add_argument(
        "--max-records-this-run",
        type=int,
        default=0,
        help="Testing hook: gracefully pause after N newly processed records",
    )


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="stage", required=True)
    for name in ("index", "split", "materialize", "audit"):
        add_common_arguments(subparsers.add_parser(name))
    status_parser = subparsers.add_parser("status")
    status_parser.add_argument("--output-dir", required=True)
    return parser


def validate_args(args):
    if args.commit_interval <= 0:
        raise ValueError("--commit-interval must be positive")
    if args.map_size_gb <= 0:
        raise ValueError("--map-size-gb must be positive")
    test_ratio = 1.0 - args.train_ratio - args.valid_ratio
    if args.train_ratio <= 0 or args.valid_ratio <= 0 or test_ratio <= 0:
        raise ValueError("train, valid, and test ratios must all be positive")
    for path in (args.sdf_path, args.nmrexp_dir):
        if not Path(path).expanduser().exists():
            raise FileNotFoundError(path)


def main(argv=None):
    global STOP_REQUESTED, STOP_SIGNAL
    STOP_REQUESTED, STOP_SIGNAL = False, None
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.stage == "status":
        return run_status(args.output_dir)
    validate_args(args)
    install_signal_handlers()
    RDLogger.DisableLog("rdApp.*")
    conn = connect_db(args.output_dir)
    if args.reset_stage:
        reset_stage(conn, args.output_dir, args.reset_stage, args.overwrite_completed)
    cfg_hash = config_hash(args)
    log(f"Stage={args.stage}; config_hash={cfg_hash}; output={Path(args.output_dir).resolve()}")
    functions = {
        "index": run_index,
        "split": run_split,
        "materialize": run_materialize,
        "audit": run_audit,
    }
    try:
        completed = functions[args.stage](conn, args, cfg_hash)
    finally:
        conn.close()
    if not completed:
        log(
            f"Stage {args.stage} paused after {STOP_SIGNAL or 'stop request'}. "
            "Submit the identical command again to resume; no new job was submitted automatically."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
