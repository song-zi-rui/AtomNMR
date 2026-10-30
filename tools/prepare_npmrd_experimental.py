#!/usr/bin/env python
"""Download NP-MRD measured assignments and deduplicate against NMRShiftDB2.

Uses official public bulk snapshots; preserves downloads and source numbering.
"""
from __future__ import annotations

import argparse
import collections
import concurrent.futures
import csv
import hashlib
import html
import io
import json
import re
import time
import math
import sqlite3
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data/np_mrd_experimental"
BASE = "https://np-mrd.org/system/downloads/current/"
VERSION = "npmrd-experimental-v1"


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def fetch(url, output_dir=DEFAULT_OUTPUT, filename=None, byte_range=None, quiet=False):
    """Cache a public source and checksum it. Validate partial downloads before rename."""
    output_dir = Path(output_dir)
    filename = filename or Path(urllib.parse.urlparse(url).path).name
    if Path(filename).name != filename or not filename:
        raise ValueError("filename must be a basename")
    raw = output_dir / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    path = raw / filename
    metadata_path = raw / (filename + ".download.json")
    if path.exists() and metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("url") == url and metadata.get("requested_range") == byte_range and metadata.get("sha256") == sha256(path):
            return path
    temporary = raw / (filename + ".part")
    for attempt in range(3):
        try:
            headers = {"User-Agent": "Mozilla/5.0 (compatible; NMRResearchDataset/1.0)", "Accept": "*/*"}
            if byte_range:
                headers["Range"] = byte_range
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=90) as response, temporary.open("wb") as target:
                size = 0
                for chunk in iter(lambda: response.read(1024 * 1024), b""):
                    target.write(chunk)
                    size += len(chunk)
                expected = response.headers.get("Content-Length")
                if expected and size != int(expected):
                    raise IOError(f"incomplete response: {size} / {expected}")
                metadata = {"url": url, "resolved_url": response.url, "downloaded_at": now(),
                            "last_modified": response.headers.get("Last-Modified"), "etag": response.headers.get("ETag"),
                            "bytes": size, "content_type": response.headers.get("Content-Type"),
                            "requested_range": byte_range, "content_range": response.headers.get("Content-Range")}
            if path.suffix == ".zip":
                with zipfile.ZipFile(temporary) as archive:
                    bad = archive.testzip()
                    if bad:
                        raise IOError(f"CRC failure: {bad}")
            metadata["sha256"] = sha256(temporary)
            temporary.replace(path)
            write_json(metadata_path, metadata)
            if not quiet:
                print(f"[fetch] {filename}: {size:,} bytes", flush=True)
            return path
        except (OSError, urllib.error.URLError, zipfile.BadZipFile) as exc:
            if attempt == 2 or isinstance(exc, urllib.error.HTTPError) and exc.code in {400, 401, 403, 404, 410}:
                raise
            time.sleep(2 ** attempt)


def strip_html(value):
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", "", value))).strip()


def page_info(text):
    title = re.search(r"<title>(.*?)</title>", text, re.S)
    conditions = {strip_html(k).rstrip(":"): strip_html(v) for k, v in re.findall(r"<tr[^>]*>\s*<th[^>]*>((?:(?!</?tr).)*?)</th>\s*<td[^>]*>((?:(?!</?tr).)*?)</td>\s*</tr>", text, re.S)}
    paths = re.findall(r'data-(?:spectra_path|path)="([^"]+)"', text)
    links = [html.unescape(s) for s in re.findall(r'href="([^"]+)"', text)]
    reference = re.search(r'id="references".*?<ol>(.*?)</ol>', text, re.S)
    return {"title": strip_html(title.group(1)) if title else None,
            "conditions_raw": conditions, "nmrml_paths": paths,
            "links": links, "references": strip_html(reference.group(1)) if reference else None}


def assignment_members(output_dir):
    with zipfile.ZipFile(Path(output_dir) / "raw/assignment_tables.zip") as archive:
        return sorted(n for n in archive.namelist() if n.endswith("_user_assignmenttable.txt"))


def fetch_sources(output_dir):
    """Resolve each public submission and its own molecular numbering.

    The generic assignmenttable contains computed reference shifts. Only the
    user_assignmenttable is considered, and its current public link must match.
    """
    output_dir = Path(output_dir)
    fetch(BASE + "assignment_tables.zip", output_dir)
    members = assignment_members(output_dir)
    by_submission = collections.defaultdict(list)
    for name in members:
        by_submission[name.split("_")[0]].append(name)

    def resolve(item):
        sid, names = item
        url = f"https://np-mrd.org/natural_products/{sid}/submission_spectrum_view?spectrum_type=simulation"
        result = {"submission_id": sid, "archive_members": names, "page_url": url}
        try:
            page = fetch(url, output_dir, f"submission_{sid}.html", quiet=True)
            info = page_info(page.read_text(encoding="utf-8"))
            result.update(info)
            verified = [name for name in names if any(Path(urllib.parse.urlparse(link).path).name == name for link in info["links"])]
            result["verified_members"] = verified
            if not verified:
                result["status"] = "archive_member_not_linked_by_current_submission"
            elif not info["nmrml_paths"]:
                result["status"] = "no_nmrml_structure"
            else:
                nmrml_url = urllib.parse.urljoin(url, html.unescape(info["nmrml_paths"][0]))
                # Spectrum arrays can be several MB; the actual atom/bond and
                # assignment blocks are at the tail. Save the exact HTTP range.
                tail = fetch(nmrml_url, output_dir, f"submission_{sid}.nmrml.tail", "bytes=-200000", quiet=True)
                text = tail.read_text(encoding="utf-8", errors="replace")
                if "<structure>" not in text or "</structure>" not in text:
                    result["status"] = "incomplete_nmrml_tail"
                else:
                    result["status"] = "resolved"
                result["nmrml_url"] = nmrml_url
                result["nmrml_tail"] = str(tail.relative_to(output_dir))
                assignment_url = next(link for link in info["links"] if Path(urllib.parse.urlparse(link).path).name == verified[0])
                live_table = fetch(urllib.parse.urljoin(url, assignment_url), output_dir, f"submission_{sid}.user_assignmenttable.txt", quiet=True)
                result["live_assignment_path"] = str(live_table.relative_to(output_dir))
                result["live_assignment_url"] = urllib.parse.urljoin(url, assignment_url)
        except Exception as exc:
            result.update(status="fetch_failed", error=f"{type(exc).__name__}: {exc}")
        return result

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(resolve, item) for item in sorted(by_submission.items(), key=lambda item: int(item[0]))]
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())
            if len(results) % 20 == 0:
                print(f"[sources] {len(results)}/{len(futures)} {dict(collections.Counter(r['status'] for r in results))}", flush=True)
                write_json(output_dir / "source_resolution.partial.json", sorted(results, key=lambda r: int(r["submission_id"])))
    results.sort(key=lambda r: int(r["submission_id"]))
    write_json(output_dir / "source_resolution.json", results)
    print(f"[sources] completed {dict(collections.Counter(r['status'] for r in results))}", flush=True)
    return results


def sdf_blocks(handle):
    block = []
    for line in handle:
        block.append(line)
        if line.strip() == "$$$$":
            yield "".join(block)
            block = []
    if any(line.strip() for line in block):
        yield "".join(block)


def properties(block):
    return {key: value.strip() for key, value in re.findall(r">\s*<([^>]+)>[^\n]*\n(.*?)(?=\n\s*\n|\n\$\$\$\$|\Z)", block, re.S)}


def load_np_structures(output_dir, wanted):
    from rdkit import Chem
    output_dir = Path(output_dir)
    path = fetch(BASE + "structures_NP0000001_NP0050000.zip", output_dir)
    structures = {}
    with zipfile.ZipFile(path) as z, z.open(z.namelist()[0]) as handle:
        for block in sdf_blocks(io.TextIOWrapper(handle, encoding="utf-8", errors="replace")):
            props = properties(block)
            np_id = props.get("DATABASE_ID") or props.get("NP-MRD_ID")
            if np_id in wanted:
                structures[np_id] = (block, props)
    missing = sorted(set(wanted) - structures.keys())
    metadata_path = fetch(BASE + "npmrd_natural_products_NP0000001_NP0050000_json.zip", output_dir)
    if missing and metadata_path.exists():
        with zipfile.ZipFile(metadata_path) as z:
            metadata = json.loads(z.read(z.namelist()[0]))["np_mrd"]["natural_product"]
        for row in metadata:
            np_id = row.get("accession")
            if np_id in missing and row.get("smiles"):
                params = Chem.SmilesParserParams()
                params.removeHs = False
                mol = Chem.MolFromSmiles(row["smiles"], params)
                if mol is not None:
                    structures[np_id] = (Chem.MolToMolBlock(mol), {"INCHI_KEY": row.get("inchikey", ""), "structure_source": "official_metadata_smiles"})
        missing = sorted(set(wanted) - structures.keys())

    def get_structure(np_id):
        url = f"https://np-mrd.org/structures/natural_products/{np_id}/download.sdf"
        try:
            path = fetch(url, output_dir, f"{np_id}.sdf", quiet=True)
            block = path.read_text(encoding="utf-8")
            mol = Chem.MolFromMolBlock(block, removeHs=False)
            if mol is None:
                return np_id, None
            return np_id, (block, properties(block))
        except Exception as exc:
            print(f"[structure] {np_id}: {exc}", flush=True)
            return np_id, None

    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for np_id, value in pool.map(get_structure, missing):
            if value:
                structures[np_id] = value
    print(f"[structure] resolved {len(structures)}/{len(wanted)} compounds", flush=True)
    return structures


KEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")


def source_keys(block, props, preserve_hs=False):
    from rdkit import Chem
    keys = set()
    for name in ("INChI key", "INCHI_KEY", "InChIKey", "INCHIKEY"):
        value = props.get(name, "").strip()
        if KEY_RE.fullmatch(value):
            keys.add(value)
    mol = Chem.MolFromMolBlock(block, removeHs=not preserve_hs)
    if mol is not None:
        try:
            key = Chem.MolToInchiKey(mol)
            if KEY_RE.fullmatch(key):
                keys.add(key)
        except Exception:
            pass
    return keys, mol


def index_baseline(source_sdf, output_dir):
    from rdkit import rdBase
    quiet = rdBase.BlockLogs()
    output_dir = Path(output_dir)
    fingerprint = sha256(source_sdf)
    db_path = output_dir / "dedup_audit.sqlite"
    connection = sqlite3.connect(db_path)
    connection.executescript("""
    CREATE TABLE IF NOT EXISTS baseline_records(record_index INTEGER PRIMARY KEY, source_id TEXT, valid_rdkit INTEGER, keys_json TEXT);
    CREATE TABLE IF NOT EXISTS baseline_keys(inchikey TEXT, record_index INTEGER, PRIMARY KEY(inchikey,record_index));
    CREATE TABLE IF NOT EXISTS run_metadata(key TEXT PRIMARY KEY, value TEXT);
    CREATE INDEX IF NOT EXISTS baseline_key_lookup ON baseline_keys(inchikey);
    """)
    signature = f"{fingerprint}:{rdBase.rdkitVersion}:removeHs=True:v1"
    cached = connection.execute("SELECT value FROM run_metadata WHERE key='baseline_index_signature'").fetchone()
    if not cached or cached[0] != signature:
        connection.execute("DELETE FROM baseline_records")
        connection.execute("DELETE FROM baseline_keys")
        with Path(source_sdf).open(encoding="utf-8", errors="replace", newline="") as handle:
            for index, block in enumerate(sdf_blocks(handle), 1):
                props = properties(block)
                keys, mol = source_keys(block, props)
                connection.execute("INSERT INTO baseline_records VALUES (?,?,?,?)", (index, props.get("nmrshiftdb2 ID"), int(mol is not None), json.dumps(sorted(keys))))
                connection.executemany("INSERT INTO baseline_keys VALUES (?,?)", [(key, index) for key in sorted(keys)])
                if index % 10000 == 0:
                    print(f"[baseline] indexed {index:,} records", flush=True)
        connection.execute("INSERT OR REPLACE INTO run_metadata VALUES ('baseline_sha256',?)", (fingerprint,))
        connection.execute("INSERT OR REPLACE INTO run_metadata VALUES ('baseline_index_signature',?)", (signature,))
        connection.commit()
    exact = collections.defaultdict(set)
    connectivity = collections.defaultdict(set)
    for key, index in connection.execute("SELECT inchikey,record_index FROM baseline_keys"):
        exact[key].add(index)
        connectivity[key[:14]].add(index)
    counts = dict(zip(("records", "valid_rdkit", "invalid_rdkit"), connection.execute("SELECT COUNT(*),SUM(valid_rdkit),SUM(1-valid_rdkit) FROM baseline_records").fetchone()))
    counts["records_without_any_key"] = connection.execute("SELECT COUNT(*) FROM baseline_records WHERE keys_json='[]'").fetchone()[0]
    counts["sha256"] = fingerprint
    del quiet
    return connection, exact, connectivity, counts


def parse_user_table(text):
    """Keep per-atom observations, including missing rows; do not fuse nearby shifts."""
    rows, skipped, errors = [], [], []
    for line_number, values in enumerate(csv.reader(io.StringIO(text)), 1):
        if not values or not any(v.strip() for v in values):
            continue
        if len(values) < 3:
            errors.append(f"line {line_number}: fewer than 3 columns")
            continue
        symbol, atom_id, raw_shift = [v.strip() for v in values[:3]]
        if raw_shift.lower() in {"", "na", "n/a", "null", "none"}:
            skipped.append({"line": line_number, "symbol": symbol, "atom_id": atom_id, "raw": values})
            continue
        try:
            atom_id = int(atom_id)
            shift = float(raw_shift)
            if not math.isfinite(shift) or atom_id < 1 or symbol not in {"H", "C"}:
                raise ValueError("invalid atom or shift")
            if symbol == "H" and not -10 <= shift <= 30 or symbol == "C" and not -50 <= shift <= 300:
                raise ValueError("shift outside broad physical sanity range")
        except ValueError as exc:
            errors.append(f"line {line_number}: {exc}")
            continue
        multiplicity = values[3].strip() if len(values) > 3 else ""
        coupling = values[4].strip() if len(values) > 4 else ""
        j_values = []
        if coupling.lower() not in {"", "na", "n/a", "null", "none"}:
            try:
                numbers = [float(s.strip()) for s in coupling.split(",")]
                if all(math.isfinite(n) and abs(n) <= 500 for n in numbers):
                    j_values = numbers
            except ValueError:
                pass
        rows.append({"source_row": line_number, "source_atom_id": atom_id, "symbol": symbol,
                     "nucleus": "13C" if symbol == "C" else "1H", "shift_ppm": shift, "shift_raw": raw_shift,
                     "multiplicity_raw": multiplicity or None, "j_values_raw": coupling or None,
                     "j_values_numeric_unverified_unit": j_values, "j_values_hz": None,
                     "coupling_unit_status": "not_explicit_in_headerless_user_table", "raw_columns": values})
    return rows, skipped, errors


def parse_nmrml_structure(text):
    import xml.etree.ElementTree as ET
    from rdkit import Chem
    atom_block = re.search(r"<atomList>(.*?)</atomList>", text, re.S)
    bond_block = re.search(r"<bondList>(.*?)</bondList>", text, re.S)
    if not atom_block or not bond_block:
        raise ValueError("source nmrML has no complete atom/bond lists")
    atoms = [ET.fromstring(tag).attrib for tag in re.findall(r"<atom\s[^>]+/>", atom_block.group(1))]
    bonds = [ET.fromstring(tag).attrib for tag in re.findall(r"<bond\s[^>]+/>", bond_block.group(1))]
    ids = [int(a["id"]) for a in atoms]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate source atom IDs")
    positions = {atom_id: i for i, atom_id in enumerate(ids)}
    mol = Chem.RWMol()
    for row in atoms:
        atom = Chem.Atom(row["elementType"])
        atom.SetNoImplicit(True)
        atom.SetFormalCharge(int(row.get("formalCharge", 0)))
        atom.SetIntProp("source_atom_id", int(row["id"]))
        mol.AddAtom(atom)
    orders = {"1": Chem.BondType.SINGLE, "2": Chem.BondType.DOUBLE, "3": Chem.BondType.TRIPLE,
              "1.5": Chem.BondType.AROMATIC, "A": Chem.BondType.AROMATIC}
    for bond in bonds:
        left, right = [int(v) for v in bond["atomRefs"].split()]
        mol.AddBond(positions[left], positions[right], orders[bond["order"]])
    result = mol.GetMol()
    Chem.SanitizeMol(result)
    return result, atoms, bonds


def map_source_atoms(source, reference):
    """Map the complete heavy graph, rejecting alternatives across stereo-aware
    symmetry classes. Explicit proton numbering is retained in its source graph;
    pro-R/pro-S labels are never invented from a heavy-atom mapping.
    """
    from rdkit import Chem
    source_heavy = Chem.RemoveHs(source)
    reference_heavy = Chem.RemoveHs(reference)
    if source_heavy.GetNumAtoms() != reference_heavy.GetNumAtoms() or source_heavy.GetNumBonds() != reference_heavy.GetNumBonds():
        raise ValueError("source/reference heavy-graph size mismatch")
    # A complete molecular graph match is required, including bond orders,
    # element, formal charge, isotope and H count, without assuming atom order.
    matches = reference_heavy.GetSubstructMatches(source_heavy, uniquify=False, useChirality=False, maxMatches=4096)
    candidates = []
    for match in matches:
        if all((a.GetAtomicNum(), a.GetFormalCharge(), a.GetIsotope(), a.GetTotalNumHs()) ==
               (reference_heavy.GetAtomWithIdx(j).GetAtomicNum(), reference_heavy.GetAtomWithIdx(j).GetFormalCharge(),
                reference_heavy.GetAtomWithIdx(j).GetIsotope(), reference_heavy.GetAtomWithIdx(j).GetTotalNumHs())
               for a, j in zip(source_heavy.GetAtoms(), match)):
            candidates.append(match)
    if not candidates:
        raise ValueError("source/reference structure or hydrogen-count mismatch")
    if len(matches) >= 4096:
        raise ValueError("too many atom mappings to resolve")
    ranks = list(Chem.CanonicalRankAtoms(reference_heavy, breakTies=False, includeChirality=True))
    if any(len({ranks[match[i]] for match in candidates}) != 1 for i in range(source_heavy.GetNumAtoms())):
        raise ValueError("source/reference mapping ambiguous beyond stereochemistry-aware symmetry")
    match = min(candidates)
    mapping = {a.GetIntProp("source_atom_id"): j for a, j in zip(source_heavy.GetAtoms(), match)}
    mapped = Chem.Mol(reference_heavy)
    for atom_id, index in mapping.items():
        mapped.GetAtomWithIdx(index).SetAtomMapNum(atom_id)
    mapped_smiles = Chem.MolToSmiles(mapped, isomericSmiles=True)
    parsed = Chem.MolFromSmiles(mapped_smiles)
    mapping = {a.GetAtomMapNum(): a.GetIdx() for a in parsed.GetAtoms()}
    for atom in parsed.GetAtoms():
        atom.SetAtomMapNum(0)
    return mapping, parsed, mapped_smiles, len(candidates)


def normalize_conditions(info):
    raw = info.get("conditions_raw", {})
    solvent = raw.get("Solvent")
    aliases = {"water": "H2O", "h2o": "H2O", "d2o": "D2O", "cdcl3": "CDCl3", "chloroform-d": "CDCl3",
               "dmso-d6": "DMSO-d6", "dimethylsulfoxide-d6": "DMSO-d6", "cd3od": "CD3OD", "methanol-d4": "CD3OD"}
    solvent_norm = aliases.get((solvent or "").lower(), solvent)
    temperature_raw = raw.get("Sample Temperature")
    temperature = None
    if temperature_raw and (number := re.fullmatch(r"\s*([+-]?\d+(?:\.\d+)?)\s*(Celsius|C|Kelvin|K)\s*", temperature_raw, re.I)):
        value = float(number[1]) + (273.15 if number[2].lower() in {"celsius", "c"} else 0)
        if 150 <= value <= 600:
            temperature = value
    frequency_raw = raw.get("Frequency")
    frequency = None
    if frequency_raw and (number := re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*MHz\s*", frequency_raw, re.I)):
        frequency = float(number[1])
    source_type = "experimental" if "experimental)" in (info.get("title") or "").lower() else "simulated_from_submitted_assignments"
    return {"solvent": solvent_norm, "solvent_raw": solvent, "temperature_k": temperature,
            "temperature_raw": temperature_raw, "display_frequency_mhz": frequency,
            "frequency_raw": frequency_raw, "acquisition_frequency_mhz": frequency if source_type == "experimental" else None,
            "frequency_role": "acquisition" if source_type == "experimental" else "simulation_display_frequency_not_verified_acquisition",
            "chemical_shift_reference": raw.get("Chemical Shift Reference"), "ph_raw": raw.get("Sample pH"),
            "concentration_raw": raw.get("Sample Concentration"), "spectrum_display_type": source_type}


def write_jsonl(path, records):
    with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n")


def build_dataset(output_dir, source_sdf):
    from rdkit import Chem, rdBase
    from rdkit.Chem import rdMolDescriptors
    output_dir = Path(output_dir)
    quiet = rdBase.BlockLogs()
    members = assignment_members(output_dir)
    wanted = {name.split("_")[1] for name in members}
    structures = load_np_structures(output_dir, wanted)
    connection, exact, connectivity, baseline = index_baseline(source_sdf, output_dir)
    baseline_ids = dict(connection.execute("SELECT record_index,source_id FROM baseline_records"))
    resolution = json.loads((output_dir / "source_resolution.json").read_text(encoding="utf-8"))
    by_sid = {row["submission_id"]: row for row in resolution}
    accepted, quarantined, all_raw = [], [], []
    failures = collections.Counter()
    raw_identity = {}
    for np_id, (block, props) in structures.items():
        keys, _ = source_keys(block, props)
        raw_identity[np_id] = keys
    with zipfile.ZipFile(output_dir / "raw/assignment_tables.zip") as archive:
        for number, member in enumerate(members, 1):
            sid, np_id, *_ = member.split("_")
            text = archive.read(member).decode("utf-8-sig", errors="replace")
            rows, missing, errors = parse_user_table(text)
            raw_record = {"record_id": member.removesuffix(".txt"), "np_mrd_id": np_id, "submission_id": sid,
                          "source_archive_member": member, "assignment_rows": rows, "missing_rows": missing,
                          "parse_errors": errors, "source_member_sha256": hashlib.sha256(archive.read(member)).hexdigest()}
            raw_keys = raw_identity.get(np_id, set())
            raw_exact_hits = sorted(set().union(*(exact.get(key, set()) for key in raw_keys)))
            raw_conn_hits = sorted(set().union(*(connectivity.get(key[:14], set()) for key in raw_keys)))
            raw_record["identity_keys_checked"] = sorted(raw_keys)
            raw_record["dedup"] = {
                "status": "unknown_identity" if not raw_keys else "exact_overlap" if raw_exact_hits else "connectivity_only_overlap" if raw_conn_hits else "novel_connectivity",
                "nmrshiftdb2_exact_record_indices": raw_exact_hits,
                "nmrshiftdb2_connectivity_record_indices": raw_conn_hits,
                "nmrshiftdb2_exact_ids": sorted({baseline_ids[i] for i in raw_exact_hits if baseline_ids[i]}),
                "nmrshiftdb2_connectivity_ids": sorted({baseline_ids[i] for i in raw_conn_hits if baseline_ids[i]}),
                "assignment_mapping_verified": False,
            }
            all_raw.append(raw_record)
            try:
                if not rows:
                    raise ValueError("no_measured_H_or_C_assignments")
                if errors:
                    raise ValueError("invalid_assignment_rows")
                info = by_sid.get(sid) or {}
                if info.get("status") == "incomplete_nmrml_tail" and info.get("nmrml_tail"):
                    tail_text = (output_dir / info["nmrml_tail"]).read_text(encoding="utf-8", errors="replace")
                    if "<structure>" in tail_text and "</structure>" in tail_text:
                        info = {**info, "status": "resolved"}
                if info.get("status") != "resolved":
                    raise ValueError(info.get("status", "source_not_resolved"))
                if member not in info.get("verified_members", []):
                    raise ValueError("archive_member_not_linked_by_current_submission")
                if np_id not in structures:
                    raise ValueError("reference_structure_unavailable")
                # Reparse cached HTML using strict row boundaries (no nested-table bleed).
                info = {**info, **page_info((output_dir / f"raw/submission_{sid}.html").read_text(encoding="utf-8"))}
                if info.get("conditions_raw", {}).get("NP-MRD ID") != np_id:
                    raise ValueError("source_page_compound_mismatch")
                if "predicted)" in (info.get("title") or "").lower():
                    raise ValueError("predicted_spectrum_excluded")
                url = next(link for link in info["links"] if Path(urllib.parse.urlparse(link).path).name == member)
                live = fetch(urllib.parse.urljoin(info["page_url"], url), output_dir, f"submission_{sid}.user_assignmenttable.txt", quiet=True)
                live_rows, _, live_errors = parse_user_table(live.read_text(encoding="utf-8-sig"))
                signature = lambda values: sorted((r["source_atom_id"], r["symbol"], r["shift_ppm"], r["multiplicity_raw"], r["j_values_raw"]) for r in values)
                if live_errors or signature(live_rows) != signature(rows):
                    raise ValueError("archive_vs_live_assignment_mismatch")
                source, atoms_raw, bonds_raw = parse_nmrml_structure((output_dir / info["nmrml_tail"]).read_text(encoding="utf-8", errors="replace"))
                block, props = structures[np_id]
                # Use the canonical heavy-atom graph for identity/deduplication.
                # Explicit hydrogens in an NP-MRD SDF are representation detail and
                # must not make an otherwise identical InChIKey look different.
                keys, reference = source_keys(block, props, preserve_hs=False)
                if reference is None or not keys:
                    raise ValueError("invalid_reference_structure_or_key")
                mapping, heavy, mapped_smiles, count = map_source_atoms(source, reference)
                computed_key = Chem.MolToInchiKey(heavy)
                source_atoms = {a.GetIntProp("source_atom_id"): a for a in source.GetAtoms()}
                unique_ids = set()
                for row in rows:
                    atom_id = row["source_atom_id"]
                    atom = source_atoms.get(atom_id)
                    if atom is None or atom.GetSymbol() != row["symbol"]:
                        raise ValueError("assignment_element_or_atom_id_mismatch")
                    if atom_id in unique_ids:
                        raise ValueError("duplicate_source_atom_assignment")
                    unique_ids.add(atom_id)
                    parent = atom if row["symbol"] == "C" else next((a for a in atom.GetNeighbors() if a.GetAtomicNum() != 1), None)
                    if parent is None:
                        raise ValueError("proton_without_heavy_atom_parent")
                    row["source_atom_index"] = atom.GetIdx()
                    row["source_heavy_atom_id"] = parent.GetIntProp("source_atom_id")
                    row["heavy_atom_index"] = mapping[parent.GetIntProp("source_atom_id")]
                exact_hits = sorted(set().union(*(exact.get(key, set()) for key in keys)))
                conn_hits = sorted(set().union(*(connectivity.get(key[:14], set()) for key in keys)))
                record = {"schema_version": VERSION, "record_id": raw_record["record_id"], "np_mrd_id": np_id,
                          "name": info["conditions_raw"].get("Compound name"), "submission_id": sid,
                          "data_origin": "user_submitted_experimental_chemical_shift_assignments",
                          "structure": {"inchikey": computed_key, "identity_keys_checked": sorted(keys), "connectivity_key": computed_key[:14],
                                        "smiles": Chem.MolToSmiles(heavy, isomericSmiles=True), "mapped_smiles_source_atom_ids": mapped_smiles,
                                        "formula": rdMolDescriptors.CalcMolFormula(heavy),
                                        "source_atoms": atoms_raw, "source_bonds": bonds_raw,
                                        "source_heavy_atom_id_to_heavy_index": mapping,
                                        "heavy_atom_index_space": "RDKit parse of mapped_smiles_source_atom_ids (not the unmapped canonical SMILES)",
                                        "mapping_isomorphisms": count,
                                        "hydrogen_assignment_scope": "explicit source H IDs retained; heavy-atom parent verified; pro-R/pro-S not inferred"},
                          "conditions": normalize_conditions(info),
                          "assignments": {nucleus: [r for r in rows if r["nucleus"] == nucleus] for nucleus in ("13C", "1H")},
                          "missing_assignment_rows": missing,
                          "provenance": {"assignment_archive": "raw/assignment_tables.zip", "assignment_member": member,
                                         "assignment_member_sha256": raw_record["source_member_sha256"], "page_url": info["page_url"],
                                         "live_assignment_url": urllib.parse.urljoin(info["page_url"], url), "nmrml_url": info["nmrml_url"],
                                         "nmrml_tail": info["nmrml_tail"], "references": info.get("references"),
                                         "spectrum_page_title": info.get("title"), "conditions_raw": info.get("conditions_raw")},
                          "dedup": {"status": "exact_overlap" if exact_hits else "connectivity_only_overlap" if conn_hits else "novel_connectivity",
                                    "nmrshiftdb2_exact_record_indices": exact_hits,
                                    "nmrshiftdb2_connectivity_record_indices": conn_hits,
                                    "nmrshiftdb2_exact_ids": sorted({baseline_ids[i] for i in exact_hits if baseline_ids[i]}),
                                    "nmrshiftdb2_connectivity_ids": sorted({baseline_ids[i] for i in conn_hits if baseline_ids[i]})}}
                accepted.append(record)
            except Exception as exc:
                reason = str(exc)
                failures[reason] += 1
                quarantined.append({**raw_record, "reason": reason, "page_url": by_sid.get(sid, {}).get("page_url")})
            if number % 25 == 0:
                print(f"[build] {number}/{len(members)} validated={len(accepted)} review={len(quarantined)}", flush=True)

    # Collapse byte-identical scientific observations/conditions from repeated
    # submissions, preserving all original record IDs as provenance aliases.
    unique = {}
    internal_duplicates = []
    for record in accepted:
        fingerprint = json.dumps({"key": record["structure"]["inchikey"], "conditions": record["conditions"],
                                  "mapped_smiles": record["structure"]["mapped_smiles_source_atom_ids"],
                                  "assignments": [(r["nucleus"], r["source_atom_id"], r["shift_ppm"], r["multiplicity_raw"], r["j_values_raw"])
                                                  for values in record["assignments"].values() for r in values]}, sort_keys=True)
        if fingerprint in unique:
            kept = unique[fingerprint]
            kept.setdefault("duplicate_sources", []).append({"record_id": record["record_id"], "provenance": record["provenance"]})
            internal_duplicates.append({"record_id": record["record_id"], "kept_record_id": kept["record_id"], "reason": "identical_assignment_and_conditions"})
        else:
            unique[fingerprint] = record
    records = sorted(unique.values(), key=lambda r: (r["np_mrd_id"], int(r["submission_id"])))
    novel_exact = [r for r in records if r["dedup"]["status"] != "exact_overlap"]
    novel_connectivity = [r for r in records if r["dedup"]["status"] == "novel_connectivity"]
    overlaps = [r for r in records if r["dedup"]["status"] != "novel_connectivity"]
    for filename, values in (("experimental_assignments.raw.jsonl", all_raw), ("experimental_assignments.validated.jsonl", records),
                             ("novel_exact.jsonl", novel_exact), ("novel_connectivity.jsonl", novel_connectivity),
                             ("overlap_with_nmrshiftdb2.jsonl", overlaps), ("review_queue.jsonl", quarantined),
                             ("internal_duplicates.jsonl", internal_duplicates)):
        write_jsonl(output_dir / filename, values)
    write_jsonl(output_dir / "downloaded_novel_exact.unverified.jsonl", [r for r in all_raw if r["dedup"]["status"] in {"connectivity_only_overlap", "novel_connectivity"}])
    write_jsonl(output_dir / "downloaded_novel_connectivity.unverified.jsonl", [r for r in all_raw if r["dedup"]["status"] == "novel_connectivity"])
    for name, values in (("novel_exact", novel_exact), ("novel_connectivity", novel_connectivity)):
        molecules = {}
        for record in values:
            key = record["structure"]["inchikey"]
            molecules.setdefault(key, {"inchikey": key, "np_mrd_ids": [], "records": []})
            if record["np_mrd_id"] not in molecules[key]["np_mrd_ids"]:
                molecules[key]["np_mrd_ids"].append(record["np_mrd_id"])
            molecules[key]["records"].append(record)
        write_jsonl(output_dir / f"{name}.molecules.jsonl", [molecules[k] for k in sorted(molecules)])
    writer = Chem.SDWriter(str(output_dir / "novel_connectivity.sdf"))
    seen = set()
    for record in novel_connectivity:
        key = record["structure"]["inchikey"]
        if key not in seen:
            seen.add(key)
            mol = Chem.MolFromSmiles(record["structure"]["smiles"])
            mol.SetProp("NP_MRD_ID", record["np_mrd_id"])
            mol.SetProp("INCHIKEY", key)
            mol.SetProp("ASSIGNMENTS_FILE", "novel_connectivity.jsonl")
            writer.write(mol)
    writer.close()
    connection.execute("CREATE TABLE IF NOT EXISTS npmrd_decisions(record_id TEXT PRIMARY KEY,np_mrd_id TEXT,status TEXT,details_json TEXT)")
    connection.execute("DELETE FROM npmrd_decisions")
    connection.executemany("INSERT INTO npmrd_decisions VALUES (?,?,?,?)", [(r["record_id"], r["np_mrd_id"], r["dedup"]["status"], json.dumps(r["dedup"])) for r in records])
    connection.executemany("INSERT INTO npmrd_decisions VALUES (?,?,?,?)", [(r["record_id"], r["np_mrd_id"], "review", json.dumps({"reason": r["reason"]})) for r in quarantined])
    connection.executemany("INSERT INTO npmrd_decisions VALUES (?,?,?,?)", [(r["record_id"], r["record_id"].split("_")[1], "internal_duplicate", json.dumps(r)) for r in internal_duplicates])
    connection.commit()
    connection.close()

    def stats(values):
        compound_nuclei = collections.defaultdict(set)
        paired_solvent = collections.defaultdict(set)
        for r in values:
            key = r["structure"]["inchikey"]
            for nucleus in ("13C", "1H"):
                if r["assignments"][nucleus]:
                    compound_nuclei[key].add(nucleus)
                    if r["conditions"]["solvent"]:
                        paired_solvent[(key, r["conditions"]["solvent"])].add(nucleus)
        return {"records": len(values), "unique_inchikeys": len(compound_nuclei),
                "unique_np_mrd_ids": len({r["np_mrd_id"] for r in values}),
                "spectra_by_nucleus": {n: sum(bool(r["assignments"][n]) for r in values) for n in ("13C", "1H")},
                "atom_assignments_by_nucleus": {n: sum(len(r["assignments"][n]) for r in values) for n in ("13C", "1H")},
                "structures_with_both_nuclei": sum(len(n) == 2 for n in compound_nuclei.values()),
                "structures_with_both_nuclei_same_known_solvent": len({k for (k, solvent), n in paired_solvent.items() if len(n) == 2}),
                "records_with_solvent": sum(bool(r["conditions"]["solvent"]) for r in values),
                "records_with_temperature_k": sum(r["conditions"]["temperature_k"] is not None for r in values),
                "solvents": dict(collections.Counter(r["conditions"]["solvent"] or "unknown" for r in values))}

    report = {"schema_version": VERSION, "finished_at": now(), "baseline": baseline,
              "official_assignment_archive_members": len(zipfile.ZipFile(output_dir / "raw/assignment_tables.zip").namelist()),
              "user_assignment_files": len(members), "user_assignment_compounds": len(wanted),
              "all_downloaded_dedup_counts": dict(collections.Counter(r["dedup"]["status"] for r in all_raw)),
              "all_downloaded_novel_exact_np_ids": len({r["np_mrd_id"] for r in all_raw if r["dedup"]["status"] in {"connectivity_only_overlap", "novel_connectivity"}}),
              "all_downloaded_novel_connectivity_np_ids": len({r["np_mrd_id"] for r in all_raw if r["dedup"]["status"] == "novel_connectivity"}),
              "validated_before_internal_dedup": len(accepted), "internal_duplicates": len(internal_duplicates),
              "review_records": len(quarantined), "review_reasons": dict(failures),
              "validated": stats(records), "novel_exact": stats(novel_exact), "novel_connectivity": stats(novel_connectivity),
              "overlap": stats(overlaps), "dedup_status_counts": dict(collections.Counter(r["dedup"]["status"] for r in records)),
              "scope": "All user assignment files in official assignment_tables.zip; not a claim of exhaustive coverage of every NP-MRD literature or raw-spectrum entry.",
              "license": "CC BY-NC 4.0", "license_url": "https://np-mrd.org/about",
              "source_sdf_unchanged": sha256(source_sdf) == baseline["sha256"]}
    assert len(records) + len(internal_duplicates) + len(quarantined) == len(members)
    assert all(not r["dedup"]["nmrshiftdb2_exact_record_indices"] for r in novel_exact)
    assert all(not r["dedup"]["nmrshiftdb2_connectivity_record_indices"] for r in novel_connectivity)
    write_json(output_dir / "report.json", report)
    write_report(output_dir, report)
    manifest = {"version": VERSION, "created_at": now(), "script_sha256": sha256(__file__), "rdkit_version": rdBase.rdkitVersion,
                "source_sdf": str(Path(source_sdf).resolve()), "source_sdf_sha256": baseline["sha256"],
                "downloads": [json.loads(p.read_text(encoding="utf-8")) for p in sorted((output_dir / "raw").glob("*.download.json"))],
                "outputs": {p.name: {"sha256": sha256(p), "bytes": p.stat().st_size} for p in sorted(output_dir.glob("*.jsonl"))}}
    write_json(output_dir / "manifest.json", manifest)
    del quiet
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def write_report(output_dir, report):
    b = report["baseline"]
    counts = report["dedup_status_counts"]
    lines = ["# NP-MRD 实验原子归属下载与 NMRShiftDB2 去重报告", "",
             f"生成时间（UTC）：{report['finished_at']}", "",
             "## 本次数据范围", "",
             "本次处理 [NP-MRD 官方下载页](https://np-mrd.org/downloads) 的 `assignment_tables.zip` 全部用户归属表。",
             "官方页面标注该压缩包发布日期为 2025-08-01。这是公开批量归属包的处理结果，不代表已经遍历网站所有文献归属、原始谱图或后续新增提交。",
             f"压缩包共 {report['official_assignment_archive_members']} 个文件；其中 `user_assignmenttable` 共 {report['user_assignment_files']} 个，涉及 {report['user_assignment_compounds']} 个 NP-MRD 编号。", "",
             "`assignmenttable` 和 `chsqc_assignmenttable` 中存在程序生成的参考值。本流程仅使用用户表，并核验其官方页面下载链接、实时用户表内容及该谱自身的原子/键表。不能把完整压缩包直接作为实验监督标签。",
             "谱图页面显示 simulated 时，保存的是用户提交的归属位移，不使用模拟曲线的逐点数值、峰强度或线宽。simulated 与 predicted 的定义见 [NP-MRD 2025 论文](https://doi.org/10.1093/nar/gkae1067)。", "",
             "## 数量", "", "| 数据范围 | 归属记录 | 独立完整结构 | ¹³C 谱 | ¹H 谱 | ¹³C 原子归属 | ¹H 原子归属 |", "|---|---:|---:|---:|---:|---:|---:|"]
    for key, title in (("validated", "核验通过、内部重复合并后"), ("novel_exact", "排除完整结构重复"), ("novel_connectivity", "排除连接层重叠")):
        s = report[key]
        lines.append(f"| {title} | {s['records']} | {s['unique_inchikeys']} | {s['spectra_by_nucleus']['13C']} | {s['spectra_by_nucleus']['1H']} | {s['atom_assignments_by_nucleus']['13C']} | {s['atom_assignments_by_nucleus']['1H']} |")
    lines += ["", f"核验通过的原始记录 {report['validated_before_internal_dedup']} 条；同结构、同条件、同原子归属的重复提交合并 {report['internal_duplicates']} 条；待核验/排除 {report['review_records']} 条。",
              f"与 NMRShiftDB2 完整结构重叠 {counts.get('exact_overlap', 0)} 条；只有连接层重叠 {counts.get('connectivity_only_overlap', 0)} 条。", "",
              "## 去重规则", "",
              "1. 按分子身份判定跨库重复，不用分子名、分子式或相近峰表代替结构。保留不同条件下的实验记录，并另提供一分子一行的聚合文件。",
              "2. 完整 InChIKey 一致归为 exact_overlap；前 14 位连接层一致但完整键不同归为 connectivity_only_overlap。后者可能涉及立体、同位素或质子化差异，不称为完全相同分子。",
              "3. 同时检查源文件中的官方 InChIKey 和 RDKit 重新计算的键，避免解析器差异导致漏查。不做中和、脱盐、取最大片段或删立体化学操作。",
              "4. 两个新增数据集均已排除相应跨库结构重叠；现有 NMRShiftDB2 不修改。内部相同实验记录保留一个主记录及全部重复来源。",
              f"本地 SDF 共 {b['records']} 条，其中 RDKit 有效 {b['valid_rdkit']}、无效 {b['invalid_rdkit']}；无法获得任何结构键 {b['records_without_any_key']} 条。无效记录若带有效官方键，仍纳入去重查询。",
              f"源 SDF SHA256：`{b['sha256']}`；结束后校验未修改：`{report['source_sdf_unchanged']}`。", "",
              "## 文件怎么用", "",
              "- `novel_exact.jsonl`：扩充训练数据时，可从这份完整结构去重结果开始。", 
              "- `novel_connectivity.jsonl`：更保守，连接层也与本地库不重叠；用于独立测试时还须与其他训练来源再核对。",
              "- `novel_exact.molecules.jsonl` / `novel_connectivity.molecules.jsonl`：每行一个完整结构，`records` 保存该分子的各次归属记录；不会把不同溶剂的 H/C 强行组合。",
              "- `novel_connectivity.sdf`：上述新增结构的查看文件；归属在配套 JSONL 中，不能直接按该 SDF 的原子行号解释源归属。",
              "- `experimental_assignments.validated.jsonl`：包括重叠结构在内的核验数据。",
              "- `experimental_assignments.raw.jsonl`：全部用户表的解析结果，保留尚未核验的归属；不宜直接训练。",
              "- `overlap_with_nmrshiftdb2.jsonl`：重叠记录和对应本地 SDF 记录序号。",
              "- `review_queue.jsonl`：来源链接失效、结构/编号不一致等未通过记录及具体原因。",
              "- `dedup_audit.sqlite`：本地结构索引与处理决定；`manifest.json`：URL、下载时间、HTTP 分段信息、SHA256 和程序版本。",
              "- `raw/`：原始官方包及核验材料。额外下载的第一分区通用峰表混有 predicted 数据，仅作格式调查材料，不进入输出数据集。", "",
              "## 原子编号和实验条件", "",
              "每条记录在 `assignments.13C` / `assignments.1H` 中逐原子保存位移、原编号、峰型和 J 值。`source_atom_id` 是谱图自身的编号；`source_atom_index` 是该原子列表的 0 基索引；`heavy_atom_index` 是输出 `structure.smiles` 对应的独立重原子图索引（须用下面的 atom-map 方式核对）。",
              "**使用 `structure.mapped_smiles_source_atom_ids` 解析结构，再按 atom-map 数字对应 `source_heavy_atom_id`。不要用重新生成的 canonical SMILES 的自然解析序号替代 atom-map。** 原始原子/键表同时保留在 `source_atoms` / `source_bonds`。",
              "H 的源编号保留到显式氢，并校验其连接重原子；没有据此推断 pro-R/pro-S。¹H 中一个重原子可对应两个不等价氢。化学位移相同的不同原子行不自动认定为结构等价，也不按相近位移合并。",
              "仅当完整重原子图、键型、元素、同位素、电荷和氢数一致，且映射歧义不超出含立体化学的对称组时，才输出核验数据。",
              "条件保留原始值；温度只有单位明确时转成 K。模拟页面的 Frequency 保存在 display_frequency_mhz，不能当作原始实验采集频率。未报告信息为 null。J 值保留原文及数字候选，但无显式单位的无表头用户表不强行标注为 Hz（coupling_unit_status=not_explicit_in_headerless_user_table）。", "",
              "## H/C 成对覆盖", ""]
    for key, title in (("validated", "全部核验数据"), ("novel_connectivity", "连接层去重后的数据")):
        s = report[key]
        lines.append(f"{title}：有 H 和 C 的结构 {s['structures_with_both_nuclei']} 个，其中存在同一已知溶剂 H/C 组合的结构 {s['structures_with_both_nuclei_same_known_solvent']} 个。不同提交之间的温度、样品状态和实测场强仍须逐项核查，不能仅凭同一分子就认定条件完全一致。")
    lines += ["", "## 待核验原因", "", "| 原因 | 记录数 |", "|---|---:|"]
    lines += [f"| {reason} | {n} |" for reason, n in sorted(report["review_reasons"].items(), key=lambda item: (-item[1], item[0]))]
    lines += ["", "## 复现", "", "在项目根目录执行（Python 3.10+，需要 RDKit）：", "", "```text",
              "python tools/prepare_npmrd_experimental.py fetch",
              "python tools/prepare_npmrd_experimental.py sources",
              "python tools/prepare_npmrd_experimental.py build --source-sdf data/nmrshiftdb2/nmrshiftdb2withsignals.sd",
              "```", "", "下载文件经 SHA256 校验后复用；需要获取未来新版时请换一个 `--output-dir`，以保留旧快照。", "",
              "NP-MRD 数据使用 [CC BY-NC 4.0](https://np-mrd.org/about) 许可，使用和发表时须注明来源并遵守非商业用途条件。本次完成数据下载、核验和去重，没有训练模型或证明性能增益。", ""]
    (Path(output_dir) / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["fetch", "sources", "build"])
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--url", action="append")
    parser.add_argument("--source-sdf", type=Path, default=ROOT / "data/nmrshiftdb2/nmrshiftdb2withsignals.sd")
    args = parser.parse_args()
    if args.command == "fetch":
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(lambda url: fetch(url, args.output_dir), args.url or [BASE + "assignment_tables.zip"]))
    elif args.command == "sources":
        fetch_sources(args.output_dir)
    elif args.command == "build":
        build_dataset(args.output_dir, args.source_sdf)


if __name__ == "__main__":
    main()
