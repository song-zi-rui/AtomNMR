from __future__ import annotations

import hashlib
import json
import math
import re
import sqlite3
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data/assignment_database_v1/assignments.sqlite"
SCHEMA = "experimental-assignment-sqlite-v1"
ASSIGNED = {"assigned", "equivalent_group", "overlap", "source_group"}
STATES = ASSIGNED | {"unassigned", "ambiguous", "solvent", "impurity", "excluded"}
METAL_ATOMIC_NUMBERS = {3, 4, 11, 12, 13, 19, 20, *range(21, 32), 37, 38,
                        *range(39, 51), 55, 56, *range(57, 85), 87, 88, *range(89, 117)}


def dumps(value):
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), allow_nan=False)


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def uid(prefix, *values):
    return prefix + hashlib.sha256(dumps(values).encode()).hexdigest()[:24]


QUEUE_ORDER = "CASE queue WHEN 'complete_remaining' THEN 0 WHEN 'both_unassigned' THEN 1 WHEN 'waiting_spectrum' THEN 2 ELSE 3 END"
ACTIVE = " NOT EXISTS(SELECT 1 FROM spectrum_quality q WHERE q.spectrum_id=spectra.id AND q.decision!='passed') AND NOT EXISTS(SELECT 1 FROM curation_state cs WHERE cs.scope='spectrum' AND cs.target_id=spectra.id AND cs.state!='active')"
COVERAGE_SQL = "CASE WHEN c_assigned+h_assigned=0 THEN 'none' WHEN c_count+h_count=c_assigned+h_assigned THEN 'full' ELSE 'partial' END"


def annotation_coverage(total, assigned):
    return {'status': 'none' if not assigned else 'full' if assigned == total else 'partial',
            'assigned_peaks': assigned, 'total_peaks': total,
            'basis': 'saved_assignment_status_not_correctness_certificate'}


def connect(path=DEFAULT_DB, timeout=60):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    c = sqlite3.connect(path, timeout=timeout)
    c.row_factory = sqlite3.Row
    c.execute("PRAGMA foreign_keys=ON")
    if c.execute("PRAGMA journal_mode").fetchone()[0] != "wal":
        c.execute("PRAGMA journal_mode=WAL")
    return c


def is_busy_error(exc):
    code = getattr(exc, "sqlite_errorcode", 0) & 255
    return code in (sqlite3.SQLITE_BUSY, sqlite3.SQLITE_LOCKED)


def initialize(c):
    c.executescript("""
    CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT NOT NULL);
    CREATE TABLE IF NOT EXISTS molecules(
      id TEXT PRIMARY KEY,smiles TEXT NOT NULL UNIQUE,inchikey TEXT,connectivity TEXT,
      formula TEXT NOT NULL,mw REAL NOT NULL,structure_json TEXT NOT NULL,
      origin TEXT NOT NULL,queue TEXT NOT NULL DEFAULT 'waiting_spectrum',
      c_count INTEGER NOT NULL DEFAULT 0,h_count INTEGER NOT NULL DEFAULT 0,
      c_assigned INTEGER NOT NULL DEFAULT 0,h_assigned INTEGER NOT NULL DEFAULT 0,
      task_done INTEGER NOT NULL DEFAULT 0,note TEXT NOT NULL DEFAULT '',
      revision INTEGER NOT NULL DEFAULT 0,created_at TEXT NOT NULL);
    CREATE INDEX IF NOT EXISTS molecule_queue ON molecules(queue,task_done,id);
    CREATE INDEX IF NOT EXISTS molecule_inchikey ON molecules(inchikey);
    CREATE INDEX IF NOT EXISTS molecule_formula ON molecules(formula,id);
    CREATE INDEX IF NOT EXISTS molecule_created ON molecules(created_at DESC,id);
    CREATE TABLE IF NOT EXISTS molecule_activity(
      molecule_id TEXT PRIMARY KEY REFERENCES molecules(id),last_processed_at TEXT NOT NULL,
      last_action TEXT NOT NULL,last_revision INTEGER NOT NULL);
    CREATE INDEX IF NOT EXISTS molecule_activity_time ON molecule_activity(last_processed_at DESC,molecule_id);
    CREATE TABLE IF NOT EXISTS spectra(
      id TEXT PRIMARY KEY,molecule_id TEXT NOT NULL REFERENCES molecules(id),
      nucleus TEXT NOT NULL CHECK(nucleus IN ('13C','1H')),source TEXT NOT NULL,
      source_key TEXT NOT NULL UNIQUE,solvent TEXT,frequency_mhz REAL,temperature_k REAL,
      conditions_json TEXT NOT NULL,peaks_json TEXT NOT NULL,source_json TEXT NOT NULL);
    CREATE INDEX IF NOT EXISTS spectrum_molecule ON spectra(molecule_id,nucleus);
    CREATE TABLE IF NOT EXISTS spectrum_quality(
      spectrum_id TEXT PRIMARY KEY REFERENCES spectra(id),policy TEXT NOT NULL,
      decision TEXT NOT NULL,details_json TEXT NOT NULL,updated_at TEXT NOT NULL);
    CREATE TABLE IF NOT EXISTS spectrum_assignment_audit(
      spectrum_id TEXT PRIMARY KEY REFERENCES spectra(id),policy TEXT NOT NULL,
      validation_json TEXT NOT NULL,identity_json TEXT NOT NULL,updated_at TEXT NOT NULL);
    CREATE TABLE IF NOT EXISTS spectrum_numbering_reviews(
      id INTEGER PRIMARY KEY,molecule_id TEXT NOT NULL REFERENCES molecules(id),
      spectrum_id TEXT NOT NULL REFERENCES spectra(id),candidate_id TEXT NOT NULL,
      decision TEXT NOT NULL CHECK(decision IN ('confirm','reject','undo')),note TEXT NOT NULL,
      flow_json TEXT NOT NULL,before_json TEXT NOT NULL,after_json TEXT NOT NULL,
      target_review_id INTEGER REFERENCES spectrum_numbering_reviews(id),created_at TEXT NOT NULL);
    CREATE INDEX IF NOT EXISTS numbering_review_spectrum ON spectrum_numbering_reviews(spectrum_id,id);
    CREATE INDEX IF NOT EXISTS spectrum_quality_decision ON spectrum_quality(decision);
    CREATE TABLE IF NOT EXISTS candidate_quality(
      smiles TEXT PRIMARY KEY,policy TEXT NOT NULL,decision TEXT NOT NULL,
      details_json TEXT NOT NULL,updated_at TEXT NOT NULL);
    CREATE TABLE IF NOT EXISTS quality_history(
      id INTEGER PRIMARY KEY,molecule_id TEXT,smiles TEXT NOT NULL,
      old_pair_json TEXT,new_pair_json TEXT,details_json TEXT NOT NULL,created_at TEXT NOT NULL);
    CREATE TABLE IF NOT EXISTS annotations(
      spectrum_id TEXT NOT NULL REFERENCES spectra(id),peak_id TEXT NOT NULL,
      status TEXT NOT NULL,atoms_json TEXT NOT NULL,note TEXT NOT NULL DEFAULT '',
      origin TEXT NOT NULL,updated_at TEXT NOT NULL,PRIMARY KEY(spectrum_id,peak_id));
    CREATE TABLE IF NOT EXISTS pairs(
      id TEXT PRIMARY KEY,molecule_id TEXT NOT NULL REFERENCES molecules(id),
      c_id TEXT NOT NULL REFERENCES spectra(id),h_id TEXT NOT NULL REFERENCES spectra(id),
      basis TEXT NOT NULL,solvent_match TEXT NOT NULL,UNIQUE(c_id,h_id));
    CREATE INDEX IF NOT EXISTS pair_molecule ON pairs(molecule_id,id);
    CREATE TABLE IF NOT EXISTS candidates(
      smiles TEXT PRIMARY KEY,exact_key TEXT,lmdb_split TEXT NOT NULL,lmdb_key TEXT NOT NULL,
      raw_smiles TEXT,solvent TEXT,c_shifts TEXT,h_shifts TEXT,
      state TEXT NOT NULL DEFAULT 'indexed',reason TEXT,molecule_id TEXT,
      source_pair_json TEXT);
    CREATE INDEX IF NOT EXISTS candidate_state ON candidates(state,smiles);
    CREATE INDEX IF NOT EXISTS candidate_molecule ON candidates(molecule_id,state);
    CREATE INDEX IF NOT EXISTS candidate_lmdb_order ON candidates(state,lmdb_split,lmdb_key);
    CREATE INDEX IF NOT EXISTS candidate_recovery_order ON candidates(state,(molecule_id IS NULL),smiles);
    CREATE TABLE IF NOT EXISTS source_rows(
      id TEXT PRIMARY KEY,raw_smiles TEXT NOT NULL,nucleus TEXT NOT NULL,
      doi_encoded TEXT,solvent TEXT,page REAL,shard TEXT NOT NULL,row_index INTEGER NOT NULL,
      raw_json TEXT NOT NULL);
    CREATE INDEX IF NOT EXISTS source_smiles ON source_rows(raw_smiles,nucleus);
    CREATE TABLE IF NOT EXISTS import_issues(
      id TEXT PRIMARY KEY,stage TEXT NOT NULL,source_key TEXT NOT NULL,reason TEXT NOT NULL);
    CREATE TABLE IF NOT EXISTS legacy_annotations(
      task_id TEXT PRIMARY KEY,legacy_id TEXT,payload_json TEXT NOT NULL,state TEXT NOT NULL);
    CREATE TABLE IF NOT EXISTS progress(stage TEXT PRIMARY KEY,cursor TEXT NOT NULL,updated_at TEXT NOT NULL);
    CREATE TABLE IF NOT EXISTS curation_state(
      scope TEXT NOT NULL,target_id TEXT NOT NULL,molecule_id TEXT NOT NULL REFERENCES molecules(id),
      state TEXT NOT NULL,reason TEXT NOT NULL,updated_at TEXT NOT NULL,PRIMARY KEY(scope,target_id));
    CREATE INDEX IF NOT EXISTS curation_molecule ON curation_state(molecule_id);
    CREATE TABLE IF NOT EXISTS curation_events(
      id INTEGER PRIMARY KEY,molecule_id TEXT NOT NULL REFERENCES molecules(id),action TEXT NOT NULL,
      details_json TEXT NOT NULL,created_at TEXT NOT NULL);
    CREATE INDEX IF NOT EXISTS curation_event_molecule ON curation_events(molecule_id,id);
    CREATE TABLE IF NOT EXISTS curation_issues(
      id INTEGER PRIMARY KEY,molecule_id TEXT NOT NULL REFERENCES molecules(id),spectrum_id TEXT REFERENCES spectra(id),
      category TEXT NOT NULL,details TEXT NOT NULL,status TEXT NOT NULL,created_at TEXT NOT NULL,updated_at TEXT NOT NULL);
    CREATE INDEX IF NOT EXISTS curation_issue_molecule ON curation_issues(molecule_id,id);
    CREATE TABLE IF NOT EXISTS spectrum_completeness(
      spectrum_id TEXT PRIMARY KEY REFERENCES spectra(id),status TEXT NOT NULL,reason TEXT NOT NULL,updated_at TEXT NOT NULL);
    """)
    c.execute(f"CREATE INDEX IF NOT EXISTS molecule_priority ON molecules(({QUEUE_ORDER}),id)")
    c.execute(f"CREATE INDEX IF NOT EXISTS molecule_pending_priority ON molecules(task_done,({QUEUE_ORDER}),id)")
    c.execute(f"CREATE INDEX IF NOT EXISTS molecule_coverage_processing ON molecules(({COVERAGE_SQL}),task_done,({QUEUE_ORDER}),id)")
    if not c.execute("SELECT 1 FROM meta WHERE key='activity_backfill_v1'").fetchone():
        # Only recover actual timestamped operations, never infer a save time
        # from creation time or revision. Note-only old saves were not timestamped.
        for row in c.execute("SELECT s.molecule_id,max(a.updated_at) stamp FROM annotations a JOIN spectra s ON s.id=a.spectrum_id WHERE a.origin='manual' GROUP BY s.molecule_id").fetchall():
            touch_activity(c, row['molecule_id'], 'historical_annotation', row['stamp'])
        for row in c.execute("SELECT molecule_id,max(created_at) stamp FROM curation_events GROUP BY molecule_id").fetchall():
            touch_activity(c, row['molecule_id'], 'historical_management', row['stamp'])
        set_meta(c, 'activity_backfill_v1', True)
    set_meta(c, "schema", SCHEMA)
    c.commit()
    from . import workload
    workload.initialize(c)
    from . import ai_database
    ai_database.initialize(c)
    c.commit()


def set_meta(c, key, value):
    c.execute("INSERT INTO meta VALUES(?,?) ON CONFLICT(key) DO UPDATE SET value=excluded.value", (key, dumps(value)))


def touch_activity(c, mid, action, stamp=None, *, details=None):
    stamp = stamp or datetime.now(timezone.utc).isoformat(timespec='microseconds')
    c.execute("""INSERT INTO molecule_activity
      SELECT id,?,?,revision FROM molecules WHERE id=?
      ON CONFLICT(molecule_id) DO UPDATE SET last_processed_at=excluded.last_processed_at,
      last_action=excluded.last_action,last_revision=excluded.last_revision
      WHERE excluded.last_processed_at>=molecule_activity.last_processed_at""", (stamp, action, mid))
    # Work history is committed/rolled back together with the corresponding edit.
    # Legacy backfills update only the pointer; workload.initialize recovers them.
    if not action.startswith('historical_') and c.execute("SELECT 1 FROM meta WHERE key='workload_started_at'").fetchone():
        from .workload import append
        append(c, mid, action, stamp, details=details)


def issue(c, stage, source_key, reason):
    c.execute("INSERT OR IGNORE INTO import_issues VALUES(?,?,?,?)", (uid("I", stage, source_key, reason), stage, source_key, reason))


def progress(c, stage, cursor):
    c.execute("INSERT INTO progress VALUES(?,?,?) ON CONFLICT(stage) DO UPDATE SET cursor=excluded.cursor,updated_at=excluded.updated_at", (stage, str(cursor), now()))
    c.commit()


def cursor(c, stage, default="0"):
    r = c.execute("SELECT cursor FROM progress WHERE stage=?", (stage,)).fetchone()
    return r[0] if r else default


def canonicalize(mol):
    """Compatibility return values backed by verified, full atom identities."""
    from tools.shiftdb_atom_identity import canonical_identity
    identity = canonical_identity(mol)
    return identity.smiles, identity.processed, identity.legacy_mapping


def structure(mol):
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=True))
    groups = defaultdict(list)
    atoms = []
    for a in mol.GetAtoms():
        groups[ranks[a.GetIdx()]].append(a.GetIdx())
        atoms.append({"index": a.GetIdx(), "symbol": a.GetSymbol(), "total_hydrogens": a.GetTotalNumHs(includeNeighbors=True),
                      "aromatic": a.GetIsAromatic(), "formal_charge": a.GetFormalCharge(),
                      "neighbors": [b.GetIdx() for b in a.GetNeighbors()]})
    return {"atoms": atoms, "equivalent_groups": [x for x in groups.values() if len(x) > 1],
            "atom_index_basis": "zero_based_canonical_smiles", "h_assignment_resolution": "attached_heavy_atom",
            "ring_count": rdMolDescriptors.CalcNumRings(mol)}


def add_molecule(c, smiles, origin):
    existing = c.execute("SELECT id FROM molecules WHERE smiles=?", (smiles,)).fetchone()
    if existing:
        return existing[0]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid structure")
    key = Chem.MolToInchiKey(mol)
    mid = uid("M", smiles)
    c.execute("INSERT INTO molecules(id,smiles,inchikey,connectivity,formula,mw,structure_json,origin,created_at) VALUES(?,?,?,?,?,?,?,?,?)",
              (mid, smiles, key or None, key.split("-")[0] if key else smiles,
               rdMolDescriptors.CalcMolFormula(mol), Descriptors.MolWt(mol), dumps(structure(mol)), origin, now()))
    return mid


def add_spectrum(c, mid, nucleus, source, source_key, conditions, peaks, provenance, assignments=None):
    sid = uid("S", source, source_key)
    c.execute("INSERT OR IGNORE INTO spectra VALUES(?,?,?,?,?,?,?,?,?,?,?)", (
        sid, mid, nucleus, source, source_key, conditions.get("solvent"), conditions.get("frequency_mhz"),
        conditions.get("temperature_k"), dumps(conditions), dumps(peaks), dumps(provenance)))
    from .peak_table import effective_peaks
    for p in effective_peaks(peaks, provenance):
        a = (assignments or {}).get(p["peak_id"])
        if a:
            c.execute("INSERT OR IGNORE INTO annotations VALUES(?,?,?,?,?,?,?)", (sid, p["peak_id"], a["status"], dumps(a["atom_indices"]), a.get("note", ""), source, now()))
    if provenance.get('atom_identity') and provenance.get('source_assignment_validation'):
        c.execute('INSERT OR IGNORE INTO spectrum_assignment_audit VALUES(?,?,?,?,?)',
                  (sid, provenance['atom_identity']['policy'], dumps(provenance['source_assignment_validation']),
                   dumps(provenance['atom_identity']), now()))
    return sid


def add_pair(c, mid, cs, hs, basis):
    a, b = [c.execute("SELECT solvent FROM spectra WHERE id=?", (s,)).fetchone()[0] for s in (cs, hs)]
    match = "same_known" if a and a == b and a.lower() != "mixed" else "same_label_only" if a and a == b else "different" if a and b else "unknown"
    c.execute("INSERT INTO pairs VALUES(?,?,?,?,?,?) ON CONFLICT(c_id,h_id) DO UPDATE SET solvent_match=excluded.solvent_match", (uid("P", cs, hs), mid, cs, hs, basis, match))


def refresh_molecule(c, mid):
    from .peak_table import effective_peaks
    totals = {"13C": [0, 0], "1H": [0, 0]}
    for r in c.execute("SELECT id,nucleus,peaks_json,source_json FROM spectra WHERE molecule_id=? AND" + ACTIVE, (mid,)):
        annotations = {a["peak_id"]: a["status"] for a in c.execute("SELECT peak_id,status FROM annotations WHERE spectrum_id=?", (r["id"],))}
        for p in effective_peaks(json.loads(r["peaks_json"]), json.loads(r['source_json'])):
            totals[r["nucleus"]][0] += 1
            totals[r["nucleus"]][1] += annotations.get(p["peak_id"]) in ASSIGNED
    cn, ca = totals["13C"]
    hn, ha = totals["1H"]
    if not cn and not hn:
        queue = "quality_hold"
    elif not cn or not hn:
        queue = "waiting_spectrum"
    elif ca == cn and ha == hn:
        queue = "complete"
    elif ca or ha:
        queue = "complete_remaining"
    else:
        queue = "both_unassigned"
    manual = c.execute("SELECT state FROM curation_state WHERE scope='molecule' AND target_id=?", (mid,)).fetchone()
    if manual and manual[0] != 'active':
        queue = 'manual_hold' if manual[0] == 'held' else 'excluded_training'
    c.execute("UPDATE molecules SET queue=?,c_count=?,h_count=?,c_assigned=?,h_assigned=? WHERE id=?", (queue, cn, hn, ca, ha, mid))


def source_assignment_targets(peak, provenance):
    """Read explicit source indices; this is evidence, not a new assignment."""
    mapping = provenance.get('original_atom_mapping', {})
    targets = set()
    if 'source_assignment_fragments' in peak:
        return sorted({a['attached_heavy_processed_idx'] for f in peak['source_assignment_fragments']
                       for a in f['assignments'] if a.get('attached_heavy_processed_idx') is not None})
    for raw in peak.get('raw_source_atom_indices', []):
        if raw is None:
            continue
        for index in re.split(r'[,;_\s]+', str(raw).strip()):
            if index in mapping:
                targets.add(mapping[index])
    return sorted(targets)


def record(c, mid, require_quality=False, include_archived=False):
    r = c.execute("SELECT * FROM molecules WHERE id=?", (mid,)).fetchone()
    if not r:
        raise KeyError(mid)
    result = dict(r)
    activity = c.execute("SELECT * FROM molecule_activity WHERE molecule_id=?", (mid,)).fetchone()
    result['activity'] = dict(activity) if activity else None
    result['navigation_rank'] = queue_rank(result['queue'])
    result["structure"] = json.loads(result.pop("structure_json"))
    mol = Chem.MolFromSmiles(result['smiles'])
    metal_ids = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in METAL_ATOMIC_NUMBERS}
    if metal_ids:
        result['structure']['metal_atoms'] = sorted(metal_ids)
        result['structure']['metal_connections'] = [
            {'atoms': [b.GetBeginAtomIdx(), b.GetEndAtomIdx()], 'bond_type': str(b.GetBondType())}
            for b in mol.GetBonds() if b.GetBeginAtomIdx() in metal_ids or b.GetEndAtomIdx() in metal_ids]
    result["spectra"] = []
    result["quality_issues"] = []
    from . import curation
    result["curation"] = curation.summary(c, mid)
    result["training_eligible"] = result["curation"]["molecule_state"] == "active"
    for s in c.execute("SELECT * FROM spectra WHERE molecule_id=? ORDER BY nucleus,source,id", (mid,)):
        quality = c.execute("SELECT details_json FROM spectrum_quality WHERE spectrum_id=?", (s["id"],)).fetchone()
        quality = json.loads(quality[0]) if quality else {"decision": "not_checked" if s["source"] == "NMRexp" else "human_adopted_reference" if s["source"] in ('BMRB_reference', 'manual_revision') else "original_experimental_source"}
        if quality["decision"] == "quarantined" or (require_quality and quality["decision"] == "not_checked"):
            result["quality_issues"].append({"spectrum_id": s["id"], "nucleus": s["nucleus"], **quality})
            if not include_archived:
                continue
        manual_state = curation.state_for(c, 'spectrum', s['id'])
        if manual_state != 'active' and not include_archived:
            continue
        item = dict(s)
        item["quality"] = quality
        item['manual_state'] = manual_state
        item['editable'] = (result['curation']['molecule_state'] == 'active'
                            and manual_state == 'active' and quality['decision'] != 'quarantined')
        completeness = c.execute("SELECT status,reason FROM spectrum_completeness WHERE spectrum_id=?", (s['id'],)).fetchone()
        item["completeness"] = dict(completeness) if completeness else {"status": "unknown", "reason": "未核实完整性；已归属峰数不代表谱完整。"}
        item["conditions"] = json.loads(item.pop("conditions_json"))
        item["peaks"] = json.loads(item.pop("peaks_json"))
        item["provenance"] = json.loads(item.pop("source_json"))
        from .spectrum_versions import annotate_role
        annotate_role(item)
        from .peak_table import effective_peaks
        item['peaks'] = effective_peaks(item['peaks'], item['provenance'])
        # The side table repairs imported identities without rewriting or
        # pretending to validate human annotations made since import.
        audit = c.execute('SELECT validation_json,identity_json FROM spectrum_assignment_audit WHERE spectrum_id=?', (s['id'],)).fetchone()
        item['source_assignment_validation'] = json.loads(audit[0]) if audit else item['provenance'].get('source_assignment_validation')
        item['atom_identity'] = json.loads(audit[1]) if audit else item['provenance'].get('atom_identity')
        item['numbering_reviews'] = []
        for reviewed in c.execute('SELECT id,candidate_id,decision,note,flow_json,target_review_id,created_at FROM spectrum_numbering_reviews WHERE spectrum_id=? ORDER BY id DESC', (s['id'],)):
            entry=dict(reviewed);entry['flow']=json.loads(entry.pop('flow_json'))
            item['numbering_reviews'].append(entry)
        annotations = {}
        for a in c.execute("SELECT * FROM annotations WHERE spectrum_id=?", (s["id"],)):
            v = dict(a)
            v["atom_indices"] = json.loads(v.pop("atoms_json"))
            annotations[v["peak_id"]] = v
        for p in item["peaks"]:
            p["assignment"] = annotations.get(p["peak_id"], {"status": "unassigned", "atom_indices": [], "note": "", "origin": None})
            v = item['source_assignment_validation']
            p['assignment']['source_identity_eligible'] = (p['assignment'].get('origin') != 'nmrshiftdb2'
                or not v or v['status'] == 'identity_validated')
            if not p['assignment']['source_identity_eligible']:
                p['untrusted_source_assignment'] = dict(p['assignment'])
                p['assignment'] = {**p['assignment'], 'status': 'unassigned', 'atom_indices': []}
            if metal_ids and item['nucleus'] == '1H' and (item['source'] == 'nmrshiftdb2' or item['provenance'].get('original_source') == 'nmrshiftdb2'):
                p['source_assignment_candidates'] = source_assignment_targets(p, item['provenance'])
                if any(result['structure']['atoms'][i]['total_hydrogens'] == 0 for i in p['source_assignment_candidates']):
                    p['source_assignment_note'] = '原库明确给出此重原子归属；当前金属键编码使推算氢数不可靠，需人工核对。'
        if any(str(p['assignment'].get('origin') or '').startswith('ai:') for p in item['peaks']):
            from .ai_database import peak_provenance
            machine=peak_provenance(c,s['id'])
            for p in item['peaks']:
                if str(p['assignment'].get('origin') or '').startswith('ai:'):
                    p['assignment']['ai']=machine.get(p['peak_id'])
        from .spectrum_revision import project
        project(item, result['structure']['atoms'])
        item['annotation_coverage'] = annotation_coverage(len(item['peaks']),
            sum(p['assignment']['status'] in ASSIGNED for p in item['peaks']))
        result["spectra"].append(item)
    result['spectra'].sort(key=lambda s: (s['nucleus'], s['manual_state'] != 'active', not bool(s.get('version_info')),
                                        -(s.get('version_info') or {}).get('number', 0), s['id']))
    active = {s["id"] for s in result["spectra"]
              if s['manual_state'] == 'active' and s['quality']['decision'] != 'quarantined'}
    result["pairs"] = [dict(x) for x in c.execute("SELECT * FROM pairs WHERE molecule_id=? ORDER BY id", (mid,)) if x["c_id"] in active and x["h_id"] in active]
    if require_quality:
        result["export_quality_filtered"] = True
        for nuc, prefix in (("13C", "c"), ("1H", "h")):
            peaks = [p for s in result["spectra"] if s["nucleus"] == nuc for p in s["peaks"]]
            result[prefix + "_count"] = len(peaks)
            result[prefix + "_assigned"] = sum(p["assignment"]["status"] in ASSIGNED for p in peaks)
    active_peaks = [p for s in result['spectra'] if s['id'] in active for p in s['peaks']]
    result['annotation_coverage'] = annotation_coverage(len(active_peaks),
        sum(p['assignment']['status'] in ASSIGNED for p in active_peaks))
    result['processing_status'] = 'finished' if result['task_done'] else 'pending'
    result["schema_version"] = SCHEMA
    result['annotation_sources']=sorted({p['assignment'].get('origin') for s in result['spectra'] for p in s['peaks'] if p['assignment'].get('origin')})
    result['ai_history']=[dict(row) for row in c.execute('SELECT id,job_id,provider,model,response_model,created_at,applied_count FROM ai_annotation_results WHERE molecule_id=? ORDER BY created_at DESC LIMIT 30',(mid,))]
    result["training_eligible"] = result["training_eligible"] and bool(active)
    result["training_eligibility_note"] = "Operational admission only; not evidence of complete spectra or correct assignments."
    return result


def save_annotations(c, mid, payload):
    try:
        c.execute("BEGIN IMMEDIATE")
        # Validate only the small fields required by a save. Do not load source
        # documents / pair records twice while holding the database writer lock.
        molecule = c.execute("SELECT revision,structure_json FROM molecules WHERE id=?", (mid,)).fetchone()
        if not molecule:
            raise KeyError(mid)
        if payload.get("revision") != molecule["revision"]:
            raise RuntimeError("记录已在其他窗口更新，请重新加载后保存。")
        from .curation import state_for
        if state_for(c, 'molecule', mid) != 'active':
            raise ValueError("分子已停用或排除训练，请在数据管理中恢复后再标注。")
        allowed = {s["id"]: {"id": s["id"], "nucleus": s["nucleus"], "source": s['source'],
                             "provenance": json.loads(s['source_json']), "peaks": json.loads(s["peaks_json"])}
                   for s in c.execute("SELECT id,nucleus,source,source_json,peaks_json FROM spectra WHERE molecule_id=? AND" + ACTIVE, (mid,))}
        from .peak_table import effective_peaks
        for s in allowed.values():
            s['peaks'] = effective_peaks(s['peaks'], s['provenance'])
        atoms = json.loads(molecule["structure_json"])["atoms"]
        has_metal = any(Chem.GetPeriodicTable().GetAtomicNumber(a['symbol']) in METAL_ATOMIC_NUMBERS for a in atoms)
        seen = set()
        for row in payload.get("assignments", []):
            sid, pid = row.get("spectrum_id"), row.get("peak_id")
            if (sid, pid) in seen or sid not in allowed or pid not in {p["peak_id"] for p in allowed[sid]["peaks"]}:
                raise ValueError("峰不属于当前分子，或提交中存在重复峰。")
            seen.add((sid, pid))
            from .spectrum_revision import group_for
            if group_for(allowed[sid]['provenance'], pid):
                raise ValueError('此峰属于位移组，组内顺序未定。请在管理当前谱中先解除该组，再编辑逐峰归属。')
            status = row.get("status")
            targets = row.get("atom_indices", [])
            if status not in STATES or not isinstance(targets, list) or any(type(i) is not int or not 0 <= i < len(atoms) for i in targets):
                raise ValueError("归属状态或原子编号无效。")
            if len(targets) != len(set(targets)):
                raise ValueError("归属原子重复。")
            if status in ASSIGNED and not targets:
                raise ValueError("已归属的峰必须选择原子。")
            if status == "assigned" and len(targets) != 1:
                raise ValueError("多个原子请选择等价组或重叠峰。")
            if status not in ASSIGNED | {"ambiguous"} and targets:
                raise ValueError("未归属或排除的峰不能保存确定原子。")
            for i in targets:
                if allowed[sid]["nucleus"] == "13C" and atoms[i]["symbol"] != "C":
                    raise ValueError("碳谱只能归属到碳原子。")
                current_peak = next(p for p in allowed[sid]['peaks'] if p['peak_id'] == pid)
                source_supported = (has_metal and (allowed[sid]['source'] == 'nmrshiftdb2' or allowed[sid]['provenance'].get('original_source') == 'nmrshiftdb2')
                                    and i in source_assignment_targets(current_peak, allowed[sid]['provenance'])
                                    and atoms[i]['symbol'] not in ('H', 'Co')
                                    and Chem.GetPeriodicTable().GetAtomicNumber(atoms[i]['symbol']) not in METAL_ATOMIC_NUMBERS)
                if allowed[sid]["nucleus"] == "1H" and not atoms[i]["total_hydrogens"] and not source_supported:
                    raise ValueError("氢谱请选择带氢的重原子环境。")
            c.execute("INSERT INTO annotations VALUES(?,?,?,?,?,?,?) ON CONFLICT(spectrum_id,peak_id) DO UPDATE SET status=excluded.status,atoms_json=excluded.atoms_json,note=excluded.note,origin=excluded.origin,updated_at=excluded.updated_at",
                      (sid, pid, status, dumps(targets), str(row.get("note", ""))[:5000], "manual", now()))
        if payload.get("task_done"):
            if not allowed:
                raise ValueError("没有通过检查的实验谱，不能标为完成。")
            # Finishing records the annotator's pass, independently of coverage.
            # Empty assignments remain unknown; do not synthesize annotations.
        c.execute("UPDATE molecules SET revision=revision+1,note=?,task_done=? WHERE id=?", (str(payload.get("note", ""))[:10000], int(bool(payload.get("task_done"))), mid))
        refresh_molecule(c, mid)
        touch_activity(c, mid, 'finish' if payload.get('task_done') else 'save', details={
            'edited_peaks': [{'spectrum_id': sid, 'peak_id': pid, 'nucleus': allowed[sid]['nucleus']}
                             for sid, pid in sorted(seen)]})
        # Build the acknowledgement within the same transaction, so it cannot
        # accidentally acknowledge a different writer's later revision. Any
        # failure here still rolls back; after commit no further DB read can fail.
        result = record(c, mid, include_archived=True)
        c.commit()
        return result
    except Exception:
        c.rollback()
        raise


def selection_where(params):
    tests, vals = [], []
    for name in ("queue", "origin"):
        if params.get(name):
            if name == 'queue' and params[name] == 'dual':
                tests.append("queue IN ('complete_remaining','both_unassigned','complete')")
            else:
                tests.append(f"{name}=?")
                vals.append(params[name])
    if params.get('coverage'):
        if params['coverage'] not in ('full', 'partial', 'none'):
            raise ValueError('归属情况筛选无效。')
        tests.append(f'({COVERAGE_SQL})=?')
        vals.append(params['coverage'])
    if params.get('annotation_source'):
        selected=params['annotation_source']
        if selected=='ai_result':
            tests.append('EXISTS(SELECT 1 FROM ai_annotation_results ar WHERE ar.molecule_id=molecules.id)')
        elif selected in ('ai','manual','shiftdb2'):
            origin="a.origin>='ai:' AND a.origin<'ai;'" if selected=='ai' else "a.origin IN ('manual','manual_v2')" if selected=='manual' else "a.origin='nmrshiftdb2'"
            tests.append("EXISTS(SELECT 1 FROM annotations a JOIN spectra s ON s.id=a.spectrum_id WHERE s.molecule_id=molecules.id AND a.status IN ('assigned','equivalent_group','overlap','source_group') AND "+origin+")")
        else:raise ValueError('归属来源筛选无效')
    q = str(params.get('q', '')).strip()
    if q:
        if re.fullmatch(r'M[0-9a-f]{24}', q, re.I):
            tests.append('id=?')
            vals.append('M' + q[1:].lower())
        elif re.fullmatch(r'[A-Z]{14}-[A-Z]{10}-[A-Z]', q, re.I):
            tests.append('inchikey=?')
            vals.append(q.upper())
        elif is_formula_query(q):
            tests.append('formula=?')
            vals.append(q)
        else:
            tests.append("(id LIKE ? ESCAPE '!' OR smiles LIKE ? ESCAPE '!' OR formula LIKE ? ESCAPE '!' OR inchikey LIKE ? ESCAPE '!')")
            literal = q.replace('!', '!!').replace('%', '!%').replace('_', '!_')
            vals.extend(['%' + literal + '%'] * 4)
    if params.get('view') == 'recent':
        tests.append('EXISTS(SELECT 1 FROM molecule_activity WHERE molecule_id=molecules.id)')
    if params.get('view') == 'finished' or params.get('finished'):
        tests.append('task_done=1')
    elif params.get("pending"):
        tests.append("task_done=0")
    return (" WHERE " + " AND ".join(tests)) if tests else "", vals


def is_formula_query(q):
    if not re.fullmatch(r'(?:[A-Z][a-z]?\d*)+', q) or not re.search(r'\d', q):
        return False
    symbols = re.findall(r'[A-Z][a-z]?', q)
    # Ring SMILES such as C1CC1 repeat an element; don't treat them as formulas.
    return len(symbols) == len(set(symbols))


def queue_rank(queue):
    return {'complete_remaining': 0, 'both_unassigned': 1, 'waiting_spectrum': 2}.get(queue, 3)


def listing(c, params):
    where, vals = selection_where(params)
    limit = min(100, max(1, int(params.get("limit", 40))))
    offset = max(0, int(params.get("offset", 0)))
    view = params.get('view', '')
    sort = 'created_at' if view == 'imported' else "coalesce(last_processed_at,'')" if view in ('recent', 'finished') else QUEUE_ORDER
    descending = view in ('recent', 'finished', 'imported')
    source = 'molecule_activity a JOIN molecules ON molecules.id=a.molecule_id' if view == 'recent' else 'molecules LEFT JOIN molecule_activity a ON a.molecule_id=molecules.id'
    order_sort = sort
    if params.get('queue') and params['queue'] != 'dual' and not params.get('q') and not descending:
        if params.get('pending'):
            # Queue + pending is a small, contiguous range in this index.
            # An ordering index over ALL pending molecules otherwise scans the
            # whole database when a queue contains fewer rows than the limit.
            source = 'molecules INDEXED BY molecule_queue LEFT JOIN molecule_activity a ON a.molecule_id=molecules.id'
            order_sort = 'id'  # Queue rank is constant; keep it in cursors.
        else:
            source = 'molecules INDEXED BY molecule_priority LEFT JOIN molecule_activity a ON a.molecule_id=molecules.id'
            where += (' AND ' if where else ' WHERE ') + f'({QUEUE_ORDER})=?'
            vals.append(queue_rank(params['queue']))
            order_sort = 'id'
    # Count molecules, not save events, across all pages of the same filters.
    # Navigation requests don't display totals and need no extra count query.
    backwards = bool(params.get('before'))
    if backwards and params.get('after'):
        raise ValueError('不能同时指定上一个和下一个游标。')
    total = c.execute('SELECT count(*) FROM ' + source + where, vals).fetchone()[0] if view == 'recent' and not (params.get('after') or backwards) else None
    if params.get('after') or backwards:
        anchor = json.loads(params.get('before') or params['after'])
        if not isinstance(anchor, list) or len(anchor) != 3 or anchor[0] != view or not isinstance(anchor[2], str):
            raise ValueError('下一分子的筛选游标无效，请重新打开当前分子。')
        value = str(anchor[1]) if descending else int(anchor[1])
        comparator = '<' if descending != backwards else '>'
        id_comparator = '<' if backwards else '>'
        where += (' AND ' if where else ' WHERE ') + f'(({sort}) {comparator} ? OR (({sort})=? AND id{id_comparator}?)) AND id!=?'
        vals.extend([value, value, anchor[2], anchor[2]])
        offset = 0
    order = f"{order_sort} {'DESC' if descending != backwards else 'ASC'}"
    if order_sort != 'id':
        order += f",id {'DESC' if backwards else 'ASC'}"
    rows = c.execute("SELECT id,formula,smiles,origin,queue,c_count,h_count,c_assigned,h_assigned,task_done,revision,created_at,last_processed_at,last_action," + sort + " AS sort_value FROM " + source + where + f" ORDER BY {order} LIMIT ? OFFSET ?", [*vals, limit, offset])
    result = []
    for row in rows:
        item = dict(row)
        item['cursor'] = [view, item.pop('sort_value'), item['id']]
        item['annotation_coverage'] = annotation_coverage(item['c_count'] + item['h_count'],
                                                        item['c_assigned'] + item['h_assigned'])
        item['processing_status'] = 'finished' if item['task_done'] else 'pending'
        result.append(item)
    return {"rows": result, "offset": offset, "limit": limit, "total": total}


def stats(c):
    result = {"molecules": c.execute("SELECT count(*) FROM molecules").fetchone()[0],
              "spectra": c.execute("SELECT count(*) FROM spectra").fetchone()[0],
              "queues": {x[0]: x[1] for x in c.execute("SELECT queue,count(*) FROM molecules GROUP BY queue")},
              "candidates": {x[0]: x[1] for x in c.execute("SELECT state,count(*) FROM candidates GROUP BY state")},
              "progress": [dict(x) for x in c.execute("SELECT * FROM progress")],
              "issues": c.execute("SELECT count(*) FROM import_issues").fetchone()[0]}
    from .ai_database import statistics as annotation_statistics
    result['annotation_origins']=annotation_statistics(c)
    result["spectrum_quality"] = {x[0]: x[1] for x in c.execute("SELECT decision,count(*) FROM spectrum_quality GROUP BY decision")}
    result['dual_assignment_progress'] = dual_assignment_progress(result['queues'])
    result['annotation_processing'] = {status: {'pending': 0, 'finished': 0} for status in ('full', 'partial', 'none')}
    for status, done, count in c.execute(f'SELECT ({COVERAGE_SQL}),task_done,count(*) FROM molecules GROUP BY ({COVERAGE_SQL}),task_done'):
        result['annotation_processing'][status]['finished' if done else 'pending'] = count
    result['spectrum_coverage'] = spectrum_coverage(c, result['queues'])
    return result


def spectrum_coverage(c, queues):
    """Count molecules by their current C/H peaks, not by spectrum versions.

    The maintained dual-spectrum queues already encode coverage. Only the much
    smaller single-spectrum/held groups need their per-nucleus counters read.
    This avoids rescanning the large structure JSON of every molecule.
    """
    counts = {kind: {status: 0 for status in ('full', 'partial', 'none')}
              for kind in ('dual', 'c_only', 'h_only', 'no_current')}
    dual_queues = {'complete': 'full', 'complete_remaining': 'partial', 'both_unassigned': 'none'}
    for queue, status in dual_queues.items():
        counts['dual'][status] = queues.get(queue, 0)
    counts['no_current']['none'] = queues.get('quality_hold', 0)
    remaining = [queue for queue in queues if queue not in (*dual_queues, 'quality_hold')]
    if remaining:
        kind = "CASE WHEN c_count>0 AND h_count>0 THEN 'dual' WHEN c_count>0 THEN 'c_only' WHEN h_count>0 THEN 'h_only' ELSE 'no_current' END"
        placeholders = ','.join('?' for _ in remaining)
        for row in c.execute(f'''SELECT ({kind}),({COVERAGE_SQL}),count(*) FROM molecules
                               INDEXED BY molecule_queue WHERE queue IN ({placeholders})
                               GROUP BY ({kind}),({COVERAGE_SQL})''', remaining):
            counts[row[0]][row[1]] += row[2]
    rows = [{'type': kind, **values, 'total': sum(values.values())} for kind, values in counts.items()]
    total = {status: sum(row[status] for row in rows) for status in ('full', 'partial', 'none', 'total')}
    return {'rows': rows, 'total': total, 'unit': 'molecule',
            'basis': 'current_retained_peak_counts',
            'note': '按分子去重，依据当前谱的峰数和已保存归属状态统计。隔离、停用的谱以及参考、历史版本不计入当前谱；分子层面停用或排除训练仍保留在存档分类中。没有当前峰的分子归入无当前可用谱，并在无归属列计数，不表示适合继续标注。'}


def dual_assignment_progress(queues):
    # refresh_molecule maintains these queues from ACTIVE spectra and peak
    # assignments in the same transaction as edits. Reuse its indexed counts;
    # never rescan 1M peak JSON documents on each UI save.
    ready = queues.get('complete', 0)
    total = ready + queues.get('complete_remaining', 0) + queues.get('both_unassigned', 0)
    return {'ready_molecules': ready, 'total_dual_molecules': total,
            'remaining_molecules': total - ready,
            'percent': ready / total * 100 if total else 0.0,
            'unit': 'molecule', 'basis': 'all_active_reported_peaks_assigned_v1',
            'note': '同一分子的当前未隔离、未停用 C/H 谱均非空，且全部纳入峰表的峰有归属；包含等价、重叠及原库多原子组。歧义、未处理及排除峰不算全归属。保守口径：多条当前谱需全部归属。不等于谱完整性、条件一致性或归属正确性认证。'}


def molecule_svg(smiles, selected=(), view='auto'):
    if view not in ('auto', 'full', 'ligands'):
        raise ValueError('结构视图无效。')
    mol = Chem.MolFromSmiles(smiles)
    metals = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in METAL_ATOMIC_NUMBERS}
    # Depiction only: never change the stored structure, atom order or labels.
    if metals and view in ('auto', 'ligands'):
        drawing = Chem.RWMol(mol)
        for b in list(mol.GetBonds()):
            if b.GetBeginAtomIdx() in metals or b.GetEndAtomIdx() in metals:
                drawing.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        mol = drawing.GetMol()
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)
        rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(1100 if metals else 720, 700 if metals else 440)
    options = drawer.drawOptions()
    options.addAtomIndices = not bool(metals)
    options.baseFontSize = 0.85
    if metals:
        options.includeRadicals = False
        options.padding = .08
        for atom in mol.GetAtoms():
            options.atomLabels[atom.GetIdx()] = atom.GetSymbol() + str(atom.GetIdx())
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=list(selected))
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
