"""Compare mapped experimental spectra and commit literature assignments to AtomNMR."""

from collections import Counter, defaultdict
from contextlib import closing
from pathlib import Path
from urllib.parse import unquote
import json

from annotation_tool_v3 import core as db
from annotation_tool_v3.peak_table import effective_peaks
from annotation_tool_v3.spectrum_revision import KEY as REVISIONS

from .dedup import normalized
from .extract import SOURCE, spectrum_payload

CHANGES = {'added', 'filled'}


def initialize(conn):
    conn.executescript("""
    CREATE TABLE IF NOT EXISTS literature_imports(
      job_id TEXT NOT NULL, source_key TEXT NOT NULL, nucleus TEXT NOT NULL,
      status TEXT NOT NULL CHECK(status IN ('added','filled','duplicate','review')),
      molecule_id TEXT REFERENCES molecules(id), spectrum_id TEXT REFERENCES spectra(id),
      new_molecule INTEGER NOT NULL, assigned_peaks INTEGER NOT NULL,
      reason TEXT NOT NULL, citation_json TEXT NOT NULL, created_at TEXT NOT NULL,
      PRIMARY KEY(job_id,source_key)
    );
    CREATE INDEX IF NOT EXISTS literature_import_molecule ON literature_imports(molecule_id);
    """)


def number(value):
    try:
        return normalized(float(value))
    except (TypeError, ValueError):
        return normalized(value)


def solvent(value):
    label = normalized(value).replace('-', '').replace(' ', '')
    return {'chloroformd': 'cdcl3', 'deuteratedchloroform': 'cdcl3',
            'dmso6': 'dmsod6', 'methanold4': 'cd3od', 'meodd4': 'cd3od',
            'acetoned6': '(cd3)2co', 'c3d6o': '(cd3)2co',
            'benzened6': 'c6d6', 'deuteriumoxide': 'd2o',
            'unknown': '', '未知': '', 'mixed': '', 'n/a': ''}.get(label, label)


def doi(provenance):
    value = unquote(normalized(provenance.get('doi') or provenance.get('doi_encoded')))
    return value.removeprefix('https://doi.org/').removeprefix('http://doi.org/').removeprefix('doi:').strip()


def condition_relation(old, incoming, old_source, incoming_source):
    left, right = solvent(old.get('solvent')), solvent(incoming.get('solvent'))
    if left and right and left != right:
        return 'different'
    for key, normalize in (('frequency_mhz', number), ('temperature_k', number),
                           ('chemical_shift_reference', normalized), ('sample', normalized)):
        a = normalize(old.get(key) or old_source.get(key))
        b = normalize(incoming.get(key) or incoming_source.get(key))
        if a and b and a != b:
            return 'different'
    same_article = any(a and a == b for a, b in ((doi(old_source), doi(incoming_source)),
                       (old_source.get('pdf_sha256'), incoming_source.get('pdf_sha256'))))
    return 'same' if left and right or same_article else 'unknown'


def shift_key(peak):
    if peak.get('range_min') is not None and peak.get('range_max') is not None:
        return tuple(number(peak[key]) for key in ('range_min', 'range_max'))
    return (number(peak.get('shift_ppm')),)


def peak_buckets(peaks):
    result = defaultdict(list)
    for peak in peaks:
        result[shift_key(peak)].append(peak)
    return result


def peak_metadata_agrees(left, right):
    def couplings(values):
        seq = values if isinstance(values, (tuple, list)) else [values]
        return sorted(number(v) for v in seq if v is not None and v != '')

    for key, normalize in (('reported_nH', number), ('multiplicity', normalized), ('j_values', couplings)):
        a, b = normalize(left.get(key)), normalize(right.get(key))
        if a and b and a != b:
            return False
    return True


def assignment_groups(provenance, peaks):
    shifts = {peak['peak_id']: shift_key(peak) for peak in peaks}
    groups = provenance.get(REVISIONS, {}).get('assignment_groups', [])
    covered, signatures = set(), Counter()
    for group in groups:
        ids = group['peak_ids']
        if not all(pid in shifts for pid in ids):
            continue
        covered.update(ids)
        signatures[(tuple(sorted(shifts[pid] for pid in ids)), tuple(sorted(group['atom_indices'])),
                    group['relation'])] += 1
    return covered, signatures


def assigned_atoms(assignment):
    if assignment and assignment['status'] in db.ASSIGNED:
        return tuple(sorted(set(assignment['atom_indices'])))
    return ()


def align(conn, spectrum, old_peaks, old_source, peaks, source, assignments):
    """Return additions for an exact observation match, or a review reason."""
    old_groups, old_group_keys = assignment_groups(old_source, old_peaks)
    new_groups, new_group_keys = assignment_groups(source, peaks)
    if old_group_keys != new_group_keys:
        return None, '可互换或歧义归属组与主库不同，须核对原表'
    old = {a['peak_id']: {**dict(a), 'atom_indices': json.loads(a['atoms_json'])}
           for a in conn.execute('SELECT * FROM annotations WHERE spectrum_id=?', (spectrum['id'],))}
    left, right = peak_buckets(old_peaks), peak_buckets(peaks)
    additions = {}
    for shift, incoming in right.items():
        existing = left[shift]
        if len(incoming) > 1:
            def signatures(rows, labels, grouped):
                return Counter((normalized(p.get('multiplicity')), number(p.get('reported_nH')),
                                db.dumps(p.get('j_values') or []), p['peak_id'] in grouped,
                                assigned_atoms(labels.get(p['peak_id']))) for p in rows)

            if signatures(existing, old, old_groups) != signatures(incoming, assignments, new_groups):
                return None, '有重复位移，无法唯一对应到主库峰'
            continue
        previous, peak = existing[0], incoming[0]
        if not peak_metadata_agrees(previous, peak):
            return None, '同位移峰的峰型、氢数或 J 值与主库不一致'
        old_id, new_id = previous['peak_id'], peak['peak_id']
        if (old_id in old_groups) != (new_id in new_groups):
            return None, '峰所在的归属组与主库不一致'
        if new_id in new_groups:
            continue
        incoming_label = assignments.get(new_id)
        targets = assigned_atoms(incoming_label)
        if not targets:
            continue
        current = old.get(old_id)
        prior_targets = assigned_atoms(current)
        if prior_targets == targets:
            continue
        if current and (current['status'] != 'unassigned' or current['note'] or current['origin'].startswith('manual')):
            return None, '论文归属与主库现有标签或人工处理状态冲突'
        additions[old_id] = (incoming_label, peak['literature'])
    return additions, ''


def candidate(conn, mid, nucleus, conditions, peaks, provenance):
    if not mid:
        return None, ''
    pattern = Counter(map(shift_key, peaks))
    matches = []
    for row in conn.execute('SELECT *, (' + db.ACTIVE + ') AS active FROM spectra WHERE molecule_id=? AND nucleus=?',
                            (mid, nucleus)):
        old_source = json.loads(row['source_json'])
        old_conditions = {**dict(row), **json.loads(row['conditions_json'])}
        relation = condition_relation(old_conditions, conditions, old_source, provenance)
        if relation == 'different':
            continue
        old_peaks = effective_peaks(json.loads(row['peaks_json']), old_source)
        if Counter(map(shift_key, old_peaks)) == pattern:
            matches.append((row, old_peaks, old_source, relation))
    if len(matches) > 1:
        return None, '多条主库谱的完整位移峰表相同，无法唯一对应'
    if not matches:
        return None, ''
    row, old_peaks, old_source, relation = matches[0]
    if not row['active'] or relation == 'unknown':
        return None, ('相同峰表的主库谱已停用或隔离' if not row['active']
                      else '完整位移峰表相同，但溶剂和文献来源不足以确认同一实验')
    return (row, old_peaks, old_source), ''


class LiteratureImporter:
    def __init__(self, path, on_change=None):
        self.path = path
        self.on_change = on_change
        with closing(db.connect(path)) as conn:
            initialize(conn)

    def ingest(self, job, table):
        table = {**table, 'mapping': {int(k): value for k, value in table['mapping'].items()}}
        entries = {}
        with closing(db.connect(self.path, timeout=15)) as conn, conn:
            conn.execute('BEGIN IMMEDIATE')
            molecule = conn.execute('SELECT id FROM molecules WHERE smiles=?', (table['smiles'],)).fetchone()
            mid = molecule['id'] if molecule else None
            pair_ids, changed = {}, False
            for nucleus, rows in table['peaks'].items():
                if not rows:
                    continue
                source_key, conditions, peaks, provenance, labels = spectrum_payload(table, nucleus, rows)
                receipt = conn.execute('SELECT * FROM literature_imports WHERE job_id=? AND source_key=?',
                                       (job, source_key)).fetchone()
                if receipt:
                    entries[nucleus] = dict(receipt)
                    if receipt['status'] != 'review' and receipt['spectrum_id']:
                        pair_ids[nucleus] = receipt['spectrum_id']
                    continue
                provenance['literature_import_job'] = job
                citation = {**provenance, 'source_key': source_key, 'conditions': conditions}
                status, sid, added, new_molecule = 'review', None, 0, 0
                existing_source = conn.execute('SELECT *, (' + db.ACTIVE + ') AS active '
                                               'FROM spectra WHERE source_key=?', (source_key,)).fetchone()
                target, reason = None, ''
                if existing_source:
                    if existing_source['molecule_id'] != mid or existing_source['nucleus'] != nucleus:
                        reason = '同一文献来源键关联了不同的结构或核种'
                    elif not existing_source['active']:
                        reason = '同一原文谱已在主库中停用或隔离'
                    else:
                        old_source = json.loads(existing_source['source_json'])
                        old_peaks = effective_peaks(json.loads(existing_source['peaks_json']), old_source)
                        if Counter(map(shift_key, old_peaks)) != Counter(map(shift_key, peaks)):
                            reason = '同一原文谱的位移峰表和当前版本不一致'
                        else:
                            target = (existing_source, old_peaks, old_source)
                else:
                    target, reason = candidate(conn, mid, nucleus, conditions, peaks, provenance)
                if target:
                    previous, old_peaks, old_source = target
                    sid = previous['id']
                    additions, reason = align(conn, previous, old_peaks, old_source, peaks, provenance, labels)
                    if additions is not None:
                        status = 'filled' if additions else 'duplicate'
                        reason = '已补入主库谱的空白归属' if additions else '主库已保存同一峰表及其归属'
                        if additions:
                            self.fill(conn, sid, additions, old_source, citation)
                            added = len(additions)
                elif not reason and not existing_source:
                    new_molecule = int(mid is None)
                    mid = mid or db.add_molecule(conn, table['smiles'], SOURCE)
                    sid = db.add_spectrum(conn, mid, nucleus, SOURCE, source_key, conditions, peaks, provenance, labels)
                    status, reason = 'added', '已新增文献实验谱'
                    added = sum(bool(assigned_atoms(a)) for a in labels.values())
                receipt = dict(job_id=job, source_key=source_key, nucleus=nucleus, status=status, molecule_id=mid,
                               spectrum_id=sid, new_molecule=new_molecule, assigned_peaks=added, reason=reason,
                               citation_json=db.dumps(citation), created_at=db.now())
                conn.execute('''INSERT INTO literature_imports VALUES(
                    :job_id,:source_key,:nucleus,:status,:molecule_id,:spectrum_id,:new_molecule,:assigned_peaks,
                    :reason,:citation_json,:created_at)''', receipt)
                entries[nucleus] = receipt
                changed |= status in CHANGES
                if status != 'review' and sid:
                    pair_ids[nucleus] = sid
            if changed:
                if len(pair_ids) == 2:
                    db.add_pair(conn, mid, pair_ids['13C'], pair_ids['1H'], 'literature_same_table')
                db.refresh_molecule(conn, mid)
                conn.execute('UPDATE molecules SET revision=revision+1 WHERE id=?', (mid,))
        if changed and self.on_change:
            self.on_change()
        return {nucleus: {k: receipt[k] for k in ('status', 'reason', 'molecule_id', 'spectrum_id')}
                for nucleus, receipt in entries.items()}

    @staticmethod
    def fill(conn, sid, additions, old_source, citation):
        provenance = dict(old_source)
        sources = provenance.setdefault('literature_assignment_sources', [])
        sources.append({**citation, 'peak_rows_by_id': {pid: row for pid, (_, row) in additions.items()}})
        conn.execute('UPDATE spectra SET source_json=? WHERE id=?', (db.dumps(provenance), sid))
        conn.executemany('''INSERT INTO annotations(spectrum_id,peak_id,status,atoms_json,note,origin,updated_at)
            VALUES(?,?,?,?,?,?,?) ON CONFLICT(spectrum_id,peak_id) DO UPDATE SET status=excluded.status,
            atoms_json=excluded.atoms_json,note=excluded.note,origin=excluded.origin,updated_at=excluded.updated_at''',
            [(sid, pid, a['status'], db.dumps(a['atom_indices']), a['note'], SOURCE, db.now())
             for pid, (a, _) in additions.items()])

    def statistics(self):
        with closing(db.connect(self.path, timeout=5)) as conn:
            conn.execute('BEGIN')
            by_nucleus = [dict(r) for r in conn.execute('''SELECT nucleus,status,count(*) AS spectra,
                sum(assigned_peaks) AS assigned_peaks FROM literature_imports GROUP BY nucleus,status''')]
            affected = conn.execute("""SELECT count(DISTINCT molecule_id) FROM literature_imports
                WHERE status IN ('added','filled')""").fetchone()[0]
            new = conn.execute('SELECT coalesce(sum(new_molecule),0) FROM literature_imports').fetchone()[0]
            batches = {r['job_id']: dict(r) for r in conn.execute('''SELECT job_id, count(*) AS compared,
                sum(status='added') AS added,sum(status='filled') AS filled,sum(status='duplicate') AS duplicates,
                sum(status='review') AS review FROM literature_imports GROUP BY job_id''')}
        return dict(by_nucleus=by_nucleus, new_molecules=new, affected_molecules=affected, batches=batches)
