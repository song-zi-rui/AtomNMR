"""Curated carbon codes and unordered peak-set/atom-set correspondence.

Raw peak tables are immutable. Edits live in provenance and curation events;
unordered groups must never become a cross-product of hard peak labels.
"""
from __future__ import annotations

import copy
import json
import re
import uuid

from . import core as db

KEY = 'spectrum_revisions'
HCOUNT = {'S': 0, 'D': 1, 'T': 2, 'Q': 3}


def edits(provenance):
    return provenance.get(KEY) or {'signal_codes': {}, 'assignment_groups': []}


def group_for(provenance, peak_id):
    return next((g for g in edits(provenance).get('assignment_groups', [])
                 if peak_id in g['peak_ids']), None)


def original_code(peak):
    return peak.get('source_signal_code') or ''


def effective_code(peak, letter):
    raw = original_code(peak)
    return re.sub(r'[SDTQ]$', letter, raw) if re.search(r'[SDTQ]$', raw) else raw + letter


def project(spectrum, atoms):
    """Attach effective metadata without rewriting raw signals or audit fields."""
    current = edits(spectrum['provenance'])
    groups = current.get('assignment_groups', [])
    spectrum['assignment_groups'] = copy.deepcopy(groups)
    for g in spectrum['assignment_groups']:
        g.pop('before_assignments', None)
        g.pop('applied_annotations', None)
    for p in spectrum['peaks']:
        override = current.get('signal_codes', {}).get(p['peak_id'])
        if override:
            p['signal_code_revision'] = copy.deepcopy(override)
            p['effective_signal_code'] = effective_code(p, override['code'])
            targets = p['assignment']['atom_indices'] if p['assignment']['status'] in db.ASSIGNED else []
            group = group_for(spectrum['provenance'], p['peak_id'])
            if group:
                targets = group['atom_indices']
            p['effective_signal_code_eligible'] = bool(targets) and all(
                type(i) is int and 0 <= i < len(atoms)
                and atoms[i]['symbol'] == 'C' and atoms[i]['total_hydrogens'] == HCOUNT[override['code']]
                for i in targets)
            p['signal_code_revision_check'] = ('consistent_with_current_atoms'
                if p['effective_signal_code_eligible'] else 'no_confirmed_atoms_or_hcount_conflict')
        group = group_for(spectrum['provenance'], p['peak_id'])
        if group:
            p['assignment_group_id'] = group['id']
            p['assignment_group'] = next(g for g in spectrum['assignment_groups'] if g['id'] == group['id'])
            # Even a later/imported hard label cannot silently override a set relation.
            p['assignment'] = {**p['assignment'], 'status': 'ambiguous', 'atom_indices': []}


def require_editable(c, mid, spec):
    from .curation import state_for
    sid = spec['id']
    state = state_for(c, 'spectrum', sid)
    q = c.execute('SELECT decision FROM spectrum_quality WHERE spectrum_id=?', (sid,)).fetchone()
    if state_for(c, 'molecule', mid) != 'active':
        raise ValueError('请先恢复分子使用状态，再修订谱。')
    if not (state == 'active' and (not q or q[0] == 'passed')
            or state == 'reference' and spec['source'] == 'manual_revision'
            and not json.loads(spec['source_json']).get('dataset_version')
            and not json.loads(spec['source_json']).get('reference_only_after_revision')):
        raise ValueError('只能编辑当前可用谱，或尚未采用的归属修订副本。历史版本请先新建副本。')
def mutate(c, mid, spec, payload, reason, atoms):
    require_editable(c, mid, spec)
    sid = spec['id']
    provenance = json.loads(spec['source_json'])
    current = copy.deepcopy(edits(provenance))
    current.setdefault('signal_codes', {})
    current.setdefault('assignment_groups', [])
    from .peak_table import effective_peaks
    peaks = {p['peak_id']: p for p in effective_peaks(json.loads(spec['peaks_json']), provenance)}
    name = payload['action']
    detail = {}
    if name == 'signal_code':
        pid, code = payload.get('peak_id'), payload.get('code')
        if spec['nucleus'] != '13C' or pid not in peaks or code not in (None, *HCOUNT):
            raise ValueError('请选择当前 C 谱的峰，以及 S / D / T / Q 或恢复原码。')
        before = current['signal_codes'].get(pid)
        if code is None:
            current['signal_codes'].pop(pid, None)
        else:
            current['signal_codes'][pid] = {'code': code, 'original_signal_code': original_code(peaks[pid]),
                                          'reason': reason, 'updated_at': db.now()}
        detail = {'peak_id': pid, 'before': before, 'after': current['signal_codes'].get(pid)}
    elif name == 'assignment_group':
        pids, targets = payload.get('peak_ids'), payload.get('atom_indices')
        if (not isinstance(pids, list) or len(pids) < 2 or any(not isinstance(p, str) for p in pids)
                or len(set(pids)) != len(pids) or any(p not in peaks for p in pids)):
            raise ValueError('请选择同一谱内至少两个不同的位移记录。')
        if (not isinstance(targets, list) or len(targets) < 2
                or any(type(i) is not int or not 0 <= i < len(atoms) for i in targets)
                or len(set(targets)) != len(targets)):
            raise ValueError('请选择至少两个不同且有效的原子编号。')
        if any(group_for({KEY: current}, pid) for pid in pids):
            raise ValueError('所选位移已属于其他位移组，请先解除原组。')
        for i in targets:
            if (spec['nucleus'] == '13C' and atoms[i]['symbol'] != 'C'
                    or spec['nucleus'] == '1H' and not atoms[i]['total_hydrogens']):
                raise ValueError('原子不适合当前核种；C 谱选碳，H 谱选带氢重原子。')
        group = {'id': 'G' + uuid.uuid4().hex[:24], 'peak_ids': pids, 'atom_indices': sorted(targets),
                 'relation': 'unordered_set_correspondence', 'within_group_order': 'unknown',
                 'chemical_equivalence_asserted': False, 'reason': reason, 'created_at': db.now(),
                 'before_assignments': [], 'applied_annotations': []}
        for pid in pids:
            old = c.execute('SELECT * FROM annotations WHERE spectrum_id=? AND peak_id=?', (sid, pid)).fetchone()
            group['before_assignments'].append({'peak_id': pid, 'annotation': dict(old) if old else None})
            c.execute('INSERT INTO annotations VALUES(?,?,?,?,?,?,?) ON CONFLICT(spectrum_id,peak_id) DO UPDATE SET status=excluded.status,atoms_json=excluded.atoms_json,note=excluded.note,origin=excluded.origin,updated_at=excluded.updated_at',
                      (sid, pid, 'ambiguous', '[]', '位移组 ' + group['id'] + '；组内顺序未定。', 'manual', db.now()))
            group['applied_annotations'].append(dict(c.execute('SELECT * FROM annotations WHERE spectrum_id=? AND peak_id=?', (sid, pid)).fetchone()))
        current['assignment_groups'].append(group)
        detail = {'group': copy.deepcopy(group), 'warning': 'set relation, not all-to-all hard assignments'}
    elif name == 'remove_assignment_group':
        group = next((g for g in current['assignment_groups'] if g['id'] == payload.get('group_id')), None)
        if group is None:
            raise ValueError('位移组不存在或已经解除。')
        before = {r['peak_id']: r['annotation'] for r in group.get('before_assignments', [])}
        preserved = 0
        for expected in group.get('applied_annotations', []):
            pid = expected['peak_id']
            actual = c.execute('SELECT * FROM annotations WHERE spectrum_id=? AND peak_id=?', (sid, pid)).fetchone()
            if actual is None or dict(actual) != expected:
                preserved += 1
                continue
            restore = before[pid]
            if restore is None:
                c.execute('DELETE FROM annotations WHERE spectrum_id=? AND peak_id=?', (sid, pid))
            else:
                c.execute('UPDATE annotations SET status=?,atoms_json=?,note=?,origin=?,updated_at=? WHERE spectrum_id=? AND peak_id=?',
                          (restore['status'], restore['atoms_json'], restore['note'], restore['origin'], restore['updated_at'], sid, pid))
        current['assignment_groups'] = [g for g in current['assignment_groups'] if g['id'] != group['id']]
        detail = {'removed_group': group, 'later_edits_preserved': preserved}
    else:
        raise ValueError('未知谱修订操作。')
    provenance[KEY] = current
    c.execute('UPDATE spectra SET source_json=? WHERE id=?', (db.dumps(provenance), sid))
    return detail
