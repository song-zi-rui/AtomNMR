"""Explicit, reversible curation. Automatic screening and human decisions stay separate."""
from __future__ import annotations

import json
import uuid

from . import core as db


def state_for(c, scope, target):
    row = c.execute("SELECT state FROM curation_state WHERE scope=? AND target_id=?", (scope, target)).fetchone()
    return row[0] if row else 'active'


def summary(c, mid):
    states = [dict(x) for x in c.execute("SELECT * FROM curation_state WHERE molecule_id=?", (mid,))]
    return {'molecule_state': next((x['state'] for x in states if x['scope'] == 'molecule'), 'active'),
            'states': states,
            'open_issues': c.execute("SELECT count(*) FROM curation_issues WHERE molecule_id=? AND status='open'", (mid,)).fetchone()[0]}


def details(c, mid):
    molecule = c.execute('SELECT revision FROM molecules WHERE id=?', (mid,)).fetchone()
    if not molecule:
        raise KeyError(mid)
    result = {**summary(c, mid), 'revision': molecule[0], 'spectra': []}
    for row in c.execute('SELECT * FROM spectra WHERE molecule_id=? ORDER BY nucleus,id', (mid,)):
        s = dict(row)
        for field, dest in [('peaks_json', 'peaks'), ('conditions_json', 'conditions'), ('source_json', 'provenance')]:
            s[dest] = json.loads(s.pop(field))
        from .peak_table import effective_peaks
        s['peaks'] = effective_peaks(s['peaks'], s['provenance'])
        s['manual_state'] = state_for(c, 'spectrum', s['id'])
        q = c.execute('SELECT details_json FROM spectrum_quality WHERE spectrum_id=?', (s['id'],)).fetchone()
        s['automatic_quality'] = json.loads(q[0]) if q else {'decision': 'not_checked'}
        from .spectrum_versions import annotate_role
        annotate_role(s)
        complete = c.execute('SELECT status,reason FROM spectrum_completeness WHERE spectrum_id=?', (s['id'],)).fetchone()
        s['completeness'] = dict(complete) if complete else {'status': 'unknown', 'reason': ''}
        s['annotations'] = [{**dict(a), 'atom_indices': json.loads(a['atoms_json'])} for a in c.execute('SELECT * FROM annotations WHERE spectrum_id=?', (s['id'],))]
        from .spectrum_revision import project
        assignments = {a['peak_id']: a for a in s['annotations']}
        for p in s['peaks']:
            p['assignment'] = assignments.get(p['peak_id'], {'status': 'unassigned', 'atom_indices': []})
        project(s, json.loads(c.execute('SELECT structure_json FROM molecules WHERE id=?', (mid,)).fetchone()[0])['atoms'])
        # Raw downloaded files remain in SQLite, not repeatedly sent to the browser.
        s['provenance'].pop('raw_star', None)
        s['provenance'].pop('raw_mol', None)
        result['spectra'].append(s)
    result['spectra'].sort(key=lambda s: (s['nucleus'], s['manual_state'] != 'active', not bool(s.get('version_info')),
                                        -(s.get('version_info') or {}).get('number', 0), s['id']))
    result['issues'] = [dict(x) for x in c.execute('SELECT * FROM curation_issues WHERE molecule_id=? ORDER BY id DESC', (mid,))]
    result['events'] = [{**dict(x), 'details': json.loads(x['details_json'])} for x in c.execute('SELECT * FROM curation_events WHERE molecule_id=? ORDER BY id DESC LIMIT 50', (mid,))]
    return result


def put_state(c, mid, scope, target, state, reason):
    c.execute('INSERT INTO curation_state VALUES(?,?,?,?,?,?) ON CONFLICT(scope,target_id) DO UPDATE SET state=excluded.state,reason=excluded.reason,updated_at=excluded.updated_at',
              (scope, target, mid, state, reason, db.now()))


def action(c, mid, payload, reference=None):
    """Every mutation is revision guarded and short; fetching occurs before this call."""
    try:
        c.execute('BEGIN IMMEDIATE')
        mol = c.execute('SELECT * FROM molecules WHERE id=?', (mid,)).fetchone()
        if not mol:
            raise KeyError(mid)
        if payload.get('revision') != mol['revision']:
            raise RuntimeError('记录已更新。请重新打开数据管理，核对后操作。')
        reason = str(payload.get('reason', '')).strip()
        name = payload.get('action')
        if not reason and name == 'disposition' and payload.get('state') == 'active':
            reason = '用户恢复人工停用 / 排除状态；保留此前操作记录。'
        if not reason or len(reason) > 5000:
            raise ValueError('请填写具体原因（1–5000 字）。')
        sid = payload.get('spectrum_id')
        spec = c.execute('SELECT * FROM spectra WHERE id=? AND molecule_id=?', (sid, mid)).fetchone() if sid else None
        if sid and not spec:
            raise ValueError('谱不属于当前分子。')
        event = {'reason': reason, 'spectrum_id': sid}
        from .spectrum_versions import final_for
        if name == 'disposition':
            scope, target = ('spectrum', sid) if sid else ('molecule', mid)
            state = payload.get('state')
            if state not in ('active', 'held', 'excluded'):
                raise ValueError('处理状态无效。')
            previous = state_for(c, scope, target)
            if previous in ('reference', 'superseded'):
                raise ValueError('参考或历史版本请使用“采用此版本”，不能直接恢复。')
            if state == 'active' and spec:
                q = c.execute("SELECT decision FROM spectrum_quality WHERE spectrum_id=?", (sid,)).fetchone()
                if q and q[0] != 'passed':
                    raise ValueError('自动隔离仍生效，不能一键清除。请先查明问题，再新建修订版本。')
            put_state(c, mid, scope, target, state, reason)
            event.update(before=previous, after=state, scope=scope)
        elif name == 'issue':
            category = payload.get('category')
            if category not in ('missing_peaks', 'wrong_assignment', 'structure', 'conditions', 'source', 'other'):
                raise ValueError('问题类型无效。')
            c.execute("INSERT INTO curation_issues(molecule_id,spectrum_id,category,details,status,created_at,updated_at) VALUES(?,?,?,?,'open',?,?)", (mid, sid, category, reason, db.now(), db.now()))
            event['category'] = category
        elif name == 'resolve_issue':
            issue = c.execute("SELECT * FROM curation_issues WHERE id=? AND molecule_id=?", (payload.get('issue_id'), mid)).fetchone()
            if not issue or issue['status'] != 'open':
                raise ValueError('没有该待处理问题。')
            c.execute("UPDATE curation_issues SET status='resolved',updated_at=? WHERE id=?", (db.now(), issue['id']))
            event['issue_id'] = issue['id']
        elif name == 'completeness':
            status = payload.get('status')
            if not spec or status not in ('unknown', 'incomplete', 'complete'):
                raise ValueError('请选择谱和完整性状态。')
            c.execute('INSERT INTO spectrum_completeness VALUES(?,?,?,?) ON CONFLICT(spectrum_id) DO UPDATE SET status=excluded.status,reason=excluded.reason,updated_at=excluded.updated_at', (sid, status, reason, db.now()))
            event['status'] = status
        elif name == 'clone' and spec and (existing_final := final_for(c, spec)) is not None:
            # The second slot is reused, never create additional intermediate drafts.
            event.update(new_spectrum_ids=[existing_final['id']], reused_existing_final=True)
        elif name == 'clone':
            if not spec:
                raise ValueError('请选择要修订的原谱。')
            # Copy, never update, the original experiment and its human labels.
            peaks = json.loads(spec['peaks_json'])
            visible = next(s for s in db.record(c, mid, include_archived=True)['spectra'] if s['id'] == sid)
            assignments = {p['peak_id']: {'status': p['assignment']['status'], 'atom_indices': p['assignment']['atom_indices'], 'note': p['assignment'].get('note', '')} for p in visible['peaks']}
            provenance = json.loads(spec['source_json'])
            provenance.pop('dataset_version', None)
            provenance.pop('reference_only_after_revision', None)
            provenance.update(parent_spectrum_id=sid, revision_reason=reason, requires_explicit_adoption=True)
            # Inherited groups retain set meaning; undo in a copy must not replay
            # snapshots belonging to the parent spectrum.
            for group in provenance.get('spectrum_revisions', {}).get('assignment_groups', []):
                group['parent_group_id'] = group['id']
                group['id'] = 'G' + uuid.uuid4().hex[:24]
                group['before_assignments'] = []
                group['applied_annotations'] = []
            quality = c.execute('SELECT details_json FROM spectrum_quality WHERE spectrum_id=?', (sid,)).fetchone()
            if quality:
                provenance['parent_automatic_quality'] = json.loads(quality[0])
            elif spec['source'] == 'NMRexp':
                provenance['parent_automatic_quality'] = {'decision': 'not_checked'}
            provenance['warnings'] = [*provenance.get('warnings', []), '这是原谱及已有修订的副本；原谱的来源或峰表问题仍需核实，复制本身不会解决这些问题。']
            new = db.add_spectrum(c, mid, spec['nucleus'], 'manual_revision', uuid.uuid4().hex, json.loads(spec['conditions_json']), peaks, provenance, assignments)
            complete = c.execute('SELECT status,reason FROM spectrum_completeness WHERE spectrum_id=?', (sid,)).fetchone()
            if complete:
                c.execute('INSERT INTO spectrum_completeness VALUES(?,?,?,?)', (new, complete[0], complete[1], db.now()))
            put_state(c, mid, 'spectrum', new, 'reference', reason)
            event['new_spectrum_ids'] = [new]
        elif name == 'promote_existing_revision':
            if not spec:
                raise ValueError('请选择已有修订的谱。')
            provenance = json.loads(spec['source_json'])
            revisions = provenance.get('spectrum_revisions') or {}
            if provenance.get('dataset_version') or not (any(x.get('active', True) for x in provenance.get('peak_table_supplements', []))
                    or revisions.get('signal_codes') or revisions.get('assignment_groups')):
                raise ValueError('该谱没有需要转换的旧修订，或已经是新版本。')
            from .spectrum_versions import create_current
            spec, version_event = create_current(c, mid, spec, reason, name)
            event.update(version_event, spectrum_id=spec['id'])
        elif name in ('supplement_peak', 'withdraw_peak_supplement', 'edit_peak', 'delete_peak',
                      'signal_code', 'assignment_group', 'remove_assignment_group'):
            if not spec:
                raise ValueError('请选择要修订的谱。')
            from .spectrum_versions import create_current
            spec, version_event = create_current(c, mid, spec, reason, name)
            event.update(version_event, spectrum_id=spec['id'])
            if name in ('signal_code', 'assignment_group', 'remove_assignment_group'):
                from .spectrum_revision import mutate
                event.update(mutate(c, mid, spec, payload, reason, json.loads(mol['structure_json'])['atoms']))
            else:
                from .peak_table import mutate
                event.update(mutate(c, mid, spec, payload, reason))
        elif name == 'import_bmrb':
            if not reference or reference['local_smiles'] != mol['smiles']:
                raise ValueError('参考文件或结构验证失败。')
            event['new_spectrum_ids'] = []
            for s in reference['spectra']:
                key = f"reference:{mid}:{reference['entry_id']}:{reference['sha256']}:{s['list_id']}:{s['nucleus']}"
                existing = c.execute('SELECT id FROM spectra WHERE source_key=?', (key,)).fetchone()
                if existing:
                    event['new_spectrum_ids'].append(existing[0])
                    continue
                new = db.add_spectrum(c, mid, s['nucleus'], 'BMRB_reference', key, s['conditions'], s['peaks'], reference['provenance'], s['assignments'])
                put_state(c, mid, 'spectrum', new, 'reference', reason)
                event['new_spectrum_ids'].append(new)
            event['entry_id'] = reference['entry_id']
        elif name == 'adopt':
            if not spec or state_for(c, 'spectrum', sid) not in ('reference', 'superseded'):
                raise ValueError('只能采用参考版本或历史版本。')
            if payload.get('acknowledge') is not True:
                raise ValueError('请确认已比较来源、结构、实验条件和归属；采用不代表补齐缺峰。')
            parent_quality = json.loads(spec['source_json']).get('parent_automatic_quality')
            if parent_quality and parent_quality.get('decision') != 'passed':
                raise ValueError('修订副本继承了未检查 / 自动隔离的来源问题。复制或改归属不能解决结构、来源、条件或峰表冲突；请引入通过结构检查的新实验参考。')
            old_id = payload.get('replace_spectrum_id')
            from .spectrum_versions import root_id
            root = root_id(c, spec)
            family_current = [s['id'] for s in c.execute('SELECT * FROM spectra WHERE molecule_id=? AND nucleus=?', (mid, spec['nucleus']))
                              if s['id'] != sid and state_for(c, 'spectrum', s['id']) == 'active' and root_id(c, s) == root]
            if family_current and (len(family_current) != 1 or old_id != family_current[0]):
                raise ValueError('该实验已有当前版本，请明确选择它作为替换目标，不能同时导出原谱和修订谱。')
            old = c.execute('SELECT * FROM spectra WHERE id=? AND molecule_id=?', (old_id, mid)).fetchone() if old_id else None
            if old_id and (not old or old_id == sid or old['nucleus'] != spec['nucleus']):
                raise ValueError('替换目标必须是本分子的另一条同核种谱。')
            q = c.execute('SELECT decision FROM spectrum_quality WHERE spectrum_id=?', (sid,)).fetchone()
            if q and q[0] != 'passed':
                raise ValueError('这个历史版本仍有自动隔离，不能采用。')
            if old:
                if state_for(c, 'spectrum', old_id) in ('reference', 'superseded'):
                    raise ValueError('替换目标已经不是当前版本，请重新选择。')
                put_state(c, mid, 'spectrum', old_id, 'superseded', reason)
            put_state(c, mid, 'spectrum', sid, 'active', reason)
            # Do not fabricate a new same-sample pair. Any replacement is independent.
            event.update(replaced_spectrum_id=old_id, adopted_spectrum_id=sid,
                         warning='人工采用实验归属表，不保证完整1D谱；不自动建立同条件H/C配对。')
        else:
            raise ValueError('未知数据管理操作。')
        c.execute('INSERT INTO curation_events(molecule_id,action,details_json,created_at) VALUES(?,?,?,?)', (mid, name, db.dumps(event), db.now()))
        c.execute('UPDATE molecules SET revision=revision+1,task_done=0 WHERE id=?', (mid,))
        db.refresh_molecule(c, mid)
        db.touch_activity(c, mid, 'curation:' + name)
        result = {'record': db.record(c, mid, include_archived=True), 'management': details(c, mid), 'event': event}
        c.commit()
        return result
    except Exception:
        c.rollback()
        raise
