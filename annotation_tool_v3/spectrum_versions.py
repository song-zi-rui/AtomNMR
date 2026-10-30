"""One immutable original and one reusable final spectrum per experiment family."""
from __future__ import annotations

import copy
import json
import uuid

from . import core as db

EDIT_ACTIONS = {'supplement_peak', 'withdraw_peak_supplement', 'edit_peak', 'delete_peak',
                'signal_code', 'assignment_group', 'remove_assignment_group'}
POLICY = 'original-and-final-spectrum-v1'


def root_id(c, spec):
    seen = set()
    while spec and spec['id'] not in seen:
        seen.add(spec['id'])
        provenance = json.loads(spec['source_json'])
        version = provenance.get('dataset_version') or {}
        if version.get('root_spectrum_id'):
            return version['root_spectrum_id']
        parent = provenance.get('parent_spectrum_id')
        if not parent:
            return spec['id']
        row = c.execute('SELECT * FROM spectra WHERE id=? AND molecule_id=? AND nucleus=?',
                        (parent, spec['molecule_id'], spec['nucleus'])).fetchone()
        if row is None:
            return spec['id']
        spec = row
    raise ValueError('谱图版本关系出现循环，请检查来源。')


def family_rows(c, spec):
    root = root_id(c, spec)
    return root, [s for s in c.execute('SELECT * FROM spectra WHERE molecule_id=? AND nucleus=?',
                                     (spec['molecule_id'], spec['nucleus'])) if root_id(c, s) == root]


def final_for(c, spec):
    from .curation import state_for
    root, family = family_rows(c, spec)
    choices = [s for s in family if s['id'] != root and s['source'] == 'manual_revision']
    active = [s for s in choices if state_for(c, 'spectrum', s['id']) == 'active']
    if len(active) > 1:
        raise ValueError('同一实验存在多个当前修订，请先整理版本。')
    return active[0] if active else max(choices, key=lambda s: (
        (json.loads(s['source_json']).get('dataset_version') or {}).get('number', 0), s['id']), default=None)


def create_current(c, mid, spec, reason, action):
    from .curation import state_for, put_state
    from .spectrum_revision import require_editable
    require_editable(c, mid, spec)
    provenance = json.loads(spec['source_json'])
    quality = c.execute('SELECT * FROM spectrum_quality WHERE spectrum_id=?', (spec['id'],)).fetchone()
    inherited_quality = json.loads(quality['details_json']) if quality else provenance.get('parent_automatic_quality')
    if inherited_quality and inherited_quality.get('decision') != 'passed':
        raise ValueError('来源仍未通过检查或已隔离；修订不会自动清除来源问题。')
    if spec['source'] == 'NMRexp' and not quality:
        raise ValueError('这条 NMRexp 来源尚未检查，不能直接生成当前数据集版本。')
    root, family = family_rows(c, spec)
    active = [s for s in family if state_for(c, 'spectrum', s['id']) == 'active']
    parent = provenance.get('parent_spectrum_id')
    if state_for(c, 'spectrum', spec['id']) == 'reference' and any(s['id'] not in (parent, root) for s in active):
        raise ValueError('这份副本已有更新的当前版本，请改为编辑当前版本。')
    if not active:
        raise ValueError('该实验家族没有当前使用版本，请先明确恢复或采用对象。')
    # Once the final slot is active, subsequent edits update it in place.
    if spec['id'] != root and spec['source'] == 'manual_revision' and state_for(c, 'spectrum', spec['id']) == 'active':
        if len(active) != 1:
            raise ValueError('同一实验的当前版本不唯一，请先整理版本。')
        version = provenance.get('dataset_version') or {}
        version.update(policy=POLICY, root_spectrum_id=root, parent_spectrum_id=root, number=1,
                       edit_count=version.get('edit_count', version.get('number', 0)) + 1,
                       last_edit_at=db.now(), reason=reason, action=action, role='final')
        provenance.update(dataset_version=version, parent_spectrum_id=root)
        c.execute('UPDATE spectra SET source_json=? WHERE id=?', (db.dumps(provenance), spec['id']))
        return c.execute('SELECT * FROM spectra WHERE id=?', (spec['id'],)).fetchone(), {
            'input_spectrum_id': spec['id'], 'current_spectrum_id': spec['id'],
            'new_spectrum_ids': [], 'reference_spectrum_ids': [root], 'version_number': 1,
            'updated_final_in_place': True, 'inherited_pair_ids': []}
    existing_final = final_for(c, spec)
    # db.add_spectrum generates exactly this ID from source and key.
    key = existing_final['source_key'] if existing_final else uuid.uuid4().hex
    new_id = existing_final['id'] if existing_final else db.uid('S', 'manual_revision', key)
    # Capture the base before reusing a reference draft or replacing an inactive final.
    visible = next(s for s in db.record(c, mid, include_archived=True)['spectra'] if s['id'] == spec['id'])
    old_annotations = c.execute('SELECT * FROM annotations WHERE spectrum_id=?', (spec['id'],)).fetchall()
    audit = c.execute('SELECT * FROM spectrum_assignment_audit WHERE spectrum_id=?', (spec['id'],)).fetchone()
    complete = c.execute('SELECT status,reason FROM spectrum_completeness WHERE spectrum_id=?', (spec['id'],)).fetchone()
    provenance = copy.deepcopy(provenance)
    provenance.pop('reference_only_after_revision', None)
    original_source = c.execute('SELECT source FROM spectra WHERE id=?', (root,)).fetchone()[0]
    provenance.update(parent_spectrum_id=root, original_source=provenance.get('original_source', original_source),
                      requires_explicit_adoption=False)
    provenance['dataset_version'] = {'policy': POLICY, 'root_spectrum_id': root,
        'parent_spectrum_id': root, 'number': 1, 'edit_count': 1, 'reason': reason, 'action': action,
        'created_at': db.now(), 'role': 'final'}
    # Group undo snapshots belong to this copied spectrum, with stable peak IDs.
    for group in (provenance.get('spectrum_revisions') or {}).get('assignment_groups', []):
        for item in group.get('before_assignments', []):
            if item.get('annotation'):
                item['annotation']['spectrum_id'] = new_id
        for item in group.get('applied_annotations', []):
            item['spectrum_id'] = new_id
    if existing_final:
        c.execute('UPDATE spectra SET solvent=?,frequency_mhz=?,temperature_k=?,conditions_json=?,peaks_json=?,source_json=? WHERE id=?',
                  (spec['solvent'], spec['frequency_mhz'], spec['temperature_k'], spec['conditions_json'],
                   spec['peaks_json'], db.dumps(provenance), new_id))
        for table in ('annotations', 'spectrum_assignment_audit', 'spectrum_quality', 'spectrum_completeness'):
            c.execute(f'DELETE FROM {table} WHERE spectrum_id=?', (new_id,))
    else:
        db.add_spectrum(c, mid, spec['nucleus'], 'manual_revision', key,
                        json.loads(spec['conditions_json']), json.loads(spec['peaks_json']), provenance)
    visible_labels = {p['peak_id']: p['assignment'] for p in visible['peaks']}
    for old in old_annotations:
        value = dict(old)
        projected = visible_labels.get(old['peak_id'])
        if projected:
            value.update(status=projected['status'], atoms_json=db.dumps(projected['atom_indices']))
        c.execute('INSERT INTO annotations VALUES(?,?,?,?,?,?,?)',
                  (new_id, old['peak_id'], value['status'], value['atoms_json'], old['note'], old['origin'], old['updated_at']))
    if audit:
        c.execute('INSERT OR REPLACE INTO spectrum_assignment_audit VALUES(?,?,?,?,?)',
                  (new_id, audit['policy'], audit['validation_json'], audit['identity_json'], audit['updated_at']))
    if inherited_quality:
        c.execute('INSERT INTO spectrum_quality VALUES(?,?,?,?,?)',
                  (new_id, 'manual-version-inherited-source-gate-v1', 'passed', db.dumps({
                      **inherited_quality, 'content_revision': 'user_experimental_revision',
                      'source_gate_inherited_from': root, 'revision_content_automatically_validated': False}), db.now()))
    if complete:
        c.execute('INSERT INTO spectrum_completeness VALUES(?,?,?,?)', (new_id, complete[0], complete[1], db.now()))
    retired = {root, spec['id'], *(s['id'] for s in active)} - {new_id}
    for old_id in retired:
        put_state(c, mid, 'spectrum', old_id, 'reference', '修订版本已作为当前数据集；原版本仅供参考。' + reason)
        old_provenance = json.loads(c.execute('SELECT source_json FROM spectra WHERE id=?', (old_id,)).fetchone()[0])
        old_provenance['reference_only_after_revision'] = new_id
        c.execute('UPDATE spectra SET source_json=? WHERE id=?', (db.dumps(old_provenance), old_id))
    put_state(c, mid, 'spectrum', new_id, 'active', reason)
    # Keep the existing experiment pairing and its limitations, never infer a new sample match.
    inherited_pairs = []
    for old_id in retired:
        for pair in c.execute('SELECT * FROM pairs WHERE molecule_id=? AND (c_id=? OR h_id=?)', (mid, old_id, old_id)).fetchall():
            other = pair['h_id'] if spec['nucleus'] == '13C' else pair['c_id']
            eligible = c.execute('SELECT 1 FROM spectra WHERE id=? AND' + db.ACTIVE, (other,)).fetchone()
            if not eligible:
                continue
            cs, hs = (new_id, other) if spec['nucleus'] == '13C' else (other, new_id)
            db.add_pair(c, mid, cs, hs, pair['basis'] + '|spectrum_revision_inherits_pair')
            inherited_pairs.append(pair['id'])
    return c.execute('SELECT * FROM spectra WHERE id=?', (new_id,)).fetchone(), {
        'input_spectrum_id': spec['id'], 'current_spectrum_id': new_id,
        'new_spectrum_ids': [] if existing_final else [new_id], 'reference_spectrum_ids': sorted(retired),
        'version_number': 1, 'updated_final_in_place': bool(existing_final), 'inherited_pair_ids': inherited_pairs}


def annotate_role(spectrum):
    version = spectrum['provenance'].get('dataset_version')
    active = spectrum['manual_state'] == 'active' and spectrum.get('quality', spectrum.get('automatic_quality', {})).get('decision') != 'quarantined'
    spectrum['dataset_role'] = 'current_dataset' if active else 'reference_only' if spectrum['manual_state'] in ('reference', 'superseded') else 'inactive'
    spectrum['version_info'] = version
    spectrum['version_label'] = '最终版' if version else '原始版' if spectrum['provenance'].get('reference_only_after_revision') else None
    spectrum['included_in_dataset_export'] = active
