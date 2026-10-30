"""Promote pre-versioning peak edits without losing their peaks or labels."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from annotation_tool_v3 import core as db, curation
from annotation_tool_v3.spectrum_versions import EDIT_ACTIONS, root_id


def candidates(c):
    ids = set()
    placeholders = ','.join('?' for _ in EDIT_ACTIONS)
    for row in c.execute('SELECT details_json FROM curation_events WHERE action IN ('+placeholders+')', sorted(EDIT_ACTIONS)):
        sid = json.loads(row[0]).get('spectrum_id')
        if sid:
            ids.add(sid)
    result = []
    for sid in sorted(ids):
        row = c.execute('SELECT * FROM spectra WHERE id=?', (sid,)).fetchone()
        if not row:
            continue
        provenance = json.loads(row['source_json'])
        revisions = provenance.get('spectrum_revisions') or {}
        if provenance.get('dataset_version') or provenance.get('reference_only_after_revision'):
            continue
        if not (any(x.get('active', True) for x in provenance.get('peak_table_supplements', []))
                or revisions.get('signal_codes') or revisions.get('assignment_groups')):
            continue
        result.append({'spectrum_id': sid, 'molecule_id': row['molecule_id'], 'root_spectrum_id': root_id(c, row)})
    return result


def snapshot(c, mid):
    result = {}
    for table in ('molecules', 'spectra', 'curation_state', 'curation_events', 'curation_issues', 'pairs'):
        column = 'id' if table == 'molecules' else 'molecule_id'
        result[table] = [dict(r) for r in c.execute(f'SELECT * FROM {table} WHERE {column}=?', (mid,))]
    for table in ('annotations', 'spectrum_quality', 'spectrum_completeness', 'spectrum_assignment_audit'):
        result[table] = [dict(r) for r in c.execute(f'SELECT * FROM {table} WHERE spectrum_id IN (SELECT id FROM spectra WHERE molecule_id=?)', (mid,))]
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', type=Path, default=db.DEFAULT_DB)
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--output-dir', type=Path, default=db.ROOT/'outputs')
    args = parser.parse_args()
    c = db.connect(args.db)
    try:
        selected = candidates(c)
        if not args.apply:
            print(json.dumps({'mode': 'dry_run', 'candidates': selected}, ensure_ascii=False, indent=2))
            return
        stamp = db.now().replace(':', '').replace('+', '_')
        args.output_dir.mkdir(parents=True, exist_ok=True)
        archive = args.output_dir / ('spectrum_versions_before_'+stamp+'.json')
        archive.write_text(json.dumps({mid: snapshot(c, mid) for mid in sorted({x['molecule_id'] for x in selected})}, ensure_ascii=False, indent=2), encoding='utf-8')
        results = []
        for item in selected:
            if sum(x['root_spectrum_id'] == item['root_spectrum_id'] for x in selected) != 1:
                results.append({**item, 'status': 'skipped', 'reason': 'multiple_edited_drafts_for_same_experiment'})
                continue
            mid = item['molecule_id']
            revision = c.execute('SELECT revision FROM molecules WHERE id=?', (mid,)).fetchone()[0]
            try:
                output = curation.action(c, mid, {'revision': revision, 'action': 'promote_existing_revision',
                    'spectrum_id': item['spectrum_id'], 'reason': '将旧补峰 / 修订转换为可标注的当前数据集版本；原谱保留为参考。'})
                results.append({**item, 'status': 'converted', **output['event']})
            except (ValueError, RuntimeError) as exc:
                results.append({**item, 'status': 'skipped', 'reason': str(exc)})
        report = {'created_at': db.now(), 'backup': str(archive), 'results': results}
        path = args.output_dir / ('spectrum_versions_upgrade_'+stamp+'.json')
        path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
        print(json.dumps({'converted': sum(x['status'] == 'converted' for x in results),
                          'skipped': sum(x['status'] == 'skipped' for x in results), 'backup': str(archive), 'report': str(path)}, ensure_ascii=False, indent=2))
    finally:
        c.close()


if __name__ == '__main__':
    main()
