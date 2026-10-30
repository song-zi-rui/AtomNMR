"""Move legacy intermediate spectra into a recovery archive, retaining two slots."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sqlite3
import sys
from contextlib import closing

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from annotation_tool_v3 import core as db
from annotation_tool_v3.spectrum_versions import root_id, family_rows, final_for, POLICY


def plan(c):
    # Curated spectra have state rows; avoid scanning the million raw spectra.
    candidates = c.execute("SELECT s.* FROM curation_state cs JOIN spectra s ON s.id=cs.target_id WHERE cs.scope='spectrum' AND s.source='manual_revision'").fetchall()
    result, seen = [], set()
    for spec in candidates:
        root = root_id(c, spec)
        if root in seen:
            continue
        seen.add(root)
        _, family = family_rows(c, spec)
        final = final_for(c, spec)
        removed = [s for s in family if s['id'] not in (root, final['id'])] if final else []
        if not removed:
            continue
        if any(s['source'] != 'manual_revision' for s in removed):
            raise ValueError('实验家族含独立来源，不能自动合并。')
        result.append({'molecule_id': spec['molecule_id'], 'root_id': root, 'final_id': final['id'],
                       'removed_ids': sorted(s['id'] for s in removed)})
    return result


def consolidate(c, archive_path):
    archive_path = Path(archive_path).resolve()
    if archive_path.exists():
        raise FileExistsError('恢复备份已存在，不能覆盖。')
    c.execute('BEGIN IMMEDIATE')
    try:
        groups = plan(c)
        if not groups:
            c.rollback()
            return {'families': [], 'removed_spectra': 0, 'archive': None}
        known = {'annotations', 'pairs', 'spectrum_quality', 'spectrum_assignment_audit',
                 'spectrum_completeness', 'spectrum_numbering_reviews', 'curation_issues'}
        for table in c.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall():
            name = table[0]
            for fk in c.execute(f'PRAGMA foreign_key_list("{name}")'):
                if fk['table'] == 'spectra' and name not in known:
                    raise ValueError('发现未支持的谱引用表，停止整理：' + name)
        ids = sorted({sid for g in groups for sid in (g['root_id'], g['final_id'], *g['removed_ids'])})
        mids = sorted({g['molecule_id'] for g in groups})
        placeholders = ','.join('?' for _ in ids)
        snapshots = []
        def snapshot(table, where, values):
            snapshots.extend((table, db.dumps(dict(row))) for row in c.execute(f'SELECT * FROM {table} WHERE {where}', values))
        snapshot('spectra', f'id IN ({placeholders})', ids)
        for table in ('annotations', 'spectrum_quality', 'spectrum_assignment_audit', 'spectrum_completeness', 'spectrum_numbering_reviews', 'curation_issues'):
            snapshot(table, f'spectrum_id IN ({placeholders})', ids)
        snapshot('pairs', f'c_id IN ({placeholders}) OR h_id IN ({placeholders})', ids + ids)
        snapshot('curation_state', f"scope='spectrum' AND target_id IN ({placeholders})", ids)
        mp = ','.join('?' for _ in mids)
        for table in ('molecules', 'molecule_activity', 'curation_events'):
            snapshot(table, ('id' if table == 'molecules' else 'molecule_id') + f' IN ({mp})', mids)
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        with closing(sqlite3.connect(archive_path)) as archive:
            with archive:
                archive.execute('CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL)')
                archive.execute('CREATE TABLE snapshots(table_name TEXT NOT NULL,row_json TEXT NOT NULL)')
                archive.executemany('INSERT INTO snapshots VALUES(?,?)', snapshots)
                archive.executemany('INSERT INTO metadata VALUES(?,?)', [('plan', db.dumps(groups)),
                    ('created_at', db.now()), ('status', 'prepared'), ('policy', POLICY)])
        # No destructive change occurs until the recoverable snapshots are committed.
        for g in groups:
            final = c.execute('SELECT source_json FROM spectra WHERE id=?', (g['final_id'],)).fetchone()
            provenance = json.loads(final[0])
            version = provenance.get('dataset_version') or {}
            version.update(policy=POLICY, root_spectrum_id=g['root_id'], parent_spectrum_id=g['root_id'],
                           edit_count=version.get('edit_count', version.get('number', 1)), number=1, role='final')
            provenance.update(dataset_version=version, parent_spectrum_id=g['root_id'])
            c.execute('UPDATE spectra SET source_json=? WHERE id=?', (db.dumps(provenance), g['final_id']))
            original = json.loads(c.execute('SELECT source_json FROM spectra WHERE id=?', (g['root_id'],)).fetchone()[0])
            if original.get('reference_only_after_revision') in g['removed_ids']:
                original['reference_only_after_revision'] = g['final_id']
                c.execute('UPDATE spectra SET source_json=? WHERE id=?', (db.dumps(original), g['root_id']))
            for sid in g['removed_ids']:
                for table in ('curation_issues', 'spectrum_numbering_reviews'):
                    c.execute(f'UPDATE {table} SET spectrum_id=? WHERE spectrum_id=?', (g['final_id'], sid))
                c.execute('DELETE FROM pairs WHERE c_id=? OR h_id=?', (sid, sid))
                for table in ('annotations', 'spectrum_quality', 'spectrum_assignment_audit', 'spectrum_completeness'):
                    c.execute(f'DELETE FROM {table} WHERE spectrum_id=?', (sid,))
                c.execute("DELETE FROM curation_state WHERE scope='spectrum' AND target_id=?", (sid,))
                c.execute('DELETE FROM spectra WHERE id=?', (sid,))
        for mid in mids:
            db.refresh_molecule(c, mid)
            c.execute('UPDATE molecules SET revision=revision+1 WHERE id=?', (mid,))
            c.execute('INSERT INTO curation_events(molecule_id,action,details_json,created_at) VALUES(?,?,?,?)',
                      (mid, 'consolidate_versions', db.dumps({'families': [g for g in groups if g['molecule_id'] == mid],
                       'recovery_archive': str(archive_path)}), db.now()))
        c.commit()
        with closing(sqlite3.connect(archive_path)) as archive:
            with archive:
                archive.execute("UPDATE metadata SET value='committed' WHERE key='status'")
        return {'families': groups, 'removed_spectra': sum(len(g['removed_ids']) for g in groups), 'archive': str(archive_path)}
    except Exception:
        c.rollback()
        raise


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--db', type=Path, default=db.DEFAULT_DB)
    p.add_argument('--apply', action='store_true')
    p.add_argument('--archive', type=Path)
    args = p.parse_args()
    if args.apply and not args.archive:
        p.error('--apply requires --archive')
    with db.connect(args.db) as c:
        result = consolidate(c, args.archive) if args.apply else plan(c)
        print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
