"""Cached ¹H peak information used to order unassigned molecules."""
from __future__ import annotations

from contextlib import closing
import argparse
import json
import math
import re
import sqlite3
import time

STRIDE = 4000
BASE = 10000
UNASSIGNED_H = "(queue='both_unassigned' OR (queue='waiting_spectrum' AND h_count>0 AND c_assigned+h_assigned=0))"
MAX_SCORE = STRIDE - 1
PEAK_TYPE = re.compile(r'(?:br(?:oad)?[. -]*)?(?:[sdtqphm]+|br|singlet|doublet|triplet|quartet|multiplet)')


def positive(value):
    return type(value) in (int, float) and math.isfinite(value) and value > 0


def spectrum_score(peaks, excluded=()):
    """Prioritize peak types plus confirmed Hz values, then field coverage."""
    count = kinds = couplings = integrals = 0
    for peak in peaks:
        if peak['peak_id'] in excluded:
            continue
        count += 1
        kind = str(peak.get('multiplicity') or '').strip().lower()
        kinds += PEAK_TYPE.fullmatch(kind) is not None
        values = peak.get('j_values')
        unit = str(peak.get('j_unit') or '').strip().lower()
        # Import places uncertain values in j_values_unconfirmed instead.
        # Manually supplemented j_values are specified in Hz by the editor.
        couplings += unit in ('', 'hz') and isinstance(values, list) and any(map(positive, values))
        integrals += positive(peak.get('reported_nH')) or positive(peak.get('integral'))
    if not count:
        return 0
    tier = 3 if kinds and couplings else 2 if kinds or couplings else 1 if integrals else 0
    return 1000 * tier + (kinds + 2 * couplings + integrals) * 999 // (4 * count)


def backfill(path, batch_size=512, progress=None):
    """Resume the one-time backfill in small transactions on the main database."""
    from . import core as db
    from .peak_table import KEY, effective_peaks

    if not db.main_database(path).is_file():
        raise FileNotFoundError(f'Main database does not exist: {path}')
    with closing(db.connect(path, timeout=2)) as conn:
        conn.execute('PRAGMA cache_size=-65536')
        conn.execute('PRAGMA synchronous=NORMAL')
        conn.execute('PRAGMA wal_autocheckpoint=8192')
        with conn:
            conn.execute('BEGIN IMMEDIATE')
            db.initialize_proton_order(conn)
        for queue in ('both_unassigned', 'waiting_spectrum'):
            stage = 'h_info_priority_v1:' + queue
            processed = changed = 0
            while True:
                try:
                    conn.execute('BEGIN IMMEDIATE')
                    state = conn.execute('SELECT cursor FROM progress WHERE stage=?', (stage,)).fetchone()
                    cursor = state[0] if state else ''
                    if cursor == 'complete':
                        conn.rollback()
                        break
                    ids = [row[0] for row in conn.execute(f'''SELECT id FROM molecules INDEXED BY molecule_queue
                        WHERE queue=? AND task_done=0 AND id>? AND {UNASSIGNED_H} ORDER BY id LIMIT ?''',
                        (queue, cursor, batch_size))]
                    scores = dict.fromkeys(ids, 0)
                    if ids:
                        params = ','.join('?' for _ in ids)
                        rows = conn.execute(f'''SELECT id,molecule_id,peaks_json,source_json FROM spectra
                            WHERE molecule_id IN ({params}) AND nucleus='1H' AND {db.ACTIVE}''', ids).fetchall()
                        excluded = {}
                        # Unassigned entries can still contain manually marked solvent/impurity peaks.
                        for begin in range(0, len(rows), 256):
                            spectra = [row['id'] for row in rows[begin:begin + 256]]
                            spec_params = ','.join('?' for _ in spectra)
                            for spec, peak in conn.execute(f'''SELECT spectrum_id,peak_id FROM annotations
                                WHERE spectrum_id IN ({spec_params}) AND status IN ('solvent','impurity','excluded')''', spectra):
                                excluded.setdefault(spec, set()).add(peak)
                        for row in rows:
                            peaks = json.loads(row['peaks_json'])
                            source = json.loads(row['source_json'])
                            if source.get(KEY):
                                peaks = effective_peaks(peaks, source)
                            mid = row['molecule_id']
                            scores[mid] = max(scores[mid], spectrum_score(peaks, excluded.get(row['id'], ())))
                        before = conn.total_changes
                        conn.executemany('UPDATE molecules SET h_info_priority=? WHERE id=? AND h_info_priority!=?',
                                         ((score, mid, score) for mid, score in scores.items()))
                        changed += conn.total_changes - before
                    conn.execute('''INSERT INTO progress(stage,cursor,updated_at) VALUES(?,?,?)
                        ON CONFLICT(stage) DO UPDATE SET cursor=excluded.cursor,updated_at=excluded.updated_at''',
                        (stage, ids[-1] if ids else 'complete', db.now()))
                    conn.commit()
                except sqlite3.OperationalError as exc:
                    conn.rollback()
                    if not db.is_busy_error(exc):
                        raise
                    time.sleep(.2)
                    continue
                if not ids:
                    if progress:
                        progress({'queue': queue, 'processed_this_run': processed, 'changed_this_run': changed, 'complete': True})
                    break
                processed += len(ids)
                if progress and processed % (batch_size * 100) == 0:
                    progress({'queue': queue, 'processed_this_run': processed, 'changed_this_run': changed})


def main():
    from pathlib import Path
    from .paths import DEFAULT_DB
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', type=Path, default=DEFAULT_DB)
    args = parser.parse_args()
    backfill(args.db, progress=lambda value: print(json.dumps(value), flush=True))


if __name__ == '__main__':
    main()
