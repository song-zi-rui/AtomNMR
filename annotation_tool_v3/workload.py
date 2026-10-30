"""Small append-only activity log; independent of the large spectral tables."""
from datetime import date, datetime, time, timedelta, timezone
import json

LOCAL_TIME = timezone(timedelta(hours=8), name='Asia/Shanghai')


def initialize(c):
    c.executescript('''
      CREATE TABLE IF NOT EXISTS work_events(
        event_key TEXT PRIMARY KEY,molecule_id TEXT NOT NULL REFERENCES molecules(id),
        action TEXT NOT NULL,created_at TEXT NOT NULL,source TEXT NOT NULL,
        details_json TEXT NOT NULL);
      CREATE INDEX IF NOT EXISTS work_events_time ON work_events(created_at,molecule_id);
    ''')
    if c.execute("SELECT 1 FROM meta WHERE key='workload_started_at'").fetchone():
        return
    # Recover timestamped evidence only. Original/imported annotations do not
    # become manual workload, and task_done alone cannot supply a finish date.
    for row in c.execute('SELECT id,molecule_id,action,created_at FROM curation_events').fetchall():
        append(c, row['molecule_id'], 'curation:' + row['action'], row['created_at'],
               key='curation:' + str(row['id']), source='recovered')
    for row in c.execute('SELECT id,molecule_id,decision,created_at FROM spectrum_numbering_reviews').fetchall():
        append(c, row['molecule_id'], 'numbering_' + row['decision'], row['created_at'],
               key='numbering:' + str(row['id']), source='recovered')
    for row in c.execute('SELECT * FROM molecule_activity').fetchall():
        stamp = normalize_stamp(row['last_processed_at'])
        # The last-activity pointer can refer to an already recovered audit row.
        if c.execute('''SELECT 1 FROM work_events WHERE molecule_id=? AND action=?
                       AND substr(created_at,1,19)=? LIMIT 1''',
                     (row['molecule_id'], row['last_action'], stamp[:19])).fetchone():
            continue
        append(c, row['molecule_id'], row['last_action'], stamp,
               key='activity:' + row['molecule_id'], source='recovered')
    started = datetime.now(timezone.utc).isoformat(timespec='microseconds')
    c.execute('INSERT INTO meta VALUES(?,?)', ('workload_started_at', json.dumps(started)))


def normalize_stamp(value):
    stamp = datetime.fromisoformat(value)
    if stamp.tzinfo is None:
        raise ValueError('工作记录时间缺少时区。')
    return stamp.astimezone(timezone.utc).isoformat(timespec='microseconds')


def append(c, mid, action, stamp, *, key=None, source='recorded', details=None):
    if key is None:
        revision = c.execute('SELECT revision FROM molecules WHERE id=?', (mid,)).fetchone()[0]
        key = f'operation:{mid}:{revision}:{action}'
    c.execute('INSERT OR IGNORE INTO work_events VALUES(?,?,?,?,?,?)',
              (key, mid, action, normalize_stamp(stamp), source,
               json.dumps(details or {}, separators=(',', ':'))))


def summarize(c, selected_date=None, days=14, *, current_time=None):
    current = current_time or datetime.now(timezone.utc)
    today = current.astimezone(LOCAL_TIME).date()
    selected = date.fromisoformat(selected_date) if selected_date else today
    days = max(1, min(90, int(days)))
    first = selected - timedelta(days=days - 1)
    start = datetime.combine(first, time.min, LOCAL_TIME).astimezone(timezone.utc)
    end = datetime.combine(selected + timedelta(days=1), time.min, LOCAL_TIME).astimezone(timezone.utc)
    keys = ('processed_molecules', 'finished_molecules', 'c_peaks', 'h_peaks',
            'save_operations', 'management_operations', 'recovered_records',
            'ai_molecules','ai_c_peaks','ai_h_peaks','ai_operations')
    buckets = {str(first + timedelta(days=i)): {k: set() if k in (*keys[:4],'ai_molecules','ai_c_peaks','ai_h_peaks') else 0 for k in keys}
               for i in range(days)}
    for row in c.execute('SELECT * FROM work_events WHERE created_at>=? AND created_at<? ORDER BY created_at',
                         (start.isoformat(timespec='microseconds'), end.isoformat(timespec='microseconds'))):
        day = str(datetime.fromisoformat(row['created_at']).astimezone(LOCAL_TIME).date())
        bucket = buckets[day]
        if row['source']=='ai':
            bucket['ai_molecules'].add(row['molecule_id']);bucket['ai_operations']+=1
            for peak in json.loads(row['details_json']).get('edited_peaks',[]):
                bucket['ai_c_peaks' if peak['nucleus']=='13C' else 'ai_h_peaks'].add((peak['spectrum_id'],peak['peak_id']))
            continue
        bucket['processed_molecules'].add(row['molecule_id'])
        if row['action'] == 'finish':
            bucket['finished_molecules'].add(row['molecule_id'])
        if row['action'] in ('save', 'finish'):
            bucket['save_operations'] += 1
        if row['action'].startswith(('curation:', 'numbering_')):
            bucket['management_operations'] += 1
        if row['source'] == 'recovered':
            bucket['recovered_records'] += 1
        for peak in json.loads(row['details_json']).get('edited_peaks', []):
            target = 'c_peaks' if peak['nucleus'] == '13C' else 'h_peaks'
            bucket[target].add((peak['spectrum_id'], peak['peak_id']))
    history = [{'date': day, **{k: len(v) if isinstance(v, set) else v for k, v in bucket.items()}}
               for day, bucket in buckets.items()]
    started = c.execute("SELECT value FROM meta WHERE key='workload_started_at'").fetchone()
    return {'timezone': 'Asia/Shanghai', 'today': str(today), 'date': str(selected),
            'day': history[-1], 'history': list(reversed(history)),
            'recording_started_at': json.loads(started[0]) if started else None,
            'generated_at': current.isoformat(timespec='seconds'),
            'definitions': {
                'processed_molecules': '当天成功保存或管理过的分子数；按分子去重。',
                'finished_molecules': '当天提交结束的分子数；重复结束计一个，之后重开仍保留当天的工作记录。',
                'c_peaks': '当天通过归属保存提交的不同碳峰数；同一谱同一峰去重，清空和峰备注修改也计入；谱表修订记入管理操作。',
                'h_peaks': '当天通过归属保存提交的不同氢峰数；同一谱同一峰去重，清空和峰备注修改也计入；谱表修订记入管理操作。',
                'save_operations': '已记录的成功保存次数，包含自动保存和结束提交。',
                'management_operations': '已记录的数据管理与编号审核次数。',
            }}
