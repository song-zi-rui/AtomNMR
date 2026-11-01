"""User-defined completion policy; keep coverage and saved annotations separate."""
import json

POLICY = 'finished_spectra_complete_20261031_v1'


def ensure_history(c):
    c.execute('''CREATE TABLE IF NOT EXISTS workflow_completeness_history(
      id INTEGER PRIMARY KEY,molecule_id TEXT NOT NULL,spectrum_id TEXT NOT NULL,
      policy TEXT NOT NULL,previous_json TEXT,new_reason TEXT NOT NULL,created_at TEXT NOT NULL)''')
    c.execute('CREATE INDEX IF NOT EXISTS workflow_completeness_spectrum ON workflow_completeness_history(spectrum_id,policy)')


def mark_complete(c, mid, spectrum_ids=None, *, reason, skip_previous_policy=False):
    """Only current usable spectra; no invented assignments or global work."""
    from . import core as db
    ensure_history(c)
    wanted=None if spectrum_ids is None else set(spectrum_ids)
    rows=c.execute('SELECT id FROM spectra WHERE molecule_id=? AND'+db.ACTIVE,(mid,)).fetchall()
    changed=0
    for row in rows:
        sid=row['id']
        if wanted is not None and sid not in wanted:continue
        if skip_previous_policy and c.execute('SELECT 1 FROM workflow_completeness_history WHERE spectrum_id=? AND policy=? LIMIT 1',(sid,POLICY)).fetchone():continue
        previous=c.execute('SELECT status,reason,updated_at FROM spectrum_completeness WHERE spectrum_id=?',(sid,)).fetchone()
        if previous and previous['status']=='complete':continue
        stamp=db.now()
        c.execute('INSERT INTO workflow_completeness_history(molecule_id,spectrum_id,policy,previous_json,new_reason,created_at) VALUES(?,?,?,?,?,?)',
                  (mid,sid,POLICY,db.dumps(dict(previous)) if previous else None,reason,stamp))
        c.execute('''INSERT INTO spectrum_completeness VALUES(?,?,?,?) ON CONFLICT(spectrum_id)
          DO UPDATE SET status=excluded.status,reason=excluded.reason,updated_at=excluded.updated_at''',(sid,'complete',reason,stamp))
        changed+=1
    return changed


def backfill(c, *, catch_up=False):
    """Run once for finished and historically human/AI-annotated molecules."""
    from . import core as db
    old=c.execute('SELECT value FROM meta WHERE key=?',(POLICY,)).fetchone()
    if old and not catch_up:return {'already_applied':True,**json.loads(old[0])}
    ensure_history(c)
    ids={r[0] for r in c.execute('SELECT id FROM molecules WHERE task_done=1')}
    ids.update(r[0] for r in c.execute('''SELECT DISTINCT s.molecule_id FROM annotations a JOIN spectra s ON s.id=a.spectrum_id
       WHERE a.origin IN ('manual','manual_v2') OR (a.origin>='ai:' AND a.origin<'ai;')'''))
    ids.update(r[0] for r in c.execute('SELECT DISTINCT molecule_id FROM ai_annotation_results'))
    spectra=changed_molecules=0
    for mid in sorted(ids):
        changed=mark_complete(c,mid,reason='按用户标注流程约定回填：此前已结束或已有人工/AI标注，当前实验谱设为完整。',skip_previous_policy=bool(old))
        if changed:
            c.execute('UPDATE molecules SET revision=revision+1 WHERE id=?',(mid,))
            spectra+=changed;changed_molecules+=1
    result={'eligible_molecules':len(ids),'changed_molecules':changed_molecules,'changed_spectra':spectra,'applied_at':db.now()}
    db.set_meta(c,POLICY,result)
    return {'already_applied':bool(old) and not spectra,**result}


def migrate(c):
    from . import core as db
    # Existing names are explicitly used in navigation queries: replace their
    # old expression together, keeping the indexes useful for the new order.
    for name,columns in (
        ('molecule_priority',f'({db.QUEUE_ORDER}),id'),
        ('molecule_pending_priority',f'task_done,({db.QUEUE_ORDER}),id'),
        ('molecule_coverage_processing',f'({db.COVERAGE_SQL}),task_done,({db.QUEUE_ORDER}),id')):
        row=c.execute("SELECT sql FROM sqlite_master WHERE type='index' AND name=?",(name,)).fetchone()
        if row and db.QUEUE_ORDER not in row[0]:c.execute('DROP INDEX '+name)
        c.execute('CREATE INDEX IF NOT EXISTS '+name+' ON molecules('+columns+')')
    return backfill(c,catch_up=True)
