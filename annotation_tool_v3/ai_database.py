"""Atomic, idempotent publication of machine annotations in the main database."""
import json
import sqlite3
import time
from contextlib import closing

from . import core as db
from .ai_consistency import review_issues


def initialize(c):
    c.executescript('''
    CREATE TABLE IF NOT EXISTS ai_annotation_results(
      id TEXT PRIMARY KEY,molecule_id TEXT NOT NULL REFERENCES molecules(id),job_id TEXT NOT NULL,
      provider TEXT NOT NULL,model TEXT NOT NULL,response_model TEXT,prompt_hash TEXT,
      input_revision INTEGER,input_hash TEXT,input_json TEXT NOT NULL,output_json TEXT NOT NULL,reference_json TEXT NOT NULL,
      created_at TEXT NOT NULL,applied_count INTEGER NOT NULL,skipped_count INTEGER NOT NULL);
    CREATE INDEX IF NOT EXISTS ai_result_model_molecule ON ai_annotation_results(model,molecule_id);
    CREATE INDEX IF NOT EXISTS ai_result_molecule ON ai_annotation_results(molecule_id,created_at);
    CREATE TABLE IF NOT EXISTS ai_peak_predictions(
      result_id TEXT NOT NULL REFERENCES ai_annotation_results(id),spectrum_id TEXT NOT NULL,
      peak_id TEXT NOT NULL,nucleus TEXT NOT NULL,status TEXT NOT NULL,atoms_json TEXT NOT NULL,
      candidates_json TEXT NOT NULL,confidence TEXT,reason TEXT,applied INTEGER NOT NULL,skip_reason TEXT,
      PRIMARY KEY(result_id,spectrum_id,peak_id));
    CREATE INDEX IF NOT EXISTS ai_prediction_peak ON ai_peak_predictions(spectrum_id,peak_id,applied);
    CREATE INDEX IF NOT EXISTS annotation_origin ON annotations(origin,spectrum_id);
    CREATE INDEX IF NOT EXISTS molecule_auto_queue ON molecules(queue,id);
    CREATE TABLE IF NOT EXISTS ai_annotation_replacements(
      result_id TEXT NOT NULL,spectrum_id TEXT NOT NULL,peak_id TEXT NOT NULL,
      before_json TEXT NOT NULL,after_json TEXT NOT NULL,created_at TEXT NOT NULL,
      PRIMARY KEY(result_id,spectrum_id,peak_id));
    ''')


def publish(path,result_id,job_id,mid,input_revision,context,output,references,config):
    from .ai_annotation import validate_output,HARD
    from .ai_protocol import requested
    import hashlib
    output=validate_output(output,context)|{k:output[k] for k in ('response_model','response_id','api_diagnostics','validation_issues','repair_targets','requested_peaks','supplied_anchor_peaks') if k in output}
    output['consistency_issues']=review_issues(output,context)
    review_only={(r['spectrum_id'],r['peak_id']) for r in output['consistency_issues']}
    for attempt in range(8):
        with closing(db.connect(path,timeout=.25)) as c:
            try:
                c.execute('BEGIN IMMEDIATE')
                previous=c.execute('SELECT applied_count,skipped_count FROM ai_annotation_results WHERE id=?',(result_id,)).fetchone()
                if previous:return {'result_id':result_id,'applied':previous[0],'skipped':previous[1],'already_published':True}
                current=db.record(c,mid,require_quality=True)
                if current['smiles']!=context['smiles']:raise ValueError('分子结构已经变化，不能应用旧结果')
                spectra={s['id']:s for s in current['spectra'] if s.get('editable',True)}
                inputs={s['spectrum_id']:s for s in context['spectra']}
                model=config['model'];provider=config.get('provider','deepseek');origin='ai:'+provider+':'+model
                c.execute('INSERT INTO ai_annotation_results VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                          (result_id,mid,job_id,provider,model,output.get('response_model'),config.get('prompt_sha256'),input_revision,
                           hashlib.sha256(db.dumps(context).encode()).hexdigest(),db.dumps(context),db.dumps(output),db.dumps(references),db.now(),0,0))
                applied,skipped,edited=0,0,[];wanted=requested(context)
                repair_refs={(r['spectrum_id'],r.get('peak_id')):r for r in references if r.get('repair')}
                for p in output['assignments']:
                    sid,pid=p['spectrum_id'],p['peak_id'];spec=spectra.get(sid)
                    peak=next((x for x in spec['peaks'] if x['peak_id']==pid),None) if spec else None
                    before=next((x for x in inputs[sid]['peaks'] if x['peak_id']==pid),None)
                    reason=None
                    if (sid,pid) not in wanted:reason='not_requested'
                    elif not current['training_eligible'] or not spec or not peak:reason='spectrum_not_current'
                    elif any(peak.get(k)!=v for k,v in before.items()):reason='measurement_changed'
                    elif any(spec.get('conditions',{}).get(k)!=v for k,v in inputs[sid]['conditions'].items()):reason='conditions_changed'
                    elif peak.get('assignment_group'):reason='manual_unordered_group'
                    old=c.execute('SELECT * FROM annotations WHERE spectrum_id=? AND peak_id=?',(sid,pid)).fetchone()
                    if not reason and config.get('repair_policy')=='shiftdb-equivalence-repair-v1' and not old:reason='source_label_changed'
                    if not reason and old:
                        from . import equivalence_repair
                        reason=(equivalence_repair.replacement_reason(peak,old,p,repair_refs.get((sid,pid)),spec,context)
                                if config.get('repair_policy')==equivalence_repair.POLICY else 'existing_label_preserved')
                    if not reason and (sid,pid) in review_only:reason='assignment_consistency_review'
                    if not reason and p['status'] not in HARD:reason='not_a_definite_assignment'
                    did_apply=not reason
                    if did_apply:
                        if old:
                            c.execute('INSERT INTO ai_annotation_replacements VALUES(?,?,?,?,?,?)',
                                      (result_id,sid,pid,db.dumps(dict(old)),db.dumps(p),db.now()))
                        c.execute('INSERT INTO annotations VALUES(?,?,?,?,?,?,?) ON CONFLICT(spectrum_id,peak_id) DO UPDATE SET status=excluded.status,atoms_json=excluded.atoms_json,note=excluded.note,origin=excluded.origin,updated_at=excluded.updated_at',
                                  (sid,pid,p['status'],db.dumps(p['atom_indices']),p['brief_reason'],origin,db.now()))
                        applied+=1;edited.append({'spectrum_id':sid,'peak_id':pid,'nucleus':inputs[sid]['nucleus']})
                    else:skipped+=1
                    c.execute('INSERT INTO ai_peak_predictions VALUES(?,?,?,?,?,?,?,?,?,?,?)',
                              (result_id,sid,pid,inputs[sid]['nucleus'],p['status'],db.dumps(p['atom_indices']),db.dumps(p['candidates']),
                               p['confidence'],p['brief_reason'],int(did_apply),reason))
                c.execute('UPDATE ai_annotation_results SET applied_count=?,skipped_count=? WHERE id=?',(applied,skipped,result_id))
                from .workflow_policy import mark_complete
                completed=mark_complete(c,mid,inputs,reason='按用户标注流程约定：AI标注结果已入库，已分析的当前实验谱设为完整。')
                if applied or completed:
                    c.execute('UPDATE molecules SET revision=revision+1 WHERE id=?',(mid,))
                if applied:
                    db.refresh_molecule(c,mid)
                from . import workload
                workload.append(c,mid,'ai_annotation',db.now(),key='ai:'+result_id,source='ai',
                                details={'model':model,'job_id':job_id,'result_id':result_id,'edited_peaks':edited})
                c.commit()
                return {'result_id':result_id,'applied':applied,'skipped':skipped,'already_published':False}
            except sqlite3.Error as exc:
                c.rollback()
                if not db.is_busy_error(exc) or attempt==7:raise
        time.sleep(min(.05*2**attempt,.5))


def peak_provenance(c,spectrum_id):
    rows=c.execute('''SELECT p.peak_id,r.id result_id,r.job_id,r.provider,r.model,r.response_model,
                      r.created_at,p.confidence,p.reason,p.status,p.atoms_json FROM ai_peak_predictions p
                      JOIN ai_annotation_results r ON r.id=p.result_id WHERE p.spectrum_id=? AND p.applied=1
                      ORDER BY r.created_at,r.id''',(spectrum_id,)).fetchall()
    return {r['peak_id']:{k:r[k] for k in r.keys() if k!='peak_id'} for r in rows}


def statistics(c):
    sources=[dict(r) for r in c.execute('''WITH per_spectrum AS (
                  SELECT spectrum_id,CASE WHEN origin LIKE 'ai:%' THEN 'ai'
                  WHEN origin IN ('manual','manual_v2') THEN 'manual'
                  WHEN origin='nmrshiftdb2' THEN 'shiftdb2' ELSE 'other' END label_source,count(*) n
                  FROM annotations WHERE status IN ('assigned','equivalent_group','overlap','source_group')
                  GROUP BY spectrum_id,label_source)
                  SELECT a.label_source source,s.nucleus,sum(a.n) peaks,count(DISTINCT s.molecule_id) molecules
                  FROM per_spectrum a JOIN spectra s ON s.id=a.spectrum_id
                  WHERE 1=1
                  AND NOT EXISTS(SELECT 1 FROM spectrum_quality q WHERE q.spectrum_id=s.id AND q.decision!='passed')
                  AND NOT EXISTS(SELECT 1 FROM curation_state cs WHERE cs.scope='spectrum' AND cs.target_id=s.id AND cs.state!='active')
                  GROUP BY a.label_source,s.nucleus''')]
    models=[dict(r) for r in c.execute('''SELECT model,provider,count(*) results,count(DISTINCT molecule_id) molecules,
                  sum(applied_count) published_peaks FROM ai_annotation_results GROUP BY provider,model''')]
    for row in models:
        row['current_ai_peaks']=c.execute('''SELECT count(*) FROM annotations a JOIN spectra s ON s.id=a.spectrum_id
            WHERE a.origin=? AND a.status IN ('assigned','equivalent_group','overlap')
            AND NOT EXISTS(SELECT 1 FROM spectrum_quality q WHERE q.spectrum_id=s.id AND q.decision!='passed')
            AND NOT EXISTS(SELECT 1 FROM curation_state cs WHERE cs.scope='spectrum' AND cs.target_id=s.id AND cs.state!='active')''',
            ('ai:'+row['provider']+':'+row['model'],)).fetchone()[0]
    return {'sources':sources,'models':models,
            'note':'来源按当前有效标签统计；模型结果数和历史写入数为累计值，之后的人工改标会降低当前 AI 标签数。'}
