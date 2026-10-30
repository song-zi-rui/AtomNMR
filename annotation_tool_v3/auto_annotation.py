"""Continuous assignment jobs: durable cursors, pause/resume and atomic publication."""
import hashlib
import json
import sqlite3
import threading
import time
import uuid
from contextlib import closing
from datetime import datetime,timezone
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED

from . import core as db, ai_annotation as ai, ai_database, ai_parallel


def initialize(c):
    c.executescript('''
    CREATE TABLE IF NOT EXISTS auto_jobs(
      id TEXT PRIMARY KEY,state TEXT NOT NULL,scope TEXT NOT NULL,config_json TEXT NOT NULL,
      phase INTEGER NOT NULL DEFAULT 0,cursor TEXT NOT NULL DEFAULT '',created_at TEXT NOT NULL,
      updated_at TEXT NOT NULL,error TEXT,skipped INTEGER NOT NULL DEFAULT 0,
      processed INTEGER NOT NULL DEFAULT 0,success INTEGER NOT NULL DEFAULT 0,
      failures INTEGER NOT NULL DEFAULT 0,applied_peaks INTEGER NOT NULL DEFAULT 0);
    CREATE TABLE IF NOT EXISTS auto_items(
      id INTEGER PRIMARY KEY,job_id TEXT NOT NULL REFERENCES auto_jobs(id),molecule_id TEXT NOT NULL,
      formula TEXT NOT NULL,revision INTEGER NOT NULL,input_json TEXT NOT NULL,reference_json TEXT NOT NULL,
      state TEXT NOT NULL,output_json TEXT,usage_json TEXT,error TEXT,diagnostics_json TEXT,
      publication_json TEXT,credential_id TEXT,elapsed REAL,created_at TEXT NOT NULL,updated_at TEXT NOT NULL,
      UNIQUE(job_id,molecule_id));
    CREATE INDEX IF NOT EXISTS auto_item_job ON auto_items(job_id,id);
    CREATE INDEX IF NOT EXISTS auto_item_job_state ON auto_items(job_id,state,id);
    ''')
    columns={r[1] for r in c.execute('PRAGMA table_info(auto_jobs)')}
    for name in ('processed','success','failures','applied_peaks'):
        if name not in columns:c.execute('ALTER TABLE auto_jobs ADD COLUMN '+name+' INTEGER NOT NULL DEFAULT 0')
    if 'request_config_json' not in {r[1] for r in c.execute('PRAGMA table_info(auto_items)')}:
        c.execute('ALTER TABLE auto_items ADD COLUMN request_config_json TEXT')
    if 'deleted_at' not in columns:
        c.execute('ALTER TABLE auto_jobs ADD COLUMN deleted_at TEXT')
    ai.response_columns(c,'auto_items')
    if 'usage_totals_json' not in columns:
        c.execute("ALTER TABLE auto_jobs ADD COLUMN usage_totals_json TEXT NOT NULL DEFAULT '{}'")
        totals={}
        for r in c.execute('SELECT job_id,usage_json FROM auto_items'):
            value=ai.usage_numbers(json.loads(r['usage_json'] or '{}'));total=totals.setdefault(r['job_id'],{})
            for k,v in value.items():total[k]=total.get(k,0)+v
        for identifier,total in totals.items():c.execute('UPDATE auto_jobs SET usage_totals_json=? WHERE id=?',(db.dumps(total),identifier))
    c.execute("UPDATE auto_jobs SET state='interrupted',error='服务重启；已保存进度，点击恢复继续。' WHERE state IN ('running','pausing')")
    c.execute("UPDATE auto_jobs SET failures=failures+(SELECT count(*) FROM auto_items i WHERE i.job_id=auto_jobs.id AND i.state='calling'),processed=processed+(SELECT count(*) FROM auto_items i WHERE i.job_id=auto_jobs.id AND i.state='calling')")
    c.execute("UPDATE auto_items SET state='error',error='请求期间服务中断，结果未知；不会自动重复计费。' WHERE state='calling'")
    c.commit()


class Service:
    def __init__(self,manager):
        self.manager=manager;self.path=manager.path;self.main_path=manager.main_path
        with closing(ai.store(self.path)) as c:initialize(c)

    def config(self):
        return ai.snapshot_config(self.manager.settings()) | self.execution()

    def execution(self,payload=None):
        with self.manager.lock, closing(ai.store(self.path)) as c:
            row=c.execute("SELECT value FROM settings WHERE key='auto_execution'").fetchone()
            saved=ai_parallel.validate(json.loads(row[0]) if row else {})
            if payload is not None:
                saved=ai_parallel.validate(saved | payload)
                if saved['parallel_strategy'] in ('deepseek_offpeak','deepseek_offpeak_only') and not ai_parallel.official(self.manager.settings()):
                    raise ValueError('谷价加速仅适用于 DeepSeek 官方端点；当前端点请选择固定并行')
                c.execute("INSERT OR REPLACE INTO settings VALUES('auto_execution',?)",(db.dumps(saved),))
                job=c.execute("SELECT id,config_json FROM auto_jobs WHERE id=? AND state='running'",(self.manager.active,)).fetchone()
                if job:
                    config=json.loads(job['config_json'])
                    if any(config.get(k,ai_parallel.DEFAULTS[k])!=v for k,v in saved.items()):
                        config['parallel_history']=[*config.get('parallel_history',[]),
                            {'changed_at':db.now(),'previous':ai_parallel.validate(config),'new':saved}]
                        config.update(saved)
                        c.execute('UPDATE auto_jobs SET config_json=?,updated_at=? WHERE id=?',(db.dumps(config),db.now(),job['id']))
                c.commit()
        return saved

    def delete(self,identifier):
        # Hide task history, retaining immutable request/source records and all labels.
        with self.manager.lock:
            with closing(ai.store(self.path)) as c:
                row=c.execute('SELECT state FROM auto_jobs WHERE id=?',(identifier,)).fetchone()
                if not row:raise KeyError(identifier)
                if self.manager.active==identifier or row['state'] in ('running','pausing'):
                    raise ValueError('请先暂停此任务，等待当前请求结束后再删除记录')
                c.execute('UPDATE auto_jobs SET deleted_at=? WHERE id=?',(db.now(),identifier));c.commit()
        return {'deleted':True,'labels_preserved':True}

    def start(self,payload):
        scope=payload.get('scope','pending')
        if scope not in ('pending','partial','unlabeled'):raise ValueError('请选择待补归属或无归属范围')
        with self.manager.lock:
            if self.manager.active:raise RuntimeError('已有运行中的自动任务；请先暂停')
            self.manager.key()
            config=self.config();identifier=uuid.uuid4().hex
            with closing(ai.store(self.path)) as c:
                c.execute('INSERT INTO auto_jobs(id,state,scope,config_json,created_at,updated_at) VALUES(?,?,?,?,?,?)',
                          (identifier,'running',scope,db.dumps(config),db.now(),db.now()));c.commit()
            self.manager.active=identifier;self.manager.cancel.clear()
        self.manager.pool.submit(self.work,identifier)
        return {'job_id':identifier}

    def pause(self,identifier):
        with self.manager.lock:
            if self.manager.active!=identifier:raise ValueError('此任务没有正在运行')
            self.manager.cancel.set()
        with closing(ai.store(self.path)) as c:
            c.execute("UPDATE auto_jobs SET state='pausing',updated_at=? WHERE id=? AND state='running'",(db.now(),identifier));c.commit()
        self.manager.cancel_requests()
        return {'state':'pausing'}

    def resume(self,identifier,use_current_limits=False):
        with self.manager.lock:
            if self.manager.active:raise RuntimeError('已有运行中的自动任务')
            self.manager.key()
            with closing(ai.store(self.path)) as c:
                row=c.execute('SELECT * FROM auto_jobs WHERE id=?',(identifier,)).fetchone()
                if not row:raise KeyError(identifier)
                if row['state']=='completed':raise ValueError('此任务已处理完所选队列；可另开任务检查新入库数据')
                if use_current_limits:
                    old=json.loads(row['config_json']);current=self.config()
                    if any(current[k]!=old[k] for k in ('model','thinking')):
                        raise ValueError('恢复原任务只调整 token 和超时；请保持模型及思考开关一致，换模型请新建任务')
                    keys=('max_tokens','timeout','reasoning_effort','wire_format','assignment_mode','policy','system_prompt','prompt_sha256','base_url','api_adapter','json_mode','stream','provider','tls_policy','ca_file','ca_sha256','tls_origin','parallelism','parallel_strategy')
                    if any(current[k]!=old.get(k) for k in keys):
                        # Preserve old request settings before revising future limits.
                        c.execute("UPDATE auto_items SET request_config_json=? WHERE job_id=? AND request_config_json IS NULL AND state!='pending'",(db.dumps(old),identifier))
                        updated=dict(old)
                        updated['limit_history']=[*old.get('limit_history',[]),{'changed_at':db.now(),
                            'old_max_tokens':old['max_tokens'],'max_tokens':current['max_tokens'],
                            'old_timeout':old['timeout'],'timeout':current['timeout']}]
                        updated['config_history']=[*old.get('config_history',[]),{'changed_at':db.now(),
                            'old':{k:old.get(k) for k in keys if k!='system_prompt'},'new':{k:current[k] for k in keys if k!='system_prompt'}}]
                        updated.update({k:current[k] for k in keys})
                        c.execute('UPDATE auto_jobs SET config_json=? WHERE id=?',(db.dumps(updated),identifier))
                c.execute("UPDATE auto_jobs SET state='running',error=NULL,updated_at=? WHERE id=?",(db.now(),identifier));c.commit()
            self.manager.active=identifier;self.manager.cancel.clear()
        self.manager.pool.submit(self.work,identifier)
        return {'job_id':identifier}

    def _next(self,job):
        if job['scope']=='retry':return None
        phases=['complete_remaining','both_unassigned'] if job['scope']=='pending' else ['complete_remaining'] if job['scope']=='partial' else ['both_unassigned']
        phase,cursor=job['phase'],job['cursor']
        model=json.loads(job['config_json'])['model']
        with closing(db.connect(self.main_path,timeout=2)) as c:
            while phase<len(phases):
                rows=c.execute('''SELECT id FROM molecules WHERE queue=? AND id>?
                    AND NOT EXISTS(SELECT 1 FROM ai_annotation_results r WHERE r.molecule_id=molecules.id AND r.model=?)
                    ORDER BY id LIMIT 32''',(phases[phase],cursor,model)).fetchall()
                if not rows:
                    phase+=1;cursor=''
                    with closing(ai.store(self.path)) as out:
                        out.execute('UPDATE auto_jobs SET phase=?,cursor=?,updated_at=? WHERE id=?',(phase,cursor,db.now(),job['id']));out.commit()
                    continue
                for row in rows:
                    if self.manager.cancel.is_set():return None
                    mid=row['id'];cursor=mid
                    try:
                        c.execute('BEGIN')
                        record=db.record(c,mid,require_quality=True)
                        context,references=ai.production_context(record,json.loads(job['config_json']))
                        c.rollback()
                    except (KeyError,ValueError):
                        c.rollback()
                        with closing(ai.store(self.path)) as out:
                            out.execute('UPDATE auto_jobs SET cursor=?,skipped=skipped+1,updated_at=? WHERE id=?',(cursor,db.now(),job['id']));out.commit()
                        continue
                    with closing(ai.store(self.path)) as out:
                        exists=out.execute('SELECT id FROM auto_items WHERE job_id=? AND molecule_id=?',(job['id'],mid)).fetchone()
                        out.execute('UPDATE auto_jobs SET cursor=?,phase=?,updated_at=? WHERE id=?',(cursor,phase,db.now(),job['id']))
                        if not exists:
                            out.execute("INSERT INTO auto_items(job_id,molecule_id,formula,revision,input_json,reference_json,state,created_at,updated_at) VALUES(?,?,?,?,?,?,'pending',?,?)",
                                        (job['id'],mid,record['formula'],record['revision'],db.dumps(context),db.dumps(references),db.now(),db.now()))
                        out.commit()
                        if not exists:return out.execute('SELECT * FROM auto_items WHERE id=last_insert_rowid()').fetchone()
        return None

    def _publish(self,item,config):
        config=json.loads(item['request_config_json']) if item['request_config_json'] else config
        started=time.perf_counter()
        result=ai_database.publish(self.main_path,'auto:'+item['job_id']+':'+str(item['id']),item['job_id'],item['molecule_id'],item['revision'],
                                   json.loads(item['input_json']),json.loads(item['output_json']),json.loads(item['reference_json']),config)
        diagnostics=json.loads(item['diagnostics_json'] or 'null') or {}
        diagnostics.setdefault('timing',{})['publish_ms']=(time.perf_counter()-started)*1000
        with closing(ai.store(self.path)) as c:
            current=c.execute('SELECT state FROM auto_items WHERE id=?',(item['id'],)).fetchone()
            if current['state']!='done':
                c.execute('UPDATE auto_jobs SET processed=processed+1,success=success+1,applied_peaks=applied_peaks+?,updated_at=? WHERE id=?',
                          (result['applied'],db.now(),item['job_id']))
            c.execute("UPDATE auto_items SET state='done',publication_json=?,diagnostics_json=?,updated_at=? WHERE id=?",(db.dumps(result),db.dumps(diagnostics),db.now(),item['id']));c.commit()
        return result

    def _execute(self,item,config):
        identifier=item['job_id']
        recovered=item['state']=='responded';key=''
        if recovered:config=json.loads(item['request_config_json'])
        else:
            key=self.manager.key();credential=self.manager.credential
            if self.manager.cancel.is_set():return None
            config=dict(config,execution_at_dispatch=ai_parallel.capacity(config),dispatched_at=db.now())
            with closing(ai.store(self.path)) as c:
                c.execute("UPDATE auto_items SET state='calling',credential_id=?,request_config_json=?,updated_at=? WHERE id=?",
                          (credential['id'] if credential else None,db.dumps(config),db.now(),item['id']));c.commit()
            self.manager.vault.mark(credential['id'] if credential else None,used=True)
        started=time.monotonic();usage,diagnostics={},None
        try:
            if recovered:
                envelope=json.loads(item['raw_response_json']);usage=envelope['usage']
                output=ai.decode_response(envelope,config,json.loads(item['input_json']),recover=True)
                if config.get('local_response_reuse'):
                    usage={};output['api_diagnostics']['response_reused_from']=config['local_response_reuse']
            else:
                def progress(value):
                    with closing(ai.store(self.path)) as c:
                        c.execute('UPDATE auto_items SET diagnostics_json=? WHERE id=?',(db.dumps({'progress':{**value,'progress_at':db.now()}}),item['id']));c.commit()
                with self.manager.request_scope(),ai.capture_response(lambda e:ai.archive_response(self.path,'auto_items',item['id'],e,'responded')),ai.capture_progress(progress):
                    output,usage=ai.complete(key,config,json.loads(item['input_json']))
            diagnostics=output.get('api_diagnostics')
            # Persist the response before touching the main database; recovery never re-calls it.
            with closing(ai.store(self.path)) as c:
                c.execute('BEGIN IMMEDIATE');ai.update_auto_usage(c,item['id'],usage)
                c.execute("UPDATE auto_items SET state='received',output_json=?,usage_json=?,diagnostics_json=?,elapsed=?,updated_at=? WHERE id=?",
                          (db.dumps(output),db.dumps(usage),db.dumps(diagnostics),time.monotonic()-started,db.now(),item['id']));c.commit()
                item=c.execute('SELECT * FROM auto_items WHERE id=?',(item['id'],)).fetchone()
            return None
        except Exception as exc:
            message=(str(exc).replace(key,'[REDACTED]') if key else str(exc))[:1500]
            failure_usage={} if config.get('local_response_reuse') else getattr(exc,'usage',usage)
            failure_diagnostics=getattr(exc,'diagnostics',diagnostics)
            if config.get('local_response_reuse') and failure_diagnostics:
                failure_diagnostics=dict(failure_diagnostics,response_reused_from=config['local_response_reuse'])
            with closing(ai.store(self.path)) as c:
                c.execute('BEGIN IMMEDIATE');ai.update_auto_usage(c,item['id'],failure_usage)
                c.execute('UPDATE auto_jobs SET processed=processed+1,failures=failures+1,updated_at=? WHERE id=?',(db.now(),identifier))
                c.execute("UPDATE auto_items SET state='error',error=?,usage_json=?,diagnostics_json=?,elapsed=?,updated_at=? WHERE id=?",
                          (message,db.dumps(failure_usage),db.dumps(failure_diagnostics),time.monotonic()-started,db.now(),item['id']));c.commit()
            return exc

    def work(self,identifier):
        final,error,streak='paused',None,0
        inflight={};exhausted=False;halt=False
        pool=ThreadPoolExecutor(max_workers=ai_parallel.MAX_WORKERS,thread_name_prefix='nmr-molecule')
        try:
            with closing(ai.store(self.path)) as c:
                job=c.execute('SELECT * FROM auto_jobs WHERE id=?',(identifier,)).fetchone()
                config=json.loads(job['config_json'])
            while True:
                with closing(ai.store(self.path)) as c:
                    config=json.loads(c.execute('SELECT config_json FROM auto_jobs WHERE id=?',(identifier,)).fetchone()[0])
                limit=ai_parallel.capacity(config)['effective_parallelism']
                while not halt and not self.manager.cancel.is_set() and len(inflight)<limit:
                    with closing(ai.store(self.path)) as c:
                        job=c.execute('SELECT * FROM auto_jobs WHERE id=?',(identifier,)).fetchone()
                        excluded=tuple(inflight.values())
                        exclusion=(' AND id NOT IN ('+','.join('?' for _ in excluded)+')') if excluded else ''
                        item=c.execute("SELECT * FROM auto_items WHERE job_id=? AND state IN ('pending','responded','received')"+exclusion+' ORDER BY id LIMIT 1',(identifier,*excluded)).fetchone()
                    if item and item['state']=='received':
                        self._publish(item,config);continue
                    if not item and not exhausted:
                        item=self._next(job)
                        if not item:exhausted=True
                    if not item:break
                    if self.manager.cancel.is_set():break
                    inflight[pool.submit(self._execute,item,config)]=item['id']
                if not inflight:
                    if limit==0 and not halt and not self.manager.cancel.is_set():
                        # Keep the job alive across peak hours, with no paid calls
                        # and no cursor advancement. Cancellation remains immediate.
                        self.manager.cancel.wait(2)
                        continue
                    if not halt:final='paused' if self.manager.cancel.is_set() else 'completed'
                    break
                done,_=wait(inflight,timeout=.2,return_when=FIRST_COMPLETED)
                for future in done:
                    item_id=inflight.pop(future)
                    try:
                        exc=future.result()
                        with closing(ai.store(self.path)) as c:
                            item=c.execute('SELECT * FROM auto_items WHERE id=?',(item_id,)).fetchone()
                        if item['state']=='received':self._publish(item,config)
                        if exc is None:streak=0
                        elif isinstance(exc,ai.APIRequestCancelled):pass
                        else:
                            streak+=1
                            if streak>=3 or ai.should_pause_batch(exc):
                                halt=True;final='attention';error=item['error']
                    except Exception as exc:
                        halt=True;final='attention';error=str(exc)[:1500]
                # On failure, stop dispatch but retain results of already-issued
                # requests. User pause cancels all local waits through Manager.
        except Exception as exc:final,error='attention',str(exc)[:1500]
        finally:
            pool.shutdown(wait=True)
            with closing(ai.store(self.path)) as c:
                c.execute('UPDATE auto_jobs SET state=?,error=?,updated_at=? WHERE id=?',(final,error,db.now(),identifier));c.commit()
            with self.manager.lock:self.manager.active=None

    def status(self,identifier=None,result_state=None):
        if result_state not in (None,'','error'):raise ValueError('结果筛选无效')
        with closing(ai.store(self.path)) as c:
            jobs=[dict(r) for r in c.execute('SELECT * FROM auto_jobs WHERE deleted_at IS NULL ORDER BY created_at DESC,id DESC LIMIT 30')]
            job=c.execute('SELECT * FROM auto_jobs WHERE id=? AND deleted_at IS NULL',(identifier,)).fetchone() if identifier else c.execute('SELECT * FROM auto_jobs WHERE deleted_at IS NULL ORDER BY created_at DESC,id DESC LIMIT 1').fetchone()
            if not job:return {'jobs':jobs,'job':None,'active_id':self.manager.active}
            value=dict(job);value['config']=json.loads(value.pop('config_json'))
            value['usage']=json.loads(value.pop('usage_totals_json'))
            value['execution']=ai_parallel.capacity(value['config'])
            items=c.execute('SELECT * FROM auto_items WHERE job_id=?'+(" AND state='error'" if result_state=='error' else '')+' ORDER BY id DESC LIMIT '+('100' if result_state=='error' else '30'),(job['id'],)).fetchall()
            value['counts']={'done':job['success'],'error':job['failures']}
            for state in ('pending','calling','responded','received'):
                value['counts'][state]=c.execute('SELECT count(*) FROM auto_items WHERE job_id=? AND state=?',(job['id'],state)).fetchone()[0]
            value['items']=[{k:r[k] for k in ('id','molecule_id','formula','state','error','elapsed','updated_at')}|{'publication':json.loads(r['publication_json'] or '{}')} for r in items]
            calls=c.execute("SELECT * FROM auto_items WHERE job_id=? AND state='calling' ORDER BY id",(job['id'],)).fetchall()
            value['current_requests']=[]
            value['current_request']=None
            for calling in calls:
                cfg=json.loads(calling['request_config_json'] or job['config_json'])
                elapsed=max(0,(datetime.now(timezone.utc)-datetime.fromisoformat(calling['updated_at'])).total_seconds())
                value['current_requests'].append({'item_id':calling['id'],'formula':calling['formula'],'model':cfg['model'],
                                          'started_at':calling['updated_at'],'elapsed_s':elapsed,'timeout_s':cfg['timeout'],
                                          'stream':cfg.get('stream',False),'progress':(json.loads(calling['diagnostics_json'] or '{}') or {}).get('progress') or {}})
            if value['current_requests']:value['current_request']=value['current_requests'][0]
        for j in jobs:j.pop('config_json',None)
        return {'jobs':jobs,'job':value,'active_id':self.manager.active}

    def retry_failed(self,identifier,item_id=None,mode='failed'):
        """New finite retry job; successful source items and the continuous cursor stay intact."""
        if item_id is not None and (type(item_id) is not int or item_id<1):raise ValueError('失败项目编号无效')
        if mode not in ('failed','local','issues','uncertain'):raise ValueError('恢复模式无效')
        if mode in ('issues','uncertain') and item_id is None:raise ValueError('请选择需要局部重试的分子')
        with self.manager.lock:
            if self.manager.active:raise RuntimeError('请先暂停当前任务并等待正在执行的请求完成')
            if mode!='local':self.manager.key()
            config=self.config();new_id=uuid.uuid4().hex
            with closing(ai.store(self.path)) as c:
                original=c.execute('SELECT * FROM auto_jobs WHERE id=?',(identifier,)).fetchone()
                if not original:raise KeyError(identifier)
                sql="SELECT * FROM auto_items WHERE job_id=? AND state="+ ("'done'" if mode in ('issues','uncertain') else "'error'")
                rows=c.execute(sql+(' AND id=?' if item_id is not None else '')+' ORDER BY id',
                               (identifier,item_id) if item_id is not None else (identifier,)).fetchall()
                if not rows:raise ValueError('没有符合条件的失败项；只支持重试此任务中的失败记录')
                if mode=='local':
                    rows=[r for r in rows if r['raw_response_json'] and json.loads(r['raw_response_json'])['diagnostics'].get('finish_reason')=='stop']
                    if not rows:raise ValueError('没有已完整保存的原始响应；旧失败项无法凭空恢复，截断结果仍需重新请求')
                    config=json.loads(original['config_json'])
                if mode=='issues' and not json.loads(rows[0]['output_json']).get('repair_targets'):raise ValueError('该结果没有待修复的格式错误峰')
                config.update(retry_of=identifier,retry_item_ids=[r['id'] for r in rows],
                              snapshot_basis='original_failed_input_and_reference',retry_mode=mode)
                c.execute("INSERT INTO auto_jobs(id,state,scope,config_json,created_at,updated_at) VALUES(?,'running','retry',?,?,?)",
                          (new_id,db.dumps(config),db.now(),db.now()))
                for r in rows:
                    context=json.loads(r['input_json']);references=json.loads(r['reference_json'])
                    if mode=='issues':
                        context['requested_peaks']=json.loads(r['output_json'])['repair_targets']
                    if mode=='uncertain':
                        wanted=ai.protocol.requested(context)
                        context['requested_peaks']=[{'spectrum_id':p['spectrum_id'],'peak_id':p['peak_id']} for p in json.loads(r['output_json'])['assignments']
                            if p['status'] in ('ambiguous','unassigned') and (p['spectrum_id'],p['peak_id']) in wanted]
                        if not context['requested_peaks']:raise ValueError('该结果没有可重分析的歧义或未归属峰')
                    if mode in ('issues','uncertain'):
                        wanted=ai.protocol.requested(context)
                        references=[ref for ref in references if ref['kind']=='peak' and (ref['spectrum_id'],ref['peak_id']) in wanted]
                    c.execute("INSERT INTO auto_items(job_id,molecule_id,formula,revision,input_json,reference_json,state,created_at,updated_at) VALUES(?,?,?,?,?,?,'pending',?,?)",
                              (new_id,r['molecule_id'],r['formula'],r['revision'],db.dumps(context),db.dumps(references),db.now(),db.now()))
                    if mode=='local':
                        request=json.loads(r['request_config_json'] or original['config_json'])
                        request['local_response_reuse']=r['id']
                        c.execute("UPDATE auto_items SET state='responded',raw_response_json=?,request_config_json=? WHERE id=last_insert_rowid()",
                                  (r['raw_response_json'],db.dumps(request)))
                c.commit()
            self.manager.active=new_id;self.manager.cancel.clear()
        self.manager.pool.submit(self.work,new_id)
        return {'job_id':new_id,'retry_of':identifier,'selected':len(rows),'mode':mode}

    def item(self,identifier):
        with closing(ai.store(self.path)) as c:
            row=c.execute('SELECT * FROM auto_items WHERE id=?',(int(identifier),)).fetchone()
        if not row:raise KeyError(identifier)
        item=dict(row)
        for key in ('input','reference','output','usage','diagnostics','publication','request_config','raw_response'):
            item[key]=json.loads(item.pop(key+'_json') or ('[]' if key=='reference' else 'null'))
        return item

    def analyze(self,identifier):
        # Explicit analysis, independent of the running queue and its write transactions.
        with closing(ai.store(self.path)) as c:
            rows=c.execute('SELECT state,output_json,reference_json,input_json FROM auto_items WHERE job_id=?',(identifier,))
            def independent():
                for row in rows:
                    value=dict(row);context=json.loads(row['input_json'])
                    anchors={(r['spectrum_id'],r['peak_id']) for r in context.get('anchors',[])};wanted=ai.protocol.requested(context)
                    refs=json.loads(row['reference_json'])
                    refs=[r for r in refs if (r['kind']=='peak' and (r['spectrum_id'],r['peak_id']) not in anchors and (r['spectrum_id'],r['peak_id']) in wanted) or (r['kind']=='group' and all((r['spectrum_id'],p) not in anchors and (r['spectrum_id'],p) in wanted for p in r['peak_ids']))]
                    value['reference_json']=db.dumps(refs);yield value
            metrics=ai.summarize(independent())
        return {'metrics':metrics,'scope':'reference_snapshots_before_AI_write','job_id':identifier,
                'note':'只对未提供给模型的参考标签计分；生产锚点及 AI 标签不作为独立标准答案。'}

    def molecule_results(self,mid):
        with closing(db.connect(self.main_path,timeout=2)) as c:
            return {'results':[dict(r) for r in c.execute('SELECT id,model,provider,response_model,created_at,applied_count FROM ai_annotation_results WHERE molecule_id=? ORDER BY created_at DESC',(mid,))]}

    def result(self,identifier):
        with closing(db.connect(self.main_path,timeout=2)) as c:
            row=c.execute('SELECT * FROM ai_annotation_results WHERE id=?',(identifier,)).fetchone()
        if not row:raise KeyError(identifier)
        value=dict(row)
        for key in ('input','output','reference'):value[key]=json.loads(value.pop(key+'_json'))
        value['publication']={'applied':value['applied_count'],'skipped':value['skipped_count']}
        return value

    def import_history(self):
        totals={'results':0,'applied':0,'skipped':0}
        with closing(ai.store(self.path)) as c:
            rows=c.execute("SELECT i.*,r.config_json FROM items i JOIN runs r ON r.id=i.run_id WHERE i.state='done' AND i.output_json IS NOT NULL ORDER BY i.id").fetchall()
        for r in rows:
            if json.loads(r['config_json']).get('purpose')=='evaluation':continue
            result=ai_database.publish(self.main_path,'legacy:'+r['run_id']+':'+str(r['id']),r['run_id'],r['molecule_id'],r['revision'],
                                      json.loads(r['input_json']),json.loads(r['output_json']),json.loads(r['reference_json']),json.loads(r['config_json']))
            totals['results']+=1;totals['applied']+=result['applied'];totals['skipped']+=result['skipped']
        return totals

    def timing(self,identifier=None):
        from .ai_metrics import timings
        with closing(ai.store(self.path)) as c:
            jobs=c.execute('SELECT id,state,scope,created_at,updated_at FROM auto_jobs WHERE deleted_at IS NULL'+(' AND id=?' if identifier else '')+' ORDER BY created_at', (identifier,) if identifier else ()).fetchall()
            result=[]
            for j in jobs:
                items=c.execute('SELECT state,elapsed,diagnostics_json FROM auto_items WHERE job_id=?',(j['id'],)).fetchall()
                result.append(dict(j)|{'all':timings(r['elapsed'] for r in items),'successful':timings(r['elapsed'] for r in items if r['state']=='done'),
                                      'failed':timings(r['elapsed'] for r in items if r['state']=='error')})
        return {'jobs':result,'note':'请求处理用时包含接口等待与本地解析；并行任务的请求耗时会重叠，合计不是批次实际用时。历史记录无法拆分首token等待。'}
