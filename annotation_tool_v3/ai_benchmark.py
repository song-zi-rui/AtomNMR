"""Frozen, blind model experiments in a separate database. No publication path."""
import copy
import hashlib
import json
import math
import random
import re
import sqlite3
import time
import uuid
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from contextlib import closing
from datetime import datetime,timezone
from pathlib import Path

from rdkit import Chem
from . import core as db,ai_annotation as ai,ai_protocol as protocol,ai_parallel
from .ai_metrics import timings


def connect(path):
    c=sqlite3.connect(path,timeout=15);c.row_factory=sqlite3.Row;c.execute('PRAGMA foreign_keys=ON');return c


def initialize(path):
    with closing(connect(path)) as c:
        c.execute('PRAGMA journal_mode=WAL')
        c.executescript('''
        CREATE TABLE IF NOT EXISTS profiles(id TEXT PRIMARY KEY,name TEXT NOT NULL,config_json TEXT NOT NULL,credential_id TEXT,updated_at TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS sample_sets(id TEXT PRIMARY KEY,name TEXT NOT NULL,created_at TEXT NOT NULL,options_json TEXT NOT NULL,summary_json TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS samples(set_id TEXT NOT NULL REFERENCES sample_sets(id),ordinal INTEGER NOT NULL,molecule_id TEXT NOT NULL,
          formula TEXT NOT NULL,revision INTEGER NOT NULL,input_json TEXT NOT NULL,reference_json TEXT NOT NULL,features_json TEXT NOT NULL,input_hash TEXT NOT NULL,
          PRIMARY KEY(set_id,ordinal),UNIQUE(set_id,molecule_id));
        CREATE TABLE IF NOT EXISTS experiments(id TEXT PRIMARY KEY,name TEXT NOT NULL,set_id TEXT NOT NULL REFERENCES sample_sets(id),state TEXT NOT NULL,
          config_json TEXT NOT NULL,created_at TEXT NOT NULL,updated_at TEXT NOT NULL,error TEXT);
        CREATE TABLE IF NOT EXISTS cases(id INTEGER PRIMARY KEY,run_id TEXT NOT NULL REFERENCES experiments(id),ordinal INTEGER NOT NULL,
          molecule_id TEXT NOT NULL,formula TEXT NOT NULL,state TEXT NOT NULL,output_json TEXT,elapsed REAL,error TEXT,UNIQUE(run_id,ordinal));
        CREATE TABLE IF NOT EXISTS attempts(id INTEGER PRIMARY KEY,case_id INTEGER NOT NULL REFERENCES cases(id),stage INTEGER NOT NULL,
          profile_json TEXT NOT NULL,input_json TEXT NOT NULL,state TEXT NOT NULL,raw_json TEXT,output_json TEXT,usage_json TEXT,
          diagnostics_json TEXT,elapsed REAL,error TEXT,created_at TEXT NOT NULL,UNIQUE(case_id,stage));
        CREATE TABLE IF NOT EXISTS sessions(id INTEGER PRIMARY KEY,run_id TEXT NOT NULL REFERENCES experiments(id),started_at TEXT NOT NULL,
          ended_at TEXT,elapsed REAL);
        CREATE INDEX IF NOT EXISTS case_run ON cases(run_id,id);
        CREATE INDEX IF NOT EXISTS attempt_case ON attempts(case_id,stage);
        ''')
        c.execute("UPDATE experiments SET state='interrupted',error='服务重启；未确认的请求不自动重复。' WHERE state IN ('running','stopping')")
        c.execute("UPDATE attempts SET state='unknown',error='请求结果未确认；不会自动重试' WHERE state='calling'")
        c.execute("UPDATE cases SET state='error',error='调用结果未确认，请另开对比轮次' WHERE state='calling' AND EXISTS(SELECT 1 FROM attempts a WHERE a.case_id=cases.id AND a.state='unknown')")
        c.commit()


def fingerprint(value):return hashlib.sha256(db.dumps(value).encode()).hexdigest()


class Service:
    def __init__(self,manager):
        self.manager=manager;self.path=Path(manager.main_path).parent/'ai_benchmarks.sqlite'
        initialize(self.path)
        with closing(connect(self.path)) as c:
            if not c.execute('SELECT 1 FROM profiles').fetchone():
                config={k:manager.settings()[k] for k in ai.DEFAULTS}
                config['assignment_mode']='blind'
                c.execute('INSERT INTO profiles VALUES(?,?,?,?,?)',('official-current','现有官方 DeepSeek',db.dumps(config),manager.credential['id'] if manager.credential else None,db.now()))
                config.update(model='Qwen3.5-122B-A10B',base_url='',api_adapter='qwen-vllm',thinking=True)
                c.execute('INSERT INTO profiles VALUES(?,?,?,?,?)',('school-qwen','学校 Qwen（待填写接口地址）',db.dumps(config),None,db.now()));c.commit()

    def profiles(self):
        rows=self.manager.vault.status()['credentials'];keys={r['id']:r for r in rows}
        with closing(connect(self.path)) as c:values=[dict(r) for r in c.execute('SELECT * FROM profiles ORDER BY id')]
        for v in values:
            v['config']=json.loads(v.pop('config_json'));key=keys.get(v['credential_id'])
            v['masked_key']=key['masked_key'] if key else None
            v['ready']=bool(v['config'].get('base_url') and key and key['enabled'])
        return {'profiles':values,'credentials':rows}

    def save_profile(self,payload):
        with self.manager.lock:
            if self.manager.active:raise ValueError('请等待当前请求结束后再修改端点')
            identifier=payload.get('id') or uuid.uuid4().hex;name=str(payload.get('name','')).strip()
            if not name or len(name)>100:raise ValueError('请填写端点名称')
            with closing(connect(self.path)) as c:
                old=c.execute('SELECT * FROM profiles WHERE id=?',(identifier,)).fetchone()
                current=json.loads(old['config_json']) if old else dict(ai.DEFAULTS)
                current.update({k:v for k,v in payload.get('config',{}).items() if k in ai.DEFAULTS})
                if old and current.get('ca_file')!=json.loads(old['config_json']).get('ca_file') and 'ca_sha256' not in payload.get('config',{}):current['ca_sha256']=''
                current['assignment_mode']='blind';config=self.validate_config(current,allow_empty_url=True)
                credential=old['credential_id'] if old else None
                selected_explicit='credential_id' in payload
                if selected_explicit:
                    credential=payload['credential_id'] or None
                    available={r['id'] for r in self.manager.vault.status()['credentials'] if r['enabled']}
                    if credential and credential not in available:raise ValueError('所选密钥不可用，请重新选择')
                secret=payload.get('api_key','')
                if not selected_explicit and self.manager.credential and credential==self.manager.credential['id'] and old and config['base_url']!=json.loads(old['config_json']).get('base_url'):
                    credential=None # A test endpoint cannot silently receive the production key.
                if secret:
                    provider='qwen' if config['api_adapter'].startswith('qwen') else 'deepseek' if config['api_adapter'].startswith('deepseek') else 'openai-compatible'
                    if self.manager.credential and credential==self.manager.credential['id']:
                        old_secret,_,_=self.manager.vault.read(credential)
                        if secret!=old_secret:credential=None # A test profile cannot rotate the production credential.
                        else:secret=''
                if secret:
                    before={r['id'] for r in self.manager.vault.status()['credentials']}
                    status=self.manager.vault.change({'action':'replace' if credential else 'create','id':credential,'name':name,
                        'provider':provider,'api_key':secret,'activate':False})
                    if not credential:credential=next(r['id'] for r in status['credentials'] if r['id'] not in before)
                    self.manager.refresh_credential()
                c.execute('INSERT OR REPLACE INTO profiles VALUES(?,?,?,?,?)',(identifier,name,db.dumps(config),credential,db.now()));c.commit()
        return {'id':identifier,**self.profiles()}

    @staticmethod
    def validate_config(config,allow_empty_url=False):
        cfg={k:config.get(k,v) for k,v in ai.DEFAULTS.items()}
        if cfg['base_url'] or not allow_empty_url:cfg['base_url']=ai.normalize_base_url(cfg['base_url'])
        if not isinstance(cfg['model'],str) or not re.fullmatch(r'[a-zA-Z0-9._:/-]{1,120}',cfg['model']):raise ValueError('模型 ID 无效')
        if cfg['api_adapter'] not in ai.API_ADAPTERS:raise ValueError('参数适配无效')
        if type(cfg['max_tokens']) is not int or not 1024<=cfg['max_tokens']<=65536:raise ValueError('输出上限应为1024–65536')
        if type(cfg['timeout']) is not int or not 30<=cfg['timeout']<=600:raise ValueError('超时应为30–600秒')
        if any(type(cfg[k]) is not bool for k in ('thinking','json_mode','stream')):raise ValueError('思考/JSON/流式开关无效')
        if cfg['reasoning_effort'] not in ('low','high','auto') or cfg['wire_format'] not in (*protocol.SUPPORTED_VERSIONS,'verbose-v1'):raise ValueError('思考强度/协议无效')
        cfg.update(ai.api_tls.validate_config(cfg))
        cfg['assignment_mode']='blind';return cfg

    def profile_snapshot(self,identifier,overrides=None):
        with closing(connect(self.path)) as c:row=c.execute('SELECT * FROM profiles WHERE id=?',(identifier,)).fetchone()
        if not row:raise ValueError('请先选择端点')
        allowed={'model','thinking','max_tokens','timeout','reasoning_effort','wire_format','json_mode','stream'}
        if any(k not in allowed for k in (overrides or {})):raise ValueError('测试覆盖参数不允许改变端点地址或鉴权')
        cfg=json.loads(row['config_json']);cfg.update(overrides or {})
        cfg=self.validate_config(cfg)
        if not row['credential_id']:raise ValueError('此端点尚未保存密钥；不会借用正式标注密钥')
        return {'id':row['id'],'name':row['name'],'credential_id':row['credential_id'],'config':ai.snapshot_config(cfg,experiment=True)}

    def test_profile(self,identifier):
        profile=self.profile_snapshot(identifier);secret,_,error=self.manager.vault.read(profile['credential_id'])
        if error or not secret:raise ValueError(error or '端点密钥不可用')
        response=ai.request_api(secret,'GET','/models',base_url=profile['config']['base_url'],**ai.api_tls.request_options(profile['config']))
        models=[x['id'] for x in response.get('data',[]) if isinstance(x.get('id'),str)]
        return {'models':models,'selected_model_found':profile['config']['model'] in models,'note':'只查询模型列表，未发起归属分析'}

    def adopt_for_production(self,identifier):
        profile=self.profile_snapshot(identifier)
        with self.manager.lock:
            if self.manager.active:raise ValueError('请先暂停AI任务并等待当前请求返回')
            key,_,error=self.manager.vault.read(profile['credential_id'])
            if error or not key:raise ValueError(error or '端点密钥不可用')
            settings={k:profile['config'][k] for k in ai.DEFAULTS};settings['assignment_mode']='assisted'
            with closing(ai.store(self.manager.path)) as c:
                c.execute('BEGIN IMMEDIATE')
                c.execute("INSERT OR REPLACE INTO settings VALUES('provider',?)",(db.dumps(settings),))
                c.execute("INSERT OR REPLACE INTO settings VALUES('active_credential',?)",(db.dumps(profile['credential_id']),));c.commit()
            self.manager.refresh_credential()
        return self.manager.settings()

    def sets(self):
        with closing(connect(self.path)) as c:rows=[dict(r) for r in c.execute('SELECT * FROM sample_sets ORDER BY created_at DESC,id')]
        for r in rows:r['options']=json.loads(r.pop('options_json'));r['summary']=json.loads(r.pop('summary_json'))
        return {'sets':rows}

    def sample_set(self,identifier):
        with closing(connect(self.path)) as c:
            row=c.execute('SELECT * FROM sample_sets WHERE id=?',(identifier,)).fetchone()
            if not row:raise KeyError(identifier)
            values=[dict(r) for r in c.execute('SELECT ordinal,molecule_id,formula,revision,features_json,input_hash FROM samples WHERE set_id=? ORDER BY ordinal',(identifier,))]
        for v in values:v['features']=json.loads(v.pop('features_json'))
        return {'id':identifier,'name':row['name'],'summary':json.loads(row['summary_json']),'samples':values}

    def build_set(self,payload):
        count=payload.get('count',100);seed=payload.get('seed',42);same_solvent=payload.get('same_solvent',False)
        if type(count) is not int or not 1<=count<=300 or type(seed) is not int:raise ValueError('样本数应为1–300，种子应为整数')
        name=str(payload.get('name') or 'ShiftDB 双谱分层样本').strip()[:120]
        # Use the queue index, avoiding a scan over the much larger unassigned corpus.
        with closing(sqlite3.connect(Path(self.manager.main_path).resolve().as_uri()+'?mode=ro',uri=True)) as c:
            c.row_factory=sqlite3.Row;c.execute('PRAGMA query_only=ON');c.execute('BEGIN')
            candidates=c.execute("SELECT id,smiles,connectivity FROM molecules INDEXED BY molecule_queue WHERE queue IN ('complete','complete_remaining') AND origin='nmrshiftdb2' AND c_assigned>0 AND h_assigned>0 ORDER BY id").fetchall()
            buckets=defaultdict(list);seen=set();rejected=Counter()
            for candidate in candidates:
                try:
                    record=db.record(c,candidate['id'],require_quality=True)
                    # Evaluate one explicitly chosen pair, not a mixture of repeated experiments.
                    context,refs=ai.prepare(record);by_id={s['spectrum_id']:s for s in context['spectra']}
                    ref_counts=Counter(r['spectrum_id'] for r in refs if r['source']=='shiftdb2' and r['status'] in ai.HARD and r['atom_indices'])
                    choices=[p for p in context['pairs'] if ref_counts[p['c_id']] and ref_counts[p['h_id']]]
                    if same_solvent:choices=[p for p in choices if p['solvent_match']=='same_known']
                    if not choices:raise ValueError('无双核种可评分谱对或不满足溶剂条件')
                    pair=max(choices,key=lambda p:(p['solvent_match']=='same_known',ref_counts[p['c_id']]+ref_counts[p['h_id']],p['c_id'],p['h_id']))
                    keep={pair['c_id'],pair['h_id']};context['spectra']=[s for s in context['spectra'] if s['spectrum_id'] in keep];context['pairs']=[pair]
                    refs=[r for r in refs if r['spectrum_id'] in keep and r['source']=='shiftdb2']
                    mol=Chem.MolFromSmiles(candidate['smiles']);connection=candidate['connectivity'] or candidate['smiles']
                    if connection in seen:continue
                    seen.add(connection);atoms=len(context['atoms']);aromatic=any(a.GetIsAromatic() for a in mol.GetAtoms());halogen=any(a.GetSymbol() in ('F','Cl','Br','I') for a in mol.GetAtoms())
                    stratum=('small' if atoms<=15 else 'medium' if atoms<=25 else 'large')+(':aromatic' if aromatic else ':aliphatic')+(':halogen' if halogen else ':other')
                    features={'heavy_atoms':atoms,'aromatic':aromatic,'halogen':halogen,'stratum':stratum,'solvent_match':pair['solvent_match'],
                              'peaks':sum(len(s['peaks']) for s in context['spectra']),'c_reference':ref_counts[pair['c_id']],'h_reference':ref_counts[pair['h_id']]}
                    buckets[stratum].append({'molecule_id':record['id'],'formula':record['formula'],'revision':record['revision'],'input':context,'references':refs,'features':features})
                except (ValueError,KeyError):rejected['not_eligible']+=1
        total=sum(map(len,buckets.values()))
        if not total:raise ValueError('没有满足双核种原库参考条件的分子；可取消同溶剂限制')
        target=min(count,total);rng=random.Random(seed)
        quotas={k:math.floor(target*len(v)/total) for k,v in buckets.items()}
        order=sorted(buckets,key=lambda k:(-(target*len(buckets[k])/total-quotas[k]),k))
        for k in order[:target-sum(quotas.values())]:quotas[k]+=1
        selected=[]
        for k in sorted(buckets):rng.shuffle(buckets[k]);selected.extend(buckets[k][:quotas[k]])
        rng.shuffle(selected);identifier=uuid.uuid4().hex
        summary={'eligible_pool':total,'selected':len(selected),'requested':count,'population_strata':{k:len(v) for k,v in buckets.items()},
                 'sample_strata':dict(Counter(r['features']['stratum'] for r in selected)),'rejected':dict(rejected),
                 'basis':'比例分层、连接层去重、双核种有可评分原库归属；代表该可评分子集，不代表全世界分子分布'}
        with closing(connect(self.path)) as c:
            c.execute('INSERT INTO sample_sets VALUES(?,?,?,?,?)',(identifier,name,db.now(),db.dumps({'count':count,'seed':seed,'same_solvent':same_solvent}),db.dumps(summary)))
            for i,r in enumerate(selected):c.execute('INSERT INTO samples VALUES(?,?,?,?,?,?,?,?,?)',(identifier,i,r['molecule_id'],r['formula'],r['revision'],db.dumps(r['input']),db.dumps(r['references']),db.dumps(r['features']),fingerprint(r['input'])))
            c.commit()
        return self.sample_set(identifier)

    def history(self):
        with closing(connect(self.path)) as c:rows=[dict(r) for r in c.execute('SELECT id,name,set_id,state,created_at,updated_at,error FROM experiments ORDER BY created_at DESC,id DESC LIMIT 100')]
        return {'runs':rows,'active_id':self.manager.active}

    def start(self,payload):
        strategy=payload.get('strategy','single')
        enabled=payload.get('routing_enabled',strategy!='single')
        if type(enabled) is not bool:raise ValueError('路由开关应为开启或关闭')
        if not enabled:strategy='single'
        elif strategy=='single':raise ValueError('路由已开启，请选择组合方式')
        if strategy not in ('single','difficulty','cascade'):raise ValueError('测试策略无效')
        count=payload.get('count',10);threshold=payload.get('atom_threshold',25) if strategy=='difficulty' else 25
        if type(count) is not int or not 1<=count<=300 or type(threshold) is not int or not 5<=threshold<=100:raise ValueError('测试数量或路由阈值无效')
        execution=ai_parallel.validate({'parallelism':payload.get('parallelism',1)})
        with self.manager.lock:
            if self.manager.active:raise ValueError('请先暂停其他AI任务并等待返回')
            primary=self.profile_snapshot(payload.get('primary'),payload.get('overrides'))
            shared={k:v for k,v in payload.get('overrides',{}).items() if k!='model'}
            secondary=self.profile_snapshot(payload.get('secondary'),payload.get('secondary_overrides',shared)) if strategy!='single' else None
            # Check decryptability without changing the production selected credential.
            for profile in (primary,secondary):
                if profile:
                    key,_,error=self.manager.vault.read(profile['credential_id'])
                    if error or not key:raise ValueError(error or '端点密钥不可用')
            with closing(connect(self.path)) as c:
                samples=c.execute('SELECT * FROM samples WHERE set_id=? ORDER BY ordinal LIMIT ?',(payload.get('set_id'),count)).fetchall()
                if not samples:raise ValueError('请先冻结测试样本')
                identifier=uuid.uuid4().hex;config={**execution,'purpose':'evaluation_only','strategy':strategy,'routing_enabled':enabled,'primary':primary,'secondary':secondary,
                    'atom_threshold':threshold,'router_version':'benchmark-route-v1','sample_ordinals':[s['ordinal'] for s in samples],
                    'input_set_hash':fingerprint([s['input_hash'] for s in samples]),'reference_set_hash':fingerprint([s['reference_json'] for s in samples])}
                c.execute('INSERT INTO experiments VALUES(?,?,?,?,?,?,?,NULL)',(identifier,str(payload.get('name') or primary['name'])[:120],payload['set_id'],'running',db.dumps(config),db.now(),db.now()))
                for s in samples:c.execute("INSERT INTO cases(run_id,ordinal,molecule_id,formula,state) VALUES(?,?,?,?,'pending')",(identifier,s['ordinal'],s['molecule_id'],s['formula']))
                c.commit()
            self.manager.active=identifier;self.manager.cancel.clear()
        self.manager.pool.submit(self.work,identifier)
        return {'run_id':identifier,'selected':len(samples),'purpose':'evaluation_only'}

    def pause(self,identifier):
        with self.manager.lock:
            if self.manager.active!=identifier:raise ValueError('此测试没有运行')
            self.manager.cancel.set()
        with closing(connect(self.path)) as c:c.execute("UPDATE experiments SET state='stopping',updated_at=? WHERE id=? AND state='running'",(db.now(),identifier));c.commit()
        self.manager.cancel_requests()
        return {'state':'stopping'}

    def resume(self,identifier):
        with self.manager.lock:
            if self.manager.active:raise ValueError('已有运行任务')
            with closing(connect(self.path)) as c:
                row=c.execute('SELECT * FROM experiments WHERE id=?',(identifier,)).fetchone()
                if not row:raise KeyError(identifier)
                if not c.execute("SELECT 1 FROM cases WHERE run_id=? AND state IN ('pending','calling')",(identifier,)).fetchone():raise ValueError('没有待继续样本，请复制配置新建轮次')
                c.execute("UPDATE experiments SET state='running',error=NULL,updated_at=? WHERE id=?",(db.now(),identifier));c.commit()
            self.manager.active=identifier;self.manager.cancel.clear()
        self.manager.pool.submit(self.work,identifier);return {'run_id':identifier}

    def call(self,case,stage,profile,context):
        with closing(connect(self.path)) as c:old=c.execute('SELECT * FROM attempts WHERE case_id=? AND stage=?',(case['id'],stage)).fetchone()
        if old:
            if old['state']=='done':return json.loads(old['output_json'])
            if old['state']=='responded':
                output=ai.decode_response(json.loads(old['raw_json']),profile['config'],context)
                with closing(connect(self.path)) as c:c.execute("UPDATE attempts SET state='done',output_json=? WHERE id=?",(db.dumps(output),old['id']));c.commit()
                return output
            raise ValueError('此阶段请求已执行或结果未知，未自动重复调用')
        key,_,error=self.manager.vault.read(profile['credential_id'])
        if error or not key:raise ValueError(error or '密钥不可用')
        with closing(connect(self.path)) as c:
            cursor=c.execute("INSERT INTO attempts(case_id,stage,profile_json,input_json,state,created_at) VALUES(?,?,?,?,'calling',?)",(case['id'],stage,db.dumps(profile),db.dumps(context),db.now()));aid=cursor.lastrowid;c.commit()
        def archive(envelope):
            with closing(connect(self.path)) as c:c.execute("UPDATE attempts SET state='responded',raw_json=?,usage_json=?,diagnostics_json=? WHERE id=?",(db.dumps(envelope),db.dumps(envelope['usage']),db.dumps(envelope['diagnostics']),aid));c.commit()
        progress_state={}
        def progress(value):
            progress_state.update(value,updated_at=db.now())
            with closing(connect(self.path)) as c:c.execute('UPDATE attempts SET diagnostics_json=? WHERE id=?',(db.dumps({'progress':progress_state}),aid));c.commit()
        started=time.perf_counter()
        try:
            with self.manager.request_scope(),ai.capture_response(archive),ai.capture_progress(progress):output,usage=ai.complete(key,profile['config'],context)
            with closing(connect(self.path)) as c:c.execute("UPDATE attempts SET state='done',output_json=?,usage_json=?,diagnostics_json=?,elapsed=? WHERE id=?",(db.dumps(output),db.dumps(usage),db.dumps(output.get('api_diagnostics')),time.perf_counter()-started,aid));c.commit()
            return output
        except Exception as exc:
            message=str(exc).replace(key,'[REDACTED]')[:1500]
            with closing(connect(self.path)) as c:
                row=c.execute('SELECT diagnostics_json,usage_json FROM attempts WHERE id=?',(aid,)).fetchone()
                diagnostics={**json.loads(row['diagnostics_json'] or '{}'),**(getattr(exc,'diagnostics',None) or {})}
                if progress_state:diagnostics['progress']=progress_state
                usage=getattr(exc,'usage',None) or json.loads(row['usage_json'] or '{}')
                c.execute("UPDATE attempts SET state='error',error=?,elapsed=?,diagnostics_json=?,usage_json=? WHERE id=?",(message,time.perf_counter()-started,db.dumps(diagnostics),db.dumps(usage),aid));c.commit()
            if isinstance(exc,ai.APIRequestCancelled):raise
            if ai.should_pause_batch(exc):
                failure=ai.APIRequestError(message,status=getattr(exc,'status',None),uncertain=getattr(exc,'uncertain',True),usage=usage)
                failure.diagnostics.update(diagnostics)
                raise failure from None
            raise ValueError(message) from None

    def _case(self,identifier,case,cfg,context):
        begin=time.perf_counter();prior_elapsed=case['elapsed'] or 0
        try:
            hard=len(context['atoms'])>cfg['atom_threshold'] or sum(a['aromatic'] for a in context['atoms'])>12
            primary=cfg['secondary'] if cfg['strategy']=='difficulty' and hard else cfg['primary']
            output=self.call(case,0,primary,context)
            if cfg['strategy']=='cascade':
                uncertain=[{'spectrum_id':p['spectrum_id'],'peak_id':p['peak_id']} for p in output['assignments'] if p['status'] in ('ambiguous','unassigned')]
                if uncertain:
                    if self.manager.cancel.is_set():
                        with closing(connect(self.path)) as c:c.execute('UPDATE cases SET elapsed=? WHERE id=?',(prior_elapsed+time.perf_counter()-begin,case['id']));c.commit()
                        return None
                    secondary_context=copy.deepcopy(context);secondary_context['requested_peaks']=uncertain
                    other=self.call(case,1,cfg['secondary'],secondary_context)
                    keys={(p['spectrum_id'],p['peak_id']) for p in uncertain};replace={(p['spectrum_id'],p['peak_id']):p for p in other['assignments'] if (p['spectrum_id'],p['peak_id']) in keys}
                    output['assignments']=[replace.get((p['spectrum_id'],p['peak_id']),p) for p in output['assignments']]
                    output['summary']+='；仅歧义/未归属峰由第二阶段独立分析，参考答案未提供。'
            with closing(connect(self.path)) as c:c.execute("UPDATE cases SET state='done',output_json=?,elapsed=?,error=NULL WHERE id=?",(db.dumps(output),prior_elapsed+time.perf_counter()-begin,case['id']));c.commit()
            return None
        except Exception as exc:
            with closing(connect(self.path)) as c:c.execute("UPDATE cases SET state='error',error=?,elapsed=? WHERE id=?",(str(exc)[:1500],prior_elapsed+time.perf_counter()-begin,case['id']));c.commit()
            return exc

    def work(self,identifier):
        started=time.perf_counter();sid=None;error=None;halt=False;inflight={}
        with closing(connect(self.path)) as c:
            sid=c.execute('INSERT INTO sessions(run_id,started_at) VALUES(?,?)',(identifier,db.now())).lastrowid;c.commit()
            run=c.execute('SELECT * FROM experiments WHERE id=?',(identifier,)).fetchone();cfg=json.loads(run['config_json'])
        pool=ThreadPoolExecutor(max_workers=ai_parallel.validate(cfg)['parallelism'],thread_name_prefix='nmr-benchmark')
        try:
            while True:
                while not halt and not self.manager.cancel.is_set() and len(inflight)<ai_parallel.validate(cfg)['parallelism']:
                    with closing(connect(self.path)) as c:
                        c.execute('BEGIN IMMEDIATE')
                        case=c.execute("SELECT * FROM cases WHERE run_id=? AND state IN ('pending','calling') ORDER BY ordinal",(identifier,)).fetchall()
                        case=next((r for r in case if r['id'] not in inflight.values()),None)
                        if not case:break
                        sample=c.execute('SELECT * FROM samples WHERE set_id=? AND ordinal=?',(run['set_id'],case['ordinal'])).fetchone()
                        c.execute("UPDATE cases SET state='calling' WHERE id=?",(case['id'],));c.commit()
                    inflight[pool.submit(self._case,identifier,case,cfg,json.loads(sample['input_json']))]=case['id']
                if not inflight:break
                done,_=wait(inflight,timeout=.2,return_when=FIRST_COMPLETED)
                for future in done:
                    inflight.pop(future)
                    exc=future.result()
                    if exc is not None and not isinstance(exc,ai.APIRequestCancelled) and ai.should_pause_batch(exc):
                        halt=True
                        kind=(getattr(exc,'diagnostics',None) or {}).get('error_kind')
                        prefix=('本地等待已停止，批次已暂停；后续分子尚未请求。' if kind=='total_timeout' else
                                '模型生成异常，批次已暂停；后续分子尚未请求。' if kind=='repetitive_output' else
                                '接口异常，批次已自动暂停；后续分子尚未请求。')
                        error=prefix+str(exc)
                with closing(connect(self.path)) as c:c.execute('UPDATE experiments SET updated_at=? WHERE id=?',(db.now(),identifier));c.commit()
        except Exception as exc:error=str(exc)[:1500]
        finally:
            pool.shutdown(wait=True)
            with closing(connect(self.path)) as c:
                c.execute('UPDATE sessions SET ended_at=?,elapsed=? WHERE id=?',(db.now(),time.perf_counter()-started,sid))
                c.execute('UPDATE experiments SET state=?,error=?,updated_at=? WHERE id=?',('attention' if error else 'paused' if self.manager.cancel.is_set() else 'completed',error,db.now(),identifier));c.commit()
            with self.manager.lock:self.manager.active=None

    def detail(self,identifier):
        with closing(connect(self.path)) as c:
            row=c.execute('SELECT * FROM experiments WHERE id=?',(identifier,)).fetchone()
            if not row:raise KeyError(identifier)
            cases=[dict(r) for r in c.execute('SELECT c.*,s.reference_json,s.input_hash FROM cases c JOIN samples s ON s.set_id=? AND s.ordinal=c.ordinal WHERE c.run_id=? ORDER BY c.ordinal',(row['set_id'],identifier))]
            attempts=[dict(r) for r in c.execute('SELECT a.* FROM attempts a JOIN cases c ON a.case_id=c.id WHERE c.run_id=? ORDER BY a.id',(identifier,))]
            active_seconds=c.execute('SELECT sum(elapsed) FROM sessions WHERE run_id=?',(identifier,)).fetchone()[0]
        output=dict(row);output['config']=json.loads(output.pop('config_json'))
        output['metrics']=ai.summarize(cases);output['counts']=dict(Counter(r['state'] for r in cases))
        usage={};reported=0;reasoning=0;reasoning_known=0
        for a in attempts:
            u=json.loads(a['usage_json'] or '{}');reported+=isinstance(u.get('total_tokens'),int)
            reason=(u.get('completion_tokens_details') or {}).get('reasoning_tokens')
            if reason is not None:reasoning+=reason;reasoning_known+=1
            for k in ('prompt_tokens','completion_tokens','total_tokens'):usage[k]=usage.get(k,0)+(u.get(k,0) or 0)
        output['usage']=usage;output['usage_reported_attempts']=reported;output['api_attempts']=len(attempts)
        output['reasoning_tokens']=reasoning;output['reasoning_reported_attempts']=reasoning_known
        output['attempt_timing']=timings(a['elapsed'] for a in attempts)
        output['timing']=timings(r['elapsed'] for r in cases);output['active_wall_seconds']=active_seconds
        output['cases']=[{k:r[k] for k in ('id','ordinal','molecule_id','formula','state','elapsed','error','input_hash')} for r in cases]
        attempt_kinds={a['case_id']:((json.loads(a['diagnostics_json'] or '{}') or {}).get('error_kind')) for a in attempts}
        for case in output['cases']:case['error_kind']=attempt_kinds.get(case['id'])
        output['error_kind']=next((r['error_kind'] for r in reversed(output['cases']) if r['state']=='error'),None)
        current=next((a for a in reversed(attempts) if a['state']=='calling'),None)
        output['current_request']=None
        if current:
            case=next(r for r in cases if r['id']==current['case_id']);profile=json.loads(current['profile_json'])
            progress=json.loads(current['diagnostics_json'] or '{}').get('progress') or {}
            now=datetime.now(timezone.utc)
            elapsed=max(0,(now-datetime.fromisoformat(current['created_at'])).total_seconds())
            output['current_request']={'case_id':case['id'],'formula':case['formula'],'molecule_id':case['molecule_id'],
                                      'stage':current['stage']+1,'model':profile['config']['model'],'elapsed_s':elapsed,
                                      'started_at':current['created_at'],'stream':profile['config'].get('stream',False),
                                      'timeout_s':profile['config']['timeout'],'progress':progress}
        output['isolated']=True;return output

    def case(self,identifier):
        with closing(connect(self.path)) as c:
            row=c.execute('SELECT c.*,e.set_id FROM cases c JOIN experiments e ON e.id=c.run_id WHERE c.id=?',(int(identifier),)).fetchone()
            if not row:raise KeyError(identifier)
            sample=c.execute('SELECT * FROM samples WHERE set_id=? AND ordinal=?',(row['set_id'],row['ordinal'])).fetchone()
            attempts=[dict(r) for r in c.execute('SELECT * FROM attempts WHERE case_id=? ORDER BY stage',(row['id'],))]
        result=dict(row);result['output']=json.loads(result.pop('output_json') or 'null');result['input']=json.loads(sample['input_json']);result['reference']=json.loads(sample['reference_json'])
        result['consistency_issues']=ai.review_issues(result['output'],result['input']) if result['output'] else []
        result['comparisons']=ai.score(result['output'],result['reference'])
        for r in attempts:
            for k in ('profile','input','raw','output','usage','diagnostics'):r[k]=json.loads(r.pop(k+'_json') or 'null')
        result['attempts']=attempts;return result

    def compare(self,ids):
        if not 2<=len(ids)<=6 or len(set(ids))!=len(ids):raise ValueError('请选择2–6轮测试对比')
        runs=[self.detail(i) for i in ids]
        by_run={}
        with closing(connect(self.path)) as c:
            for run in runs:
                values=[dict(r) for r in c.execute('SELECT c.*,s.reference_json,s.input_hash FROM cases c JOIN samples s ON s.set_id=? AND s.ordinal=c.ordinal WHERE c.run_id=?',(run['set_id'],run['id']))]
                by_run[run['id']]={r['input_hash']+':'+fingerprint(json.loads(r['reference_json'])):r for r in values}
        common=set.intersection(*[set(v) for v in by_run.values()])
        if not common:raise ValueError('没有相同的冻结输入；请使用同一样本集')
        paired=[]
        with closing(connect(self.path)) as c:
            for run in runs:
                rows=[v for k,v in by_run[run['id']].items() if k in common];case_ids={r['id'] for r in rows}
                attempts=[dict(r) for r in c.execute('SELECT a.* FROM attempts a JOIN cases c ON c.id=a.case_id WHERE c.run_id=?',(run['id'],)) if r['case_id'] in case_ids]
                usage={k:sum((json.loads(a['usage_json'] or '{}').get(k,0) or 0) for a in attempts) for k in ('prompt_tokens','completion_tokens','total_tokens')}
                paired.append({'run_id':run['id'],'metrics':ai.summarize(rows),'timing':timings(r['elapsed'] for r in rows),'n':len(rows),
                               'usage':usage,'api_attempts':len(attempts),'usage_reported_attempts':sum(isinstance(json.loads(a['usage_json'] or '{}').get('total_tokens'),int) for a in attempts)})
        return {'runs':runs,'paired':paired,'same_inputs':all(r['config']['input_set_hash']==runs[0]['config']['input_set_hash'] for r in runs),
                'same_references':all(r['config'].get('reference_set_hash')==runs[0]['config'].get('reference_set_hash') for r in runs),
                'common_samples':len(common),'note':'失败计入准确率分母；不同样本或未完成轮次需看共同输入和完成状态。'}

    def export(self,ids):
        comparison=self.compare(ids)
        return {'schema':'evaluation-only-v1','comparison':comparison,
                'results':{r['id']:[self.case(c['id']) for c in r['cases']] for r in comparison['runs']}}
