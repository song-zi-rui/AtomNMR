"""One explicitly requested, read-only-database API diagnostic. No DB publication.

Only final answer content is retained; private reasoning text is counted and
immediately discarded. No automatic retries or app configuration mutations.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import socket
import sqlite3
import sys
import threading
import time
from urllib.parse import urlsplit
from urllib.request import Request, build_opener, HTTPSHandler, ProxyHandler
from contextlib import closing

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from annotation_tool_v3 import ai_annotation as ai, api_tls, local_secrets


def save(path,value):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8');temp.replace(path)


def choose(c):
    # Includes both continuous production and the legacy annotation runner.
    candidates=[]
    for table,parent,fk in [('auto_items','auto_jobs','job_id'),('items','runs','run_id')]:
        for r in c.execute(f"SELECT a.*,j.config_json FROM {table} a JOIN {parent} j ON j.id=a.{fk} WHERE a.state='done' ORDER BY a.elapsed DESC"):
            cfg=json.loads((r['request_config_json'] if 'request_config_json' in r.keys() else None) or r['config_json'])
            if urlsplit(cfg.get('base_url','')).hostname=='api.deepseek.com':
                candidates.append((r['elapsed'],table,dict(r),cfg));break
    if not candidates:raise ValueError('No completed official annotation found')
    return max(candidates,key=lambda x:x[0])


def compare(context,official,group):
    rows=[]
    old={(a['spectrum_id'],a['peak_id']):a for a in official.get('assignments',[])}
    new={(a['spectrum_id'],a['peak_id']):a for a in group.get('assignments',[])}
    for s in context['spectra']:
        for peak in s['peaks']:
            key=(s['spectrum_id'],peak['peak_id']);a=old.get(key,{});b=new.get(key,{})
            hard=lambda x:x.get('status') in ai.HARD
            match=hard(a) and hard(b) and set(a.get('atom_indices',[]))==set(b.get('atom_indices',[]))
            rows.append({'nucleus':s.get('nucleus'),'spectrum_id':key[0],'peak_id':key[1],
                'shift_ppm':peak.get('shift_ppm'),'official_status':a.get('status'),
                'official_atoms':a.get('atom_indices',[]),'official_candidates':a.get('candidates',[]),
                'official_reason':a.get('brief_reason',''),'group_status':b.get('status'),
                'group_atoms':b.get('atom_indices',[]),'group_candidates':b.get('candidates',[]),
                'group_reason':b.get('brief_reason',''),'both_hard':hard(a) and hard(b),'same_hard_atom_set':match})
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--max-seconds',type=int,default=2400)
    parser.add_argument('--source-item',type=int)
    parser.add_argument('--thinking',choices=['original','off'],default='original')
    parser.add_argument('--max-tokens',type=int)
    args=parser.parse_args();out=args.output.resolve()
    if out.exists():raise ValueError('Output directory already exists; refusing duplicate request')
    out.mkdir(parents=True);(out/'official').mkdir();(out/'group').mkdir()
    source=ROOT/'data/assignment_database_v1/databases/ai_annotations.sqlite'
    with closing(sqlite3.connect(source.resolve().as_uri()+'?mode=ro',uri=True)) as c:
        c.row_factory=sqlite3.Row
        if args.source_item is None:
            elapsed,table,item,cfg=choose(c)
        else:
            record=c.execute("SELECT a.*,j.config_json FROM auto_items a JOIN auto_jobs j ON j.id=a.job_id WHERE a.id=? AND a.state='done'",(args.source_item,)).fetchone()
            if not record:raise ValueError('Completed source item not found')
            item=dict(record);table='auto_items';elapsed=item['elapsed']
            cfg=json.loads(item.get('request_config_json') or item['config_json'])
            if urlsplit(cfg.get('base_url','')).hostname!='api.deepseek.com':raise ValueError('Source is not an official result')
        group_job=c.execute("SELECT config_json FROM auto_jobs WHERE config_json LIKE '%10.8.130.88%' ORDER BY created_at DESC LIMIT 1").fetchone()
        endpoint=json.loads(group_job[0])
        credential=c.execute("SELECT encrypted FROM api_credentials WHERE name='deepseek-DigitalChemistry' AND enabled=1").fetchone()
        if not credential:raise ValueError('Enabled group credential not found')
        key=local_secrets.unprotect(credential[0])
    context=json.loads(item['input_json']);official=json.loads(item['output_json'])
    official_diag=json.loads(item['diagnostics_json'] or '{}')
    version=cfg['wire_format'];wire,mapping=ai.protocol.encode(context,version=version)
    prompt=cfg.get('system_prompt',ai.protocol.system_prompt(version))
    thinking=bool(cfg['thinking']) if args.thinking=='original' else False
    token_limit=args.max_tokens or cfg['max_tokens']
    actual=(official_diag.get('reasoning_effort') or ai.protocol.effort(cfg,context)) if thinking else 'none'
    payload={'model':endpoint['model'],'messages':[{'role':'system','content':prompt},{'role':'user','content':ai.db.dumps(wire)}],
        'max_tokens':token_limit,'stream':True,'stream_options':{'include_usage':True},
        'chat_template_kwargs':{'thinking':thinking}}
    if thinking:payload['chat_template_kwargs']['reasoning_effort']=25 if actual=='low' else 50
    if cfg.get('json_mode',True):payload['response_format']={'type':'json_object'}
    meta={'source_table':table,'source_item':item['id'],'source_job':item.get('job_id') or item.get('run_id'),
        'molecule_id':item['molecule_id'],'formula':item['formula'],'official_elapsed_s':elapsed,
        'official_diagnostics':official_diag,'started_at':datetime.now(timezone.utc).isoformat(),
        'max_seconds':args.max_seconds,'database_mode':'read-only; no writes or result publication',
        'explicit_overrides':{'thinking':args.thinking,'max_tokens':args.max_tokens},
        'controlled_differences':['group endpoint/model identifier','native thinking parameter adapter','streaming instrumentation','total wait budget'],
        'reasoning_retention':'counts and timing only; no private reasoning body'}
    save(out/'selection.json',meta);save(out/'input.json',context)
    save(out/'official/config.json',cfg);save(out/'official/assignments.json',official)
    save(out/'official/diagnostics.json',official_diag)
    if item.get('raw_response_json'):save(out/'official/response.json',json.loads(item['raw_response_json']))
    save(out/'group/request.json',payload)
    group_cfg={**cfg,**{k:endpoint[k] for k in ['base_url','model','api_adapter','tls_policy','ca_file','ca_sha256','tls_origin']},
               'stream':True,'timeout':args.max_seconds,'reasoning_effort':actual,'thinking':thinking,'max_tokens':token_limit}
    save(out/'group/config.json',group_cfg)
    status={'state':'connecting','elapsed_s':0,'events':0,'reasoning_characters':0,'answer_characters':0,
        'ttft_s':None,'first_reasoning_s':None,'first_answer_s':None,'first_finish_s':None,
        'finish_reason':None,'done_marker':False,'answer_json_complete_s':None,'usage':{},'model':None,
        'last_sse_event_s':None,'last_text_s':None,'local_deadline_trigger_s':None}
    save(out/'group/progress.json',status)
    data=ai.db.dumps(payload).encode();context_tls=api_tls.context_for(endpoint['base_url'],{k:endpoint.get(k,v) for k,v in api_tls.TLS_DEFAULTS.items()})
    opener=build_opener(ai.NoRedirect,ProxyHandler({}),HTTPSHandler(context=context_tls))
    req=Request(endpoint['base_url'].rstrip('/')+'/chat/completions',data=data,
                headers={'Authorization':'Bearer '+key,'Content-Type':'application/json','Accept':'text/event-stream'},method='POST')
    started=time.perf_counter();answer=[];last_save=0;pending=[];event_rows=[];response=None;ended=False
    full_response={};failure=None;timed_out=threading.Event()
    def deadline():
        status['local_deadline_trigger_s']=time.perf_counter()-started
        timed_out.set()
        if response is not None:
            try:response.fp.raw._sock.shutdown(socket.SHUT_RDWR)
            except (OSError,AttributeError):pass
    watchdog=threading.Timer(args.max_seconds,deadline);watchdog.daemon=True
    def checkpoint(force=False):
        nonlocal last_save
        now=time.perf_counter();status['elapsed_s']=now-started
        if force or now-last_save>=5:
            last_save=now
            # Responses from this endpoint use a separate reasoning channel.
            # Never persist an unexpected inline <think> body.
            final=''.join(answer)
            if '<think>' in final:
                end=final.rfind('</think>');final=final[end+8:] if end>=0 else ''
            (out/'group/answer.txt').write_text(final,encoding='utf-8')
            save(out/'group/progress.json',status)
    def event():
        nonlocal ended,full_response
        if not pending:return
        raw='\n'.join(pending);pending.clear();t=time.perf_counter()-started
        if raw.strip()=='[DONE]':status['done_marker']=True;ended=True;event_rows.append({'elapsed_s':t,'type':'done'});return
        value=json.loads(raw);status['events']+=1;status['last_sse_event_s']=t
        if value.get('error'):raise ValueError('Remote stream error; response body omitted')
        if value.get('model'):status['model']=value['model']
        if value.get('id'):status['response_id']=value['id']
        if isinstance(value.get('usage'),dict):status['usage'].update(value['usage'])
        added_answer=0;added_reasoning=0
        for choice in value.get('choices',[]):
            if choice.get('index',0)!=0:continue
            delta=choice.get('delta') or {};a=delta.get('content');r=delta.get('reasoning_content') or delta.get('reasoning')
            if isinstance(r,str):
                added_reasoning+=len(r)
                if r and status['first_reasoning_s'] is None:status['first_reasoning_s']=t
            # Count-only private reasoning. Neither r nor raw events are saved.
            if isinstance(a,str):
                answer.append(a);added_answer+=len(a)
                if a and status['first_answer_s'] is None:status['first_answer_s']=t
            if choice.get('finish_reason') is not None:
                status['finish_reason']=choice['finish_reason']
                if status['first_finish_s'] is None:status['first_finish_s']=t
            if added_answer or added_reasoning:
                if status['ttft_s'] is None:status['ttft_s']=t
        status['reasoning_characters']+=added_reasoning;status['answer_characters']+=added_answer
        if added_reasoning or added_answer:status['last_text_s']=t
        status['state']='answering' if status['answer_characters'] else 'reasoning' if status['reasoning_characters'] else 'waiting'
        if added_answer and status['answer_json_complete_s'] is None:
            try:
                parsed=json.loads(''.join(answer))
                if isinstance(parsed,dict):status['answer_json_complete_s']=t
            except ValueError:pass
        event_rows.append({'elapsed_s':round(t,6),'answer_added':added_answer,'reasoning_added':added_reasoning,
            'finish_reason':status['finish_reason'],'usage_received':bool(value.get('usage'))})
        checkpoint()
    watchdog.start()
    try:
        response=opener.open(req,timeout=args.max_seconds)
        status['http_status']=response.status;status['headers_s']=time.perf_counter()-started
        status['content_type']=response.headers.get('Content-Type')
        if 'text/event-stream' not in status['content_type'].lower():raise ValueError('Endpoint did not return SSE; request not retried')
        for line in iter(lambda:response.readline(1024*1024+1),b''):
            if timed_out.is_set():raise TimeoutError('standalone total deadline reached')
            if len(line)>1024*1024:raise ValueError('Stream line too large')
            s=line.decode('utf-8').rstrip('\r\n')
            if s.startswith('data:'):pending.append(s[5:].lstrip())
            elif not s:event()
            if ended:break
            checkpoint()
        event();status['state']='complete' if status['done_marker'] and status['finish_reason']=='stop' else 'stream_incomplete'
    except Exception as exc:
        # Never store credential-bearing exception bodies or stream contents.
        failure={'type':type(exc).__name__,'message':str(exc).replace(key,'[REDACTED]')[:600]}
        status['state']='timeout' if timed_out.is_set() or isinstance(exc,TimeoutError) else 'error'
        if timed_out.is_set():
            status['error']={'type':'LocalTotalTimeout','message':'达到本地总时限，已主动停止接收；远端完成情况未知。'}
            status['secondary_shutdown_exception']=failure
        else:status['error']=failure
    finally:
        watchdog.cancel()
        if response is not None:response.close()
        checkpoint(True)
        with (out/'group/events.jsonl').open('w',encoding='utf-8') as f:
            for row in event_rows:f.write(json.dumps(row)+'\n')
    final=(out/'group/answer.txt').read_text('utf-8');parsed_group={};decode_error=None
    envelope={'content':final,'usage':status['usage'],'diagnostics':{'finish_reason':status['finish_reason'],
        'response_model':status['model'],'wire_format':version,'reasoning_effort':actual,'stream_transport':dict(status)}}
    save(out/'group/response.json',envelope)
    if status['finish_reason']=='stop':
        try:
            parsed_group=ai.decode_response(envelope,group_cfg,context)
            save(out/'group/assignments.json',parsed_group)
        except Exception as exc:decode_error={'type':type(exc).__name__,'message':str(exc)[:800]}
    rows=compare(context,official,parsed_group)
    with (out/'peak_comparison.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else []);writer.writeheader();writer.writerows(rows)
    hard=lambda result:sum(a.get('status') in ai.HARD for a in result.get('assignments',[]))
    summary={'selection':meta,'group':status,'group_decode_error':decode_error,
        'official_assignment_count':len(official.get('assignments',[])),'official_hard_peaks':hard(official),
        'group_assignment_count':len(parsed_group.get('assignments',[])),'group_hard_peaks':hard(parsed_group),
        'both_hard_peaks':sum(r['both_hard'] for r in rows),'same_hard_atom_set':sum(r['same_hard_atom_set'] for r in rows),
        'interpretation':'Agreement with the official model is not chemical accuracy; one molecule is not a model benchmark.',
        'database_written':False,'api_calls':1,'reasoning_body_saved':False}
    save(out/'summary.json',summary)
    print(json.dumps({'state':status['state'],'elapsed_s':status['elapsed_s'],'summary':str(out/'summary.json')},ensure_ascii=True),flush=True)


if __name__=='__main__':main()
