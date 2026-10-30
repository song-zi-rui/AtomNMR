"""Compare saved pre/post optimization requests without API calls or key access."""
import argparse
from collections import Counter,defaultdict
from contextlib import closing
from datetime import datetime,timedelta,timezone
import json
from pathlib import Path
import sqlite3
import statistics

ROOT=Path(__file__).resolve().parents[1]


def load(path):
    with closing(sqlite3.connect(path.resolve().as_uri()+'?mode=ro',uri=True)) as c:
        c.row_factory=sqlite3.Row;c.execute('PRAGMA query_only=ON');c.execute('BEGIN')
        rows=[dict(r) for r in c.execute('SELECT i.*,j.config_json FROM auto_items i JOIN auto_jobs j ON j.id=i.job_id ORDER BY i.id')]
    for r in rows:
        r['config']=json.loads(r['request_config_json'] or r['config_json'])
        for key in ('input','output','usage','diagnostics','reference'):
            r[key]=json.loads(r[key+'_json'] or ('[]' if key=='reference' else '{}')) or {}
    return rows


def tokens(usage):
    if not isinstance(usage.get('total_tokens'),int):return None
    reason=(usage.get('completion_tokens_details') or {}).get('reasoning_tokens')
    return {'input':usage.get('prompt_tokens'), 'reasoning':reason,
            'final':usage['completion_tokens']-reason if reason is not None else None,
            'total':usage['total_tokens']}


def summarize(rows):
    reported=[r for r in rows if tokens(r['usage']) is not None]
    total={k:sum(tokens(r['usage'])[k] for r in reported if tokens(r['usage'])[k] is not None) for k in ('input','reasoning','final','total')}
    values={k:[tokens(r['usage'])[k] for r in reported if tokens(r['usage'])[k] is not None] for k in total}
    mean={k:statistics.mean(v) if v else None for k,v in values.items()}
    states=Counter(r['state'] for r in rows)
    statuses=Counter(p['status'] for r in rows if r['output'] for p in r['output']['assignments'])
    elapsed=[r['elapsed'] for r in reported if r['elapsed'] is not None]
    return {'requests':len(rows),'reported_usage':len(reported),'states':dict(states),'tokens':total,'mean_tokens':mean,
            'mean_elapsed_seconds':statistics.mean(elapsed) if elapsed else None,
            'mean_atoms':statistics.mean(len(r['input']['atoms']) for r in rows) if rows else None,
            'mean_peaks':statistics.mean(sum(len(s['peaks']) for s in r['input']['spectra']) for r in rows) if rows else None,
            'reference_labels':sum(len(r['reference']) for r in rows),'assignment_states':dict(statuses)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--db',type=Path,default=ROOT/'data/assignment_database_v1/ai_annotations.sqlite')
    p.add_argument('--output',type=Path,default=ROOT/'docs/ai_token_audit_20261028/optimization_comparison.json')
    args=p.parse_args();rows=load(args.db)
    old=[r for r in rows if r['config'].get('wire_format','verbose-v1')=='verbose-v1']
    new=[r for r in rows if r['config'].get('wire_format')=='compact-v2']
    fresh=[r for r in new if not r['config'].get('retry_of')]
    retried=[r for r in new if r['config'].get('retry_of') and not r['config'].get('local_response_reuse')]
    previous={}
    for r in old:previous.setdefault(r['molecule_id'],r)
    pairs=[]
    for r in retried:
        before=previous.get(r['molecule_id'])
        if not before:continue
        pairs.append({'molecule_id':r['molecule_id'],'formula':r['formula'],'old_id':before['id'],'new_id':r['id'],
                      'identical_context':before['input']==r['input'],'old_state':before['state'],'new_state':r['state'],
                      'old_finish':before['diagnostics'].get('finish_reason'),'new_finish':r['diagnostics'].get('finish_reason'),
                      'old_cap':before['config']['max_tokens'],'new_cap':r['config']['max_tokens'],
                      'old':tokens(before['usage']),'new':tokens(r['usage'])})
    known=[r for r in pairs if r['identical_context'] and r['old'] and r['new']]
    paired={k:{'old':sum(r['old'][k] for r in known),'new':sum(r['new'][k] for r in known)} for k in ('input','reasoning','final','total')}
    for v in paired.values():v['change_percent']=(v['new']/v['old']-1)*100
    actual=defaultdict(list)
    for r in fresh:actual[r['diagnostics'].get('reasoning_effort','unreported')].append(r)
    result={'generated_at':datetime.now(timezone(timedelta(hours=8))).isoformat(),
            'old':summarize(old),'new_retry':summarize(retried),'new_fresh':summarize(fresh),
            'actual_effort':{k:summarize(v) for k,v in actual.items()},'same_input_pairs':pairs,
            'paired_reported_count':len(known),'paired_tokens':paired,
            'limits':['Retry pairs compare old failures with new completed responses; caps differ, so they are not a controlled total-cost experiment.',
                      'Fresh old/new requests concern different molecules; group averages do not establish causal savings.',
                      'No API calls or key access. Unknown usage is excluded, never treated as zero.',
                      'Successful database publication is not a measurement of assignment accuracy.']}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:result[k] for k in ('paired_reported_count','paired_tokens','old','new_fresh')},ensure_ascii=True,indent=2))


if __name__=='__main__':main()
