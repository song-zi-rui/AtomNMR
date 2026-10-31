"""Compare archived requests without API calls or changes to annotations."""
import csv
import hashlib
import json
import statistics
from collections import Counter
from pathlib import Path
import sqlite3
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from annotation_tool_v3.ai_consistency import review_issues

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/assignment_database_v1/audits/api_comparison_20261030'
JOBS={'official':['14a491b7acbf486192ab21f71313dccf','3b0cdd7e57f24815860cd2d9a87b393f'],
      'group':['ea109a75a890480aa7723279dc166b74']}
HARD={'assigned','equivalent_group','overlap'}
def load(x,default=None):return json.loads(x) if x else default
def h(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
def stats(values):
    x=sorted(v for v in values if v is not None)
    return {'n':len(x),'mean':statistics.mean(x) if x else None,'median':statistics.median(x) if x else None,
            'p90':x[max(0,__import__('math').ceil(.9*len(x))-1)] if x else None}
def key(r):return r['spectrum_id'],r['peak_id']
def write_csv(path,rows):
    if not rows:return
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()})

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    c=sqlite3.connect((ROOT/'data/assignment_database_v1/databases/ai_annotations.sqlite').as_uri()+'?mode=ro',uri=True);c.row_factory=sqlite3.Row
    result={};records={};details=[]
    for side,ids in JOBS.items():
        values=[]
        for job_id in ids:
            job=c.execute('SELECT * FROM auto_jobs WHERE id=?',(job_id,)).fetchone()
            for row in c.execute('SELECT * FROM auto_items WHERE job_id=?',(job_id,)):
                cfg=load(row['request_config_json'],load(job['config_json']));ctx=load(row['input_json']);out=load(row['output_json'],{})
                diag=load(row['diagnostics_json'],{}) or {};usage=load(row['usage_json'],{}) or {}
                wanted={key(p) for p in ctx.get('requested_peaks',[])} if 'requested_peaks' in ctx else { (s['spectrum_id'],p['peak_id']) for s in ctx['spectra'] for p in s['peaks']}
                predictions={key(p):p for p in out.get('assignments',[]) if key(p) in wanted}
                chem={k:v for k,v in ctx.items() if k not in ('anchors','requested_peaks')}
                counts=Counter(p['status'] for p in predictions.values());hard=sum(counts[s] for s in HARD)
                nuclei={s['spectrum_id']:s['nucleus'] for s in ctx['spectra']}
                transport=diag.get('stream_transport') or {};tokens_details=usage.get('completion_tokens_details') or {}
                record={'side':side,'job_id':job_id,'item_id':row['id'],'molecule_id':row['molecule_id'],'formula':row['formula'],'state':row['state'],
                    'atoms':len(ctx['atoms']),'all_peaks':sum(len(s['peaks']) for s in ctx['spectra']),'targets':len(wanted),'anchors':len(ctx.get('anchors',[])),
                    'returned_targets':len(predictions),'hard_targets':hard,'coverage':hard/len(wanted) if wanted else None,
                    'statuses':dict(counts),'elapsed':row['elapsed'],'summary_chars':len(out.get('summary','')),
                    'mean_reason_chars':statistics.mean(len(p.get('brief_reason','')) for p in predictions.values()) if predictions else None,
                    'input_tokens':usage.get('prompt_tokens'),'completion_tokens':usage.get('completion_tokens'),'reasoning_tokens':tokens_details.get('reasoning_tokens'),
                    'ttft_ms':(diag.get('timing') or {}).get('ttft_ms'),'api_wall_ms':(diag.get('timing') or {}).get('api_wall_ms'),
                    'publication_ms':(diag.get('timing') or {}).get('publish_ms'),'reasoning_characters':transport.get('reasoning_characters'),
                    'consistency_flagged_peaks':len({(p['spectrum_id'],p['peak_id']) for p in review_issues(out,ctx) if (p['spectrum_id'],p['peak_id']) in wanted}),
                    'nucleus_counts':{n:{'targets':sum(nuclei[k[0]]==n for k in wanted),'hard':sum(nuclei[k[0]]==n and p['status'] in HARD for k,p in predictions.items())} for n in ['13C','1H']},
                    'final_tokens':usage.get('completion_tokens',0)-tokens_details.get('reasoning_tokens',0) if 'completion_tokens' in usage and 'reasoning_tokens' in tokens_details else None,
                    'dispatch_concurrency':(cfg.get('execution_at_dispatch') or {}).get('effective_parallelism'),
                    'applied':(load(row['publication_json'],{}) or {}).get('applied'),
                    'chem_hash':h(chem),'input_hash':h(ctx),'reference_count':len(load(row['reference_json'],[])),
                    'response_model':out.get('response_model'),'wire':cfg.get('wire_format'),'thinking':cfg.get('thinking'),
                    'max_tokens':cfg.get('max_tokens'),'timeout':cfg.get('timeout'),'effort':diag.get('reasoning_effort'),
                    'prompt_hash':cfg.get('prompt_sha256'),'error':row['error']}
                record['_pred']=predictions;record['_wanted']=wanted;record['_ctx']=ctx;record['_out']=out
                values.append(record);details.append({k:v for k,v in record.items() if not k.startswith('_')})
        records[side]=values;done=[v for v in values if v['state']=='done'];status=Counter()
        for v in done:status.update(v['statuses'])
        denominator=sum(v['targets'] for v in done)
        result[side]={'requests':len(values),'states':dict(Counter(v['state'] for v in values)),
            'successful_elapsed':stats(v['elapsed'] for v in done),'successful_targets':stats(v['targets'] for v in done),
            'successful_atoms':stats(v['atoms'] for v in done),'successful_all_peaks':stats(v['all_peaks'] for v in done),
            'status_counts':dict(status),'requested_peaks':denominator,'returned_peaks':sum(v['returned_targets'] for v in done),
            'hard_coverage':sum(v['hard_targets'] for v in done)/denominator if denominator else None,
            'all_targets_definite_molecules':sum(v['hard_targets']==v['targets'] for v in done),
            'no_definite_molecules':sum(v['hard_targets']==0 for v in done),
            'completion_tokens':stats(v['completion_tokens'] for v in done),'reasoning_tokens':stats(v['reasoning_tokens'] for v in done),
            'final_tokens':stats(v['final_tokens'] for v in done),
            'nuclei':{n:{k:sum(v['nucleus_counts'][n][k] for v in done) for k in ['targets','hard']} for n in ['13C','1H']},
            'summary_chars':stats(v['summary_chars'] for v in done),'mean_reason_chars':stats(v['mean_reason_chars'] for v in done),
            'api_wall_ms':stats(v['api_wall_ms'] for v in done),'publication_ms':stats(v['publication_ms'] for v in done),
            'configurations':dict(Counter((str((v['wire'],v['thinking'],v['max_tokens'],v['timeout'],v['effort']))) for v in values)),
            'response_models':dict(Counter(str(v['response_model']) for v in done)),
            'flags':sum(v['consistency_flagged_peaks'] for v in done),'applied':sum(v['applied'] or 0 for v in done),
            'dispatch_concurrency_counts':dict(Counter(str(v['dispatch_concurrency']) for v in done)),
            'existing_ai_labels_in_input':sum(any(p.get('assignment') for s in v['_ctx']['spectra'] for p in s['peaks']) for v in values)}
    pairs=[];peaks=[]
    off={r['molecule_id']:r for r in records['official'] if r['state']=='done'}
    for b in records['group']:
        if b['state']!='done' or b['molecule_id'] not in off:continue
        a=off[b['molecule_id']];common=a['_wanted']&b['_wanted'];same=a['chem_hash']==b['chem_hash']
        counts=Counter()
        for k in common:
            p=a['_pred'].get(k,{});q=b['_pred'].get(k,{})
            ah=p.get('status') in HARD;bh=q.get('status') in HARD
            counts['official_definite']+=ah;counts['group_definite']+=bh
            if ah and bh:counts['both_definite']+=1;counts['same_atom_set']+=set(p['atom_indices'])==set(q['atom_indices'])
            peaks.append({'molecule_id':a['molecule_id'],'formula':a['formula'],'spectrum_id':k[0],'peak_id':k[1],
                'same_chemical_input':same,'official_status':p.get('status'),'official_atoms':p.get('atom_indices'),
                'group_status':q.get('status'),'group_atoms':q.get('atom_indices'),
                'official_reason':p.get('brief_reason'),'group_reason':q.get('brief_reason')})
        pairs.append({'molecule_id':a['molecule_id'],'formula':a['formula'],'same_chemical_input':same,'exact_input':a['input_hash']==b['input_hash'],
            'official_targets':a['targets'],'group_targets':b['targets'],'common_targets':len(common),
            'official_seconds':a['elapsed'],'group_seconds':b['elapsed'],'official_hard_coverage':a['coverage'],'group_hard_coverage':b['coverage'],**counts})
    result['paired']={'common_molecules':len(pairs),'same_chemical_input':sum(p['same_chemical_input'] for p in pairs),
        'exact_same_input':sum(p['exact_input'] for p in pairs),'same_input_pairs': [p for p in pairs if p['exact_input']],
        'same_chemical_common_peak_counts':{k:sum(p.get(k,0) for p in pairs if p['same_chemical_input']) for k in ['common_targets','official_definite','group_definite','both_definite','same_atom_set']}}
    write_csv(OUT/'requests.csv',details);write_csv(OUT/'paired_molecules.csv',pairs);write_csv(OUT/'paired_peaks.csv',peaks)
    (OUT/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False))

if __name__=='__main__':main()
