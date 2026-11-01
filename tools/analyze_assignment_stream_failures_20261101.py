"""Read-only task diagnosis and offline recovery previews; no API or publication."""
from __future__ import annotations

import copy
from collections import Counter
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import sqlite3
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from annotation_tool_v3 import ai_annotation as ai
from annotation_tool_v3.ai_json_recovery import close_missing_root_object

CONFIG_KEYS=('model','api_adapter','base_url','thinking','reasoning_effort','stream','wire_format',
             'json_mode','max_tokens','timeout','parallelism','provider_parallelism_cap')
JOBS=('3b0cdd7e57f24815860cd2d9a87b393f','ee1e4f1c6b2d49ef8ba0fd977be33c7b')


def main():
    source=ROOT/'data/assignment_database_v1/databases/ai_annotations.sqlite'
    output=ROOT/'outputs/stream_failure_analysis_20261101'
    output.mkdir(exist_ok=True)
    findings={'captured_at_utc':datetime.now(timezone.utc).isoformat(),'database':str(source),
        'mode':'read-only SQLite transaction; offline decoder only; no API calls or database updates',
        'jobs':[],'json_previews':[]}
    c=sqlite3.connect(source.resolve().as_uri()+'?mode=ro',uri=True,timeout=5)
    c.row_factory=sqlite3.Row;c.execute('BEGIN')
    try:
        for jid in JOBS:
            job=dict(c.execute('SELECT id,name,state,processed,success,failures,updated_at FROM auto_jobs WHERE id=?',(jid,)).fetchone())
            rows=c.execute('SELECT id,state,error,elapsed,created_at,updated_at,diagnostics_json,request_config_json FROM auto_items WHERE job_id=? ORDER BY id',(jid,)).fetchall()
            cohorts=Counter();errors=[]
            for r in rows:
                config=json.loads(r['request_config_json'] or '{}');diag=json.loads(r['diagnostics_json'] or '{}')
                effective=(config.get('execution_at_dispatch') or {}).get('effective_parallelism')
                cohorts[(r['state'],config.get('model'),config.get('stream'),config.get('parallelism'),effective)]+=1
                if r['state']=='error':
                    errors.append({'id':r['id'],'error':r['error'],'elapsed_s':r['elapsed'],
                        'created_at_utc':r['created_at'],'updated_at_utc':r['updated_at'],'diagnostics':diag,
                        'request_config':{k:config.get(k) for k in CONFIG_KEYS},'effective_parallelism':effective})
            job['cohorts']=[dict(zip(('state','model','stream','configured_parallelism','effective_parallelism'),key),count=n) for key,n in cohorts.items()]
            job['errors']=errors;findings['jobs'].append(job)
        for iid in (2761,2634):
            row=c.execute('SELECT id,molecule_id,formula,input_json,request_config_json,raw_response_json FROM auto_items WHERE id=?',(iid,)).fetchone()
            config=json.loads(row['request_config_json']);context=json.loads(row['input_json']);envelope=json.loads(row['raw_response_json'])
            original=envelope['content'];candidate=close_missing_root_object(original)
            assert candidate==original+'}', 'Only one root brace may be appended to these observed examples'
            decoded=ai.decode_response(copy.deepcopy(envelope),config,context,recover=True)
            assert len(decoded['assignments'])==len(context['requested_peaks']) and decoded['missing_rows']==0
            assert envelope['content']==original
            preview={'item_id':iid,'molecule_id':row['molecule_id'],'formula':row['formula'],
                'applied_to_database':False,'remote_api_calls':0,'original_content':original,
                'original_sha256':hashlib.sha256(original.encode()).hexdigest(),'proposed_append':'}',
                'finish_reason':envelope['diagnostics']['finish_reason'],
                'done_marker_received':envelope['diagnostics']['stream_transport']['done_marker_received'],
                'requested_peaks':context['requested_peaks'],'decoded_preview':decoded,
                'note':'Syntax and requested row/atom/nucleus checks passed; this does not establish chemical assignment correctness or authorize publication.'}
            (output/f'item_{iid}_local_recovery_preview.json').write_text(json.dumps(preview,ensure_ascii=False,indent=2),encoding='utf-8')
            findings['json_previews'].append({k:preview[k] for k in ('item_id','molecule_id','formula','applied_to_database','remote_api_calls','proposed_append','finish_reason','done_marker_received')})
    finally:c.close()
    (output/'findings.json').write_text(json.dumps(findings,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'findings':str(output/'findings.json'),'jobs':[{k:j[k] for k in ('id','state','processed','success','failures','cohorts')} for j in findings['jobs']],
                      'local_recovery_previews':findings['json_previews']},ensure_ascii=True))


if __name__=='__main__':main()
