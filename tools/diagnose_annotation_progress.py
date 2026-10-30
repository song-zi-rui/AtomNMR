"""Read-only diagnosis of saved labels, queues and repair snapshots."""
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone, timedelta
import json
from pathlib import Path
import sqlite3

ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / 'data/assignment_database_v1/assignments.sqlite'
OUT = ROOT / 'outputs/annotation_progress_diagnosis_20261012.json'


def connect(path):
    c = sqlite3.connect(path.resolve().as_uri() + '?mode=ro', uri=True, timeout=10)
    c.row_factory = sqlite3.Row
    c.execute('PRAGMA query_only=ON')
    return c


def rows(c, sql, args=()):
    return [dict(r) for r in c.execute(sql, args)]


def day(stamp):
    return datetime.fromisoformat(stamp).astimezone(timezone(timedelta(hours=8))).date().isoformat()


def main():
    result = {'database': str(DB), 'read_only': True, 'timezone': 'Asia/Shanghai'}
    with connect(DB) as c:
        result['databases'] = rows(c, 'PRAGMA database_list')
        result['queues'] = rows(c, 'SELECT queue,task_done,count(*) n FROM molecules GROUP BY queue,task_done')
        result['annotation_origins'] = rows(c, 'SELECT origin,status,count(*) n,min(updated_at) oldest,max(updated_at) newest FROM annotations GROUP BY origin,status')
        manual = rows(c, '''SELECT a.*,s.molecule_id,s.nucleus,s.source,m.formula,m.queue,m.task_done,
                           m.c_count,m.h_count,m.c_assigned,m.h_assigned,m.revision
                           FROM annotations a JOIN spectra s ON s.id=a.spectrum_id
                           JOIN molecules m ON m.id=s.molecule_id WHERE a.origin!='nmrshiftdb2' ''')
        result['non_source_annotations'] = manual
        by_day = defaultdict(lambda: {'labels': 0, 'molecules': set(), 'spectra': set(), 'status': Counter()})
        for r in manual:
            d = by_day[day(r['updated_at'])]
            d['labels'] += 1; d['molecules'].add(r['molecule_id']); d['spectra'].add(r['spectrum_id']); d['status'][r['status']] += 1
        result['non_source_by_local_day'] = {k: {**v,'molecules':len(v['molecules']),'spectra':len(v['spectra'])} for k,v in sorted(by_day.items())}
        result['manual_molecule_queues'] = dict(Counter({mid: next(r['queue'] for r in manual if r['molecule_id']==mid) for mid in {r['molecule_id'] for r in manual}}.values()))
        activity = rows(c, '''SELECT a.*,m.formula,m.queue,m.task_done,m.c_count,m.h_count,m.c_assigned,m.h_assigned
                              FROM molecule_activity a JOIN molecules m ON m.id=a.molecule_id ORDER BY last_processed_at''')
        result['activity'] = activity
        result['activity_by_local_day'] = dict(Counter(day(r['last_processed_at']) for r in activity))
        result['activity_actions'] = dict(Counter(r['last_action'] for r in activity))
        result['remaining'] = rows(c, '''SELECT m.id,m.formula,m.queue,m.task_done,m.c_count,m.h_count,m.c_assigned,m.h_assigned,m.revision,
                                       a.last_processed_at,a.last_action FROM molecules m LEFT JOIN molecule_activity a ON a.molecule_id=m.id
                                       WHERE m.queue='complete_remaining' ORDER BY m.id''')
        result['numbering_reviews'] = rows(c, 'SELECT id,molecule_id,spectrum_id,decision,note,created_at FROM spectrum_numbering_reviews')
        current_annotations = {(r['spectrum_id'],r['peak_id']):r for r in manual}
        remaining_ids = {r['id'] for r in result['remaining']}
        result['remaining_audits'] = []
        for mid in remaining_ids:
            result['remaining_audits'] += rows(c, '''SELECT s.id,s.molecule_id,s.nucleus,a.validation_json
                FROM spectra s LEFT JOIN spectrum_assignment_audit a ON a.spectrum_id=s.id WHERE s.molecule_id=?''',(mid,))
    print(json.dumps({k:v for k,v in result.items() if k not in ('non_source_annotations','activity','remaining','remaining_audits')}, ensure_ascii=False), flush=True)
    archives = []
    earliest_molecule = {}
    for path in sorted((ROOT/'data/shiftdb_atom_identity_audit').glob('*/repair_rollback_*.sqlite'), key=lambda p:p.name):
        detail = {'path':str(path), 'saved_manual_labels':0, 'missing_now':[], 'changed_now':[], 'manual_molecules':set(), 'queue_changes':Counter()}
        with connect(path) as b:
            for r in b.execute('SELECT spectrum_id,molecule_id,annotations_json,molecule_json FROM snapshots'):
                before = json.loads(r['molecule_json'])
                earliest_molecule.setdefault(r['molecule_id'],before)
                for a in json.loads(r['annotations_json']):
                    if a['origin']=='nmrshiftdb2': continue
                    detail['saved_manual_labels'] += 1
                    detail['manual_molecules'].add(r['molecule_id'])
                    current = current_annotations.get((a['spectrum_id'],a['peak_id']))
                    if current is None: detail['missing_now'].append(a)
                    elif any(current[k]!=a[k] for k in ('status','atoms_json','note','origin')):
                        detail['changed_now'].append({'before':a,'after':current})
        detail['manual_molecules'] = len(detail['manual_molecules'])
        archives.append(detail)
        print(json.dumps({'archive':path.name,'saved_manual_labels':detail['saved_manual_labels'],'missing_now':len(detail['missing_now']),'changed_now':len(detail['changed_now'])}, ensure_ascii=False), flush=True)
    result['repair_archives'] = archives
    result['remaining_before_repair'] = {mid: earliest_molecule[mid] for mid in remaining_ids if mid in earliest_molecule}
    result['remaining_before_queue_counts'] = dict(Counter(v['queue'] for v in result['remaining_before_repair'].values()))
    reopened = {mid for mid,v in result['remaining_before_repair'].items() if v['queue']=='complete'}
    result['reopened_molecules'] = sorted(reopened)
    result['reopened_with_manual_labels'] = sorted(reopened & {r['molecule_id'] for r in manual})
    result['reopened_audit_statuses'] = {}
    for mid in sorted(reopened):
        statuses = {json.loads(r['validation_json'])['status'] for r in result['remaining_audits']
                    if r['molecule_id']==mid and r['validation_json']}
        result['reopened_audit_statuses'][mid] = sorted(statuses-{'identity_validated'})
    result['reopened_issue_combinations'] = dict(Counter('|'.join(v) for v in result['reopened_audit_statuses'].values()))
    with connect(DB) as c:
        missing_peak_ids=[]
        for sid in {r['spectrum_id'] for r in manual}:
            p=c.execute('SELECT peaks_json FROM spectra WHERE id=?',(sid,)).fetchone()
            ids={v['peak_id'] for v in json.loads(p[0])}
            missing_peak_ids += [r for r in manual if r['spectrum_id']==sid and r['peak_id'] not in ids]
        result['manual_labels_missing_from_current_peak_table']=missing_peak_ids
        transitions=[]
        for mid,before in earliest_molecule.items():
            after=c.execute('SELECT id,formula,queue,c_count,h_count,c_assigned,h_assigned,task_done FROM molecules WHERE id=?',(mid,)).fetchone()
            if after is not None and before['queue']!=after['queue']:
                transitions.append({'molecule_id':mid,'formula':after['formula'],'before':{k:before[k] for k in ('queue','c_count','h_count','c_assigned','h_assigned','task_done')},'after':dict(after)})
        result['queue_transitions_since_first_repair']=transitions
        result['queue_transition_counts']=dict(Counter(r['before']['queue']+' -> '+r['after']['queue'] for r in transitions))
        delta=Counter()
        for r in transitions:
            delta[r['before']['queue']]+=1
            delta[r['after']['queue']]-=1
        current_counts=Counter()
        for r in result['queues']:current_counts[r['queue']]+=r['n']
        result['reconstructed_before_queue_counts']={k:current_counts[k]+delta[k] for k in current_counts.keys()|delta.keys()}
    result['latest_saved_manual_at_local']=max(datetime.fromisoformat(r['updated_at']).astimezone(timezone(timedelta(hours=8))).isoformat() for r in manual)
    result['unique_manual_molecules']=len({r['molecule_id'] for r in manual})
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('remaining_before_queue_counts','reopened_with_manual_labels','reopened_issue_combinations','manual_labels_missing_from_current_peak_table','latest_saved_manual_at_local','unique_manual_molecules','queue_transition_counts','reconstructed_before_queue_counts')},ensure_ascii=False),flush=True)
    print(json.dumps({'report':str(OUT)},ensure_ascii=False),flush=True)


if __name__=='__main__':
    main()
