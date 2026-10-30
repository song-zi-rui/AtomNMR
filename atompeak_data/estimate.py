"""Read-only preflight: serialize in memory; never create training LMDB files."""
from __future__ import annotations

import json
import math
import pickle
import sqlite3
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from annotation_tool_v3 import core as db
from .export import prepare_record
from .schema import split_for_connectivity
from .tokenization import ROOT, Vocabulary


def estimate_training_dataset(database,*,molecule_ids=None,seed=42,valid_ratio=.1,test_ratio=.1,
        max_seq_len=1024,max_target_len=512,use_solvent=False,use_temperature=False,
        use_field_strength=False,require_supervision=False,progress=None):
    split_for_connectivity('validate',seed,valid_ratio,test_ratio)
    if not 8 <= min(max_seq_len,max_target_len) <= max(max_seq_len,max_target_len) <= 8192:
        raise ValueError('Training sequence lengths must be between 8 and 8192')
    options=dict(max_seq_len=max_seq_len,max_target_len=max_target_len,use_formula=True,
        use_h_category=True,use_h_nh=True,use_h_jvalue=True,use_range_shift=False,
        use_solvent=use_solvent,use_temperature=use_temperature,use_field_strength=use_field_strength)
    source=Vocabulary(ROOT/'dict/bart/condition_dict.txt');target=Vocabulary(ROOT/'dict/bart/smiles_dict.txt')
    started=time.monotonic();c=sqlite3.connect(f'file:{Path(database).resolve().as_posix()}?mode=ro',uri=True,timeout=60)
    c.row_factory=sqlite3.Row;c.execute('BEGIN')
    where=''
    if molecule_ids is not None:
        c.execute('CREATE TEMP TABLE selection(id TEXT PRIMARY KEY)')
        c.executemany('INSERT OR IGNORE INTO selection VALUES(?)',((mid,) for mid in molecule_ids))
        where=' WHERE id IN (SELECT id FROM selection)'
    candidate_counts=[dict(r) for r in c.execute("SELECT CASE WHEN c_count>0 AND h_count>0 THEN 'dual' WHEN c_count>0 THEN 'carbon_only' WHEN h_count>0 THEN 'proton_only' ELSE 'none' END mode,CASE WHEN c_assigned+h_assigned=0 THEN 'none' WHEN c_assigned+h_assigned=c_count+h_count THEN 'full' ELSE 'partial' END coverage,count(*) count FROM molecules"+where+' GROUP BY 1,2')]
    total=sum(v['count'] for v in candidate_counts)
    accepted=Counter();counts=Counter();skipped=Counter();splits=Counter();payloads=Counter();pages=Counter();index_bytes=0;processed=0
    try:
        for row in c.execute('SELECT id FROM molecules'+where+' ORDER BY id'):
            processed+=1
            try:
                raw=db.record(c,row['id'],require_quality=True)
                value,labelled=prepare_record(raw,source,target,options,require_supervision)
                split=split_for_connectivity(value['connectivity_key'],seed,valid_ratio,test_ratio);value['split']=split
                size=len(pickle.dumps(value,protocol=4));payloads[split]+=size
                # LMDB overflow pages plus a conservative B-tree overhead.
                pages[split]+=math.ceil((size+16)/4096)
                splits[split]+=1
                mode='dual' if all(value['branch_present'].values()) else 'carbon_only' if value['branch_present']['carbon'] else 'proton_only'
                complete=True
                hard=weak=sets=0
                for branch in ('carbon','proton'):
                    labels=value['alignment_supervision'][branch]
                    hard+=sum(any(r) for r in labels['pair_known_mask']);weak+=sum(labels['allowed_row_mask']);sets+=len(labels['set_peak_membership'])
                    for i,status in enumerate(labels['assignment_status_by_peak']):
                        if status not in db.ASSIGNED or not (any(labels['pair_known_mask'][i]) or labels['allowed_row_mask'][i]):complete=False
                coverage='full' if complete else 'partial'
                accepted[mode+':'+coverage]+=1
                counts['with_hard_labels' if hard else 'weak_or_set_only' if weak or sets else 'without_supervision']+=1
                counts['hard_peaks']+=hard;counts['weak_peaks']+=weak;counts['unordered_sets']+=sets
                index_bytes+=len((json.dumps({'id':value['id'],'split':split,'connectivity_key':value['connectivity_key'],
                    'mode':mode,'has_supervision':labelled,'revision':raw['revision'],
                    'selected_spectrum_ids':value['provenance']['selected_spectrum_ids']},ensure_ascii=False)+'\n').encode())
            except (ValueError,KeyError,TypeError,RuntimeError) as exc:
                skipped[str(exc)]+=1
            if processed%1000==0 and progress:
                progress({'processed':processed,'total':total,'eligible':sum(splits.values()),'elapsed_seconds':round(time.monotonic()-started,1)})
    finally:c.close()
    elapsed=time.monotonic()-started
    used={s:4096*(pages[s]+math.ceil(splits[s]*48/4096/.7)+4) for s in ('train','valid','test')}
    reserved={s:32*1024**2 for s in used}
    for s in reserved:
        while reserved[s]<used[s]*1.05:reserved[s]*=2
    vocab_bytes=sum((ROOT/'dict/bart'/name).stat().st_size for name in ('condition_dict.txt','smiles_dict.txt'))
    extra=index_bytes+vocab_bytes+1024**2
    return {'mode':'read_only_preflight_no_lmdb_written','created_at':datetime.now(timezone.utc).isoformat(),
        'database':str(Path(database).resolve()),'require_supervision':require_supervision,'tokenization':options,
        'selected_candidates':total,'candidate_counts':candidate_counts,'eligible_records':sum(splits.values()),
        'eligible_counts':dict(accepted),'supervision_counts':dict(counts),'split_counts':dict(splits),
        'skipped_records':sum(skipped.values()),'skipped_reasons':dict(skipped),
        'serialized_payload_bytes':sum(payloads.values()),'estimated_used_bytes':sum(used.values())+extra,
        'estimated_windows_file_bytes':sum(reserved.values())+extra,
        'estimated_per_split_bytes':{'used':used,'windows_allocated':reserved},
        'preflight_seconds':round(elapsed,2),'estimated_export_seconds_range':[round(elapsed*1.2),round(elapsed*2.5)],
        'notes':['All candidates checked in a read-only SQLite snapshot; records serialized in memory, not exported.',
                 'Counts are exact for this preflight snapshot; output size and time are estimates, not promises.',
                 'Full coverage uses saved status, including source_group; this does not certify spectral completeness or chemical correctness.',
                 'Windows LMDB files grow by map-size doubling; estimated allocated bytes can exceed serialized payload.',
                 'Export time excludes ZIP compression; simultaneous annotation, antivirus and storage load can change timing.']}
