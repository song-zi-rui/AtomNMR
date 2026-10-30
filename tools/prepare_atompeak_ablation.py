"""Prepare train-cohort/solvent ablations from one frozen partial-v1 LMDB export.

Only TRAIN membership and the solvent input switch vary. VALID/TEST observations,
labels, IDs and original structure splits are preserved. No live annotation DB is read.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import shutil
import sys
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
COHORTS=('dual_complete','dual_all','mixed_complete','mixed_all')
ARMS=tuple(c+'_'+s for c in COHORTS for s in ('nosolv','solv'))


def classify(record):
    present=[b for b in ('carbon','proton') if record['branch_present'][b]]
    dual=len(present)==2
    complete=bool(present) and all(record['alignment_supervision'][b]['label_complete_by_branch'] for b in present)
    return dual,complete


def retained(record,cohort):
    dual,complete=classify(record)
    return (not cohort.startswith('dual') or dual) and (not cohort.endswith('complete') or complete)


def write_json(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8')


def inspect(source):
    import lmdb
    manifest=json.loads((source/'manifest.json').read_text(encoding='utf-8'))
    if manifest.get('schema_version')!='atompeak-partial-v1' or manifest.get('state')!='complete':
        raise ValueError('Source must be a completed atompeak-partial-v1 export')
    if not manifest.get('require_usable_assignment_supervision'):
        raise ValueError('This cohort experiment expects the annotated source export; do not mix in a different unlabeled pool')
    for name,expected in manifest['vocab_sha256'].items():
        if hashlib.sha256((source/'dict/bart'/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('Source vocabulary differs from its manifest')
    groups={};summary={}; ids={};duplicate_ids=set()
    for split in ('train','valid','test'):
        counts=Counter();ids[split]=[]
        with lmdb.open(str(source/'lmdb_dataset'/f'{split}.lmdb'),subdir=False,readonly=True,lock=False) as env:
            with env.begin() as tx:
                for _,blob in tx.cursor():
                    r=pickle.loads(blob)
                    if r.get('split')!=split:raise ValueError('Stored split differs from file')
                    if r['id'] in duplicate_ids:raise ValueError('Duplicate molecule ID in source')
                    duplicate_ids.add(r['id']);ids[split].append(r['id'])
                    groups.setdefault(r['connectivity_key'],set()).add(split)
                    for cohort in COHORTS:
                        if retained(r,cohort):counts[cohort]+=1
        summary[split]={'records':len(ids[split]),'cohorts':dict(counts),
                        'ordered_ids_sha256':hashlib.sha256('\n'.join(ids[split]).encode()).hexdigest()}
    if any(len(v)>1 for v in groups.values()):raise ValueError('Connectivity leakage across original splits')
    return manifest,summary


def make_arm(source,root,arm,manifest):
    import lmdb
    cohort,solvent=arm.rsplit('_',1)
    target=root/arm
    if target.exists():raise ValueError(f'Output already exists; preserved: {target}')
    target.mkdir(parents=True);(target/'lmdb_dataset').mkdir();(target/'dict/bart').mkdir(parents=True)
    for name in manifest['vocab_sha256']:shutil.copyfile(source/'dict/bart'/name,target/'dict/bart'/name)
    out=dict(manifest)
    out.update(state='running',ablation={'cohort':cohort,'solvent_enabled':solvent=='solv',
        'source_manifest':str((source/'manifest.json').resolve()),
        'source_fingerprint':manifest['source_snapshot_fingerprint'],
        'changes_train_membership_only':True,'valid_test_preserved':True,
        'complete_definition':'all_observed_peaks_have_hard_assignments_in_every_present_branch_not_spectral_completeness',
        'interpretation':'Adding molecules changes structures and label counts; not a same-molecule partial-label masking experiment'},
        tokenization={**manifest['tokenization'],'use_solvent':solvent=='solv'})
    write_json(target/'manifest.json',out)
    counts=Counter();classes=Counter();index=[];fingerprints=[]
    try:
        for split in ('train','valid','test'):
            src=source/'lmdb_dataset'/f'{split}.lmdb';dst=target/'lmdb_dataset'/f'{split}.lmdb'
            with lmdb.open(str(src),subdir=False,readonly=True,lock=False,readahead=True) as e:
                with lmdb.open(str(dst),subdir=False,map_size=max(src.stat().st_size*2,16*1024**2)) as writer:
                    with e.begin() as tx,writer.begin(write=True) as tw:
                        for _,blob in tx.cursor():
                            r=pickle.loads(blob)
                            if split=='train' and not retained(r,cohort):continue
                            # Copy original bytes, including labels and source revision, verbatim.
                            tw.put(f'{counts[split]:012d}'.encode(),blob);counts[split]+=1
                            fingerprints.append((r['id'],split,counts[split]-1))
                            mode='dual' if all(r['branch_present'].values()) else 'carbon_only' if r['branch_present']['carbon'] else 'proton_only'
                            classes[mode]+=1
                            labels=r['alignment_supervision']
                            has_labels=any(any(any(row) for row in labels[b]['pair_known_mask']) or any(labels[b]['allowed_row_mask']) or labels[b]['set_peak_membership'] for b in ('carbon','proton'))
                            index.append({'id':r['id'],'split':split,'connectivity_key':r['connectivity_key'],'mode':mode,
                                'has_supervision':bool(has_labels),'revision':r['provenance']['molecule_revision'],
                                'selected_spectrum_ids':r['provenance']['selected_spectrum_ids']})
                            for b in ('carbon','proton'):
                                classes[b+'_hard_peaks']+=sum(any(row) for row in labels[b]['pair_known_mask'])
                                classes[b+'_weak_peaks']+=sum(labels[b]['allowed_row_mask'])
                                classes[b+'_unordered_sets']+=len(labels[b]['set_peak_membership'])
            if not counts[split]:raise ValueError(f'Empty split: {split}')
        # Reproduce the original export fingerprint convention in global ID order.
        handles={};transactions={};h=hashlib.sha256()
        try:
            for split in ('train','valid','test'):
                handles[split]=lmdb.open(str(target/'lmdb_dataset'/f'{split}.lmdb'),subdir=False,readonly=True,lock=False)
                transactions[split]=handles[split].begin()
            for mid,split,key in sorted(fingerprints):h.update(mid.encode()+transactions[split].get(f'{key:012d}'.encode()))
        finally:
            for tx in transactions.values():tx.abort()
            for e in handles.values():e.close()
        index.sort(key=lambda r:r['id'])
        (target/'records_index.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in index),encoding='utf-8')
        out.update(state='complete',records=sum(counts.values()),split_counts=dict(counts),classes=dict(classes),
                   source_snapshot_fingerprint=h.hexdigest(),selected_molecules=sum(counts.values()),
                   filters={'ablation_train_cohort':cohort},skipped_reasons={},warnings=[])
        write_json(target/'manifest.json',out)
        return {'arm':arm,'path':str(target.resolve()),'split_counts':dict(counts),'solvent':solvent=='solv'}
    except Exception as exc:
        out.update(state='failed',error=str(exc));write_json(target/'manifest.json',out);raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output-root',type=Path)
    p.add_argument('--arms',nargs='+',choices=ARMS,default=list(ARMS))
    p.add_argument('--inspect-only',action='store_true')
    p.add_argument('--report',type=Path)
    args=p.parse_args()
    manifest,summary=inspect(args.source)
    result={'source':str(args.source.resolve()),'source_fingerprint':manifest['source_snapshot_fingerprint'],
            'source_summary':summary,'arms':[]}
    if not args.inspect_only:
        if args.output_root is None:p.error('--output-root is required unless --inspect-only')
        if len(set(args.arms))!=len(args.arms):p.error('Duplicate arms')
        for arm in args.arms:
            if (args.output_root/arm).exists():raise ValueError('An output arm exists; choose another root or select only new arms')
        for arm in args.arms:
            value=make_arm(args.source,args.output_root,arm,manifest);result['arms'].append(value)
            print(json.dumps(value,ensure_ascii=False),flush=True)
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True);write_json(args.report,result)
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
