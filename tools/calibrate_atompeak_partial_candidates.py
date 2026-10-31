"""Fit only score weights on complete frozen VALID predictions, without a model."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pickle
from pathlib import Path
import sys

os.environ['NMRPEAK_LIGHT_IMPORT']='1'
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from nmrpeak.atompeak.candidate_policy_partial import FEATURES,GENERATION_ONLY,NORMALIZATION,normalized_features,validate_candidate
from nmrpeak.atompeak.generation_metrics_partial import parse_structure


def weight_grid(divisions):
    if not 1<=divisions<=20:raise ValueError('Grid divisions must be between 1 and 20')
    yield GENERATION_ONLY.copy()
    for g in range(divisions+1):
        for a in range(divisions-g+1):
            for u in range(divisions-g-a+1):
                values=(g,a,u,divisions-g-a-u)
                if values==(divisions,0,0,0):continue
                yield dict(zip(FEATURES,(v/divisions for v in values)))


def calibrate(predictions,data,output,divisions=4):
    import lmdb
    if output.exists():raise ValueError('Choose a new profile file')
    run_path=predictions.with_suffix('.run.json');run=json.loads(run_path.read_text(encoding='utf-8'))
    manifest=json.loads((data/'manifest.json').read_text(encoding='utf-8'))
    if run.get('split')!='valid':raise ValueError('Calibration is allowed on VALID only')
    if run['dataset_fingerprint']!=manifest['source_snapshot_fingerprint']:raise ValueError('Validation snapshot differs')
    mode=run['generation_policy']['input_mode']
    with lmdb.open(str(data/'lmdb_dataset/valid.lmdb'),subdir=False,readonly=True,lock=False) as env:
        with env.begin() as tx:records={r['id']:r for _,blob in tx.cursor() if (r:=pickle.loads(blob))}
    expected={mid for mid,r in records.items() if mode=='as_stored' or r['branch_present']['carbon' if mode=='carbon_only' else 'proton']}
    seen=set();prepared=[]
    with predictions.open(encoding='utf-8') as handle:
        for line in handle:
            if not line.strip():continue
            row=json.loads(line);mid=row['id']
            if row.get('split')!='valid' or mid not in expected or mid in seen:raise ValueError('Duplicate, wrong split or unknown validation ID: '+mid)
            seen.add(mid);record=records[mid];truth=parse_structure(record['smiles_canonical'])
            reported=parse_structure(row['target_for_evaluation_only'])
            if reported is None or reported['canonical']!=truth['canonical']:raise ValueError('Validation target differs: '+mid)
            pool=row['generation_diagnostics']['candidate_pool']
            keys=[]
            for candidate in pool:
                checked,reason=validate_candidate(candidate['smiles'],record['molecular_formula'])
                if reason or checked['canonical_smiles']!=candidate['canonical_smiles']:raise ValueError('Candidate failed strict revalidation')
                if not candidate.get('score_features_complete'):raise ValueError('Run VALID inference with --collect-scores first')
                if any(k not in candidate for k in ('generation_score','alignment_compatibility','uot_cost','shift_error')):raise ValueError('Missing score feature')
                keys.append(candidate['canonical_smiles'])
            if len(keys)!=len(set(keys)):raise ValueError('Duplicate candidate in strict pool')
            prepared.append((normalized_features(pool),np.asarray(keys),np.asarray([c['generation_score'] for c in pool]),truth['canonical']))
    if seen!=expected or not seen:raise ValueError('Complete nonempty eligible VALID set required; missing '+str(len(expected-seen)))
    results=[];best=None
    for weights in weight_grid(divisions):
        hits={k:0 for k in (1,3,10)}
        for matrix,keys,generation,truth in prepared:
            scores=matrix@np.asarray([weights[k] for k in FEATURES])
            order=np.lexsort((keys,-generation,-scores))
            for k in hits:hits[k]+=int(truth in keys[order[:k]])
        item={'weights':weights,**{'top'+str(k):v for k,v in hits.items()}}
        results.append(item)
        # Keep the simpler generation-only policy on a tie; don't force extra terms.
        if best is None or item['top1']>best['top1']:best=item
    profile={'schema':'atompeak-partial-calibrated-candidates-v1','calibrated_split':'valid','complete_validation':True,
        'normalization':NORMALIZATION,'weights':best['weights'],'samples':len(seen),
        'checkpoint_sha256':run['checkpoint_sha256'],'dataset_fingerprint':run['dataset_fingerprint'],
        'generation_policy':run['generation_policy'],'validation_file':str(predictions.resolve()),
        'validation_sha256':hashlib.sha256(predictions.read_bytes()).hexdigest(),
        'baseline':results[0],'selected':best,'grid_results':results,
        'note':'Selected on complete VALID only. These are calibration scores, not an independent test improvement claim. '
               'UOT already incorporates shift and compatibility terms; zero auxiliary weights are allowed.'}
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x',encoding='utf-8') as handle:json.dump(profile,handle,indent=2)
    return profile


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--predictions',type=Path,required=True);p.add_argument('--data',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--grid-divisions',type=int,default=4)
    args=p.parse_args();profile=calibrate(args.predictions,args.data,args.output,args.grid_divisions)
    print(json.dumps({'weights':profile['weights'],'samples':profile['samples'],'baseline':profile['baseline'],
                      'selected':profile['selected'],'output':str(args.output.resolve())}))


if __name__=='__main__':main()
