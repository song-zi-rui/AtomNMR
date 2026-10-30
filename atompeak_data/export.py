"""Snapshot export to split, directly trainable LMDB files."""
from __future__ import annotations

import hashlib
import json
import pickle
import shutil
import sqlite3
import time
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .schema import SCHEMA_VERSION, build_record, split_for_connectivity
from .tokenization import ROOT, Vocabulary, tokenize_spectrum


def has_supervision(record):
    return any(any(any(row) for row in record['alignment_supervision'][b]['pair_known_mask'])
        or any(record['alignment_supervision'][b]['allowed_row_mask'])
        or record['alignment_supervision'][b]['set_peak_membership'] for b in ('carbon','proton'))


def prepare_record(raw, source_vocab, target_vocab, tokenization, require_supervision=False):
    value=build_record(raw)
    labelled=bool(has_supervision(value))
    if require_supervision and not labelled:
        raise ValueError('Selected current spectra contain no usable assignment supervision')
    tokenize_spectrum(value,source_vocab,**{k:v for k,v in tokenization.items()
        if k not in ('max_seq_len','max_target_len')},max_length=tokenization['max_seq_len'])
    if len(target_vocab.encode_smiles(value['smiles_canonical']))>tokenization['max_target_len']:
        raise ValueError('SMILES target exceeds configured maximum; not truncated')
    return value,labelled


def export_training_dataset(database, folder, *, filters=None, molecule_id=None,
                            molecule_ids=None, seed=42, valid_ratio=.1, test_ratio=.1,
                            max_seq_len=1024, max_target_len=512, archive=False, progress=None,
                            use_solvent=False, use_temperature=False, use_field_strength=False,
                            require_supervision=False, export_options=None):
    from annotation_tool_v3 import core as db
    from annotation_tool_v3.export_policy import Policy
    policy=Policy(export_options)
    import lmdb
    split_for_connectivity('validate', seed, valid_ratio, test_ratio)
    if molecule_ids is not None and not molecule_ids:
        raise ValueError('Empty molecule selection')
    if not 8 <= min(max_seq_len,max_target_len) <= max(max_seq_len,max_target_len) <= 8192:
        raise ValueError('Training sequence lengths must be between 8 and 8192')
    folder = Path(folder).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    if any(folder.iterdir()):
        raise ValueError('Training output directory must be empty; choose a new dataset version')
    vocab_dir = folder / 'dict/bart'
    vocab_dir.mkdir(parents=True)
    for name in ('condition_dict.txt', 'smiles_dict.txt'):
        shutil.copyfile(ROOT / 'dict/bart' / name, vocab_dir / name)
    source_vocab, target_vocab = Vocabulary(vocab_dir/'condition_dict.txt'), Vocabulary(vocab_dir/'smiles_dict.txt')
    c = sqlite3.connect(f'file:{Path(database).resolve().as_posix()}?mode=ro', uri=True, timeout=60)
    c.row_factory = sqlite3.Row
    c.execute('BEGIN')
    where, values = db.selection_where(filters or {})
    if molecule_id:
        where, values = ' WHERE id=?', [molecule_id]
    elif molecule_ids is not None:
        c.execute('CREATE TEMP TABLE export_selection(id TEXT PRIMARY KEY)')
        c.executemany('INSERT OR IGNORE INTO export_selection VALUES(?)',((mid,) for mid in molecule_ids))
        where, values = ' WHERE id IN (SELECT id FROM export_selection)', []
    total = c.execute('SELECT count(*) FROM molecules' + where, values).fetchone()[0]
    manifest = {'schema_version': SCHEMA_VERSION, 'state': 'running',
        'created_at': datetime.now(timezone.utc).isoformat(), 'database': str(Path(database).resolve()),
        'split_policy': {'method':'sha256_connectivity_group_v1','seed':seed,'valid_ratio':valid_ratio,'test_ratio':test_ratio},
        'grouping_policy':'individual_atoms_v1', 'initialization':'random_from_scratch',
        'tokenization': {'max_seq_len':max_seq_len,'max_target_len':max_target_len,'use_formula':True,
                         'use_h_category':True,'use_h_nh':True,'use_h_jvalue':True,'use_range_shift':False,
                         'use_solvent':use_solvent,'use_temperature':use_temperature,'use_field_strength':use_field_strength},
        'quality_export_enforced':True,'source_candidates_not_fabricated':True,
        'unknown_labels_are_not_negative':True, 'full_assignment_does_not_certify_spectral_completeness':True,
        'filters':filters or {}, 'selection_molecule_id':molecule_id,
        'require_usable_assignment_supervision':require_supervision}
    manifest_path = folder/'manifest.json'
    manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    lmdb_dir = folder/'lmdb_dataset';lmdb_dir.mkdir()
    envs = {s:lmdb.open(str(lmdb_dir/f'{s}.lmdb'),subdir=False,map_size=32*1024**2) for s in ('train','valid','test')}
    batches = {s:[] for s in envs};counts=Counter();classes=Counter();reasons=Counter();scanned=0
    fingerprint=hashlib.sha256();started=time.monotonic()
    def flush(split):
        if not batches[split]:return
        while True:
            try:
                with envs[split].begin(write=True) as txn:
                    for key,value in batches[split]:txn.put(key,value,overwrite=False)
                batches[split].clear();return
            except lmdb.MapFullError:envs[split].set_mapsize(envs[split].info()['map_size']*2)
    def publish():
        info={'processed':scanned,'total':total,'exported':sum(counts.values()),'splits':dict(counts),
              'elapsed_seconds':round(time.monotonic()-started,1)}
        if progress:progress(info)
    try:
        with (folder/'records_index.jsonl').open('w',encoding='utf-8') as index, (folder/'skipped_records.jsonl').open('w',encoding='utf-8') as skipped:
            for row in c.execute('SELECT id FROM molecules'+where+' ORDER BY id',values):
                scanned+=1
                try:
                    raw=db.record(c,row['id'],require_quality=True)
                    if not policy.accept(raw):raise ValueError('export_policy_filtered')
                    record,has_labels=prepare_record(raw,source_vocab,target_vocab,manifest['tokenization'],require_supervision)
                    split=split_for_connectivity(record['connectivity_key'],seed,valid_ratio,test_ratio)
                    record['split']=split
                    payload=pickle.dumps(record,protocol=4)
                    fingerprint.update(record['id'].encode()+payload)
                    batches[split].append((f'{counts[split]:012d}'.encode(),payload));counts[split]+=1
                    mode='dual' if all(record['branch_present'].values()) else 'carbon_only' if record['branch_present']['carbon'] else 'proton_only'
                    classes[mode]+=1
                    for branch in ('carbon','proton'):
                        labels=record['alignment_supervision'][branch]
                        classes[branch+'_hard_peaks']+=sum(any(r) for r in labels['pair_known_mask'])
                        classes[branch+'_weak_peaks']+=sum(labels['allowed_row_mask'])
                        classes[branch+'_unordered_sets']+=len(labels['set_peak_membership'])
                    index.write(json.dumps({'id':record['id'],'split':split,'connectivity_key':record['connectivity_key'],
                        'mode':mode,'has_supervision':bool(has_labels),'revision':raw['revision'],
                        'selected_spectrum_ids':record['provenance']['selected_spectrum_ids']},ensure_ascii=False)+'\n')
                    if len(batches[split])>=200:flush(split)
                except (ValueError,KeyError,TypeError,RuntimeError) as exc:
                    reason=str(exc);reasons[reason]+=1
                    skipped.write(json.dumps({'id':row['id'],'reason':reason},ensure_ascii=False)+'\n')
                if scanned%200==0:publish()
            for split in envs:flush(split);envs[split].sync()
        manifest.update(state='complete',records=sum(counts.values()),selected_molecules=total,
                        split_counts={s:counts[s] for s in envs},classes=dict(classes),skipped_reasons=dict(reasons),
                        source_snapshot_fingerprint=fingerprint.hexdigest(),elapsed_seconds=round(time.monotonic()-started,2),
                        export_selection=policy.report(),
                        vocab_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in vocab_dir.iterdir()},
                        warnings=['A split is empty; use a larger selection for training.'] if any(not counts[s] for s in envs) else [])
    except Exception as exc:
        manifest.update(state='failed',error=str(exc));raise
    finally:
        for env in envs.values():env.close()
        c.close()
        manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    publish()
    target=manifest_path
    if archive:
        target=folder/'atompeak_dataset.zip'
        with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1) as z:
            for file in sorted(folder.rglob('*')):
                if file.is_file() and file!=target and not file.name.endswith('-lock'):
                    z.write(file,file.relative_to(folder))
    return {'file':str(target),'dataset_dir':str(folder),'manifest':str(manifest_path),
            'count':manifest['records'],'split_counts':manifest['split_counts'],'skipped':total-manifest['records'],'selection':policy.report()}
