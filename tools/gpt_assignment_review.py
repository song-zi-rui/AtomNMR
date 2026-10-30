"""Persist explicit GPT review decisions; never generate chemical labels automatically."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rdkit import Chem

SOURCE = ROOT / 'data/shift_assignment_assistant/new_targets_200.with_suggestions.v2.jsonl'
DATA = ROOT / 'data/atom_assignment_annotation_v2'
OUT = ROOT / 'data/gpt_assignment_review/batch_20260916'
IDS = ['NMRX-0015','NMRX-0163','NMRX-0037','NMRX-0073','NMRX-0088','NMRX-0137','NMRX-0082','NMRX-0130','NMRX-0022','NMRX-0034']

def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')

def digest(value):
    return hashlib.sha256(json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def atomic(path, text):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    with temp.open('w',encoding='utf-8',newline='\n') as f:
        f.write(text); f.flush(); os.fsync(f.fileno())
    os.replace(temp,path)
    if path.read_text(encoding='utf-8') != text:
        raise IOError(f'readback mismatch: {path}')

def save(path,value):
    atomic(path,json.dumps(value,ensure_ascii=False,indent=2)+'\n')

def rows(path):
    return [json.loads(s) for s in path.read_text(encoding='utf-8').splitlines() if s.strip()]

def samples():
    return {r['sample_id']:r for r in rows(SOURCE)}

def task_map():
    return {r['sample_id']:f'TASK-{i:04d}' for i,r in enumerate(json.loads((DATA/'annotation_sessions.json').read_text(encoding='utf-8')),1) if r['role']=='new'}

def molecule(sample):
    mol=Chem.MolFromSmiles(sample['structure']['smiles_canonical'])
    atoms=sample['structure']['atoms']
    assert mol is not None and mol.GetNumAtoms()==len(atoms)
    for row in atoms:
        a=mol.GetAtomWithIdx(row['index'])
        assert a.GetSymbol()==row['symbol'] and a.GetTotalNumHs()==row['total_hydrogens']
        assert sorted(n.GetIdx() for n in a.GetNeighbors())==sorted(row['neighbors'])
        a.SetAtomMapNum(a.GetIdx()+1)
    return mol

def prepare():
    from annotation_tool_v2.review_profile import build_review_profile
    pool=samples(); tasks=task_map()
    for sid in IDS:
        s=pool[sid]; mol=molecule(s); aa=s['assignment_assistant']; best=aa['solutions'][0]
        packet={'sample_id':sid,'task_id':tasks[sid],'structure':s['structure'], 'conditions':s['conditions'],
                'spectra':s['spectra'],'source':s['source'],'provenance':s['provenance'],
                'mapped_smiles':Chem.MolToSmiles(mol),'review_profile':build_review_profile(s),
                'atom_environments':[{'index':a.GetIdx(),'symbol':a.GetSymbol(),'H':a.GetTotalNumHs(),
                    'aromatic':a.GetIsAromatic(),'hybridization':str(a.GetHybridization()),
                    'neighbors':[(n.GetIdx(),n.GetSymbol(),str(mol.GetBondBetweenAtoms(a.GetIdx(),n.GetIdx()).GetBondType())) for n in a.GetNeighbors()]} for a in mol.GetAtoms()],
                'predictions':[{'atoms':g['atom_indices'],'prediction':g.get('predicted_shift_ppm'),
                    'sigma':g.get('prediction_sigma_ppm'),'members':[{'atom':m['atom_index'],'prediction':m['predicted_shift_ppm'],
                       'environment':m['environment']} for m in g.get('member_predictions',[])]} for g in aa['groups']],
                'best_mapping':best,'alternative_solutions':aa['solutions'][1:3],'ambiguities':aa.get('ambiguity_groups',[])}
        save(OUT/'packets'/f'{sid}.json',packet)
    save(OUT/'selection.json',{'sample_ids':IDS,'task_ids':{i:tasks[i] for i in IDS},'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'created_at':now()})
    print(json.dumps({i:tasks[i] for i in IDS},ensure_ascii=False))

def validate(sample,decision):
    mol=molecule(sample); assignments=decision['assignments']; seen_c=set()
    if not decision.get('summary') or not decision.get('limitations'):
        raise ValueError('summary and limitations required')
    for nucleus in ['13C','1H']:
        peaks={p['peak_id']:p for p in sample['spectra'][nucleus]['peaks']}
        if set(assignments[nucleus])!=set(peaks):
            raise ValueError(f'{nucleus} missing or extra peak decisions')
        for pid,d in assignments[nucleus].items():
            if d['status'] not in {'assigned','equivalent_group','overlap','unresolved','impurity','unobserved'}:
                raise ValueError(f'{pid}: invalid status')
            if d['confidence'] not in {'certain','probable','ambiguous'} or not d.get('rationale'):
                raise ValueError(f'{pid}: confidence/explanation missing')
            assigned=d['status'] in {'assigned','equivalent_group','overlap'}
            indices=d['atom_indices']
            if assigned != bool(indices) or len(indices)!=len(set(indices)):
                raise ValueError(f'{pid}: invalid assigned indices')
            for group in [indices]+d.get('candidate_atom_sets',[]):
                if not isinstance(group,list) or len(group)!=len(set(group)):
                    raise ValueError(f'{pid}: invalid candidate set')
                for i in group:
                    if type(i) is not int or not 0<=i<mol.GetNumAtoms(): raise ValueError('invalid atom index')
                    a=mol.GetAtomWithIdx(i)
                    if nucleus=='13C' and a.GetSymbol()!='C': raise ValueError('C assignment to non-C')
                    if nucleus=='1H' and not a.GetTotalNumHs(): raise ValueError('H assignment to atom without H')
            if nucleus=='13C':
                if seen_c.intersection(indices): raise ValueError('duplicate carbon occupancy')
                seen_c.update(indices)
            if d['training_eligible'] and (not assigned or d['confidence']!='certain'):
                raise ValueError('hard labels must be resolved and certain by declared evidence')
    return mol

def export():
    records=[json.loads(p.read_text(encoding='utf-8')) for p in sorted((OUT/'records').glob('*.json'))]
    all_labels=[]; hard=[]; unresolved=0; total=0
    for r in records:
        for nucleus,table in r['review']['assignments'].items():
            peaks={p['peak_id']:p for p in r['sample']['spectra'][nucleus]['peaks']}
            for pid,d in table.items():
                total+=1
                if d['status']=='unresolved': unresolved+=1
                if not d['atom_indices']: continue
                item={'sample_id':r['sample_id'],'task_id':r['task_id'],'nucleus':nucleus,'peak_id':pid,
                      'peak':peaks[pid],'structure':r['portable_structure'],'conditions':r['sample']['conditions'],
                      'source':r['sample']['source'],'assignment':d,'atom_map_numbers':[i+1 for i in d['atom_indices']],
                      'label_origin':'gpt_reviewed_pseudo_label','experimental_assignment_verified':False,
                      'split_group':r['portable_structure']['connectivity_key']}
                all_labels.append(item)
                if d['training_eligible']: hard.append(item)
    for name,values in [('annotations.jsonl',records),('training_resolved_all.jsonl',all_labels),('training_hard_labels.jsonl',hard)]:
        atomic(OUT/name,''.join(json.dumps(v,ensure_ascii=False,separators=(',',':'))+'\n' for v in values))
    progress={'saved_molecules':len(records),'reviewed_peaks':total,'resolved_peaks':len(all_labels),'unresolved_peaks':unresolved,
              'hard_label_peaks':len(hard),'completed_sample_ids':[r['sample_id'] for r in records],'updated_at':now()}
    save(OUT/'progress.json',progress)
    report=['# 首批 10 个分子的 GPT 归属审核', '',
            '这些是 GPT 审核的推定标签，不是实验验证的金标准。certain 是证据等级，不是校准后的正确率。', '',
            '每个分子提交后立即保存并回读，汇总可从 records 重建。氢归属指所连重原子；下表原子号为界面 0-based 编号，训练 mapped SMILES 的 atom-map = 此编号 + 1。', '',
            '本批优先选择峰数较少的分子，不能据此估计全库准确率。没有重新取得二维相关谱；source_checked 表示检查了所保存的来源文本，不代表原论文全部实验已复核。', '',
            f"已保存 {len(records)} 个分子，审核 {total} 条峰记录；已给出归属 {len(all_labels)} 条，保留未解决 {unresolved} 条；严格伪标签 {len(hard)} 条，其余已归属项仍属 probable。", '',
            '| 分子 | 界面任务 | C 已分配/总数 | H 已分配/总数 | 严格伪标签 |',
            '|---|---|---:|---:|---:|']
    for r in records:
        ds=r['review']['assignments']; counts=[]
        for n in ['13C','1H']:
            counts.append(f"{sum(bool(d['atom_indices']) for d in ds[n].values())}/{len(ds[n])}")
        report.append(f"| {r['sample_id']} | {r['task_id']} | {counts[0]} | {counts[1]} | {sum(d['training_eligible'] for t in ds.values() for d in t.values())} |")
    for r in records:
        report.extend(['',f"## {r['sample_id']} / {r['task_id']}",'',r['review']['summary'],'',
                       '限制：'+'；'.join(r['review']['limitations']),'',
                       '| 核种/峰 | ppm | 分配原子 | 备选原子组 | 置信等级 | 原因 |',
                       '|---|---:|---|---|---|---|'])
        for n,ds in r['review']['assignments'].items():
            ps={p['peak_id']:p for p in r['sample']['spectra'][n]['peaks']}
            for pid,d in ds.items():
                report.append(f"| {n}/{pid} | {ps[pid]['shift_ppm']} | {d['atom_indices'] or '未定'} | {d.get('candidate_atom_sets',[])} | {d['confidence']} | {d['rationale'].replace('|','/')} |")
    atomic(OUT/'REVIEW_REPORT.md','\n'.join(report)+'\n')
    return progress

def commit(path):
    decision=json.loads(path.read_text(encoding='utf-8')); sid=decision['sample_id']
    if sid not in IDS: raise ValueError('sample not in pilot selection')
    s=samples()[sid]; mol=validate(s,decision); task=task_map()[sid]
    mapped=Chem.MolToSmiles(mol)
    for atom in mol.GetAtoms(): atom.SetAtomMapNum(0)
    key=Chem.MolToInchiKey(mol)
    record={'schema_version':'gpt-assignment-review-v1','sample_id':sid,'task_id':task,'saved_at':now(),
            'reviewer':{'type':'GPT','execution':'current Codex conversation','independent_experimental_verification':False},
            'input_sample_sha256':digest(s),'decision_sha256':digest(decision),'sample':s,'review':decision,
            'portable_structure':{'mapped_smiles':mapped,'inchikey':key,'connectivity_key':key[:14],
              'atom_index_convention':'UI 0-based; SMILES atom-map=index+1; parse using atom-map, not SMILES position',
              'h_assignment_convention':'heavy-atom parent, not explicit proton identity','atoms':s['structure']['atoms']}}
    dest=OUT/'records'/f'{sid}.json'
    if dest.exists():
        previous=json.loads(dest.read_text(encoding='utf-8'))
        if previous['decision_sha256']!=record['decision_sha256']:
            raise ValueError('existing review differs; preserve it and use an explicit versioned batch')
        record=previous
    else:
        save(dest,record)
        with (OUT/'events.jsonl').open('a',encoding='utf-8') as f:
            f.write(json.dumps({'sample_id':sid,'saved_at':record['saved_at'],'sha256':digest(record)},ensure_ascii=False)+'\n');f.flush();os.fsync(f.fileno())
    progress=export()
    ui_path=DATA/'annotations_v2'/f'{task}.json'
    if not ui_path.exists():
        table={n:{pid:{'status':d['status'],'atom_indices':d['atom_indices'],'confidence':d['confidence'],
             'note':d['rationale']+' 候选：'+json.dumps(d.get('candidate_atom_sets',[])),
             'candidate_atom_sets':d.get('candidate_atom_sets',[]),'label_origin':'gpt_reviewed_pseudo_label'} for pid,d in ds.items()} for n,ds in decision['assignments'].items()}
        unclear=any(d['status']=='unresolved' for ds in decision['assignments'].values() for d in ds.values())
        ann={'schema_version':'atom-assignment-annotation-v2','task_id':task,'workflow_status':'needs_review' if unclear else 'completed',
             'assignments':table,'checks':{'structure_checked':True,'source_checked':True,'ambiguity_checked':True},
             'global_note':'GPT审核标签（非实验金标准）。'+decision['summary']+' 限制：'+'；'.join(decision['limitations']),
             'reviewer':'GPT','server_revision':1,'server_updated_at':now(),'gpt_record_path':str(dest.relative_to(ROOT))}
        ann['content_sha256']=digest(ann)
        # Exclusive creation avoids overwriting an annotation saved concurrently.
        ui_path.parent.mkdir(parents=True,exist_ok=True)
        try:
            with ui_path.open('x',encoding='utf-8') as f:
                json.dump(ann,f,ensure_ascii=False,indent=2);f.flush();os.fsync(f.fileno())
            print(f'UI saved: {task}')
        except FileExistsError:
            print(f'UI preserved existing: {task}')
    print(json.dumps({'committed':sid,'task_id':task,**progress},ensure_ascii=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('command',choices=['prepare','commit','export']);parser.add_argument('--decision',type=Path)
    args=parser.parse_args()
    if args.command=='prepare': prepare()
    elif args.command=='commit': commit(args.decision)
    else: print(json.dumps(export(),ensure_ascii=False))
