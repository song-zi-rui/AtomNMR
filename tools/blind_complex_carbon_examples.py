"""Five complex molecules: forced assignments locked before reference reveal."""
import argparse
import hashlib
import importlib.util
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('blind_helpers',ROOT/'tools/blind_close_carbon_examples.py')
helper=importlib.util.module_from_spec(spec); spec.loader.exec_module(helper)
SOURCE=helper.SOURCE
OUT=ROOT/'data/blind_complex_carbon_examples'

def now(): return datetime.now(timezone.utc).isoformat()
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp')
    with tmp.open('w',encoding='utf-8',newline='\n') as f:
        json.dump(obj,f,ensure_ascii=False,indent=2); f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.replace(tmp,path)
    assert json.loads(path.read_text(encoding='utf-8'))==json.loads(json.dumps(obj))
def read(path): return json.loads(path.read_text(encoding='utf-8'))

def select():
    if (OUT/'selection.json').exists(): raise FileExistsError('Preserve existing blind packets')
    RDLogger.DisableLog('rdApp.*'); candidates=[]; seen=set(); scaffolds=set()
    for row,mol in enumerate(Chem.SDMolSupplier(str(SOURCE),removeHs=False),1):
        if mol is None: continue
        heavy=Chem.RemoveHs(mol)
        nc=sum(a.GetSymbol()=='C' for a in heavy.GetAtoms())
        if not 16<=nc<=24 or '.' in Chem.MolToSmiles(heavy): continue
        if any(a.GetSymbol() not in {'C','N','O','S','F','Cl','Br'} for a in heavy.GetAtoms()): continue
        rings=heavy.GetRingInfo().AtomRings()
        if not rings or not (len(rings)>=3 or max(map(len,rings))>=8): continue
        centers=Chem.FindMolChiralCenters(heavy,includeUnassigned=True)
        if len(centers)<2: continue
        key=Chem.MolToInchiKey(heavy); scaffold=Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(heavy))
        if key in seen or scaffold in scaffolds: continue
        for sp in helper.carbon_spectra(mol):
            method=helper.indexed_metadata(mol,'Assignment Method',sp.spectrum_id)
            program=helper.indexed_metadata(mol,'Program',sp.spectrum_id)
            if 'predict' in str(method).lower() or program: continue
            shifts=sorted(set(s['shift_ppm'] for s in sp.signals))
            close=[(a,b,round(b-a,6)) for a,b in zip(shifts,shifts[1:]) if 0<b-a<=0.200001]
            if not close or len(shifts)!=nc: continue
            # Reference used ONLY for eligibility; never printed or placed in blind packet.
            reference={}
            for s in sp.signals: reference.setdefault(s['shift_ppm'],set()).update(s['atom_indices'])
            ci={a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol()=='C'}
            if any(len(v)!=1 for v in reference.values()) or set.union(*reference.values())!=ci: continue
            cid=f'CX{len(candidates)+1:02d}'
            mapped=Chem.Mol(mol)
            for a in mapped.GetAtoms(): a.SetAtomMapNum(a.GetIdx()+1)
            blind={'example_id':cid,'source_record':row,'database_id':mol.GetProp('nmrshiftdb2 ID') if mol.HasProp('nmrshiftdb2 ID') else '',
                   'spectrum_field':sp.source_field,'smiles':Chem.MolToSmiles(heavy),'mapped_smiles':Chem.MolToSmiles(mapped),
                   'carbon_count':nc,'ring_count':len(rings),'stereocenter_count':len(centers),
                   'atoms':[{'atom':a.GetIdx()+1,'element':a.GetSymbol(),'aromatic':a.GetIsAromatic(),'H':a.GetTotalNumHs(includeNeighbors=True),
                     'neighbors':[[n.GetIdx()+1,str(mol.GetBondBetweenAtoms(a.GetIdx(),n.GetIdx()).GetBondType())] for n in a.GetNeighbors()]} for a in mol.GetAtoms() if a.GetSymbol()!='H'],
                   'carbon_shifts':shifts,'close_pairs':close,'solvent':sp.solvent_raw,'temperature_k':sp.temperature_k,'field_strength_raw_mhz':sp.field_strength_mhz,
                   'assignment_method':method,'protocol':'Only structure, sorted carbon shifts, conditions. No reference assignments, model predictions, HOSE DB, compound-specific external lookup or other spectra.'}
            save(OUT/'blind'/(cid+'.json'),blind)
            save(OUT/'hidden_reference'/(cid+'.json'),{'signals':list(sp.signals)})
            drawing=Chem.RemoveHs(mapped)
            Draw.MolToFile(drawing,str(OUT/'blind'/(cid+'.png')),size=(1200,900))
            candidates.append({k:blind[k] for k in ['example_id','database_id','smiles','carbon_count','ring_count','stereocenter_count','close_pairs']})
            seen.add(key); scaffolds.add(scaffold); break
        if len(candidates)>=10: break
    save(OUT/'selection.json',{'created_at':now(),'source_sha256':sha(SOURCE),'criteria':'16-24 C; >=3 rings or ring size>=8; >=2 potential stereocenters; full one-to-one 13C reference; pair gap<=0.20; first10 distinct isomeric Murcko scaffolds','candidates':candidates})
    print(json.dumps(candidates,ensure_ascii=False,indent=2))

def lock(cid):
    path=OUT/'decisions'/(cid+'.json'); d=read(path); b=read(OUT/'blind'/(cid+'.json'))
    assert d['example_id']==cid and d['summary']
    assert sorted(x['ppm'] for x in d['assignments'])==b['carbon_shifts']
    assert len(d['assignments'])==len(b['carbon_shifts'])
    carbon={a['atom'] for a in b['atoms'] if a['element']=='C'}
    predicted=[x['atom'] for x in d['assignments']]
    assert len(predicted)==len(set(predicted)) and set(predicted)==carbon
    for a in d['assignments']: assert a['reason'] and a['confidence'] in ['high','medium','low']
    target=OUT/'locked'/(cid+'.json')
    if target.exists(): assert read(target)['decision_sha256']==sha(path)
    else: save(target,{'locked_at':now(),'decision_sha256':sha(path),'blind_sha256':sha(OUT/'blind'/(cid+'.json')),'decision':d})
    print(json.dumps({'locked':cid,'carbon_assignments':len(predicted),'sha256':sha(path)}))

def reveal(ids):
    assert len(ids)==5 and len(set(ids))==5
    # Validate EVERY lock before reading ANY answer.
    for cid in ids:
        l=read(OUT/'locked'/(cid+'.json'))
        assert sha(OUT/'decisions'/(cid+'.json'))==l['decision_sha256']
        assert sha(OUT/'blind'/(cid+'.json'))==l['blind_sha256']
    assert sha(SOURCE)==read(OUT/'selection.json')['source_sha256']
    results=[]
    for cid in ids:
        d=read(OUT/'locked'/(cid+'.json'))['decision']; b=read(OUT/'blind'/(cid+'.json'))
        ref=read(OUT/'hidden_reference'/(cid+'.json')); by={}
        for s in ref['signals']: by.setdefault(s['shift_ppm'],set()).update(i+1 for i in s['atom_indices'])
        rows=[{**a,'reference_atom':next(iter(by[a['ppm']])),'correct':a['atom'] in by[a['ppm']]} for a in d['assignments']]
        near={p for a,b,g in b['close_pairs'] for p in [a,b]}
        results.append({'example_id':cid,'database_id':b['database_id'],'total':len(rows),'correct':sum(r['correct'] for r in rows),
                        'near_total':sum(r['ppm'] in near for r in rows),'near_correct':sum(r['correct'] and r['ppm'] in near for r in rows),'rows':rows})
    save(OUT/'comparison.json',{'revealed_at':now(),'results':results})
    print(json.dumps(results,ensure_ascii=False,indent=2))

def report():
    comparison=read(OUT/'comparison.json'); results=comparison['results']
    assert len(results)==5
    total=sum(r['total'] for r in results); correct=sum(r['correct'] for r in results)
    near_total=sum(r['near_total'] for r in results); near_correct=sum(r['near_correct'] for r in results)
    confidence={c:{'total':sum(a['confidence']==c for r in results for a in r['rows']),
                    'correct':sum(a['confidence']==c and a['correct'] for r in results for a in r['rows'])} for c in ['high','medium','low']}
    titles={'CX01':'含环氧基的大环内酯及羟基酰基侧链','CX02':'六环含氮笼状骨架及甲酯','CX03':'联芳双香豆素样骨架的脱氧糖苷','CX04':'含邻二醇的大环二烯内酯及甲基丙烯酸酯','CX05':'三环桥连萜类乙酸酯'}
    diagnoses={
      'CX01':'17/21一致。35.9和36.0的C1/C7分反；81.0和81.5的内酯O-CH与侧链酯化O-CH也分反。其余官能团分组基本正确，错误主要来自相近环境的精细排序。不能把低把握的“稍低场倾向”当成已验证规律。',
      'CX02':'11/20一致。46.4–52.0区间发生五个碳的连锁错配，涉及N-CH2、N-CH、甲酯OCH3、桥头季碳和羰基alpha-CH；此外120.3/120.8以及135.2/136.3互换。并不只是一个0.2 ppm的局部错误。',
      'CX03':'17/24一致。糖碳和主要C-O/羰基分配与库一致，但111.8–115.1区间七个碳全部错配，其中包含CH和季碳；分错后会影响该区间整套一对一排列。结构未指定糖中心构型，不因糖排序这次命中就视为确定的实验归属。',
      'CX04':'5/19一致，是本次表现最差的一例。四个含氧CH全部错配，多个甲基、烯CH和烯季碳也错配。我把83.1给内酯O-CH C10，库中却是醇CH C2；这一事实否定了本次简单按“内酯O-CH更低场”排序的预测，但不构成可推广的反向规则。',
      'CX05':'13/17一致。47.96→C5、48.09→C3这一对0.13 ppm近邻峰分配正确，错误发生在17.43–29.17的三个甲基与一个CH2的排列。说明近邻峰不一定分错，较远的峰也不保证分对。溶剂/温度缺失、数值与通常经验不完全贴合均在揭盲前记录，未因此剔除本例。'}
    lines=['# 五个复杂分子：强制完整13C归属盲测','',
           '## 先看结果','',f'全部 {total} 个碳均给出唯一预测：与NMRShiftDB2一致 {correct}/{total}（{correct/total:.1%}），不一致 {total-correct} 个。没有未定项，没有用候选命中代替正确。', '',
           f'相邻不同位移差 ≤0.20 ppm 的信号共有 {near_total} 个，按信号去重计数，命中 {near_correct}/{near_total}（{near_correct/near_total:.1%}）；五个分子均非整分子全对。', '',
           '## 方法','',
           '先筛16–24个碳、有至少三个RDKit感知环或至少一个8元以上环、至少两个潜在立体中心的分子；要求存在≤0.20 ppm的不同碳信号，并有完整一对一数据库归属。先取得10个不同异构Murcko骨架候选，按出现顺序取前5个，不根据预测成绩替换分子。部分潜在立体中心没有指定构型，不能把数量称为已确认手性中心数。', '',
           '筛选程序内部用数据库标签检查完整性，但不给GPT输出标签或其顺序。GPT只读取结构、排序后的碳位移和条件，不调用HOSE/GNN，不查具体分子的文献答案，不读取原始归属。信号附加代码、H谱、DEPT/APT、HSQC/HMBC及原始谱图未用于本轮。', '',
           '每完成一个分子，就以临时文件写入、flush/fsync、替换、回读及SHA256锁定。五份预测全部锁定后，才读取任何一个数据库答案。比较程序验证每个碳恰好使用一次、每个信号恰好分配一次。', '',
           '这是本次GPT的强制预测表现，不是经验化学家的能力上限，也不能证明这些峰在所有一维证据下都不可区分。小样本、非随机选择、信息不完整，不能外推全库准确率。会话内先锁定后揭盲，不是外部预注册测试，也不能排除通用模型预训练中的分子知识。', '',
           'NMRShiftDB2仅作为比较参照，未复核原论文归属证据；未发现预测程序字段不等于已独立证明实验来源。原始SDF和训练集均未修改，强制预测不作为实验金标准导入标注系统。', '',
           '## 各分子结果','', '| 分子 | 数据库ID | 骨架 | 碳数 | 与库一致 | 比例 | 近邻信号命中 |', '|---|---|---|---:|---:|---:|---:|']
    for r in results:
        lines.append(f"| {r['example_id']} | {r['database_id']} | {titles[r['example_id']]} | {r['total']} | {r['correct']}/{r['total']} | {r['correct']/r['total']:.1%} | {r['near_correct']}/{r['near_total']} |")
    lines += ['', '## 事前把握程度','', '| 事前标记 | 数量 | 与库一致 | 比例 |', '|---|---:|---:|---:|']
    for c,v in confidence.items(): lines.append(f"| {c} | {v['total']} | {v['correct']} | {v['correct']/v['total']:.1%} |")
    lines += ['', '这些等级没有经过概率校准。没有因为揭盲后猜对/猜错而改动事前等级。', '',
              '## 近邻峰逐项对照','', '| 分子 | 位移 ppm | 揭盲前唯一预测 | 数据库 | 结果 |', '|---|---:|---|---|---|']
    exported=[]; checks=[]
    for r in results:
        cid=r['example_id']; b=read(OUT/'blind'/(cid+'.json')); near={p for x,y,g in b['close_pairs'] for p in [x,y]}
        for a in r['rows']:
            if a['ppm'] in near: lines.append(f"| {cid} | {a['ppm']} | C{a['atom']} | C{a['reference_atom']} | {'一致' if a['correct'] else '错误'} |")
    for r in results:
        cid=r['example_id']; b=read(OUT/'blind'/(cid+'.json')); l=read(OUT/'locked'/(cid+'.json'))
        assert sha(OUT/'decisions'/(cid+'.json'))==l['decision_sha256']
        assert sha(OUT/'blind'/(cid+'.json'))==l['blind_sha256']
        assert sorted(a['ppm'] for a in r['rows'])==b['carbon_shifts']
        assert len(set(a['atom'] for a in r['rows']))==b['carbon_count']
        assert len(set(a['reference_atom'] for a in r['rows']))==b['carbon_count']
        mol=Chem.MolFromSmiles(b['mapped_smiles'])
        for a in mol.GetAtoms():
            if a.GetSymbol()=='H': a.SetAtomMapNum(0)
        mol=Chem.RemoveHs(mol)
        rdDepictor.SetPreferCoordGen(True); rdDepictor.Compute2DCoords(mol)
        image_path=OUT/(cid+'_numbered.png'); Draw.MolToFile(mol,str(image_path),size=(1600,1100))
        lines += ['', f"## {cid}：{titles[cid]}",'',
                  f"数据库ID {b['database_id']}；源记录 {b['source_record']}；{b['spectrum_field']}；溶剂 {b['solvent']}；温度 {b['temperature_k']} K。原子编号为原SDF索引+1，图中C:5即原子5，不是第5个峰。",'',
                  f"![编号结构]({image_path.as_posix()})",'',
                  '揭盲前分析：'+l['decision']['summary'],'',
                  f"预测锁定时间：{l['locked_at']}；SHA256：`{l['decision_sha256']}`。",'',
                  '| ppm | 原预测 | 数据库 | 结果 | 事前把握 | 原预测理由（错误理由也保留） |',
                  '|---:|---|---|---|---|---|']
        for a in r['rows']:
            lines.append(f"| {a['ppm']} | C{a['atom']} | C{a['reference_atom']} | {'一致' if a['correct'] else '**错误**'} | {a['confidence']} | {a['reason']} |")
        lines += ['', '**揭盲后分析：** '+diagnoses[cid],'']
        if cid=='CX05': lines += ['原预测措辞说明：48.09的理由中“含氧CH2号原子”指“2号原子上的含氧CH（C2）”，不是说C2为CH₂。仅在报告澄清文字，不改变锁定预测。','']
        exported.append({'schema_version':'forced-blind-carbon-assignment-v1','input':b,'locked_prediction':l,'comparison':r,
                         'label_origin':'GPT_forced_prediction_not_experimental_gold','post_reveal_analysis':diagnoses[cid]})
        checks.append({'example_id':cid,'lock_valid':True,'full_bijection':True})
    assert sha(SOURCE)==read(OUT/'selection.json')['source_sha256']
    (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    with (OUT/'records.jsonl').open('w',encoding='utf-8') as f:
        for e in exported: f.write(json.dumps(e,ensure_ascii=False)+'\n')
        f.flush(); os.fsync(f.fileno())
    summary={'total':total,'correct':correct,'accuracy':correct/total,'near_total':near_total,'near_correct':near_correct,
             'near_accuracy':near_correct/near_total,'confidence':confidence,'source_unchanged':True,'checks':checks,
             'per_molecule':[{k:v for k,v in r.items() if k!='rows'} for r in results]}
    save(OUT/'validation_summary.json',summary)
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('command',choices=['select','lock','reveal','report']); p.add_argument('ids',nargs='*'); args=p.parse_args()
    if args.command=='select': select()
    elif args.command=='lock': lock(args.ids[0])
    elif args.command=='reveal': reveal(args.ids)
    else: report()
