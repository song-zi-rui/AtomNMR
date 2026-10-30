"""Select close 13C signals without disclosing assignments; reveal only after lock."""
import hashlib
import json
import sys
import importlib.util
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
spec=importlib.util.spec_from_file_location('blind_chemistry',ROOT/'nmrpeak/shift_assignment/chemistry.py')
chemistry=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=chemistry
spec.loader.exec_module(chemistry)
carbon_spectra=chemistry.carbon_spectra
indexed_metadata=chemistry.indexed_metadata
OUT=ROOT/'data/blind_close_carbon_examples'
SOURCE=ROOT/'data/nmrshiftdb2/nmrshiftdb2withsignals.sd'

def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

def select():
    if (OUT/'selection.json').exists():
        raise FileExistsError('Preserve the blind selection; use a new output directory for another run.')
    RDLogger.DisableLog('rdApp.*')
    candidates=[]; seen=set()
    for row,mol in enumerate(Chem.SDMolSupplier(str(SOURCE),removeHs=False),1):
        if mol is None: continue
        heavy=Chem.RemoveHs(mol)
        if heavy.GetNumHeavyAtoms()>17 or '.' in Chem.MolToSmiles(heavy): continue
        if any(a.GetSymbol() not in {'C','N','O','S','F','Cl','Br','H'} for a in mol.GetAtoms()): continue
        nc=sum(a.GetSymbol()=='C' for a in mol.GetAtoms())
        if not 6<=nc<=12: continue
        key=Chem.MolToInchiKey(heavy)
        if key in seen: continue
        for sp in carbon_spectra(mol):
            method=indexed_metadata(mol,'Assignment Method',sp.spectrum_id)
            program=indexed_metadata(mol,'Program',sp.spectrum_id)
            if 'predict' in str(method).lower() or program: continue
            shifts=sorted(set(s['shift_ppm'] for s in sp.signals))
            close=[(a,b,round(b-a,6)) for a,b in zip(shifts,shifts[1:]) if 0<b-a<=0.200001]
            if not close or len(shifts)<5: continue
            # Eligibility may inspect labels internally, but no label or label-order leaks to packets.
            assigned={i for s in sp.signals for i in s['atom_indices']}
            if not assigned or any(i>=mol.GetNumAtoms() for i in assigned): continue
            if any(mol.GetAtomWithIdx(i).GetSymbol()!='C' for i in assigned): continue
            if len(assigned)!=nc: continue
            cid=f'EX{len(candidates)+1:02d}'
            mapped=Chem.Mol(mol)
            for a in mapped.GetAtoms(): a.SetAtomMapNum(a.GetIdx()+1)
            atoms=[{'atom':a.GetIdx()+1,'element':a.GetSymbol(),'aromatic':a.GetIsAromatic(),
                    'H':a.GetTotalNumHs(includeNeighbors=True),'neighbors':[[n.GetIdx()+1,str(mol.GetBondBetweenAtoms(a.GetIdx(),n.GetIdx()).GetBondType())] for n in a.GetNeighbors()]} for a in mol.GetAtoms() if a.GetSymbol()!='H']
            blind={'example_id':cid,'source_record':row,'database_id':mol.GetProp('nmrshiftdb2 ID') if mol.HasProp('nmrshiftdb2 ID') else '',
                   'spectrum_field':sp.source_field,'smiles':Chem.MolToSmiles(heavy),'mapped_smiles':Chem.MolToSmiles(mapped),
                   'atoms':atoms,'carbon_shifts':shifts,'close_pairs':close,'solvent':sp.solvent_raw,
                   'temperature_k':sp.temperature_k,'field_strength_mhz':sp.field_strength_mhz,
                   'assignment_method_metadata':method,'program_metadata':program,
                   'blind_protocol':'Only structure, sorted 13C values and conditions. No assignment, model, HOSE DB or other spectra used. Extra signal codes withheld: acquisition origin unverified.'}
            save(OUT/'blind'/(cid+'.json'),blind)
            save(OUT/'hidden_reference'/(cid+'.json'),{'example_id':cid,'signals':list(sp.signals),'source_record':row,'spectrum_field':sp.source_field})
            candidates.append({k:blind[k] for k in ['example_id','source_record','database_id','smiles','carbon_shifts','close_pairs','solvent']})
            seen.add(key); break
        if len(candidates)>=15: break
    save(OUT/'selection.json',{'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'candidates':candidates})
    print(json.dumps(candidates,ensure_ascii=False,indent=2))

def reveal():
    decision=json.loads((OUT/'independent_decisions.json').read_text(encoding='utf-8'))
    if not (OUT/'decision_lock.json').exists():
        save(OUT/'decision_lock.json',{'sha256':hashlib.sha256((OUT/'independent_decisions.json').read_bytes()).hexdigest()})
    lock=json.loads((OUT/'decision_lock.json').read_text())
    assert lock['sha256']==hashlib.sha256((OUT/'independent_decisions.json').read_bytes()).hexdigest()
    comparisons=[]
    for d in decision['examples']:
        cid=d['example_id']; ref=json.loads((OUT/'hidden_reference'/(cid+'.json')).read_text(encoding='utf-8'))
        blind=json.loads((OUT/'blind'/(cid+'.json')).read_text(encoding='utf-8'))
        by_shift={}
        for s in ref['signals']: by_shift.setdefault(str(s['shift_ppm']),[]).extend(i+1 for i in s['atom_indices'])
        results=[]
        for item in d['assignments']:
            gold=sorted(set(by_shift[str(float(item['ppm']))]))
            pred=item['assigned_atoms']
            results.append({**item,'database_atoms':gold,'exact_match':sorted(pred)==gold if pred else None,
                            'candidate_contains_reference':any(sorted(c)==gold for c in item.get('candidate_atom_sets',[]))})
        comparisons.append({'blind':blind,'independent_review':d,'comparison':results})
        mol=Chem.MolFromSmiles(blind['mapped_smiles'])
        mol=Chem.RemoveHs(mol)
        Draw.MolToFile(mol,str(OUT/(cid+'_numbered.png')),size=(950,650))
    save(OUT/'comparison.json',{'locked_decision_sha256':lock['sha256'],'examples':comparisons})
    print(json.dumps([{'example_id':e['blind']['example_id'],'comparisons':e['comparison']} for e in comparisons],ensure_ascii=False,indent=2))

def report():
    result=json.loads((OUT/'comparison.json').read_text(encoding='utf-8'))
    selection=json.loads((OUT/'selection.json').read_text(encoding='utf-8'))
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==selection['source_sha256']
    assert hashlib.sha256((OUT/'independent_decisions.json').read_bytes()).hexdigest()==result['locked_decision_sha256']
    lines=['# NMRShiftDB2 近邻碳峰：三例隐藏归属分析与揭盲对照','',
           '## 方法和限制','',
           '先按 0 < 相邻不同位移差 ≤ 0.20 ppm 筛选，保留中小分子、可解析结构和完整碳归属，获得前15个候选后选择3个教学例子；不是全库计数或随机准确率评估。', '',
           '独立分析只读取结构、排序后的13C位移及条件；不读取数据库归属，不调用HOSE库或用NMRShiftDB2训练的GNN，不查询具体分子的文献答案。未使用H谱、原始谱图、DEPT/APT或二维相关谱。信号附加代码因采集/标记来源未经核验也未作为独立证据。', '',
           '这是同一会话内先锁定判断再揭盲的演示，不是外部预注册盲测，也不能排除通用模型预训练中见过类似知识。', '',
           '原子编号统一使用原SDF索引+1，图中C:13即原子13，不是碳谱第13个峰。相同位移的多条源记录合并为一个信号用于比较。谱仪频率保留原值，不保证该字段在每条源记录中都指碳核观察频率。', '',
           'NMRShiftDB2归属作为本次参照，未重新核查原论文的实验归属证据；与其一致不等于本次独立证明。', '',
           f"独立判断锁定SHA256：`{result['locked_decision_sha256']}`", '',
           '## 汇总','', '| 实例/数据库ID | 谱字段 | 总信号 | 高把握 | 暂定 | 未定 | 已给出归属与库一致 |',
           '|---|---|---:|---:|---:|---:|---:|']
    aggregate={'signals':0,'high':0,'probable':0,'unresolved':0,'assigned_agree':0,'unresolved_candidate_contains':0}
    for e in result['examples']:
        cs=e['comparison']; b=e['blind']; count={c:sum(x['confidence']==c for x in cs) for c in ['high','probable','unresolved']}
        agree=sum(x['exact_match'] is True for x in cs)
        lines.append(f"| {b['example_id']} / {b['database_id']} | {b['spectrum_field']} | {len(cs)} | {count['high']} | {count['probable']} | {count['unresolved']} | {agree}/{len(cs)-count['unresolved']} |")
        aggregate['signals']+=len(cs)
        aggregate['assigned_agree']+=agree
        for k,v in count.items(): aggregate[k]+=v
        aggregate['unresolved_candidate_contains']+=sum(x['confidence']=='unresolved' and x['candidate_contains_reference'] for x in cs)
        assert {float(a['ppm']) for a in cs}==set(b['carbon_shifts'])
        mol=Chem.MolFromSmiles(b['mapped_smiles']); carbons={a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetSymbol()=='C'}
        occupied=set()
        for a in cs:
            assert set(a['database_atoms'])<=carbons
            assert not occupied.intersection(a['assigned_atoms'])
            occupied.update(a['assigned_atoms'])
            for group in [a['assigned_atoms']]+a['candidate_atom_sets']: assert set(group)<=carbons
    lines += ['', '25组信号中，6组高把握、11组暂定、8组不作唯一分配。已分配17组与数据库均一致；8组未定的候选也包含库中答案，但候选命中不能算唯一归属正确。这是68%的选择性分配覆盖率，不是100%全谱准确率。', '']
    titles={'EX11':'对甲基苯甲醛：0.153 ppm 的芳香CH对', 'EX07':'三甲氧基芳香丙烯基化合物：0.18 ppm 的OCH3对', 'EX02':'S-丁基半胱氨酸亚砜骨架：0.20 ppm 的异类碳对'}
    lessons={'EX11':'129.710暂归C2/C3，129.863暂归C1/C4，揭盲后相符。但位移顺序的微小经验倾向并未消除两套配对可能，不能因揭盲正确就回改为事前高把握。',
             'EX07':'数据库为55.97→C11、56.36→C13、56.54→C15。独立分析只确定这三峰属于三个甲氧基，不确定它们的连接位置。即使知道三个芳香C-O的归属，也没有直接建立每个OCH3与其芳环连接位点的实验关联。HMBC等额外证据可能帮助。',
             'EX02':'数据库为50.9→C5（S-CH2-CH(NH2)一侧）、51.1→C9（CH(NH2)）、52.0→C4（丁基一侧S-CH2）。近邻对实际上跨越CH2/CH类别。如果额外测得可靠DEPT/APT，便有可能先识别CH的51.1，而无需靠0.2 ppm的经验排序；本轮没有使用该信息。'}
    for e in result['examples']:
        b=e['blind']; cid=b['example_id']
        lines += [f"## {titles[cid]}",'',f"数据库ID：{b['database_id']}；SDF第{b['source_record']}条；{b['spectrum_field']}；溶剂：{b['solvent']}；温度：{b['temperature_k']} K。",'',
                  f"![原子编号结构]({(OUT/(cid+'_numbered.png')).as_posix()})",'',e['independent_review']['summary'],'',
                  '| ppm | 揭盲前分配 | 事前把握 | 未排除候选 | NMRShiftDB2 | 比较 | 分配依据 |',
                  '|---:|---|---|---|---|---|---|']
        for a in e['comparison']:
            state='一致' if a['exact_match'] is True else ('不一致' if a['exact_match'] is False else '未作唯一分配')
            lines.append(f"| {a['ppm']} | {a['assigned_atoms'] or '未定'} | {a['confidence']} | {a['candidate_atom_sets']} | {a['database_atoms']} | {state} | {a['rationale']} |")
        lines += ['', '揭盲后解读：'+lessons[cid],'']
    lines += ['## 对后续标注的含义','',
              '不以0.2 ppm作为自动合并或放弃归属阈值。区分环境类别判断与精确位置判断，保存整组互换约束；数据库答案只用于对照，不能反向润色已经锁定的独立推理。', '',
              '本次没有修改原始SDF、已有训练数据或标注系统。完整逐峰结果在comparison.json，事前判断在independent_decisions.json；盲包和原归属分别存放。']
    (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    save(OUT/'validation.json',{'source_unchanged':True,'independent_decision_lock_verified':True,'complete_peak_coverage':True,'atom_indices_valid':True,'no_duplicate_assigned_carbon':True,'summary':aggregate})
    print(json.dumps(aggregate,ensure_ascii=False))

if __name__=='__main__':
    {'select':select,'reveal':reveal,'report':report}[sys.argv[1]]()
