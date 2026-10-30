"""Read-only chemical composition and exclusion-policy analysis of ShiftDB2.

The analysis simulates exclusions; it does not change spectra or annotations.
One consistent SQLite snapshot is used, then released before RDKit analysis.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sqlite3
import sys
import time

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from annotation_tool_v3.core import DEFAULT_DB, METAL_ATOMIC_NUMBERS

ROOT = Path(__file__).resolve().parents[1]
FEATURES = {
    'has_n':'含 N','n_ge_2':'至少 2 个 N','n_ge_3':'至少 3 个 N',
    'has_o':'含 O','has_s':'含 S','has_p':'含 P','has_halogen':'含卤素',
    'has_metal':'含金属','charged':'结构中有带电原子','multicomponent':'多片段结构',
    'stereo_specified':'指定了原子或双键立体信息','has_ring':'含环',
    'has_aromatic_ring':'含芳香环','heteroaromatic':'含芳香杂原子',
    'guanidine':'胍基（SMARTS）','amide':'酰胺（SMARTS）',
    'carboxylic_acid':'羧酸／羧酸根（SMARTS）','ester':'酯（SMARTS）',
    'carbonyl':'含 C=O','nitrile':'腈基','nitro':'硝基','sulfonyl':'磺酰基',
    'phenol':'酚羟基','alcohol':'脂肪族醇羟基',
}
PATTERNS = {k:Chem.MolFromSmarts(v) for k,v in {
    'guanidine':'[#6](=[#7])([#7])[#7]', 'amide':'[NX3][CX3](=[OX1])',
    'carboxylic_acid':'[CX3](=[OX1])[OX2H1,OX1-]', 'ester':'[CX3](=[OX1])[OX2][#6]',
    'carbonyl':'[#6]=[#8]', 'nitrile':'[#6]#[#7]', 'nitro':'[N+](=O)[O-]',
    'sulfonyl':'[SX4](=[OX1])(=[OX1])','phenol':'[OX2H1][c]',
    'alcohol':'[OX2H1][CX4]',
}.items()}


def pct(n,d):
    return 100*n/d if d else 0.


def chemical_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {'descriptor_status':'unparseable','numeric':{},'flags':{},'bins':{'elements':'未解析'},'scaffold':None}
    atoms=list(mol.GetAtoms()); elements=Counter(a.GetSymbol() for a in atoms if a.GetAtomicNum()!=1)
    rings=rdMolDescriptors.CalcNumRings(mol); aromatic=rdMolDescriptors.CalcNumAromaticRings(mol)
    heavy=mol.GetNumHeavyAtoms(); mw=Descriptors.MolWt(mol)
    flags={'has_n':elements['N']>0,'n_ge_2':elements['N']>=2,'n_ge_3':elements['N']>=3,
        'has_o':elements['O']>0,'has_s':elements['S']>0,'has_p':elements['P']>0,
        'has_halogen':any(elements[e] for e in ('F','Cl','Br','I')),
        'has_metal':any(a.GetAtomicNum() in METAL_ATOMIC_NUMBERS for a in atoms),
        'charged':any(a.GetFormalCharge()!=0 for a in atoms),'multicomponent':len(Chem.GetMolFrags(mol))>1,
        'stereo_specified':any(a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED for a in atoms)
            or any(b.GetStereo() not in (Chem.BondStereo.STEREONONE,Chem.BondStereo.STEREOANY) for b in mol.GetBonds()),
        'has_ring':rings>0,'has_aromatic_ring':aromatic>0,
        'heteroaromatic':any(a.GetIsAromatic() and a.GetAtomicNum() not in (1,6) for a in atoms)}
    flags.update({k:mol.HasSubstructMatch(p) for k,p in PATTERNS.items()})
    bins={'elements':''.join(e for e in ('C','N','O','S','P','F','Cl','Br','I','B','Si') if elements[e])
              + ''.join(sorted(e for e in elements if e not in ('C','N','O','S','P','F','Cl','Br','I','B','Si'))),
          'heavy_atoms':'1–10' if heavy<=10 else '11–20' if heavy<=20 else '21–30' if heavy<=30 else '31–50' if heavy<=50 else '>50',
          'mw':'<150' if mw<150 else '150–250' if mw<250 else '250–400' if mw<400 else '400–600' if mw<600 else '≥600',
          'rings':'0' if rings==0 else '1' if rings==1 else '2' if rings==2 else '3' if rings==3 else '≥4',
          'nitrogen_count':'0' if not elements['N'] else '1' if elements['N']==1 else '2' if elements['N']==2 else '≥3'}
    # Scaffold identity here deliberately ignores stereo. Clear it on a copy
    # before scaffold extraction, which otherwise can retain stale bond stereo.
    scaffold_mol=Chem.Mol(mol);Chem.RemoveStereochemistry(scaffold_mol)
    try:
        scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=scaffold_mol,includeChirality=False)
        scaffold_error=None
    except Exception as exc:
        scaffold=None;scaffold_error=str(exc)
    return {'descriptor_status':'parsed','numeric':{'mw':mw,'heavy_atoms':heavy,'carbon_atoms':elements['C'],
            'nitrogen_atoms':elements['N'],'rings':rings,'aromatic_rings':aromatic,
            'fraction_csp3':rdMolDescriptors.CalcFractionCSP3(mol)},'flags':flags,'bins':bins,
            'scaffold':scaffold,'scaffold_error':scaffold_error}


def is_flagged(s,mode):
    status=s['status']
    if status=='identity_validated':return False
    return mode=='all_flags' or status!='needs_attention'


def group_metrics(rows, memberships):
    total=len(rows); flagged=sum(r['affected'] for r in rows); removed=sum(r['lost'] for r in rows)
    kept=total-removed; result=[]
    for key,label,indices in memberships:
        group=[rows[i] for i in indices]; n=len(group)
        a=sum(r['affected'] for r in group); lost=sum(r['lost'] for r in group)
        rate=a/n if n else 0; overall=flagged/total if total else 0
        result.append({'key':key,'label':label,'baseline_n':n,'affected_n':a,'fully_lost_n':lost,
            'affected_pct':pct(a,n),'affected_rate_ratio':rate/overall if overall else None,
            'within_type_loss_pct':pct(lost,n),'baseline_share_pct':pct(n,total),
            'flagged_share_pct':pct(a,flagged),'retained_share_pct':pct(n-lost,kept),
            'share_change_pp':pct(n-lost,kept)-pct(n,total),'rare_group':n<100})
    return result


def numeric_summary(rows,key,selector=lambda r:True):
    values=sorted(r['features']['numeric'][key] for r in rows if selector(r) and key in r['features']['numeric'])
    if not values:return {'n':0}
    def quantile(q):
        pos=(len(values)-1)*q; lo=int(pos);hi=min(lo+1,len(values)-1)
        return values[lo]+(values[hi]-values[lo])*(pos-lo)
    return {'n':len(values),'mean':sum(values)/len(values),'p25':quantile(.25),'median':quantile(.5),'p75':quantile(.75)}


def load_snapshot(database):
    conn=sqlite3.connect(f'file:{database.resolve().as_posix()}?mode=ro',uri=True);conn.row_factory=sqlite3.Row
    conn.execute('BEGIN')
    states={(r['scope'],r['target_id']):r['state'] for r in conn.execute('SELECT scope,target_id,state FROM curation_state')}
    total=conn.execute('SELECT count(*) FROM molecules').fetchone()[0]
    query="""SELECT s.id,s.molecule_id,s.nucleus,s.solvent,m.smiles,m.formula,m.mw,
        a.validation_json,coalesce(q.decision,'original_experimental_source') quality,
        s.source_key
        FROM spectra s JOIN molecules m ON m.id=s.molecule_id
        LEFT JOIN spectrum_assignment_audit a ON a.spectrum_id=s.id
        LEFT JOIN spectrum_quality q ON q.spectrum_id=s.id
        WHERE s.source_key>='shiftdb:' AND s.source_key<'shiftdb;' AND s.source='nmrshiftdb2'"""
    spectra=[]
    for row in conn.execute(query):
        r=dict(row);v=json.loads(r.pop('validation_json'))
        r['status']=v['status'];r['issues']=v.get('issue_counts',{})
        r['source_fragment_count']=v.get('fragment_count')
        r.update(source_id=None,record_index=None,field=r['source_key'].split(':',2)[-1])
        if r['status']!='identity_validated':
            raw=conn.execute('SELECT source_json FROM spectra WHERE id=?',(r['id'],)).fetchone()[0]
            source=json.loads(raw)
            r.update(source_id=source.get('nmrshiftdb2_id'),record_index=source.get('source_record_index'))
        r['active']=r['quality']!='quarantined' and states.get(('spectrum',r['id']),'active')=='active' \
            and states.get(('molecule',r['molecule_id']),'active')=='active'
        spectra.append(r)
    conn.close();return spectra,total


def policy_analysis(spectra, features, mode, unit):
    if unit=='spectrum':
        rows=[{'id':s['id'],'affected':is_flagged(s,mode),'lost':is_flagged(s,mode),
               'features':features[s['molecule_id']]} for s in spectra]
    else:
        groups=defaultdict(list)
        for s in spectra:groups[s['molecule_id']].append(s)
        rows=[{'id':mid,'affected':any(is_flagged(s,mode) for s in items),
               'lost':all(is_flagged(s,mode) for s in items),'features':features[mid]} for mid,items in groups.items()]
    memberships=[(key,label,[i for i,r in enumerate(rows) if r['features']['flags'].get(key)]) for key,label in FEATURES.items()]
    table=group_metrics(rows,memberships); categorical={}
    for name in ('elements','heavy_atoms','mw','rings','nitrogen_count'):
        cats=sorted({r['features']['bins'].get(name,'未解析') for r in rows})
        ms=[(v,v,[i for i,r in enumerate(rows) if r['features']['bins'].get(name,'未解析')==v]) for v in cats]
        categorical[name]=group_metrics(rows,ms)
    return {'unit':unit,'mode':mode,'baseline_n':len(rows),'affected_n':sum(r['affected'] for r in rows),
        'fully_lost_n':sum(r['lost'] for r in rows),'retained_n':sum(not r['lost'] for r in rows),
        'affected_pct':pct(sum(r['affected'] for r in rows),len(rows)),
        'fully_lost_pct':pct(sum(r['lost'] for r in rows),len(rows)),
        'features':table,'categorical':categorical,
        'total_variation_by_bin':{k:sum(abs(r['share_change_pp']) for r in v)/200 for k,v in categorical.items()},
        'numeric':{k:{'baseline':numeric_summary(rows,k),'affected':numeric_summary(rows,k,lambda r:r['affected']),
                      'retained':numeric_summary(rows,k,lambda r:not r['lost'])}
                   for k in ('mw','heavy_atoms','carbon_atoms','nitrogen_atoms','rings','fraction_csp3')}}


def paired_availability(spectra,mode):
    groups=defaultdict(list)
    for s in spectra:groups[s['molecule_id']].append(s)
    counts=Counter()
    for items in groups.values():
        before={s['nucleus'] for s in items};kept=[s for s in items if not is_flagged(s,mode)]
        after={s['nucleus'] for s in kept}
        if before=={'13C','1H'}:counts['baseline_both_nuclei_molecules']+=1
        if after=={'13C','1H'}:counts['retained_both_nuclei_molecules']+=1
        def same_solvent(values):
            sets={n:{s['solvent'] for s in values if s['nucleus']==n and s['solvent'] and s['solvent'].lower() not in ('unknown','unreported','unspecified')} for n in ('13C','1H')}
            return bool(sets['13C']&sets['1H'])
        if same_solvent(items):counts['baseline_both_same_known_solvent_molecules']+=1
        if same_solvent(kept):counts['retained_both_same_known_solvent_molecules']+=1
    return dict(counts)


def adjusted_nitrogen(spectra,features,mode):
    groups=defaultdict(list)
    for s in spectra:groups[s['molecule_id']].append(s)
    strata=defaultdict(lambda:[0,0,0,0])
    for mid,items in groups.items():
        f=features[mid]
        if f['descriptor_status']!='parsed':continue
        exposed=f['flags']['has_n'];bad=any(is_flagged(s,mode) for s in items)
        key=(f['bins']['heavy_atoms'], '1' if len(items)==1 else '2' if len(items)==2 else '≥3',
             'C+H' if {s['nucleus'] for s in items}=={'13C','1H'} else items[0]['nucleus'])
        strata[key][0 if exposed and bad else 1 if exposed else 2 if bad else 3]+=1
    numerator=denominator=0.;included=[]
    for key,(a,b,c,d) in strata.items():
        if not (a+b and c+d):continue
        n=a+b+c+d;numerator+=a*d/n;denominator+=b*c/n
        included.append({'stratum':key,'N_flagged':a,'N_not_flagged':b,'noN_flagged':c,'noN_not_flagged':d})
    return {'description':'按重原子数、每分子谱数、核种组合分层的 Mantel–Haenszel 比值比；相关性而非因果',
            'common_odds_ratio':numerator/denominator if denominator else None,'strata':included}


def scaffold_summary(spectra,features,mode):
    molecules=defaultdict(list)
    for s in spectra:molecules[s['molecule_id']].append(s)
    scaffolds=defaultdict(lambda:[0,0,0])
    for mid,items in molecules.items():
        scaffold=features[mid]['scaffold']
        if not scaffold:continue
        v=scaffolds[scaffold];v[0]+=1;v[1]+=any(is_flagged(s,mode) for s in items);v[2]+=all(is_flagged(s,mode) for s in items)
    return {'distinct_nonempty_scaffolds':len(scaffolds),'lost_scaffolds_when_only_spectra_removed':sum(v[2]==v[0] for v in scaffolds.values()),
        'lost_scaffolds_when_any_flagged_molecule_removed':sum(v[1]==v[0] for v in scaffolds.values()),
        'affected_scaffolds':sorted([{'smiles':k,'baseline_molecules':v[0],'affected_molecules':v[1],
            'fully_lost_molecules':v[2]} for k,v in scaffolds.items() if v[1]],key=lambda r:(-r['affected_molecules'],-r['baseline_molecules']))}


def warning_density(spectra,features):
    # Any-warning-per-spectrum rises with the number of reported fragments,
    # even if per-fragment errors are random. Report both denominators.
    groups=defaultdict(lambda:Counter())
    for s in spectra:
        if s['nucleus']!='13C' or s['status'] not in ('identity_validated','needs_attention'):continue
        f=features[s['molecule_id']]
        keys=[('heavy_atoms',f['bins'].get('heavy_atoms','未解析')),('rings',f['bins'].get('rings','未解析')),
              ('nitrogen', '含 N' if f['flags'].get('has_n') else '不含 N')]
        for key in keys:
            group=groups[key];group['spectra']+=1;group['flagged_spectra']+=s['status']=='needs_attention'
            group['warning_occurrences']+=s['issues'].get('carbon_hcount_code_conflict',0)
            if s.get('source_fragment_count') is not None:
                group['reported_fragments']+=s['source_fragment_count'];group['spectra_with_fragment_count']+=1
    return [{'group':key[0],'label':key[1],**dict(v),
             'warning_occurrences_per_100_fragments':100*v['warning_occurrences']/v['reported_fragments'] if v['reported_fragments'] else None}
            for key,v in sorted(groups.items())]


def write_report(folder,result):
    analyses=result['analyses'];strict=analyses['all_flags_molecule'];hard=analyses['hard_only_molecule']
    lines=['# ShiftDB2 归属问题筛除的分布偏向分析','',f"分析时间：{result['created_at']}",'',
        '本分析只读当前 SQLite，模拟筛除，不删除数据、不变更标注。下文“硬问题策略”包括硬冲突、关联 H 暂缓和映射失败；“全部标记策略”再加编码警告。','',
        '## 结论先说清楚','',
        '全删会带来非随机的局部偏向，但常见元素比例与分子量等全局指标变化较小。含 N 本身不是风险因素；金属结构、胍基这一小类，以及较大／多环／含羟基结构需要单独关注。',
        '本次 20 个映射失败分子全部含金属：这是当前结构处理兼容性问题，不能当成 20 个源库归属错误。含胍基分子 145 个，4 个有硬冲突、另外 3 个只有编码警告。两者应分开处置。',
        '建议修正此前“716 条全部筛掉”的默认建议：它可以作为保守对照版本，不宜直接作为唯一正式数据集。优先保留结构／峰、屏蔽不可信归属，并把映射失败与编码警告分别管理。', '',
        '## 数据口径','',f"当前整库有 {result['total_current_molecules']:,} 个分子；主要分析对象是与归属问题直接相关的 ShiftDB2 子集，不用大量未标注 NMRexp 稀释分母。",
        f"ShiftDB2 谱 {result['source_spectra']:,} 条，涉及 {strict['baseline_n']:,} 个唯一分子；状态计数 {result['spectrum_status_counts']}。",
        f"现行隔离／人工停用后活跃谱 {result['active_spectra']:,} 条。谱级重复与分子级统计分开；官能团为程序规则识别，可重叠。",'',
        '## 两种筛选政策的损失','', '| 政策 | 去掉的谱 | 涉及分子 | 只删问题谱后完全失去 ShiftDB2 谱的分子 |', '|---|---:|---:|---:|']
    for mode,label in [('hard_only','仅硬问题／无法映射'),('all_flags','所有标记含警告')]:
        a=analyses[mode+'_molecule'];s=analyses[mode+'_spectrum']
        lines.append(f"| {label} | {s['fully_lost_n']:,} | {a['affected_n']:,} | {a['fully_lost_n']:,}（{a['fully_lost_pct']:.3f}%） |")
    lines+=['','涉及一个问题谱就删整个分子，损失是“涉及分子”；只删该谱，损失是“完全失去谱的分子”。这里不把 NMRexp 替代谱或人工重标救回计入，属于保守模拟。','',
        '## 按分子类型：全部标记策略','', '| 类别 | 原有分子 | 有标记分子 | 该类被标记率 | 相对全体被标记率 | 只删谱后该类完全丢失率 | 筛后组成变化 / 百分点 |', '|---|---:|---:|---:|---:|---:|---:|']
    for r in strict['features']:
        lines.append(f"| {r['label']} | {r['baseline_n']:,} | {r['affected_n']:,} | {r['affected_pct']:.2f}% | {r['affected_rate_ratio']:.2f}× | {r['within_type_loss_pct']:.2f}% | {r['share_change_pp']:+.3f} |")
    lines+=['','相对率不是错误概率，也不是因果。含 N 等类别与官能团交叉；样本少于 100 的组不据倍率单独下结论。SMARTS 官能团识别是一套固定子结构规则，不是化合物类别人工判定。','',
        '## 分子大小、环与元素组成','']
    for name in ('mw','heavy_atoms','rings','nitrogen_count','elements'):
        lines += [f'### {name}','','| 分组 | 原有 | 有标记 | 被标记率 | 完全丢失 | 组成变化 / 百分点 |','|---|---:|---:|---:|---:|---:|']
        for r in sorted(strict['categorical'][name],key=lambda r:-r['baseline_n']):
            lines.append(f"| {r['label']} | {r['baseline_n']:,} | {r['affected_n']:,} | {r['affected_pct']:.2f}% | {r['fully_lost_n']:,} | {r['share_change_pp']:+.3f} |")
    lines += ['','## 按片段数校正警告机会','','大分子有更多信号，只要一个片段出现编码警告就标记整谱，谱级标记率会自然增加。下面同时给出每 100 个报告片段的警告发生次数；片段不是独立统计样本，也不是确认的错误原子。',
              '', '| 分组 | 类别 | C 谱数 | 有警告 C 谱 | 报告片段数 | 警告次数 / 每 100 片段 |', '|---|---|---:|---:|---:|---:|']
    for r in result['warning_density']:
        lines.append(f"| {r['group']} | {r['label']} | {r['spectra']:,} | {r['flagged_spectra']:,} | {r['reported_fragments']:,} | {r['warning_occurrences_per_100_fragments']:.3f} |")
    lines += ['','## 平均与中位数变化','','| 特征 | 原始中位数 | 有标记中位数 | 筛后中位数 | 原始均值 | 筛后均值 |','|---|---:|---:|---:|---:|---:|']
    for key,values in strict['numeric'].items():
        if all(v.get('n') for v in values.values()):lines.append(f"| {key} | {values['baseline']['median']:.3f} | {values['affected']['median']:.3f} | {values['retained']['median']:.3f} | {values['baseline']['mean']:.3f} | {values['retained']['mean']:.3f} |")
    lines+=['','## 双谱、溶剂、环骨架与分层','', '双谱数量仅表示有两个核种，不保证已完整归属、同实验或温度／场强相同；相同溶剂依据当前非 unknown 的规范化标签，未再人工核实原文，也未重新清理所有历史别名。',
        '```json',json.dumps(result['pairing'],ensure_ascii=False,indent=2),'```','',
        f"按重原子大小、谱数和核种组合分层后的含 N 比值比：全标记 {result['adjusted_nitrogen']['all_flags']['common_odds_ratio']}；硬问题 {result['adjusted_nitrogen']['hard_only']['common_odds_ratio']}。仍非因果证明。",'',
        '环骨架统计（不含无环分子的空骨架）：','```json',json.dumps({k:{a:b for a,b in v.items() if a!='affected_scaffolds'} for k,v in result['scaffolds'].items()},ensure_ascii=False,indent=2),'```','',
        '溶剂按谱计数：','', '| 溶剂 | 谱数 | 标记谱 | 标记率 |','|---|---:|---:|---:|']
    for r in result['solvents'][:25]:lines.append(f"| {r['solvent']} | {r['n']:,} | {r['flagged']:,} | {pct(r['flagged'],r['n']):.2f}% |")
    lines+=['','## 解释边界与建议','',
        '- 这是相对本地 ShiftDB2 的筛选偏向，不是与“真实世界分子分布”的符合程度。当前没有定义用途对应的真实世界参照库，不能从此统计宣称具有真实分布代表性。',
        '- 编码警告不是已证明错位。不要因为某类警告富集就说该类型结构必然解析错误。',
        '- 优先保留结构与谱峰，仅关闭不可靠原子归属的监督 mask；可以继续用于没有原子级归属监督的谱到结构训练，但需独立满足谱—结构质量检查。',
        '- 若必须剔除，优先按谱剔除硬问题，不扩大到整个分子；警告另存一层，针对高损失类型优先人工检查。',
        '- 不用复制、过采样未核实标签来补平分布；用人工核实或可靠实验替代谱恢复缺少的结构类型。',
        '- 分类及统计的详细结果为 analysis.json；具体受影响分子与来源谱为 affected_molecules.jsonl。',
        '- 当前修复保留人工标签；此处的盲筛模拟会把人工重标过的问题源谱也算作损失，不能等同于现行训练 mask 的实际损失。']
    (folder/'analysis.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def analyze(database,folder,snapshot=None):
    start=time.monotonic();folder.mkdir(parents=True,exist_ok=True)
    print('Reading consistent, read-only ShiftDB2 snapshot',flush=True)
    if snapshot:
        saved=json.loads(snapshot.read_text(encoding='utf-8'));spectra=saved['spectra'];total=saved['total_current_molecules']
    else:
        spectra,total=load_snapshot(database)
        (folder/'source_snapshot.json').write_text(json.dumps({'total_current_molecules':total,'database':str(database.resolve()),
            'read_at':datetime.now(timezone.utc).isoformat(),'spectra':spectra},ensure_ascii=False),encoding='utf-8')
    print(f'Loaded {len(spectra):,} spectra',flush=True)
    molecules={s['molecule_id']:s for s in spectra};features={};RDLogger.DisableLog('rdApp.warning');RDLogger.DisableLog('rdApp.error')
    for i,(mid,s) in enumerate(molecules.items(),1):
        features[mid]=chemical_features(s['smiles'])
        if i%5000==0:print(f'Calculated descriptors {i:,}/{len(molecules):,}, {time.monotonic()-start:.1f}s',flush=True)
    analyses={mode+'_'+unit:policy_analysis(spectra,features,mode,unit) for mode in ('all_flags','hard_only') for unit in ('spectrum','molecule')}
    active=[s for s in spectra if s['active']]
    solvents=defaultdict(lambda:[0,0,0])
    for s in spectra:
        v=solvents[s['solvent'] or 'unknown'];v[0]+=1;v[1]+=is_flagged(s,'all_flags');v[2]+=is_flagged(s,'hard_only')
    groups=defaultdict(list)
    for s in spectra:groups[s['molecule_id']].append(s)
    with (folder/'affected_molecules.jsonl').open('w',encoding='utf-8') as handle:
        for mid,items in groups.items():
            if any(is_flagged(s,'all_flags') for s in items):
                handle.write(json.dumps({'molecule_id':mid,'smiles':items[0]['smiles'],'formula':items[0]['formula'],
                    'features':features[mid],'spectra':[ {k:s[k] for k in ('id','nucleus','solvent','source_id','record_index','field','status','issues','active')} for s in items]},ensure_ascii=False)+'\n')
    result={'created_at':datetime.now(timezone.utc).isoformat(),'database':str(database.resolve()),'simulation_only':True,
        'total_current_molecules':total,'source_spectra':len(spectra),'source_unique_molecules':len(molecules),
        'active_spectra':len(active),'descriptor_failures':sum(f['descriptor_status']!='parsed' for f in features.values()),
        'spectrum_status_counts':dict(Counter(s['status'] for s in spectra)),'analyses':analyses,
        'active_subset':{mode+'_'+unit:policy_analysis(active,features,mode,unit) for mode in ('all_flags','hard_only') for unit in ('spectrum','molecule')},
        'pairing':{mode:paired_availability(spectra,mode) for mode in ('all_flags','hard_only')},
        'adjusted_nitrogen':{mode:adjusted_nitrogen(spectra,features,mode) for mode in ('all_flags','hard_only')},
        'scaffolds':{mode:scaffold_summary(spectra,features,mode) for mode in ('all_flags','hard_only')},
        'scaffold_descriptor_failures':sum(bool(f.get('scaffold_error')) for f in features.values()),
        'warning_density':warning_density(spectra,features),
        'nucleus_status_counts':{n:dict(Counter(s['status'] for s in spectra if s['nucleus']==n)) for n in ('13C','1H')},
        'solvents':sorted([{'solvent':k,'n':v[0],'flagged':v[1],'hard_problem':v[2]} for k,v in solvents.items()],key=lambda r:-r['n']),
        'elapsed_seconds':round(time.monotonic()-start,2)}
    (folder/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');write_report(folder,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('analyses','active_subset','solvents','scaffolds','adjusted_nitrogen')},ensure_ascii=False,indent=2),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--db',type=Path,default=DEFAULT_DB)
    p.add_argument('--output-dir',type=Path,default=ROOT/'data/shiftdb_atom_identity_audit/exclusion_bias_20261002')
    p.add_argument('--snapshot',type=Path,help='Explicitly reuse a saved read-only snapshot, not current production data')
    args=p.parse_args();analyze(args.db,args.output_dir,args.snapshot)
