"""Audit stored structure predictions against a frozen split, without model calls."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pickle
import sys
from collections import Counter,defaultdict
from functools import lru_cache
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from rdkit import Chem,RDLogger
from rdkit.Chem import Descriptors,rdMolDescriptors

MODE_NAMES={'dual':'H/C 双谱','carbon_only':'仅 C 谱','proton_only':'仅 H 谱'}
ERROR_NAMES={'correct':'严格正确','stereo_difference':'仅立体标记有差异','same_formula_wrong_graph':'同分子式、连接结构不同',
             'wrong_formula':'分子式不符合输入','invalid':'无合法首选候选'}


@lru_cache(maxsize=80000)
def molecule(smiles):
    m=Chem.MolFromSmiles(smiles or '')
    if m is None or not m.GetNumAtoms() or any(a.GetAtomicNum()==0 for a in m.GetAtoms()):return None
    copy=Chem.Mol(m);Chem.RemoveStereochemistry(copy)
    # Remove only stereo labels; isotope and charge distinctions are retained.
    stereo=any(a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED for a in m.GetAtoms()) or any(
        b.GetStereo() not in (Chem.BondStereo.STEREONONE,Chem.BondStereo.STEREOANY) for b in m.GetBonds())
    return {'canonical':Chem.MolToSmiles(m,isomericSmiles=True),'without_stereo':Chem.MolToSmiles(copy,isomericSmiles=True),
            'formula':rdMolDescriptors.CalcMolFormula(m),'heavy_atoms':m.GetNumHeavyAtoms(),
            'mw':Descriptors.MolWt(m),'rings':rdMolDescriptors.CalcNumRings(m),'stereo_specified':stereo}


def buckets(n):
    return '<=10' if n<=10 else '11-20' if n<=20 else '21-30' if n<=30 else '31-40' if n<=40 else '>40'


def proton_feature_coverage(data):
    """Count fields actually recognized by the frozen input tokenizer."""
    import lmdb
    from atompeak_data.schema import finite,j_values
    from atompeak_data.tokenization import Vocabulary
    from atompeak_data.vocab_rules_partial import CATEGORICAL_RULES
    vocabulary=Vocabulary(data/'dict/bart/condition_dict.txt')
    result={}
    for split in ('train','valid','test'):
        groups=defaultdict(Counter)
        with lmdb.open(str(data/'lmdb_dataset'/f'{split}.lmdb'),subdir=False,readonly=True,lock=False) as env:
            with env.begin() as tx:
                for _,blob in tx.cursor():
                    r=pickle.loads(blob)
                    if not r['branch_present']['proton']:continue
                    features=Counter(peaks=len(r['spectra']['proton']['peaks']))
                    for peak in r['spectra']['proton']['peaks']:
                        category=str(peak.get('category') or 'unk').replace(' ','_')
                        nh=finite(peak.get('nH'))
                        features['known_multiplicity']+=category in CATEGORICAL_RULES['hnmr_category'] and 'hnmr_category_'+category in vocabulary.indices
                        features['known_integral']+=nh is not None and nh.is_integer() and 1<=nh<=50 and f'hnmr_nH_{int(nh)}' in vocabulary.indices
                        features['with_J']+=bool(j_values(peak.get('j_values')))
                    mode='dual' if r['branch_present']['carbon'] else 'proton_only'
                    for group in ('all_H',mode):
                        groups[group].update(features);groups[group]['molecules']+=1
                        for key in ('known_multiplicity','known_integral','with_J'):
                            groups[group]['molecules_any_'+key]+=features[key]>0
        result[split]={}
        for label,counts in groups.items():
            values=dict(counts)
            for key in ('known_multiplicity','known_integral','with_J'):
                values[key+'_peak_rate']=counts[key]/counts['peaks'] if counts['peaks'] else 0
            result[split][label]=values
    return result


def audit(predictions,data,split):
    import lmdb
    from atompeak_data.tokenization import Vocabulary,condition_tokens
    with lmdb.open(str(data/'lmdb_dataset'/f'{split}.lmdb'),subdir=False,readonly=True,lock=False) as env:
        with env.begin() as tx:
            records={r['id']:r for _,blob in tx.cursor() if (r:=pickle.loads(blob))}
    source=Vocabulary(data/'dict/bart/condition_dict.txt')
    groups=defaultdict(Counter);errors=Counter();counts=Counter();rank_hist=Counter();valid_hist=Counter();details=[]
    seen=set();string_mismatches=0;flag_mismatches=0;canonical_mismatches=0
    import re
    for row in predictions:
        mid=row['id']
        if mid in seen or mid not in records:raise ValueError('Duplicate or unknown prediction ID: '+mid)
        seen.add(mid);r=records[mid];truth=molecule(r['smiles_canonical']);reported=molecule(row.get('target_for_evaluation_only'))
        if truth is None or reported is None or reported['canonical']!=truth['canonical']:
            raise ValueError('Prediction target differs from frozen split: '+mid)
        if r['molecular_formula']!=truth['formula']:raise ValueError('Stored input formula differs from structure: '+mid)
        if row.get('input_mode')!='as_stored':raise ValueError('This analysis expects as_stored input; masked-input experiments require separate group definitions')
        mode='dual' if all(r['branch_present'].values()) else 'carbon_only' if r['branch_present']['carbon'] else 'proton_only'
        candidates=[molecule(b.get('smiles')) for b in row['beams']]
        for saved,parsed in zip(row['beams'],candidates):
            flag_mismatches+=bool(saved.get('valid'))!=bool(parsed)
            canonical_mismatches+=bool(parsed and saved.get('canonical_smiles')!=parsed['canonical'])
        exact=[i+1 for i,b in enumerate(candidates) if b and b['canonical']==truth['canonical']]
        no_stereo=[i+1 for i,b in enumerate(candidates) if b and b['without_stereo']==truth['without_stereo']]
        first=candidates[0] if candidates else None
        error=('correct' if first and first['canonical']==truth['canonical'] else 'invalid' if first is None else
               'stereo_difference' if first['without_stereo']==truth['without_stereo'] else
               'wrong_formula' if first['formula']!=truth['formula'] else 'same_formula_wrong_graph')
        errors[error]+=1;rank_hist[str(min(exact) if exact else 0)]+=1
        valid_unique=len({b['canonical'] for b in candidates if b});valid_hist[str(valid_unique)]+=1
        formula_candidates=[b for b in candidates if b and b['formula']==r['molecular_formula']]
        complete=all(r['alignment_supervision'][b]['label_complete_by_branch'] for b in ('carbon','proton') if r['branch_present'][b])
        tokens=[]
        for b in ('carbon','proton'):
            if r['branch_present'][b]:tokens.append(condition_tokens(r,b,True,False,False)[1])
        specific=any(not t.endswith(('_missing','_other')) for t in tokens)
        reported_solvent=any(not t.endswith('_missing') for t in tokens)
        solvent='有具体溶剂 token' if specific else '仅 other 溶剂 token' if reported_solvent else '所有可用谱均缺溶剂'
        formula_tokens=[]
        for element,count in re.findall(r'([A-Z][a-z]?|[+-])(\d*)',r['molecular_formula']):
            formula_tokens.append('formula_charge_'+element if element in '+-' else 'formula_ele_'+element)
            if count and int(count)>1:formula_tokens.append('formula_cnt_'+count)
        unknown_formula=[t for t in formula_tokens if t not in source.indices]
        names=[('overall','all'),('mode',mode),('heavy_atoms',buckets(truth['heavy_atoms'])),
               ('coverage','完整确定标签' if complete else '部分或弱标签'),('solvent',solvent),
               ('origin',str(r.get('origin'))),('target_stereo','原结构有立体标记' if truth['stereo_specified'] else '原结构无立体标记'),
               ('mode_and_size',mode+'|'+buckets(truth['heavy_atoms'])),('mode_and_solvent',mode+'|'+solvent),
               ('formula_token_support','含 UNK 分子式 token' if unknown_formula else '分子式 token 齐全')]
        for kind,label in names:
            g=groups[(kind,label)];g['n']+=1;g['sum_heavy_atoms']+=truth['heavy_atoms']
            g['sum_valid_unique']+=valid_unique;g['sum_emitted']+=len(candidates)
            for k in (1,3,10):
                g['top'+str(k)]+=any(i<=k for i in exact)
                g['no_stereo_top'+str(k)]+=any(i<=k for i in no_stereo)
                g['formula_filtered_top'+str(k)]+=any(b['canonical']==truth['canonical'] for b in formula_candidates[:k])
            g['formula_top1']+=bool(first and first['formula']==r['molecular_formula'])
            g['formula_any']+=bool(formula_candidates);g['valid_top1']+=bool(first)
        details.append({'id':mid,'mode':mode,'formula':r['molecular_formula'],'target':r['smiles_canonical'],
            'top1_smiles':row['beams'][0]['smiles'] if row['beams'] else None,
            'top1_formula':first['formula'] if first else None,'target_heavy_atoms':truth['heavy_atoms'],
            'predicted_heavy_atoms':first['heavy_atoms'] if first else None,'target_mw':truth['mw'],
            'error_category':error,'strict_hit_rank':min(exact) if exact else None,'without_stereo_hit_rank':min(no_stereo) if no_stereo else None,
            'emitted_candidates':len(candidates),'valid_unique_candidates':valid_unique,'formula_compatible_candidates':len(formula_candidates),
            'formula_filtered_top1_smiles':formula_candidates[0]['canonical'] if formula_candidates else None,
            'formula_filtered_top1_correct':bool(formula_candidates and formula_candidates[0]['canonical']==truth['canonical']),
            'label_coverage':'complete' if complete else 'partial_or_weak','solvent_group':solvent,'solvent_tokens':tokens,
            'unknown_formula_tokens':unknown_formula})
    grouped=[]
    for (kind,label),values in groups.items():
        d={'kind':kind,'label':label,**dict(values)}
        for key in ('top1','top3','top10','no_stereo_top1','no_stereo_top3','no_stereo_top10',
                    'formula_filtered_top1','formula_filtered_top3','formula_filtered_top10','formula_top1','formula_any','valid_top1'):
            d[key+'_rate']=d[key]/d['n']
        d['mean_heavy_atoms']=d['sum_heavy_atoms']/d['n'];d['mean_valid_unique']=d['sum_valid_unique']/d['n'];d['mean_emitted']=d['sum_emitted']/d['n']
        grouped.append(d)
    return {'processed':len(predictions),'split_records':len(records),'missing_split_ids':sorted(set(records)-seen),
        'unknown_or_duplicate_ids':0,'valid_flag_mismatches':flag_mismatches,'stored_canonical_string_mismatches':canonical_mismatches,
        'groups':grouped,'error_categories':dict(errors),'valid_unique_histogram':dict(sorted(valid_hist.items(),key=lambda x:int(x[0]))),
        'strict_hit_rank_histogram':dict(sorted(rank_hist.items(),key=lambda x:int(x[0]))),
        'note':'Formula filtering is a post-hoc diagnostic using an input field, not a newly validated model or an untouched-test improvement claim.'},details


def plot(summary,out):
    os.environ.setdefault('MPLCONFIGDIR',str(out/'matplotlib_cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    font=FontProperties(fname='C:/Windows/Fonts/msyh.ttc') if Path('C:/Windows/Fonts/msyh.ttc').exists() else FontProperties()
    plt.rcParams.update({'svg.fonttype':'none','axes.unicode_minus':False,'font.size':10})
    def group(kind,label):return next(x for x in summary['groups'] if x['kind']==kind and x['label']==label)
    fig,axes=plt.subplots(2,2,figsize=(13.4,9.6),layout='constrained')
    modes=['dual','carbon_only','proton_only'];ax=axes[0,0]
    for i,k in enumerate((1,3,10)):
        bars=ax.bar([j+(i-1)*.24 for j in range(3)],[100*group('mode',m)['top'+str(k)+'_rate'] for m in modes],width=.23,label='Top-'+str(k),color=['#196f72','#67a7a0','#d69b67'][i])
        ax.bar_label(bars,fmt='%.1f',padding=3,fontsize=9)
    ax.set_xticks(range(3),[MODE_NAMES[m]+'\nn='+str(group('mode',m)['n']) for m in modes],fontproperties=font)
    ax.set_ylabel('严格结构命中率 (%)',fontproperties=font);ax.set_title('双谱与单谱的表现差异',fontproperties=font,fontsize=13)
    ax.set_ylim(0,53);ax.legend(frameon=False)
    ax=axes[0,1];order=['wrong_formula','same_formula_wrong_graph','correct','invalid','stereo_difference']
    vals=[summary['error_categories'].get(k,0) for k in order]
    bars=ax.barh(range(len(order)),[100*v/summary['processed'] for v in vals],color=['#b76750','#d49b6b','#388680','#939e9b','#8598b6'])
    ax.set_yticks(range(len(order)),[ERROR_NAMES[k] for k in order],fontproperties=font);ax.invert_yaxis()
    ax.bar_label(bars,labels=[f'{v:,} ({100*v/summary["processed"]:.2f}%)' for v in vals],padding=4,fontsize=9)
    ax.set_xlim(0,88);ax.set_title('首选结果的互斥错误分类',fontproperties=font,fontsize=13)
    ax=axes[1,0];bins=['<=10','11-20','21-30','31-40','>40']
    for k,c in [(1,'#196f72'),(10,'#d69b67')]:
        ax.plot(range(5),[100*group('heavy_atoms',b)['top'+str(k)+'_rate'] for b in bins],marker='o',label='Top-'+str(k),color=c,lw=2)
    ax.set_xticks(range(5),[b+'\nn='+str(group('heavy_atoms',b)['n']) for b in bins]);ax.set_xlabel('真实结构重原子数',fontproperties=font)
    ax.set_ylabel('严格结构命中率 (%)',fontproperties=font);ax.set_title('分子变大后，命中率快速下降',fontproperties=font,fontsize=13);ax.legend(frameon=False)
    ax=axes[1,1];hist=summary['valid_unique_histogram']
    ax.bar(range(11),[hist.get(str(i),0) for i in range(11)],color='#568e92')
    ax.set_xticks(range(11));ax.set_xlabel('每个分子的合法不同候选数',fontproperties=font);ax.set_ylabel('测试分子数',fontproperties=font)
    mean=group('overall','all')['mean_valid_unique']
    ax.set_title(f'请求 beam=10，实际平均 {mean:.2f} 个合法不同候选',fontproperties=font,fontsize=13)
    for ax in axes.flat:ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    fig.suptitle('AtomPeak 测试集逐样本分析（4,277 个分子）',fontproperties=font,fontsize=17)
    fig.savefig(out/'prediction_diagnostics.svg');fig.savefig(out/'prediction_diagnostics.png',dpi=170);plt.close(fig)


def make_report(summary,details,out):
    def g(kind,label):return next(x for x in summary['groups'] if x['kind']==kind and x['label']==label)
    all_=g('overall','all');n=all_['n']
    lines=['# AtomPeak 测试结果逐分子分析','', '日期：2026-10-26。读取当前上传的 test_best_generation_v1.jsonl、汇总 JSON 和冻结测试集，不重新训练或调用模型。',
        '', '最突出的现象是：双谱优于单谱，H 单谱很弱；大分子生成明显失败；模型虽然输入了分子式，首选结果仍大量违反分子式。只忽略立体标记，命中率变化很小，主要问题是组成与连接结构。',
        '', '## 1. 数据和汇总核对','',
        f'- {n:,} 条预测与测试集逐 ID 对应，没有重复、缺失或未知 ID；评估目标与冻结结构逐条一致。',
        f'- 重新解析候选：valid 标志不一致 {summary["valid_flag_mismatches"]} 条，存储规范 SMILES 字符串与本地计算不一致 {summary["stored_canonical_string_mismatches"]} 条。',
        '- 独立复算 Top-1/3/10 与上传汇总完全一致：314 / 542 / 749。',
        '- 当前代码先把合法结构排在非法结构前。因此 valid_top1=4160/4277 实际是“至少生成一个合法候选”的覆盖情况，不是原始 beam 首选合法率，更不是正确率。',
        '- 这些都是结构生成评价，没有从该文件计算出原子归属准确率。文件名为 best，但本身没有 checkpoint 哈希；正式归档时应保留运行命令及断点身份。',
        '', '![预测诊断](prediction_diagnostics.png)',
        '', '## 2. 单谱与双谱拆分','', '| 输入 | 分子数 | Top-1 | Top-3 | Top-10 | 首选分子式符合 |','|---|---:|---:|---:|---:|---:|']
    for mode in ('dual','carbon_only','proton_only'):
        x=g('mode',mode);lines.append(f'| {MODE_NAMES[mode]} | {x["n"]:,} | {x["top1_rate"]:.2%} | {x["top3_rate"]:.2%} | {x["top10_rate"]:.2%} | {x["formula_top1_rate"]:.2%} |')
    lines += ['', '整体 Top-1=7.34% 不能作为双谱性能。双谱子集 Top-1=19.08%、Top-10=42.60%，而 H 单谱 Top-1 只有 2.30%。这些是不同测试分子的分组结果；结构大小、化学类型和来源可能不同，不能把差异直接当作“给同一个分子额外输入 H/C 后提升多少”的因果结论。']
    if 'proton_feature_coverage' in summary:
        lines += ['', '### H 谱实际提供了多少附加信息','',
            '以下统计冻结 LMDB 中 tokenizer 实际读取的 spectra.proton.peaks。多重性必须被规则和词表识别；积分必须可编码为 1–50 的整数；J 值必须通过当前数值解析规则。各比例分母都是相应组的 H 峰数，不是分子数。',
            '', '| split / H 输入组 | 分子数 | H 峰数 | 有效多重性峰 | 有效积分峰 | 含 J 值峰 |','|---|---:|---:|---:|---:|---:|']
        for split in ('train','valid','test'):
            for mode in ('all_H','dual','proton_only'):
                x=summary['proton_feature_coverage'][split].get(mode)
                if x:
                    cells=[f'{x[key]:,} ({x[key+"_peak_rate"]:.2%})' for key in ('known_multiplicity','known_integral','with_J')]
                    lines.append(f'| {split} / {"全部 H" if mode=="all_H" else MODE_NAMES[mode]} | {x["molecules"]:,} | {x["peaks"]:,} | '+ ' | '.join(cells)+' |')
        lines += ['', '开启 H 谱多重性、积分、J 值输入开关，不等于每条 H 谱都具备这些信息。缺失多重性和积分进入 hnmr_unk；缺失 J 值不产生 J token。不能把峰的原子归属数量直接当作实验积分来填补：归属并不保证谱完整，也不保证积分已测量。',
            '', '本轮还按 selected_spectrum_ids 回查当前 SQLite 中 3 条训练 H 单谱（M001915436b8562d64f970a6b、M00223d9dde39582b63fd068e、M003056ca609be036361bd5bf）：来源均为 nmrshiftdb2，其所有峰的 multiplicity、reported_nH、j_values 均为空。observed_peak 会直接读取这三个来源字段。这 3 条没有显示出“来源已有这些值、导出时丢失”的迹象；该回查是定位样例，不代表已审计全部原始 SDF。',
            '', 'H 单谱结果还受到分子大小差异影响：其平均重原子数为 24.61，双谱为 14.42，C 单谱为 15.16。该审计提示应优先补充实验信息完整的 H 谱并单独检查 H 输入质量，不能仅据此确定某一个因素造成了全部性能差距。']
    lines += ['', '## 3. 错在哪里','', '| 首选结果类别 | 数量 | 占全部测试分子 |','|---|---:|---:|']
    for key in ('correct','wrong_formula','same_formula_wrong_graph','stereo_difference','invalid'):
        count=summary['error_categories'].get(key,0);lines.append(f'| {ERROR_NAMES[key]} | {count:,} | {count/n:.2%} |')
    lines += ['', '“仅立体标记有差异”包含漏标、多标或相反标记。分析仅删除立体标记，保留元素、同位素和电荷。忽略立体标记后 Top-1 从 314 增至 331（7.74%），Top-10 从 749 增至 775（18.12%）；不能用立体评估过严解释大多数失败。',
        '', f'首选分子式符合的只有 {all_["formula_top1"]:,}/{n:,}（{all_["formula_top1_rate"]:.2%}）。有任意分子式相符候选的样本是 {all_["formula_any"]:,}（{all_["formula_any_rate"]:.2%}）。分子式来自模型已知输入，不是额外泄漏真实结构。',
        '', '分子式信息当前通过 token 提供，没有在解码或候选选择时硬性检查。这说明本轮模型没有稳定遵守这个输入约束；不等于分子式没有被传入，也不能据此断言全部问题都是训练实现错误。',
        '', '### 只做输入分子式过滤的诊断','',
        '离线保留与已知输入分子式相符的候选，保持原候选顺序，不调用模型、不生成新结构。无匹配候选的分子仍留在分母内，计为未命中。',
        '', '| 输入 | 原 Top-1 | 分子式过滤后 Top-1 |','|---|---:|---:|']
    for mode in ('all','dual','carbon_only','proton_only'):
        x=all_ if mode=='all' else g('mode',mode)
        lines.append(f'| {"全部" if mode=="all" else MODE_NAMES[mode]} | {x["top1_rate"]:.2%} | {x["formula_filtered_top1_rate"]:.2%} |')
    lines += ['', '全部从 314 增至 435（10.17%），双谱从 129 增至 182（26.92%）。这是测试集上的事后诊断，提示利用已知输入约束有实际空间；不是新训练模型的独立验证提升。本次没有修改生产推理代码。今后应在 valid 固定规则，再作正式评估。',
        '', '## 4. 分子大小','', '| 重原子数 | 分子数 | Top-1 | Top-10 | 首选分子式符合 |','|---|---:|---:|---:|---:|']
    for b in ('<=10','11-20','21-30','31-40','>40'):
        x=g('heavy_atoms',b);lines.append(f'| {b} | {x["n"]:,} | {x["top1_rate"]:.2%} | {x["top10_rate"]:.2%} | {x["formula_top1_rate"]:.2%} |')
    lines += ['', '超过 20 个重原子的 1,097 个测试分子中，Top-1 仅 11 个，Top-10 仅 24 个。大分子的组成和连接组合更难，这是当前生成能力明显不足的区域。该文件没有完整 token ID、EOS 或执行耗时，不能仅凭候选字符串进一步确认是否触及最大生成长度。',
        '', '### 按模态和大小交叉检查','', '| 输入 | 重原子数 | n | Top-1 | Top-10 |','|---|---|---:|---:|---:|']
    for mode in ('dual','carbon_only','proton_only'):
        for b in ('<=10','11-20','21-30','31-40','>40'):
            xs=[x for x in summary['groups'] if x['kind']=='mode_and_size' and x['label']==mode+'|'+b]
            if xs:
                x=xs[0];lines.append(f'| {MODE_NAMES[mode]} | {b} | {x["n"]} | {x["top1_rate"]:.2%} | {x["top10_rate"]:.2%} |')
    lines += ['', '## 5. 候选池与排序','',
        f'- 平均返回 {all_["mean_emitted"]:.2f} 个候选，其中平均 {all_["mean_valid_unique"]:.2f} 个合法不同结构；beam=10 并不保证留下 10 个不同合法结构。',
        '- 只有 1,298 个分子留下 10 个合法不同候选，117 个没有任何合法候选。',
        '- 3,528 个分子的正确结构没有出现在当前候选池，固定池内重排无法补救。',
        '- 749 个候选池含正确答案，其中 435 个没排首位。即使固定池内排序完美，Top-1 理论上限也只是 17.51%，不是预期可达到值。',
        '- 扩大 beam 可作为后续验证集实验，但会增加计算量，也不保证候选多样性、Top-1 或正确结构覆盖一定改善。',
        '', '## 6. 部分标签、溶剂与来源','', '| 测试样本标签覆盖 | n | Top-1 | Top-10 |','|---|---:|---:|---:|']
    for label in ('完整确定标签','部分或弱标签'):
        x=g('coverage',label);lines.append(f'| {label} | {x["n"]} | {x["top1_rate"]:.2%} | {x["top10_rate"]:.2%} |')
    lines += ['', '该分组按测试记录的标签覆盖划分。推理没有读取这些标签。完整与部分测试子集的 Top-10 相近，不能据此认为部分归属训练造成性能低，也不能证明加入部分归属训练有益；训练数据贡献仍需固定数据的消融。',
        '', '| 模型可见溶剂类别 | n | 平均重原子数 | Top-1 | Top-10 |','|---|---:|---:|---:|---:|']
    for x in [x for x in summary['groups'] if x['kind']=='solvent']:
        lines.append(f'| {x["label"]} | {x["n"]} | {x["mean_heavy_atoms"]:.1f} | {x["top1_rate"]:.2%} | {x["top10_rate"]:.2%} |')
    lines += ['', '具体溶剂按实际 token 是否属于已识别类别判定；other 单独列出。溶剂已知/缺失组不是随机对照，化学空间、模态和来源可能不同，因此不能由此推断溶剂 token 的因果增益。',
        '', '| 原分子来源字段 | n | Top-1 | Top-10 |','|---|---:|---:|---:|']
    for x in [x for x in summary['groups'] if x['kind']=='origin']:
        lines.append(f'| {x["label"]} | {x["n"]} | {x["top1_rate"]:.2%} | {x["top10_rate"]:.2%} |')
    lines += ['', 'origin 是分子来源，不是逐峰人工/原库标签来源；NMRexp 子集仅 15 个样本，不适合据此作稳定的泛化结论。',
        '', '## 7. 下一步优先级','',
        '1. 把这轮结果完整保存，主指标按整体和三种模态分别报告；不要把 7.34% 当作双谱性能，也不要把 97.26% 当作准确率。',
        '2. 先在 valid 检查分子式约束、候选多样性，以及 10/20/50 beam 的覆盖与开销，再决定是否修改推理。当前诊断不改变正式测试原始结果。',
        '3. 检查 H 单谱：区分完整实验信息与类别/积分/J 缺失；比较同一双谱遮蔽另一支后的预测，避免不同分子之间作因果比较。',
        '4. 对生成候选内已有答案的样本分析排序；UOT 权重仍只在 valid 上定。组成约束与排序都不能替代缺失正确结构的候选生成。',
        '5. 再开展单谱/部分标签/溶剂消融，并单独保留小分子与大分子成绩。数据扩充应覆盖 H 谱信息及较大分子的化学环境，不能只重复增加易识别小分子。',
        '', '## 8. 可追踪例子（并非随机代表性抽样）','']
    for key in ('wrong_formula','same_formula_wrong_graph','stereo_difference','invalid','correct'):
        cases=sorted((x for x in details if x['error_category']==key),key=lambda x:(x['target_heavy_atoms'],x['id']))[:3]
        lines += ['### '+ERROR_NAMES[key],'']
        for x in cases:
            lines += [f'- `{x["id"]}`，{MODE_NAMES[x["mode"]]}，输入分子式 `{x["formula"]}`。',
                      f'  - 真实结构：`{x["target"]}`',f'  - 首选结构：`{x["top1_smiles"]}`；其分子式：`{x["top1_formula"]}`。']
        lines.append('')
    lines += ['## 9. 产物与边界','',
        '- summary.json：全部分组、命中率和错误分类；per_molecule.jsonl：逐分子的模式、错误类别、公式相容性、大小和命中位置。',
        '- prediction_diagnostics.svg：矢量图；prediction_diagnostics.png：预览。',
        '- 这是对已上传测试输出的分析，没有重新训练、重跑推理或覆写测试文件；没有据测试集挑选模型/权重。',
        '- 原始 beam 第一名和被规范化去重掉的条目没有保留，不能从当前文件还原其合法性或原始候选重复率。',
        '- 排名率按当前输出顺序和严格规范 SMILES 定义复算；仅去立体标记的结果是额外诊断，不替换严格成绩。']
    (out/'prediction_analysis.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--predictions',type=Path,required=True);p.add_argument('--data',type=Path,required=True)
    p.add_argument('--split',choices=('valid','test'),default='test');p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise ValueError('Choose a new analysis folder')
    RDLogger.DisableLog('rdApp.error');RDLogger.DisableLog('rdApp.warning')
    predictions=[json.loads(s) for s in args.predictions.read_text(encoding='utf-8').splitlines() if s.strip()]
    summary,details=audit(predictions,args.data,args.split)
    summary['proton_feature_coverage']=proton_feature_coverage(args.data)
    overall=next(x for x in summary['groups'] if x['kind']=='overall')
    metrics_path=args.predictions.with_suffix('.metrics.json')
    if metrics_path.exists():
        original=json.loads(metrics_path.read_text(encoding='utf-8'))
        for k in ('processed','top1','top3','top10','valid_top1'):
            value=summary['processed'] if k=='processed' else overall[k]
            if original[k]!=value:raise ValueError('Recomputed metric differs: '+k)
        summary['original_metrics']=original
    summary.update(prediction_file=str(args.predictions.resolve()),data=str(args.data.resolve()),
                   predictions_sha256=hashlib.sha256(args.predictions.read_bytes()).hexdigest(),
                   split=args.split,rdkit_version=__import__('rdkit').__version__)
    args.output.mkdir(parents=True)
    (args.output/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    with (args.output/'per_molecule.jsonl').open('w',encoding='utf-8') as f:
        for row in details:f.write(json.dumps(row,ensure_ascii=False)+'\n')
    plot(summary,args.output);make_report(summary,details,args.output)
    print(json.dumps({'report':str((args.output/'prediction_analysis.md').resolve()),'overall':overall,
                      'by_mode':[x for x in summary['groups'] if x['kind']=='mode'],'errors':summary['error_categories']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
