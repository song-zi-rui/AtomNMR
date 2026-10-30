"""Audit the uploaded 2026-10-26 run; the narrative is specific to this snapshot.

Checkpoint, metrics and dataset inputs are read-only. A separate output folder
contains the report and figures; this is not a general-purpose training evaluator.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import heapq
import json
import math
import os
import pickle
import sys
from collections import Counter, defaultdict
from contextlib import ExitStack
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def audit_dataset(folder):
    import lmdb
    from atompeak_data.tokenization import condition_tokens
    manifest = json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
    counts, connectivity = Counter(), defaultdict(set)
    for line in (folder/'records_index.jsonl').read_text(encoding='utf-8').splitlines():
        r = json.loads(line)
        counts[(r['split'],r['mode'],r['has_supervision'])] += 1
        connectivity[r['connectivity_key']].add(r['split'])
    split_stats = {}
    fingerprint = hashlib.sha256()
    # The exporter hashes id+original-pickle in global molecule-ID order.
    with ExitStack() as stack:
        cursors, queue, previous = {}, [], {}
        for split in ('train','valid','test'):
            split_stats[split] = {'records':0,'modes':Counter(),'origins':Counter(),
                'branches':{b:Counter() for b in ('carbon','proton')},
                'solvent_tokens':{b:Counter() for b in ('carbon','proton')}}
            env = stack.enter_context(lmdb.open(str(folder/'lmdb_dataset'/f'{split}.lmdb'),
                subdir=False,readonly=True,lock=False,readahead=True))
            txn = stack.enter_context(env.begin())
            cursors[split] = iter(txn.cursor())
            value = next(cursors[split],None)
            if value:
                r=pickle.loads(value[1]);heapq.heappush(queue,(r['id'],split,value[1],r))
        while queue:
            mid,split,blob,r = heapq.heappop(queue)
            if split in previous and mid <= previous[split]:
                raise ValueError('LMDB order is not the original ID traversal')
            previous[split] = mid
            fingerprint.update(mid.encode()+blob)
            t=split_stats[split];t['records']+=1;t['origins'][str(r.get('origin'))]+=1
            mode='dual' if all(r['branch_present'].values()) else 'carbon_only' if r['branch_present']['carbon'] else 'proton_only'
            t['modes'][mode]+=1
            for b in ('carbon','proton'):
                if not r['branch_present'][b]:continue
                labels=r['alignment_supervision'][b];stats=t['branches'][b]
                stats['records']+=1;stats['peaks']+=len(r['spectra'][b]['peaks'])
                stats['hard_peaks']+=sum(any(row) for row in labels['pair_known_mask'])
                stats['weak_peaks']+=sum(labels['allowed_row_mask'])
                stats['shift_targets']+=sum(labels['shift_known_mask'])
                stats['sets']+=len(labels['set_peak_membership'])
                stats['mass_valid_records']+=bool(labels['mass_supervision_valid'])
                stats['spectrum_complete_records']+=bool(labels.get('spectrum_complete'))
                stats['label_complete_records']+=bool(labels.get('label_complete_by_branch'))
                token=condition_tokens(r,b,True,False,False)[1]
                t['solvent_tokens'][b][token]+=1
            value=next(cursors[split],None)
            if value:
                nxt=pickle.loads(value[1]);heapq.heappush(queue,(nxt['id'],split,value[1],nxt))
    actual = fingerprint.hexdigest()
    return {'manifest':manifest,'splits':split_stats,'connectivity_groups':len(connectivity),
            'cross_split_connectivity_groups':sum(len(v)>1 for v in connectivity.values()),
            'index_counts':[{'split':s,'mode':m,'has_supervision':label,'count':n} for (s,m,label),n in sorted(counts.items())],
            'all_index_rows_have_supervision':all(label for s,m,label in counts),
            'recomputed_payload_fingerprint':actual,
            'payload_fingerprint_matches_manifest':actual==manifest['source_snapshot_fingerprint'],
            'vocabulary_checksums':{name:digest(folder/'dict/bart'/name)==value for name,value in manifest['vocab_sha256'].items()}}


def checkpoint_metadata(folder):
    import torch
    torch.set_num_threads(1)
    output={}
    for name in ('checkpoint_best.pt','checkpoint_prev.pt','checkpoint_last.pt'):
        path=folder/name
        if not path.exists():continue
        # Trusted project-generated artifacts; mmap avoids allocating optimizer tensor copies.
        saved=torch.load(path,map_location='cpu',mmap=True,weights_only=False)
        keys=('checkpoint_version','epoch','batches_completed','epoch_completed','step','best','dataset_fingerprint',
              'model_config','training_signature','checkpoint_reason','tokenization','initialization')
        data={k:saved.get(k) for k in keys}
        data.update(file_bytes=path.stat().st_size,epoch_display=saved['epoch']+1,
                    optimizer_lrs=[g['lr'] for g in saved.get('optimizer',{}).get('param_groups',[])],
                    scheduler=saved.get('scheduler'),rng_world_size=len(saved.get('rng',[])),
                    source_vocab_size=len(saved['source_tokens']),target_vocab_size=len(saved['target_tokens']),
                    model_state_keys=len(saved['model']))
        output[name]=data
        del saved;gc.collect()
    return {'reader_torch_version':torch.__version__,'checkpoints':output}


def plots(rows,out):
    os.environ.setdefault('MPLCONFIGDIR',str(out/'matplotlib_cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    font_path=Path('C:/Windows/Fonts/msyh.ttc')
    prop=FontProperties(fname=str(font_path)) if font_path.exists() else FontProperties()
    plt.rcParams.update({'svg.fonttype':'none','axes.unicode_minus':False,'font.size':10})
    epoch=[r['epoch']+1 for r in rows]
    best=min(rows,key=lambda r:r['validation']['generation_loss'])['epoch']+1
    fig,axes=plt.subplots(2,2,figsize=(13.4,8.6),layout='constrained')
    specs=[('验证集生成损失与总损失',[('generation_loss','生成交叉熵','#117b78'),('loss','总损失','#879299')],'损失'),
           ('Token 准确率（给定真实前文）',[('token_accuracy','Token 准确率','#3a70a3')],'准确率 (%)'),
           ('确定归属的配对损失',[('carbon_assignment','C 配对 BCE','#117b78'),('proton_assignment','H 配对 BCE','#d77c49')],'损失'),
           ('归属分支综合损失',[('carbon_alignment_total','C 综合归属损失','#117b78'),('proton_alignment_total','H 综合归属损失','#d77c49')],'损失')]
    for ax,(title,lines,ylabel) in zip(axes.flat,specs):
        for key,label,color in lines:
            y=[r['validation'][key]*(100 if key=='token_accuracy' else 1) for r in rows]
            ax.plot(epoch,y,label=label,color=color,lw=2)
        ax.axvline(best,color='#ac4155',ls='--',lw=1.3,label=f'生成损失最小：第 {best} 轮')
        ax.axvspan(best,epoch[-1],color='#bb6e5c',alpha=.055)
        ax.set_title(title,fontproperties=prop,fontsize=13,pad=12)
        ax.set_xlabel('训练轮次（从 1 开始）',fontproperties=prop)
        ax.set_ylabel(ylabel,fontproperties=prop);ax.grid(alpha=.18)
        ax.spines[['right','top']].set_visible(False);ax.legend(prop=prop,frameon=False,fontsize=9)
    fig.suptitle('AtomPeak 本轮训练：验证表现较早触底，后期损失回升',fontproperties=prop,fontsize=17)
    fig.savefig(out/'validation_curves.svg');fig.savefig(out/'validation_curves.png',dpi=170);plt.close(fig)


def report(rows,config,metadata,audit,out):
    best=min(rows,key=lambda r:r['validation']['generation_loss']);last=rows[-1]
    bv,lv=best['validation'],last['validation'];manifest=audit['manifest']
    lines=['# AtomPeak 正式训练结果分析','',
        '分析日期：2026-10-26。原始配置、日志和三个 checkpoint 仅以只读方式打开；分析产物保存在独立目录。',
        '', '本轮正式训练已完成，按当前选模规则应优先评估第 12 轮的 `checkpoint_best.pt`。第 50 轮虽然 token accuracy 更高，但生成交叉熵和归属损失均明显变差，不能据此认为最后模型更好。',
        '', '## 1. 训练和断点完整性','',
        f'- 50 条验证记录覆盖 epoch=0…49，全部 `interrupted=false`；每次验证均为 2,196 个 batch。指标未发现 NaN/Inf。',
        f'- `tiny=false`，BART 768 维、6 层编码器/6 层解码器、12 个注意力头；单卡、batch=2、累积=8，有效 batch=16。',
        f'- 数据指纹与本地 `data/atompeak_annotated_v1` 对应。重新按导出顺序计算全部 LMDB 原始字节指纹，匹配结果：{audit["payload_fingerprint_matches_manifest"]}；两份词表校验：{all(audit["vocabulary_checksums"].values())}。',
        '- 最后断点 epoch_completed=true、step=107350，与总训练预算一致；最终学习率为 0，属于训练计划结束，不能沿用原计划直接追加更新。',
        '- 没有本轮 `.out/.err` 日志。无法重建训练集损失曲线、耗时、抢占/OOM 历史；文件复制时间不能当作训练耗时。保存文件存在且可读，不等于已经做过所有权重张量的有限性或 GPU 推理验证。',
        '', '| 文件 | 轮次（从 1 开始） | 更新步数 | 轮次完成 | 用途 |','|---|---:|---:|---|---|']
    for name,m in metadata['checkpoints'].items():
        lines.append(f'| {name} | {m["epoch_display"]} | {m["step"]:,} | {m["epoch_completed"]} | '+('按验证生成损失选出的最佳模型' if name=='checkpoint_best.pt' else '最近的滚动断点；不是额外的早期候选')+' |')
    lines += ['', 'epoch 在文件中从 0 开始：`epoch=11` 对应第 12 轮。三个文件均约 1.29 GB，包含优化器状态。',
        '', '## 2. 验证曲线','', '![验证曲线](validation_curves.png)',
        '', '| 指标 | 第 1 轮 | 第 12 轮 | 第 50 轮 |','|---|---:|---:|---:|']
    for key,name in [('generation_loss','生成交叉熵'),('loss','总验证损失'),('carbon_assignment','C 归属配对 BCE'),('proton_assignment','H 归属配对 BCE'),('carbon_alignment_total','C 综合归属损失'),('proton_alignment_total','H 综合归属损失')]:
        lines.append(f'| {name} | {rows[0]["validation"][key]:.6f} | {bv[key]:.6f} | {lv[key]:.6f} |')
    lines += [f'| Token accuracy | {rows[0]["validation"]["token_accuracy"]:.2%} | {bv["token_accuracy"]:.2%} | {lv["token_accuracy"]:.2%} |',
        '', f'第 50 轮生成交叉熵比第 12 轮高 {(lv["generation_loss"]/bv["generation_loss"]-1)*100:.2f}%。第 10 轮生成损失 0.504888、第 12 轮 0.504848，两者只差约 0.000041（0.008%）；更合理的观察是约第 10–14 轮已进入最优平台，而不是第 12 轮具有压倒性优势。',
        '', 'Token accuracy 是 teacher-forcing 下的局部 token 正确比例。交叉熵还受给正确答案分配的概率影响，且当前交叉熵按 batch 平均、accuracy 按 token 计数加权，因此两者可以反向变化。模型错误时更自信、不同长度/难度样本的影响都可能参与，现有汇总不能确定是哪一种。86.35% 不能当作整分子 Top-1，也不能当作原子归属正确率。',
        '', '验证损失持续回升符合泛化退化、过拟合或校准恶化的特征。但缺训练集日志，尚不能严格证实“训练集持续变好、验证集变差”的泛化差距，更不能排除标签噪声、模态分布和多任务优化的影响。',
        '', '## 3. 不同任务的最佳时期','', '| 验证指标 | 最低值 | 轮次 |','|---|---:|---:|']
    for key,name in [('loss','总验证损失'),('generation_loss','生成损失'),('carbon_assignment','C 配对损失'),('proton_assignment','H 配对损失'),('carbon_alignment_total','C 综合归属损失'),('proton_alignment_total','H 综合归属损失')]:
        r=min(rows,key=lambda r:r['validation'][key]);lines.append(f'| {name} | {r["validation"][key]:.6f} | {r["epoch"]+1} |')
    lines += ['', '当前 checkpoint_best 只由 generation_loss 决定，并不是按 C/H 归属表现选择。H 归属分支较早恶化；第 4、5、8 轮已经不是目前三个文件中可恢复的独立模型，不能靠 checkpoint_prev 恢复。后续若以原子归属为主要成果，需要单独保存归属指标最佳模型。',
        '', '## 4. 这一轮实际使用的数据','',
        f'共 {manifest["records"]:,} 条分子记录，连接层分组数 {audit["connectivity_groups"]:,}，跨 split 的相同连接层组数 {audit["cross_split_connectivity_groups"]}。这项检查说明结构连接层未跨划分，不保证文献、作者或近邻骨架也全部隔离。',
        '', '| 划分 | 分子数 | 双谱 | 仅 C | 仅 H |','|---|---:|---:|---:|---:|']
    for split,t in audit['splits'].items():
        lines.append(f'| {split} | {t["records"]:,} | {t["modes"]["dual"]:,} | {t["modes"]["carbon_only"]:,} | {t["modes"]["proton_only"]:,} |')
    lines += ['', '全库双谱 7,127（16.57%）、仅 C 25,308（58.85%）、仅 H 10,571（24.58%）。这一轮并非全部样本都同时输入 C/H，也没有把数据库几十万条无归属分子加入训练。所有导出记录的 has_supervision=true；`labelled_fraction=0.5` 在没有无监督样本时不会凭空制造一半无监督数据。',
        '', '训练采样器对有监督样本按双谱/C单谱/H单谱 = 50%/25%/25% 的期望比例有放回抽样，再对双谱施加 20% 随机隐藏一支。该比例是采样策略，不是原库真实比例。验证保持原始模态分布，报告应按双谱、C单谱、H单谱进一步拆分。',
        '', '按分子 origin 汇总，43,006 条中 42,854 条来自 nmrshiftdb2、152 条来自 NMRexp；该字段不是逐峰标签来源，不能据此得到“原库标注 vs 手工标注”的逐峰成绩。当前训练导出缺逐峰 origin，后续分来源评估应使用冻结版本的来源映射或另行升级导出字段，不能拿今天的主库标签无条件替换当时快照。',
        '', '### 标签覆盖和辅助损失','', '| 划分 | C 确定标签峰 | H 确定标签峰 | C 弱标签峰 | H 弱标签峰 | C 无序组 | H 无序组 |','|---|---:|---:|---:|---:|---:|---:|']
    for split,t in audit['splits'].items():
        c,h=t['branches']['carbon'],t['branches']['proton'];lines.append(f'| {split} | {c["hard_peaks"]:,} | {h["hard_peaks"]:,} | {c["weak_peaks"]:,} | {h["weak_peaks"]:,} | {c["sets"]} | {h["sets"]} |')
    lines += ['', '验证集中有 26,523 个 C 确定标签峰、10,673 个 H 确定标签峰，与日志的平均每 batch 已知峰数乘以 2,196 一致。每一核种只有 1 个无序组；测试集两种核种均无无序组。因此极小的 sets 损失不能支撑“组归属能力很好”的结论，这个测试划分也无法给出原生无序组泛化成绩。',
        '', '全部三个划分的 mass_supervision_valid 都为 0，spectrum_complete 也都为 0。因此 carbon_mass/proton_mass 始终为 0 的原因是质量监督没有启用，不是模型实现了完美峰质量/积分匹配。程序要求谱完整性明确且标签满足条件才启用该项；不要为追求非零损失把未知完整性改成 complete。',
        '', 'carbon_shift/proton_shift 是带预测不确定性的归一化高斯 NLL，公式约为 `0.5*((预测-真实)/scale)^2*exp(-2*log_sigma)+log_sigma`，不是 ppm 单位 MAE。NLL 可以为负，现有日志不能直接还原为 C/H 的 ppm 误差；应另算 MAE/RMSE 和区间覆盖率。',
        '', '### 溶剂条件','', '本次 use_solvent=true，但 temperature/field_strength=false；分子式、H 峰型、H 积分和 J token 开启。开启某字段不意味着其原始数值普遍齐全。',
        '', '| 训练分支 | 当前存在的谱 | missing token | 非 missing token |','|---|---:|---:|---:|']
    for b,title in [('carbon','C'),('proton','H')]:
        t=audit['splits']['train'];n=t['branches'][b]['records'];missing=sum(v for k,v in t['solvent_tokens'][b].items() if k.endswith('_missing'))
        lines.append(f'| {title} | {n:,} | {missing:,}（{missing/n:.2%}） | {n-missing:,}（{(n-missing)/n:.2%}） |')
    lines += ['', '非 missing 中包括 other 及混合/未充分规范化的溶剂字符串，不能全当作准确的具体溶剂类别。当前没有无溶剂模型的固定 split 对照，不能用这轮训练证明加入溶剂带来提升。',
        '', '## 5. 建议下一步','',
        '1. 先保留这次结果。不要因为第 50 轮 token accuracy 更高，就覆盖第 12 轮 best，也不建议直接续跑到更多轮次。',
        '2. 使用同一个验证子集、相同生成参数比较 best 与 last，先运行生成基线 `--uot-weight 0`。只有实际生成结果才能比较整分子 Top-1/3/10；先不使用 test 调参。',
        '3. 在验证集确定是否需要 UOT 重排及其权重。生成损失最佳不保证重排最佳；C/H 归属分支本身也有后期退化迹象。',
        '4. 固定选模和重排规则后，在完整 test 上评估一次，按双谱/C单谱/H单谱、溶剂已知/缺失拆分。记录输入包含分子式，避免与未提供分子式的模型作不公平比较。',
        '5. 为归属头增加单独评估：C/H 分开，确定逐峰标签、弱标签和无序组分别统计；给出集合精确匹配、Top-1/3、可评估峰数及 ppm 误差。给定真实结构做归属诊断，与未知结构的谱到分子生成是两种实验，应分开报告。',
        '6. 下一轮再考虑提前停止、分别保存生成/归属最佳模型，以及正则化、学习率和损失权重的受控实验。目前没有实施这些训练代码变更。',
        '', '当前没有本轮的自由生成 JSONL 或 test 指标文件，所以本报告不提供未经计算的整分子准确率、归属准确率或模型提升结论。',
        '', '## 6. 可直接在训练机运行的初步对照命令','',
        '以下复用原训练数据目录，输出选新的文件名；200 条 valid 仅用于流程与候选错误检查，最终选择应扩大到完整 valid。', '', '```bash',
        'python infer_atompeak_partial.py \\',
        '  --data /data/dataset/atompeak_partial_v1 \\',
        '  --checkpoint results/atompeak_partial_scratch_v1/checkpoint_best.pt \\',
        '  --output results/atompeak_partial_scratch_v1/valid_best200_generation.jsonl \\',
        '  --split valid --max-samples 200 --beam-size 10 --uot-weight 0', '',
        'python infer_atompeak_partial.py \\',
        '  --data /data/dataset/atompeak_partial_v1 \\',
        '  --checkpoint results/atompeak_partial_scratch_v1/checkpoint_last.pt \\',
        '  --output results/atompeak_partial_scratch_v1/valid_last200_generation.jsonl \\',
        '  --split valid --max-samples 200 --beam-size 10 --uot-weight 0',
        '```', '',
        '每条命令生成一个逐样本 JSONL 和同名 `.metrics.json`。Top-K 比较带立体信息的规范 SMILES；当前候选先将有效结构放在无效结构前，因此 valid_top1_rate 是该排序后第一候选的合法性，不是未过滤的原始 beam 首选合法率。候选去重/过滤和分母口径应在正式评估中注明。',
        '', '## 7. 分析范围与产物','',
        '- 读取 run_config、50 行验证记录、三个 checkpoint 元数据，以及匹配数据导出的全部 LMDB 和索引。',
        '- `analysis.json`：数值摘要和最优轮次；`checkpoint_metadata.json`：三个断点的核查信息；`dataset_audit.json`：数据组成、监督与溶剂 token 覆盖。',
        '- `validation_curves.svg`：矢量图；`validation_curves.png`：预览图。',
        '- 本机只用 CPU 读取元数据，不重新训练，不运行收费接口，不修改原始权重或数据。',
        '- 元数据读取环境与训练环境不同：本地 PyTorch '+metadata['reader_torch_version']+'，原训练 '+config.get('torch_version','未知')+'。本报告没有在不同运行时重放模型数值结果。',
        '- 没有读取主标注库来替换本次训练的快照标签。']
    (out/'training_analysis.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--data',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise ValueError('Choose a new analysis folder')
    args.output.mkdir(parents=True)
    runtime=ROOT/'tmp/atompeak_partial_runtime'
    if runtime.exists():sys.path.insert(0,str(runtime))
    rows=[json.loads(s) for s in (args.run/'metrics.jsonl').read_text(encoding='utf-8').splitlines() if s.strip()]
    config=json.loads((args.run/'run_config.json').read_text(encoding='utf-8'))
    if [r['epoch'] for r in rows]!=list(range(config['args']['epochs'])):
        raise ValueError('Expected complete ordered epoch records for this report')
    if any(not math.isfinite(v) for r in rows for v in r['validation'].values() if isinstance(v,float)):
        raise ValueError('Found nonfinite validation values')
    metadata=checkpoint_metadata(args.run);write_json(args.output/'checkpoint_metadata.json',metadata)
    print('checkpoint metadata read',flush=True)
    audit=audit_dataset(args.data);write_json(args.output/'dataset_audit.json',audit)
    if not audit['payload_fingerprint_matches_manifest']:raise ValueError('LMDB payload fingerprint differs')
    if any(m['dataset_fingerprint']!=audit['manifest']['source_snapshot_fingerprint'] for m in metadata['checkpoints'].values()):
        raise ValueError('Checkpoint belongs to a different dataset')
    print('dataset and fingerprints verified',flush=True)
    extrema={}
    for key in rows[0]['validation']:
        if key in ('batches','interrupted'):continue
        row=(max if key=='token_accuracy' else min)(rows,key=lambda r:r['validation'][key])
        extrema[key]={'epoch_display':row['epoch']+1,'step':row['step'],'value':row['validation'][key]}
    summary={'run':str(args.run.resolve()),'data':str(args.data.resolve()),'epochs':len(rows),
             'last_step':rows[-1]['step'],'extrema':extrema,'rows':rows,
             'input_checksums':{name:digest(args.run/name) for name in ('run_config.json','metrics.jsonl')},
             'config':config,'model_predictions_evaluated':False}
    write_json(args.output/'analysis.json',summary)
    plots(rows,args.output);report(rows,config,metadata,audit,args.output)
    print(json.dumps({'report':str((args.output/'training_analysis.md').resolve()),'best_generation':extrema['generation_loss'],
                      'last_generation_loss':rows[-1]['validation']['generation_loss'],
                      'dataset_fingerprint_verified':True},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
