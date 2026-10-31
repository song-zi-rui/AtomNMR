"""Render the standalone API diagnostic files into a readable local report."""
import argparse
from collections import Counter
import csv
from datetime import datetime,timezone
import html
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def main():
    p=argparse.ArgumentParser();p.add_argument('folder',type=Path);a=p.parse_args();folder=a.folder
    read=lambda name:json.loads((folder/name).read_text('utf-8'))
    selection=read('selection.json');ctx=read('input.json');official=read('official/assignments.json');diag=read('official/diagnostics.json')
    progress=read('group/progress.json');group=read('group/assignments.json') if (folder/'group/assignments.json').exists() else {}
    summary=read('summary.json') if (folder/'summary.json').exists() else {}
    group_cfg=read('group/config.json')
    settings_note=('本次组内开启思考' if group_cfg.get('thinking') else '本次组内关闭思考')+'，生成上限 '+str(group_cfg['max_tokens'])+' token，独立总时限 '+str(group_cfg['timeout'])+' 秒。'
    esc=lambda v:html.escape(str(v if v is not None else '—'))
    sec=lambda v:'未记录' if v is None else f'{v:.2f} 秒'
    number=lambda v:'未回报' if v is None else f'{v:,}'
    atoms=lambda v:'、'.join(ctx['atoms'][i]['element']+str(i) for i in v) if v else '—'
    statuses={'assigned':'单原子归属','equivalent_group':'等价组','overlap':'重叠归属','ambiguous':'歧义','unassigned':'未归属'}
    count=lambda output:'；'.join(statuses.get(k,k)+' '+str(v) for k,v in Counter(x['status'] for x in output.get('assignments',[])).items()) or '尚无可校验归属'
    usage=progress.get('usage',{});reason=(usage.get('completion_tokens_details') or {}).get('reasoning_tokens')
    running=progress['state'] in ('connecting','waiting','reasoning','answering')
    observed=(datetime.now(timezone.utc)-datetime.fromisoformat(selection['started_at'])).total_seconds() if running else progress['elapsed_s']
    metrics=[('总耗时',sec(selection['official_elapsed_s']),sec(observed)),
        ('首个文字输出','历史非流式，未记录',sec(progress['ttft_s'])),
        ('开始输出归属答案','历史未记录',sec(progress['first_answer_s'])),
        ('收到 finish_reason','stop',str(progress['finish_reason'])),
        ('收到 [DONE]','历史非流式',str(progress['done_marker'])),
        ('生成 token 总数',number(diag.get('completion_tokens')),number(usage.get('completion_tokens'))),
        ('其中思考 token',number(diag.get('reasoning_tokens')),number(reason)),
        ('归属状态分布',count(official),count(group))]
    rows=[];om={(x['spectrum_id'],x['peak_id']):x for x in official.get('assignments',[])};gm={(x['spectrum_id'],x['peak_id']):x for x in group.get('assignments',[])}
    for s in ctx['spectra']:
        for peak in s['peaks']:
            key=(s['spectrum_id'],peak['peak_id']);old=om.get(key,{});new=gm.get(key,{})
            def cell(x):
                candidate=' / '.join(atoms(v) for v in x.get('candidates',[]))
                return '<strong>'+esc(atoms(x.get('atom_indices',[])))+'</strong><div>'+esc(statuses.get(x.get('status'),x.get('status','尚未取得结果')))+'</div>'+('<div>候选：'+esc(candidate)+'</div>' if candidate else '')+'<p>'+esc(x.get('brief_reason',''))+'</p>'
            rows.append('<tr><td>'+esc(s['nucleus'])+'<br>'+esc(peak['peak_id'])+'</td><td>'+esc(peak.get('shift_ppm'))+'</td><td>'+cell(old)+'</td><td>'+cell(new)+'</td></tr>')
    structure=''
    mol=Chem.MolFromSmiles(ctx['smiles'])
    if mol and len(ctx['atoms'])==mol.GetNumAtoms() and all(mol.GetAtomWithIdx(x['index']).GetSymbol()==x['element'] and sorted(n.GetIdx() for n in mol.GetAtomWithIdx(x['index']).GetNeighbors())==sorted(x['neighbors']) for x in ctx['atoms']):
        draw=rdMolDraw2D.MolDraw2DSVG(880,440);draw.drawOptions().addAtomIndices=True
        draw.DrawMolecule(mol);draw.FinishDrawing();structure=draw.GetDrawingText()
        (folder/'structure.svg').write_text(structure,encoding='utf-8')
    atom_rows=''.join('<tr>'+''.join('<td>'+esc(x.get(k))+'</td>' for k in ['index','element','hydrogens','aromatic','hybridization','neighbors'])+'</tr>' for x in ctx['atoms'])
    conditions=''.join('<p><b>'+esc(s['nucleus'])+'</b>：'+esc(json.dumps(s['conditions'],ensure_ascii=False))+'</p>' for s in ctx['spectra'])
    state=progress['state'];notes='本次组内结果未写入标注系统。两模型的一致率不是化学正确率；单个分子也不能代表整体模型性能。'
    extra=''
    probe=read('network_probe.json') if (folder/'network_probe.json').exists() else {}
    if probe:
        extra+='<p>只读 /models 检查：'+esc('成功' if probe['ok'] else '失败')+'，耗时 '+sec(probe['elapsed_s'])+'。这只证明网关可达，不代表生成请求正常。</p>'
    if (folder/'local_runtime.json').exists():
        local=json.loads((folder/'local_runtime.json').read_text('utf-8-sig'))
        extra+='<p>本地测试进程累计 CPU 时间 '+sec(local['cpu_seconds'])+'；内存约 '+str(round(local['working_set_bytes']/1024**2))+' MB。</p>'
    if running:extra+='<p>此报告为运行中快照；最近进展记录在第 '+sec(progress['elapsed_s'])+'。正在等待的空白时段不会增加生成字符数。</p>'
    agreement=''
    if summary:
        agreement=f"双方均给确定归属的峰 {summary['both_hard_peaks']} 个，其中原子集合相同 {summary['same_hard_atom_set']} 个。图等价未作额外放宽；请结合下表人工检查。"
    body=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><title>官方与组内 API 单分子对照</title><style>
body{{font:16px/1.7 "Microsoft YaHei",sans-serif;color:#193931;background:#f4f7f5;margin:0}}main{{max-width:1300px;margin:auto;padding:32px}}section{{background:white;border:1px solid #d9e4df;border-radius:12px;padding:24px;margin:20px 0}}h1{{font-size:28px}}h2{{font-size:20px}}table{{border-collapse:collapse;width:100%;font-size:16px}}td,th{{text-align:left;vertical-align:top;padding:12px;border-bottom:1px solid #e4eae6}}th{{background:#edf4ef}}td p{{margin:6px 0}}code{{overflow-wrap:anywhere}}pre{{white-space:pre-wrap;word-break:break-all;font-size:14px}}svg{{width:100%;height:auto;max-height:440px}}.note{{background:#fff4da;padding:14px;border-radius:8px}}details{{margin:16px 0}}.table{{overflow:auto}}</style></head><body><main>
<h1>{esc(selection['formula'])} · 官方与组内 API 对照</h1><p><code>{esc(selection['molecule_id'])}</code> · 运行状态：<b>{esc(state)}</b></p><p class="note">{esc(notes)}</p>
<section><h2>运行过程</h2><p>沿用官方冻结的分子资料、提示词和V2格式，组内使用流式接收。{esc(settings_note)}没有自动重试。</p><table><tr><th>指标</th><th>官方历史结果</th><th>组内本次结果</th></tr>{''.join('<tr><td>'+esc(k)+'</td><td>'+esc(v)+'</td><td>'+esc(w)+'</td></tr>' for k,v,w in metrics)}</table><p>组内已接收思考 {progress['reasoning_characters']:,} 字符、归属答案 {progress['answer_characters']:,} 字符；字符不等于 token。思考正文没有保存。</p>{extra}<p>{esc(agreement)}</p></section>
<section><h2>结构与实验资料</h2>{structure}<code>{esc(ctx['smiles'])}</code>{conditions}<details><summary>原子编号与环境</summary><div class="table"><table><tr><th>编号</th><th>元素</th><th>连接 H 数</th><th>芳香</th><th>杂化</th><th>相邻原子</th></tr>{atom_rows}</table></div></details></section>
<section><h2>逐峰归属对照</h2><p>下面的“依据”是模型输出的简要化学解释；官方结果是对照参考，不视为标准答案。</p><div class="table"><table><tr><th>核种 / 峰</th><th>δ / ppm</th><th>官方：原子与依据</th><th>组内：原子与依据</th></tr>{''.join(rows)}</table></div></section>
<section><h2>完整归属输出</h2><details><summary>官方的原始最终答案</summary><pre>{esc(read('official/response.json').get('content','') if (folder/'official/response.json').exists() else '')}</pre></details><details><summary>组内的最终答案通道（可能不完整）</summary><pre>{esc((folder/'group/answer.txt').read_text('utf-8') if (folder/'group/answer.txt').exists() else '')}</pre></details><details><summary>组内过程统计</summary><pre>{esc(json.dumps(progress,ensure_ascii=False,indent=2))}</pre></details></section>
</main></body></html>'''
    (folder/'comparison.html').write_text(body,encoding='utf-8')
    lines=['# 官方与组内 API 单分子对照',f"\n分子：{selection['formula']}，`{selection['molecule_id']}`。状态：{state}。",'\n'+notes,'\n| 指标 | 官方 | 组内 |','| --- | --- | --- |']
    lines += [f'| {k} | {v} | {w} |' for k,v,w in metrics]
    lines += ['\n'+agreement,'\n文件说明：','- comparison.html：结构、原子环境、逐峰归属及简要依据。','- peak_comparison.csv：逐峰表格。','- input.json、group/request.json：原始资料与本次请求参数（不含密钥）。','- official/：官方历史答案与统计。','- group/answer.txt：组内归属答案通道。','- group/progress.json：耗时、思考长度、结束状态与用量。','- group/events.jsonl：每个事件的时间与字符增量；不含思考正文。','- summary.json：最终汇总。','\n原始标注数据库只读；本次结果没有入库，也没有更改当前模型配置。']
    (folder/'report.md').write_text('\n'.join(lines),encoding='utf-8')
    print(folder/'comparison.html')


if __name__=='__main__':main()
