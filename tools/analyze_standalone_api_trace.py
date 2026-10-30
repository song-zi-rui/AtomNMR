"""Offline analysis of the requested one-molecule API trace. No network/DB calls."""
import csv
import json
from pathlib import Path
import sys
from html import escape

folder=Path(sys.argv[1]).resolve()
read=lambda name:json.loads((folder/name).read_text('utf-8-sig'))
summary=read('summary.json');selection=read('selection.json');progress=read('group/progress.json')
events=[json.loads(x) for x in (folder/'group/events.jsonl').read_text('utf-8').splitlines() if x]
text_events=[x for x in events if x.get('reasoning_added',0) or x.get('answer_added',0)]
duration=progress['elapsed_s'];last=text_events[-1]['elapsed_s'] if text_events else None
gaps=[];previous=0
for event in text_events:
    t=event['elapsed_s']
    if t-previous>=10:gaps.append({'start_s':previous,'end_s':t,'duration_s':t-previous})
    previous=t
if duration-previous>=10:gaps.append({'start_s':previous,'end_s':duration,'duration_s':duration-previous})
with (folder/'output_gaps.csv').open('w',encoding='utf-8-sig',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=['start_s','end_s','duration_s']);writer.writeheader();writer.writerows(gaps)
bins=[{'minute':i+1,'reasoning_characters':0,'answer_characters':0} for i in range(int(duration//60)+1)]
for event in text_events:
    i=min(int(event['elapsed_s']//60),len(bins)-1)
    bins[i]['reasoning_characters']+=event.get('reasoning_added',0);bins[i]['answer_characters']+=event.get('answer_added',0)
with (folder/'output_by_minute.csv').open('w',encoding='utf-8-sig',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(bins[0]));writer.writeheader();writer.writerows(bins)
maximum=max((x['reasoning_characters']+x['answer_characters'] for x in bins),default=1) or 1
svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="330" viewBox="0 0 1200 330">',
     '<rect width="1200" height="330" fill="#ffffff"/>',
     '<g font-family="Microsoft YaHei, sans-serif" fill="#193931"><text x="52" y="35" font-size="22">组内 API：每分钟收到的文字字符数</text>',
     '<text x="52" y="61" font-size="14">蓝色为思考长度统计，橙色为归属答案；不含思考正文，字符数不是 token 数。</text>']
for y in [0,.5,1]:
    yp=255-160*y;svg.append(f'<path d="M55 {yp} H1160" stroke="#e0e8e3"/><text x="8" y="{yp+5}" font-size="12">{round(maximum*y)}</text>')
width=1100/len(bins)
for i,row in enumerate(bins):
    x=55+i*width;r=160*row['reasoning_characters']/maximum;a=160*row['answer_characters']/maximum
    svg.append(f'<rect x="{x+2}" y="{255-r}" width="{max(1,width-4)}" height="{r}" fill="#2375a6"/>')
    svg.append(f'<rect x="{x+2}" y="{255-r-a}" width="{max(1,width-4)}" height="{a}" fill="#d88a32"/>')
    if i%5==0:svg.append(f'<text x="{x}" y="280" font-size="12">{i}分</text>')
svg.append('</g></svg>');(folder/'timeline.svg').write_text(''.join(svg),encoding='utf-8')
probe=read('network_probe.json');controls=read('comparison_controls.json')
analysis={'inference_calls':1,'diagnostic_get_requests':1,'database_written':False,
          'total_s':duration,'last_text_s':last,'no_text_tail_s':duration-last if last is not None else duration,
          'largest_no_text_gap_s':max((g['duration_s'] for g in gaps),default=0),
          'local_deadline_reached':duration>=selection['max_seconds'],
          'response_complete':bool(progress['finish_reason']=='stop' and progress['done_marker']),
          'answer_characters':progress['answer_characters'],'reasoning_characters':progress['reasoning_characters'],
          'usage_reported':bool(progress['usage']),'accuracy':None,'accuracy_reason':'No completed group answer or independent reference labels',
          'termination_note':'The local 2400-second watchdog shut down the socket. The resulting BrokenPipeError is a shutdown artifact, not evidence of the original stall cause.',
          'network_probe':probe,'controls':controls}
(folder/'analysis.json').write_text(json.dumps(analysis,ensure_ascii=False,indent=2),encoding='utf-8')
summary.update(inference_calls=1,diagnostic_get_requests=1,api_calls_note='api_calls counts annotation inference only; one separate read-only GET /models probe was also sent')
(folder/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
minutes=lambda seconds:f'{seconds/60:.2f} 分钟'
report=f'''# 单分子 API 对照：结果与排障结论

分子：**{selection['formula']}**，`{selection['molecule_id']}`。选择范围是正式自动标注和早期标注的成功官方记录；本次没有把新结果或测试任务写入系统。

本次组内请求未取得完整归属。它很快开始返回思考文字，之后长时间没有新的文字，直到独立测试的40分钟上限。不能把它概括为“分子复杂，所以一直在正常思考”。

| 项目 | 官方历史记录 | 组内本次请求 |
| --- | --- | --- |
| 总耗时 | {selection['official_elapsed_s']:.2f} 秒 | {duration:.2f} 秒 |
| 首个文字 | 旧非流式记录未保存 | {progress['ttft_s']:.2f} 秒 |
| 完整归属 | 23 个峰，正常结束 | 未取得 |
| 思考信息 | 55,100 token | {progress['reasoning_characters']:,} 字符，用量未回报 |
| 归属答案通道 | 完整答案已保存 | {progress['answer_characters']} 字符 |
| 结束原因 / DONE | stop / 非流式 | {progress['finish_reason']} / {progress['done_marker']} |

最后一次收到文字约在第 **{minutes(last or 0)}**，之后到本地结束约 **{minutes(analysis['no_text_tail_s'])}** 没有新文字。精确时间见 `output_gaps.csv`，逐事件计数见 `group/events.jsonl`。

![文字接收时间分布]({(folder/'timeline.svg').as_posix()})

## 这次能确认什么

1. 不是初始连接失败：HTTP 200，约1.31秒开始收到文字；没有收到401/429或模型名错误。
2. 不是在本地解析归属或写数据库时卡住：答案尚未开始，测试不执行数据库写入。本地进程在长时间等待中累计CPU约4.5秒，内存约57MB；未见本机算力瓶颈。
3. 停顿期间只读查询 `/models` 用时 **{probe['elapsed_s']:.3f} 秒**，模型名存在。说明当时网关可达；不能据此证明GPU推理正常，也不能完全排除单条流连接的问题。
4. 不能归因于V3严格schema：这次与官方一样使用V2和json_object，输入及提示词校验一致。
5. 本次没有任何归属答案和结束标记，因此不是已有完整答案却只漏了DONE的证据。也没有收到length，无法把它诊断为已确认的token上限错误。

## 还不能确认什么

客户端无法区分：服务端推理挂起、调度暂停、模型进程异常，或该条流式转发停止。TCP保持ESTABLISHED不代表模型还在计算。没有最终用量，不能把用量/费用填成零；本地关闭接收也不能保证远端停止。

原始记录中的BrokenPipeError/WinError 10058发生在本地40分钟定时器主动关闭连接后，是结束测试时的伴随错误，不能倒推它导致了前面的28分钟停滞。实际生成请求没有收到服务端提供的明确错误原因。

这是一例诊断，不代表所有组内调用都失败，也不证明模型本身的化学能力差。该分子没有独立原库/人工标准标签；官方答案是对照结果，不是标准答案。组内未完成，所以不存在可计算的本次化学正确率。

## 对照的边界

两次分子、完整C/H资料、提示词、V2格式及65,536生成上限相同。差别包括接口、模型部署、原生思考参数、流式/非流式及等待机制。官方旧记录的timeout字段为200秒、实际耗时315秒，不能把那个字段当作本轮同样的严格总时限。

本次复现现有组内适配的 `thinking=true, reasoning_effort=50`。数值50是思考强度，不能当作50个token；官方与开源端同名high也不保证相同数值预算。参见 [vLLM DeepSeek-V4.1说明](https://recipes.vllm.ai/deepseek-ai/DeepSeek-V4.1-Flash)。没有确认组内权重版本、硬件和采样默认值，不把本例写成严格模型能力排名。

## 接下来怎么用这个结果

不要把等待上限继续无条件加大，作为正式批量标注的唯一解决方案。先对“开始响应后长时间无新数据”建立单独状态，记录最后文字时间、答案开始时间与结束标记。组内开启思考仍不适合直接作为已经验证稳定的24小时批量默认项。

后续如需比较化学效果，应另做关闭思考或较低强度的对照并明确标为不同配置；本次只发送了一次分子生成请求，没有自动重试。另有一次不生成答案的GET /models检查。

## 看哪些文件

- `comparison.html`：结构、原子编号、全部实验条件与官方逐峰归属，组内栏明确显示未取得结果。
- `input.json`、`group/request.json`、`comparison_controls.json`：冻结输入、实际请求和一致性验证，不含API Key。
- `official/response.json`：官方完整最终答案与回报用量。
- `group/answer.txt`：本次答案通道；空文件表示没有收到答案，不是程序遗漏导出。
- `group/progress.json`、`group/events.jsonl`、`network_probe.json`、`analysis.json`：详细诊断。
- `peak_comparison.csv`：全部23个峰的对照表；组内空白不计为错误归属。

思考只保留长度与时间统计，没有保存正文。
'''
(folder/'analysis.md').write_text(report,encoding='utf-8')
page=folder/'comparison.html'
if page.exists():
    value=page.read_text('utf-8')
    if 'id="timelineSection"' not in value:
        value=value.replace('<section><h2>结构与实验资料</h2>',
            '<section id="timelineSection"><h2>接收时间分布</h2><img src="timeline.svg" style="width:100%" alt="每分钟收到的文字"/><p>完整结论见同目录 analysis.md；长时间无新文字不能等同持续正常思考。</p></section><section><h2>结构与实验资料</h2>')
    page.write_text(value,encoding='utf-8')
print(json.dumps(analysis,ensure_ascii=True))
