"""Read diagnostic artifacts and write a follow-up report; no network or DB."""
from pathlib import Path
import csv
import hashlib
import json
import zipfile

ROOT=Path(__file__).resolve().parents[1]
folder=ROOT/'outputs/transport_root_cause_20261031'
old=ROOT/'outputs/official_slowest_group_comparison_20261030'
read=lambda path:json.loads(path.read_text('utf-8-sig'))
original=read(old/'summary.json');controls=[]
for name in ['group_thinking_off_control','group_off_65536_control']:
    path=folder/name;value=read(path/'summary.json');cfg=read(path/'group/config.json');progress=read(path/'group/progress.json')
    controls.append({'name':name,'thinking':cfg['thinking'],'max_tokens':cfg['max_tokens'],
        'elapsed_s':progress['elapsed_s'],'state':progress['state'],'finish_reason':progress['finish_reason'],
        'done_marker':progress['done_marker'],'prompt_tokens':progress['usage'].get('prompt_tokens'),
        'completion_tokens':progress['usage'].get('completion_tokens'),'reasoning_tokens':(progress['usage'].get('completion_tokens_details') or {}).get('reasoning_tokens'),
        **{k:value[k] for k in ['group_assignment_count','group_hard_peaks','both_hard_peaks','same_hard_atom_set','group_decode_error']}})
    # Verify this is the same frozen scientific input and original prompt.
    current=read(path/'group/request.json');prior=read(old/'group/request.json')
    assert current['messages']==prior['messages']
    assert value['database_written'] is False
with (folder/'comparison.csv').open('w',encoding='utf-8-sig',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(controls[0]));writer.writeheader();writer.writerows(controls)
result={'date':'2026-10-31','controls':controls,'local_10058':read(folder/'exact_10058_reproduction.json'),
    'current_gateway_checks':read(folder/'gateway_checks.json'),'new_inference_requests':2,'database_written':False,
    'confirmed':'Original WinError 10058 follows the local deadline shutdown; same-input thinking-off controls complete in ~26–29 seconds.',
    'unresolved':'Cannot distinguish long-reasoning backend failure, scheduling, per-connection network loss or streaming gateway stalls from client-only evidence.',
    'causal_limit':'The stalled thinking-on run was on Oct 30; controls were on Oct 31. Different server load/weights/time remain confounders. No independent gold labels.'}
(folder/'findings.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
report='''# 错误继续排查：10058含义与同分子对照

日期：2026-10-31。分子 C18H16INO3，M004facedebcf2db5c5244284。所有测试均不写入标注数据库、任务列表或当前API设置。

## 先回答10058是什么意思

这次的 BrokenPipeError / WinError 10058 表示：程序已经关闭了这条连接的收发方向，另一个仍在读取的操作又尝试使用它。它是本地40分钟定时器停止等待时出现的伴随异常，不能解释前面为什么没有新数据。

顺序是：1.31秒开始收到文字 → 716.08秒最后一次收到文字 → 后续1684.02秒无新文字 → 2400秒定时器主动调用shutdown → 读取操作抛出10058。

本机已用两个本地连接复现完全相同的BrokenPipeError/10058，未调用组内API。TLS控制还显示：不同读取状态下，本地主动终止也可能表现为TimeoutError。不能只凭异常名称判断远端故障。

参考：[Microsoft Winsock错误码](https://learn.microsoft.com/en-us/windows/win32/winsock/windows-sockets-error-codes-2)、[shutdown说明](https://learn.microsoft.com/en-us/windows/win32/api/winsock2/nf-winsock2-shutdown)。

## 新增两次独立标注对照

保留同一个冻结输入、原提示词、完整C/H谱、V2、json_object、相同模型标识与TLS配置。数据与提示词已逐项校验一致。今天只新增这两次生成请求，均未自动重试。

| 组内配置 | 输出上限 | 耗时 | 输出结果 | 结束原因 |
| --- | ---: | ---: | --- | --- |
| 昨日开启思考，强度50 | 65536 | 2400.10秒 | 无归属答案 | 本地时限，未收到结束标记 |
| 今日关闭思考 | 8192 | 28.94秒 | 23峰条目完整返回 | stop + DONE |
| 今日关闭思考 | 65536 | 26.03秒 | 23峰条目完整返回 | stop + DONE |

第一条关闭思考对照同时降低了输出上限，存在混杂因素，所以又保留原65536上限做了第二条控制。第二条说明65536并非一个必然让请求失败的值。

第二条关闭思考对照回报输入2625、生成673、思考0 token。完整返回的含义是格式和目标条目齐全，不等于每峰归属正确。

## 当前能把原因缩小到哪里

- 本地40分钟关闭操作导致10058：已确认，并已本地复现。
- 模型名、基本密钥、连接证书、客户端完全无法接收V2：不符合这两次成功返回的观察。
- 65536这个上限本身必然导致失败：被同预算关闭思考的成功对照排除。
- 问题更集中于组内长思考或长时间生成/流式连接的稳定性：与现有证据一致，优先查这一层。
- 具体是GPU推理、调度/KV资源、网关转发，还是单条网络连接静默丢失：仍不能从客户端日志确定。

原思考请求发生在10月30日，今天对照在10月31日。服务器负载、路由、权重更新和采样随机性未完全控制，不能断言思考开关就是唯一原因，也不能据此断言组内硬件损坏。

另外一个日志盲区：原测试主要记录已解析SSE事件和文字增量，没有分别记录原始字节、完整行、心跳时间。如果网关只送出半行或缓存未刷出，客户端也可能表现为“没有新的可解析文字”。下次长请求应增加这几个时间点，不应把最后一次频道为reasoning一直显示成“模型还在思考”。

## 网关的只读检查

今天模型列表0.13秒返回HTTP 200，目标模型存在。公开的/health/liveliness和/health/readiness返回404，说明当前入口未暴露这两条路由，不能把404解释为模型坏了。未尝试管理接口或绕过权限。

LiteLLM将网关探针和真实模型生成健康检查分开；网关活着不能证明GPU模型可以完成长请求。参考：[LiteLLM健康检查说明](https://docs.litellm.ai/docs/proxy/health)。本次只调用不生成答案的探针，未调用会对各模型发起推理的/health。

## 返回速度与归属质量分开看

| 关闭思考配置 | 返回峰条目 | 确定归属峰 | 与官方原子集合相同 | 未归属峰 |
| --- | ---: | ---: | ---: | ---: |
| 8192 | 23 | 22 | 12/22 | 1 |
| 65536 | 23 | 20 | 10/20 | 3 |

本分子没有独立人工/原库标准归属，这里的相同只表示与官方输出一致，不是正确率。65536对照将C2同时分到C001和C014，触发同一碳对应多个C峰的复核提示；8192对照在C7上也有类似提示。存在慢交换/多构象例外，因此保留为检查提示，不直接修改标签。

不能仅因为关闭思考能快速返回，就直接把组内结果无筛选地用于训练。还需在固定、有标准归属的样本上评估准确率、歧义处理及重复原子问题。

## 对本地程序的处理

后续独立测试已把本地总时限作为主错误，将关闭连接时的原始异常放进secondary_shutdown_exception，并记录定时器触发时间。旧40分钟测试的原始文件未改写。

这只修正诊断表达，不会让已经停住的远端模型继续生成。正式系统后续应分别显示：等待首字、接收思考、接收答案、长时间无新数据；并保留首字时间、最后字节/心跳时间、最后文字时间、finish_reason和DONE。静默检测可以避免长时间空等，不能当作服务端根因修复；结果未知的请求仍不应自动重发。

## 材料

- group_off_65536_control/comparison.html：保留原预算的同分子逐峰对照，建议先看这份。
- group_thinking_off_control/comparison.html：8192对照。
- 两个子目录的group/answer.txt、group/assignments.json、group/request.json、group/events.jsonl：完整最终答案、实际参数、时间事件；不含API Key，不保存内部思考正文。
- exact_10058_reproduction.json：完全相同错误的本地复现。
- local_tls_controls.json：正常结束、主动停止和远端提前结束的控制记录。
- gateway_checks.json：只读连通与健康路由检查。
- comparison.csv、findings.json：汇总。
'''
(folder/'analysis.md').write_text(report,encoding='utf-8')
manifest={x.relative_to(folder).as_posix():hashlib.sha256(x.read_bytes()).hexdigest() for x in sorted(folder.rglob('*')) if x.is_file() and x.name!='manifest.json'}
(folder/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
with zipfile.ZipFile(folder.with_suffix('.zip'),'w',compression=zipfile.ZIP_DEFLATED) as z:
    for x in sorted(folder.rglob('*')):
        if x.is_file():z.write(x,x.relative_to(folder))
print(json.dumps({'report':str(folder/'analysis.md'),'archive':str(folder.with_suffix('.zip'))}))
