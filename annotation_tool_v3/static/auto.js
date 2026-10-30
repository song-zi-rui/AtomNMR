'use strict';
const $=id=>document.getElementById(id),esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const state={job:null,status:null,poll:null,replace:null,keys:[],selectedKey:null,selectedItem:null,busy:false,activeId:null};
const names={running:'持续标注中',pausing:'正在停止本地等待',paused:'已暂停',interrupted:'服务中断 · 可恢复',attention:'需要处理错误',completed:'所选队列已处理完',pending:'等待调用',calling:'分析中',responded:'已收响应 · 待校验',received:'等待入库',done:'已存数据库',error:'失败'};
async function api(path,body){const r=await fetch(path,body===undefined?{}:{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});const x=await r.json();if(!r.ok)throw Error(x.error||'请求失败');return x;}
function message(text){const el=$('vault').open?$('vaultMessage'):$('message');el.textContent=text;el.hidden=!text;}
function action(id,fn){$(id).onclick=async()=>{const button=$(id);button.disabled=true;message('');try{await fn();}catch(e){message(e.message);}finally{button.disabled=(id==='pause'&&state.status==='pausing')||(id==='start'&&state.busy);}};}
const date=v=>v?new Date(v).toLocaleString('zh-CN'):'—';
async function keys(){const x=await api('/api/ai/credentials');state.keys=x.credentials;state.selectedKey=x.active_id;
 $('credential').innerHTML='<option value="">请选择已保存密钥</option>'+x.credentials.filter(k=>k.enabled).map(k=>`<option value="${k.id}">${esc(k.name)} · ${esc(k.masked_key)}</option>`).join('');$('credential').value=x.active_id||'';
 const current=x.credentials.find(k=>k.id===x.active_id);$('credentialHint').textContent=current?(current.enabled?'当前：':'已停用：')+current.name+' · '+current.masked_key:'首次使用请在密钥管理中添加；保存后自动沿用。';
 $('credentialRows').innerHTML=x.credentials.map(k=>`<tr><td><strong>${esc(k.name)}</strong><br><code>${esc(k.masked_key)}</code></td><td>${k.active?'当前 · ':''}${k.enabled?'可用':'已停用'}</td><td>${k.verification==='ok'?'验证成功':k.verification==='failed'?'验证失败':'尚未验证'}<br><small>验证 ${date(k.verified_at)}<br>使用 ${date(k.last_used_at)}</small></td><td>${k.enabled&&!k.active?`<button data-key="${k.id}" data-action="activate">设为当前</button>`:''}<button data-key="${k.id}" data-action="replace">替换内容</button><button data-key="${k.id}" data-action="rename">改名</button><button data-key="${k.id}" data-action="${k.enabled?'disable':'enable'}">${k.enabled?'停用':'启用'}</button><button data-key="${k.id}" data-action="remove">移除</button></td></tr>`).join('')||'<tr><td colspan="4">还没有保存密钥。</td></tr>';
}
async function settings(){const s=await api('/api/ai/settings');$('model').value=s.model;$('thinking').checked=s.thinking;$('tokens').value=s.max_tokens;$('timeout').value=s.timeout;$('effort').value=s.reasoning_effort||'high';$('wireFormat').value=s.wire_format||'compact-v2';$('assignmentMode').value=s.assignment_mode||'assisted';$('baseUrl').value=s.base_url||'https://api.deepseek.com';$('apiAdapter').value=s.api_adapter||'deepseek';$('jsonMode').checked=s.json_mode!==false;$('tlsPolicy').value=s.tls_policy||'system';$('caFile').value=s.ca_file||'';if(s.key_error)message(s.key_error);updateExecutionHint();}
async function saveConfig(){const result=await api('/api/ai/settings',{model:$('model').value.trim(),thinking:$('thinking').checked,max_tokens:Number($('tokens').value),timeout:Number($('timeout').value),reasoning_effort:$('effort').value,wire_format:$('wireFormat').value,assignment_mode:$('assignmentMode').value,base_url:$('baseUrl').value.trim(),api_adapter:$('apiAdapter').value,json_mode:$('jsonMode').checked,tls_policy:$('tlsPolicy').value,ca_file:$('caFile').value.trim()});await api('/api/auto/execution',{parallelism:Number($('parallelism').value),parallel_strategy:$('parallelStrategy').value});return result;}
async function loadProfiles(){const r=await api('/api/benchmark/profiles');state.profiles=r.profiles;$('endpointProfile').innerHTML='<option value="">选择已配置端点</option>'+r.profiles.filter(p=>p.ready).map(p=>`<option value="${p.id}">${esc(p.name)} · ${esc(p.config.model)}</option>`).join('');$('modelChoice').innerHTML='<option value="">沿用当前配置</option>'+r.profiles.filter(p=>p.ready).map(p=>`<option value="${p.id}">${esc(p.name)}</option>`).join('');}
action('adoptEndpoint',async()=>{if(!$('endpointProfile').value)throw Error('请先在测试台配置可用端点');await api('/api/auto/endpoint',{profile_id:$('endpointProfile').value});await Promise.all([keys(),settings()]);await execution();message('已采用端点配置与对应密钥；未启动任务，也未导入任何测试结果。');});
action('saveSettings',async()=>{const r=await saveConfig();message(r.saved_for_next_run?'下一轮参数已保存；当前请求仍按原设置等待，不会启动或恢复任务。密钥切换需等当前请求结束。':'配置已保存；未启动或恢复任务。');});
async function execution(){const r=await api('/api/auto/execution');$('parallelism').value=r.parallelism;$('parallelStrategy').value=r.parallel_strategy;updateExecutionHint();}
function updateExecutionHint(){
 const official=$('apiAdapter').value==='deepseek'&&/^https:\/\/api\.deepseek\.com(?::443)?(?:\/v1)?\/?$/.test($('baseUrl').value.trim())&&['deepseek-flash','deepseek-v4-pro','deepseek-v4-flash','deepseek-v4-flash-vision-exp'].includes($('model').value.trim());
 for(const value of ['deepseek_offpeak','deepseek_offpeak_only']){const option=$('parallelStrategy').querySelector?.('[value="'+value+'"]');if(option)option.disabled=!official;}
 if(!official)$('parallelStrategy').value='fixed';
 const n=Number($('parallelism').value)||1,scheduled=$('parallelStrategy').value==='deepseek_offpeak';
 $('parallelHint').textContent=$('parallelStrategy').value==='deepseek_offpeak_only'?`峰价停止领取新分子；临近谷价结束会按历史耗时提前收尾，已发请求继续完成。谷价按 ${n} 并行自动继续。本程序最多 64 并行；已发出的请求继续完成。北京时间工作日 09–12、14–18 点为峰价，周末及法定节假日为谷价。`:scheduled?`峰价 ${n} 个、谷价 16 个同时分析。北京时间工作日 09–12、14–18 点为峰价，周末及法定节假日为谷价。`:`同时分析 ${n} 个分子，每个分子联合读取 C/H 谱。`;
 $('modelHint').textContent='当前：'+$('model').value+(official?' · DeepSeek官方':'');
}
$('parallelism').oninput=updateExecutionHint;
async function saveExecution(){updateExecutionHint();try{await api('/api/auto/execution',{parallelism:Number($('parallelism').value),parallel_strategy:$('parallelStrategy').value});await refresh();}catch(e){message(e.message);await execution();}}
$('parallelism').onchange=saveExecution;$('parallelStrategy').onchange=saveExecution;
$('modelChoice').onchange=async()=>{if(!$('modelChoice').value)return;try{await api('/api/auto/endpoint',{profile_id:$('modelChoice').value});await Promise.all([keys(),settings()]);await execution();message('已切换模型配置；未启动标注。');}catch(e){message(e.message);$('modelChoice').value='';}};
$('timingPanel').ontoggle=()=>{if($('timingPanel').open)loadTiming().catch(e=>message(e.message));};

async function openVault(){await keys();$('vault').showModal();}
function resetForm(){state.replace=null;$('credentialName').value='';$('credentialSecret').value='';$('credentialFormTitle').textContent='添加密钥';}
action('manageKeys',openVault);action('keys',openVault);$('closeVault').onclick=()=>$('vault').close();$('newCredential').onclick=resetForm;
action('saveCredential',async()=>{await api('/api/ai/credentials',{action:state.replace?'replace':'create',id:state.replace,name:$('credentialName').value.trim(),api_key:$('credentialSecret').value.trim(),activate:true});resetForm();await keys();message('已加密保存并设为当前，后续自动沿用。');});
action('testKey',async()=>{try{await api('/api/ai/test',{});message('当前密钥验证成功。');}finally{await keys();}});
$('credentialRows').onclick=async e=>{const b=e.target.closest('[data-key]');if(!b)return;const key=state.keys.find(x=>x.id===b.dataset.key),kind=b.dataset.action;if(kind==='replace'){state.replace=key.id;$('credentialName').value=key.name;$('credentialSecret').value='';$('credentialFormTitle').textContent='替换：'+key.name;$('credentialSecret').focus();return;}
 let name;if(kind==='rename'){name=prompt('新的密钥名称',key.name);if(!name)return;}if(kind==='remove'&&!confirm('仅移除此系统中保存的“'+key.name+'”？历史标注来源保留。'))return;
 b.disabled=true;try{await api('/api/ai/credentials',{action:kind,id:key.id,name});await keys();}catch(e){message(e.message);}finally{b.disabled=false;}};
$('credential').onchange=async()=>{if(!$('credential').value)return;try{await api('/api/ai/credentials',{action:'activate',id:$('credential').value});await keys();}catch(e){message(e.message);$('credential').value=state.selectedKey||'';}};
action('start',async()=>{await saveConfig();const r=await api('/api/auto/jobs',{scope:$('scope').value});state.job=r.job_id;$('analysis').hidden=true;await refresh();});
action('pause',async()=>{await api('/api/auto/jobs/'+state.job+'/pause',{});await refresh();});
action('resume',async()=>{await saveConfig();await api('/api/auto/jobs/'+state.job+'/resume',{use_current_limits:true});await refresh();});
async function retry(job,itemId,mode='failed'){if(mode!=='local')await saveConfig();const body=itemId?{item_id:itemId}:{};if(mode!=='failed')body.mode=mode;const r=await api('/api/auto/jobs/'+job+'/retry',body);state.job=r.job_id;state.selectedItem=null;$('detail').hidden=true;$('analysis').hidden=true;await refresh();message(mode==='local'?'已建立本地恢复任务，不调用 API；原响应和原任务保留。':'已建立独立重试任务，仅包含 '+r.selected+' 个分子；原队列进度保持不变。');}
action('retryFailed',()=>retry(state.job));
action('retryOne',async()=>{if(!state.selectedItem)throw Error('请先打开一个失败分子');await retry(state.selectedItem.job_id,state.selectedItem.id);});
action('recoverOne',async()=>{if(!state.selectedItem)return;await retry(state.selectedItem.job_id,state.selectedItem.id,'local');});
action('retryIssues',async()=>{if(!state.selectedItem)return;await retry(state.selectedItem.job_id,state.selectedItem.id,'issues');});
action('retryUncertain',async()=>{if(!state.selectedItem)return;await retry(state.selectedItem.job_id,state.selectedItem.id,'uncertain');});
async function loadTiming(){const r=await api('/api/auto/timing');$('timingNote').textContent=r.note;const f=n=>n==null?'—':n<60?n.toFixed(1)+'秒':Math.floor(n/60)+'分'+(n%60).toFixed(0)+'秒';$('timingRows').innerHTML=r.jobs.map(j=>`<tr><td>${j.scope==='retry'?'有限重试':'持续队列'}<br><code>${j.id.slice(0,8)}</code></td><td>${j.all.n}</td><td>${f(j.all.total_s)}</td><td>${f(j.all.mean_s)} / ${f(j.all.median_s)}</td><td>${f(j.all.p90_s)}</td><td>${f(j.all.max_s)}</td></tr>`).join('');}
action('timingLoad',loadTiming);
function showPauseNotice(job){
 const q=job.current_request,waiting=job.state==='pausing';$('pauseNotice').hidden=!q&&!waiting;if(!q&&!waiting)return;
 let text=q?`${q.formula} · 已等待 ${Math.floor(q.elapsed_s/60)}分${Math.floor(q.elapsed_s%60)}秒；本次等待时限 ${q.timeout_s}秒。`:'正在保存当前请求结果。';
 if(q?.stream){const p=q.progress||{};text+=p.ttft_ms!=null?`已收到 ${p.chunks||0}个数据块。`:'正在等待接口首个输出。';
  if(p.reasoning_characters>0&&!p.content_characters)text+='模型正在思考，归属答案尚未开始返回；思考输出仍在接收。';}
 else if(q)text+='当前为非流式调用，完整结果返回前不会显示输出。';
 if(waiting)text+='正在停止本地等待，不再领取新分子。远端是否停止计算取决于API；未确认的当前分子会留在中断记录中，不自动重试。';
 const requests=job.current_requests||[];if(requests.length>1)text+='\n其他正在分析：'+requests.slice(1).map(r=>r.formula+'（'+Math.floor(r.elapsed_s)+'秒）').join('、');
 $('pauseNotice').textContent=text;
}
async function refresh(){clearTimeout(state.poll);const params=new URLSearchParams();if(state.job)params.set('id',state.job);if($('resultFilter').value)params.set('result_state',$('resultFilter').value);const query=params.toString();const r=await api('/api/auto/jobs'+(query?'?'+query:''));const j=r.job;state.activeId=r.active_id;state.busy=!!r.active_id||['running','pausing'].includes(j?.state);
 $('jobs').innerHTML=r.jobs.map(x=>`<div class="job-row ${x.id===state.job?'selected':''}"><button data-job="${x.id}">${esc(names[x.state]||x.state)}<small>${date(x.created_at)} · ${x.id.slice(0,8)}</small></button><button data-delete="${x.id}" class="delete-job" aria-label="删除任务记录 ${x.id.slice(0,8)}" title="删除任务记录，保留已入库归属" ${x.id===r.active_id||['running','pausing'].includes(x.state)?'disabled':''}>×</button></div>`).join('')||'<p class="hint">尚无任务记录。</p>';
 if(!j){state.job=null;$('jobTitle').textContent='标注进度';$('status').textContent='尚未运行';$('counters').innerHTML='';$('items').innerHTML='<tr><td colspan="5">尚无处理记录。</td></tr>';$('parallelStatus').textContent='';$('jobConfig').textContent='';$('jobError').textContent='';$('usageCards').innerHTML='';$('usageHint').textContent='';$('pauseNotice').hidden=true;for(const id of ['pause','resume','retryFailed','analyze'])$(id).hidden=true;$('start').disabled=state.busy;return;}
state.job=j.id;state.status=j.state;$('status').textContent=names[j.state]||j.state;$('jobTitle').textContent='任务 '+j.id.slice(0,8);$('jobConfig').textContent=j.config.model+' · 实际任务输出上限 '+j.config.max_tokens+' · '+(j.config.thinking?'思考 '+(j.config.reasoning_effort||'high'):'非思考模式')+' · '+(j.config.wire_format||'verbose-v1')+' · '+(j.config.assignment_mode==='assisted'?'生产补全':'盲标')+' · '+(j.scope==='retry'?(j.config.retry_mode==='local'?'本地恢复 ':'有限重试 ')+j.config.retry_item_ids.length+' 个分子':'持续处理');
 const ex=j.execution||{},n=ex.effective_parallelism??j.config.parallelism??1;$('parallelStatus').textContent=['paused','attention','interrupted'].includes(j.state)?'恢复时应用上方并行设置，已完成的分子不重复分析。':j.state==='completed'?'本任务已完成。':`同时分析 ${j.counts?.calling||0} / ${n} 个分子`+(ex.schedule_active?' · '+(ex.draining_before_peak?'谷价收尾：不再派发新分子，完成在途请求':ex.tariff==='offpeak'?'谷价运行':(n===0?'峰价等待，谷价自动继续':'峰价运行')):'');
 if(ex.duration_estimate&&$('parallelStrategy').value==='deepseek_offpeak_only'){const d=ex.duration_estimate;$('parallelHint').textContent=`仅谷价运行 · 已统计 ${d.samples} 次同配置成功调用`+(d.mean_seconds!=null?`，平均 ${Math.round(d.mean_seconds)} 秒/分子`:'')+`；结束前约 ${Math.ceil(d.guard_seconds/60)} 分钟停止派发，已有请求继续完成。`;}
 const c=j.counts||{},metrics=[['已处理分子',j.processed],['成功入库结果',c.done||0],['新增 AI 归属峰',j.applied_peaks],['失败分子',c.error||0]];
 $('counters').innerHTML=metrics.map(([label,value])=>`<div class="metriccard"><span>${label}</span><strong>${Number(value||0).toLocaleString()}</strong></div>`).join('');
 const u=j.usage||{},um=[['输入 token',u.input],['内部思考 token',u.reasoning],['最终输出 token',u.final],['总 token',u.total]];
 $('usageCards').innerHTML=um.map(([label,value])=>`<div class="metriccard"><span>${label}</span><strong>${Number(value||0).toLocaleString()}</strong></div>`).join('');
 $('usageHint').textContent='用量已回报 '+(u.reported||0)+' 次；未回报用量不计作零费用。缓存命中输入 '+Number(u.cache_hit||0).toLocaleString()+' token。'+((u.reasoning_reported||0)<(u.reported||0)?'部分历史响应未单独回报思考用量，分项仅显示已知部分。':'');
 showPauseNotice(j);const running=['running','pausing'].includes(j.state);$('start').disabled=state.busy;$('pause').hidden=!running;$('pause').disabled=j.state==='pausing';$('resume').hidden=running||j.state==='completed';$('analyze').hidden=!j.processed;$('jobError').textContent=j.error||'';
 $('retryFailed').hidden=running||!c.error;
 const requests=new Map((j.current_requests??(j.current_request?[j.current_request]:[])).map(q=>[String(q.item_id),q]));
 $('items').innerHTML=j.items.map(i=>{const request=requests.get(String(i.id)),elapsed=i.elapsed??(i.state==='calling'?request?.elapsed_s:null),time=Number.isFinite(elapsed)?elapsed.toFixed(1)+' s'+(i.elapsed==null?'（进行中）':''):'—';return `<tr data-item="${i.id}"><td>${esc(i.formula)}</td><td><code>${esc(i.molecule_id)}</code></td><td>${esc(names[i.state]||i.state)}${i.error?'<br>'+esc(i.error):''}</td><td>${i.publication.applied??'—'}</td><td>${time}</td></tr>`;}).join('')||'<tr><td colspan="5">此筛选暂无处理记录。</td></tr>';
 if(running)state.poll=setTimeout(()=>refresh().catch(e=>message(e.message)),2500);
}
$('jobs').onclick=async e=>{const del=e.target.closest('[data-delete]');if(del){if(!confirm('从任务列表删除此记录？已入库归属和来源保留，不会重新标注。'))return;del.disabled=true;try{await api('/api/auto/jobs/'+del.dataset.delete+'/delete',{});if(state.job===del.dataset.delete)state.job=null;state.selectedItem=null;$('detail').hidden=true;$('analysis').hidden=true;await refresh();message('已删除任务记录，已入库归属保留。');}catch(error){message(error.message);}finally{del.disabled=false;}return;}const b=e.target.closest('[data-job]');if(b){state.job=b.dataset.job;state.selectedItem=null;$('detail').hidden=true;$('analysis').hidden=true;refresh().catch(e=>message(e.message));}};

action('refresh',refresh);
$('resultFilter').onchange=()=>refresh().catch(e=>message(e.message));
action('analyze',async()=>{const r=await api('/api/auto/jobs/'+state.job+'/analyze',{}),pct=x=>x==null?'—':(x*100).toFixed(1)+'%';$('analysis').hidden=false;$('analysisNote').textContent=r.note;$('metrics').innerHTML=r.metrics.map(m=>`<tr><td>${m.source==='manual'?'人工':'ShiftDB2'}</td><td>${m.nucleus==='all'?'H + C':m.nucleus}</td><td>${m.reference_peaks}</td><td>${pct(m.accuracy)}</td><td>${pct(m.coverage)}</td><td>${pct(m.top3_rate)}</td></tr>`).join('');});
$('items').onclick=async e=>{const b=e.target.closest('[data-item]');if(!b)return;try{const r=await api('/api/auto/items/'+b.dataset.item);state.selectedItem={id:r.id,job_id:r.job_id};$('detail').hidden=false;$('detailTitle').textContent=r.formula+' · AI 归属';$('openMolecule').href='/?molecule_id='+r.molecule_id;$('structure').src='/api/molecule/'+r.molecule_id+'/svg';$('summary').textContent=r.output?.summary||r.error||'尚未返回';$('warnings').textContent=(r.output?.warnings||[]).join('\n')+((r.output?.consistency_issues||[]).length?'\n本地一致性检查：存在原子复用或分组问题，相关AI建议保留供查看，不自动补入确定归属。':'');
 const usage=r.usage||{},diag=r.diagnostics||{},reasoning=usage.completion_tokens_details?.reasoning_tokens;
 $('detailUsage').textContent='本次上限 '+(diag.max_tokens??r.request_config?.max_tokens??'未记录')+'；实际强度 '+(diag.reasoning_effort||'历史未单列')+'；终止原因 '+(diag.finish_reason||'未返回')+'；输入 '+(usage.prompt_tokens??'未返回')+'，生成 '+(usage.completion_tokens??'未返回')+' token'+(reasoning!=null?'，其中思考 '+reasoning+'、其余输出 '+(usage.completion_tokens-reasoning):'')+'。'+(diag.wire_characters!=null?' 紧凑输入 '+diag.wire_characters+' / 原表示 '+diag.input_characters+' 字符。':'')+(diag.response_reused_from?' 本地复用响应，未重新调用 API。':'');
 $('retryOne').hidden=r.state!=='error';
 $('recoverOne').hidden=r.state!=='error'||r.raw_response?.diagnostics?.finish_reason!=='stop';
 $('retryIssues').hidden=r.state!=='done'||!r.output?.repair_targets?.length;
 const wanted=r.output?.requested_peaks||r.input.requested_peaks;
 $('retryUncertain').hidden=r.state!=='done'||!(r.output?.assignments||[]).some(p=>['ambiguous','unassigned'].includes(p.status)&&(!wanted||wanted.some(t=>t.spectrum_id===p.spectrum_id&&t.peak_id===p.peak_id)));
 $('rawJSON').textContent=r.raw_response?.content||'此历史记录没有保存原始最终响应。';
 const actions={ambiguous_atoms_moved_to_candidates:'歧义原子移到候选字段',candidates_limited_to_top3:'保留前三组有序候选',quarantined_invalid_peak:'该峰未通过校验，已留空',missing_peak:'遗漏峰，已留空',ignored_not_requested:'未请求的峰，未采用'};
 $('validationIssues').textContent=(r.output?.validation_issues||[]).map(x=>x.peak_id+'：'+(actions[x.action]||x.action)).join('；');
 const atoms=r.input.atoms,targets=x=>(x||[]).map(i=>atoms[i].element+i).join(', ')||'—';
 $('predictions').innerHTML='<div class="tablewrap"><table><thead><tr><th>谱峰</th><th>归属状态</th><th>AI 原子建议</th><th>AI 简要依据</th></tr></thead><tbody>'+(r.output?.assignments||[]).map(p=>`<tr><td>${esc(p.peak_id)}<br><small>${esc(p.spectrum_id)}</small></td><td>${esc(p.status)}</td><td>${esc(targets(p.atom_indices))}${p.candidates?.length?'<br><small>候选 '+esc(p.candidates.map(targets).join(' / '))+'</small>':''}</td><td>${esc(p.brief_reason)}</td></tr>`).join('')+'</tbody></table></div>';
 $('detail').scrollIntoView({behavior:'smooth'});
 }catch(e){message(e.message);}};
Promise.all([keys(),settings(),refresh(),loadProfiles()]).then(async()=>{await execution();if(location.hash==='#keys')return openVault();}).catch(e=>message(e.message));
const requestedMolecule=new URLSearchParams(location.search).get('molecule_id');
if(requestedMolecule){(async()=>{const history=await api('/api/auto/results?molecule_id='+encodeURIComponent(requestedMolecule));if(!history.results.length)return;const r=await api('/api/auto/results/'+encodeURIComponent(history.results[0].id));$('detail').hidden=false;$('detailTitle').textContent=r.provider+' / '+r.model+' · 已入库结果';$('openMolecule').href='/?molecule_id='+r.molecule_id;$('structure').src='/api/molecule/'+r.molecule_id+'/svg';$('summary').textContent=r.output.summary;$('warnings').textContent=(r.output.warnings||[]).join('\n');$('predictions').innerHTML='<p class="hint">完整结果已归档。补入 '+r.applied_count+' 个空缺归属；已有标签保留。此记录的入库时间：'+date(r.created_at)+'</p><pre>'+esc(JSON.stringify(r.output.assignments,null,2))+'</pre>';$('detail').scrollIntoView({behavior:'smooth'});})().catch(e=>message(e.message));}
