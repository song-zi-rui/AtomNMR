'use strict';
const $=id=>document.getElementById(id),esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const params=new URLSearchParams(location.search);let inherited={};let job=null,poll=null;
try{inherited=JSON.parse(params.get('filters')||'{}');}catch{}
const mid=params.get('molecule_id');if(mid)$('query').value=mid;else $('query').value=inherited.q||'';
$('origin').value=inherited.origin||'';
if(Object.keys(inherited).some(k=>!['q','origin'].includes(k))){$('inherited').hidden=false;$('inherited').textContent='已继承标注页的队列及来源等筛选条件。';}
async function api(path,body){const r=await fetch(path,body===undefined?{}:{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});const x=await r.json();if(!r.ok)throw Error(x.error||'请求失败');return x;}
function notice(text){$('message').textContent=text;$('message').hidden=!text;}
function options(){return {exclude_partial:$('excludePartial').checked,exclude_unassigned:$('excludeNone').checked,equivalence:$('equivalence').value,require_dual:$('requireDual').checked,require_same_solvent:$('sameSolvent').checked};}
function payload(){const filters={...inherited,q:$('query').value.trim(),origin:$('origin').value};const p={format:$('format').value,filters,export_options:options()};
 if(p.format==='atompeak')p.training_options={seed:Number($('seed').value),valid_ratio:Number($('validRatio').value),test_ratio:Number($('testRatio').value),max_seq_len:Number($('maxSeq').value),use_solvent:$('useSolvent').checked};return p;}
function describe(){const o=options(),parts=[];if(o.exclude_partial)parts.push('筛除部分归属');if(o.exclude_unassigned)parts.push('筛除无确定归属');if(o.equivalence!=='none')parts.push('筛除等价原子疑似问题');if(o.require_dual)parts.push('仅双谱');if(o.require_same_solvent)parts.push('已知同溶剂谱对');$('summary').textContent=parts.length?parts.join(' · '):'保留所有归属完整度，不额外筛除等价原子疑似问题。';$('trainingOptions').hidden=$('format').value!=='atompeak';}
for(const id of ['excludePartial','excludeNone','equivalence','requireDual','sameSolvent','format'])$(id).onchange=describe;
$('allPreset').onclick=()=>{$('excludePartial').checked=false;$('excludeNone').checked=false;$('requireDual').checked=false;$('sameSolvent').checked=false;$('equivalence').value='none';describe();};
$('strictPreset').onclick=()=>{$('excludePartial').checked=true;$('excludeNone').checked=true;$('requireDual').checked=true;$('equivalence').value='carbon_missing';describe();};
const reasonNames={coverage_full:'完整归属',coverage_partial:'部分归属',coverage_none:'无确定归属',no_eligible_spectra:'无当前可用谱',missing_paired_nucleus:'缺少H谱或C谱',no_same_known_solvent_pair:'没有已知同溶剂谱对',suspected_source_equivalence:'ShiftDB2等价原子疑似问题'};
function busy(value){for(const id of ['preview','run','origin','query','format','excludePartial','excludeNone','equivalence','requireDual','sameSolvent','allPreset','strictPreset','seed','maxSeq','validRatio','testRatio','useSolvent'])$(id).disabled=value;}
async function refresh(){clearTimeout(poll);if(!job)return;try{const r=await api('/api/export/'+job);if(r.state==='running'){$('progress').textContent=r.total!=null?`已检查 ${Number(r.processed).toLocaleString()} / ${Number(r.total).toLocaleString()} 个分子，筛选保留 ${Number(r.exported).toLocaleString()} 个。`:'正在建立只读数据快照…';poll=setTimeout(refresh,2000);return;}
 if(r.state==='error')throw Error(r.error);busy(false);try{localStorage.removeItem('atomnmr.export.pending');}catch{}
 $('progress').textContent=r.preview?`预览完成：可保留 ${r.count.toLocaleString()} 个分子。`:`导出完成：${r.count.toLocaleString()} 个分子。`;
 $('downloads').innerHTML=r.preview?'':`<a href="/api/export/${job}/download">下载数据文件${$('format').value==='atompeak'?' / 训练包':''}</a><a href="/api/export/${job}/manifest">下载筛选与校验清单</a>`;
 const s=r.selection;if(s){$('selection').hidden=false;$('counts').textContent=`检查 ${s.scanned.toLocaleString()} 个，保留 ${s.kept.toLocaleString()} 个，筛除 ${s.excluded.toLocaleString()} 个。`;$('reasonRows').innerHTML=Object.entries(s.reasons).map(([k,v])=>`<tr><td>${esc(reasonNames[k]||k)}</td><td>${Number(v).toLocaleString()}</td></tr>`).join('')||'<tr><td colspan="2">未命中筛除条件。</td></tr>';}
 $('note').textContent=(r.note||'筛选设置和计数已写入manifest，未修改原数据库。')+' 下方结果对应任务保存的筛选条件；改变上方设置后需重新预览或导出。';$('examples').hidden=!r.examples?.length;$('exampleRows').innerHTML=(r.examples||[]).map(x=>`<p><a href="/?molecule_id=${encodeURIComponent(x.id)}" target="_blank" rel="noopener">${esc(x.formula)} · ${esc(x.id)}</a></p>`).join('');
 }catch(e){busy(false);notice(e.message);}}
async function start(preview){notice('');busy(true);$('downloads').innerHTML='';$('selection').hidden=true;try{const r=await api('/api/export'+(preview?'/preview':''),payload());job=r.job_id;try{localStorage.setItem('atomnmr.export.pending',job);}catch{}await refresh();}catch(e){busy(false);notice(e.message);}}
$('preview').onclick=()=>start(true);$('run').onclick=()=>start(false);describe();
try{job=localStorage.getItem('atomnmr.export.pending');}catch{}if(job){busy(true);refresh();}
