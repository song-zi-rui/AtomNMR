"use strict";
const $ = selector => document.querySelector(selector);
const $$ = selector => [...document.querySelectorAll(selector)];
const labels = {complete_remaining:"部分归属",both_unassigned:"无归属",waiting_spectrum:"缺少实验谱",complete:"全归属",assigned:"唯一归属",equivalent_group:"等价组",overlap:"重叠",source_group:"原库多原子归属",ambiguous:"歧义",unassigned:"未归属",solvent:"溶剂峰",impurity:"杂质",excluded:"不计入"};
const coverageLabels={full:'全归属',partial:'部分归属',none:'无归属'};
function coverageLabel(r){const total=(r.c_count||0)+(r.h_count||0),assigned=(r.c_assigned||0)+(r.h_assigned||0);return coverageLabels[r.annotation_coverage?.status||(assigned?(assigned===total?'full':'partial'):'none')];}
labels.quality_hold="谱图已隔离";
const state = {record:null,nucleus:"13C",sid:null,pid:null,rows:[],offset:0,dirty:false,generation:0,saving:null,timer:null,loading:false,listRequest:0,navigating:null,anchor:null,anchorFilters:null,backStack:[]};
const escape = value => String(value??"").replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
async function api(path,body,options={}){
  const controller=new AbortController();let timedOut=false;
  const abort=()=>controller.abort();
  if(options.signal?.aborted)abort();else options.signal?.addEventListener("abort",abort,{once:true});
  const timer=setTimeout(()=>{timedOut=true;abort();},options.timeoutMs??(body===undefined?15000:path.endsWith("/management")?120000:30000));
  try{
    const response=await fetch(path,{signal:controller.signal,...(body===undefined?{}:{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(body)})});
    const value=await response.json();
    if(!response.ok){const error=new Error(value.error||response.statusText);error.status=response.status;error.retryable=value.retryable===true;throw error;}
    return value;
  }catch(error){
    if(error.status)throw error;
    if(options.signal?.aborted&&!timedOut)throw error;
    const failure=new Error(timedOut?(body===undefined?"请求超时，请稍后重试。":"提交响应超时，操作结果尚未确认。"):error.message||"连接已中断。");
    failure.deliveryUnknown=body!==undefined;throw failure;
  }finally{clearTimeout(timer);options.signal?.removeEventListener("abort",abort);}
}
function toast(text){
  const notification=$('#toast'),dialog=$$('dialog[open]').at(-1),host=dialog||document.body;
  // Native modal dialogs are in the top layer; body z-index cannot cover them.
  if(host?.appendChild&&notification.parentElement!==host)host.appendChild(notification);
  if(notification.style)notification.style.top=dialog?`${dialog.scrollTop+12}px`:'';
  if(dialog?.id==='managementDialog')$('#manageMessage').textContent=text;
  notification.textContent=text;notification.classList.add('visible');
  clearTimeout(toast.timer);toast.timer=setTimeout(()=>notification.classList.remove('visible'),5000);
}
function filters(){const q=$("#search").value.trim();if(q&&$("#searchScope").value!=="filtered")return {q,queue:"",origin:"",pending:"",view:""};return {q,queue:$("#queue").value,origin:$("#origin").value,pending:$("#processing").value==='pending'?"1":"",view:$("#view").value,...($('#processing').value==='finished'?{finished:'1'}:{}),...($('#coverage').value?{coverage:$('#coverage').value}:{}),...($('#annotationSource')?.value?{annotation_source:$('#annotationSource').value}:{})};}
function dateLabel(value){return value?new Date(value).toLocaleString("zh-CN",{hour12:false}):"时间未记录";}
function activityLabel(action){return {save:"保存标注 / 备注",finish:"结束标注",historical_annotation:"历史人工标注",historical_management:"历史数据管理"}[action]||(action?.startsWith("curation:")?"数据管理操作":"");}
let listController;
async function loadList(reset=false){
  if(reset)state.offset=0;
  const request=++state.listRequest,selection=filters(),offset=state.offset;
  listController?.abort();listController=new AbortController();
  $("#searchHint").textContent=selection.q&&$("#searchScope").value!=="filtered"?"全库检索中：忽略下方视图、队列、来源及完成状态。":"正在读取当前筛选…";
  let result;
  try{result=await api("/api/molecules?"+new URLSearchParams({...selection,offset,limit:40}),undefined,{signal:listController.signal});}
  catch(error){if(request!==state.listRequest)return;$("#searchHint").textContent="检索失败，请重试："+error.message;throw error;}
  if(request!==state.listRequest)return;
  state.rows=result.rows;$("#pageNo").textContent=String(1+offset/40);$("#prevPage").disabled=offset===0;$("#nextPage").disabled=result.rows.length<40;
  const hints={recent:"按最近保存标注、备注或数据管理操作排序。",finished:"已结束分子；结束时当前实验谱自动设为完整，归属覆盖另行记录。",imported:"按首次入库时间倒序。"};
  $("#searchHint").textContent=(selection.q&&$("#searchScope").value!=="filtered"?"全库检索：已忽略下方筛选。":hints[selection.view]||"结果受当前队列、来源与完成状态筛选。")+(result.total!=null?` 共 ${result.total.toLocaleString()} 个分子 ·`:"")+` 本页 ${result.rows.length} 条。`;
  $("#list").innerHTML=result.rows.length?result.rows.map(r=>`<button class="listitem ${state.record?.id===r.id?"active":""}" data-mid="${r.id}"><strong>${escape(r.formula)}</strong><span class="done-tag">${r.task_done?'本轮已结束':'尚未结束'}</span><small>C ${r.c_assigned}/${r.c_count} · H ${r.h_assigned}/${r.h_count}</small><small>${escape(coverageLabel(r))} · ${escape(labels[r.queue])===escape(coverageLabel(r))?'':escape(labels[r.queue])+' · '}${escape(r.origin)}</small><small>${selection.view==="imported"?"入库 "+dateLabel(r.created_at):r.last_processed_at?escape(activityLabel(r.last_action))+" · "+dateLabel(r.last_processed_at):""}</small><code>${escape(r.id)}</code><code>${escape(r.smiles)}</code></button>`).join(""):"<p class='hint'>当前筛选没有结果。输入编号时可选“全库检索”；查看已保存的分子可切换“最近处理”。</p>";
}
function recordCursor(record,selection){const view=selection.view||"";return [view,view==="imported"?record.created_at:["recent","finished"].includes(view)?record.activity?.last_processed_at||"":record.navigation_rank??({both_unassigned:0,complete_remaining:1,waiting_spectrum:2}[record.queue]??3),record.id];}
async function goNext(done=false,direction=1){
  if(state.navigating)return state.navigating;
  if(!state.record||state.loading)return;
  state.navigating=(async()=>{
    const selection=filters(),signature=JSON.stringify(selection),mid=state.record.id,originOffset=state.offset;
    const anchor=state.anchorFilters===signature&&state.anchor?.[2]===mid?state.anchor:recordCursor(state.record,selection);
    const previous=direction<0?state.backStack.at(-1):null;
    $("#prevTask").disabled=true;$("#nextTask").disabled=true;$("#completeNext").disabled=true;
    try{
      await save(done);
      if(done&&!state.record.task_done){toast("保存期间有新修改，尚未结束此分子；请确认后再按完成。");return;}
      if(signature!==JSON.stringify(filters())){toast("已保存；筛选已改变，请从新列表选择分子。");return;}
      const result=previous?{rows:[previous]}:await api("/api/molecules?"+new URLSearchParams({...selection,[direction<0?"before":"after"]:JSON.stringify(anchor),limit:1}));
      if(signature!==JSON.stringify(filters())||state.record.id!==mid)return;
      const next=result.rows[0];
      if(!next){toast(direction<0?"已保存；没有可返回的上一个分子。":done?"此分子已结束；当前筛选顺序下已没有下一个。":"已保存；当前筛选顺序下已没有下一个。");return;}
      // The server seeks past the captured pre-save position, even if finishing
      // moved the molecule to another queue. Never wrap to the first list row.
      if(previous&&previous.signature===signature)state.offset=previous.offset;
      else if(!state.rows.some(r=>r.id===next.id))state.offset=direction<0?Math.max(0,state.offset-40):state.rows.length===40?state.offset+40:0;
      await openRecord(next.id,direction>0,originOffset);
      if(previous)state.backStack.pop();
      if(next.cursor&&(!previous||previous.signature===signature)){state.anchor=next.cursor;state.anchorFilters=signature;}
    }finally{$("#prevTask").disabled=false;$("#nextTask").disabled=false;$("#completeNext").disabled=false;}
  })().finally(()=>{state.navigating=null;if(window.renderCurationSummary&&state.record)window.renderCurationSummary(state.record);});
  return state.navigating;
}
function goPrevious(){return goNext(false,-1);}
function rememberRecord(record,offset=state.offset){if(!record)return;state.backStack.push({id:record.id,cursor:state.anchor,signature:state.anchorFilters,offset});if(state.backStack.length>1000)state.backStack.shift();}
function renderAssignmentProgress(p){
  if(!p){$("#assignmentProgress").textContent="双谱归属进度暂不可用";return;}
  const expanded=$("#assignmentProgress").querySelector?.('details')?.open;
  $("#assignmentProgress").innerHTML=`<details class="header-disclosure progress-help" ${expanded?"open":""}><summary><strong>${p.ready_molecules.toLocaleString()}</strong> C/H 双谱全归属 <span>· ${p.percent.toFixed(2)}%</span></summary><div class="header-popover"><h2>双谱归属进度</h2><p><b>${p.ready_molecules.toLocaleString()}</b> / ${p.total_dual_molecules.toLocaleString()} 个双谱分子 · ${p.percent.toFixed(2)}%</p><progress aria-label="双谱全归属进度" max="${Math.max(1,p.total_dual_molecules)}" value="${p.ready_molecules}"></progress><p>还剩 ${p.remaining_molecules.toLocaleString()} 个双谱分子未达到全归属。</p><p>按分子去重；多条当前谱全部有归属才计入；含原库归属，不仅是人工新增。未处理、歧义及排除峰不算全归属；不代表谱完整、归属正确或 H/C 条件一致。</p></div></details>`;
}
let statsRequest;
function loadStats(){if(!statsRequest)statsRequest=loadStatsOnce().finally(()=>{statsRequest=null;});return statsRequest;}
async function loadStatsOnce(){
  // A cold database can take longer to count; this background request must not
  // share the short deadline of interactive searches or block annotation.
  const s=await api("/api/stats",undefined,{timeoutMs:60000});
  const remaining=(s.candidates.indexed||0)+(s.candidates.awaiting_source||0);
  const job=s.migration_job||{};
  const heading={running:"后台迁移进行中",generating_report:"后台正在生成迁移报告",failed:"后台迁移已中断",unconfirmed:"后台迁移状态未确认（未收到近期运行信号）",finished:remaining?"本轮处理已结束，仍有候选未导入":"后台迁移已结束"}[job.state]||"尚有候选未导入（未确认后台任务在运行）";
  const phase={audit_complements:"检查原补谱",rematch_missing:"重新搜索缺失谱",audit_imported:"检查已入库 NMRexp",recover_remaining:"回源并导入剩余候选",integrity_report:"校验数据库并生成统计"}[job.phase]||"";
  const progress=job.state==="running"&&phase?`${phase}：${(job.processed||0).toLocaleString()} / ${(job.total||0).toLocaleString()}。${job.phase_eta_seconds!=null?"本阶段预计还需约 "+Math.ceil(job.phase_eta_seconds/60)+" 分钟（不含后续阶段）。":""}`:"";
  $("#migrationStatus").hidden=!remaining&&!["running","failed","unconfirmed","generating_report"].includes(job.state);
  $("#migrationStatus").textContent=`${heading}。${progress}剩余 ${remaining.toLocaleString()} 个候选待回源 / 导入。${job.state==="failed"?"失败原因："+(job.error||"未知")+"。":""}已隔离和未完成新版检查的 NMRexp 谱不进入导出；自动检查通过不等于已核对原论文。`;
}
function spectrum(){return state.record?.spectra.find(s=>s.id===state.sid);}
function peak(){return spectrum()?.peaks.find(p=>p.peak_id===state.pid);}
function canEditMolecule(){return (state.record?.curation?.molecule_state||'active')==='active';}
function canEditSpectrum(s=spectrum()){return canEditMolecule()&&!!s&&s.editable!==false&&(!s.manual_state||s.manual_state==='active')&&s.quality?.decision!=='quarantined';}
function atomEligible(a){return !!a&&(state.nucleus==='13C'?a.symbol==='C':a.total_hydrogens>0||(peak()?.source_assignment_candidates||[]).includes(a.index));}
function spectrumStateLabel(s){return s.quality?.decision==='quarantined'?'自动隔离 · 只读':({held:'暂时停用 · 不导出',excluded:'排除训练 · 不导出',reference:(s.version_label||'参考版本')+' · 只读 / 不导出',superseded:(s.version_label||'参考版本')+' · 只读 / 不导出'}[s.manual_state]||(!canEditMolecule()?'分子已停用 · 只读':s.version_info?'当前数据集 · 最终版':'当前数据集 · 原始谱'));}
function spectrumOrder(a,b){return Number(!canEditSpectrum(a))-Number(!canEditSpectrum(b))||Number(!a.version_info)-Number(!b.version_info)||(b.version_info?.number||0)-(a.version_info?.number||0);}
function renderStructureNotice(){
  const r=state.record,metals=r?.structure.metal_atoms||[];
  $('#structureView').hidden=!metals.length;$('#complexNotice').hidden=!metals.length;
  if(!metals.length)return;
  const atom=i=>(r.structure.atoms[i]?.symbol||'?')+i;
  const links=(r.structure.metal_connections||[]).map(b=>b.atoms.map(atom).join(' ↔ ')+' ('+b.bond_type+')');
  const nc=r.spectra.filter(s=>s.nucleus==='13C').map(s=>s.peaks.length).join(' / ')||'未收录';
  const nh=r.spectra.filter(s=>s.nucleus==='1H').map(s=>s.peaks.length).join(' / ')||'未收录';
  $('#complexNotice').innerHTML='<p>含金属 '+metals.map(atom).join('、')+'；'+r.structure.atoms.filter(a=>a.symbol==='C').length+' 个碳原子。原始记录 C 位移项：'+nc+'；H 位移项：'+nh+'。峰数少可能来自等价、重叠或缺项，不能自动补齐。</p><p>配位键编码可能影响推算氢数和分子式，须查原文确认。原库明确指向的 H 归属会另行提示，不自动替你填写。</p><p>配体展开只在绘图副本中隐藏金属连接，原子编号不变；不表示游离配体、三维构型或重新确定配位结构。可切换“完整连接”对照。</p><details><summary>原结构中的金属连接</summary><p>'+escape(links.join('；'))+'</p></details>';
}
function dirty(p=null){state.dirty=true;state.generation++;if(p)p._editGeneration=state.generation;state.record.task_done=0;$("#saveState").textContent=state.saving?"有新修改 · 等待上一笔保存":"有修改 · 即将保存";clearTimeout(state.timer);state.timer=setTimeout(()=>save().catch(e=>toast(e.message)),900);}
function confirmsSave(record,id,payload){
  if(record.id!==id||record.revision!==payload.revision+1||record.note!==payload.note||Boolean(record.task_done)!==Boolean(payload.task_done))return false;
  return payload.assignments.every(saved=>{
    const current=record.spectra.find(s=>s.id===saved.spectrum_id)?.peaks.find(p=>p.peak_id===saved.peak_id)?.assignment;
    return current&&current.status===saved.status&&(current.note||"")===saved.note&&JSON.stringify([...current.atom_indices].sort((a,b)=>a-b))===JSON.stringify([...saved.atom_indices].sort((a,b)=>a-b));
  });
}
async function performSave(done){
  const generation=state.generation;
  const record=state.record;
  // Freeze the request, including atom arrays, before any retry / later edit.
  const payload=JSON.parse(JSON.stringify({revision:record.revision,note:record.note||"",task_done:done,assignments:record.spectra.flatMap(s=>s.peaks.filter(p=>p._dirty).map(p=>({spectrum_id:s.id,peak_id:p.peak_id,status:p.assignment.status,atom_indices:p.assignment.atom_indices,note:p.assignment.note||""})))}));
  $("#save").disabled=true;$("#complete").disabled=true;
  $("#saveState").textContent="正在保存…";
  try{
    let result;
    for(let attempt=0;attempt<8;attempt++){
      try{result=await api("/api/molecule/"+record.id,payload);break;}
      catch(error){
        if(error.deliveryUnknown){
          $("#saveState").textContent="响应中断 · 正在核对数据库中的保存结果…";
          try{const confirmed=await api("/api/molecule/"+record.id);if(confirmsSave(confirmed,record.id,payload)){result=confirmed;break;}}catch{}
          throw error;
        }
        // Retry only an explicit, rolled-back lock conflict. A dropped network
        // reply may already have committed and must not be replayed blindly.
        if(error.status!==503||!error.retryable||attempt===7)throw error;
        $("#saveState").textContent=`等待后台写入 · 即将重试 ${attempt+1}/7（修改仍保留）`;
        await new Promise(resolve=>setTimeout(resolve,300*Math.min(attempt+1,4)));
      }
    }
    record.revision=result.revision;record.queue=result.queue;
    if(generation===state.generation){
      state.record=result;state.dirty=false;
      $("#saveState").textContent="已保存";
    }else{
      // Do not clear or replace edits made while the request was in flight.
      record.task_done=0;
      for(const s of record.spectra)for(const p of s.peaks){
        if(p._dirty&&(p._editGeneration??0)<=generation)p._dirty=false;
      }
      $("#saveState").textContent="上一笔已保存 · 即将保存新修改";
      clearTimeout(state.timer);
      state.timer=setTimeout(()=>save().catch(e=>toast(e.message)),300);
    }
    renderHeader();
    if(done){
      toast(state.dirty?"上一笔已保存；新增修改仍待保存，分子尚未结束。":"已结束此分子；歧义记录原样保留。");
    }
    if(done||["recent","finished"].includes(filters().view)){
      // A list refresh error must not turn a successful save into a failure.
      void loadList().catch(error=>toast("归属已保存，但列表刷新失败："+error.message));
    }
  }catch(error){$("#saveState").textContent=error.deliveryUnknown?"保存结果未确认 · 修改仍保留，请勿关闭或刷新":"保存失败 · 修改仍在当前页面，请勿关闭或刷新";throw error;}
  finally{$("#save").disabled=false;$("#complete").disabled=false;}
}
async function save(done=false){
  clearTimeout(state.timer);
  if(state.saving){await state.saving;return save(done);}
  if(!state.record||(!state.dirty&&(!done||state.record.task_done)))return;
  if(!canEditMolecule())throw new Error('当前分子为只读存档；可浏览谱峰，恢复分子后再保存标注。');
  // Serialize the whole acknowledgement handling, not just the HTTP request.
  state.saving=performSave(done).finally(()=>{state.saving=null;});
  await state.saving;
  // Navigation/export await save(): they must not leave edits made in flight.
  if(state.dirty)return save();
}
async function openRecord(id,remember=true,historyOffset=state.offset){if(state.loading)return;state.loading=true;try{await save();const loaded=await api("/api/molecule/"+id);if(remember&&state.record?.id!==id)rememberRecord(state.record,historyOffset);state.record=loaded;const selection=filters();state.anchor=recordCursor(state.record,selection);state.anchorFilters=JSON.stringify(selection);state.pid=null;$("#atomNumber").value="";state.dirty=false;const choices=[...state.record.spectra].sort(spectrumOrder);const current=choices.filter(canEditSpectrum);const pending=current.find(s=>s.peaks.some(p=>p.assignment.status==="unassigned"))||current[0]||choices[0];state.nucleus=pending?.nucleus||"13C";state.sid=pending?.id;$("#empty").hidden=true;$("#workspace").hidden=false;$("#recordNote").value=state.record.note||"";$("#saveState").textContent="已加载";renderHeader();renderStructure();renderSpectra();void loadList().catch(e=>toast("分子已加载，列表刷新失败："+e.message));}finally{state.loading=false;}}
function renderHeader(){
  const r=state.record;if(!r)return;
  $("#formula").textContent=r.formula;$("#recordId").textContent=r.id;
  $("#recordBadges").innerHTML=`<span class="badge">${escape(r.origin)}</span><span class="badge" title="按已保存的归属状态统计，不代表正确性认证">${escape(coverageLabel(r))}</span><span class="badge">${r.task_done?'本轮标注已结束':'尚未结束本轮标注'}</span>`;
  if(r.ai_history?.length)$('#recordBadges').innerHTML+=`<a class="badge" href="/auto.html?molecule_id=${r.id}" target="_blank" rel="noopener">查看 AI 归属记录 · ${escape([...new Set(r.ai_history.map(x=>x.model))].join(' / '))}</a>`;
  if(window.renderCurationSummary)window.renderCurationSummary(r);
}
function renderStructure(){const r=state.record;if(!r)return;const selected=peak()?.assignment.atom_indices||[];$("#molecule").src="/api/molecule/"+r.id+"/svg?"+new URLSearchParams({atoms:selected.join(","),view:$("#structureView").value||"auto"});renderStructureNotice();$("#smiles").textContent=r.smiles;$("#atoms").innerHTML=r.structure.atoms.filter(a=>a.symbol!=="H").map(a=>{const eligible=atomEligible(a);return `<button data-atom="${a.index}" class="${selected.includes(a.index)?"selected":""} ${eligible?"":"disabled"}" ${!eligible||!canEditSpectrum()?"disabled":""} title="${escape(a.symbol)}${a.index} · 连氢数 ${a.total_hydrogens} · 邻接 ${a.neighbors.join(",")}">${escape(a.symbol)}${a.index}</button>`;}).join("");const groups=r.structure.equivalent_groups.filter(g=>g.every(i=>state.nucleus==="13C"?r.structure.atoms[i].symbol==="C":r.structure.atoms[i].total_hydrogens>0));$("#groups").innerHTML=groups.length?"结构对称组提示（不自动证明化学等价）：<br>"+groups.map(g=>`<button data-group="${g.join(",")}" ${!canEditSpectrum()?"disabled":""}>${g.map(i=>escape(r.structure.atoms[i].symbol)+i).join(" / ")}</button>`).join(""):"该核种没有多原子的结构对称组提示。";renderEditor();}
let selectedSpectra={recordId:null,ids:{}};
function renderSpectra(){
  const r=state.record;if(!r)return;
  if(selectedSpectra.recordId!==r.id)selectedSpectra={recordId:r.id,ids:{}};
  $$('[data-nuc]').forEach(b=>{b.classList.toggle("active",b.dataset.nuc===state.nucleus);b.setAttribute("aria-selected",String(b.dataset.nuc===state.nucleus));});
  const choices=r.spectra.filter(s=>s.nucleus===state.nucleus).sort(spectrumOrder);
  if(!choices.some(s=>s.id===state.sid))state.sid=choices.find(s=>s.id===selectedSpectra.ids[state.nucleus])?.id||choices[0]?.id;
  if(state.sid)selectedSpectra.ids[state.nucleus]=state.sid;
  $("#spectrum").innerHTML=choices.map((s,i)=>`<option value="${s.id}">${i+1}. ${escape(s.source)} · ${escape(s.solvent||"溶剂未知")} · ${s.frequency_mhz??"频率未知"}${s.frequency_mhz?" MHz":""} · ${s.peaks.length} 信号 · ${escape(spectrumStateLabel(s))}</option>`).join("")||"<option>没有这类实验谱</option>";
  $("#spectrum").value=state.sid||"";
  const s=spectrum();
  if(!s){$("#peaks").innerHTML="";$("#sourceDetails").textContent="";$("#solventReference").textContent="";state.pid=null;window.renderNumberingNotice?.();window.renderSpectrumCuration?.();renderStructure();return;}
  if(window.NMRReference)window.NMRReference.render($("#solventReference"),{conditions:{by_nucleus:{[state.nucleus]:{solvent_raw:s.conditions.solvent_raw||s.solvent}}},spectra:{[state.nucleus]:{peaks:s.peaks}}},state.nucleus);
  renderPeaks();window.renderSpectrumCuration?.();
  const source=s.provenance;const doi=source.doi;
  $("#sourceDetails").innerHTML=(doi?`<a href="https://doi.org/${encodeURIComponent(doi).replace(/%2F/g,"/")}" target="_blank" rel="noopener">${escape(doi)}</a>`:"")+`<p>${escape(source.source_shard||source.source_file||"")} · ${source.source_row_index!==undefined?"原表行 "+source.source_row_index+"（从 0 开始）":"源记录 "+escape(source.source_record_index||"")}</p>`+(source.original_row?`<p>页码口径：筛选后 NMR 页面索引，不是原 SI 页码。原值 ${escape(source.original_row.Page_in_file_para)}</p><pre>${escape(source.original_row.NMR_shift_text||"")}</pre>`:`<p>ShiftDB2 ID：${escape(source.nmrshiftdb2_id||"")} · ${escape(source.spectrum_field||"")}</p><pre>${escape(s.peaks.map(p=>p.raw_signal||"").join(" | "))}</pre>`);
  renderStructure();
}
function spectrumForCopy(nucleus){
  const choices=(state.record?.spectra.filter(s=>s.nucleus===nucleus)||[]).sort(spectrumOrder);
  const remembered=selectedSpectra.recordId===state.record?.id?selectedSpectra.ids[nucleus]:null;
  return choices.find(s=>s.id===state.sid)||choices.find(s=>s.id===remembered)||choices[0];
}
function spectrumListText(s){
  if(!s?.peaks.length)return '未提供';
  return s.peaks.map(p=>{
    const low=p.range_min??p.rangeMin,high=p.range_max??p.rangeMax;
    const shift=low!=null&&high!=null&&Number(low)!==Number(high)?`${Math.max(Number(low),Number(high))}-${Math.min(Number(low),Number(high))}`:String(p.shift_ppm);
    if(s.nucleus==='13C')return shift;
    const details=[];
    if(p.multiplicity)details.push(p.multiplicity);
    if(p.j_values?.length)details.push('J = '+p.j_values.join(', ')+' Hz');
    if(p.reported_nH!=null)details.push(String(p.reported_nH).replace(/H$/i,'')+'H');
    return shift+(details.length?' ('+details.join(', ')+')':'');
  }).join(', ');
}
function spectraCopyText(){return '碳谱列表为：'+spectrumListText(spectrumForCopy('13C'))+'\n氢谱列表为：'+spectrumListText(spectrumForCopy('1H'));}
async function copySpectra(){
  if(!state.record)return toast('请先选择一个分子。');
  const text=spectraCopyText();
  try{
    let copied=false;
    if(typeof navigator!=='undefined'&&navigator.clipboard?.writeText){try{await navigator.clipboard.writeText(text);copied=true;}catch{}}
    if(!copied){
      const focused=document.activeElement,field=document.createElement('textarea');
      field.value=text;field.readOnly=true;field.style.cssText='position:fixed;left:-9999px;top:0';document.body.appendChild(field);
      try{field.select();copied=document.execCommand('copy');}finally{field.remove();focused?.focus?.();}
    }
    if(!copied)throw new Error('未取得剪贴板权限，请检查浏览器的复制权限。');
    toast('已复制 C/H 谱列表。');
  }catch(e){toast('复制失败：'+e.message);}
}
function renderPeaks(){
  const s=spectrum();if(!s)return;
  window.renderNumberingNotice?.();
  $('#assignmentWarning')?.remove?.();
  const validation=s.source_assignment_validation;
  if(validation&&validation.status!=='identity_validated'){
    const reason=validation.status==='needs_attention'?'原库归属存在氢数编码等警告，需核实':validation.status==='record_assignment_conflict'?'同一源记录的实验 C 谱原子身份冲突，本 H 谱源归属暂不采信':'原库归属与原始结构编号不一致';
    $('#spectrumTools').insertAdjacentHTML('afterend',`<p id="assignmentWarning" role="status" style="color:#9b401b;background:#fff3df;padding:12px;border-radius:8px">${escape(reason)}。峰位与原文保留；不按位移猜测归属。人工标注保持不变。</p>`);
  }
  $("#peaks").innerHTML=s.peaks.map(p=>{
    const a=p.assignment;
    const annotationOrigin=String(a.origin||'').startsWith('ai:')?'AI · '+(a.ai?.response_model||a.ai?.model||String(a.origin).slice(12)):a.origin==='nmrshiftdb2'?'原库归属':a.origin?'人工标注':'';
    const atoms=a.atom_indices.map(i=>state.record.structure.atoms[i]?.symbol+i).join(", ");
    const detail=[p.multiplicity,p.reported_nH!=null?p.reported_nH+"H":"",p.j_values?.length?"J "+p.j_values.join(", ")+" Hz":""].filter(Boolean).join(" · ");
    const fragments=p.source_assignment_fragments||[];
    const references=fragments.flatMap(f=>f.assignments.map(x=>`${x.source_assignment_id} → ${x.display_atom_label||'无法映射'}${x.transfer_status==='mapped'?'':'（冲突）'}`));
    const warning=validation&&validation.status!=='identity_validated'?`<small title="${escape(references.join('; '))}" style="color:#9b401b">源归属待核实${a.origin==='manual'?' · 当前为人工标注':''}</small>`:'';
    return `<tr data-peak="${p.peak_id}" tabindex="0" class="${p.peak_id===state.pid?"active":""}"><td>${escape(p.peak_id)}</td><td class="shift">${p.shift_ppm}<small>${p.range_min!==undefined&&p.range_min!==p.range_max?`${p.range_min}–${p.range_max}`:""}</small></td><td class="detail">${escape(detail||"未报告")}<small>${p.source_signal_code?"原库编码 "+escape(p.source_signal_code):""}</small></td><td class="target">${escape(atoms||"—")}<small title="${escape(a.ai?'任务 '+a.ai.job_id+' · '+a.ai.created_at:'')}">${escape(annotationOrigin)}</small>${warning}${p.source_assignment_note?`<small title="${escape(p.source_assignment_note)}">原库指向 ${escape((p.source_assignment_candidates||[]).map(i=>state.record.structure.atoms[i].symbol+i).join(" / "))} · 氢数需核对</small>`:""}</td><td class="${a.status}">${escape(labels[a.status]||a.status)}</td></tr>`;
  }).join("");
}
function selectPeak(pid){$("#atomNumber").value="";state.pid=pid;renderPeaks();renderStructure();}
function renderEditor(){const p=peak();$("#editor").hidden=!p;if(!p)return;$("#peakTitle").textContent=`${state.nucleus} · ${p.peak_id} · ${p.shift_ppm} ppm`;$("#selectedAtoms").textContent=p.assignment.atom_indices.length?p.assignment.atom_indices.map(i=>state.record.structure.atoms[i].symbol+i).join("、"):"尚未选择原子";$("#status").value=p.assignment.status;$("#peakNote").value=p.assignment.note||"";$("#status").disabled=!canEditSpectrum();$("#peakNote").readOnly=!canEditSpectrum();$("#clearAtoms").disabled=!canEditSpectrum();$("#atomNumber").disabled=!canEditSpectrum();$("#atomApply").disabled=!canEditSpectrum();}
function changeAtoms(indices,group=false){
  const p=peak();
  if(!p)return toast("请先选一个谱峰。");
  if(p.assignment_group)return toast('此峰属于位移组，组内顺序未定；请在管理当前谱中解除分组后再编辑。');
  if(!canEditSpectrum())return toast("当前谱为只读存档，请恢复后再修改归属。");
  for(const i of indices){
    const a=state.record.structure.atoms[i];
    if(!a)return toast("没有原子编号 "+i+"，请按结构图中的编号输入。");
    if(!atomEligible(a))return toast("这个原子不适合作为当前核种的归属目标。");
  }
  let set=new Set(p.assignment.atom_indices);
  // Extending an explicitly selected equivalent group no longer describes
  // just that group. Pure removals keep the remaining group interpretation.
  const extendsEquivalentGroup=!group&&p.assignment.status==="equivalent_group"&&indices.some(i=>!set.has(i));
  const deselectGroup=group&&p.assignment.status==='equivalent_group'&&set.size===indices.length&&indices.every(i=>set.has(i));
  if(group)set=deselectGroup?new Set():new Set(indices);
  else indices.forEach(i=>set.has(i)?set.delete(i):set.add(i));
  p.assignment.atom_indices=[...set].sort((a,b)=>a-b);
  if(group)p.assignment.status=set.size?"equivalent_group":"unassigned";
  else if(!set.size)p.assignment.status="unassigned";
  else if(set.size===1)p.assignment.status="assigned";
  else if(extendsEquivalentGroup||!["equivalent_group","overlap","ambiguous"].includes(p.assignment.status))p.assignment.status="overlap";
  p.assignment.origin="manual";p._dirty=true;dirty(p);renderPeaks();renderStructure();return true;
}
function stepPeak(direction){
  const s=spectrum();if(!s?.peaks.length)return;
  const index=s.peaks.findIndex(p=>p.peak_id===state.pid);
  if(direction>0&&state.nucleus==='13C'&&index===s.peaks.length-1){
    const choices=state.record.spectra.filter(h=>h.nucleus==='1H'&&h.peaks.length&&canEditSpectrum(h)).sort(spectrumOrder);
    const remembered=selectedSpectra.recordId===state.record.id?selectedSpectra.ids['1H']:null;
    const next=choices.find(h=>h.id===remembered)||choices[0];
    if(!next)return toast('当前分子没有可用的 H 谱峰，仍停留在 C 谱最后一个峰。');
    state.nucleus='1H';state.sid=next.id;state.pid=next.peaks[0].peak_id;
    $('#atomNumber').value='';renderSpectra();
    $('#peaks').querySelector('tr[data-peak]')?.scrollIntoView?.({block:'center',inline:'nearest'});
    return;
  }
  selectPeak(s.peaks[Math.min(s.peaks.length-1,Math.max(0,index+direction))].peak_id);
}
function clearCurrentAssignment(){
  const p=peak();if(!p)return toast('请先选择需要清空归属的峰。');
  if(!canEditSpectrum())return toast('当前谱为只读存档，不能清空归属。');
  if(p.assignment_group)return toast('此峰属于组内未定的位移组，请在管理当前谱中解除分组。');
  if(p.assignment.status==='unassigned'&&!p.assignment.atom_indices.length)return;
  $('#atomNumber').value='';
  p.assignment={status:'unassigned',atom_indices:[],note:p.assignment.note||'',origin:'manual'};
  p._dirty=true;dirty(p);renderPeaks();renderStructure();
}
$("#collapse").onclick=()=>{$(".shell").classList.toggle("collapsed");localStorage.setItem("nmr-db-sidebar",$(".shell").classList.contains("collapsed")?"closed":"open");};if(localStorage.getItem("nmr-db-sidebar")==="closed")$(".shell").classList.add("collapsed");
$("#list").onclick=e=>{const b=e.target.closest("[data-mid]");if(b&&!state.navigating)openRecord(b.dataset.mid).catch(x=>toast(x.message));};
let searchTimer;$("#search").oninput=()=>{clearTimeout(searchTimer);searchTimer=setTimeout(()=>loadList(true).catch(e=>toast(e.message)),300);};
$("#search").onkeydown=e=>{if(e.key==="Enter"){e.preventDefault();clearTimeout(searchTimer);loadList(true).catch(x=>toast(x.message));}};
["#queue","#origin","#coverage","#searchScope","#annotationSource"].forEach(id=>{if($(id))$(id).onchange=()=>loadList(true).catch(e=>toast(e.message));});
$('#processing').onchange=()=>{if($('#view').value==='finished'&&$('#processing').value!=='finished')$('#view').value='';loadList(true).catch(e=>toast(e.message));};
$("#view").onchange=()=>{clearTimeout(searchTimer);$("#search").value="";$("#queue").value="";$("#origin").value="";$('#coverage').value='';$('#processing').value=$('#view').value==='finished'?'finished':'all';loadList(true).catch(e=>toast(e.message));};
$("#refreshList").onclick=()=>loadList(true).catch(e=>toast(e.message));
$("#prevPage").onclick=()=>{state.offset=Math.max(0,state.offset-40);loadList().catch(e=>toast(e.message));};$("#nextPage").onclick=()=>{state.offset+=40;loadList().catch(e=>toast(e.message));};
$$('[data-nuc]').forEach(b=>b.onclick=()=>{state.nucleus=b.dataset.nuc;state.pid=null;$("#atomNumber").value="";renderSpectra();});$("#spectrum").onchange=()=>{state.sid=$("#spectrum").value;state.pid=null;$("#atomNumber").value="";renderSpectra();};
$("#peaks").onclick=e=>{const row=e.target.closest("[data-peak]");if(row)selectPeak(row.dataset.peak);};$("#peaks").onkeydown=e=>{if(e.key==="Enter"){e.preventDefault();if($("#atomNumber").value){e.stopPropagation();commitAtomNumber();return;}const row=e.target.closest("[data-peak]");if(row)selectPeak(row.dataset.peak);}};
$("#atoms").onclick=e=>{const b=e.target.closest("[data-atom]");if(b)changeAtoms([Number(b.dataset.atom)]);};$("#groups").onclick=e=>{const b=e.target.closest("[data-group]");if(b)changeAtoms(b.dataset.group.split(",").map(Number),true);};
$("#status").onchange=()=>{const p=peak();if(!p||!canEditSpectrum()||p.assignment_group)return;p.assignment.status=$("#status").value;if(["unassigned","solvent","impurity","excluded"].includes(p.assignment.status))p.assignment.atom_indices=[];p.assignment.origin="manual";p._dirty=true;dirty(p);renderPeaks();renderStructure();};
$("#peakNote").oninput=()=>{const p=peak();if(!p||!canEditSpectrum()||p.assignment_group)return;p.assignment.note=$("#peakNote").value;p._dirty=true;dirty(p);};$("#recordNote").oninput=()=>{if(state.record&&canEditMolecule()){state.record.note=$("#recordNote").value;dirty();}};
$("#closeEditor").onclick=()=>{state.pid=null;renderPeaks();renderStructure();};$("#clearAtoms").onclick=clearCurrentAssignment;
$("#prevPeak").onclick=()=>stepPeak(-1);$("#nextPeak").onclick=()=>stepPeak(1);$("#save").onclick=()=>save().then(()=>toast("已保存")).catch(e=>toast(e.message));$("#complete").onclick=()=>save(true).catch(e=>toast(e.message));
$("#nextTask").onclick=()=>goNext().catch(e=>toast(e.message));
$("#prevTask").onclick=()=>goPrevious().catch(e=>toast(e.message));
$("#completeNext").onclick=()=>goNext(true).catch(e=>toast(e.message));
function handleShortcut(e){
  const command=e.ctrlKey||e.metaKey,key=e.key.toLowerCase();
  const saveKey=command&&key==="s",finishKey=command&&e.key==="Enter",nextKey=e.altKey&&e.key==="ArrowRight",prevKey=e.altKey&&e.key==="ArrowLeft";
  if(saveKey||finishKey||nextKey||prevKey)e.preventDefault();
  if(e.repeat||e.isComposing||$$('dialog[open]').length)return;
  if(saveKey){if(!canEditMolecule())return toast('当前分子为只读存档，可查看谱峰；恢复后再标注。');save().catch(x=>toast(x.message));return;}
  if(prevKey){goPrevious().catch(x=>toast(x.message));return;}
  if(finishKey||nextKey){goNext(finishKey).catch(x=>toast(x.message));return;}
  if(command||e.altKey||["INPUT","TEXTAREA","SELECT"].includes(document.activeElement?.tagName)||document.activeElement?.isContentEditable)return;
  if(key==='e'){e.preventDefault();clearCurrentAssignment();return;}
  if(/^[0-9]$/.test(e.key)){
    if(!peak())return toast('先用 J / K 选择一个峰，再输入原子编号。');
    if(!canEditSpectrum())return toast('当前谱为只读存档，可浏览；恢复后才能改归属。');
    e.preventDefault();$('#atomNumber').value+=e.key;return;
  }
  if(e.key==='Enter'&&$('#atomNumber').value){e.preventDefault();commitAtomNumber();return;}
  if(e.key==='Backspace'&&$('#atomNumber').value){e.preventDefault();$('#atomNumber').value=$('#atomNumber').value.slice(0,-1);return;}
  if(e.key==='Escape'&&$('#atomNumber').value){e.preventDefault();$('#atomNumber').value='';return;}
  if((key==='c'||key==='h')&&state.record){state.nucleus=key==='c'?'13C':'1H';state.pid=null;$('#atomNumber').value='';renderSpectra();return;}
  if(key==="j")stepPeak(-1);if(key==="k")stepPeak(1);if(e.key==="Escape"){state.pid=null;renderPeaks();renderStructure();}
}

function commitAtomNumber(){
  const text=$('#atomNumber').value.trim();
  if(!/^\d+(?:[\s,，]+\d+)*$/.test(text))return toast('请输入原子编号，例如 12；也可输入 12,13。');
  const indices=[...new Set(text.split(/[\s,，]+/).map(Number))];
  if(changeAtoms(indices))$('#atomNumber').value='';
}
$('#atomNumber').onkeydown=e=>{if(e.key==='Enter'){e.preventDefault();e.stopPropagation();commitAtomNumber();}else if(e.key==='Escape'){e.preventDefault();e.stopPropagation();e.target.value='';e.target.blur();}};
$('#atomApply').onclick=commitAtomNumber;
$('#copySpectra').onclick=copySpectra;
$('#structureView').onchange=renderStructure;
$('#structureEnlarge').onclick=()=>{$('#moleculeLarge').src=$('#molecule').src;$('#structureDialog').showModal();};
$('#structureClose').onclick=()=>$('#structureDialog').close();
document.addEventListener("keydown",handleShortcut);
window.addEventListener("beforeunload",e=>{if(state.dirty||state.saving){e.preventDefault();e.returnValue="";}});
$("#exportOpen").onclick=()=>$("#exportDialog").showModal();$("#exportRun").onclick=async()=>{const button=$("#exportRun");button.disabled=true;try{await save();const scope=$("#exportScope").value;if(scope==="current"&&!state.record)throw new Error("请先选择一个分子。");const payload={format:$("#exportFormat").value,molecule_id:scope==="current"?state.record.id:null,filters:scope==="filtered"?filters():{}};const job=await api("/api/export",payload);$("#exportResult").textContent="正在生成一致性快照，可继续标注…";const poll=async()=>{try{const result=await api("/api/export/"+job.job_id);if(result.state==="running"){setTimeout(poll,1800);return;}if(result.state==="error")throw new Error(result.error);$("#exportResult").innerHTML=`已导出 ${result.count.toLocaleString()} 个分子。<br><a href="/api/export/${job.job_id}/download">下载文件</a><br><small>${escape(result.file)}</small>`;button.disabled=false;}catch(e){$("#exportResult").textContent=e.message;button.disabled=false;}};poll();}catch(e){$("#exportResult").textContent=e.message;button.disabled=false;}};
Promise.all([loadList(),loadStats()]).then(()=>{const mid=new URLSearchParams(location.search).get('molecule_id');if(mid&&/^M[a-f0-9]{24}$/.test(mid))return openRecord(mid);}).catch(e=>toast(e.message));setInterval(()=>loadStats().catch(()=>{}),60000);
