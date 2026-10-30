"use strict";
// Separate panel: original annotation save serialization remains unchanged.
labels.manual_hold="人工暂时停用";
labels.excluded_training="排除训练";
const managementLabels={active:"当前数据集（导出）",held:"暂时停用",excluded:"排除训练",reference:"参考版本（不导出）",superseded:"历史参考（不导出）",unknown:"完整性未核实",incomplete:"峰表不完整",complete:"已核实完整",not_checked:"未做自动检查",passed:"通过自动检查",quarantined:"自动隔离"};
const managementActions={disposition:"更改使用状态",issue:"记录问题",resolve_issue:"标记问题已处理",completeness:"记录谱完整性",clone:"新建归属修订副本",import_bmrb:"引入 BMRB 参考",adopt:"采用版本"};
const sourceNames={nmrshiftdb2:'ShiftDB2 原库',NMRexp:'NMRexp',manual_revision:'人工修订谱',BMRB_reference:'BMRB 实验参考'};
let management=null,managementBusy=false;
$("#queue").insertAdjacentHTML("beforeend",'<option value="quality_hold">自动隔离后无可用谱</option><option value="manual_hold">人工暂时停用</option><option value="excluded_training">排除训练</option>');
$('#spectrumTools').insertAdjacentHTML('afterbegin','<button id="manageOpen">管理当前谱</button><span id="curationSummary" hidden></span><button id="restoreMolecule" hidden>恢复分子</button>');
$('#spectrumTools').insertAdjacentHTML('afterend','<section id="spectrumUseState" class="spectrum-use-state" hidden><span id="spectrumUseText"></span><button id="restoreSpectrum" hidden>恢复当前谱</button></section>');
document.body.insertAdjacentHTML("beforeend",`
<dialog id="managementDialog" aria-labelledby="managementTitle">
 <div class="management-top"><div class="management-head"><h2 id="managementTitle">谱图管理</h2><button id="manageClose" aria-label="关闭谱图管理">×</button></div>
 <div class="management-context"><p id="manageRecord"></p><label for="manageTarget">操作对象<select id="manageTarget"></select></label></div></div>
 <nav class="management-nav" role="tablist" aria-label="谱图管理功能">
  <button id="manageTabPeaks" role="tab" data-management-tab="peaks" aria-controls="managePeaksPanel" aria-selected="true">实验峰表</button>
  <button id="manageTabRevision" role="tab" data-management-tab="revision" aria-controls="manageRevisionPanel" aria-selected="false">编码与分组</button>
  <button id="manageTabStatus" role="tab" data-management-tab="status" aria-controls="manageStatusPanel" aria-selected="false">状态与完整性</button>
  <button id="manageTabReference" role="tab" data-management-tab="reference" aria-controls="manageReferencePanel" aria-selected="false">参考与版本</button>
  <button id="manageTabSource" role="tab" data-management-tab="source" aria-controls="manageSourcePanel" aria-selected="false">来源与历史</button>
 </nav>
 <div id="managementBody" class="management-body">
 <section id="managePeaksPanel" role="tabpanel" aria-labelledby="manageTabPeaks" data-management-panel="peaks"></section>
 <section id="manageRevisionPanel" role="tabpanel" aria-labelledby="manageTabRevision" data-management-panel="revision" hidden></section>
 <section id="manageStatusPanel" role="tabpanel" aria-labelledby="manageTabStatus" data-management-panel="status" hidden>
 <section class="management-section"><h3>使用状态</h3><div class="manage-actions"><button data-disposition="held">暂时停用</button><button data-disposition="excluded">排除训练（保留存档）</button><button data-disposition="active">恢复使用（原因可选）</button><button id="manageClone">打开 / 建立最终版</button></div><p>每条实验谱保留原始版和一份最终版。首次修改建立最终版，后续修改持续更新同一份最终版；原始版可只读对照。</p></section>
 <section class="management-section"><h3>记录问题与完整性</h3><div class="manage-fields"><label>问题类型<select id="manageCategory"><option value="missing_peaks">实验峰表缺项</option><option value="wrong_assignment">原子归属可疑</option><option value="structure">结构 / 立体信息问题</option><option value="conditions">实验条件问题</option><option value="source">来源或匹配问题</option><option value="other">其他</option></select></label><button id="manageIssue">记录问题</button><label>当前谱完整性<select id="manageCompleteness"><option value="unknown">尚未核实</option><option value="incomplete">已确认不完整</option><option value="complete">已核实完整</option></select></label><button id="manageSetCompleteness">保存完整性</button></div><p>记录问题、标记不完整不会自动停用；如不宜进入训练，请另外选择“排除训练”。峰数少于原子数也可能是等价或重叠，不能单凭数量判为缺峰。</p><div id="manageIssues"></div></section>
 </section>
 <section id="manageReferencePanel" role="tabpanel" aria-labelledby="manageTabReference" data-management-panel="reference" hidden>
 <section class="management-section"><h3>引入 BMRB 实验参考</h3><div class="manage-fields"><label>官方条目编号<input id="manageEntry" placeholder="例如 bmse000194" autocomplete="off"></label><button id="manageImport">下载并存为参考（不替换）</button></div><p>仅导入可确定结构与样品关联的实验原子位移表。保存原始 STAR / MOL、来源和原子映射。不导入预测值，不推算缺失峰，不自动确定当前分子的立体结构。</p></section>
 <section class="management-section"><h3>采用所选参考 / 恢复历史版本</h3><label>替换哪一条同核种谱<select id="manageReplace"></select></label><div id="manageCompare"></div><label class="inline"><input id="manageAcknowledge" type="checkbox">我已比较结构、峰表、归属和条件；理解采用不代表补齐缺峰，也不证明 H/C 为同一样品。</label><button id="manageAdopt">采用此版本（旧谱保留）</button><p>只替换明确选中的谱；另一核种不动，不自动建立新的同条件双谱。恢复旧版本也用此操作，并指定需要撤下的新版本。</p></section>
 </section>
 <section id="manageSourcePanel" role="tabpanel" aria-labelledby="manageTabSource" data-management-panel="source" hidden>
 <p class="notice">原始版用于只读对照，最终版用于标注与导出。后续修改直接保存到最终版；操作说明保留在下方记录中。</p>
 <details class="management-section"><summary>查看当前峰表、来源与条件</summary><div id="manageInspection"></div></details>
 <details class="management-section"><summary>操作记录（最近 50 条）</summary><div id="manageHistory"></div></details>
 </section>
 </div>
 <div class="management-footer">
 <label for="manageReason">本次操作原因 / 原文位置 <span>恢复使用可留空</span><textarea id="manageReason" rows="2" placeholder="填写实验记录、文献页码或修改原因，再点击对应操作的保存按钮。" aria-describedby="manageReasonHint"></textarea></label>
 <div class="management-feedback"><small id="manageReasonHint">各项操作单独保存，共用此处说明。</small><p id="manageMessage" role="status" aria-live="polite"></p></div>
 </div>
</dialog>`);
let managementTab='peaks';
function setManagementTab(tab){
 const tabs=$$('[data-management-tab]');
 if(!tabs.some(b=>b.dataset.managementTab===tab))return;
 managementTab=tab;
 for(const button of tabs){const selected=button.dataset.managementTab===tab;button.setAttribute('aria-selected',String(selected));button.tabIndex=selected?0:-1;}
 for(const panel of $$('[data-management-panel]'))panel.hidden=panel.dataset.managementPanel!==tab;
 $('#managementBody').scrollTop=0;
}

window.renderCurationSummary=r=>{
 const status=r.curation?.molecule_state||"active";
 const held=status!=="active";
 $(".workgrid").inert=false;
 $('#recordNote').readOnly=held;
 $('#restoreMolecule').hidden=!held;
 $("#save").disabled=held;$("#complete").disabled=held;$("#completeNext").disabled=held||!!state.navigating;
 window.renderSpectrumCuration();
};
window.renderSpectrumCuration=()=>{
 const s=spectrum();
 const summary=[];
 if(!canEditMolecule())summary.push('分子'+(managementLabels[state.record?.curation?.molecule_state]||'已停用'));
 if(s&&s.manual_state&&s.manual_state!=='active')summary.push(managementLabels[s.manual_state]||s.manual_state);
 if(s?.completeness?.status&&s.completeness.status!=='unknown')summary.push(managementLabels[s.completeness.status]);
 $('#curationSummary').textContent=summary.join(' · ');$('#curationSummary').hidden=!summary.length;
 $('#manageOpen').textContent=s?'管理当前谱':'数据管理';
 $('#spectrumUseState').hidden=!s||canEditSpectrum(s);
 $('#spectrumUseText').textContent=s?(spectrumStateLabel(s)+' · 只读'):'';
 $('#restoreSpectrum').hidden=!s||!canEditMolecule()||!['held','excluded'].includes(s.manual_state)||s.quality?.decision==='quarantined';
};
function managedSpectrum(){return management?.spectra.find(s=>s.id===$("#manageTarget").value);}
function atomText(a){return a?.atom_indices?.length?a.atom_indices.map(i=>(state.record.structure.atoms[i]?.symbol||"?")+i).join(" / "):"未归属";}
function manageInspect(){
 const s=managedSpectrum();
 $("#manageAcknowledge").checked=false;
 const versions=management.spectra.filter(x=>s&&x.nucleus===s.nucleus&&x.id!==s.id&&!["reference","superseded"].includes(x.manual_state));
 $("#manageReplace").innerHTML='<option value="">新增为独立实验记录，不替换</option>'+versions.map(x=>`<option value="${escape(x.id)}">${escape(x.source)} · ${x.peaks.length} 项 · ${escape(x.id)} · ${escape(managementLabels[x.manual_state])}</option>`).join("");
 // Require an explicit replacement target; never silently pick a spectrum.
 if(!s){$("#manageInspection").textContent="对整个分子操作：影响该分子所有当前谱的标注与导出。";$("#manageCompare").textContent="先选择参考谱或历史版本。";}
 else{
  const annotations=new Map(s.annotations.map(a=>[a.peak_id,a]));
  $("#manageCompleteness").value=s.completeness.status;
  $("#manageInspection").innerHTML=`<p><strong>${escape(managementLabels[s.manual_state])}</strong> · ${escape(managementLabels[s.automatic_quality.decision]||s.automatic_quality.decision)} · ${escape(managementLabels[s.completeness.status])}<br>${escape(s.source)} · ${escape(s.nucleus)} · ${s.peaks.length} 个记录信号 · 溶剂 ${escape(s.solvent||"未知")} · 温度 ${escape(s.temperature_k??"未知")} K</p><div class="management-warnings">${(s.provenance.warnings||[]).map(escape).join("<br>")||"这是原始 / 历史归属的只读视图。已有标签不保证正确，谱表完整性需另行判断。"}${s.automatic_quality.reasons?.length?"<br>自动隔离原因："+s.automatic_quality.reasons.map(escape).join("；"):""}</div><div class="management-table"><table><thead><tr><th>δ / ppm</th><th>当前结构原子</th><th>状态</th><th>来源原子</th></tr></thead><tbody>${s.peaks.map(p=>`<tr><td>${escape(p.shift_ppm)}</td><td>${escape(atomText(annotations.get(p.peak_id)))}</td><td>${escape(labels[annotations.get(p.peak_id)?.status]||"未归属")}</td><td>${escape((p.source_atom_ids||[]).join(" / "))}</td></tr>`).join("")}</tbody></table></div><details><summary>来源、原子映射及实验条件</summary><pre>${escape(JSON.stringify({conditions:s.conditions,provenance:s.provenance,automatic_quality:s.automatic_quality,completeness:s.completeness},null,2))}</pre></details>`;
  manageCompare();
 }
 $("#manageClone").disabled=!s;$("#manageSetCompleteness").disabled=!s;
 if(s)$('#manageInspection').insertAdjacentHTML('afterbegin',`<p class="management-warnings"><strong>${s.dataset_role==='current_dataset'?'当前数据集版本：保存标注与导出使用此谱':'参考 / 停用版本：不进入数据集导出'}</strong>${s.version_info?'<br>最终版 · 原始谱 '+escape(s.version_info.root_spectrum_id):''}</p>`);
 $("#manageAdopt").disabled=!s||!["reference","superseded"].includes(s.manual_state);
 $$('[data-disposition]').forEach(b=>b.disabled=!!s&&["reference","superseded"].includes(s.manual_state));
}
function manageCompare(){
 const s=managedSpectrum(),old=management.spectra.find(x=>x.id===$("#manageReplace").value);
 if(!s)return;
 if(!old){$("#manageCompare").textContent="未选择替换对象：采用后将新增一条独立实验记录。若只是参考，请勿采用。";return;}
 const ref=new Map(s.annotations.map(a=>[a.peak_id,a])),prev=new Map(old.annotations.map(a=>[a.peak_id,a]));
 $("#manageCompare").innerHTML=`<p>原版本 ${old.peaks.length} 项 → 所选版本 ${s.peaks.length} 项。当前结构含 ${state.record.structure.atoms.filter(a=>a.symbol==='C').length} 个 C 原子。仅对照完全相同的位移；这不是自动配峰。</p><p>溶剂：${escape(old.solvent||'未知')} → ${escape(s.solvent||'未知')}；温度：${escape(old.temperature_k??'未知')} → ${escape(s.temperature_k??'未知')} K。另一核种的条件不会跟着改变。</p><div class="management-table"><table><thead><tr><th>所选版本 δ</th><th>原版本归属</th><th>所选版本归属</th></tr></thead><tbody>${s.peaks.map(p=>{const matches=old.peaks.filter(x=>x.shift_ppm===p.shift_ppm);const before=matches.length===1?atomText(prev.get(matches[0].peak_id)):matches.length?'多个同位移记录，不能唯一对照':'原版本无同位移项';const after=atomText(ref.get(p.peak_id));return `<tr class="${before===after?'':'manage-difference'}"><td>${escape(p.shift_ppm)}</td><td>${escape(before)}</td><td>${escape(after)}</td></tr>`;}).join('')}</tbody></table></div>`;
}
function renderManagement(selected=''){
 const r=state.record;
 $("#manageRecord").textContent=`${r.formula} · ${r.id}`;
 $("#manageTarget").innerHTML='<option value="">整个分子（影响全部谱）</option>'+management.spectra.map(s=>`<option value="${escape(s.id)}">${escape(s.nucleus)} · ${escape(s.version_label||s.source)} · ${s.peaks.length} 项 · ${escape(managementLabels[s.manual_state])} · ${escape(s.id)}</option>`).join('');
 $("#manageTarget").value=selected;
 $("#manageIssues").innerHTML=management.issues.map(i=>`<div class="manage-issue">${i.status==='open'?'待处理':'已处理'} · ${escape(i.category)} · ${escape(i.details)} ${i.spectrum_id?'<small>'+escape(i.spectrum_id)+'</small>':''}${i.status==='open'?`<button data-resolve="${i.id}">标记已处理</button>`:''}</div>`).join('')||'<p>尚未记录问题。</p>';
 $("#manageHistory").innerHTML=management.events.map(e=>`<p>${escape(e.created_at)} · ${escape(managementActions[e.action]||e.action)}<br>${escape(JSON.stringify(e.details))}</p>`).join('')||'暂无管理操作。';
 manageInspect();
}
async function managementAction(payload,quick=false){
 if(managementBusy)return;
 const reason=quick?'':$("#manageReason").value.trim();
 if(!reason&&!(payload.action==='disposition'&&payload.state==='active')){
  $('#manageReason').setAttribute?.('aria-invalid','true');$('#manageReason').focus?.();
  return toast("请在窗口底部填写处理原因或原文位置，再提交当前操作。");
 }
 const mid=state.record.id;
 managementBusy=true;
 $("#managementDialog").querySelectorAll('button,input,select,textarea').forEach(x=>x.disabled=true);
 $("#manageMessage").textContent=payload.action==='import_bmrb'?"正在下载官方参考并核对结构；不会替换现有谱…":"正在保存处理决定…";
 try{
  const result=await api('/api/molecule/'+mid+'/management',{...payload,revision:management.revision,reason,spectrum_id:Object.hasOwn(payload,'spectrum_id')?payload.spectrum_id:$("#manageTarget").value||null});
  state.record=result.record;state.dirty=false;state.pid=null;management=result.management;
  if(result.event.current_spectrum_id)state.sid=result.event.current_spectrum_id;
  else if(['supplement_peak','withdraw_peak_supplement'].includes(payload.action))state.sid=result.event.spectrum_id;
  if(result.event.added_peak_id)state.pid=result.event.added_peak_id;
  if(!state.record.spectra.some(s=>s.id===state.sid))state.sid=state.record.spectra[0]?.id;
  state.nucleus=spectrum()?.nucleus||'13C';
  $("#saveState").textContent='已保存';
  renderHeader();renderSpectra();
  const selected=result.event.current_spectrum_id||result.event.new_spectrum_ids?.[0]||payload.spectrum_id||$("#manageTarget").value;
  $("#managementDialog").querySelectorAll('button,input,select,textarea').forEach(x=>x.disabled=false);
  renderManagement(selected);
  $("#manageMessage").textContent=result.event.current_spectrum_id?'最终版已保存，可关闭本窗口标注。后续修改持续更新此最终版；原始版用于对照。':result.event.reused_existing_final?'已打开已有最终版。':payload.action==='import_bmrb'?'已存为参考。请查看来源和对照，再决定是否采用；原谱未修改。':'处理已保存。';
  if(payload.action==='disposition'&&payload.state==='active'){
   const remaining=state.record.spectra.filter(s=>['held','excluded'].includes(s.manual_state));
   toast(payload.spectrum_id?'当前谱已恢复；旧原因保留在操作记录。':remaining.length?'分子已恢复；仍有 '+remaining.length+' 条谱单独停用，可在该谱上点击“恢复当前谱”。':'分子已恢复，可继续标注；旧原因保留在操作记录。');
  }
  refreshProgressSoon();
  Promise.all([loadList(),loadStats()]).catch(e=>toast('处理已保存；列表刷新失败：'+e.message));
  return result;
 }catch(e){$("#manageMessage").textContent='操作未确认完成：'+e.message+' 如连接中断，请关闭后重新打开管理面板确认记录，不要重复提交。';if(quick)toast('恢复结果未确认：'+e.message);}
 finally{managementBusy=false;$("#managementDialog").querySelectorAll('button,input,select,textarea').forEach(x=>x.disabled=false);manageInspect();}
}
$("#manageOpen").onclick=async()=>{try{if(!state.record)return;await save();management=await api('/api/molecule/'+state.record.id+'/management');if(management.revision!==state.record.revision){await openRecord(state.record.id);management=await api('/api/molecule/'+state.record.id+'/management');}renderManagement(state.sid||'');setManagementTab(state.sid?'peaks':'status');$("#manageReason").value='';$('#manageReason').removeAttribute?.('aria-invalid');$("#manageMessage").textContent='修改后点击对应操作按钮保存，原因可在底部随时填写。';$("#managementDialog").showModal();}catch(e){toast(e.message);}};
$("#manageClose").onclick=()=>{if(!managementBusy)$("#managementDialog").close();};
async function restoreTarget(spectrum_id=null){
 if(managementBusy||state.loading||state.navigating||!state.record)return;
 const mid=state.record.id;state.loading=true;
 try{await save();management=await api('/api/molecule/'+mid+'/management');await managementAction({action:'disposition',state:'active',spectrum_id},true);}catch(e){toast(e.message);}finally{state.loading=false;}
}
$('#restoreMolecule').onclick=()=>restoreTarget();
$('#restoreSpectrum').onclick=()=>restoreTarget(state.sid);
$("#managementDialog").addEventListener('cancel',e=>{if(managementBusy)e.preventDefault();});
$('.management-nav').onclick=e=>{const button=e.target.closest('[data-management-tab]');if(button)setManagementTab(button.dataset.managementTab);};
$('.management-nav').onkeydown=e=>{
 if(!['ArrowLeft','ArrowRight','Home','End'].includes(e.key))return;
 const tabs=$$('[data-management-tab]');if(!tabs.length)return;
 const index=tabs.findIndex(b=>b.dataset.managementTab===managementTab);
 const next=e.key==='Home'?0:e.key==='End'?tabs.length-1:(index+(e.key==='ArrowRight'?1:-1)+tabs.length)%tabs.length;
 e.preventDefault();setManagementTab(tabs[next].dataset.managementTab);tabs[next].focus();
};
$('#manageReason').oninput=()=>$('#manageReason').removeAttribute?.('aria-invalid');
$("#manageTarget").onchange=()=>{manageInspect();if(!managedSpectrum())setManagementTab('status');};$("#manageReplace").onchange=()=>{ $("#manageAcknowledge").checked=false;manageCompare();};
$$('[data-disposition]').forEach(b=>b.onclick=()=>managementAction({action:'disposition',state:b.dataset.disposition}));
$("#manageClone").onclick=()=>managementAction({action:'clone'});
$("#manageIssue").onclick=()=>managementAction({action:'issue',category:$("#manageCategory").value});
$("#manageSetCompleteness").onclick=()=>managementAction({action:'completeness',status:$("#manageCompleteness").value});
$("#manageImport").onclick=()=>managementAction({action:'import_bmrb',entry_id:$("#manageEntry").value.trim()});
$("#manageAdopt").onclick=()=>managementAction({action:'adopt',replace_spectrum_id:$("#manageReplace").value||null,acknowledge:$("#manageAcknowledge").checked});
$("#manageIssues").onclick=e=>{const b=e.target.closest('[data-resolve]');if(b)managementAction({action:'resolve_issue',issue_id:Number(b.dataset.resolve)});};

// Extend source display without changing the annotation save state machine.
const originalRenderPeaks=renderPeaks;
renderPeaks=function(){
 originalRenderPeaks();
 const s=spectrum();if(!s)return;
 for(const row of $$('#peaks tr[data-peak]')){
  const p=s.peaks.find(x=>x.peak_id===row.dataset.peak);
  const label=row.querySelector('.target small');
  if(label&&p)label.textContent=String(p.assignment.origin||'').startsWith('ai:')?'AI · '+(p.assignment.ai?.response_model||p.assignment.ai?.model||p.assignment.origin):p.assignment.origin==='manual'?'人工标注':sourceNames[p.assignment.origin]||p.assignment.origin||'';
 }
};
const originalRenderSpectra=renderSpectra;
renderSpectra=function(){
 originalRenderSpectra();
 window.renderSpectrumCuration();
 const s=spectrum();if(!s)return;
 if(['BMRB_reference','manual_revision'].includes(s.source)){
  const source=s.provenance;
  $('#sourceDetails').innerHTML=`<p>${escape(sourceNames[s.source])} · ${escape(source.entry_id||'')} · ${escape(source.title||'')}</p>`+(source.doi?`<p><a href="https://doi.org/${encodeURIComponent(source.doi)}" target="_blank" rel="noopener">${escape(source.doi)}</a></p>`:'')+`<p>${(source.warnings||[]).map(escape).join('<br>')}</p><pre>${escape(JSON.stringify({reference_smiles:source.reference_smiles,parent_spectrum_id:source.parent_spectrum_id,named_atom_mapping:source.named_atom_mapping,sha256:source.sha256,conditions:s.conditions},null,2))}</pre>`;
 }
};
