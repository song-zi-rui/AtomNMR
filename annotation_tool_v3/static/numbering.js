"use strict";
// Numbering review is intentionally separate from normal annotation autosave.
$('#spectrumTools').insertAdjacentHTML('afterend','<div id="numberingNotice" class="numbering-notice" hidden><button id="numberingOpen">原子编号检查 / 恢复候选</button><span id="numberingStatus"></span></div>');
document.body.insertAdjacentHTML('beforeend',`<dialog id="numberingDialog" class="numbering-dialog">
 <button id="numberingClose" class="dialogclose" aria-label="关闭编号检查">×</button>
 <h2>源编号 → 原始原子 → 当前图中编号</h2>
 <p id="numberingSummary" role="status"></p>
 <label>查看编号解释<select id="numberingOption"></select></label>
 <p id="numberingChecks"></p>
 <div class="numbering-table"><table><thead><tr><th>δ / ppm</th><th>源 ID</th><th>原始映射</th><th>所选解释对应原子<br>原始 idx → 当前图</th><th>原库编码 / 警告</th></tr></thead><tbody id="numberingRows"></tbody></table></div>
 <details><summary>完整原子身份、InChI 和检查细节</summary><pre id="numberingDetails"></pre></details>
 <p class="numbering-caution">检查通过只证明结构兼容，不证明历史实验归属正确。请核对原条目或论文原子标签；不能仅凭位移合理就确认。采用只影响当前谱，人工标注不会被覆盖，S/D 编码警告不会自动修好。</p>
 <label>原文位置、编号依据或拒绝原因<textarea id="numberingNote" rows="3" placeholder="例如：核对原条目的结构原子标签及全部谱峰，说明为何采用这套编号（至少 8 字）"></textarea></label>
 <label><input type="checkbox" id="numberingSourceConfirmed">已核对独立来源，确认该谱采用所选编号，而不只是化学位移看起来合理</label>
 <label><input type="checkbox" id="numberingMapConfirmed">已逐峰核对映射；理解同位移多原子不自动代表等价</label>
 <div class="numbering-actions"><button id="numberingAccept">确认采用当前编号解释</button><button id="numberingReject">拒绝此解释，保留原记录</button></div>
 <h3>编号操作记录</h3><div id="numberingHistory"></div>
</dialog>`);
let numberingPreview=null,numberingBusy=false;
const numberingReasons={source_inchi_missing_or_not_exact_equal:'原始 InChI 缺失或不完全一致',candidate_same_as_original_order:'与原始顺序相同，不是新的恢复方案',candidate_structural_checks_failed:'元素、H 容量或原库编码检查未通过',candidate_has_empty_assignments:'存在空归属'};
const numberingIssueNames={carbon_assignment_to_noncarbon:'C 谱指向非 C',carbon_hcount_code_conflict:'S/D/T/Q 与结构氢数不一致',proton_assignment_without_h_capacity:'H 谱指向无氢容量的原子',source_atom_out_of_range:'源编号越界'};
function numberingOptions(){return numberingPreview?[numberingPreview.default,...numberingPreview.candidates]:[];}
function selectedNumbering(){return numberingOptions().find(o=>o.candidate_id===$('#numberingOption').value);}
window.renderNumberingNotice=()=>{
 const s=spectrum();$('#numberingNotice').hidden=!s||s.source!=='nmrshiftdb2';
 if(!s)return;
 const v=s.source_assignment_validation;
 $('#numberingStatus').textContent=v?.status==='identity_validated'?'源原子身份检查通过；不等于实验归属已验证。':v?.status==='needs_attention'?'编号与编码分开检查；警告不一定是编号错位。':'源编号有冲突或尚未确认，可查看恢复候选。';
 $('#numberingOpen').disabled=numberingBusy;
};
function renderNumberingOption(){
 const o=selectedNumbering();if(!o)return;
 $('#numberingSourceConfirmed').checked=false;$('#numberingMapConfirmed').checked=false;
 const issues=o.validation.issues||[];
 const warnings=issues.map(i=>(numberingIssueNames[i.code]||i.code)+' · 源 ID '+(i.source_assignment_id??'—'));
 $('#numberingChecks').textContent=(o.adoptable?'可以人工确认采用，但不是自动修复。':'当前解释不能直接采用。')+(warnings.length?' '+warnings.join('；'):' 元素及已有编码的结构兼容检查通过。')+(o.blocked_reasons?.length?' 阻止原因：'+o.blocked_reasons.map(x=>numberingReasons[x]||x).join('；'):'');
 const original=new Map(numberingPreview.default.peaks.flatMap(p=>p.source_assignment_fragments.flatMap(f=>f.assignments.map(a=>[f.source_fragment_ordinal+':'+a.source_assignment_id,a]))));
 $('#numberingRows').innerHTML=o.peaks.flatMap(p=>p.source_assignment_fragments.flatMap(f=>f.assignments.map(a=>{
  const before=original.get(f.source_fragment_ordinal+':'+a.source_assignment_id);
  const different=before?.display_atom_label!==a.display_atom_label;
  const problems=f.issues.map(i=>numberingIssueNames[i.code]||i.code).join('；');
  return `<tr class="${different?'numbering-difference':''}"><td>${escape(f.shift_text)}</td><td>${escape(a.source_assignment_id)}</td><td>${escape(before?.original_atom_idx??'—')} → ${escape(before?.display_atom_label||'无法映射')}</td><td>${escape(a.original_atom_idx??'—')} → ${escape(a.display_atom_label||'无法映射')}</td><td>${escape(f.source_signal_code||'未报告')}<br><small>${escape(problems)}</small></td></tr>`;
 }))).join('');
 $('#numberingDetails').textContent=JSON.stringify({source:numberingPreview.source,original_identity:numberingPreview.atom_identity,basis:o.basis,generated_inchi:o.generated_inchi,aux_info:o.aux_info,candidate_identity:o.identity,validation:o.validation},null,2);
 updateNumberingAccept();
}
function updateNumberingAccept(){
 $('#numberingAccept').disabled=numberingBusy||!canEditSpectrum()||!selectedNumbering()?.adoptable||!$('#numberingSourceConfirmed').checked||!$('#numberingMapConfirmed').checked||$('#numberingNote').value.trim().length<8;
 $('#numberingReject').disabled=numberingBusy||!canEditSpectrum();
}
function renderNumbering(){
 const p=numberingPreview;
 $('#numberingSummary').textContent=`${state.record.formula} · ${p.nucleus} · 源谱 ${p.spectrum_id}。先按原始可见 molfile 编号检查；${p.fallback_attempted?'有异常，已尝试 InChI 编号解释。':'原始检查通过，不需要切换编号。'}${p.fallback_error?' 候选生成失败：'+p.fallback_error:''}`;
 $('#numberingOption').innerHTML=numberingOptions().map((o,i)=>`<option value="${escape(o.candidate_id)}">${i===0?'原始 molfile 编号（默认）':'InChI 规范编号（恢复假设）'} · ${o.adoptable?'待人工确认':'不能直接采用 / 无需变更'}</option>`).join('');
 const undone=new Set(p.history.filter(h=>h.decision==='undo').map(h=>h.target_review_id));
 $('#numberingHistory').innerHTML=p.history.map(h=>`<p>${escape(dateLabel(h.created_at))} · ${escape({confirm:'采用',reject:'拒绝',undo:'撤销'}[h.decision])}<br>${escape(h.note)} ${h.decision==='confirm'&&!undone.has(h.id)?`<button data-numbering-undo="${h.id}" ${!canEditSpectrum()?'disabled':''}>撤销这次采用（保留后续人工修改）</button>`:''}</p>`).join('')||'<p>尚未进行编号审核。</p>';
 renderNumberingOption();
}
async function openNumbering(){
 if(numberingBusy||!state.record||spectrum()?.source!=='nmrshiftdb2')return;
 numberingBusy=true;window.renderNumberingNotice();
 try{
  await save();const mid=state.record.id,sid=state.sid;
  const result=await api(`/api/molecule/${mid}/numbering?`+new URLSearchParams({spectrum_id:sid}));
  if(state.record.id!==mid||state.sid!==sid)return;
  numberingPreview=result;$('#numberingNote').value='';renderNumbering();$('#numberingDialog').showModal();
 }catch(e){toast(e.message);}finally{numberingBusy=false;window.renderNumberingNotice();updateNumberingAccept();}
}
async function submitNumbering(decision,reviewId=null){
 if(numberingBusy||!numberingPreview)return;
 const note=$('#numberingNote').value.trim();if(note.length<8)return toast('请填写至少 8 字的依据或原因。');
 if(decision==='confirm'&&(!selectedNumbering()?.adoptable||!$('#numberingSourceConfirmed').checked||!$('#numberingMapConfirmed').checked))return toast('请完成独立来源与逐峰对应确认。');
 numberingBusy=true;updateNumberingAccept();
 try{
  const p=numberingPreview;
  const result=await api(`/api/molecule/${p.molecule_id}/numbering`,{revision:p.revision,spectrum_id:p.spectrum_id,candidate_id:selectedNumbering()?.candidate_id,decision,review_id:reviewId,note,confirmed_source:$('#numberingSourceConfirmed').checked,confirmed_mapping:$('#numberingMapConfirmed').checked});
  if(state.record.id===p.molecule_id){state.record=result;state.dirty=false;state.pid=null;renderHeader();renderSpectra();refreshProgressSoon();}
  $('#numberingDialog').close();toast(`编号操作已保存；更改 ${result.numbering_result.changed_peaks} 个展示峰，保留 ${result.numbering_result.preserved_manual_or_later_edits} 个已有人工或后续修改。`);
  void loadList().catch(()=>{});
 }catch(e){toast(e.message);}finally{numberingBusy=false;updateNumberingAccept();window.renderNumberingNotice();}
}
$('#numberingOpen').onclick=openNumbering;
$('#numberingClose').onclick=()=>{if(!numberingBusy)$('#numberingDialog').close();};
$('#numberingDialog').addEventListener('cancel',e=>{if(numberingBusy)e.preventDefault();});
$('#numberingOption').onchange=renderNumberingOption;
$('#numberingSourceConfirmed').onchange=updateNumberingAccept;
$('#numberingMapConfirmed').onchange=updateNumberingAccept;
$('#numberingNote').oninput=updateNumberingAccept;
$('#numberingAccept').onclick=()=>submitNumbering('confirm');
$('#numberingReject').onclick=()=>submitNumbering('reject');
$('#numberingHistory').onclick=e=>{const button=e.target.closest('[data-numbering-undo]');if(button)submitNumbering('undo',Number(button.dataset.numberingUndo));};
