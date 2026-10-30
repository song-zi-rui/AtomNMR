"use strict";
// Edits use the existing revision-guarded management endpoint and save history.
Object.assign(managementActions,{signal_code:'修订碳氢数编码',assignment_group:'建立位移组—原子组对应',remove_assignment_group:'解除位移组'});
$('#manageRevisionPanel').insertAdjacentHTML('beforeend',`
<section class="management-section" id="spectrumRevisionSection">
 <h3>修订 S / D / T / Q 编码</h3>
 <p>保存编码会生成当前数据集的修订版本，原版本保留为参考。S / D / T / Q 分别对应碳上 0 / 1 / 2 / 3 个氢。</p>
 <div class="management-table" id="revisionCodeTable"></div>
 <h3>位移组 ↔ 原子组</h3>
 <p>勾选至少两个位移和两个原子，建立整体对应，组内顺序留空。所有原始位移分别保留，不取平均。相近位移本身不能证明原子化学等价。</p>
 <div class="revision-group-grid"><fieldset><legend>位移组</legend><div id="revisionGroupPeaks"></div></fieldset><fieldset><legend>原子组（当前图编号）</legend><div id="revisionGroupAtoms"></div></fieldset></div>
 <button id="revisionCreateGroup">保存组间对应（组内未定）</button>
 <p>建立分组会撤下组内逐峰确定标签，原值记录在操作历史中。解除分组时可恢复本次分组前的标签；不会覆盖后续修改。</p>
 <div id="revisionExistingGroups"></div>
</section>`);
const revisionStyle=document.createElement('style');
revisionStyle.textContent=`.revision-group-grid{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin:12px 0}.revision-group-grid fieldset{min-width:0;max-height:280px;overflow:auto;border:1px solid #d7e4de;border-radius:7px}.revision-group-grid label{display:flex;align-items:center;gap:9px;margin:8px 0;font-size:16px}.revision-group-grid input[type=checkbox]{width:auto;flex:0 0 auto}#revisionCodeTable select{width:85px;display:inline-block;margin:0}#revisionCodeTable button{white-space:nowrap;margin:2px}#revisionExistingGroups p{background:#eef5f2;padding:12px;line-height:1.8}.revision-code-change{color:#795112;font-weight:600}@media(max-width:650px){.revision-group-grid{grid-template-columns:1fr}}`;
document.head.appendChild(revisionStyle);
function revisionEditable(s){return !!s&&canEditMolecule()&&((s.manual_state==='active'&&s.automatic_quality?.decision!=='quarantined')||(s.source==='manual_revision'&&s.manual_state==='reference'&&!s.provenance?.dataset_version&&!s.provenance?.reference_only_after_revision));}
function renderSpectrumRevision(){
 const s=managedSpectrum(),editable=revisionEditable(s);
 $('#spectrumRevisionSection').hidden=!s;
 if(!s)return;
 $('#revisionCodeTable').innerHTML=s.nucleus!=='13C'?'<p>S/D/T/Q 编码修订仅用于 C 谱。</p>':`<table><thead><tr><th>位移 / ppm</th><th>原库编码</th><th>当前修订 / 检查</th><th>操作</th></tr></thead><tbody>${s.peaks.map(p=>{
  const value=p.signal_code_revision?.code||(p.source_signal_code||'').match(/[SDTQ]$/)?.[0]||'S';
  const note=p.signal_code_revision?(p.effective_signal_code_eligible?'与当前归属原子的氢数一致':'尚无确定原子或氢数仍不一致；暂不作为模型编码'):'未修订';
  return `<tr><td>${escape(p.peak_id)} · ${escape(p.shift_ppm)}</td><td>${escape(p.source_signal_code||'未报告')}</td><td>${escape(p.effective_signal_code||p.source_signal_code||'未报告')}<br><small>${escape(note)}</small></td><td><select aria-label="${escape(p.peak_id)} 修订编码" data-revision-code="${escape(p.peak_id)}" ${editable?'':'disabled'}>${['S','D','T','Q'].map(c=>`<option ${c===value?'selected':''}>${c}</option>`).join('')}</select><button data-save-revision-code="${escape(p.peak_id)}" ${editable?'':'disabled'}>保存编码</button><button data-reset-revision-code="${escape(p.peak_id)}" ${editable&&p.signal_code_revision?'':'disabled'}>恢复原码</button></td></tr>`;
 }).join('')}</tbody></table>`;
 const grouped=new Set((s.assignment_groups||[]).flatMap(g=>g.peak_ids));
 $('#revisionGroupPeaks').innerHTML=s.peaks.map(p=>`<label><input type="checkbox" data-revision-peak="${escape(p.peak_id)}" ${editable&&!grouped.has(p.peak_id)?'':'disabled'}>${escape(p.peak_id)} · ${escape(p.shift_ppm)} ppm${grouped.has(p.peak_id)?'（已有组）':''}</label>`).join('');
 $('#revisionGroupAtoms').innerHTML=state.record.structure.atoms.filter(a=>s.nucleus==='13C'?a.symbol==='C':a.total_hydrogens>0).map(a=>`<label><input type="checkbox" data-revision-atom="${a.index}" ${editable?'':'disabled'}>${escape(a.symbol)}${a.index} · 连氢 ${a.total_hydrogens}</label>`).join('');
 $('#revisionCreateGroup').disabled=!editable;
 $('#revisionExistingGroups').innerHTML=(s.assignment_groups||[]).map(g=>`<p>位移组 {${g.peak_ids.map(id=>escape(s.peaks.find(p=>p.peak_id===id)?.shift_ppm)+' ppm').join('、')}} ↔ 原子组 {${g.atom_indices.map(i=>escape(state.record.structure.atoms[i].symbol)+i).join('、')}}<br>组内顺序未定 · ${escape(g.reason)}<br><button data-remove-revision-group="${escape(g.id)}" ${editable?'':'disabled'}>解除此组并恢复此前标签</button></p>`).join('')||'<p>尚未建立位移组。</p>';
}
const inspectBeforeRevision=manageInspect;
manageInspect=function(){inspectBeforeRevision();renderSpectrumRevision();};
$('#revisionCodeTable').onclick=e=>{
 const saveButton=e.target.closest('[data-save-revision-code]'),resetButton=e.target.closest('[data-reset-revision-code]');
 if(!saveButton&&!resetButton)return;
 const pid=saveButton?.dataset.saveRevisionCode||resetButton.dataset.resetRevisionCode;
 const selector=$$('[data-revision-code]').find(x=>x.dataset.revisionCode===pid);
 managementAction({action:'signal_code',peak_id:pid,code:saveButton?selector.value:null});
};
$('#revisionCreateGroup').onclick=()=>{
 const peakIds=$$('[data-revision-peak]:checked').map(x=>x.dataset.revisionPeak),atoms=$$('[data-revision-atom]:checked').map(x=>Number(x.dataset.revisionAtom));
 if(peakIds.length<2||atoms.length<2)return toast('请至少勾选两个位移和两个原子。');
 managementAction({action:'assignment_group',peak_ids:peakIds,atom_indices:atoms});
};
$('#revisionExistingGroups').onclick=e=>{const button=e.target.closest('[data-remove-revision-group]');if(button)managementAction({action:'remove_assignment_group',group_id:button.dataset.removeRevisionGroup});};
const peaksBeforeRevision=renderPeaks;
renderPeaks=function(){
 peaksBeforeRevision();const s=spectrum();if(!s)return;
 const revised=s.peaks.filter(p=>p.signal_code_revision);
 const warning=$('#assignmentWarning');
 if(warning&&revised.length)warning.textContent+=' 当前已有 '+revised.length+' 个峰保存了编码修订；采用值见峰表“修订编码”，原库警告继续保留作来源记录。';
 for(const row of $$('#peaks tr[data-peak]')){
  const p=s.peaks.find(x=>x.peak_id===row.dataset.peak);if(!p)continue;
  if(p.signal_code_revision){const detail=row.querySelector('.detail small');if(detail){detail.textContent='修订编码 '+p.effective_signal_code+'（原库 '+(p.source_signal_code||'未报告')+'）';detail.classList.add('revision-code-change');}}
  if(p.assignment_group){const target=row.querySelector('.target'),status=row.querySelector('td:last-child');if(target){target.textContent='原子组 '+p.assignment_group.atom_indices.map(i=>state.record.structure.atoms[i].symbol+i).join(' / ');target.title='共同位移组：'+p.assignment_group.peak_ids.map(id=>s.peaks.find(x=>x.peak_id===id)?.shift_ppm).join(' / ')+' ppm；组内顺序未定。';}if(status)status.textContent='位移组对应 · 组内未定';}
 }
};
const editorBeforeRevision=renderEditor;
renderEditor=function(){
 editorBeforeRevision();const p=peak();if(!p?.assignment_group)return;
 $('#selectedAtoms').textContent='位移组整体对应原子组；组内顺序未定。请在管理当前谱中查看或解除。';
 ['#status','#clearAtoms','#atomNumber','#atomApply'].forEach(id=>$(id).disabled=true);$('#peakNote').readOnly=true;
};
