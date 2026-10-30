"use strict";
Object.assign(managementActions,{supplement_peak:'补充实验峰表',withdraw_peak_supplement:'撤销峰表修订',edit_peak:'修改实验峰',delete_peak:'删除实验峰'});
$('#managePeaksPanel').insertAdjacentHTML('beforeend',`
<section class="management-section" id="peakSupplementSection">
 <h3>实验峰表增删改</h3>
 <p>新增峰先留空归属；保存后可回主界面标注。操作原因填写在窗口底部。</p>
 <label>补充目标<select id="supplementTarget"><option value="">新增实验峰</option></select></label>
 <label>已有峰操作<select id="supplementMode"><option value="fill">补充缺失字段</option><option value="edit">修改现有数值</option></select></label>
 <div class="supplement-fields">
  <label>化学位移 / ppm<input id="supplementShift" type="number" step="any" placeholder="新增峰必填，如 177.2"></label>
  <label>峰型<input id="supplementMultiplicity" placeholder="如 s、dd、m、br s"></label>
  <label>J 值 / Hz<input id="supplementJ" placeholder="如 7.2, 1.5"></label>
  <label>文献明确报告的 H 数<input id="supplementH" type="number" step="any" min="0" placeholder="未报告可留空"></label>
  <label>相对积分<input id="supplementIntegral" type="number" step="any" min="0" placeholder="与文献 H 数分开记录"></label>
  <label>线宽 / Hz<input id="supplementWidth" type="number" step="any" min="0"></label>
  <label>位移范围下限 / ppm<input id="supplementLow" type="number" step="any"></label>
  <label>位移范围上限 / ppm<input id="supplementHigh" type="number" step="any"></label>
 </div>
 <label class="inline"><input id="supplementSameShift" type="checkbox">如已有相同位移，仍保留为独立实验信号</label>
 <button id="supplementSave">保存峰表修订并作为当前数据集</button><button id="supplementDelete">删除所选峰（原版本保留参考）</button>
 <p>补缺模式会锁定已有数值；修改模式允许更正位移及其他字段。修改和删除保存到同一份最终版，原始版保留供对照。</p>
 <div id="supplementHistory"></div>
</section>`);
const supplementFields={shift_ppm:'#supplementShift',multiplicity:'#supplementMultiplicity',j_values:'#supplementJ',reported_nH:'#supplementH',integral:'#supplementIntegral',line_width_hz:'#supplementWidth',range_min:'#supplementLow',range_max:'#supplementHigh'};
const supplementStyle=document.createElement('style');
supplementStyle.textContent='.supplement-fields{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px 16px;margin:12px 0}.supplement-fields label{margin:0}.supplement-fields input:disabled{background:#f1f4f3;color:#63716c}#supplementHistory p{padding:10px;background:#edf5f1;overflow-wrap:anywhere}#supplementHistory button{margin-left:8px}@media(max-width:650px){.supplement-fields{grid-template-columns:1fr}}';
document.head.appendChild(supplementStyle);
let supplementContext=null;
function resetSupplementFields(){
  for(const selector of Object.values(supplementFields))$(selector).value='';
  $('#supplementSameShift').checked=false;
}
function renderPeakSupplement(reset=false){
  const s=managedSpectrum();$('#peakSupplementSection').hidden=!s;if(!s)return;
  const context=state.record.id+':'+s.id;
  if(context!==supplementContext){reset=true;supplementContext=context;$('#supplementTarget').value='';$('#supplementMode').value='fill';}
  const selected=$('#supplementTarget').value;
  $('#supplementTarget').innerHTML='<option value="">新增实验峰</option>'+s.peaks.map(p=>`<option value="${escape(p.peak_id)}">${escape(p.peak_id)} · ${escape(p.shift_ppm)} ppm：补充缺失字段</option>`).join('');
  $('#supplementTarget').value=s.peaks.some(p=>p.peak_id===selected)?selected:'';
  if(reset)resetSupplementFields();
  const peak=s.peaks.find(p=>p.peak_id===$('#supplementTarget').value),editable=revisionEditable(s),editing=!!peak&&$('#supplementMode').value==='edit';
  $('#supplementMode').disabled=!editable||!peak;
  for(const [key,selector] of Object.entries(supplementFields)){
    const present=peak&&peak[key]!=null&&peak[key]!==''&&!(Array.isArray(peak[key])&&!peak[key].length);
    $(selector).disabled=!editable||!!present&&!editing;
    if(present)$(selector).value=Array.isArray(peak[key])?peak[key].join(', '):peak[key];
  }
  $('#supplementSameShift').disabled=!editable||!!peak&&!editing;
  $('#supplementSave').disabled=!editable;
  $('#supplementDelete').disabled=!editable||!peak;
  $('#supplementHistory').innerHTML=(s.provenance.peak_table_supplements||[]).map(x=>`<p>${x.active!==false?'当前修订':'已撤销'} · ${escape(x.peak_id)} · ${escape({add:'新增实验峰',fields:'补充字段',edit:'修改数值',delete:'删除峰'}[x.kind])}<br>${escape(JSON.stringify(x.fields))}<br>${escape(x.reason)}${x.active!==false?`<button data-withdraw-supplement="${escape(x.id)}" ${editable?'':'disabled'}>撤销这次修订（更新最终版）</button>`:''}</p>`).join('');
}
const inspectBeforePeakSupplement=manageInspect;
manageInspect=function(){inspectBeforePeakSupplement();renderPeakSupplement();};
$('#supplementTarget').onchange=()=>renderPeakSupplement(true);
$('#supplementMode').onchange=()=>renderPeakSupplement(true);
$('#supplementSave').onclick=async()=>{
  const fields={};
  const selected=managedSpectrum()?.peaks.find(p=>p.peak_id===$('#supplementTarget').value);
  const editing=!!selected&&$('#supplementMode').value==='edit';
  for(const [key,selector] of Object.entries(supplementFields)){
    if($(selector).disabled)continue;
    const text=String($(selector).value).trim();
    if(!text){if(editing&&key!=='shift_ppm'&&selected[key]!=null)fields[key]=null;continue;}
    if(key==='multiplicity')fields[key]=text;
    else if(key==='j_values')fields[key]=text.split(/[,，;；\s]+/).filter(Boolean).map(Number);
    else fields[key]=Number(text);
  }
  if(editing)for(const key of Object.keys(fields))if(JSON.stringify(fields[key])===JSON.stringify(selected[key]??null))delete fields[key];
  if(!Object.keys(fields).length)return toast('请填写需要补充的实验峰表数值。');
  if(Object.values(fields).some(v=>typeof v==='number'&&!Number.isFinite(v)||Array.isArray(v)&&v.some(n=>!Number.isFinite(n))))return toast('请检查数值格式，J 值使用逗号分隔。');
  const result=await managementAction({action:editing?'edit_peak':'supplement_peak',peak_id:$('#supplementTarget').value||null,fields,allow_same_shift:$('#supplementSameShift').checked});
  if(result)renderPeakSupplement(true);
};
$('#supplementDelete').onclick=()=>{const pid=$('#supplementTarget').value;if(pid)managementAction({action:'delete_peak',peak_id:pid});};
$('#supplementHistory').onclick=e=>{const b=e.target.closest('[data-withdraw-supplement]');if(b)managementAction({action:'withdraw_peak_supplement',supplement_id:b.dataset.withdrawSupplement});};
const renderBeforePeakSupplement=renderPeaks;
renderPeaks=function(){
  renderBeforePeakSupplement();const s=spectrum();if(!s)return;
  for(const row of $$('#peaks tr[data-peak]')){
    const p=s.peaks.find(p=>p.peak_id===row.dataset.peak);
    if(p?.experimental_peak_added){const cell=row.querySelector('td');if(cell){cell.innerHTML=escape(p.peak_id)+'<small style="display:block;color:#356957">实验补录</small>';cell.title='人工补录的实验峰；来源与撤销操作见管理当前谱。';}}
  }
};
