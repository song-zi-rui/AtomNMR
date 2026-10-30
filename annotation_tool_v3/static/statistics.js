"use strict";
(() => {
  const $=id=>document.getElementById(id),n=value=>Number(value||0).toLocaleString('zh-CN');
  const text=value=>String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const stamp=value=>new Date(value).toLocaleString('zh-CN',{timeZone:'Asia/Shanghai',hour12:false});
  let today=null,lastDaily=null,lastTotal=null,sequence=0;
  async function get(path){const controller=new AbortController(),timer=setTimeout(()=>controller.abort(),60000);try{const response=await fetch(path,{signal:controller.signal});const value=await response.json();if(!response.ok)throw Error(value.error||'请求失败');return value;}finally{clearTimeout(timer);}}
  const card=(title,value,unit='')=>`<article class="card"><p>${text(title)}</p><strong>${n(value)}<span class="unit">${text(unit)}</span></strong></article>`;
  function renderDaily(data){
    const d=data.day;today=data.today;$('day').value=data.date;$('day').max=today;
    $('dailyTitle').textContent=data.date===today?'今日工作量':data.date+' 工作量';
    $('dailyCards').innerHTML=card('处理分子',d.processed_molecules,'个')+card('结束分子',d.finished_molecules,'个')+card('修改 C 峰',d.c_peaks,'个')+card('修改 H 峰',d.h_peaks,'个')+card('保存记录',d.save_operations,'次')+card('管理操作',d.management_operations,'次');
    $('autoDailyCards').innerHTML=card('AI 结果入库分子',d.ai_molecules,'个')+card('AI 新增 C 归属',d.ai_c_peaks,'峰')+card('AI 新增 H 归属',d.ai_h_peaks,'峰')+card('AI 入库结果',d.ai_operations,'条');
    $('historyNote').textContent=`完整逐次记录从 ${data.recording_started_at?stamp(data.recording_started_at):'本次升级'} 开始。更早日期按已有最后处理时间和管理记录恢复，可能遗漏早先保存及结束记录；逐峰工作量从启用后开始统计。${d.recovered_records?' 当天含 '+n(d.recovered_records)+' 条恢复记录。':''}`;
    const max=Math.max(1,...data.history.map(x=>x.processed_molecules));
    $('historyRows').innerHTML=data.history.map(r=>`<tr class="${r.date===data.date?'selected-day':''}"><td><button class="date-link" data-date="${r.date}">${r.date}${r.date===today?' · 今天':''}</button></td><td>${n(r.processed_molecules)}<div class="mini-bar" style="width:${(r.processed_molecules/max*100).toFixed(1)}%"></div></td><td>${n(r.finished_molecules)}</td><td>${n(r.c_peaks)}</td><td>${n(r.h_peaks)}</td><td>${n(r.save_operations)}</td><td>${n(r.management_operations)}</td></tr>`).join('');
  }
  function renderTotal(s){
    const quality=s.spectrum_quality||{},held=quality.quarantined||0,p=s.annotation_processing||{},progress=s.dual_assignment_progress||{};
    const origins=s.annotation_origins||{sources:[],models:[]},sourceNames={shiftdb2:'ShiftDB2 原库',manual:'人工标注',ai:'AI 自动标注',other:'其他来源'};
    $('annotationSources').innerHTML=origins.sources.map(r=>`<tr><td>${sourceNames[r.source]||text(r.source)}</td><td>${r.nucleus}</td><td>${n(r.peaks)}</td><td>${n(r.molecules)}</td></tr>`).join('')||'<tr><td colspan="4">暂无来源统计</td></tr>';
    $('annotationModels').innerHTML=origins.models.map(r=>`<tr><td>${text(r.model)}<br><small>${text(r.provider)}</small></td><td>${n(r.results)}</td><td>${n(r.molecules)}</td><td>${n(r.current_ai_peaks)}</td><td>${n(r.published_peaks)}</td></tr>`).join('')||'<tr><td colspan="5">尚无 AI 结果入库</td></tr>';
    const finished=Object.values(p).reduce((sum,v)=>sum+v.finished,0);
    $('totalCards').innerHTML=card('存档分子总量',s.molecules,'个')+card('存档实验谱总量',s.spectra,'条')+card('未自动隔离的存档谱',s.spectra-held,'条')+card('自动隔离谱',held,'条')+card('本轮已结束的分子',finished,'个')+card('C/H 双谱全归属分子',progress.ready_molecules,'个');
    $('coverageRows').innerHTML=[['full','全归属'],['partial','部分归属'],['none','无归属']].map(([key,label])=>{const v=p[key]||{pending:0,finished:0};return `<tr><td>${label}</td><td>${n(v.pending)}</td><td>${n(v.finished)}</td><td>${n(v.pending+v.finished)}</td></tr>`;}).join('');
    const matrix=s.spectrum_coverage;
    if(matrix){
      const names={dual:'H/C 双谱',c_only:'只有 C 谱',h_only:'只有 H 谱',no_current:'无当前可用谱'};
      $('spectrumCoverageRows').innerHTML=matrix.rows.map(row=>`<tr class="${row.type==='no_current'?'no-current':''}"><th scope="row">${names[row.type]||text(row.type)}</th>${['full','partial','none','total'].map(key=>`<td>${n(row[key])}</td>`).join('')}</tr>`).join('');
      $('spectrumCoverageTotal').innerHTML=`<tr><th scope="row">合计</th>${['full','partial','none','total'].map(key=>`<td>${n(matrix.total[key])}</td>`).join('')}</tr>`;
      const absent=matrix.rows.find(row=>row.type==='no_current')?.total||0;
      $('spectrumCoverageNote').textContent=`全归属：当前全部峰有归属；部分归属：有归属峰但尚未覆盖全部峰；完全未归属：没有确定归属的峰。${absent?' 另有 '+n(absent)+' 个分子无当前可用谱，单独列出并计入无归属列，以与存档总量核对。':''} 当前谱不含已隔离、已停用、参考或历史版本；分子层面停用的记录仍按保留谱分类。`;
    }else{
      $('spectrumCoverageRows').innerHTML='<tr><td colspan="5">谱类型统计暂不可用，请刷新重试。</td></tr>';
      $('spectrumCoverageTotal').innerHTML='';
    }
    $('dualProgress').innerHTML=`<div class="progress-number">${n(progress.ready_molecules)} <small>/ ${n(progress.total_dual_molecules)} 个</small></div><progress aria-label="C/H 双谱全归属比例" max="${Math.max(1,progress.total_dual_molecules||0)}" value="${progress.ready_molecules||0}"></progress><p>${Number(progress.percent||0).toFixed(2)}% 已全归属 · 还剩 ${n(progress.remaining_molecules)} 个双谱分子</p>`;
  }
  async function refresh(totals=true){
    const request=++sequence;$('refresh').disabled=true;$('updated').textContent='正在更新…';
    const params=new URLSearchParams({days:14});if($('day').value)params.set('date',$('day').value);
    const jobs=[get('/api/workload?'+params).then(data=>{if(request===sequence){lastDaily=data;renderDaily(data);}})];
    if(totals)jobs.push(get('/api/stats').then(data=>{if(request===sequence){lastTotal=data;renderTotal(data);}}));
    const results=await Promise.allSettled(jobs);if(request!==sequence)return;
    const error=results.find(x=>x.status==='rejected');$('updated').classList.toggle('error',!!error);
    $('updated').textContent=error?'部分统计读取失败，请刷新重试。':`更新于 ${stamp(new Date())} · 北京时间`;
    if(!lastDaily&&results[0].status==='rejected')$('dailyCards').textContent='工作记录暂时无法读取，请刷新重试。';
    if(!lastTotal&&totals&&results[1].status==='rejected')$('totalCards').textContent='数据库总量暂时无法读取，请刷新重试。';
    $('refresh').disabled=false;
  }
  $('refresh').onclick=()=>refresh();$('day').onchange=()=>{if($('day').value)refresh(false);};
  $('today').onclick=()=>{$('day').value='';refresh(false);};
  $('historyRows').onclick=event=>{const button=event.target.closest('[data-date]');if(button){$('day').value=button.dataset.date;refresh(false);}};
  refresh();setInterval(()=>{if(!document.hidden){if($('day').value===today)$('day').value='';refresh();}},60000);
})();
