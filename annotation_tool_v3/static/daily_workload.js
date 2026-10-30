"use strict";
(() => {
  const host=document.querySelector('#stats');
  let pending=null;
  async function refresh(){
    if(pending)return pending;
    pending=(async()=>{
      const controller=new AbortController(),timer=setTimeout(()=>controller.abort(),15000);
      try{
        const response=await fetch('/api/workload?days=1',{signal:controller.signal});
        const result=await response.json();if(!response.ok)throw new Error(result.error||'读取失败');
        const d=result.day,n=v=>Number(v||0).toLocaleString('zh-CN');
        const metrics=[['processed_molecules','处理分子'],['finished_molecules','结束分子'],['c_peaks','修改 C 峰'],['h_peaks','修改 H 峰'],['save_operations','保存记录']];
        host.innerHTML=`<span title="按北京时间 00:00–24:00 统计">今日 · ${result.date.slice(5)}</span>`+metrics.map(([key,label])=>`<span title="${result.definitions[key]}"><strong>${n(d[key])}</strong>${label}</span>`).join('')+'<a href="/statistics.html" target="_blank" rel="noopener">查看全部统计 ↗</a>';
        host.title=d.recovered_records?'含升级前可恢复记录；完整逐次统计从启用时开始。':'分子和峰均按日去重；保存记录包含自动保存。';
      }catch(error){host.textContent='今日工作量暂时无法读取，可在数据统计页刷新。';}
      finally{clearTimeout(timer);}
    })().finally(()=>{pending=null;});
    return pending;
  }
  window.AnnotationWorkload={refresh};
  refresh();setInterval(refresh,30000);
  window.addEventListener('focus',refresh);
  document.addEventListener('visibilitychange',()=>{if(!document.hidden)refresh();});
})();
