'use strict';
const $=id=>document.getElementById(id), esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
let timer,busy=false;
async function api(path,body){const r=await fetch(path,body?{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)}:{});const v=await r.json();if(!r.ok)throw Error(v.error||'操作失败');return v;}
function message(text){$('message').hidden=!text;$('message').textContent=text;}
function payload(){return {webdav_url:$('url').value.trim(),username:$('username').value.trim(),password:$('password').value,enabled:$('enabled').checked,interval_hours:Number($('hours').value),cleanup_uploaded:$('cleanup').checked,github_repo:$('repo').value.trim(),github_enabled:$('github').checked};}
async function save(){await api('/api/backups/settings',payload());$('password').value='';}
const phases={queued:'等待开始',snapshot:'正在生成数据库快照',checksum:'正在计算校验值',upload:'正在上传',github:'正在备份程序',done:'完成'};
async function load(fill=false){clearTimeout(timer);const r=await api('/api/backups'),s=r.settings;
 if(fill){$('url').value=s.webdav_url;$('username').value=s.username;$('enabled').checked=s.enabled;$('hours').value=s.interval_hours;$('cleanup').checked=s.cleanup_uploaded;$('repo').value=s.github_repo;$('github').checked=s.github_enabled;}
 $('connection').textContent=s.password_saved?'配置已保存（可检查连接）':'尚未配置';$('note').textContent=r.note;
 const active=r.jobs.find(j=>j.state==='running');$('run').disabled=r.active||busy;
 $('progress').textContent=active?(phases[active.phase]||active.phase)+(active.current_file?' · '+active.current_file:'')+(active.phase==='snapshot'&&active.pages_total?' · '+Math.floor(100*active.pages_done/active.pages_total)+'%':''):'';
 $('history').innerHTML=r.jobs.map(j=>`<tr><td>${esc(new Date(j.started_at).toLocaleString('zh-CN'))}<br><small>${esc(j.snapshot_id)}</small></td><td>${j.state==='running'?esc(phases[j.phase]||j.phase):j.state==='done'?'已完成':j.state==='interrupted'?'服务中断':'未全部完成'}<br>${j.kind==='code'?'仅程序代码':j.cloud_complete?'云端完成':'云端未确认完成'}${j.github_complete?' · GitHub完成':''}${j.error?'<br><span class="error">'+esc(j.error)+'</span>':''}</td><td>${((j.uploaded_bytes||0)/1024**3).toFixed(2)} GB</td><td>${!r.active&&j.can_retry?'<button data-retry="'+esc(j.id)+'">重新上传此快照</button>':''}</td></tr>`).join('')||'<tr><td colspan="4">尚无备份记录</td></tr>';
 timer=setTimeout(()=>load().catch(e=>message(e.message)),r.active?2500:15000);
}
function action(id,fn){$(id).onclick=async()=>{if(busy)return;busy=true;$(id).disabled=true;try{message('');await fn();await load();}catch(e){message(e.message);}finally{busy=false;$(id).disabled=false;}};}
action('save',async()=>{await save();message('设置已保存。尚未配置云端时不会自动上传。');});
action('test',async()=>{await save();const r=await api('/api/backups/test',{});message(r.message);});
action('run',async()=>{await save();await api('/api/backups/run',{});message('已开始后台备份，可以关闭本页面。');});
action('code',async()=>{await save();await api('/api/backups/code',{});message('已开始单独备份程序代码；不生成数据库快照。');});
action('refresh',()=>load());
$('history').onclick=async e=>{const b=e.target.closest('[data-retry]');if(!b||busy)return;b.disabled=true;try{await api('/api/backups/run',{retry_id:b.dataset.retry});await load();}catch(x){message(x.message);}finally{b.disabled=false;}};
load(true).catch(e=>message(e.message));
