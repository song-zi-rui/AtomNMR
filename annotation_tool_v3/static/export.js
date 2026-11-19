'use strict';
const $=id=>document.getElementById(id),esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const params=new URLSearchParams(location.search);
let query='';
try{query=JSON.parse(params.get('filters')||'{}')?.q||'';}catch{}
$('query').value=params.get('molecule_id')||query;
let job=null,poll=null,working=false,catalog={};
const sourceNames={IBM_DFT:'IBM', 'CASCADE-DFT8K':'CASCADE DFT8K', 'CASPRE-JEOL':'CASPRE · JEOL',
  nmrshiftdb2:'NMRShiftDB2',manual_revision:'人工修订谱',literature_pdf:'文献实验谱'};
const coverageNames={all:'全部归属情况',full:'全部归属',partial:'部分归属',none:'无归属'};
const nucleiNames={all:'全部核种组合',dual:'全双谱',proton_only:'仅 H 谱',carbon_only:'仅 C 谱',none:'无当前可用谱'};
async function api(path,body){
 const r=await fetch(path,body===undefined?{}:{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
 const x=await r.json();if(!r.ok)throw Error(x.error||'请求失败');return x;
}
function notice(text){$('message').textContent=text;$('message').hidden=!text;}
function updateSources(){
 const s=$('spectrumSource'),previous=s.value,available=catalog[$('spectrumType').value]||[];
 s.innerHTML='<option value="">当前类型的全部来源</option>'+available.map(v=>'<option value="'+esc(v)+'">'+esc(sourceNames[v]||v)+'</option>').join('');
 s.value=available.includes(previous)?previous:'';
}
function options(){return {
 spectrum_type:$('spectrumType').value,spectrum_source:$('spectrumSource').value,
 coverage:$('coverage').value,nuclei:$('nuclei').value,
 equivalence:$('equivalence').disabled?'none':$('equivalence').value,
 require_same_solvent:!$('sameSolvent').disabled&&$('sameSolvent').checked
};}
function payload(){
 const p={format:$('format').value,filters:{q:$('query').value.trim()},export_options:options()};
 if(p.format==='atompeak')p.training_options={seed:Number($('seed').value),valid_ratio:Number($('validRatio').value),
  test_ratio:Number($('testRatio').value),max_seq_len:Number($('maxSeq').value),use_solvent:$('useSolvent').checked};
 return p;
}
const controls=['preview','run','spectrumType','spectrumSource','query','format','coverage','nuclei','equivalence',
 'sameSolvent','seed','maxSeq','validRatio','testRatio','useSolvent'];
function describe(){
 const experimental=$('spectrumType').value==='experimental',single=['carbon_only','proton_only'].includes($('nuclei').value);
 if(!experimental||single)$('sameSolvent').checked=false;
 $('sameSolvent').disabled=working||!experimental||single;
 $('equivalence').disabled=working||!experimental;
 const o=options(),training=$('format').value==='atompeak';
 const parts=[experimental?'实验谱':'模拟谱',sourceNames[o.spectrum_source]||o.spectrum_source||'全部来源',
  coverageNames[o.coverage],nucleiNames[o.nuclei]];
 if($('equivalence').value!=='none'&&experimental)parts.push('筛除 ShiftDB2 等价原子疑似问题');
 if($('sameSolvent').checked)parts.push('已知同溶剂谱对');
 $('summary').textContent=parts.join(' · ');
 $('trainingOptions').hidden=!training;
 $('coverageHelp').textContent=(training?'训练包每个分子各选最多一条 H 谱和 C 谱，按实际选入的谱判定。':'只按本次所选类型和来源内的谱判定。')+
  ' 全部归属要求各峰有确定的原子集合。歧义峰和弱标签不计入确定归属。';
}
for(const id of ['spectrumSource','coverage','nuclei','equivalence','sameSolvent','format'])$(id).onchange=describe;
$('spectrumType').onchange=()=>{updateSources();describe();};
const reasonNames={coverage_full:'实际为全部归属',coverage_partial:'实际为部分归属',coverage_none:'实际为无归属',
 no_eligible_spectra:'无符合条件的当前可用谱',missing_paired_nucleus:'缺少H谱或C谱',no_same_known_solvent_pair:'没有已知同溶剂谱对',
 suspected_source_equivalence:'ShiftDB2 等价原子疑似问题',
 ...Object.fromEntries(Object.entries(nucleiNames).map(([k,v])=>['nuclei_'+k,'实际为'+v]))};
function busy(value){working=value;for(const id of controls)$(id).disabled=value;describe();}
function showSourceCounts(counts){
 $('sourceBreakdown').hidden=!counts||!Object.keys(counts).length;
 $('sourceRows').innerHTML=Object.entries(counts||{}).map(([name,s])=>
  '<tr><td>'+esc(name.split(' + ').map(v=>sourceNames[v]||v).join(' + '))+'</td><td>'+Number(s.train).toLocaleString()+
  '</td><td>'+Number(s.valid).toLocaleString()+'</td><td>'+Number(s.test).toLocaleString()+'</td></tr>').join('');
}
function clearSavedJob(){try{localStorage.removeItem('atomnmr.export.pending');}catch{}}
async function refresh(){
 clearTimeout(poll);if(!job)return;
 try{
  const r=await api('/api/export/'+job);
  if(r.state==='running'){
   $('progress').textContent=r.total!=null?'已检查 '+Number(r.processed).toLocaleString()+' / '+Number(r.total).toLocaleString()+' 个分子。':'正在建立只读数据快照…';
   poll=setTimeout(refresh,2000);return;
  }
  if(r.state==='error')throw Error(r.error);
  busy(false);clearSavedJob();
  $('progress').textContent=(r.preview?'预览完成：可保留 ':'导出完成：')+r.count.toLocaleString()+' 个分子。';
  $('downloads').innerHTML=r.preview?'':'<a href="/api/export/'+job+'/download">下载数据文件</a><a href="/api/export/'+job+'/manifest">下载筛选与校验清单</a>';
  showSourceCounts(r.source_split_counts);
  const s=r.selection;
  if(s){
   $('selection').hidden=false;
   $('counts').textContent='检查 '+s.scanned.toLocaleString()+' 个，筛选保留 '+s.kept.toLocaleString()+' 个，筛除 '+s.excluded.toLocaleString()+' 个。';
   $('reasonRows').innerHTML=Object.entries(s.reasons).map(([k,v])=>'<tr><td>'+esc(reasonNames[k]||k)+'</td><td>'+Number(v).toLocaleString()+'</td></tr>').join('')||'<tr><td colspan="2">未命中筛除条件。</td></tr>';
  }
  $('note').textContent=r.note||'谱类型、谱来源、归属和核种组合已写入导出清单。';
  $('examples').hidden=!r.examples?.length;
  $('exampleRows').innerHTML=(r.examples||[]).map(x=>'<p><a href="/?molecule_id='+encodeURIComponent(x.id)+
   '" target="_blank" rel="noopener">'+esc(x.formula)+' · '+esc(x.id)+'</a>　'+esc([coverageNames[x.coverage],nucleiNames[x.nuclei],
    ...(x.sources||[]).map(v=>sourceNames[v]||v)].filter(Boolean).join(' · '))+'</p>').join('');
 }catch(e){clearSavedJob();busy(false);notice(e.message);}
}
async function start(preview){
 const request=payload();
 notice('');busy(true);$('downloads').innerHTML='';$('selection').hidden=true;showSourceCounts(null);
 try{
  const r=await api('/api/export'+(preview?'/preview':''),request);job=r.job_id;
  try{localStorage.setItem('atomnmr.export.pending',job);}catch{}
  await refresh();
 }catch(e){busy(false);notice(e.message);}
}
$('preview').onclick=()=>start(true);$('run').onclick=()=>start(false);describe();
api('/api/export/sources').then(s=>{catalog=s;updateSources();describe();}).catch(e=>notice('来源列表暂时读取失败：'+e.message));
try{job=localStorage.getItem('atomnmr.export.pending');}catch{}if(job){busy(true);refresh();}
