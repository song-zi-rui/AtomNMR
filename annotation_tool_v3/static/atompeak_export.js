"use strict";
// An explicit new export profile; portable JSONL/Parquet/LMDB remain available.
$('#exportFormat').insertAdjacentHTML('afterbegin','<option value="atompeak">AtomPeak 从头训练包 · 单/双谱与部分归属</option>');
$('#exportScope').closest('label').insertAdjacentHTML('afterend',`
<section id="atompeakExportOptions" hidden>
 <p>生成可直接训练的 ZIP：train / valid / test、配套词表、逐峰监督掩码和导出统计。同一分子的各谱与修订版按连接层固定分组，未归属位置保留未知。</p>
 <label>划分随机种子<input id="exportTrainingSeed" type="number" value="42" step="1"></label>
 <label>验证集比例<input id="exportValidRatio" type="number" value="0.1" min="0" max="0.9" step="0.01"></label>
 <label>测试集比例<input id="exportTestRatio" type="number" value="0.1" min="0" max="0.9" step="0.01"></label>
 <label class="inline"><input id="exportUseSolvent" type="checkbox">输入各核种溶剂，缺失保留未知</label>
 <p>全库导出包含通过来源检查的未归属谱，可用于结构生成训练。只有确定标签参与硬归属监督；歧义和位移组另存。全库数据量较大，界面会显示处理进度。</p>
</section>`);
$('#exportFormat').onchange=()=>{
 const training=$('#exportFormat').value==='atompeak';$('#atompeakExportOptions').hidden=!training;
 if(training)$('#exportScope').value='all';
};
const portableExport=$('#exportRun').onclick;
$('#exportRun').onclick=async()=>{
 if($('#exportFormat').value!=='atompeak')return portableExport();
 const button=$('#exportRun');button.disabled=true;
 try{
  await save();const scope=$('#exportScope').value;
  if(scope==='current'&&!state.record)throw new Error('请先选择一个分子。');
  const payload={format:'atompeak',molecule_id:scope==='current'?state.record.id:null,filters:scope==='filtered'?filters():{},
   training_options:{seed:Number($('#exportTrainingSeed').value),valid_ratio:Number($('#exportValidRatio').value),
    test_ratio:Number($('#exportTestRatio').value),use_solvent:$('#exportUseSolvent').checked}};
  const job=await api('/api/export',payload);$('#exportResult').textContent='正在建立训练数据快照，可继续标注…';
  const poll=async()=>{
   try{
    const result=await api('/api/export/'+job.job_id);
    if(result.state==='running'){
     $('#exportResult').textContent=result.total!==undefined?`已处理 ${result.processed.toLocaleString()} / ${result.total.toLocaleString()} 个分子，已导出 ${result.exported.toLocaleString()} 个。`:'正在生成训练数据快照…';
     setTimeout(poll,1800);return;
    }
    if(result.state==='error')throw new Error(result.error);
    $('#exportResult').innerHTML=`已导出 ${result.count.toLocaleString()} 个分子。train ${result.split_counts.train.toLocaleString()} · valid ${result.split_counts.valid.toLocaleString()} · test ${result.split_counts.test.toLocaleString()}。<br><a href="/api/export/${job.job_id}/download">下载 AtomPeak 训练包</a><br><small>${escape(result.dataset_dir)}</small>`;
    button.disabled=false;
   }catch(error){$('#exportResult').textContent=error.message;button.disabled=false;}
  };poll();
 }catch(error){$('#exportResult').textContent=error.message;button.disabled=false;}
};
