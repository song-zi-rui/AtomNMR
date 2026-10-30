"use strict";
$('#structureEnlarge').insertAdjacentHTML('beforebegin','<button id="copyMolecule" title="复制当前显示的分子图，包含原子编号和高亮，可直接粘贴图片。">复制分子图</button>');

async function moleculePng(src){
  const response=await fetch(src);
  if(!response.ok)throw new Error('分子图读取失败。');
  const url=URL.createObjectURL(await response.blob());
  try{
    const img=new Image();
    await new Promise((resolve,reject)=>{img.onload=resolve;img.onerror=()=>reject(new Error('分子图加载失败。'));img.src=url;});
    const canvas=document.createElement('canvas');
    canvas.width=Math.min(3600,Math.max(1200,img.naturalWidth*2));
    canvas.height=Math.max(1,Math.round(canvas.width*img.naturalHeight/img.naturalWidth));
    const context=canvas.getContext('2d');
    if(!context)throw new Error('浏览器无法生成图片。');
    context.fillStyle='#fff';context.fillRect(0,0,canvas.width,canvas.height);
    context.drawImage(img,0,0,canvas.width,canvas.height);
    return await new Promise((resolve,reject)=>canvas.toBlob(blob=>blob?resolve(blob):reject(new Error('PNG 图片生成失败。')),'image/png'));
  }finally{URL.revokeObjectURL(url);}
}

async function copyMolecule(){
  if(!state.record||!$('#molecule').src)return toast('请先打开一个分子。');
  if(!navigator.clipboard?.write||typeof ClipboardItem==='undefined')return toast('当前浏览器不支持图片复制，请使用 Edge / Chrome。');
  const button=$('#copyMolecule');button.disabled=true;
  try{
    // Start the clipboard operation during the click's user activation.
    const png=moleculePng($('#molecule').src);png.catch(()=>{});
    await navigator.clipboard.write([new ClipboardItem({'image/png':png})]);
    toast('已复制分子图，可直接粘贴到文档或聊天中。');
  }catch(e){toast('分子图复制失败：'+e.message);}
  finally{button.disabled=false;}
}
$('#copyMolecule').onclick=copyMolecule;
