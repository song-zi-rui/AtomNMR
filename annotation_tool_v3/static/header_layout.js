"use strict";
// Keep sticky panels below the real header, including wrapping or font zoom.
(() => {
  const header=document.querySelector('.app-header');
  const update=()=>document.documentElement.style.setProperty('--app-header-height',Math.ceil(header.getBoundingClientRect().height)+'px');
  update();
  if(window.ResizeObserver)new ResizeObserver(update).observe(header);
  else window.addEventListener('resize',update);
  const opened=()=>[...header.querySelectorAll('.header-disclosure[open]')];
  document.addEventListener('click',event=>{
    for(const item of opened())if(!item.contains(event.target))item.open=false;
  });
  document.addEventListener('keydown',event=>{
    if(event.key!=='Escape'||!opened().length)return;
    const items=opened();items.forEach(item=>item.open=false);
    items[0].querySelector('summary').focus();
    event.preventDefault();event.stopImmediatePropagation();
  },true);
})();
