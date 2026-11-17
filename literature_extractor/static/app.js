const element = id => document.getElementById(id);
const state = {table: null, nucleus: '1H', row: null, atom: null};
const relation = {interchangeable: '可互换', overlapping: '重叠区间', shared_shift: '共享位移', unassigned: '未给归属'};
const integrated = location.pathname.startsWith('/literature/');
const jobStorage = integrated ? 'atomnmr-literature-job' : 'literature-job';
const importLabels = {added: '新增实验谱', filled: '补入归属', duplicate: '重复跳过', review: '待复核'};
if (integrated) {
  element('main-nav').hidden = false;
  element('intro').textContent = '提取论文的 ¹H / ¹³C 位移和归属，与当前数据库比对后自动入库。';
}
let activeJob = null;
let poll = null;

function text(tag, value, css) {
  const node = document.createElement(tag);
  node.textContent = value;
  if (css) node.className = css;
  return node;
}

function link(url, title) {
  const node = text('a', title);
  node.href = url;
  node.target = '_blank';
  node.rel = 'noopener';
  return node;
}

function notice(value, error = false) {
  element('notice').textContent = value;
  element('notice').classList.toggle('error', error);
}

async function request(url, payload, signal) {
  const options = payload === undefined ? {} : {
    method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(payload)
  };
  const response = await fetch(integrated ? '/literature' + url : url, {...options, signal});
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || `请求失败 (${response.status})`);
  return data;
}

function busy(active) {
  for (const id of ['folder','start','choose','limit','latest','group-preview','group-apply']) element(id).disabled = active;
  element('pause').hidden = !active || !activeJob;
}

async function track(id) {
  clearTimeout(poll);
  activeJob = id;
  try {
    const job = await request(`/api/jobs/${id}`);
    if (activeJob !== id) return;
    element('folder').value = job.folder;
    const running = job.state === 'running';
    busy(running);
    element('progress').hidden = false;
    element('resume').hidden = !['paused','error'].includes(job.state);
    element('pause').disabled = !!job.pause;
    element('pause').textContent = job.pause ? '正在暂停…' : '暂停';
    const bar = element('progress').querySelector('progress');
    const tally = [`${job.saved ?? job.result?.saved ?? 0}${job.limit ? ` / ${job.limit}` : ''} 个分子`,
      job.skipped ? `跳过 ${job.skipped} 个重复或已处理目录` : '',
      job.duplicate_spectra ? `跳过 ${job.duplicate_spectra} 条重复谱` : ''].filter(Boolean).join(' · ');
    if (job.total) {bar.max = job.total; bar.value = running ? Math.max(0, (job.index || 1) - 1) : (job.index || (job.state === 'complete' && !job.limit ? job.total : 0));}
    else bar.removeAttribute('value');
    element('progress').querySelector('span').textContent = [tally, running ? job.filename || job.current || job.reason : ''].filter(Boolean).join(' · ');
    if (running) {
      notice(job.pause ? '当前文件夹解析结束后暂停。' : '');
      poll = setTimeout(() => track(id), 700);
    } else {
      if (document.activeElement !== element('limit')) element('limit').value = job.limit ?? '';
      if (job.result) showResults(job.result);
      notice([job.reason, ...(job.errors || [])].filter(Boolean).join(' · '), job.state === 'error');
    }
  } catch (error) {
    notice(error.message, true);
    busy(false);
    element('resume').hidden = true;
    sessionStorage.removeItem(jobStorage);
  }
}

function showResults(result) {
  element('results').hidden = false;
  const review = result.tables.filter(table => result.integrated
    ? !Object.keys(table.main_import || {}).length || Object.values(table.main_import).some(item => item.status === 'review')
    : !table.molecule_id).length;
  const tableUnit = result.tables.some(table => table.comparison) ? '组归属表' : '张归属表';
  element('counts').textContent = `${result.pdf_count} 份 PDF · ${result.tables.length} ${tableUnit} · ${result.peak_count} 行峰表 · ${result.saved} 个分子${result.integrated ? '已入库或补归属' : '已保存'}${review ? ` · ${review} 张表待复核` : ''}`;
  const totals = {};
  for (const table of result.tables) for (const item of Object.values(table.main_import || {})) {
    totals[item.status] = (totals[item.status] || 0) + 1;
  }
  element('destination').textContent = result.integrated
    ? Object.entries(importLabels).map(([key, label]) => `${label} ${totals[key] || 0} 条`).join(' · ')
    : result.database || result.output_folder;
  element('database').hidden = !result.database_url;
  if (result.database_url) element('database').href = result.database_url;
  element('report').href = result.report_url;
  const records = result.files.map(file => `${file.name} · ${file.status}${file.issue ? ` · ${file.issue}` : ''}`);
  element('file-list').replaceChildren(...records.concat(result.warnings).map(value => text('li', value)));
  element('files').open = !result.tables.length || result.files.some(file => file.issue);
  notice(result.integrated ? '主库比对结果已保存。'
    : result.database ? '提取完成，数据库已保存。' : '提取完成，报告已保存。请查看文件处理记录和待复核原因。');
  element('molecules').replaceChildren();
  for (const table of result.tables) {
    const button = text('button', table.name);
    button.type = 'button';
    button.addEventListener('click', () => showTable(table, button));
    element('molecules').append(button);
  }
  element('selected').hidden = !result.tables.length;
  element('molecules').firstElementChild?.click();
}

function showTable(table, button) {
  state.table = table;
  state.nucleus = table.peaks['1H'].length ? '1H' : '13C';
  state.row = null;
  state.atom = null;
  element('filter').value = '';
  [...element('molecules').children].forEach(item => item.setAttribute('aria-pressed', item === button));
  element('chemical').textContent = table.name;
  const imports = Object.entries(table.main_import || {});
  const status = table.main_import ? imports.map(([n, item]) => `${n} ${importLabels[item.status]}`).join(' · ') || '待复核'
    : table.molecule_id ? '已保存' : '待复核';
  const reference = [table.formula, `PDF 第 ${table.page} 页${table.table ? ` · 表 ${table.table}` : ''}`, table.conditions?.solvent, status];
  element('reference').textContent = reference.filter(Boolean).join(' · ');
  const sources = [link(table.pdf_url, '打开原文归属表 ↗')];
  const mainId = imports.map(([, item]) => item.molecule_id).find(id => /^M[a-f0-9]{24}$/.test(id || ''));
  if (mainId) sources.push(link('/?molecule_id=' + mainId, '打开数据库记录 ↗'));
  if (table.supplement_url) sources.push(link(table.supplement_url, '打开补充谱图 ↗'));
  element('sources').replaceChildren(...sources);
  const issues = [...table.issues, ...imports.filter(([, item]) => item.status === 'review').map(([n, item]) => n + ' · ' + item.reason)];
  element('problems').replaceChildren(...issues.map(reason => text('li', reason)));
  element('problems').hidden = !issues.length;
  const missing = table.unreported_proton_sites || [];
  element('missing').textContent = `原文未报告位移 ${missing.join('、')}`;
  element('missing').hidden = !missing.length;
  const tabs = Object.entries(table.peaks).map(([key, rows]) => {
    const tab = text('button', `${key === '1H' ? '¹H' : '¹³C'} · ${rows.length}`);
    tab.type = 'button';
    tab.dataset.nucleus = key;
    tab.onclick = () => {state.nucleus = key; state.row = null; renderPeaks();};
    return tab;
  });
  element('nuclei').replaceChildren(...tabs);
  const structure = element('structure');
  structure.replaceChildren();
  if (table.diagram) {
    const canvas = text('div', '', 'drawing');
    canvas.innerHTML = table.diagram; // Markup comes from the RDKit drawing, never from a PDF.
    const overlay = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    overlay.setAttribute('viewBox', '0 0 580 380');
    overlay.classList.add('hotspots');
    overlay.setAttribute('aria-label', '点击原子筛选峰表');
    for (const [number, [x, y]] of Object.entries(table.points)) {
      const circle = document.createElementNS(overlay.namespaceURI, 'circle');
      circle.setAttribute('cx', x); circle.setAttribute('cy', y); circle.setAttribute('r', '11');
      circle.setAttribute('tabindex', '0'); circle.setAttribute('role', 'button');
      circle.setAttribute('aria-label', `原子 ${table.atom_display?.[number] || number}`);
      circle.dataset.atom = number;
      circle.onclick = () => {
        state.atom = state.atom === number ? null : number;
        state.row = null;
        element('filter').value = '';
        renderPeaks();
      };
      circle.onkeydown = event => {
        if (event.key === ' ' || event.key === 'Enter') {event.preventDefault(); circle.onclick();}
      };
      overlay.append(circle);
    }
    canvas.append(overlay);
    structure.append(canvas);
  } else {
    structure.append(text('p', '结构待核对，请打开原文查看编号图。'));
  }
  renderPeaks();
}

function highlight(peak, rows) {
  const overlay = element('structure').querySelector('.hotspots');
  if (!overlay) return;
  const selected = !peak ? [] : peak.group ? rows.filter(row => row.group === peak.group) : [peak];
  const atomNumbers = row => row.paper_atom_numbers?.map(String) || row.source_labels.map(label => label.match(/\d+/)[0]);
  const numbers = new Set(selected.flatMap(atomNumbers));
  for (const circle of overlay.children) {
    const active = numbers.has(circle.dataset.atom);
    circle.setAttribute('aria-pressed', circle.dataset.atom === state.atom);
    circle.classList.toggle('highlight', active);
    circle.classList.toggle('candidate', active && !!peak?.group);
  }
}

function renderPeaks() {
  [...element('nuclei').children].forEach(tab => tab.setAttribute('aria-pressed', tab.dataset.nucleus === state.nucleus));
  const rows = state.table.peaks[state.nucleus];
  const proton = state.nucleus === '1H';
  const comparison = state.table.comparison === 'synthetic';
  if (comparison) {
    const source = element('sources').firstElementChild;
    const page = state.table.component_pages?.[state.nucleus]?.[0];
    if (source && page) source.href = `${state.table.pdf_url.split('#')[0]}#page=${page}`;
  }
  const columns = [comparison ? '合成 δ / ppm' : 'δ / ppm', '论文归属'];
  if (proton) columns.push('峰型', '氢数', 'J / Hz');
  if (comparison) columns.push('天然样品原文');
  const heading = document.createElement('tr');
  columns.forEach(column => heading.append(text('th', column)));
  element('peak-head').replaceChildren(heading);
  const normalize = value => value.toLowerCase().replace(/[’']/g, '′').replace(/\s+/g, '');
  const query = normalize(element('filter').value);
  const filtered = rows.filter(row => {
    const members = row.group ? rows.filter(peak => peak.group === row.group) : [row];
    const labels = members.flatMap(peak => peak.source_labels);
    const numbers = members.flatMap(peak => peak.paper_atom_numbers?.map(String) || []);
    return (!state.atom || numbers.includes(state.atom)) &&
      [row.shift_text, row.source_assignment_text, ...labels].some(value => normalize(value).includes(query));
  });
  const content = filtered.map(row => {
    const line = document.createElement('tr');
    if (row === state.row) line.classList.add('current');
    const assignment = document.createElement('td');
    const pick = text('button', row.source_assignment_text || '—', 'pick');
    pick.type = 'button';
    pick.addEventListener('click', () => {
      const selected = state.row !== row;
      state.row = selected ? row : null;
      element('peaks').querySelector('.current')?.classList.remove('current');
      line.classList.toggle('current', selected);
      highlight(state.row, rows);
    });
    assignment.append(pick);
    if (relation[row.relation]) {
      const badge = text('span', row.group_reason ? '编号映射待定' : relation[row.relation], `relation${row.group ? ' swap' : ''}`);
      if (row.group_reason) badge.title = row.group_reason;
      assignment.append(badge);
    }
    line.append(text('td', row.shift_text), assignment);
    if (proton) line.append(text('td', row.multiplicity ?? '—'), text('td', row.reported_nH ?? '—'), text('td', row.j_values_hz.join(', ') || '—'));
    if (comparison) {
      const reference = text('td', row.natural_reference_cell || '—');
      if (row.source_page) reference.title = `补充材料第 ${row.source_page} 页`;
      line.append(reference);
    }
    return line;
  });
  if (!content.length) {
    const line = document.createElement('tr');
    const empty = text('td', '没有符合条件的峰。');
    empty.colSpan = columns.length;
    line.append(empty);
    content.push(line);
  }
  element('peaks').replaceChildren(...content);
  highlight(filtered.includes(state.row) ? state.row : null, rows);
}

element('filter').oninput = () => {state.atom = null; renderPeaks();};
const picker = element('folder-picker');
let folderView = null;
let folderRequest = null;

async function browseFolder(folder) {
  folderRequest?.abort();
  const controller = new AbortController();
  folderRequest = controller;
  const timeout = setTimeout(() => controller.abort('timeout'), 9000);
  const status = element('folder-status');
  const list = element('folder-list');
  element('folder-address').value = folder;
  element('folder-select').disabled = true;
  element('folder-up').disabled = true;
  folderView = null;
  list.replaceChildren();
  status.textContent = '正在读取文件夹…';
  status.classList.remove('error');
  try {
    const view = await request('/api/folders', {folder}, controller.signal);
    if (folderRequest !== controller) return;
    folderView = view;
    element('folder-address').value = view.folder;
    element('folder-up').disabled = !view.parent;
    element('folder-select').disabled = false;
    const buttons = view.children.map(child => {
      const button = text('button', child.name, 'folder-entry');
      button.type = 'button';
      button.setAttribute('aria-label', child.name);
      button.onclick = () => browseFolder(child.path);
      return button;
    });
    list.replaceChildren(...buttons);
    status.textContent = buttons.length ? '' : '此文件夹没有子文件夹，可以直接使用。';
  } catch (error) {
    if (folderRequest !== controller || !picker.open) return;
    status.textContent = controller.signal.reason === 'timeout'
      ? '读取超时。可以输入其他路径后重试。' : error.message;
    status.classList.add('error');
  } finally {
    clearTimeout(timeout);
  }
}

element('choose').onclick = () => {
  picker.showModal();
  browseFolder(element('folder').value.trim());
};
element('folder-navigation').onsubmit = event => {
  event.preventDefault();
  browseFolder(element('folder-address').value.trim());
};
element('folder-address').oninput = () => {
  element('folder-select').disabled = !folderView || element('folder-address').value !== folderView.folder;
};
element('folder-up').onclick = () => {
  if (folderView?.parent) browseFolder(folderView.parent);
};
element('folder-cancel').onclick = () => picker.close();
element('folder-select').onclick = () => {
  if (!folderView) return;
  element('folder').value = folderView.folder;
  resetRoot();
  notice('');
  picker.close();
};
picker.onclose = () => {
  if (picker.open) return;
  folderRequest?.abort();
  folderRequest = null;
};
function resetRoot() {
  activeJob = null;
  clearTimeout(poll);
  sessionStorage.removeItem(jobStorage);
  history.replaceState(null, '', location.pathname);
  element('progress').hidden = true;
  element('results').hidden = true;
  element('resume').hidden = true;
  element('group-apply').hidden = true;
  element('group-list').replaceChildren();
  element('group-status').textContent = '';
}

function remember(id) {
  sessionStorage.setItem(jobStorage, id);
  history.replaceState(null, '', `?job=${id}`);
  activeJob = id;
}

async function launch(resume) {
  const folder = element('folder').value.trim();
  if (!element('extract').reportValidity()) return;
  element('results').hidden = true;
  element('resume').hidden = true;
  notice('');
  busy(true);
  element('progress').hidden = false;
  element('progress').querySelector('span').textContent = '正在查找待处理文件夹';
  try {
    const limit = element('limit').value ? Number(element('limit').value) : null;
    const {id} = await request('/api/extract', {folder, limit, resume});
    remember(id);
    await track(id);
  } catch (error) {notice(error.message, true); busy(false); element('resume').hidden = !resume;}
}
element('folder').oninput = resetRoot;
element('extract').onsubmit = event => {event.preventDefault(); launch(null);};
element('resume').onclick = () => launch(activeJob);
element('pause').onclick = async () => {
  element('pause').disabled = true;
  try { await request('/api/pause', {id:activeJob}); }
  catch (error) { notice(error.message, true); element('pause').disabled = false; }
};
element('latest').onclick = async () => {
  try {
    const {id} = await request('/api/latest', {folder:element('folder').value.trim()});
    if (!id) return notice('该目录暂无批量提取任务。');
    remember(id);
    element('results').hidden = true;
    track(id);
  } catch (error) {notice(error.message, true);}
};

async function group(apply) {
  const folder = element('folder').value.trim();
  if (!element('folder').reportValidity()) return;
  const status = element('group-status');
  element('group-preview').disabled = true;
  element('group-apply').hidden = true;
  status.textContent = '正在核对 PDF 首页和文件名…';
  try {
    const {groups, pending, moved, duplicates} = await request('/api/group', {folder, apply});
    if (folder !== element('folder').value.trim()) return;
    element('group-list').replaceChildren(...[
      ...groups.map(row => `${row.status} · ${row.folder} · ${row.evidence.join('、')} · ${row.files.join(' + ')}`),
      ...pending.map(row => `待手动分组 · ${row.file} · ${row.reason}`)
    ].map(x => text('li', x)));
    status.textContent = [apply ? `移动 ${moved} 份 PDF` : `${groups.length} 篇论文可整理`,
      apply && duplicates ? `清理 ${duplicates} 份内容相同的散放文件` : '',
      `${pending.length} 份 PDF 待确认`].filter(Boolean).join(' · ');
    element('group-apply').hidden = !groups.some(row => row.status === '可分组');
  } catch (error) {status.textContent = error.message;}
  finally {element('group-preview').disabled = !element('pause').hidden;}
}
element('group-preview').onclick = () => group(false);
element('group-apply').onclick = () => group(true);

let found = [];
element('search-papers').onsubmit = async event => {
  event.preventDefault();
  const categories = [...event.currentTarget.querySelectorAll('input[type=checkbox]:checked')].map(input => input.value);
  const status = element('search-status');
  if (!categories.length) {status.textContent = '请选择至少一个检索类别。'; return;}
  element('search-start').disabled = true;
  element('search-export').hidden = true;
  element('search-results').hidden = true;
  status.textContent = '正在检索 Crossref 和 Europe PMC…';
  try {
    const {results, warnings} = await request('/api/search', {name:element('compound-name').value,
      limit:Number(element('search-limit').value), categories});
    found = results;
    element('paper-rows').replaceChildren(...results.map(paper => {
      const row = document.createElement('tr');
      const title = document.createElement('td');
      title.append(link(paper.url, paper.title));
      title.append(text('div', [paper.journal,paper.year].filter(Boolean).join(' · '), 'paper-meta'));
      if (paper.open_url) title.append(link(paper.open_url, '开放全文 ↗'));
      row.append(title, text('td', paper.categories.join('、')), text('td', paper.evidence));
      return row;
    }));
    element('search-results').hidden = !results.length;
    element('search-export').hidden = !results.length;
    status.textContent = `${results.length} 篇不同论文。打开正文和补充材料后可确认是否有归属表。${warnings.length ? '部分检索请求未成功。' : ''}`;
    status.title = warnings.join('\n');
  } catch (error) {status.textContent = error.message;}
  finally {element('search-start').disabled = false;}
};
element('search-export').onclick = () => {
  const fields = ['title','url','open_url','doi','journal','year','categories','evidence','provider'];
  const cell = v => `"${String(v ?? '').replace(/^[=+@-]/, m => `'${m}`).replaceAll('"','""')}"`;
  const lines = [fields.join(','), ...found.map(row => fields.map(k => cell(k === 'categories' ? row[k].join(' / ') : row[k])).join(','))];
  const href = URL.createObjectURL(new Blob(['\ufeff'+lines.join('\r\n')], {type:'text/csv;charset=utf-8'}));
  const file = document.createElement('a');
  file.href = href; file.download = 'nmr_literature_search.csv'; file.click();
  setTimeout(() => URL.revokeObjectURL(href), 500);
};
const lastJob = new URLSearchParams(location.search).get('job') || sessionStorage.getItem(jobStorage);
if (/^[a-f0-9]{12}$/.test(lastJob || '')) {
  busy(true);
  track(lastJob);
}
