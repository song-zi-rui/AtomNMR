const slots = ['group', 'first', 'second'];
const labels = ['组内 API', '官方第 1 次', '官方第 2 次'];
const settings = ['model', 'api_adapter', 'base_url', 'thinking', 'reasoning_effort',
  'max_tokens', 'timeout', 'wire_format', 'json_mode', 'stream', 'prompt_sha256'];
const statusLabels = {
  assigned: '确定归属', equivalent_group: '等价原子组', overlap: '重叠峰',
  ambiguous: '有歧义', unassigned: '未归属', solvent: '溶剂峰',
  impurity: '杂质峰', excluded: '排除'
};
const peakKey = (spectrum, peak) => JSON.stringify([spectrum, peak]);
const atomSet = atoms => [...new Set(atoms || [])].sort((a, b) => a - b);
const primaryKey = p => JSON.stringify(p ? [p.status, atomSet(p.atom_indices)] : null);
const candidateKey = p => JSON.stringify((p?.candidates || []).map(atomSet));

export function comparePredictions(first, second) {
  if (!first || !second) return 'missing';
  if (primaryKey(first) !== primaryKey(second)) return 'primary';
  return candidateKey(first) === candidateKey(second) ? 'same' : 'candidates';
}

export function alignRuns(runs) {
  const baseline = runs[0];
  if (!baseline.cases.length || !baseline.config.input_set_hash || !baseline.config.reference_set_hash ||
      runs.some(r => r.config.input_set_hash !== baseline.config.input_set_hash ||
        r.config.reference_set_hash !== baseline.config.reference_set_hash ||
        r.cases.length !== baseline.cases.length)) {
    throw Error('这三轮的冻结输入或参考答案不完全相同，请选择用同一份样本、同一分子数量运行的记录。');
  }
  return baseline.cases.map(c => runs.map(r => {
    const match = r.cases.find(other => other.molecule_id === c.molecule_id && other.input_hash === c.input_hash);
    if (!match) throw Error('无法对齐分子 ' + c.molecule_id);
    return match;
  }));
}

export function buildRows(caseGroups) {
  return caseGroups.flatMap(details => {
    const baseline = details[0];
    const isRef = r => r.source === 'shiftdb2';
    const byPeak = (items, filter = () => true) => new Map((items || []).filter(filter)
      .map(p => [peakKey(p.spectrum_id, p.peak_id), p]));
    const references = byPeak(baseline.reference, r => isRef(r) && r.kind === 'peak');
    const predictions = details.map(d => byPeak(d.output?.assignments));
    const scores = details.map(d => byPeak(d.comparisons, r => isRef(r) && r.eligible));
    return baseline.input.spectra.flatMap(s => s.peaks.map(peak => {
      const key = peakKey(s.spectrum_id, peak.peak_id);
      const answers = predictions.map(map => map.get(key));
      const marks = scores.map(map => map.get(key));
      const referenceGroups = baseline.reference.filter(r => isRef(r) && r.kind === 'group' &&
        r.spectrum_id === s.spectrum_id && r.peak_ids.includes(peak.peak_id));
      return {molecule: baseline.molecule_id, formula: baseline.formula, nucleus: s.nucleus,
        spectrum: s.spectrum_id, peak, reference: references.get(key), referenceGroups, answers, marks,
        scored: marks.every(Boolean), change: comparePredictions(answers[1], answers[2])};
    }));
  });
}

export function summarizeRows(rows) {
  const scored = rows.filter(r => r.scored);
  return {
    scored: scored.length,
    hits: slots.map((_, i) => scored.filter(r => r.marks[i].correct).length),
    top3: slots.map((_, i) => scored.filter(r => r.marks[i].top3).length),
    changes: Object.fromEntries(['same', 'primary', 'candidates', 'missing']
      .map(k => [k, scored.filter(r => r.change === k).length])),
  };
}

const $ = id => document.getElementById(id);
const elem = (tag, content, cls) => {
  const el = document.createElement(tag);
  if (content != null) el.textContent = content;
  if (cls) el.className = cls;
  return el;
};
const selected = () => slots.map(k => $(k).value);
const textSet = (items, nucleus) => atomSet(items).map(i => nucleus === '13C' ? 'C' + i : '原子 ' + i).join('、');
const hasPpm = val => typeof val === 'number' && Number.isFinite(val);
const percent = (n, total) => total ? (100 * n / total).toFixed(1) + '%' : '—';
const number = n => typeof n === 'number' ? n.toLocaleString('zh-CN') : '—';
let rows = [], timer = null;

async function api(path) {
  const response = await fetch('/api/benchmark' + path);
  if (!response.ok) throw Error('无法读取测试记录（HTTP ' + response.status + '）。');
  return response.json();
}

function answerCell(p, mark, nucleus, showReasons) {
  const td = elem('td');
  if (!p) {
    td.textContent = '未返回';
    return td;
  }
  td.append(elem('strong', textSet(p.atom_indices, nucleus) || statusLabels[p.status] || p.status));
  if (p.atom_indices?.length) td.append(elem('span', statusLabels[p.status] || p.status, 'repeat-secondary'));
  if (p.candidates?.length) {
    td.append(elem('span', '候选：' + p.candidates.map(c => '[' + textSet(c, nucleus) + ']').join(' → '), 'repeat-secondary'));
  }
  if (mark) {
    const msg = mark.correct ? '符合参考' : mark.attempted ? '与参考不符' : '未给确定归属';
    td.append(elem('div', msg + (!mark.correct && mark.top3 ? ' · 前三候选含参考' : ''), 'repeat-mark ' + (mark.correct ? 'good' : mark.attempted ? 'bad' : 'neutral')));
  }
  if (showReasons && p.brief_reason) td.append(elem('small', p.brief_reason, 'repeat-reason'));
  return td;
}

function renderRows() {
  const [mol, nu, filter] = ['molecule', 'nucleus', 'filter'].map(id => $(id).value);
  const visible = rows.filter(r => (!mol || r.molecule === mol) && (!nu || r.nucleus === nu) &&
    (filter === 'all' || filter === 'reference' && r.scored ||
      filter === 'changed' && ['primary', 'candidates'].includes(r.change) ||
      filter === 'primary' && r.change === 'primary'));
  $('rowCount').textContent = '显示 ' + visible.length + ' 个峰 · 其中 ' + visible.filter(r => r.scored).length + ' 个有可评分参考。¹H 的原子编号表示连接的重原子。';
  const fragment = document.createDocumentFragment();
  for (const r of visible) {
    const tr = elem('tr');
    tr.dataset.change = r.change;
    const peak = elem('th');
    peak.scope = 'row';
    peak.title = r.spectrum;
    peak.append(elem('div', r.formula), elem('div', r.nucleus + ' · ' + r.peak.peak_id +
      (hasPpm(r.peak.shift_ppm) ? ' · ' + r.peak.shift_ppm.toFixed(3) + ' ppm' : ''), 'repeat-shift'));
    tr.append(peak);
    const ref = answerCell(r.reference, null, r.nucleus, false);
    if (!r.reference) ref.textContent = '无逐峰参考';
    for (const group of r.referenceGroups) {
      ref.append(elem('span', '组参考：' + textSet(group.atom_indices, r.nucleus) + '（' + group.peak_ids.join('、') + '）', 'repeat-secondary'));
    }
    tr.append(ref, ...r.answers.map((p, i) => answerCell(p, r.marks[i], r.nucleus, $('reasons').checked)));
    const changes = {same:'一致', primary:'主归属或状态变化', candidates:'仅候选变化', missing:'结果不全'};
    tr.append(elem('td', changes[r.change], 'repeat-change'));
    fragment.append(tr);
  }
  if (!visible.length) {
    const empty = elem('td', '没有符合当前筛选的峰。');
    empty.colSpan = 6;
    const tr = elem('tr');
    tr.append(empty);
    fragment.append(tr);
  }
  $('peakRows').replaceChildren(fragment);
}

function card(title, value, detail) {
  const div = elem('div', null, 'repeat-metric');
  div.append(elem('span', title), elem('strong', value), elem('small', detail));
  return div;
}

function renderSummary(runs, caseGroups) {
  const first = runs[1].config, second = runs[2].config;
  const differing = settings.filter(k => first.primary.config[k] !== second.primary.config[k]);
  for (const k of ['parallelism', 'strategy', 'routing_enabled']) {
    if (first[k] !== second[k]) differing.push(k);
  }
  const wires = caseGroups.map(cs => cs.slice(1).map(c => c.attempts[0]?.diagnostics?.wire_sha256));
  const identicalWires = wires.filter(([a, b]) => a && a === b).length;
  const comparableWires = wires.filter(([a, b]) => a && b).length;
  $('controlsCheck').textContent = '冻结输入与 ShiftDB2 参考相同。' +
    (differing.length ? ' 两次官方的已保存设置有差异：' + differing.join('、') + '。' : ' 两次官方的已保存模型和请求设置一致。') +
    ' 谱数据请求指纹相同：' + identicalWires + '/' + comparableWires + ' 个可核对分子。';
  const s = summarizeRows(rows);
  const n = s.scored;
  $('metrics').replaceChildren(...runs.map((r, i) => {
    const time = r.cases.every(c => typeof c.elapsed === 'number') ? r.cases.reduce((a, c) => a + c.elapsed, 0).toFixed(1) + ' 秒' : '用时未齐';
    const usage = r.usage_reported_attempts === r.api_attempts && r.api_attempts ? number(r.usage.total_tokens) : '—';
    return card(labels[i], s.hits[i] + '/' + n + ' · ' + percent(s.hits[i], n),
      r.config.primary.config.model + ' · Top 3 ' + s.top3[i] + '/' + n + ' · ' + time + ' · ' + usage + ' tokens');
  }), card('两次官方一致', s.changes.same + '/' + (n - s.changes.missing),
    '主归属或状态变 ' + s.changes.primary + ' · 仅候选变 ' + s.changes.candidates +
    (s.changes.missing ? ' · 无法比较 ' + s.changes.missing : '')));
}

async function loadRuns() {
  clearTimeout(timer);
  const ids = selected();
  if (ids.some(id => !id) || new Set(ids).size !== 3) throw Error('请选择三个不同的测试记录。');
  $('loadRuns').disabled = true;
  $('results').hidden = true;
  $('loadStatus').textContent = '正在读取三轮结果…';
  try {
    const runs = await Promise.all(ids.map(id => api('/runs/' + encodeURIComponent(id))));
    const aligned = alignRuns(runs);
    const url = new URL(location.href);
    slots.forEach((key, i) => url.searchParams.set(key, ids[i]));
    history.replaceState(null, '', url);
    if (runs.some(r => r.state === 'running')) {
      $('loadStatus').textContent = '测试仍在运行：' + runs.map((r, i) => labels[i] + ' ' + (r.counts.done || 0) + '/' + r.cases.length).join(' · ') + '。完成后自动显示结果。';
      timer = setTimeout(() => {
        if (selected().every((id, i) => id === ids[i])) loadRuns().catch(showError);
      }, 12000);
      return;
    }
    const details = await Promise.all(aligned.map(cases => Promise.all(cases.map(c => api('/cases/' + c.id)))));
    rows = buildRows(details);
    const options = details.map(([d]) => new Option(d.formula + ' · ' + d.molecule_id, d.molecule_id));
    $('molecule').replaceChildren(new Option('全部分子', ''), ...options);
    renderSummary(runs, details);
    renderRows();
    $('results').hidden = false;
    $('loadStatus').textContent = runs.every(r => r.state === 'completed')
      ? '三轮测试均已完成。' : '有轮次尚未完成；没有返回的归属会标记为“未返回”。';
  } finally {
    $('loadRuns').disabled = false;
  }
}

function showError(error) {
  $('loadStatus').textContent = error.message;
}

async function start() {
  $('runForm').addEventListener('submit', event => {
    event.preventDefault();
    loadRuns().catch(showError);
  });
  for (const id of ['molecule', 'nucleus', 'filter', 'reasons']) $(id).addEventListener('change', renderRows);
  const {runs} = await api('/runs');
  const requested = new URLSearchParams(location.search);
  for (const id of slots) {
    const select = $(id);
    select.replaceChildren(new Option('请选择一轮', ''),
      ...runs.map(r => new Option(r.name + ' · ' + String(r.created_at).slice(0, 10), r.id)));
    const wanted = requested.get(id);
    if (wanted && !runs.some(r => r.id === wanted)) select.add(new Option('链接指定的轮次 · ' + wanted, wanted));
    select.value = wanted || '';
  }
  if (selected().every(Boolean)) await loadRuns();
  else $('loadStatus').textContent = '请选择三轮共用相同冻结样本的记录。';
}

if (typeof document !== 'undefined') start().catch(showError);
