"""Read-only, snapshot-consistent census of molecule assignment coverage.

Run from the repository root with ``python -m tools.summarize_annotation_database``.
The application database is opened mode=ro and is never initialized or refreshed.
"""
from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import time
from collections import Counter, defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

from annotation_tool_v3 import core as db


def log(message):
    print(message, flush=True)


def coverage(total, assigned):
    return 'none' if assigned == 0 else 'full' if total == assigned else 'partial'


def kind(cn, hn):
    return 'dual' if cn and hn else 'carbon_only' if cn else 'proton_only' if hn else 'no_current_spectrum'


def row_key(row):
    return (row['origin'], kind(row['c_count'], row['h_count']),
            coverage(row['c_count'] + row['h_count'], row['c_assigned'] + row['h_assigned']),
            int(row['task_done']))


def summarize(rows):
    counts = Counter()
    for row in rows:
        counts['molecules'] += 1
        counts[row['coverage']] += 1
        counts[row['kind']] += 1
        counts['finished'] += row['task_done'] == 1
        counts['manual_annotation_history'] += row['manual_annotation_history']
        counts['human_participation'] += row['human_participation']
    return {key: counts[key] for key in ('molecules', 'full', 'partial', 'none', 'dual',
            'carbon_only', 'proton_only', 'no_current_spectrum', 'finished',
            'manual_annotation_history', 'human_participation')}


def write_markdown(report, output):
    kinds = ('dual', 'carbon_only', 'proton_only', 'no_current_spectrum')
    names = {'dual': 'H/C 双谱', 'carbon_only': '只有 C 谱', 'proton_only': '只有 H 谱',
             'no_current_spectrum': '当前无可用谱（隔离 / 停用）'}
    cross = defaultdict(Counter)
    processing = defaultdict(Counter)
    for row in report['cross_tabulation']:
        cross[row['kind']][row['coverage']] += row['molecules']
        processing[row['task_done']][row['coverage']] += row['molecules']
    scopes = report['scope_summaries']
    manual = scopes['manual_annotation_history']
    human = scopes['human_participation']
    native = report['native_original']
    fmt = lambda n: f'{n:,}'
    lines = [
        '# 当前实验 NMR 数据库归属统计', '',
        f"统计快照：{report['snapshot_at']}（北京时间）。",
        f"读取文件：`{report['database']}`。所有查询在同一个只读事务中执行。", '',
        '## 计数口径', '',
        '- 按数据库分子 ID 去重；同分子的多条实验谱、原始版、修订版只计一个分子。结构身份依据存储的规范异构 SMILES，不按分子式合并。',
        '- “全归属”指当前纳入峰表的所有峰均有确定归属状态，包含唯一归属、等价组、重叠和原库多原子组。它没有核实原始峰表是否漏峰，也没有证明标签化学正确。',
        '- “部分归属”指至少一个峰已有上述归属，但仍有峰处于未处理、歧义、位移组顺序不明或其他非确定状态。“完全未归属”指当前没有这样的确定峰归属。',
        '- 当前谱采用标注界面口径：排除自动隔离和人工停用、参考及历史谱。整分子的排除训练状态另行保留，统计不等于训练集准入审查。',
        '- “已结束”是用户提交的任务状态，允许保留部分归属，因此单独统计。', '',
        '## 1. 全库归属与单 / 双谱', '',
        f"存档 **{fmt(report['molecules_total'])}** 个分子；其中 **{fmt(report['current_coverage']['full']+report['current_coverage']['partial'])}** 个至少有一部分归属（含原库自带标签）。", '',
        '| 当前谱类型 | 全归属 | 部分归属 | 完全未归属 | 合计 |',
        '| --- | ---: | ---: | ---: | ---: |',
    ]
    for key in kinds:
        values = cross[key]
        lines.append(f"| {names[key]} | {fmt(values['full'])} | {fmt(values['partial'])} | {fmt(values['none'])} | {fmt(sum(values.values()))} |")
    lines += [f"| 合计 | {fmt(report['current_coverage']['full'])} | {fmt(report['current_coverage']['partial'])} | {fmt(report['current_coverage']['none'])} | {fmt(report['molecules_total'])} |", '',
        '全库有标签分子的数量不能当作本研究人工标注的数量。', '',
        '## 2. 后续人工参与', '',
        f"有过人工峰标注记录的分子 **{fmt(manual['molecules'])}** 个。这里包括人工填写、歧义 / 清空处理以及归属修订历史。另有 **{fmt(scopes['human_management_or_save_only']['molecules'])}** 个仅保存备注或进行了管理操作；广义人工参与合计 **{fmt(human['molecules'])}** 个。", '',
        '| 有人工峰标注历史的分子来源 | 全归属 | 部分归属 | 当前无确定归属 | 合计 |',
        '| --- | ---: | ---: | ---: | ---: |',
    ]
    for key, label in (('manual_nmrshiftdb2', 'ShiftDB2 来源'), ('manual_NMRexp', 'NMRexp 来源'), ('manual_annotation_history', '合计')):
        values = scopes[key]
        lines.append(f"| {label} | {fmt(values['full'])} | {fmt(values['partial'])} | {fmt(values['none'])} | {fmt(values['molecules'])} |")
    lines += ['',
        f"人工标注历史分子中，当前仍有双谱的 **{fmt(manual['dual'])}** 个，只有 C 谱 **{fmt(manual['carbon_only'])}** 个，只有 H 谱 **{fmt(manual['proton_only'])}** 个，当前无可用谱 **{fmt(manual['no_current_spectrum'])}** 个。当前仍有人工峰处理记录的分子为 **{fmt(report['manual_active_molecules'])}** 个；273 的历史口径包含后来停用、撤下等记录。", '',
        f"当前标为“已结束”的任务共有 **{fmt(manual['finished'])}** 个，归属情况如下：", '',
        '| 任务状态 | 全归属 | 部分归属 | 完全未归属 | 合计 |',
        '| --- | ---: | ---: | ---: | ---: |',
    ]
    for done, label in ((1, '已结束'), (0, '尚未结束')):
        values = processing[done]
        lines.append(f"| {label} | {fmt(values['full'])} | {fmt(values['partial'])} | {fmt(values['none'])} | {fmt(sum(values.values()))} |")
    lines += ['', '## 3. ShiftDB2 原始双谱', '',
        f"主库包含 **{fmt(native['structures'])}** 个 ShiftDB2 来源结构，保留 **{fmt(native['spectra'])}** 条原始 H/C 实验谱。", '',
        f"在原始 ShiftDB2 谱中同时有非空 C 谱和 H 谱的去重结构有 **{fmt(native['dual_structures'])}** 个。其中 **{fmt(native['same_record_dual_structures'])}** 个的两类谱可追到同一条源记录；其余 **{fmt(native['dual_structures']-native['same_record_dual_structures'])}** 个是在相同结构下汇集了不同源记录的 H/C 谱。", '',
        '| 原始双谱分子后续是否有人工参与 | 全归属 | 部分归属 | 当前无确定归属 | 合计 |',
        '| --- | ---: | ---: | ---: | ---: |',
    ]
    for key, label in (('native_dual_without_human_participation', '未记录后续人工参与'),
                       ('native_dual_with_human_participation', '记录过后续人工参与'),
                       ('original_native_dual', '合计')):
        values = scopes[key]
        lines.append(f"| {label} | {fmt(values['full'])} | {fmt(values['partial'])} | {fmt(values['none'])} | {fmt(values['molecules'])} |")
    lines += ['',
        f"原始双谱集合与广义人工参与集合重合 **{fmt(scopes['native_dual_with_human_participation']['molecules'])}** 个，不能将两个集合的总数直接相加。原始双谱中有 **{fmt(scopes['original_native_dual']['manual_annotation_history'])}** 个留有实际人工峰标注历史。", '',
        f"上述原始双谱分子目前仍有当前 H/C 双谱的有 **{fmt(scopes['original_native_dual']['dual'])}** 个；另 **{fmt(scopes['original_native_dual']['carbon_only'])}** 个只剩 C 谱、**{fmt(scopes['original_native_dual']['proton_only'])}** 个只剩 H 谱。当前两条原始谱本身都仍有效的结构是 **{fmt(native['currently_both_original_spectra_active'])}** 个，差异来自人工修订版本。", '',
        f"同源双谱不直接证明实验条件一致：现有原始配对表中，明确标记已知相同溶剂的原始双谱结构有 **{fmt(native['same_record_known_same_solvent_structures'])}** 个；未知溶剂未自动当作一致。这个数量受已有配对表覆盖范围限制，不是对全部原始谱重新穷举条件匹配的结果。", '',
        '## 4. 是否存在单谱', '',
        f"存在。目前只有 C 谱 **{fmt(report['current_spectrum_kinds']['carbon_only'])}** 个，只有 H 谱 **{fmt(report['current_spectrum_kinds']['proton_only'])}** 个，共 **{fmt(report['current_spectrum_kinds']['carbon_only']+report['current_spectrum_kinds']['proton_only'])}** 个。", '',
        f"若连隔离、参考和历史版本也算入存档，则当前存档仍仅有一类谱的结构为：C 谱 **{fmt(report['raw_spectrum_kinds']['carbon_only'])}** 个，H 谱 **{fmt(report['raw_spectrum_kinds']['proton_only'])}** 个，共 **{fmt(report['raw_spectrum_kinds']['carbon_only']+report['raw_spectrum_kinds']['proton_only'])}** 个。当前可用单谱比存档单谱多，是因为部分分子的另一核种被停用或隔离。", '',
        '## 5. 文件与校验', '',
        '- `statistics.json`：精确统计、口径、来源交叉表和原子身份门控核对。',
        '- `annotated_molecules.csv`：全部至少有部分归属的分子清单。',
        '- `human_participation.csv`：人工参与分子，列出来源、归属程度、单 / 双谱、最近操作、是否原始 ShiftDB2 双谱。',
        '- `original_shiftdb_dual.csv`：原始 ShiftDB2 双谱结构清单。',
        '- `current_single_spectrum_molecules.csv`：当前单谱分子清单。', '',
        f"核对了全部 {fmt(sum(report['source_identity_audit'].values()))} 条谱的已有原子身份审计状态，并对 {fmt(sum(v for k,v in report['source_identity_audit'].items() if k!='identity_validated'))} 条非通过谱涉及的分子及所有人工参与分子，使用当前界面的有效峰投影重新计数。本次缓存汇总与有效投影的分类差异为 **{fmt(len(report['cached_vs_effective_changed_molecules']))}** 个分子。", '',
        '统计过程没有改写分子、谱表、归属、任务状态或数据库结构。', '',
    ]
    (output / 'statistics.md').write_text('\n'.join(lines), encoding='utf-8')


def run(path, output):
    started = time.perf_counter()
    c = sqlite3.connect(f'file:{path.resolve().as_posix()}?mode=ro', uri=True, timeout=30)
    c.row_factory = sqlite3.Row
    c.execute('BEGIN')
    snapshot = datetime.now(timezone(timedelta(hours=8))).isoformat(timespec='seconds')
    # Establish the snapshot before reading any other table.
    total_molecules = c.execute('SELECT count(*) FROM molecules').fetchone()[0]
    log(f'Snapshot {snapshot}: {total_molecules:,} molecules')
    brief = 'id,formula,smiles,origin,queue,task_done,c_count,h_count,c_assigned,h_assigned,revision'
    aggregate = Counter()
    sql = f'''SELECT origin,c_count>0 AS cn,h_count>0 AS hn,({db.COVERAGE_SQL}) AS coverage,
              task_done,count(*) AS n FROM molecules GROUP BY origin,cn,hn,coverage,task_done'''
    for r in c.execute(sql):
        aggregate[(r['origin'], kind(r['cn'], r['hn']), r['coverage'], r['task_done'])] += r['n']
    cached_aggregate = aggregate.copy()
    log('Read current spectrum/coverage census')

    annotation_rows = c.execute('''SELECT s.molecule_id,a.spectrum_id,a.peak_id,a.origin,a.status,a.updated_at
          FROM annotations a JOIN spectra s ON s.id=a.spectrum_id
          WHERE a.origin IN ('manual','manual_revision')''').fetchall()
    manual_ids = {r['molecule_id'] for r in annotation_rows}
    manual_counts = Counter((r['origin'], r['status']) for r in annotation_rows)
    activity = {r['molecule_id']: dict(r) for r in c.execute('SELECT * FROM molecule_activity')}
    event_ids = {r[0] for r in c.execute('SELECT DISTINCT molecule_id FROM curation_events')}
    numbering_ids = {r[0] for r in c.execute('SELECT DISTINCT molecule_id FROM spectrum_numbering_reviews')}
    human_ids = manual_ids | set(activity) | event_ids | numbering_ids
    management_only = human_ids - manual_ids
    spectrum_states = {r['target_id']: r['state'] for r in c.execute("SELECT target_id,state FROM curation_state WHERE scope='spectrum'")}
    molecule_states = {r['target_id']: r['state'] for r in c.execute("SELECT target_id,state FROM curation_state WHERE scope='molecule'")}
    qualities = {r['spectrum_id']: r['decision'] for r in c.execute('SELECT spectrum_id,decision FROM spectrum_quality')}
    audit_statuses = Counter()
    blocked_sids = set()
    for r in c.execute('SELECT spectrum_id,validation_json FROM spectrum_assignment_audit'):
        status = json.loads(r['validation_json'])['status']
        audit_statuses[status] += 1
        if status != 'identity_validated':
            blocked_sids.add(r['spectrum_id'])
    log(f'Human activity: {len(human_ids)} molecules; {len(blocked_sids)} spectra require source-identity gating')

    native_nuclei = defaultdict(set)
    native_current_nuclei = defaultdict(set)
    native_spectra = {}
    native_source_rows = defaultdict(set)
    # The stable source-key range uses the UNIQUE source_key index and avoids a
    # full scan of >1M NMRexp provenance documents.
    for r in c.execute('''SELECT id,molecule_id,nucleus,source,source_key,peaks_json
           FROM spectra WHERE source_key>='shiftdb:' AND source_key<'shiftdb;' '''):
        if r['source'] != 'nmrshiftdb2':
            continue
        nonempty = bool(json.loads(r['peaks_json']))
        native_spectra[r['id']] = (r['molecule_id'], r['nucleus'], nonempty)
        if not nonempty:
            continue
        native_nuclei[r['molecule_id']].add(r['nucleus'])
        native_source_rows[(r['molecule_id'], r['source_key'].rsplit(':Spectrum', 1)[0])].add(r['nucleus'])
        if spectrum_states.get(r['id'], 'active') == 'active' and qualities.get(r['id'], 'passed') == 'passed':
            native_current_nuclei[r['molecule_id']].add(r['nucleus'])
    native_dual = {mid for mid, nuclei in native_nuclei.items() if len(nuclei) == 2}
    native_same_record = {mid for (mid, _), nuclei in native_source_rows.items() if len(nuclei) == 2}
    native_current_dual = {mid for mid, nuclei in native_current_nuclei.items() if len(nuclei) == 2}
    native_pair_mids, native_solvent_mids, current_native_pair_mids = set(), set(), set()
    for r in c.execute("SELECT * FROM pairs WHERE basis='same_shiftdb_record_sample_unknown'"):
        if r['c_id'] not in native_spectra or r['h_id'] not in native_spectra:
            continue
        native_pair_mids.add(r['molecule_id'])
        if r['solvent_match'] == 'same_known':
            native_solvent_mids.add(r['molecule_id'])
        if all(spectrum_states.get(sid, 'active') == 'active' and qualities.get(sid, 'passed') == 'passed' for sid in (r['c_id'], r['h_id'])):
            current_native_pair_mids.add(r['molecule_id'])
    log(f'Original ShiftDB2 structures: {len(native_nuclei):,}; original dual: {len(native_dual):,}')

    records = {r['id']: dict(r) for r in c.execute(f'''SELECT {brief} FROM molecules
              WHERE ({db.COVERAGE_SQL}) IN ('full','partial')''')}
    for mid in (set(native_nuclei) | human_ids) - records.keys():
        row = c.execute(f'SELECT {brief} FROM molecules WHERE id=?', (mid,)).fetchone()
        if row:
            records[mid] = dict(row)
    log(f'Read {len(records):,} labelled/native/human molecule summaries')

    affected = human_ids.copy()
    for sid in blocked_sids:
        row = c.execute('SELECT molecule_id FROM spectra WHERE id=?', (sid,)).fetchone()
        if row:
            affected.add(row[0])
    changed = []
    manual_active_ids = set()
    active_groups = Counter()
    effective_records = {}
    for index, mid in enumerate(sorted(affected), 1):
        if mid not in records:
            records[mid] = dict(c.execute(f'SELECT {brief} FROM molecules WHERE id=?', (mid,)).fetchone())
        row = records[mid]
        old = row_key(row)
        value = db.record(c, mid)
        effective_records[mid] = value
        for nuc, prefix in (('13C', 'c'), ('1H', 'h')):
            spectra = [s for s in value['spectra'] if s['nucleus'] == nuc]
            peaks = [p for s in spectra for p in s['peaks']]
            row[prefix + '_count'] = len(peaks)
            row[prefix + '_assigned'] = sum(p['assignment']['status'] in db.ASSIGNED for p in peaks)
            for spec in spectra:
                active_groups['spectra_with_unordered_groups'] += bool(spec.get('assignment_groups'))
                active_groups['unordered_groups'] += len(spec.get('assignment_groups', []))
                for p in spec['peaks']:
                    if p['assignment'].get('origin') in ('manual', 'manual_revision'):
                        manual_active_ids.add(mid)
        new = row_key(row)
        if old != new:
            aggregate[old] -= 1
            aggregate[new] += 1
            changed.append({'id': mid, 'formula': row['formula'], 'old': list(old), 'effective': list(new)})
        if index % 100 == 0:
            log(f'Validated effective display coverage {index}/{len(affected)}')

    for mid, row in records.items():
        row['coverage'] = coverage(row['c_count'] + row['h_count'], row['c_assigned'] + row['h_assigned'])
        row['kind'] = kind(row['c_count'], row['h_count'])
        row['manual_annotation_history'] = mid in manual_ids
        row['human_participation'] = mid in human_ids
        row['manual_labels_current'] = mid in manual_active_ids
        row['original_native_dual'] = mid in native_dual
        row['original_native_same_record_dual'] = mid in native_same_record
        row['molecule_use_state'] = molecule_states.get(mid, 'active')
        row['last_processed_at'] = activity.get(mid, {}).get('last_processed_at', '')
        row['last_action'] = activity.get(mid, {}).get('last_action', '')

    raw_kinds = Counter()
    for row in c.execute('''SELECT has_c,has_h,count(*) n FROM (
          SELECT molecule_id,max(nucleus='13C') has_c,max(nucleus='1H') has_h
          FROM spectra INDEXED BY spectrum_molecule GROUP BY molecule_id)
          GROUP BY has_c,has_h'''):
        raw_kinds[kind(row['has_c'], row['has_h'])] += row['n']
    log('Validated raw stored spectrum nuclei')
    all_native = [records[mid] for mid in native_nuclei]
    manual = [records[mid] for mid in manual_ids]
    human = [records[mid] for mid in human_ids]
    scopes = {
        'human_participation': human,
        'manual_annotation_history': manual,
        'human_management_or_save_only': [records[mid] for mid in management_only],
        'original_native_dual': [records[mid] for mid in native_dual],
        'native_dual_without_human_participation': [records[mid] for mid in native_dual-human_ids],
        'native_dual_with_human_participation': [records[mid] for mid in native_dual&human_ids],
        'all_native_structures': all_native,
    }
    for origin in ('nmrshiftdb2', 'NMRexp'):
        scopes['manual_' + origin] = [r for r in manual if r['origin'] == origin]
        scopes['human_' + origin] = [r for r in human if r['origin'] == origin]
    cross = []
    for (origin, nuclei, cov, done), n in sorted(aggregate.items()):
        if n:
            cross.append({'origin': origin, 'kind': nuclei, 'coverage': cov, 'task_done': done, 'molecules': n})
    total_cov = Counter()
    total_kinds = Counter()
    for r in cross:
        total_cov[r['coverage']] += r['molecules']
        total_kinds[r['kind']] += r['molecules']
    labelled = [r for r in records.values() if r['coverage'] != 'none']
    report = {
        'snapshot_at': snapshot, 'database': str(path.resolve()),
        'definitions': {
            'unit': 'Distinct molecules.id (canonical isomeric SMILES), not spectrum records or atom counts',
            'current_spectra': 'core.record default: active manual state and not automatically quarantined; NMRexp not_checked is retained as in UI',
            'full': 'Every reported current peak has assigned/equivalent_group/overlap/source_group status; does not certify physical spectral completeness or correctness',
            'partial': 'At least one such assigned peak, with at least one current peak not assigned; unordered groups and ambiguity do not count as hard assignments',
            'none': 'No current peak has one of these assignment statuses',
            'human_participation': 'Union of manual/manual_revision annotations, molecule_activity, curation_events and numbering reviews; includes note-only and status operations',
            'manual_annotation_history': 'At least one saved manual/manual_revision peak annotation, including ambiguous/cleared/impurity states and archived versions',
            'original_native_dual': 'Both nonempty C/H source spectra originally stored as source=nmrshiftdb2, independent of current reference/exclusion states',
            'raw_nuclei': 'All stored spectra, including quarantined, reference and historical versions',
            'correction': 'Cached summary counts corrected by core.record projections for all human-touched molecules and all source-identity-blocked spectra',
        },
        'molecules_total': total_molecules,
        'current_coverage': dict(total_cov), 'current_spectrum_kinds': dict(total_kinds),
        'raw_spectrum_kinds': dict(raw_kinds),
        'labelled_molecules': summarize(labelled),
        'cross_tabulation': cross,
        'scope_summaries': {name: summarize(rows) for name, rows in scopes.items()},
        'native_original': {
            'structures': len(native_nuclei), 'spectra': len(native_spectra), 'dual_structures': len(native_dual),
            'same_record_dual_structures': len(native_same_record),
            'same_record_paired_structures': len(native_pair_mids),
            'same_record_known_same_solvent_structures': len(native_solvent_mids),
            'currently_both_original_spectra_active': len(native_current_dual),
            'currently_same_record_pair_active': len(current_native_pair_mids),
        },
        'manual_annotation_rows': [{'origin': o, 'status': s, 'rows': n} for (o, s), n in sorted(manual_counts.items())],
        'manual_active_molecules': len(manual_active_ids),
        'source_identity_audit': dict(audit_statuses), 'unordered_groups_in_checked_molecules': dict(active_groups),
        'cached_vs_effective_changed_molecules': changed,
        'cached_coverage': dict(Counter({cov: sum(n for (o,k,s,d),n in cached_aggregate.items() if s==cov) for cov in ('full','partial','none')})),
        'query_elapsed_seconds': round(time.perf_counter() - started, 3),
    }
    assert sum(total_cov.values()) == total_molecules
    assert sum(total_kinds.values()) == total_molecules
    assert total_cov['full'] + total_cov['partial'] == len(labelled)
    assert len(human_ids) == len(manual_ids) + len(management_only)
    assert sum(raw_kinds.values()) <= total_molecules
    c.rollback()
    c.close()
    output.mkdir(parents=True, exist_ok=True)
    (output / 'statistics.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    columns = ['id','formula','smiles','origin','kind','coverage','c_count','c_assigned','h_count','h_assigned',
               'task_done','manual_annotation_history','manual_labels_current','human_participation',
               'original_native_dual','original_native_same_record_dual','molecule_use_state','last_processed_at','last_action']
    for filename, values in (('annotated_molecules.csv', labelled), ('human_participation.csv', human),
                             ('original_shiftdb_dual.csv', scopes['original_native_dual']),
                             ('current_single_spectrum_molecules.csv',
                              [r for r in records.values() if r['kind'] in ('carbon_only','proton_only')])):
        with (output / filename).open('w', encoding='utf-8-sig', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(sorted(values, key=lambda r: r['id']))
    write_markdown(report, output)
    log(json.dumps({key: report[key] for key in ('snapshot_at','current_coverage','current_spectrum_kinds','raw_spectrum_kinds','scope_summaries','native_original','cached_coverage')}, ensure_ascii=False, indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', type=Path, default=db.DEFAULT_DB)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    run(args.db, args.output_dir)


if __name__ == '__main__':
    main()
