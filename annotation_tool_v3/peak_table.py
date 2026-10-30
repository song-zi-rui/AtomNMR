"""Experimental peak supplements projected over immutable source peak tables."""
from __future__ import annotations

import copy
import json
import math
import uuid

from . import core as db

KEY = 'peak_table_supplements'
FIELDS = {'shift_ppm', 'range_min', 'range_max', 'multiplicity', 'j_values',
          'reported_nH', 'integral', 'line_width_hz'}


def effective_peaks(peaks, provenance):
    result = copy.deepcopy(peaks)
    by_id = {p['peak_id']: p for p in result}
    for item in provenance.get(KEY, []):
        if not item.get('active', True):
            continue
        pid = item['peak_id']
        if item['kind'] == 'delete':
            if pid not in by_id:
                raise ValueError('删除修订引用的峰已不存在。')
            result = [p for p in result if p['peak_id'] != pid]
            by_id.pop(pid)
        elif item['kind'] == 'add':
            if pid in by_id:
                raise ValueError('补充峰编号与源峰重复。')
            peak = {'peak_id': pid, **copy.deepcopy(item['fields']),
                    'experimental_peak_added': True, 'peak_supplement': copy.deepcopy(item)}
            result.append(peak)
            by_id[pid] = peak
        else:
            if pid not in by_id:
                raise ValueError('补充字段引用的源峰已不存在，请核实来源。')
            by_id[pid].update(copy.deepcopy(item['fields']))
            by_id[pid].setdefault('field_supplements', []).append(copy.deepcopy(item))
    return result


def normalized_fields(value, allow_clear=False):
    if not isinstance(value, dict) or not value or set(value) - FIELDS:
        raise ValueError('请选择要补充的峰表字段；原子归属请在标注界面填写。')
    result = {}
    for key, raw in value.items():
        if raw is None and allow_clear and key != 'shift_ppm':
            result[key] = None
            continue
        if raw is None or raw == '' or raw == []:
            continue
        if key == 'multiplicity':
            if not isinstance(raw, str) or not raw.strip() or len(raw.strip()) > 80:
                raise ValueError('峰型应为不超过 80 字的实验记录。')
            result[key] = raw.strip()
            continue
        values = raw if key == 'j_values' else [raw]
        if not isinstance(values, list) or not values or len(values) > 100:
            raise ValueError('J 值应为 Hz 单位的数值列表。')
        numbers = []
        for v in values:
            if isinstance(v, bool):
                raise ValueError('峰表数值不能是布尔值。')
            try:
                number = float(v)
            except (ValueError, TypeError):
                raise ValueError('峰表数值格式不正确。') from None
            if not math.isfinite(number) or (key in {'j_values', 'reported_nH', 'integral', 'line_width_hz'} and number < 0):
                raise ValueError('数值必须有限；J、氢数、积分及线宽不能为负。')
            numbers.append(number)
        result[key] = numbers if key == 'j_values' else numbers[0]
    if not result:
        raise ValueError('请至少填写一个缺失值。')
    return result


def mutate(c, mid, spec, payload, reason):
    from .spectrum_revision import require_editable, group_for
    require_editable(c, mid, spec)
    provenance = json.loads(spec['source_json'])
    items = copy.deepcopy(provenance.get(KEY, []))
    peaks = effective_peaks(json.loads(spec['peaks_json']), provenance)
    by_id = {p['peak_id']: p for p in peaks}
    name = payload['action']
    if name in ('supplement_peak', 'edit_peak', 'delete_peak'):
        pid = payload.get('peak_id')
        kind = 'delete' if name == 'delete_peak' else 'edit' if name == 'edit_peak' else 'fields' if pid else 'add'
        if kind in ('edit', 'delete') and pid not in by_id:
            raise ValueError('请选择当前谱中要修改或删除的峰。')
        fields = {} if kind == 'delete' else normalized_fields(payload.get('fields'), allow_clear=kind == 'edit')
        if kind == 'delete':
            if group_for(provenance, pid):
                raise ValueError('请先解除该峰所在的位移组，再删除峰。')
            combined = by_id[pid]
        elif pid:
            if pid not in by_id:
                raise ValueError('补充目标不属于当前谱。')
            if 'shift_ppm' in fields and kind != 'edit':
                raise ValueError('补充缺项不能修改已有位移；缺失信号请选择新增实验峰。')
            if kind != 'edit' and any(by_id[pid].get(k) not in (None, '', []) for k in fields):
                raise ValueError('所选字段已有数值，补充操作不能覆盖；如要改正先撤销相应补充。')
            combined = {**by_id[pid], **fields}
            if 'shift_ppm' in fields and any(p['peak_id'] != pid and p['shift_ppm'] == fields['shift_ppm'] for p in peaks) and payload.get('allow_same_shift') is not True:
                raise ValueError('修改后与其他峰同位移，请明确勾选保留独立信号。')
        else:
            if 'shift_ppm' not in fields:
                raise ValueError('新增实验峰必须填写化学位移 ppm。')
            if any(p['shift_ppm'] == fields['shift_ppm'] for p in peaks) and payload.get('allow_same_shift') is not True:
                raise ValueError('当前谱已有同位移项。如确为另一个实验信号，请明确勾选保留同位移的独立记录。')
            pid = ('C' if spec['nucleus'] == '13C' else 'H') + 'M' + uuid.uuid4().hex[:10]
            combined = fields
        low, high = combined.get('range_min'), combined.get('range_max')
        if ((low is not None and low > combined['shift_ppm'])
                or (high is not None and high < combined['shift_ppm'])):
            raise ValueError('范围下限应不大于中心位移，范围上限应不小于中心位移。')
        item = {'id': 'PS' + uuid.uuid4().hex[:24], 'kind': kind, 'peak_id': pid,
                'fields': fields, 'reason': reason, 'source_type': 'manual_experimental_supplement',
                'active': True, 'created_at': db.now()}
        items.append(item)
        detail = {'supplement': copy.deepcopy(item), 'added_peak_id': pid if kind == 'add' else None}
    else:
        item = next((x for x in items if x['id'] == payload.get('supplement_id') and x.get('active', True)), None)
        if item is None:
            raise ValueError('补充记录不存在或已撤销。')
        pid = item['peak_id']
        if item['kind'] == 'add':
            annotation = c.execute('SELECT * FROM annotations WHERE spectrum_id=? AND peak_id=?', (spec['id'], pid)).fetchone()
            if annotation and (annotation['status'] != 'unassigned' or json.loads(annotation['atoms_json']) or annotation['note']):
                raise ValueError('新增峰已有归属、处理状态或备注，请先清空后再撤销该峰。')
            if group_for(provenance, pid) or pid in (provenance.get('spectrum_revisions') or {}).get('signal_codes', {}):
                raise ValueError('请先解除该峰的位移组或编码修订，再撤销新增峰。')
            if any(x['kind'] != 'add' and x['peak_id'] == pid and x.get('active', True) for x in items):
                raise ValueError('请先撤销该峰后续的字段补充，再撤销新增峰。')
        item.update(active=False, withdrawn_at=db.now(), withdrawal_reason=reason)
        detail = {'withdrawn_supplement': copy.deepcopy(item)}
    provenance[KEY] = items
    c.execute('UPDATE spectra SET source_json=? WHERE id=?', (db.dumps(provenance), spec['id']))
    # Adding/withdrawing an observed signal changes what still needs checking.
    if item['kind'] in ('add', 'delete') or 'shift_ppm' in item['fields']:
        c.execute('INSERT INTO spectrum_completeness VALUES(?,?,?,?) ON CONFLICT(spectrum_id) DO UPDATE SET status=excluded.status,reason=excluded.reason,updated_at=excluded.updated_at',
                  (spec['id'], 'unknown', '实验峰表已补充或撤销补充，请重新核实完整性。', db.now()))
    return detail
