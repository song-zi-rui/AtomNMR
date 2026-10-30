"""Bounded official BMRB import, named-atom graph mapping, reference-only admission.

This deliberately supports a conservative subset of NMR-STAR. Unsupported or
ambiguous sample/component associations fail closed instead of flattening tables.
Assigned shifts are NOT advertised as a complete experimental 1D peak list.
"""
from __future__ import annotations

import hashlib
import math
import re
import urllib.request
from collections import defaultdict

from rdkit import Chem
from . import core as db


def tokens(text):
    i, n = 0, len(text)
    while i < n:
        if text[i].isspace():
            i += 1
            continue
        if text[i] == '#':
            end = text.find('\n', i)
            i = n if end < 0 else end + 1
            continue
        if text[i] == ';' and (i == 0 or text[i - 1] == '\n'):
            end = re.search(r'\n;(?=\s|$)', text[i + 1:])
            if not end:
                raise ValueError('STAR 多行文本未闭合。')
            stop = i + 1 + end.start()
            yield text[i + 1:stop], True
            i = stop + 2
        elif text[i] in "\"'":
            quote, start = text[i], i + 1
            i += 1
            while i < n and not (text[i] == quote and (i + 1 == n or text[i + 1].isspace())):
                i += 1
            if i == n:
                raise ValueError('STAR 引号未闭合。')
            yield text[start:i], True
            i += 1
        else:
            start = i
            while i < n and not text[i].isspace():
                i += 1
            yield text[start:i], False


def parse_star(text):
    ts, frames, current, i = list(tokens(text)), [], None, 0
    while i < len(ts):
        value, quoted = ts[i]
        i += 1
        if quoted:
            raise ValueError('无法解析 STAR 文件结构。')
        if value.startswith('data_'):
            continue
        if value.startswith('save_'):
            current = None if value == 'save_' else {'name': value[5:], 'tags': {}, 'loops': []}
            if current is not None:
                frames.append(current)
        elif current is None:
            raise ValueError('STAR 数据不在 saveframe 内。')
        elif value == 'loop_':
            headers, rows = [], []
            while i < len(ts) and not ts[i][1] and ts[i][0].startswith('_'):
                headers.append(ts[i][0]); i += 1
            if not headers:
                raise ValueError('STAR 表头为空。')
            values = []
            while i < len(ts) and ts[i] != ('stop_', False):
                if not ts[i][1] and (ts[i][0].startswith(('_', 'save_', 'data_')) or ts[i][0] == 'loop_'):
                    raise ValueError('STAR 表格未正确结束。')
                values.append(ts[i][0]); i += 1
            if i == len(ts) or len(values) % len(headers):
                raise ValueError('STAR 表格列数不一致。')
            i += 1
            for start in range(0, len(values), len(headers)):
                rows.append(dict(zip(headers, values[start:start + len(headers)])))
            current['loops'].append(rows)
        elif value.startswith('_'):
            if i == len(ts) or value in current['tags']:
                raise ValueError('STAR 重复或缺失字段。')
            current['tags'][value] = ts[i][0]; i += 1
        else:
            raise ValueError('不支持的 STAR 控制字段：' + value)
    return frames


def table(frames, prefix):
    return [{k.removeprefix(prefix + '.'): v for k, v in row.items()} for frame in frames
            for loop in frame['loops'] for row in loop if any(k.startswith(prefix + '.') for k in row)]


def plain_smiles(mol):
    m = Chem.RemoveHs(Chem.Mol(mol))
    Chem.RemoveStereochemistry(m)
    return Chem.MolToSmiles(m)


def build_reference(star, molblock, local_smiles, entry_id, metadata=None):
    frames = parse_star(star)
    tags = {k: v for f in frames for k, v in f['tags'].items() if k.startswith('_Entry.')}
    if tags.get('_Entry.ID') != entry_id or tags.get('_Entry.Experimental_method') != 'NMR':
        raise ValueError('参考不是指定条目的实验 NMR 数据。')
    atoms, bonds = table(frames, '_Chem_comp_atom'), table(frames, '_Chem_comp_bond')
    if not atoms or len({a['Comp_ID'] for a in atoms}) != 1:
        raise ValueError('暂不支持多化学组件表或缺少命名原子表的条目。')
    rw, names = Chem.RWMol(), {}
    for row in atoms:
        name = row['Atom_ID']
        if name in names:
            raise ValueError('BMRB 原子名重复。')
        a = Chem.Atom(row['Type_symbol'])
        charge = row.get('Charge', '.')
        if charge not in ('.', '?'):
            a.SetFormalCharge(int(charge))
        a.SetProp('_reference_name', name)
        names[name] = rw.AddAtom(a)
    orders = {'SING': Chem.BondType.SINGLE, 'DOUB': Chem.BondType.DOUBLE, 'TRIP': Chem.BondType.TRIPLE}
    for b in bonds:
        if b['Value_order'].upper() not in orders:
            raise ValueError('参考含暂不支持的键类型，需要手工核对。')
        rw.AddBond(names[b['Atom_ID_1']], names[b['Atom_ID_2']], orders[b['Value_order'].upper()])
    named = rw.GetMol()
    Chem.SanitizeMol(named)
    mol = Chem.MolFromMolBlock(molblock, removeHs=False)
    local = Chem.MolFromSmiles(local_smiles)
    if mol is None or local is None or plain_smiles(mol) != plain_smiles(named) or plain_smiles(named) != plain_smiles(local):
        raise ValueError('完整连接图、组成、电荷或同位素不一致，不能按分子式导入参考。')
    ref_heavy = Chem.RemoveHs(mol)
    if not ref_heavy.HasSubstructMatch(local, useChirality=True):
        raise ValueError('当前明确的立体结构与 BMRB 不一致或参考未明确，不能导入。')
    heavy = Chem.RemoveHs(named)
    matches = local.GetSubstructMatches(heavy, uniquify=False, useChirality=False, maxMatches=4097)
    if not matches or len(matches) >= 4097:
        raise ValueError('原子映射过多或无法建立，不能安全导入。')
    mapped = {a.GetProp('_reference_name'): sorted({m[a.GetIdx()] for m in matches}) for a in heavy.GetAtoms()}
    parent = {}
    for a in named.GetAtoms():
        if a.GetAtomicNum() == 1 and a.GetDegree() == 1:
            parent[a.GetProp('_reference_name')] = a.GetNeighbors()[0].GetProp('_reference_name')
    warnings = ['这是实验原子位移归属表，不是完整 1D 峰表；未出现的位移不会被推算补入。',
                'H 按连接重原子表示；同一重原子上同位移的 H 合并，不把原子数冒充实验积分。',
                '结构映射按命名原子和键建立，不按分子式或原子名数字猜编号。']
    stereo_diff = Chem.MolToSmiles(ref_heavy) != Chem.MolToSmiles(local)
    if stereo_diff:
        warnings.append('BMRB 立体信息比当前结构更具体：仅建立连接图映射，不改写当前结构或声称已确定其立体化学。')
    samples = table(frames, '_Sample_component')
    condition_rows = table(frames, '_Sample_condition_variable')
    spectra = []
    for frame in frames:
        list_id = frame['tags'].get('_Assigned_chem_shift_list.ID')
        if not list_id:
            continue
        shifts = table([frame], '_Atom_chem_shift')
        if len({r.get('ID') for r in shifts}) != len(shifts):
            raise ValueError('位移源行编号重复，不能保证逐行追溯。')
        experiments = table([frame], '_Chem_shift_experiment')
        sample_ids = {e.get('Sample_ID') for e in experiments}
        if len(sample_ids) != 1 or sample_ids & {None, '.', '?'}:
            raise ValueError('归属表对应样品不唯一，暂不能安全引入。')
        sample_id = next(iter(sample_ids))
        cond_id = frame['tags'].get('_Assigned_chem_shift_list.Sample_condition_list_ID')
        solvents = [r['Mol_common_name'] for r in samples if r['Sample_ID'] == sample_id and r.get('Type', '').lower() == 'solvent']
        conditions = {'solvent': solvents[0] if len(solvents) == 1 else None,
                      'sample_id': sample_id, 'sample_condition_list_id': cond_id,
                      'sample_components': [r for r in samples if r['Sample_ID'] == sample_id],
                      'frequency_mhz': None, 'temperature_k': None}
        for kind, key, unit in [('temperature', 'temperature_k', 'K'), ('pH', 'ph', 'pH')]:
            rows = [r for r in condition_rows if r['Sample_condition_list_ID'] == cond_id and r['Type'] == kind and r['Val_units'] == unit]
            if len(rows) == 1:
                number = float(rows[0]['Val'])
                if not math.isfinite(number):
                    raise ValueError('非有限实验条件。')
                conditions[key] = number
        for nucleus, symbol, isotope in [('13C', 'C', '13'), ('1H', 'H', '1')]:
            selected = [r for r in shifts if r['Atom_type'] == symbol and r['Atom_isotope_number'] == isotope]
            if not selected:
                continue
            evidence = [e for e in experiments if re.fullmatch(r'1D\s+' + nucleus, e.get('Experiment_name', ''), re.I)]
            if not evidence:
                raise ValueError('该归属表缺少相应核种的一维实验关联。')
            grouped = defaultdict(list)
            for row in selected:
                if row.get('Comp_ID') != atoms[0]['Comp_ID'] or row.get('Assigned_chem_shift_list_ID') != list_id:
                    raise ValueError('位移表的化学组件或归属表 ID 不一致。')
                atom_id = row['Atom_ID']
                if atom_id not in names or named.GetAtomWithIdx(names[atom_id]).GetSymbol() != symbol:
                    raise ValueError('位移表原子核种与命名结构原子不一致。')
                name = parent.get(row['Atom_ID']) if symbol == 'H' else row['Atom_ID']
                if name not in mapped:
                    raise ValueError('位移原子不在命名结构图中。')
                value = float(row['Val'])
                if not math.isfinite(value):
                    raise ValueError('非有限化学位移。')
                grouped[(name, value)].append(row)
            peaks, assignments = [], {}
            for (name, value), rows in grouped.items():
                pid = 'B' + str(len(peaks) + 1)
                targets = mapped[name]
                unique = len(targets) == 1 and all(r.get('Ambiguity_code') == '1' for r in rows)
                peaks.append({'peak_id': pid, 'shift_ppm': value, 'source_atom_ids': [r['Atom_ID'] for r in rows],
                              'source_row_ids': [r['ID'] for r in rows], 'raw_rows': rows, 'reported_nH': None})
                assignments[pid] = {'status': 'assigned' if unique else 'ambiguous', 'atom_indices': targets,
                                    'note': 'BMRB 命名原子映射；多候选保留歧义，不随机选择。'}
            spectra.append({'nucleus': nucleus, 'list_id': list_id, 'conditions': {**conditions, 'experiments': evidence},
                            'peaks': peaks, 'assignments': assignments})
    if not spectra:
        raise ValueError('没有可关联的实验归属表。')
    digest = hashlib.sha256((star + '\0' + molblock).encode()).hexdigest()
    provenance = {'entry_id': entry_id, 'doi': tags.get('_Entry.DOI'), 'title': tags.get('_Entry.Title'),
                  'source_url': f'https://bmrb.io/metabolomics/mol_summary/show_data.php?id={entry_id}',
                  'source_kind': 'experimental_assigned_shift_table', 'sha256': digest,
                  'download': metadata or {}, 'raw_star': star, 'raw_mol': molblock,
                  'reference_smiles': Chem.MolToSmiles(ref_heavy), 'local_smiles': local_smiles,
                  'named_atom_mapping': mapped, 'warnings': warnings,
                  'mapping_scope': 'connectivity_only', 'stereo_more_specific': stereo_diff}
    return {'entry_id': entry_id, 'local_smiles': local_smiles, 'sha256': digest, 'spectra': spectra, 'provenance': provenance}


def fetch_reference(entry_id, local_smiles):
    if not isinstance(entry_id, str) or not re.fullmatch(r'bmse\d{6}', entry_id):
        raise ValueError('请输入 BMRB Metabolomics 编号，例如 bmse000194。')
    content, metadata = {}, {}
    for suffix in ('str', 'mol'):
        url = f'https://bmrb.io/ftp/pub/bmrb/metabolomics/entry_directories/{entry_id}/{entry_id}.{suffix}'
        try:
            with urllib.request.urlopen(url, timeout=25) as response:
                raw = response.read(8 * 1024 * 1024 + 1)
                if len(raw) > 8 * 1024 * 1024:
                    raise ValueError('参考文件过大，暂停导入。')
                metadata[suffix] = {'url': url, 'sha256': hashlib.sha256(raw).hexdigest(),
                                    'last_modified': response.headers.get('Last-Modified'), 'fetched_at': db.now()}
                content[suffix] = raw.decode('utf-8')
        except (OSError, UnicodeError) as exc:
            raise ValueError(f'BMRB 下载失败，原库未修改：{exc}') from exc
    return build_reference(content['str'], content['mol'], local_smiles, entry_id, metadata)
