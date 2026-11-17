"""Folder extraction and portable AtomNMR database output."""

from collections import defaultdict
from contextlib import closing
from hashlib import file_digest
from importlib.util import find_spec
from pathlib import Path
import json
import os
import re
import subprocess

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

from annotation_tool_v3 import ai_database, core as db
from annotation_tool_v3.spectrum_revision import KEY as REVISIONS
from .pdf_tables import extract_document, probe_document

HERE = Path(__file__).resolve().parent
SAMPLE = json.loads((HERE / 'mauritines.json').read_text(encoding='utf-8'))
COMPARISON_SAMPLE = json.loads((HERE / 'lyngbyabellins.json').read_text(encoding='utf-8'))
SOURCE = 'literature_pdf'


def digest(path):
    with path.open('rb') as handle:
        return file_digest(handle, 'sha256').hexdigest()


def source_paths(folder):
    return sorted((p for p in folder.iterdir() if p.is_file() and not p.is_symlink()
                   and p.suffix.lower() in ('.pdf', '.sdf')), key=lambda p: p.name.casefold())


def pdf_helper(path, mode):
    configured = os.environ.get('ATOMNMR_PDF_PYTHON')
    bundled = Path.home() / '.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
    python = Path(configured) if configured else bundled
    if not python.is_file():
        raise RuntimeError('找不到 PDF 解析器。请在启动环境运行 python -m pip install pdfplumber。')
    result = subprocess.run([str(python), str(HERE / 'pdf_tables.py'), str(path), mode],
                            capture_output=True, encoding='utf-8', errors='replace', timeout=180,
                            creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
    if result.returncode != 0:
        message = result.stderr.strip().splitlines()[-1:] or ['PDF 读取失败']
        raise ValueError(message[0])
    return json.loads(result.stdout)


def read_pdf(path, checksum, profiles):
    with path.open('rb') as handle:
        if b'%PDF-' not in handle.read(1024):
            raise ValueError('文件没有 PDF 标头，可能下载成了网页或错误响应')
        handle.seek(max(0, path.stat().st_size - 65536))
        if not re.search(rb'(?m)^%%EOF[ \t\r]*$', handle.read().replace(b'\r', b'\n')):
            raise ValueError('PDF 文件尾缺少 %%EOF，疑似下载被截断。请重新下载完整 PDF。')
    return extract_document(path, profiles) if find_spec('pdfplumber') else pdf_helper(path, checksum)


def inspect_pdf(path):
    return probe_document(path) if find_spec('pdfplumber') else pdf_helper(path, 'probe')


def read_sdf(files):
    structures, issues = defaultdict(list), []
    for path in files:
        try:
            with path.open('rb') as handle:
                for number, mol in enumerate(Chem.ForwardSDMolSupplier(handle), 1):
                    if mol is None:
                        issues.append(f'{path.name} 第 {number} 条结构无法读取')
                        continue
                    props = {key.lower(): value.strip() for key, value in mol.GetPropsAsDict(
                        autoConvertStrings=False).items() if value.strip()}
                    scope = props.get('pdf') or props.get('doi', '').lower()
                    table = props.get('table', '').upper()
                    if not scope or not table or not any(a.GetAtomMapNum() for a in mol.GetAtoms()):
                        issues.append(f'{path.name} 第 {number} 条需要 PDF 或 DOI、Table 和带论文编号的原子映射')
                        continue
                    structures[scope, table].append({
                        'name': mol.GetProp('_Name') if mol.HasProp('_Name') else f'Table {table}',
                        'paper_atom_mapped_smiles': Chem.MolToSmiles(mol),
                        'source': {'mapping_file': path.name, 'mapping_record': number},
                    })
        except (OSError, ValueError) as exc:
            issues.append(f'{path.name} · {exc}')
    return structures, issues


def number(label):
    return int(re.match(r'[HC](\d+)', label)[1])


def molecular_mapping(raw, expected):
    mapped = Chem.MolFromSmiles(raw)
    if mapped is None:
        raise ValueError('无法解析编号结构')
    paper = {}
    for atom in mapped.GetAtoms():
        tag = atom.GetAtomMapNum()
        if tag:
            if tag in paper:
                raise ValueError(f'结构含重复的论文编号 {tag}')
            paper[tag] = atom.GetIdx()
        atom.SetAtomMapNum(0)
    smiles, processed, old_to_new = db.canonicalize(mapped)
    key = Chem.MolToInchiKey(processed)
    if expected and key != expected:
        raise ValueError('结构与已核对的 InChIKey 不一致')
    return smiles, processed, {tag: old_to_new[idx] for tag, idx in paper.items() if idx in old_to_new}


def diagram(mol, mapping, special_labels):
    rdDepictor.Compute2DCoords(mol)
    drawing = rdMolDraw2D.MolDraw2DSVG(580, 380)
    options = drawing.drawOptions()
    options.padding = .08
    for tag, idx in mapping.items():
        mol.GetAtomWithIdx(idx).SetProp('atomNote', special_labels.get(tag, str(tag)))
    drawing.DrawMolecule(mol)
    points = {tag: [drawing.GetDrawCoords(idx).x, drawing.GetDrawCoords(idx).y] for tag, idx in mapping.items()}
    drawing.FinishDrawing()
    content = drawing.GetDrawingText()
    return content[content.index('<svg '):], points


def check_assignments(table, mapping, mol, overrides):
    for nucleus, peaks in table['peaks'].items():
        for peak in peaks:
            peak['paper_atom_numbers'] = []
            for label in peak['source_labels']:
                tag = overrides.get(label, number(label))
                if label.endswith('-NH') and label not in overrides:
                    carbon_index = mapping.get(number(label))
                    neighbors = [] if carbon_index is None else [a for a in mol.GetAtomWithIdx(carbon_index).GetNeighbors()
                        if a.GetSymbol() == 'N' and a.GetTotalNumHs(includeNeighbors=True) > 0]
                    if len(neighbors) != 1:
                        raise ValueError(f'{label} 未能唯一对应邻接的带氢氮原子')
                    index = neighbors[0].GetIdx()
                    tag = next((k for k, idx in mapping.items() if idx == index), max(mapping, default=0) + 1)
                    mapping[tag] = index
                    overrides[label] = tag
                index = mapping.get(tag)
                if index is None:
                    raise ValueError(f'{label} 在编号结构中没有对应原子')
                atom = mol.GetAtomWithIdx(index)
                if (nucleus == '13C' and atom.GetSymbol() != 'C' or label.endswith('-NH') and atom.GetSymbol() != 'N' or
                        nucleus == '1H' and atom.GetTotalNumHs(includeNeighbors=True) == 0):
                    raise ValueError(f'{label} 与结构中的原子类型不符')
                peak['paper_atom_numbers'].append(tag)


def preserve_numbering_uncertainty(table, site_sets):
    for peaks in table['peaks'].values():
        for tags in site_sets:
            group = [peak for peak in peaks if set(peak['paper_atom_numbers']).intersection(tags)]
            if len(group) == len(tags) and all(not peak.get('group') for peak in group):
                for peak in group:
                    peak.update(relation='interchangeable', group='structure-' + '-'.join(map(str, tags)),
                                group_reason='论文保留各编号，规范分子图无法唯一映射这些甲基的编号')


def prepare_table(table, source, profile, sdf, supplement, *, draw=True):
    table.update(pdf=source['name'], file_id=source['id'], doi=source.get('doi'), sha256=source['sha256'])
    entry = profile or sdf
    table['name'] = entry['name'] if entry else table['caption'] or f"表 {table['table'] or '?'}"
    if not entry:
        table['issues'].append('需要带论文原子编号的结构 SDF，已保留原文峰表')
        return
    table['metadata'] = {**entry['source'], 'pdf_sha256': source['sha256'], 'pdf': source['name'],
                         'pdf_page': table['page'], 'table': table['table'], 'doi': table['doi']}
    table['conditions'] = entry.get('conditions', {})
    table['unreported_proton_sites'] = entry.get('unreported_proton_sites', [])
    if profile and supplement is not None and 'supplement_pdf_page' in profile:
        table.update(supplement_file_id=supplement['id'], supplement_page=profile['supplement_pdf_page'])
    if profile and 'spectrum_page' in profile:
        table.update(supplement_file_id=source['id'], supplement_page=profile['spectrum_page'])
        if {nucleus: len(peaks) for nucleus, peaks in table['peaks'].items()} != profile['table_rows']:
            table['issues'].append('对照表提取数量与已核对行数不一致，请检查续页')
    try:
        smiles, mol, mapping = molecular_mapping(entry['paper_atom_mapped_smiles'], entry.get('computed_inchikey'))
        overrides = dict(entry.get('label_overrides', {}))
        check_assignments(table, mapping, mol, overrides)
        preserve_numbering_uncertainty(table, entry.get('permutable_structure_sites', []))
        display = {tag: label.removeprefix('H') for label, tag in overrides.items()}
        table['smiles'] = smiles
        table['formula'] = rdMolDescriptors.CalcMolFormula(mol)
        table['mapping'] = mapping
        table['atom_display'] = display
        if draw:
            table['diagram'], table['points'] = diagram(mol, mapping, display)
    except ValueError as exc:
        table['issues'].append(str(exc))


def spectrum_payload(table, nucleus, rows):
    mapping = table['mapping']
    provenance = {**table['metadata'], 'experiment': 'experimental',
                  'paper_atom_number_to_atom_index': mapping, 'footnotes_text': table['footnotes_text'],
                  'unreported_proton_sites': table['unreported_proton_sites'] if nucleus == '1H' else []}
    if 'component_pages' in table:
        provenance['pdf_pages'] = table['component_pages'].get(nucleus, [])
    fields = table['conditions']
    conditions = {'solvent': fields.get('solvent'), 'frequency_mhz': fields.get('frequencies_mhz', {}).get(nucleus),
                  'temperature_k': fields.get('temperature_k'), 'instrument': fields.get('instrument'),
                  'chemical_shift_reference': fields.get('chemical_shift_reference'),
                  'frequency_source': fields.get('frequency_source')}
    peaks, assignments, grouped = [], {}, defaultdict(list)
    for row in rows:
        peak_id = f"{nucleus[-1]}{row['row']:03}"
        bounds = row['range_ppm']
        targets = sorted({mapping[tag] for tag in row['paper_atom_numbers']})
        entry = {'peak_id': peak_id, 'shift_ppm': sum(bounds)/2 if bounds else row['shift_ppm'],
                 'multiplicity': row['multiplicity'], 'reported_nH': row['reported_nH'],
                 'j_values': row['j_values_hz'], 'literature': row}
        if bounds:
            entry.update(range_min=bounds[0], range_max=bounds[1], shift_origin='midpoint_of_reported_range')
        peaks.append(entry)
        if row.get('group'):
            grouped[row['group']].append((peak_id, targets, row.get('group_reason', '论文表或脚注的可互换归属')))
            assignments[peak_id] = {'status': 'ambiguous', 'atom_indices': [], 'note': '论文归属组，组内对应未定'}
        elif targets:
            status = 'assigned' if row['relation'] == 'single_site' else 'source_group'
            assignments[peak_id] = {'status': status, 'atom_indices': targets, 'note': row['source_assignment_text']}
    source_key = f"{table['sha256']}:{table['page']}:{table['table']}:{table['ordinal']}:{nucleus}"
    groups = [{'id': db.uid('G', source_key, key), 'peak_ids': [p for p, _, _ in entries],
               'atom_indices': sorted({i for _, atoms, _ in entries for i in atoms}),
               'relation': 'unordered_set_correspondence', 'within_group_order': 'unknown',
               'chemical_equivalence_asserted': False, 'reason': entries[0][2]}
              for key, entries in grouped.items()]
    provenance[REVISIONS] = {'signal_codes': {}, 'assignment_groups': groups}
    return source_key, conditions, peaks, provenance, assignments


def save_table(conn, table):
    mid = db.add_molecule(conn, table['smiles'], SOURCE)
    ids = {}
    for nucleus, rows in table['peaks'].items():
        if not rows:
            continue
        ids[nucleus] = db.add_spectrum(conn, mid, nucleus, SOURCE, *spectrum_payload(table, nucleus, rows))
    if len(ids) == 2:
        db.add_pair(conn, mid, ids['13C'], ids['1H'], 'literature_same_table')
    db.refresh_molecule(conn, mid)
    table['molecule_id'] = mid
    table['spectrum_ids'] = ids


def save_results(output, tables, report, *, portable_tables=None):
    export_tables = tables if portable_tables is None else portable_tables
    output.mkdir(parents=True, exist_ok=False)
    if any(table.get('smiles') and not table['issues'] for table in export_tables):
        pending = output / 'assignments.tmp.sqlite'
        with closing(db.connect(pending)) as conn:
            db.initialize(conn)
            ai_database.initialize(conn)
            with conn:
                for table in export_tables:
                    if table.get('smiles') and not table['issues']:
                        save_table(conn, table)
            if conn.execute('PRAGMA foreign_key_check').fetchone() or conn.execute('PRAGMA quick_check').fetchone()[0] != 'ok':
                raise RuntimeError('SQLite 校验失败')
            conn.execute('PRAGMA wal_checkpoint(TRUNCATE)')
        pending.replace(output / 'assignments.sqlite')
        report['database'] = str(output / 'assignments.sqlite')
    report['saved'] = len({table['molecule_id'] for table in export_tables if 'molecule_id' in table})
    report['saved_spectra'] = sum(len(t.get('spectrum_ids', {})) for t in export_tables)
    report['output_folder'] = str(output)
    payload = {**report, 'schema': db.SCHEMA,
               'tables': [{k: v for k, v in table.items() if k not in ('diagram', 'points')} for table in tables]}
    (output / 'extraction_report.json').write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf-8')
    return report


def collect(folder, progress, *, checksums=None, read_document=read_pdf, draw=True):
    paths = source_paths(folder)
    pdfs = [p for p in paths if p.suffix.lower() == '.pdf']
    if not pdfs:
        raise ValueError('所选文件夹里没有 PDF。请选包含论文 PDF 的这一层文件夹。')
    structures, issues = read_sdf(p for p in paths if p.suffix.lower() == '.sdf')
    files = [{'id': i, 'name': p.name, 'sha256': checksums[p.name] if checksums else digest(p)} for i, p in enumerate(pdfs)]
    supplement = next((file for file in files if file['sha256'] == SAMPLE['supplement_sha256']), None)
    tables = []
    seen = set()
    for i, (path, file) in enumerate(zip(pdfs, files), 1):
        # When an SDF binds a PDF by name, each filename may carry different assignments.
        if file['sha256'] in seen and not structures:
            file['status'] = '内容重复，已跳过'
            continue
        seen.add(file['sha256'])
        progress(i, len(pdfs), file['name'])
        if file is supplement:
            file['status'] = '补充谱图'
            continue
        known = file['sha256'] == SAMPLE['pdf_sha256']
        comparison_si = file['sha256'] == COMPARISON_SAMPLE['supplement_sha256']
        comparison_main = file['sha256'] == COMPARISON_SAMPLE['pdf_sha256']
        profiles = {str(row['pdf_page']): row for row in SAMPLE['compounds']} if known else {}
        try:
            document = read_document(path, file['sha256'], profiles)
        except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
            file.update(status='读取失败', issue=str(exc))
            continue
        file.update(doi=SAMPLE['doi'] if known else COMPARISON_SAMPLE['doi'] if comparison_si else document['doi'],
                    pages=document['pages'])
        file['status'] = ('找到归属表' if document['tables'] else
                          '正文实验记录（归属表在补充材料）' if comparison_main else
                          '需 OCR' if not document['text_pages'] else '待复核版式')
        for index, table in enumerate(document['tables'], 1):
            table['ordinal'] = index
            if table.get('comparison') and not table['table']:
                table['table'] = table['caption'].upper()
            candidates = structures.get((file['name'], table['table']), []) or structures.get((file.get('doi'), table['table']), [])
            profile = profiles.get(str(table['page']))
            if profile and str(profile['table']) != table['table']:
                table['issues'].append('提取表号与已核对的原表不一致')
                profile = None
            if comparison_si and table.get('comparison'):
                profile = next((entry for entry in COMPARISON_SAMPLE['compounds'] if entry['caption'] == table['caption']), None)
            if not profile and len(candidates) > 1:
                table['issues'].append('多个 SDF 结构匹配到同一论文表')
            prepare_table(table, file, profile, candidates[0] if len(candidates) == 1 else None, supplement, draw=draw)
            tables.append(table)
    return {'folder': str(folder), 'files': files, 'warnings': issues, 'tables': tables}


def run(folder, output, progress):
    report = collect(folder, progress)
    tables = report['tables']
    report.update(pdf_count=len(report['files']), peak_count=sum(len(v) for t in tables for v in t['peaks'].values()))
    report = save_results(output, tables, report)
    return {**report, 'tables': tables}
