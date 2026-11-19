"""Attach the reviewed IBM H/C batch to AtomNMR without editing experimental spectra."""

import argparse
from collections import Counter
from contextlib import closing
import hashlib
import json
from pathlib import Path
import sqlite3
import sys

import numpy as np
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

from audit_ibm_computational_nmr import MD5, OUT, SOURCE, frame_labels, geometry_atom_order

sys.path.insert(0, str(OUT.parents[1]))
from annotation_tool_v3 import computed_spectra, core
from external_spectrum_atoms import target_mappings

DOI = '10.5281/zenodo.16417648'
NUCLEI = {'c': '13C', 'h': '1H'}
# The IBM batch audit confirmed that this shared standard InChIKey represents
# different localized H/bond orders. Keep the relation as reference only.
TAUTOMER_VARIANT = 1369


def read_jsonl(name):
    with (OUT / name).open(encoding='utf-8') as f:
        return [json.loads(line) for line in f if line.strip()]


def requested_batch():
    targets = {r['source_id']: r for r in read_jsonl('ibm_novel_dual_candidates.jsonl')}
    known = [r for r in read_jsonl('ibm_exact_overlap_assignments.jsonl')
             if r['main_assignment_coverage'] in {'unassigned', 'no_spectrum'}]
    all_rows = {r['source_id']: r for r in read_jsonl('ibm_all_records_audit.jsonl')}
    for r in known:
        sid = r['source_id']
        if sid in targets:
            raise ValueError(f'IBM record {sid} appears in both selection lists')
        targets[sid] = {**all_rows[sid], 'expected_id': r['main_molecule_id'],
                        'requested_group': r['main_assignment_coverage']}
    return targets


def resolve(c, identity):
    key = identity['inchikey']
    exact = c.execute('SELECT id,smiles,structure_json,c_count,h_count FROM molecules WHERE inchikey=?',
                      (key,)).fetchall()
    expected_id = identity.get('expected_id')
    if len(exact) > 1:
        raise ValueError(f'IBM {identity["source_id"]} has more than one exact AtomNMR target')
    if exact:
        if expected_id and exact[0]['id'] != expected_id:
            raise ValueError(f'IBM {identity["source_id"]} now maps to a different AtomNMR ID')
        return exact[0]
    if expected_id:
        raise ValueError(f'IBM {identity["source_id"]} has lost its exact AtomNMR match')
    prefix = key[:14]
    if c.execute('SELECT 1 FROM molecules WHERE inchikey>=? AND inchikey<? LIMIT 1',
                 (prefix+'-', prefix+'.')).fetchone():
        raise ValueError(f'IBM {identity["source_id"]} now has a connectivity-only AtomNMR match')
    return None


def source_order(row, identity):
    if identity.get('source_order_molblock'):
        mol = Chem.MolFromMolBlock(identity['source_order_molblock'], removeHs=False)
        if mol is None or [a.GetSymbol() for a in mol.GetAtoms()] != row['atoms']:
            raise ValueError('Validated IBM atom order could not be loaded')
        return 'validated_audit_geometry', mol
    template = Chem.MolFromSmiles(row['smiles'])
    status, ordered = geometry_atom_order(row, template, allow_unspecified_stereo=True)
    if ordered is not None:
        return status, ordered
    for frame in row['frames'].values():
        if frame is None or not frame.get('xyz'):
            continue
        status, ordered = geometry_atom_order({**row, 'xyz': frame['xyz']}, template,
                                              allow_unspecified_stereo=True)
        if ordered is not None:
            return f'{status}_from_trajectory', ordered
    return status, None


def make_peaks(n, atom_shifts, groups, mappings, attached):
    by_atom = {p['atom_index']: p for p in atom_shifts}
    rows = []
    for group in groups:
        source_shifts = [by_atom[i] for i in group]
        positions = [
            tuple(sorted({mapping[attached[i] if n == 'h' else i] for i in group}))
            for mapping in mappings
        ]
        targets = sorted({i for position in positions for i in position})
        status = ('unassigned' if not targets else 'ambiguous' if len(set(positions)) > 1
                  else 'assigned' if len(targets) == 1 else 'equivalent_group')
        rows.append({
            'shift_ppm': round(float(np.mean([s['ppm'] for s in source_shifts])), 4),
            'source_atom_ids': group, 'source_atom_shifts': source_shifts,
            'assignment': {'status': status, 'atom_indices': targets, 'origin': computed_spectra.IBM, 'note': ''},
        })
    rows.sort(key=lambda r: (-r['shift_ppm'], r['source_atom_ids']))
    return [{**row, 'peak_id': f'{n.upper()}{i}'} for i, row in enumerate(rows, 1)]


def prepare(c, identities):
    plans = []
    for batch in pq.ParquetFile(SOURCE).iter_batches(batch_size=32):
        for raw in batch.to_pylist():
            sid = raw['id']
            if sid not in identities:
                continue
            identity = identities[sid]
            source_graph = Chem.MolFromSmiles(raw['smiles'])
            if Chem.MolToInchiKey(source_graph) != identity['inchikey']:
                raise ValueError(f'IBM {sid} differs from the audited structure')
            matched = resolve(c, identity)
            smiles = matched['smiles'] if matched else identity['canonical_smiles']
            target = Chem.MolFromSmiles(smiles)
            if Chem.MolToInchiKey(target) != identity['inchikey']:
                raise ValueError(f'IBM {sid} and the current structure do not share a complete InChIKey')
            same_graph = Chem.MolToSmiles(source_graph) == Chem.MolToSmiles(target)
            if not same_graph and sid != TAUTOMER_VARIANT:
                raise ValueError(f'Unexpected molecular graph variant for IBM {sid}')
            structure = json.loads(matched['structure_json']) if matched else core.structure(target)
            status, source = source_order(raw, identity)
            mappings, attached = target_mappings(source, target, structure) if same_graph else ([], {})
            labels, groups_valid, _, groups = frame_labels(raw)
            peaks = {n: make_peaks(n, labels[n], groups[n] if groups_valid[n] else
                                  [[p['atom_index']] for p in labels[n]], mappings, attached) for n in NUCLEI}
            note = ('已对齐当前结构的原子编号；对称结构中的多种有效对应显示为歧义。' if mappings
                    else '已保存 IBM 的位移与源原子编号；与当前结构的编号对应待核对。')
            if 'declared_stereo_only' in status:
                note += '源结构未指定所有立体位置，计算构象的立体信息仅供参考。'
            if not same_graph:
                note = 'IBM 与主库记录为不同互变异构体，仅作为同一 InChIKey 下的对照谱。主图显示主库结构，IBM 谱保留其原始原子编号。'
            plans.append({
                'identity': identity, 'smiles': smiles, 'planned_id': matched['id'] if matched else None,
                'structure': structure, 'peaks': peaks,
                'provenance': {'source_record_id': sid, 'source_url': f'https://doi.org/{DOI}',
                               'source_md5': MD5, 'source_smiles': raw['smiles'],
                               'structure_match': 'exact_graph' if same_graph else 'tautomer_variant',
                               'spectrum_kind': 'computed', 'frame_count': raw['number_of_frames'],
                               'frame_ids': raw['frame_ids'], 'chemical_shift_reference': 'TMS',
                               'atom_mapping_status': status if mappings else 'source_indices_only',
                               'atom_mapping_note': note, 'source_atom_index_basis': 'zero_based_IBM_XYZ_order',
                               'source_atom_symbols': raw['atoms']},
            })
    if len(plans) != len(identities):
        raise ValueError('Some requested IBM records were missing')
    return plans


def apply(c, plans):
    counts = Counter()
    c.execute('BEGIN IMMEDIATE')
    try:
        computed_spectra.initialize(c)
        for plan in plans:
            identity = plan['identity']
            target = resolve(c, identity)
            if (target and (target['smiles'] != plan['smiles']
                            or json.loads(target['structure_json']) != plan['structure'])) or (
                    plan['planned_id'] and (not target or target['id'] != plan['planned_id'])):
                raise ValueError('An AtomNMR structure changed during preparation; please rerun the import')
            if target:
                mid = target['id']
                counts['matched_existing'] += 1
            else:
                mid = core.add_molecule(c, plan['smiles'], computed_spectra.IBM)
                counts['new_molecules'] += 1
            source_id = str(identity['source_id'])
            retained = 0
            for n, nucleus in NUCLEI.items():
                existing = c.execute('SELECT id,molecule_id,peaks_json,source_json FROM computed_spectra WHERE source=? AND source_record_id=? AND nucleus=?',
                                     (computed_spectra.IBM, source_id, nucleus)).fetchone()
                peaks_json = core.dumps(plan['peaks'][n])
                source_json = core.dumps(plan['provenance'])
                if existing:
                    if existing['molecule_id'] != mid:
                        raise ValueError(f'IBM record {source_id} already belongs to a different molecule')
                    if (existing['peaks_json'], existing['source_json']) != (peaks_json, source_json):
                        c.execute('UPDATE computed_spectra SET peaks_json=?,source_json=? WHERE id=?',
                                  (peaks_json, source_json, existing['id']))
                        counts['updated_existing_spectra'] += 1
                    retained += 1
                    continue
                sid = core.uid('S', computed_spectra.IBM, source_id, nucleus)
                c.execute('INSERT INTO computed_spectra VALUES(?,?,?,?,?,?,?)',
                          (sid, mid, nucleus, computed_spectra.IBM, source_id, peaks_json, source_json))
            if retained != len(NUCLEI):
                # Preserve the experimental counts and status of matched records.
                # Only a genuinely experiment-free record needs a new queue label.
                if not target or not (target['c_count'] or target['h_count']):
                    core.refresh_molecule(c, mid)
            counts['inserted_spectra'] += len(NUCLEI)-retained
            counts['existing_spectra'] += retained
        c.commit()
    except Exception:
        c.rollback()
        raise
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', type=Path, default=core.DEFAULT_DB)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    if hashlib.md5(SOURCE.read_bytes()).hexdigest() != MD5:
        raise ValueError('IBM Parquet differs from the downloaded, verified version')
    database = core.main_database(args.db).resolve()
    with closing(sqlite3.connect(database.as_uri()+'?mode=ro', uri=True)) as c:
        c.row_factory = sqlite3.Row
        c.execute('BEGIN')
        plans = prepare(c, requested_batch())
    result = {
        'records': len(plans),
        'requested_groups': dict(Counter(p['identity'].get('requested_group', 'new') for p in plans)),
        'atom_mapping': dict(Counter(p['provenance']['atom_mapping_status'] for p in plans)),
        'source_indices_only_records': [p['identity']['source_id'] for p in plans
                                        if p['provenance']['atom_mapping_status']=='source_indices_only'],
        'tautomer_variant_records': [p['identity']['source_id'] for p in plans
                                     if p['provenance']['structure_match']=='tautomer_variant'],
    }
    if args.apply:
        with closing(core.connect(database)) as c:
            result['database'] = dict(apply(c, plans))
        (OUT / 'atomnmr_import_summary.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
