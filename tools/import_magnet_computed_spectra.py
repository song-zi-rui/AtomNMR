"""Import verified MagNET H/C structures and DFT stick spectra into AtomNMR."""

import argparse
from collections import Counter, defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack, closing
from pathlib import Path
import json
import math
import sqlite3

import h5py
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

from scan_magnet_structures import (
    METHODS, OUT, SOURCE, declared_geometry, file_for, recover,
)
from annotation_tool_v3 import computed_spectra, core
from annotation_tool_v3.export_policy import coverage, current_spectra
from external_spectrum_atoms import target_mappings

DATA_ORDER = ('concentrate-pcm', 'pepper-reference', 'pepper-pbe0', 'shake')
SHIELDINGS = {
    'concentrate-pcm': ('b3lyp_pcm', 'shielding_b3lyp_pcm', 'shielding_b3lyp_pcm'),
    'pepper-reference': ('pepper_gas', 'shielding_wp04', 'shielding_wb97xd'),
    'shake': ('pbe0_gas', 'shielding_stationary', 'shielding_stationary'),
    'pepper-pbe0': ('pbe0_gas', 'shielding_pbe0', 'shielding_pbe0'),
}
NUCLEI = ('1H', '13C')
IMPORTABLE = {'new', 'unassigned', 'no_spectrum'}
GROUP_SPREAD = {'1H': 0.03, '13C': 0.25}


def key_for(row):
    return f"{row['method']}:{row['source_id']}"


def classify(c, row, minimum):
    inchikey, smiles = row['inchikey'], row['smiles']
    layer = inchikey[:14]
    matches = c.execute(
        'SELECT id,inchikey,smiles,structure_json,revision FROM molecules '
        'WHERE inchikey>=? AND inchikey<?', (layer+'-', layer+'.')
    ).fetchall()
    if not matches:
        # Some older records lack an InChIKey. Check their graph separately.
        for other in c.execute('SELECT smiles FROM molecules WHERE inchikey IS NULL '
                               'AND formula=?', (row['formula'],)):
            m = Chem.MolFromSmiles(other['smiles'])
            if m is not None and Chem.MolToInchiKey(m)[:14] == layer:
                return 'unkeyed_possible_match', None
        if row['heavy_atoms'] < minimum:
            return 'novel_below_minimum_heavy_atoms', None
        return 'new', None
    targets = [m for m in matches if m['inchikey'] == inchikey and
               (m['smiles'] == smiles or Chem.MolToSmiles(Chem.MolFromSmiles(m['smiles'])) == smiles)]
    if len(targets) != 1:
        return ('ambiguous_exact' if targets else 'connectivity_or_tautomer_only'), None
    target = targets[0]
    if c.execute('SELECT 1 FROM computed_spectra WHERE molecule_id=? AND source=? LIMIT 1',
                 (target['id'], computed_spectra.MAGNET)).fetchone():
        return 'already_has_magnet', target
    record = core.record(c, target['id'], require_quality=True, include_computed=True)
    usable = [s for s in current_spectra(record) if s['peaks']]
    status = ('inactive_molecule' if record['curation']['molecule_state'] != 'active' else
              'no_spectrum' if not usable else 'unassigned' if coverage(usable) == 'none' else
              'has_valid_assignment')
    return status, target


class ShieldingFile:
    def __init__(self, method, files):
        self.method = method
        path = file_for(method)
        self.group = files[path][METHODS[method][1]]
        g = self.group
        self.coord_name = 'stationary_coords' if method == 'shake' else 'coordinates'
        _, *self.shielding_names = SHIELDINGS[method]
        if method == 'shake':
            self.start, self.end = g['atom_start'][:], g['atom_end'][:]
        else:
            offsets = np.r_[0, np.cumsum(g['n_atoms'][:], dtype=np.int64)]
            self.start, self.end = offsets[:-1], offsets[1:]
        self.scales = {name: float(g[name].attrs['scale'])
                       for name in {self.coord_name, *self.shielding_names}}

    def load(self, row_numbers):
        # Read adjacent molecules once to avoid inflating the same HDF5 chunk
        # separately for every atom or every candidate.
        origin, end = self.start[row_numbers[0]], self.end[row_numbers[-1]]
        span = slice(origin, end)
        z = self.group['atomic_numbers'][span]
        coords = self.group[self.coord_name][span] * self.scales[self.coord_name]
        fields = {name: self.group[name][span] * self.scales[name] for name in set(self.shielding_names)}
        for i in row_numbers:
            sl = slice(self.start[i]-origin, self.end[i]-origin)
            yield z[sl], coords[sl], [fields[name][sl] for name in self.shielding_names]


def make_peaks(z, values, nucleus, mappings, attached, calibration, orbits):
    atomic_number = 1 if nucleus == '1H' else 6
    label = 'H' if atomic_number == 1 else 'C'
    lo, hi = (-20, 40) if atomic_number == 1 else (-50, 350)
    by_orbit = defaultdict(list)
    for i in np.flatnonzero(z == atomic_number):
        shielding = round(float(values[i]), 4)
        delta = calibration['intercept'] + calibration['slope'] * shielding
        if not math.isfinite(delta) or not lo <= delta <= hi or shielding == 0:
            return None
        idx = attached[int(i)] if atomic_number == 1 else int(i)
        target_orbit = orbits[idx]
        by_orbit[target_orbit].append((delta, int(i), idx, shielding))
    peaks = []
    for orbit, members in by_orbit.items():
        members.sort()
        groups = []
        for item in members:
            if not groups or item[0] - groups[-1][0][0] > GROUP_SPREAD[nucleus]:
                groups.append([])
            groups[-1].append(item)
        for group in groups:
            positions = {tuple(sorted({m[item[2]] for item in group})) for m in mappings}
            targets = list(orbit)
            status = ('ambiguous' if positions != {orbit} else
                      'equivalent_group' if len(targets) > 1 else 'assigned')
            values = [{'atom_index': item[1], 'ppm': round(item[0], 4), 'shielding_ppm': item[3]}
                      for item in sorted(group, key=lambda item: item[1])]
            peak = {
                'shift_ppm': round(sum(item[0] for item in group) / len(group), 4),
                'source_atom_ids': [item['atom_index'] for item in values],
                'source_atom_shifts': values,
                'assignment': {'status': status, 'atom_indices': targets,
                               'origin': computed_spectra.MAGNET, 'note': ''},
            }
            if atomic_number == 1:
                peak['reported_nH'] = len(group)
            peaks.append(peak)
    peaks.sort(key=lambda p: (-p['shift_ppm'], p['source_atom_ids']))
    return [{**peak, 'peak_id': f'{label}{i}'} for i, peak in enumerate(peaks, 1)]


def novel_mapping(source, heavy, z, structure):
    # Canonical SMILES was just generated from this geometry graph. RDKit
    # supplies the exact heavy-atom traversal order used to write the SMILES.
    order = json.loads(heavy.GetProp('_smilesAtomOutputOrder'))
    heavy_indices = np.flatnonzero(z != 1).tolist()
    if len(order) != len(heavy_indices) or len(order) != len(structure['atoms']):
        return [], {}, {}
    mapping = {heavy_indices[i]: j for j, i in enumerate(order)}
    attached = {a.GetIdx(): a.GetNeighbors()[0].GetIdx() for a in source.GetAtoms()
                if a.GetAtomicNum() == 1 and a.GetDegree() == 1}
    equivalent = {i: tuple(sorted(g)) for g in structure['equivalent_groups'] for i in g}
    orbits = {i: equivalent.get(j, (j,)) for i, j in mapping.items()}
    return [mapping], attached, orbits


def prepare(row, data, target, models, revision, hash_by_dataset):
    method = row['method']
    z, coords, shieldings = data
    if method.startswith('pepper'):
        recovered = declared_geometry(z, coords, row['smiles'])
    else:
        recovered = recover(z, coords, alternate=True)
    if not recovered or recovered[2] != row['smiles']:
        return 'geometry_graph_or_stereo_unconfirmed', None
    source, heavy, _ = recovered
    mol = Chem.MolFromSmiles(target['smiles'] if target else row['smiles'])
    structure = json.loads(target['structure_json']) if target else core.structure(mol)
    if target:
        mapping, attached = target_mappings(source, mol, structure)
        orbits = {i: tuple(sorted({m[i] for m in mapping})) for i in mapping[0]} if mapping else {}
    else:
        mapping, attached, orbits = novel_mapping(source, heavy, z, structure)
    if not mapping or len(attached) != int((z == 1).sum()) or any(i not in orbits for i in attached.values()):
        return 'unresolved_atom_mapping', None
    model = SHIELDINGS[method][0]
    cal = models[model]
    peaks = {n: make_peaks(z, values, n, mapping, attached, cal[n], orbits)
             for n, values in zip(NUCLEI, shieldings)}
    if not all(peaks.values()):
        return 'invalid_shielding_or_shift', None
    if any(not any(p['assignment']['status'] != 'ambiguous' for p in pp) for pp in peaks.values()):
        return 'no_certain_assignment_in_one_nucleus', None
    dataset = METHODS[method][0]
    source_url = f'https://huggingface.co/ekwan16/MagNET/tree/{revision}/data/{dataset}'
    provenance = {
        'source_dataset': dataset, 'source_record_id': key_for(row), 'source_row_zero_based': row['source_row'],
        'source_molecule_id': row['source_id'], 'sigma_shake_id': row['shake_id'],
        'source_url': source_url, 'source_sha256': hash_by_dataset[dataset],
        'source_smiles': row['smiles'], 'structure_match': 'exact_graph', 'spectrum_kind': 'computed',
        'label_origin': 'per_atom_DFT_shielding', 'source_atom_index_basis': 'zero_based_hdf5_molecule_atom_order',
        'source_atomic_numbers': z.tolist(),
        'atom_mapping_status': 'validated_declared_smiles' if method.startswith('pepper') else 'validated_geometry_reconstruction',
        'atom_mapping_note': '保留源文件逐原子计算值。氢谱在结构图上指向所连重原子，近似重合的同类位移可合并显示。',
        'calibration_id': 'magnet-delta22-chloroform-v1', 'calibration_model': model,
        'calculation_environment': 'PCM chloroform' if method == 'concentrate-pcm' else 'gas',
        'chemical_shift_reference': 'Delta22 实验氯仿标定',
        'shift_model': {n: {k: cal[n][k] for k in ('intercept', 'slope')} for n in NUCLEI},
        'stationary_geometry_only': True,
    }
    return 'ready', (Descriptors.MolWt(mol), core.dumps(structure),
                     {n: core.dumps(p) for n, p in peaks.items()}, provenance)


def attach(c, row, target, plan, status):
    mw, structure_json, peaks, provenance = plan
    mid = target['id'] if target else core.uid('M', row['smiles'])
    if not target:
        c.execute('INSERT INTO molecules(id,smiles,inchikey,connectivity,formula,mw,structure_json,origin,created_at,queue) '
                  "VALUES(?,?,?,?,?,?,?,?,?,'simulation_only')",
                  (mid, row['smiles'], row['inchikey'], row['inchikey'][:14], row['formula'],
                   mw, structure_json, computed_spectra.MAGNET, core.now()))
    provenance['import_category'] = status
    cited = core.dumps(provenance)
    for nucleus, values in peaks.items():
        c.execute('INSERT INTO computed_spectra VALUES(?,?,?,?,?,?,?)',
                  (core.uid('S', computed_spectra.MAGNET, key_for(row), nucleus), mid, nucleus,
                   computed_spectra.MAGNET, key_for(row), values, cited))
    if status == 'no_spectrum':
        core.refresh_molecule(c, mid)
    return mid


def prior_import(c, row):
    stored = c.execute('SELECT molecule_id,source_json,nucleus FROM computed_spectra '
                       'WHERE source=? AND source_record_id=?', (computed_spectra.MAGNET, key_for(row))).fetchall()
    if not stored:
        return None
    if {s['nucleus'] for s in stored} != set(NUCLEI) or len({s['molecule_id'] for s in stored}) != 1:
        raise ValueError(f'Previous partial import for {key_for(row)}')
    return json.loads(stored[0]['source_json'])['import_category'], stored[0]['molecule_id']


def groups(stage, method, *, limit, overlaps_only, novel_heavy, new_quota):
    filter_sql, params = '', ()
    if overlaps_only:
        filter_sql = 'AND original_layer=1'
    elif novel_heavy is not None:
        filter_sql = 'AND original_layer=0 AND heavy_atoms=?'
        params = (novel_heavy,)
    base = 'candidates'
    if new_quota is not None:
        # InChIKey ordering samples a stratum independently of the source's
        # molecule order; geometry/shieldings are subsequently read in row order.
        stage.execute('DROP TABLE IF EXISTS temp.magnet_panel')
        stage.execute("""
            CREATE TEMP TABLE magnet_panel AS SELECT b.* FROM candidates b
            WHERE b.method=? AND b.original_layer=0 AND b.heavy_atoms=?
            AND NOT EXISTS (SELECT 1 FROM import_result r
                            WHERE r.method=b.method AND r.source_row=b.source_row)
            ORDER BY b.inchikey,b.source_row LIMIT ?
        """, (method, novel_heavy, int(new_quota * 1.7) + 64))
        stage.execute('CREATE INDEX temp.magnet_panel_row ON magnet_panel(source_row)')
        base = 'temp.magnet_panel'
    sql = f'SELECT * FROM {base} WHERE method=? AND source_row>? {filter_sql} ORDER BY source_row LIMIT ?'
    offset, remaining = -1, limit
    while remaining is None or remaining > 0:
        batch = stage.execute(sql, (method, offset, *params,
                                    600 if remaining is None else min(600, remaining))).fetchall()
        if not batch:
            break
        offset = batch[-1]['source_row']
        if remaining is not None:
            remaining -= len(batch)
        start = 0
        while start < len(batch):
            end = start + 1
            block = batch[start]['source_row'] // 8192
            while end < len(batch) and batch[end]['source_row'] // 8192 == block:
                end += 1
            yield batch[start:end]
            start = end


def init_worker(models, revision, hashes):
    global CALIBRATION, REVISION, HASHES, FILES, READERS
    CALIBRATION, REVISION, HASHES = models, revision, hashes
    FILES, READERS = {}, {}
    RDLogger.DisableLog('rdApp.*')


def prepare_batch(method, selected):
    if method not in READERS:
        path = file_for(method)
        if path not in FILES:
            FILES[path] = h5py.File(path)
        READERS[method] = ShieldingFile(method, FILES)
    data = READERS[method].load([x[0]['source_row'] for x in selected])
    return [prepare(r, source, target, CALIBRATION, REVISION, HASHES)
            for (r, _, target), source in zip(selected, data, strict=True)]


def commit_batch(stage, con, method, decisions, selected, future, report, minimum, apply, new_quota):
    result = future.result() if future else []
    plans = []
    for (row, category, target), (state, plan) in zip(selected, result, strict=True):
        if plan:
            plans.append((row, category, target, plan))
        else:
            decisions.append((method, row['source_row'], state, target['id'] if target else None))
    if apply and plans:
        con.execute('BEGIN IMMEDIATE')
        try:
            already_new = report[(method, 'new')]
            for row, category, target, plan in plans:
                if category == 'new' and new_quota is not None and already_new >= new_quota:
                    continue
                check, latest = classify(con, row, minimum)
                same_target = target is None or (latest and
                    all(target[k] == latest[k] for k in ('id', 'revision', 'structure_json')))
                if check == category and same_target:
                    mid = attach(con, row, target, plan, category)
                    decisions.append((method, row['source_row'], category, mid))
                    already_new += category == 'new'
                elif check not in IMPORTABLE:
                    decisions.append((method, row['source_row'], check, latest['id'] if latest else None))
                else:
                    report[(method, 'changed_during_preparation')] += 1
            con.commit()
        except Exception:
            con.rollback()
            raise
    if not apply:
        decisions.extend((method, row['source_row'], category, target['id'] if target else None)
                         for row, category, target, _ in plans)
    else:
        stage.executemany('INSERT INTO import_result VALUES(?,?,?,?)', decisions)
        stage.commit()
    report.update((m, state) for m, _, state, _ in decisions)


def run(database, stage_path, calibration, *, apply, limit, overlaps_only, methods, workers,
        novel_heavy=None, new_quota=None):
    report = Counter()
    manifest = json.loads((SOURCE / 'manifest.json').read_text(encoding='utf-8'))
    hashes = {Path(f['file']).stem: f['sha256'] for f in manifest['files'] if f['file'].endswith('.hdf5')}
    with ExitStack() as stack:
        stage = stack.enter_context(closing(sqlite3.connect(stage_path, timeout=60)))
        stage.row_factory = sqlite3.Row
        _, indexed_db, minimum, _, revision = json.loads(
            stage.execute("SELECT value FROM settings WHERE key='scan'").fetchone()[0])
        if Path(indexed_db) != database or revision != manifest['revision']:
            raise ValueError('Database or source revision differs from scanned inputs')
        if apply:
            con = stack.enter_context(closing(core.connect(database)))
            con.execute('PRAGMA cache_size=-196608')
            computed_spectra.initialize(con)
            con.commit()
            stage.execute("""CREATE TABLE IF NOT EXISTS import_result(
                method TEXT, source_row INTEGER, category TEXT, molecule_id TEXT,
                PRIMARY KEY(method,source_row))""")
            stage.commit()
        else:
            con = stack.enter_context(closing(sqlite3.connect(database.as_uri()+'?mode=ro',uri=True)))
            con.row_factory = sqlite3.Row
        pool = stack.enter_context(ProcessPoolExecutor(
            max_workers=workers, initializer=init_worker,
            initargs=(calibration['models'], manifest['revision'], hashes)))
        for method in methods:
            if apply:
                with h5py.File(file_for(method)) as f:
                    group = f[METHODS[method][1]]
                    size = len(group['atom_start' if method == 'shake' else 'n_atoms'])
                scanned = stage.execute(
                    "SELECT sum(json_extract(counts,'$.input')) FROM scan WHERE method=?", (method,)).fetchone()[0]
                if size != scanned:
                    raise ValueError(f'Finish scanning {method} before importing it ({scanned}/{size})')
            pending = deque()
            for i, rows in enumerate(groups(stage, method, limit=limit, overlaps_only=overlaps_only,
                                            novel_heavy=novel_heavy, new_quota=new_quota), 1):
                if apply:
                    done = {r[0] for r in stage.execute(
                        'SELECT source_row FROM import_result WHERE method=? AND source_row BETWEEN ? AND ?',
                        (method, rows[0]['source_row'], rows[-1]['source_row']))}
                    rows = [r for r in rows if r['source_row'] not in done]
                if not rows:
                    continue
                selected, decisions = [], []
                for row in rows:
                    existing = prior_import(con, row)
                    if existing:
                        decisions.append((method, row['source_row'], *existing))
                        continue
                    category, target = classify(con, row, minimum)
                    if category in IMPORTABLE:
                        selected.append((dict(row), category, dict(target) if target else None))
                    else:
                        decisions.append((method, row['source_row'], category, target['id'] if target else None))
                future = pool.submit(prepare_batch, method, selected) if selected else None
                pending.append((decisions, selected, future))
                if len(pending) >= 2 * workers:
                    commit_batch(stage, con, method, *pending.popleft(), report, minimum, apply, new_quota)
                    if new_quota is not None and report[(method, 'new')] >= new_quota:
                        break
                if i % 25 == 0:
                    print(method, i, json.dumps({s: n for (m, s), n in report.items() if m == method}),
                          flush=True)
            while pending:
                commit_batch(stage, con, method, *pending.popleft(), report, minimum, apply, new_quota)
        if apply:
            report = Counter({(m, state): n for m, state, n in stage.execute(
                'SELECT method,category,count(*) FROM import_result GROUP BY method,category')})
    return {m: {s: n for (method, s), n in report.items() if m == method} for m in DATA_ORDER}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', type=Path, default=core.DEFAULT_DB)
    parser.add_argument('--stage', type=Path, default=OUT / 'magnet_structure_index.sqlite')
    parser.add_argument('--calibration', type=Path, default=OUT / 'shielding_calibration.json')
    parser.add_argument('--method', action='append', choices=DATA_ORDER)
    parser.add_argument('--limit', type=int)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--overlaps-only', action='store_true')
    parser.add_argument('--novel-heavy-atoms', type=int)
    parser.add_argument('--new-quota', type=int)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    if args.overlaps_only and args.novel_heavy_atoms is not None:
        parser.error('Pick overlaps or one novel atom-count stratum')
    if args.new_quota is not None and (not args.apply or args.new_quota < 1 or
                                      args.novel_heavy_atoms is None or len(args.method or []) != 1):
        parser.error('Apply a positive quota to one method and one novel atom-count stratum')
    RDLogger.DisableLog('rdApp.*')
    calibration = json.loads(args.calibration.read_text(encoding='utf-8'))
    report = run(core.main_database(args.db).resolve(), args.stage, calibration, apply=args.apply,
                 limit=args.limit, overlaps_only=args.overlaps_only,
                 methods=tuple(m for m in DATA_ORDER if not args.method or m in args.method), workers=args.workers,
                 novel_heavy=args.novel_heavy_atoms, new_quota=args.new_quota)
    if args.apply:
        (args.stage.parent / 'magnet_import_summary.json').write_text(json.dumps(report, indent=2),
                                                                      encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
