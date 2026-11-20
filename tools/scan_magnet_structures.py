"""Build a resumable structure index from the MagNET shielding files."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import closing
import hashlib
import itertools
import json
from pathlib import Path
import pickle
import re
import sqlite3
import sys

import h5py
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from annotation_tool_v3.paths import DEFAULT_DB

SOURCE = ROOT / 'data/nmr_atom_assignment_sources/magnet'
OUT = ROOT / 'outputs/magnet_atomnmr_import'
METHODS = {
    'pepper-reference': ('sigma-pepper', 'wp04_wb97xd_pcSseg2'),
    'pepper-pbe0': ('sigma-pepper', 'pbe0_pcSseg1'),
    'concentrate-pcm': ('sigma-concentrate', '/'),
    'shake': ('sigma-shake', '/'),
}
ELEMENTS = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 16: 'S', 17: 'Cl'}
ELEMENT_CODES = {symbol: 1 << (7 * i) for i, symbol in enumerate(ELEMENTS.values())}
NUMBER_CODES = np.array([ELEMENT_CODES.get(ELEMENTS.get(z), 0) for z in range(18)], dtype=np.int64)
FORMULA_PART = re.compile(r'([A-Z][a-z]?)([0-9]*)')
FORMAT_VERSION = 1


def file_for(method):
    name, _ = METHODS[method]
    return SOURCE / name / f'{name}.hdf5'


def verify_sources():
    manifest = json.loads((SOURCE / 'manifest.json').read_text(encoding='utf-8'))
    for entry in manifest['files']:
        path = SOURCE / entry['file']
        if path.suffix == '.hdf5':
            with path.open('rb') as f:
                if path.stat().st_size != entry['bytes'] or hashlib.file_digest(f, 'sha256').hexdigest() != entry['sha256']:
                    raise ValueError(f'Invalid MagNET snapshot: {path}')
    return manifest


def formula_code(formula):
    parts = FORMULA_PART.findall(formula)
    if ''.join(symbol + n for symbol, n in parts) != formula:
        return None
    if any(symbol not in ELEMENT_CODES or int(n or 1) >= 128 for symbol, n in parts):
        return None
    return sum(ELEMENT_CODES[symbol] * int(n or 1) for symbol, n in parts)


def xyz_molecule(z, coordinates):
    xyz = '\n'.join(f'{ELEMENTS[int(n)]} {x:.4f} {y:.4f} {z:.4f}'
                    for n, (x, y, z) in zip(z, coordinates))
    return Chem.MolFromXYZBlock(f'{len(z)}\n\n{xyz}\n')


def finish_geometry(mol, z):
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    heavy = Chem.RemoveHs(mol)
    if [a.GetAtomicNum() for a in mol.GetAtoms()] != z.tolist():
        return None
    return mol, heavy, Chem.MolToSmiles(heavy)


def recover(z, coordinates, *, alternate=False):
    try:
        mol = xyz_molecule(z, coordinates)
        if mol is None:
            return None
        rdDetermineBonds.DetermineBonds(mol, charge=0, useVdw=alternate,
                                       covFactor=1.25, maxIterations=2000)
        return finish_geometry(mol, z)
    except (RuntimeError, ValueError):
        return None


def neighbor_signature(mol):
    """Exact structures must share their element-labeled unweighted heavy graph."""
    labels = [a.GetAtomicNum() for a in mol.GetAtoms()]
    neighbors = [[] for _ in labels]
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if labels[a] > 1 and labels[b] > 1:
            neighbors[a].append(labels[b])
            neighbors[b].append(labels[a])
    return tuple(sorted((z, tuple(sorted(n))) for z, n in zip(labels, neighbors) if z > 1))


def match_graph(z, xyz, formula):
    try:
        mol = xyz_molecule(z, xyz)
        if mol is None:
            return 'geometry_recovery_failed', None
        rdDetermineBonds.DetermineConnectivity(mol)
        if (formula, neighbor_signature(mol)) not in NEIGHBORS:
            return 'outside_novel_size_and_absent_neighbor_signature', None
        rdDetermineBonds.DetermineBondOrders(mol, charge=0, maxIterations=2000)
        recovered = finish_geometry(mol, z)
        return ('candidate' if recovered else 'geometry_recovery_failed'), recovered
    except (RuntimeError, ValueError):
        return 'geometry_recovery_failed', None


def neighbor_index(db, minimum, out, *, lower=0, extend=None):
    RDLogger.DisableLog('rdApp.*')
    neighbors = set()
    if extend:
        with extend.open('rb') as f:
            neighbors = pickle.load(f)
    with closing(sqlite3.connect(db.as_uri() + '?mode=ro', uri=True)) as c:
        formulas = {}
        for (formula,) in c.execute('SELECT DISTINCT formula FROM molecules'):
            code = formula_code(formula)
            if code is not None and lower <= sum((code >> (7*i)) & 127 for i in range(1,len(ELEMENTS))) < minimum:
                formulas[formula] = code
        if lower:
            rows = itertools.chain.from_iterable(
                c.execute('SELECT formula,smiles FROM molecules WHERE formula=?', (f,))
                for f in sorted(formulas)
            )
        else:
            rows = c.execute('SELECT formula,smiles FROM molecules')
        for formula, smiles in rows:
            if formula in formulas and (mol := Chem.MolFromSmiles(smiles)) is not None:
                neighbors.add((formulas[formula], neighbor_signature(mol)))
    with out.open('wb') as f:
        pickle.dump(neighbors, f)
    print(f'Indexed {len(neighbors)} original heavy-neighbor signatures', flush=True)


def declared_geometry(z, xyz, smiles):
    """Use declared bond orders and verify they agree with the XYZ atom order."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        return None
    inferred = recover(z, xyz)
    if not inferred:
        return None
    source, heavy, source_smiles = inferred
    if Chem.MolToSmiles(mol) == source_smiles:
        return source, heavy, source_smiles
    # The generated GDB SMILES may omit a stereocenter. Keep this as a
    # mismatch for exact import, even when the connectivity is consistent.
    return None


def read_index(method):
    with h5py.File(file_for(method)) as f:
        group = f[METHODS[method][1]]
        if method == 'shake':
            start, stop = group['atom_start'][:], group['atom_end'][:]
            ids = group['molecule_id'][:]
            failed = np.isin(ids, group['remove_molecule_ids'][:])
        else:
            offsets = np.r_[0, np.cumsum(group['n_atoms'][:], dtype=np.int64)]
            start, stop = offsets[:-1], offsets[1:]
            failed = np.zeros(len(start), dtype=bool)
            failed[group['failed_indices'][:]] = True
        return start, stop, failed


def initialize_worker(db, minimum, neighbors):
    global FORMULAS, LAYERS, MINIMUM, NEIGHBORS
    RDLogger.DisableLog('rdApp.*')
    MINIMUM = minimum
    NEIGHBORS = None
    if neighbors:
        with neighbors.open('rb') as f:
            NEIGHBORS = pickle.load(f)
    with closing(sqlite3.connect(db.as_uri() + '?mode=ro', uri=True)) as c:
        FORMULAS = {code for (f,) in c.execute('SELECT DISTINCT formula FROM molecules')
                    if (code := formula_code(f)) is not None}
        LAYERS = {key[:14] for (key,) in c.execute('SELECT inchikey FROM molecules WHERE inchikey IS NOT NULL')}


def scan_block(job):
    method, first, start, stop, failed = job
    counts = Counter(input=len(start), excluded_by_source=int(failed.sum()))
    rows = []
    with h5py.File(file_for(method)) as f:
        group = f[METHODS[method][1]]
        first_atom, end_atom = int(start[0]), int(stop[-1])
        base = slice(first_atom, end_atom)
        starts, stops = start - first_atom, stop - first_atom
        z = group['atomic_numbers'][base]
        inventory = np.add.reduceat(NUMBER_CODES[z], starts)
        heavy_count = np.add.reduceat(z != 1, starts)
        dual = ((inventory & 127) > 0) & (((inventory >> 7) & 127) > 0)
        counts['no_dual_nuclei'] = int((~dual & ~failed).sum())
        considered = ~failed & dual
        possible = np.fromiter((value in FORMULAS for value in inventory), dtype=bool, count=len(inventory))
        wanted = considered & (possible | (heavy_count >= MINIMUM))
        counts['outside_novel_size_and_absent_formula'] = int((considered & ~wanted).sum())
        if not wanted.any():
            return method, first, counts, rows
        coord = group['stationary_coords' if method == 'shake' else 'coordinates']
        ids = group['molecule_id'][first:first + len(start)] if method == 'shake' else None
        linked = group['sigma_shake_id'][first:first + len(start)] if method == 'concentrate-pcm' else None
        if method.startswith('pepper'):
            offsets = group['smiles_offsets'][first:first + len(start) + 1]
            buf = group['smiles_buffer'][offsets[0]:offsets[-1]].tobytes()
            offsets -= offsets[0]
        else:
            xyz = coord[base] * float(coord.attrs['scale'])
        for i in np.flatnonzero(wanted):
            n = int(heavy_count[i])
            atoms = slice(starts[i], stops[i])
            external = first + int(i) if ids is None else int(ids[i])
            if method.startswith('pepper'):
                declared = buf[offsets[i]:offsets[i + 1]].decode()
                mol = Chem.MolFromSmiles(declared)
                if mol is None or len(Chem.GetMolFrags(mol)) != 1:
                    counts['invalid_declared_graph'] += 1
                    continue
                smiles = Chem.MolToSmiles(mol)
            else:
                if NEIGHBORS is not None and n < MINIMUM:
                    state, recovered = match_graph(z[atoms], xyz[atoms], inventory[i])
                else:
                    recovered = recover(z[atoms], xyz[atoms])
                    state = 'geometry_recovery_failed'
                if not recovered:
                    counts[state] += 1
                    continue
                _, mol, smiles = recovered
            key = Chem.MolToInchiKey(mol)
            if len(key) != 27:
                counts['invalid_inchikey'] += 1
                continue
            in_original = key[:14] in LAYERS
            if not in_original and n < MINIMUM:
                counts['novel_below_minimum_heavy_atoms'] += 1
                continue
            formula = rdMolDescriptors.CalcMolFormula(mol)
            if formula_code(formula) != inventory[i]:
                counts['identity_failure'] += 1
                continue
            counts['overlap_connectivity_candidate' if in_original else 'novel_size_candidate'] += 1
            rows.append((method, first + int(i), external, None if linked is None else int(linked[i]),
                         key, smiles, n, formula, int(in_original)))
    return method, first, counts, rows


def scan(db, stage, minimum, screening_floor, workers, block_size, *, neighbors=None, baseline=None):
    manifest = verify_sources()
    with sqlite3.connect(stage) as c:
        c.execute('PRAGMA journal_mode=WAL')
        c.executescript("""
        CREATE TABLE IF NOT EXISTS settings(key TEXT PRIMARY KEY, value TEXT);
        CREATE TABLE IF NOT EXISTS scan(method TEXT, first INTEGER, counts TEXT, PRIMARY KEY(method, first));
        CREATE TABLE IF NOT EXISTS candidates(
            method TEXT, source_row INTEGER, source_id INTEGER, shake_id INTEGER,
            inchikey TEXT, smiles TEXT, heavy_atoms INTEGER, formula TEXT, original_layer INTEGER,
            PRIMARY KEY(method, source_row));
        """)
        settings = json.dumps([FORMAT_VERSION, str(db), minimum, block_size, manifest['revision']])
        prior = c.execute("SELECT value FROM settings WHERE key='scan'").fetchone()
        if prior and prior[0] != settings:
            raise ValueError('The staging file was created with different inputs')
        c.execute("INSERT OR IGNORE INTO settings VALUES('scan',?)", (settings,))
        c.commit()
        done = set(c.execute('SELECT method,first FROM scan'))
        jobs = []
        for method in METHODS:
            start, stop, failed = read_index(method)
            for first in range(0, len(start), block_size):
                if (method, first) not in done:
                    sl = slice(first, first + block_size)
                    jobs.append((method, first, start[sl], stop[sl], failed[sl]))
        # Finish the inexpensive SMILES-backed overlap work, then coordinates.
        jobs.sort(key=lambda job: (job[0] == 'shake', job[1]))
        with ProcessPoolExecutor(max_workers=workers, initializer=initialize_worker,
                                 initargs=(baseline or db, screening_floor, neighbors)) as workers:
            futures = [workers.submit(scan_block, job) for job in jobs]
            for j, future in enumerate(as_completed(futures), 1):
                method, first, counts, rows = future.result()
                c.executemany('INSERT INTO candidates VALUES(?,?,?,?,?,?,?,?,?)', rows)
                c.execute('INSERT INTO scan VALUES(?,?,?)', (method, first, json.dumps(counts)))
                c.commit()
                if j % 10 == 0 or j == len(jobs):
                    print(f'Indexed {j}/{len(jobs)} blocks this run', flush=True)
        c.execute('CREATE INDEX IF NOT EXISTS candidates_key ON candidates(original_layer,inchikey,smiles)')
        counts_by_method = {}
        for method, records in itertools.groupby(c.execute('SELECT method,counts FROM scan ORDER BY method'), lambda r: r[0]):
            counts = Counter()
            for _, stats in records:
                counts.update(json.loads(stats))
            counts_by_method[method] = counts
        print(json.dumps(counts_by_method, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', type=Path, default=DEFAULT_DB)
    parser.add_argument('--stage', type=Path, default=OUT / 'magnet_structure_index.sqlite')
    parser.add_argument('--minimum-novel-heavy-atoms', type=int, default=16)
    parser.add_argument('--screening-floor', type=int)
    parser.add_argument('--workers', type=int, default=6)
    parser.add_argument('--block-size', type=int, default=12_000)
    parser.add_argument('--neighbor-index', type=Path)
    parser.add_argument('--build-neighbor-index', action='store_true')
    parser.add_argument('--extend-neighbor-index', type=Path)
    parser.add_argument('--neighbor-index-lower', type=int, default=0)
    parser.add_argument('--baseline', type=Path)
    args = parser.parse_args()
    args.stage.parent.mkdir(parents=True, exist_ok=True)
    baseline = args.baseline.resolve() if args.baseline else None
    if args.build_neighbor_index:
        neighbor_index(baseline or args.db.resolve(), args.minimum_novel_heavy_atoms, args.neighbor_index,
                       lower=args.neighbor_index_lower, extend=args.extend_neighbor_index)
    else:
        floor = args.screening_floor or args.minimum_novel_heavy_atoms
        if floor < args.minimum_novel_heavy_atoms:
            parser.error('Screening floor must be at least the novel admission threshold')
        scan(args.db.resolve(), args.stage, args.minimum_novel_heavy_atoms, floor, args.workers, args.block_size,
             baseline=baseline, neighbors=args.neighbor_index)


if __name__ == '__main__':
    main()
