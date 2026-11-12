"""Compare the local 2DNMRGym and NP-MRD files with the live annotation database."""
import collections
import concurrent.futures
import csv
import gzip
import io
import json
from pathlib import Path
import re
import sqlite3
import sys
import urllib.request
import zipfile

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from annotation_tool_v3.core import ACTIVE
from annotation_tool_v3.migrate import solvent
from annotation_tool_v3.peak_table import effective_peaks

DATA = ROOT / 'data'
DB = DATA / 'assignment_database_v1/databases/assignments.sqlite'
OUT = ROOT / 'outputs/experimental_2d_nmr_overlap_20261112'
GITHUB = 'https://raw.githubusercontent.com/siriusxiao62/2dNMR/main/data_csv/2dnmr/'
GYM_INDEX = GITHUB + 'nmr_smile_solventclass_filtered4_3dgnn.csv'
GYM_TEST = 'https://drive.google.com/uc?export=download&id=1ChbiPb19OZ2JveQdcDaw4KlpQP8GoCCm'
NP_CSV = 'https://np-mrd.org/system/downloads/current/smiles_{}.csv.gz'
NP_LOCAL = DATA / 'np_mrd_experimental/raw/npmrd_natural_products_NP0000001_NP0050000_json.zip'
NP_PEAKS = DATA / 'np_mrd_experimental/raw/NP-MRD_nmr_peak_lists_NP0000001_NP0050000.zip'
NP_RAW = NP_LOCAL.parent
BONDS = [Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE, Chem.BondType.AROMATIC]
DIRS = [Chem.BondDir.NONE, Chem.BondDir.ENDUPRIGHT, Chem.BondDir.ENDDOWNRIGHT]
CHIRAL = [Chem.ChiralType.CHI_UNSPECIFIED, Chem.ChiralType.CHI_TETRAHEDRAL_CW,
          Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_OTHER]
HYBRID = [Chem.HybridizationType.S, Chem.HybridizationType.SP, Chem.HybridizationType.SP2,
          Chem.HybridizationType.SP3, Chem.HybridizationType.SP3D,
          Chem.HybridizationType.SP3D2, Chem.HybridizationType.UNSPECIFIED]
NP_NAME = re.compile(r'(\d+)_(NP\d+)_(\d+)_(assignmenttable|user_assignmenttable|chsqc_assignmenttable)\.txt')


def get(url, cache):
    path = OUT / 'source_index' / cache
    if not path.exists():
        request = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(request, timeout=18) as response:
            payload = response.read()
        if payload.lstrip().startswith((b'<!DOCTYPE', b'<html')):
            raise ValueError(f'Expected data: {url}')
        path.write_bytes(payload)
    return path.read_bytes()


def csv_bytes(data):
    return list(csv.DictReader(io.StringIO(data.decode('utf-8-sig'), newline='')))


def connect():
    connection = sqlite3.connect(DB.as_uri() + '?mode=ro', uri=True)
    connection.row_factory = sqlite3.Row
    connection.execute('PRAGMA query_only=ON')
    return connection


def identity(mol):
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    key = Chem.MolToInchiKey(mol)
    return {'smiles': smiles, 'inchikey': key, 'connectivity': key[:14],
            'formula': rdMolDescriptors.CalcMolFormula(mol), 'heavy_atoms': mol.GetNumHeavyAtoms(),
            'mw': round(Descriptors.MolWt(mol), 2)}


def graph_equal(graph, mol):
    explicit = Chem.AddHs(mol)
    features = [[a.GetAtomicNum()-1, CHIRAL.index(a.GetChiralTag()), HYBRID.index(a.GetHybridization())]
                for a in explicit.GetAtoms()]
    actual = {(a, b): tuple(v) for (a, b), v in zip(zip(*graph['edge_index']), graph['edge_attr'])}
    expected = {}
    for bond in explicit.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        code = (BONDS.index(bond.GetBondType()), DIRS.index(bond.GetBondDir()))
        expected[a, b] = expected[b, a] = code
    return features == graph['x'] and expected == actual


def graph_mol(graph):
    mol = Chem.RWMol()
    for num, chiral, _ in graph['x']:
        atom = Chem.Atom(int(num) + 1)
        atom.SetNoImplicit(True)
        atom.SetChiralTag(CHIRAL[int(chiral)])
        mol.AddAtom(atom)
    seen = set()
    for (a, b), (bond, direction) in zip(zip(*graph['edge_index']), graph['edge_attr']):
        edge = frozenset((a, b))
        if edge not in seen:
            mol.AddBond(a, b, BONDS[bond])
            mol.GetBondBetweenAtoms(a, b).SetBondDir(DIRS[direction])
            seen.add(edge)
    Chem.SanitizeMol(mol)
    return Chem.RemoveHs(mol)


def gym_data():
    indexes = collections.defaultdict(list)
    for url, cache in [(GYM_INDEX, 'gym_author_metadata.csv'), (GYM_TEST, 'transpeaknet_test_smiles.csv')]:
        for r in csv_bytes(get(url, cache)):
            indexes[Path(r['File_name']).stem].append(r)
    data, peaks = [], []
    for line in (DATA/'eval_data.txt').read_text(encoding='utf-8-sig').splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        filename, graph = row['filename'], row['graph_data']
        sources = []
        for entry in indexes[filename]:
            mol = Chem.MolFromSmiles(entry['SMILES'])
            if mol is not None and graph_equal(graph, mol):
                sources.append((mol, entry))
        keys = {Chem.MolToInchiKey(mol) for mol, _ in sources}
        conditions = {r['Solvent'] for _, r in sources if r.get('Solvent') and
                      r.get('solvent_class') == str(graph['solvent_class']) and r['Solvent'].lower() != 'unknown'}
        known_solvent = next(iter(conditions)) if len(conditions) == 1 else ''
        basis = 'author_smiles_graph_verified' if len(keys) == 1 else 'reconstructed_graph_only'
        if len(keys) == 1:
            mol = sources[0][0]
        else:
            try:
                mol = graph_mol(graph)
            except ValueError:
                mol = None
        info = identity(mol) if mol is not None else {
            'heavy_atoms': sum(a[0] != 0 for a in graph['x']), 'mw': None, 'formula': '', 'smiles': '',
            'inchikey': '', 'connectivity': ''}
        carbons = sorted({a for a, b in zip(*graph['edge_index'])
                          if graph['x'][a][0] == 5 and graph['x'][b][0] == 0})
        if len(carbons) != len(row['c_peaks']) or len(carbons) != len(row['h_peaks']):
            raise ValueError(f'Carbon index mismatch for {filename}')
        pairs = set()
        for idx, (c,), hs in zip(carbons, row['c_peaks'], row['h_peaks']):
            c_ppm = round(c*200, 4)
            for h in sorted(set(hs)):
                h_ppm = round(h*10, 4)
                pairs.add((h_ppm, c_ppm))
                peaks.append({'source_id': filename, 'source_carbon_index_0based': idx,
                              'h_ppm': h_ppm, 'c_ppm': c_ppm, 'basis': 'expert_HSQC'})
        data.append({'source': '2DNMRGym', 'source_id': filename, **info, 'identity_basis': basis,
                     'source_smiles': sources[0][1]['SMILES'] if len(keys) == 1 else '',
                     'solvent_raw': known_solvent, 'solvent_class': graph['solvent_class'],
                     'hsqc_coordinate_pairs': len(pairs), '_pairs': sorted(pairs)})
    return data, peaks


def np_data():
    grouped, files = collections.defaultdict(list), {}
    for path in sorted((DATA/'assignment_tables').glob('*.txt')):
        match = NP_NAME.fullmatch(path.name)
        if not match:
            raise ValueError(f'Unexpected filename: {path.name}')
        submission, accession, batch, kind = match.groups()
        grouped[accession].append((submission, batch, kind, path))
        files[submission, accession, batch, kind] = path
    with zipfile.ZipFile(NP_LOCAL) as z:
        metadata = json.loads(z.read(z.namelist()[0]))['np_mrd']['natural_product']
    refs = {r['accession']: r for r in metadata if r['accession'] in grouped}
    aliases = collections.defaultdict(list)
    for ref in metadata:
        secondary = (ref.get('secondary_accessions') or {}).get('accession') or []
        for alias in [secondary] if isinstance(secondary, str) else secondary:
            if alias in grouped:
                aliases[alias].append(ref)
    for alias, candidates in aliases.items():
        if alias not in refs and len(candidates) == 1:
            refs[alias] = candidates[0]
    missing = set(grouped) - set(refs)
    partitions = {max(0, (int(a[2:])-1)//50_000) for a in missing}

    def partition(part):
        start = 1 if not part else part*50_000 + 1
        end = (part+1)*50_000
        name = f'NP{start:07}_NP{end:07}'
        return csv_bytes(gzip.decompress(get(NP_CSV.format(name), f'np_smiles_{name}.csv.gz')))

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for rows in pool.map(partition, sorted(partitions)):
            for r in rows:
                accession = next((v for v in r.values() if isinstance(v, str) and re.fullmatch(r'NP\d+', v)), '')
                if accession in missing:
                    fields = {k.strip().lower(): v for k, v in r.items()}
                    refs[accession] = {'accession': accession, **fields}
    for accession in set(grouped) - set(refs):
        path = NP_RAW / f'{accession}.sdf'
        if path.exists():
            block = path.read_text(encoding='utf-8-sig')
            mol = Chem.MolFromMolBlock(block)
            if mol is not None:
                refs[accession] = {'accession': accession, 'smiles': Chem.MolToSmiles(mol)}
    backfilled = collections.defaultdict(list)
    with zipfile.ZipFile(NP_PEAKS) as z:
        for n in z.namelist():
            accession = Path(n).name.split('_', 1)[0]
            if '_nmrtwod_' in n and accession in grouped:
                backfilled[accession].append(Path(n).name)
    data, audit = [], []
    for accession, members in sorted(grouped.items()):
        ref = refs.get(accession, {})
        mol = Chem.MolFromSmiles(ref['smiles']) if ref.get('smiles') else None
        info = identity(mol) if mol is not None else {}
        official_key = ref.get('inchikey') or ref.get('inchi_key') or ''
        simulated_matches = 0
        for sub, batch, kind, path in members:
            if kind != 'chsqc_assignmenttable':
                continue
            companion = files.get((sub, accession, batch, 'assignmenttable'))
            sim = {(r[0], r[1]): float(r[2]) for r in read_table(companion)
                   if len(r) > 2 and r[2] not in ('', 'NA') and r[0] in ('H', 'C')}
            coords = [r for r in read_table(path) if len(r) >= 6 and r[0] == 'H' and r[2] == 'C']
            reproduced = sum(sim.get(('H', r[1])) == float(r[4]) and
                             sim.get(('C', r[3])) == float(r[5]) for r in coords)
            all_match = bool(coords) and reproduced == len(coords)
            simulated_matches += all_match
            audit.append({'source_id': accession, 'submission_id': sub, 'file': path.name,
                          'reference_table': companion.name if companion else '', 'crosspeak_rows': len(coords),
                          'exact_coordinate_rows': reproduced, 'fully_reproduced': all_match})
        counts = collections.Counter(t[2] for t in members)
        data.append({'source': 'NP-MRD', 'source_id': accession, **info,
                     'current_accession': ref.get('accession', ''),
                     'identity_basis': 'np_mrd_structure' if info else 'unresolved',
                     'official_inchikey': official_key, 'name': ref.get('name') or ref.get('natural_products_name') or '',
                     'table_files': len(members), 'np_hsqc_files': counts['chsqc_assignmenttable'],
                     'np_reproduced_hsqc_files': simulated_matches,
                     'np_user_files': counts['user_assignmenttable'],
                     'backfilled_2d_files': len(backfilled[accession]),
                     'backfilled_2d_file_names': '|'.join(backfilled[accession]),
                     'source_url': 'https://np-mrd.org/natural_products/' + accession})
    return data, audit


def read_table(path):
    if not path:
        return []
    return list(csv.reader(io.StringIO(path.read_text(encoding='utf-8-sig'))))


def size_bin(n):
    if n is None:
        return 'unknown'
    return next((label for bound, label in [(15,'<=15'),(25,'16-25'),(35,'26-35'),(45,'36-45')]
                 if n <= bound), '>=46')


def db_index(c, wanted):
    matches = collections.defaultdict(list)
    totals = collections.Counter()
    for r in c.execute('SELECT id,inchikey,connectivity,formula,mw,queue FROM molecules'):
        formula = r['formula']
        n = sum(int(count or 1) for element, count in re.findall(r'([A-Z][a-z]?)(\d*)', formula) if element != 'H')
        totals[size_bin(n), r['queue']] += 1
        if r['connectivity'] in wanted:
            matches[r['connectivity']].append(dict(r))
    return matches, totals


def overlap(row, index):
    layer = row.get('connectivity', '')
    matches = index.get(layer, [])
    key = row.get('inchikey', '')
    exact = [m for m in matches if m['inchikey'] == key]
    trusted = row['identity_basis'] in ('author_smiles_graph_verified', 'np_mrd_structure')
    key_conflict = bool(row.get('official_inchikey') and row['official_inchikey'] != key)
    if not layer:
        status = 'unresolved'
    elif not trusted:
        status = 'graph_connectivity_hit' if matches else 'graph_no_connectivity_hit'
    elif key_conflict:
        status = 'source_key_conflict'
    elif exact:
        status = 'exact'
    else:
        status = 'connectivity_only' if matches else 'new_connectivity'
    row.update(match=status, size_bin=size_bin(row.get('heavy_atoms')),
               db_exact_ids='|'.join(m['id'] for m in exact) if status == 'exact' else '',
               db_connectivity_ids='|'.join(m['id'] for m in matches),
               db_exact_queues='|'.join(m['queue'] for m in exact) if status == 'exact' else '')


def norm_solvent(raw):
    if not raw:
        return None
    if raw.strip().upper() == 'C5D5N':
        return 'pyridine-d5'
    return solvent(raw)


def near_peak(peaks, ppm, tolerance):
    for p in peaks:
        lo, hi = p.get('range_min'), p.get('range_max')
        lo = p['shift_ppm'] if lo is None else lo
        hi = p['shift_ppm'] if hi is None else hi
        if lo-tolerance <= ppm <= hi+tolerance:
            return True
    return False


def gym_pair_candidates(c, gym):
    rows = []
    for r in gym:
        for mid in filter(None, r['db_exact_ids'].split('|')):
            spectra = {}
            for s in c.execute(f'SELECT * FROM spectra WHERE molecule_id=? AND {ACTIVE}', (mid,)):
                source = json.loads(s['source_json'])
                peaks = effective_peaks(json.loads(s['peaks_json']), source)
                drop = {a[0] for a in c.execute("SELECT peak_id FROM annotations WHERE spectrum_id=? AND status IN ('solvent','impurity','excluded')", (s['id'],))}
                spectra[s['id']] = (s, source, [p for p in peaks if p['peak_id'] not in drop])
            for p in c.execute('SELECT * FROM pairs WHERE molecule_id=?', (mid,)):
                if p['h_id'] not in spectra or p['c_id'] not in spectra:
                    continue
                h, h_src, h_peaks = spectra[p['h_id']]
                carbon, c_src, c_peaks = spectra[p['c_id']]
                sol = [norm_solvent(v) for v in (r['solvent_raw'], h['solvent'], carbon['solvent'])]
                known = set(filter(None, sol))
                condition = 'solvent_conflict' if len(known) > 1 else 'same_known_solvent' if all(sol) else 'solvent_unknown'
                hits = sum(near_peak(h_peaks, hp, .04) and near_peak(c_peaks, cp, .4) for hp, cp in r['_pairs'])
                num = len(r['_pairs'])
                rows.append({'source_id': r['source_id'], 'heavy_atoms': r['heavy_atoms'], 'molecule_id': mid,
                             'queue': r['db_exact_queues'], 'pair_id': p['id'], 'h_spectrum_id': h['id'],
                             'c_spectrum_id': carbon['id'], 'h_doi': h_src.get('doi'), 'c_doi': c_src.get('doi'),
                             'h_solvent': h['solvent'], 'c_solvent': carbon['solvent'], 'hsqc_solvent': r['solvent_raw'],
                             'condition': condition, 'pair_basis': p['basis'], 'hsqc_pairs': num, 'both_axes_hits': hits,
                             'coverage': round(hits/num, 4) if num else 0,
                             'screen_pass': num >= 5 and hits/num >= .90 and condition == 'same_known_solvent'})
    return rows


def write_csv(name, rows):
    cols = list(dict.fromkeys(k for r in rows for k in r if not k.startswith('_')))
    with (OUT / name).open('w', encoding='utf-8-sig', newline='') as out:
        writer = csv.DictWriter(out, cols, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main():
    RDLogger.DisableLog('rdApp.*')
    (OUT/'source_index').mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        gym_job, np_job = pool.submit(gym_data), pool.submit(np_data)
        gym, peaks = gym_job.result()
        print('2DNMRGym:', len(gym), collections.Counter(r['identity_basis'] for r in gym), flush=True)
        np, audit = np_job.result()
        print('NP-MRD:', len(np), collections.Counter(r['identity_basis'] for r in np), flush=True)
    with connect() as c:
        c.execute('BEGIN')
        index, totals = db_index(c, {r.get('connectivity') for r in gym+np if r.get('connectivity')})
        for r in gym+np:
            overlap(r, index)
        align = gym_pair_candidates(c, gym)
    write_csv('2dnmrgym_overlap.csv', gym)
    write_csv('2dnmrgym_expert_hsqc.csv', peaks)
    write_csv('2dnmrgym_existing_pairs.csv', align)
    write_csv('npmrd_overlap.csv', np)
    write_csv('npmrd_hsqc_provenance.csv', audit)
    gym_by_key = collections.defaultdict(list)
    for r in gym:
        if r['identity_basis'] == 'author_smiles_graph_verified':
            gym_by_key[r['inchikey']].append(r['source_id'])
    shared = [dict(inchikey=r['inchikey'], np_id=r['source_id'], heavy_atoms=r['heavy_atoms'],
                   gym_ids='|'.join(gym_by_key[r['inchikey']]), db_match=r['match'])
              for r in np if r['match'] not in ('source_key_conflict','unresolved') and r['inchikey'] in gym_by_key]
    write_csv('cross_source_exact_overlap.csv', shared)
    write_csv('database_size_distribution.csv', [dict(size_bin=b, queue=q, count=n) for (b,q),n in sorted(totals.items())])
    write_csv('large_molecule_candidates.csv', sorted((r for r in gym+np if (r.get('heavy_atoms') or 0) >= 26),
                                                    key=lambda r: (r['source'], r['match'], -r['heavy_atoms'])))
    print('Match:', collections.Counter(r['match'] for r in gym), collections.Counter(r['match'] for r in np), flush=True)
    print('Paired screens:', len(align), sum(r['screen_pass'] for r in align), flush=True)


if __name__ == '__main__':
    main()
