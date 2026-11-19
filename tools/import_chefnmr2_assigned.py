"""Import verified published atom assignments for ChefNMR2 benchmark molecules."""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from contextlib import closing
import hashlib
import json
from pathlib import Path
import re
import sqlite3
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

from annotation_tool_v3 import core
from annotation_tool_v3.export_policy import hard_assignment


DATA = core.ROOT / 'chefnmr2_atompeak/data'
TAXANE_DOI = '10.1038/s41586-025-09090-z'
TAXANE_HASH = 'cead16f1647f8b01503d956828ad58dead1b3628bddd2db2635fcf3a2505f879'
TAXANE_URL = 'https://pmc-oa-opendata.s3.amazonaws.com/PMC12240809.1/41586_2025_9090_MOESM1_ESM.docx'
RAPAMYCIN_URL = 'https://www.jstage.jst.go.jp/article/antibiotics1968/44/6/44_6_688/_pdf'
RAPAMYCIN_HASH = '9e1f13b9fbddb6ac1510328446d6afd7c989097d4a239501c6f6e68537689131'
RAPAMYCIN_ERRATUM = {
    'citation': 'J Antibiot (Tokyo) 1991 Oct;44(10):C-3',
    'index_url': 'https://pubmed.ncbi.nlm.nih.gov/1712766/',
    'verification': 'Indexed erratum located; its full text could not be verified.',
}
# Compound table, first shift column, and O-acetyl attachment sites, as drawn.
TAXANES = {
    'taxane-001': ('S3', 1, ()),
    'taxane-002': ('S4', 1, (5, 9, 10, 13)),
    'taxane-003': ('S5', 1, (5, 9, 10, 13)),
    'taxane-004': ('S6', 1, (5, 9, 10, 13)),
    'taxane-005': ('S10', 1, (5, 13)),
    'taxane-006': ('S4', 3, (5, 9, 10, 13)),
    'taxane-007': ('S8', 1, ()),
}
# Table 1, major form in CDCl3. The remaining columns in the publication are
# the minor form and older assignments, which are not joined to benchmark peaks.
RAPAMYCIN_MAJOR = """1 169.2
2 51.3 5.29
3 27.0 2.34 1.76
4 20.6 1.78 1.47
5 25.3 1.75 1.48
6 44.2 3.59 3.44
8 166.8
9 192.5
10 98.5
11 33.7 1.98
12 27.3 1.60 1.60
13 31.3 1.62 1.33
14 67.2 3.86
15 38.8 1.85 1.52
16 84.4 3.67
17 135.5
18 129.6 5.97
19 126.4 6.39
20 133.6 6.32
21 130.1 6.15
22 140.2 5.54
23 35.2 2.32
24 40.2 1.50 1.20
25 41.4 2.74
26 215.6
27 84.9 3.71
28 77.3 4.17
29 136.1
30 126.8 5.42
31 46.6 3.33
32 208.2
33 40.7 2.74 2.60
34 75.7 5.17
35 33.1 1.98
36 38.4 1.22 1.12
37 33.2 1.39
38 34.2 2.10 0.68
39 84.4 2.93
40 73.9 3.37
41 31.3 1.99 1.33
42 31.7 1.70 1.00
43 16.2 0.95
44 10.2 1.65
45 21.5 1.05
46 13.8 1.00
47 13.0 1.74
48 16.0 1.11
49 15.9 0.92
50 55.8 3.13
51 59.5 3.34
52 56.5 3.41"""
SOURCE_NOTES = {
    'known-002': 'Original supplement Table 1 numbers signals for comparison, with no atom index. The characterization text lists shifts and peak types.',
    'known-003': 'The numbered shifts in doi:10.1039/D1RA07915H differ from the benchmark spectrum. No matching source table was verified.',
    'known-004': 'The original supplement contains plotted spectra and shift picks without an atom-assignment table.',
    'known-009': 'The original supplement contains spectra of intermediates and derivatives. No matching atom table verified.',
}
FILE_NOTES = {
    'novel_natural_products': 'The ChefNMR2 paper releases the experimental peak lists and structure without an explicit peak-to-atom table.',
    'specteach': 'The benchmark uses automatically picked peaks from experimental spectra. The JSONL has no explicit peak-to-atom labels.',
    'spectranp+': 'The ChefNMR2 Mnova simulation release has no source atom-to-peak mapping.',
    'uspto': 'The released simulated test peak lists have no atom-to-peak mapping.',
}


def table_text(node):
    return ''.join(n.text or '' for n in node.iter() if n.tag.endswith('}t')).strip()


def read_taxane_tables(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != TAXANE_HASH:
        raise ValueError('Unexpected version of the publisher supplement')
    with ZipFile(path) as archive:
        body = ET.fromstring(archive.read('word/document.xml')).find('.//{*}body')
    tables, caption = {}, None
    for node in body:
        if node.tag.endswith('}p'):
            if match := re.match(r'Table (S\d+)\.', table_text(node)):
                caption = match[1]
        elif node.tag.endswith('}tbl') and caption:
            tables[caption] = [[table_text(c) for c in tr.findall('{*}tc')] for tr in node.findall('{*}tr')]
            caption = None
    return tables


class NumberedStructure:
    """Connectivity transcribed from the numbered published structure drawings."""
    def __init__(self):
        self.graph, self.positions = Chem.RWMol(), {}

    def atom(self, label, element='C'):
        if label not in self.positions:
            self.positions[label] = self.graph.AddAtom(Chem.Atom(element))
        return self.positions[label]

    def bond(self, left, right, double=False):
        self.graph.AddBond(self.atom(left), self.atom(right), Chem.BondType.DOUBLE if double else Chem.BondType.SINGLE)

    def oxygen(self, site, double=False):
        label = f'O{site}'
        self.atom(label, 'O')
        self.bond(site, label, double)
        return label

    def maps(self, mol):
        query = self.graph.GetMol()
        Chem.SanitizeMol(query)
        if Chem.MolToSmiles(mol, isomericSmiles=False) != Chem.MolToSmiles(query, isomericSmiles=False):
            raise ValueError('Numbered drawing and benchmark graph differ')
        matches = mol.GetSubstructMatches(query, uniquify=False, maxMatches=1000)
        if not 0 < len(matches) < 1000:
            raise ValueError('Failed to enumerate numbering correspondences')
        return [{label: match[idx] for label, idx in self.positions.items()} for match in matches]


def taxane_structure(source_id, acetates):
    s = NumberedStructure()
    epoxy = not acetates
    for n in range(1, 14):
        s.bond(n, n + 1, n == 11)
    bridge = (1, 11) if source_id == 'taxane-004' else (11, 15)
    for left, right in [(14, 1), (1, 15), (15, 16), (15, 17), (12, 18), (8, 19), (3, 8), bridge]:
        s.bond(left, right)
    s.bond(4, 20, not epoxy)
    if epoxy:
        s.bond(s.oxygen(4), 20)
        s.oxygen(5)
    for site in acetates:
        o = s.oxygen(site)
        acyl = f'{site}-acetyl-carbonyl'
        s.bond(o, acyl)
        s.oxygen(acyl, double=True)
        s.bond(acyl, f'{site}-acetyl-methyl')
    extra = {'taxane-003': (1, False), 'taxane-004': (15, False), 'taxane-007': (13, True)}
    if source_id in extra:
        s.oxygen(*extra[source_id])
    return s


def published_positions(label, nucleus, ppm, acetates):
    if match := re.match(r'^(\d+)(?:[ab](?:-\([EZ]\))?)?$', label):
        return [int(match[1])]
    if match := re.fullmatch(r'C(\d+)-OH', label):
        return [f'O{match[1]}']
    if match := re.fullmatch(r'(?:C(\d+))?-OCOMe', label):
        group = 'methyl' if nucleus == '1H' or ppm < 50 else 'carbonyl'
        sites = [int(match[1])] if match[1] else acetates
        return [f'{n}-acetyl-{group}' for n in sites]
    raise ValueError(f'Unknown published atom label: {label}')


def taxane_observations(rows, left, acetates):
    out = {'13C': defaultdict(list), '1H': defaultdict(list)}
    for number, row in enumerate(rows[2:], 3):
        label = row[0]
        for nucleus, col in [('13C', left), ('1H', left + 1)]:
            raw = row[col]
            if not (match := re.match(r'^(~?)\s*(\d+(?:\.\d+)?)', raw)):
                continue
            ppm = float(match[2])
            desc = re.search(r'\(([^;)]+)', raw)
            mult = desc[1].strip().removeprefix('br, ').removeprefix('br ') if desc else None
            if mult == 'overlapping':
                mult = 'm'
            j = raw.partition(';')[2].partition(')')[0]
            j = [float(n) for n in re.findall(r'\d+(?:\.\d+)?', j)]
            out[nucleus][ppm].append({'row': number, 'label': label, 'reported': raw,
                                     'positions': published_positions(label, nucleus, ppm, acetates),
                                     'approximate': bool(match[1]), 'multiplicity': mult, 'j_hz': j})
    return out


def resolve_peak(observations, maps, conflicted):
    if not observations:
        return {'status': 'unassigned', 'atom_indices': []}
    targets = {tuple(sorted({m[p] for row in observations for p in row['positions']})) for m in maps}
    atoms = sorted({i for t in targets for i in t})
    uncertain = len(targets) != 1 or any(row['approximate'] or len(row['positions']) != 1 or
                 any(p in conflicted for p in row['positions']) for row in observations)
    return {'status': 'ambiguous' if uncertain else 'assigned' if len(atoms) == 1 else 'overlap',
            'atom_indices': atoms} if atoms else {'status': 'unassigned', 'atom_indices': []}


def taxane_peaks(row, mol, tables):
    table, column, acetates = TAXANES[row['id']]
    observed = taxane_observations(tables[table], column, acetates)
    maps = taxane_structure(row['id'], acetates).maps(mol)
    carbon_sites = defaultdict(set)
    for ppm, entries in observed['13C'].items():
        for item in entries:
            if len(item['positions']) == 1:
                carbon_sites[item['positions'][0]].add(ppm)
    conflicts = {site for site, shifts in carbon_sites.items() if len(shifts) > 1}
    duplicate_h = Counter(h[0] for h in row['h'])
    spectra = {}
    for nucleus, prefix, peaks in [('13C', 'C', [(p, None, None) for p in row['c']]),
                                    ('1H', 'H', row['h'])]:
        specs = []
        for i, (ppm, nh, mult) in enumerate(peaks, 1):
            entries = observed[nucleus].get(ppm, [])
            if nucleus == '1H' and duplicate_h[ppm] > 1:
                entries = [x for x in entries if x['multiplicity'] == mult]
            if not entries:
                raise ValueError(f'Benchmark peak not matched in table {table}: {nucleus} {ppm}')
            result = resolve_peak(entries, maps, conflicts if nucleus == '13C' else set())
            for idx in result['atom_indices']:
                atom = mol.GetAtomWithIdx(idx)
                if (nucleus == '13C' and atom.GetAtomicNum() != 6) or (nucleus == '1H' and not atom.GetTotalNumHs()):
                    raise ValueError(f'Published atom incompatible with {nucleus}: {row["id"]} {ppm}')
            if nucleus == '1H' and hard_assignment({'assignment': result}):
                accessible_h = sum(mol.GetAtomWithIdx(i).GetTotalNumHs() for i in result['atom_indices'])
                if accessible_h < nh:
                    raise ValueError(f'Published proton assignment cannot explain integral: {row["id"]} {ppm}')
            peak = {'peak_id': f'{prefix}{i}', 'shift_ppm': ppm, 'assignment': result,
                    'published_assignment_rows': entries}
            if nucleus == '1H':
                peak.update(reported_nH=nh, multiplicity=mult)
                if len(entries) == 1 and result['status'] == 'assigned' and entries[0]['j_hz']:
                    peak['j_values'] = entries[0]['j_hz']
            specs.append(peak)
        if not any(hard_assignment(p) for p in specs):
            raise ValueError(f'No determined assignments in verified table {table} for {nucleus}')
        spectra[nucleus] = specs
    provenance = {
        'doi': TAXANE_DOI, 'table': f'Supplementary Table {table}',
        'pdf': 'McClune et al. — Supplementary Information', 'document_url': TAXANE_URL,
        'document_sha256': TAXANE_HASH, 'benchmark_id': row['id'],
        'benchmark_file': 'taxanes.jsonl', 'compound_name': row['name'],
        'source_structure_numbering': 'Numbered skeletal drawing above the assignment table',
        'source_atom_map_candidates': maps, 'original_text': {'13C': row['c_text'], '1H': row['h_text']},
    }
    return spectra, provenance


def rapamycin_structure():
    # The macrolide numbering follows the Rapamycin drawing on page 688.
    s = NumberedStructure()
    s.atom(7, 'N')
    s.bond(1, 2)
    s.bond(2, 7)
    for n in range(2, 7):
        s.bond(n, n + 1)
    for n in range(7, 42):
        s.bond(n, n + 1, n in (17, 19, 21, 29))
    for left,right in [(42,37),(11,43),(17,44),(23,45),(25,46),(29,47),(31,48),(35,49)]:
        s.bond(left, right)
    for n in (1,8,9,26,32):
        s.oxygen(n, double=True)
    for n in (10,28,40):
        s.oxygen(n)
    for label,left,right in [('ester_oxygen',1,34), ('ring_oxygen',10,14),
                             ('methoxy16',16,50),('methoxy27',27,51),('methoxy39',39,52)]:
        s.atom(label, 'O')
        s.bond(left, label)
        s.bond(label, right)
    return s


def rapamycin_peaks(row, mol, document):
    if hashlib.sha256(document.read_bytes()).hexdigest() != RAPAMYCIN_HASH:
        raise ValueError('Unexpected version of the rapamycin assignment paper')
    maps = rapamycin_structure().maps(mol)
    reference = {'13C': defaultdict(list), '1H': defaultdict(list)}
    for line in RAPAMYCIN_MAJOR.splitlines():
        site, carbon, *hydrogens = line.split()
        reference['13C'][float(carbon)].append(int(site))
        for ppm in hydrogens:
            reference['1H'][float(ppm)].append(int(site))
    spectra = {}
    for nucleus, prefix, peaks in [('13C', 'C', [(p, None, None) for p in row['c']]),
                                    ('1H', 'H', row['h'])]:
        rows = []
        for i, (ppm, nh, mult) in enumerate(peaks, 1):
            # Admit exact reported positions only. For a shared value, keep
            # candidate positions rather than assuming the benchmark signal is
            # the superposition of the full published list.
            sites = reference[nucleus].get(ppm, [])
            entries = [{'label': n, 'positions': [n], 'reported': str(ppm), 'approximate': False}
                       for n in sorted(set(sites))]
            label = resolve_peak(entries, maps, set())
            if len(set(sites)) > 1:
                label['status'] = 'ambiguous'
            if nucleus == '1H' and label['status'] == 'assigned':
                h_count = mol.GetAtomWithIdx(label['atom_indices'][0]).GetTotalNumHs()
                if not h_count:
                    raise ValueError('Rapamycin proton site has no attached hydrogen')
                if nh > h_count:
                    label['status'] = 'ambiguous'
            peak = {'peak_id': f'{prefix}{i}', 'shift_ppm': ppm, 'assignment': label,
                    'published_assignment_rows': entries}
            if nucleus == '1H':
                peak.update(reported_nH=nh, multiplicity=mult)
            rows.append(peak)
        if any(hard_assignment(p) for p in rows):
            spectra[nucleus] = rows
    if not spectra:
        raise ValueError('Rapamycin benchmark has no exact shifts matching the numbered table')
    provenance = {
        'doi': '10.7164/antibiotics.44.688', 'table': 'Table 1, major form', 'pdf_pages': [1, 2],
        'pdf': 'McAlpine et al. — Revised NMR Assignments for Rapamycin',
        'document_url': RAPAMYCIN_URL, 'document_sha256': RAPAMYCIN_HASH,
        'benchmark_id': row['id'], 'benchmark_file': 'known_natural_products.jsonl',
        'compound_name': row['name'], 'source_structure_numbering': 'Numbered structure on page 688',
        'source_atom_map_candidates': maps, 'matching_rule': 'exact_major_form_table_shift_only',
        'original_text': {'13C': row['c_text'], '1H': row['h_text']},
        'erratum': RAPAMYCIN_ERRATUM,
    }
    return spectra, provenance


def benchmark_rows():
    for path in sorted(DATA.glob('*.jsonl')):
        if path.stem not in ('nmrexp', 'nmrshiftdb2'):
            with path.open(encoding='utf-8') as f:
                for line in f:
                    if line.strip():
                        yield path.stem, json.loads(line)


def scan(c, prepared):
    decisions, plans = [], []
    for file, row in benchmark_rows():
        mol = Chem.MolFromSmiles(row['smiles'])
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True) if mol else None
        key = Chem.MolToInchiKey(mol) if mol else ''
        prefix = key[:14]
        hits = list(c.execute('SELECT id,inchikey FROM molecules WHERE inchikey>=? AND inchikey<?',
                              (prefix+'-', prefix+'.'))) if key else []
        exact = [mid for mid, k in hits if k == key]
        if not exact and smiles:
            exact = [r[0] for r in c.execute('SELECT id FROM molecules WHERE smiles=?', (smiles,))]
        info = {'file': file, 'id': row['id'], 'name': row.get('name'), 'formula': row['formula'],
                'inchikey': key, 'smiles': smiles, 'existing_ids': [mid for mid, _ in hits]}
        if not key or len(key) != 27 or rdMolDescriptors.CalcMolFormula(mol) != row['formula']:
            status = 'invalid_identity'
        elif exact:
            status = 'existing_exact'
        elif hits:
            status = 'existing_connectivity'
        elif row['id'] not in prepared:
            status = 'novel_without_verified_atom_assignments'
            info['source_check'] = SOURCE_NOTES.get(row['id'], FILE_NOTES.get(file))
        else:
            status = 'prepared'
            plans.append((info, prepared[row['id']]))
        decisions.append({**info, 'decision': status})
    return decisions, plans


def save(c, plans):
    entries = []
    for info, (spectra, provenance, conditions) in plans:
        mid = core.add_molecule(c, info['smiles'], 'literature_pdf')
        sids = {}
        counts = {}
        for nucleus, items in spectra.items():
            labels = {p['peak_id']: p['assignment'] for p in items}
            source_peaks = [{k: v for k, v in p.items() if k != 'assignment'} for p in items]
            source_key = f'chefnmr2:{info["id"]}:{nucleus}'
            sid = core.add_spectrum(c, mid, nucleus, 'literature_pdf', source_key, conditions[nucleus],
                                    source_peaks, {**provenance, 'benchmark_release': 'ChefNMR2'}, labels)
            sids[nucleus] = sid
            counts[nucleus] = dict(Counter(p['assignment']['status'] for p in items))
        if set(sids) == {'13C', '1H'}:
            core.add_pair(c, mid, sids['13C'], sids['1H'], 'same_publication_assignment_table_sample_unknown')
        core.refresh_molecule(c, mid)
        entries.append({**info, 'molecule_id': mid, 'spectrum_ids': sids, 'peak_assignments': counts,
                        'source': {k: provenance[k] for k in ('doi', 'table', 'document_url', 'erratum')
                                   if k in provenance}})
    return entries


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--supplement', required=True, type=Path)
    parser.add_argument('--rapamycin-paper', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--apply', action='store_true')
    options = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    taxane_tables = read_taxane_tables(options.supplement)
    prepared = {}
    for file, row in benchmark_rows():
        if file == 'taxanes':
            mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(row['smiles']), isomericSmiles=True))
            spectra, provenance = taxane_peaks(row, mol, taxane_tables)
            prepared[row['id']] = (spectra, provenance, {n: {'solvent': 'CDCl3', 'solvent_raw': 'CDCl3'} for n in spectra})
        if row['id'] == 'known-007' and options.rapamycin_paper:
            mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(row['smiles']), isomericSmiles=True))
            spectra, provenance = rapamycin_peaks(row, mol, options.rapamycin_paper)
            prepared[row['id']] = (spectra, provenance, {n: {'solvent': 'CDCl3', 'solvent_raw': 'CDCl3'} for n in spectra})
    database = core.DEFAULT_DB.resolve()
    connection = (core.connect(database) if options.apply else sqlite3.connect(database.as_uri()+'?mode=ro', uri=True))
    with closing(connection) as c:
        c.execute('BEGIN IMMEDIATE' if options.apply else 'BEGIN')
        decisions, plans = scan(c, prepared)
        admitted = save(c, plans) if options.apply else []
        if options.apply:
            c.commit()
    counts = defaultdict(Counter)
    for row in decisions:
        counts[row['file']][row['decision']] += 1
    report = {'by_file': counts, 'screened': len(decisions), 'eligible': len(plans),
              'applied': options.apply, 'admitted': admitted}
    options.output.mkdir(parents=True, exist_ok=True)
    (options.output / 'summary.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    with (options.output / 'screened.jsonl').open('w', encoding='utf-8') as f:
        for row in decisions:
            f.write(json.dumps(row, ensure_ascii=False) + '\n')
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
