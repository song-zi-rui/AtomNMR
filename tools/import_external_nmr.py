"""Import validated CASCADE DFT8K and CASPRE CDCl3/JEOL spectra to AtomNMR."""

import argparse
from collections import Counter, defaultdict
from contextlib import closing
import hashlib
import json
from pathlib import Path
import sqlite3
import sys

from rdkit import Chem, RDLogger

from screen_balanced_nmr_sources import CASPRE, DFT, OUT, ROOT, read_caspre, read_dft, write_jsonl
from audit_ibm_computational_nmr import unkeyed_main_rows

sys.path.insert(0, str(ROOT))
from annotation_tool_v3 import computed_spectra as sources, core as db
from annotation_tool_v3.export_policy import coverage, current_spectra
from external_spectrum_atoms import target_mappings

URLS = {
    sources.DFT8K: 'https://github.com/patonlab/CASCADE/tree/master/data/DFT8K',
    sources.CASPRE: 'https://caspre.ca/carbonpreds/training_data',
}


def checksums():
    manifest = json.loads((DFT / 'manifest.json').read_text(encoding='utf-8'))
    expected = {DFT / f['name']: f['sha256'] for f in manifest['files']}
    expected[CASPRE] = 'cbdfd82253285316fd5e9a8fc032a890807d4cff728dcf08097106b6d010d13a'
    for path, digest in expected.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f'Checksum mismatch: {path}')
    return {p.name: digest for p, digest in expected.items()}


def existing_spectrum(c, row):
    table = 'computed_spectra' if row['source'] == sources.DFT8K else 'spectra'
    field = 'source_record_id' if table == 'computed_spectra' else 'source_key'
    return list(c.execute(f'SELECT molecule_id,nucleus FROM {table} WHERE source=? AND {field}=?',
                          (row['source'], source_key(row))))


def source_key(row):
    if row['source'] == sources.DFT8K:
        return row['source_id']
    return f"{row['source_sdf_record_1based']}:{row['source_id']}"


def candidates(c, row):
    layer = row['connectivity']
    return c.execute('SELECT id,smiles,inchikey,revision,structure_json FROM molecules WHERE inchikey>=? AND inchikey<? ORDER BY id',
                     (layer+'-', layer+'.')).fetchall()


def classify(c, row, matches, unkeyed, cached):
    exact = [m for m in matches if m['inchikey'] == row['inchikey']]
    same_graph = [m for m in exact if Chem.MolToSmiles(Chem.MolFromSmiles(m['smiles'])) == row['canonical_smiles']]
    if len(same_graph) != 1:
        if exact:
            return ('ambiguous_exact' if same_graph else 'exact_key_different_graph'), None
        if matches:
            return 'connectivity_only', None
        smiles_no_stereo = Chem.MolToSmiles(Chem.RemoveHs(row['source_mol']), isomericSmiles=False)
        return ('unkeyed_possible_match' if smiles_no_stereo in unkeyed else 'new'), None
    target = same_graph[0]
    mid = target['id']
    if mid not in cached:
        record = db.record(c, mid, require_quality=True, include_computed=True)
        usable = [s for s in current_spectra(record) if s['peaks']]
        cached[mid] = ('inactive_molecule' if record['curation']['molecule_state'] != 'active' else
                       'no_spectrum' if not usable else
                       'unassigned' if coverage(usable) == 'none' else 'has_valid_assignment')
    return cached[mid], target


def make_spectrum(shifts, nucleus, mappings, attached):
    by_shift = defaultdict(list)
    for shift in shifts:
        by_shift[shift['ppm']].append(shift)
    peaks = []
    prefix = 'C' if nucleus == '13C' else 'H'
    for i, (ppm, values) in enumerate(sorted(by_shift.items(), reverse=True), 1):
        source_ids = [s['atom_index'] for s in values]
        mapped = {tuple(sorted({m[attached[j] if nucleus == '1H' else j] for j in source_ids}))
                  for m in mappings}
        targets = sorted({j for group in mapped for j in group})
        status = ('unassigned' if not targets else 'ambiguous' if len(mapped) > 1
                  else 'assigned' if len(targets) == 1 else 'overlap')
        peaks.append({'peak_id': f'{prefix}{i}', 'shift_ppm': ppm,
                      'source_atom_ids': source_ids, 'source_atom_shifts': values,
                      'assignment': {'status': status, 'atom_indices': targets}})
    return peaks


def prepare(c, batches, unkeyed, file_hashes):
    plans, statuses, records = [], {}, []
    cached = {}
    for rows, source_info in batches:
        source = rows[0]['source']
        counts = Counter()
        used_keys = set()
        for row in rows:
            key = source_key(row)
            want = {'13C', '1H'} if source == sources.DFT8K else {'13C'}
            entry = {'source': source, 'source_id': row['source_id'], 'source_key': key,
                     'inchikey': row['inchikey'], 'formula': row['formula']}
            source_in_smiles = Chem.MolFromSmiles(row['canonical_smiles'])
            if source_in_smiles is None or Chem.MolToInchiKey(source_in_smiles) != row['inchikey']:
                counts['identity_roundtrip_failure'] += 1
                records.append({**entry, 'decision': 'identity_roundtrip_failure'})
                continue
            already = existing_spectrum(c, row)
            if already:
                if {x['nucleus'] for x in already} != want or len({x['molecule_id'] for x in already}) != 1:
                    raise ValueError(f'Inconsistent prior import for {source} {key}')
                counts['already_imported'] += 1
                continue
            if key in used_keys:
                raise ValueError(f'Repeated source key: {source} {key}')
            used_keys.add(key)
            matches = candidates(c, row)
            status, target = classify(c, row, matches, unkeyed, cached)
            if target:
                entry['molecule_id'] = target['id']
            if status not in ('new', 'unassigned', 'no_spectrum'):
                counts[status] += 1
                records.append({**entry, 'decision': status})
                continue
            smiles = target['smiles'] if target else row['canonical_smiles']
            mol = Chem.MolFromSmiles(smiles)
            struct = json.loads(target['structure_json']) if target else db.structure(mol)
            mappings, attached = target_mappings(row['source_mol'], mol, struct)
            peaks = {n: make_spectrum(row['c_shifts' if n == '13C' else 'h_shifts'], n, mappings, attached) for n in sorted(want)}
            mid = target['id'] if target else db.uid('M', smiles)
            provenance = {
                'source_record_id': row['source_id'], 'source_sdf_record_1based': row['source_sdf_record_1based'],
                'source_url': URLS[source], 'source_smiles': row['canonical_smiles'], 'source_files_sha256':
                    {k: v for k, v in file_hashes.items() if k.startswith('DFT') == (source == sources.DFT8K)},
                'label_origin': row['label_origin'], 'source_atom_index_basis': row['atom_index_basis'],
                'structure_match': 'exact_graph',
                'atom_mapping_status': 'validated_source_sdf' if mappings else 'source_indices_only',
            }
            if source == sources.CASPRE:
                provenance.update(compound_name=row['compound_name'], chemical_shift_reference='TMS')
            plans.append((row, mid, smiles, peaks, provenance, status, signatures(matches)))
            record = {'spectra': len(peaks), 'peaks': sum(map(len, peaks.values())),
                      'ambiguous_peaks': sum(p['assignment']['status']=='ambiguous' for ps in peaks.values() for p in ps),
                      'molecule_id': mid, 'decision': status, 'atom_mapping_status': provenance['atom_mapping_status']}
            records.append({**entry, **record})
            counts[status] += 1
        statuses[source] = {**source_info, 'validated': len(rows), 'decisions': dict(counts)}
    return plans, statuses, records


def signatures(matches):
    return [(r['id'], r['revision'], r['smiles'], r['inchikey']) for r in matches]


def apply(c, plans):
    c.execute('BEGIN IMMEDIATE')
    try:
        sources.initialize(c)
        for row, _, _, _, _, _, prior in plans:
            if signatures(candidates(c, row)) != prior or existing_spectrum(c, row):
                raise ValueError('The live data changed during preparation. Rerun the import.')
        refresh = set()
        for row, mid, smiles, spectra, provenance, status, _ in plans:
            source = row['source']
            if status == 'new':
                if db.add_molecule(c, smiles, source) != mid:
                    raise ValueError(f'Unexpected identifier for {source} {source_key(row)}')
                refresh.add(mid)
            for nucleus, peaks in spectra.items():
                if source == sources.DFT8K:
                    for p in peaks:
                        p['assignment'].update(origin=source, note='')
                    sid = db.uid('S', source, source_key(row), nucleus)
                    c.execute('INSERT INTO computed_spectra VALUES(?,?,?,?,?,?,?)',
                              (sid, mid, nucleus, source, source_key(row), db.dumps(peaks), db.dumps(provenance)))
                else:
                    labels = {p['peak_id']: p.pop('assignment') for p in peaks}
                    db.add_spectrum(c, mid, nucleus, source, source_key(row),
                                    {'solvent': row['solvent'], 'solvent_raw': row['solvent']}, peaks, provenance, labels)
                    refresh.add(mid)
            if status == 'no_spectrum':
                refresh.add(mid)
        for mid in sorted(refresh):
            db.refresh_molecule(c, mid)
        c.commit()
    except Exception:
        c.rollback()
        raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--db', type=Path, default=db.DEFAULT_DB)
    p.add_argument('--apply', action='store_true')
    args = p.parse_args()
    RDLogger.DisableLog('rdApp.*')
    file_hashes = checksums()
    database = db.main_database(args.db).resolve()
    batches = [read_dft(include_source_structure=True), read_caspre(include_source_structure=True)]
    with closing(sqlite3.connect(database.as_uri() + '?mode=ro', uri=True)) as c:
        c.row_factory = sqlite3.Row
        c.execute('BEGIN')
        unkeyed, _ = unkeyed_main_rows(c)
        plans, statuses, records = prepare(c, batches, unkeyed, file_hashes)
    accepted = [r for r in records if 'spectra' in r]
    report = {'sources': statuses, 'input_sha256': file_hashes,
              'source_records_to_import': len(accepted), 'spectra_to_import': sum(r['spectra'] for r in accepted),
              'unique_molecules': len({r['molecule_id'] for r in accepted}),
              'new_molecules': len({r['molecule_id'] for r in accepted if r['decision'] == 'new'}),
              'manifest': str(OUT / 'external_nmr_import_records.jsonl')}
    if args.apply and plans:
        with closing(db.connect(database)) as c:
            apply(c, plans)
        write_jsonl(OUT / 'external_nmr_import_records.jsonl', records)
        report['applied'] = True
        (OUT / 'external_nmr_import_summary.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
