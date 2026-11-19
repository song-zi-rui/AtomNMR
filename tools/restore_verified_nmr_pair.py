"""Restore an experimental pair after a full simulated pair and primary-source review."""

from argparse import ArgumentParser
from contextlib import closing
import json
from pathlib import Path
import re

from annotation_tool_v3 import computed_spectra, core, quality
from annotation_tool_v3.export_policy import coverage


def prepare(c, evidence):
    mid = evidence['molecule_id']
    record = core.record(c, mid, include_archived=True, include_computed=True)
    if record['curation']['molecule_state'] != 'active':
        raise ValueError('The molecule is manually held or excluded')
    by_id = {s['id']: s for s in record['spectra']}
    spectra = {nuc: by_id[sid] for nuc, sid in evidence['spectrum_ids'].items()}
    carbon, proton = spectra['13C'], spectra['1H']
    experimental_pair = c.execute('SELECT id FROM pairs WHERE molecule_id=? AND c_id=? AND h_id=?',
                                  (mid, carbon['id'], proton['id'])).fetchone()
    if not experimental_pair:
        raise ValueError('The reviewed experimental pair no longer exists')
    old = carbon['quality']
    if all(s['quality'].get('source_review') == evidence and s['quality']['decision'] == 'passed'
           for s in spectra.values()):
        return None
    simulations = [s for s in record['spectra'] if s.get('spectrum_kind') == 'computed'
                   and s['provenance'].get('structure_match') != 'tautomer_variant']
    supporting = next((p for p in computed_spectra.pairs(simulations)
                       if by_id[p['c_id']]['source'] == evidence['computed_source']
                       and by_id[p['c_id']]['source_record_id'] == evidence['computed_source_record_id']
                       and all(coverage([by_id[p[key]]]) == 'full' for key in ('c_id', 'h_id'))), None)
    if not supporting:
        raise ValueError('No matching simulated record with full assignments for both nuclei')
    if (sorted(evidence['spectrum_ids']) != ['13C', '1H']
        or not re.fullmatch('[0-9a-f]{64}', evidence['document_sha256'])
        or evidence['pages'][1] != evidence['pages'][0] + 1 or not evidence['compound']):
        raise ValueError('The verified supplement reference is incomplete')
    for nucleus, s in spectra.items():
        if (s['quality'] != old or s['quality']['decision'] != 'quarantined'
            or s['quality']['reasons'] != ['structure_text_on_different_pages']
            or s['source'] != 'NMRexp' or s['nucleus'] != nucleus or s['manual_state'] != 'active'
            or s['provenance'].get('doi') != evidence['doi']
            or s['solvent'] != evidence['solvent'] or s['frequency_mhz'] != evidence['frequency_mhz'][nucleus]):
            raise ValueError('The spectrum differs from the reviewed NMRexp record or has another hold')
        observed = [p['shift_ppm'] for p in s['peaks']]
        verified = evidence['verified_shifts_ppm'][nucleus]
        match = quality.peak_match(observed, verified, 1e-6)
        if match['f1'] != 1.0:
            raise ValueError(f'The stored {nucleus} peaks differ from the supplement')
    candidate = c.execute('SELECT * FROM candidates WHERE molecule_id=? AND state=?',
                          (mid,'quality_held')).fetchall()
    if len(candidate) != 1:
        raise ValueError('Expected a single held source candidate')
    candidate = candidate[0]
    paired = json.loads(candidate['source_pair_json'])
    if (paired.get('c_source_id') != old['source_ids']['13C']
        or paired.get('h_source_id') != old['source_ids']['1H']):
        raise ValueError('The held candidate points to different NMRexp source rows')
    result = {**old, 'decision': 'passed', 'reasons': [], 'source_document_verified': True,
              'source_anchor_check': 'published_supporting_information_same_compound',
              'source_review': evidence, 'supporting_computed_pair': supporting,
              'prior_quarantine_reasons': old['reasons']}
    return candidate, old, result


def run(c, evidence, apply=False):
    c.execute('BEGIN IMMEDIATE' if apply else 'BEGIN')
    try:
        prepared = prepare(c, evidence)
        if prepared is None:
            return {'state': 'already_restored', 'molecule_id': evidence['molecule_id']}
        candidate, previous, result = prepared
        if apply:
            mid, ids = evidence['molecule_id'], evidence['spectrum_ids']
            quality.write_result(c, mid, candidate, result, ids.values())
            c.execute("UPDATE candidates SET state='imported',reason=NULL WHERE smiles=?", (candidate['smiles'],))
            c.execute('''INSERT INTO quality_history(molecule_id,smiles,old_pair_json,new_pair_json,details_json,created_at)
                         VALUES(?,?,?,?,?,?)''',
                      (mid, candidate['smiles'], candidate['source_pair_json'], candidate['source_pair_json'],
                       core.dumps({**result, 'action': 'verified_experimental_source_restore', 'previous_quality': previous}), core.now()))
            c.execute('INSERT INTO curation_events(molecule_id,action,details_json,created_at) VALUES(?,?,?,?)',
                      (mid, 'verified_experimental_source_restore', core.dumps({
                          'compound': evidence['compound'], 'doi': evidence['doi'], 'pages': evidence['pages'],
                          'spectrum_ids': ids,
                      }), core.now()))
            c.commit()
        return {'state': 'restored' if apply else 'eligible', 'molecule_id': evidence['molecule_id'],
                'spectrum_ids': evidence['spectrum_ids'], 'supporting_computed_pair': result['supporting_computed_pair']}
    finally:
        if c.in_transaction:
            c.rollback()


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('evidence', type=Path)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    evidence = json.loads(args.evidence.read_text(encoding='utf-8'))
    with closing(core.connect()) as connection:
        print(json.dumps(run(connection, evidence, args.apply), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
