"""Count current hard assignments and the contribution of the three NMR imports."""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from contextlib import closing
from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import sqlite3

from rdkit import RDLogger

from annotation_tool_v3 import computed_spectra, core
from annotation_tool_v3.export_policy import current_spectra, hard_assignment


CASPRE = computed_spectra.CASPRE
DFT8K = computed_spectra.DFT8K
IBM = computed_spectra.IBM
ADDED = {CASPRE, DFT8K, IBM}
NUCLEI = ('13C', '1H')


def read_jsonl(path):
    with path.open(encoding='utf-8') as stream:
        return [json.loads(line) for line in stream if line.strip()]


def import_categories():
    folders = core.ROOT / 'outputs'
    rows = read_jsonl(folders / 'nmr_balance_selection_20261119/external_nmr_import_records.jsonl')
    values = {(r['source'], r['source_key']): r['decision'] for r in rows if 'spectra' in r}
    folder = folders / 'ibm_computational_nmr_audit_20261119'
    for row in read_jsonl(folder / 'ibm_novel_dual_candidates.jsonl'):
        values[IBM, str(row['source_id'])] = 'new'
    for row in read_jsonl(folder / 'ibm_exact_overlap_assignments.jsonl'):
        if row['main_assignment_coverage'] in ('unassigned', 'no_spectrum'):
            values[IBM, str(row['source_id'])] = row['main_assignment_coverage']
    return values


def source_key(s):
    source = s['provenance'].get('original_source') or s['source']
    p = s['provenance']
    key = (f"{p.get('source_sdf_record_1based')}:{p.get('source_record_id')}" if source == CASPRE else
           str(s.get('source_record_id') or p.get('source_record_id')))
    return source, key


def profile(spectra):
    return [(s, sum(hard_assignment(p) for p in s['peaks'])) for s in spectra if s['peaks']]


def accumulate(counts, entries):
    labelled = [(s, n) for s, n in entries if n]
    if not labelled:
        return
    nuclei = {s['nucleus'] for s, _ in entries}
    by_label = {s['nucleus'] for s, _ in labelled}
    counts['molecules'] += 1
    counts['full_molecules' if all(n == len(s['peaks']) for s, n in entries) else 'partial_molecules'] += 1
    counts['labelled_in_both_nuclei'] += len(by_label) == 2
    counts['full_dual_molecules'] += len(nuclei) == 2 and all(n == len(s['peaks']) for s, n in entries)
    counts['observed_dual' if len(nuclei) == 2 else 'observed_c_only' if '13C' in nuclei else 'observed_h_only'] += 1
    for s, n in labelled:
        nuc = s['nucleus']
        counts[f'{nuc}_spectra'] += 1
        counts[f'{nuc}_full_spectra' if n == len(s['peaks']) else f'{nuc}_partial_spectra'] += 1
        counts[f'{nuc}_peaks'] += n


def summarize(c, kinds):
    panels = defaultdict(Counter)
    details = defaultdict(lambda: defaultdict(Counter))
    gains = Counter()
    assignments = defaultdict(Counter)
    heavy_atoms = defaultdict(Counter)
    mids = c.execute('''SELECT id FROM molecules WHERE c_assigned+h_assigned>0
                       UNION SELECT molecule_id FROM computed_spectra
                       UNION SELECT molecule_id FROM spectra WHERE source=?''', (CASPRE,))
    for i, (mid,) in enumerate(mids, 1):
        r = core.record(c, mid, require_quality=True, include_computed=True)
        if r['curation']['molecule_state'] != 'active':
            continue
        spectra = current_spectra(r)
        exp = profile(s for s in spectra if s.get('spectrum_kind') != 'computed')
        comp = profile(s for s in spectra if s.get('spectrum_kind') == 'computed'
                       and s['provenance'].get('structure_match') != 'tautomer_variant')
        baseline = [(s, n) for s, n in exp if source_key(s)[0] not in ADDED]
        for name, entries in [('experimental_prior_sources', baseline), ('experimental', exp), ('computed', comp)]:
            accumulate(panels[name], entries)
        pre, ex, sim = (any(n for _, n in entries) for entries in (baseline, exp, comp))
        if ex or sim:
            gains['distinct_labelled_molecules'] += 1
            gains['labelled_in_both_domains'] += ex and sim
            gains['net_first_labelled_molecules'] += not pre
            gains['first_labelled_via_experiment'] += not pre and ex
            gains['first_labelled_via_simulation_only'] += not ex and sim
        for s, n in exp+comp:
            if not n:
                continue
            for peak in s['peaks']:
                if hard_assignment(peak):
                    origin = peak['assignment'].get('origin') or 'unknown'
                    assignments[s['nucleus']]['AI' if origin.startswith('ai:') else origin] += 1
        sources = defaultdict(list)
        groups = defaultdict(list)
        per_record = defaultdict(list)
        for s, n in exp+comp:
            src, key = source_key(s)
            if src not in ADDED:
                continue
            sources[src].append((s, n))
            groups[src, kinds.get((src, key), 'unknown')].append((s, n))
            per_record[src, key].append((s, n))
        for src, entries in sources.items():
            accumulate(panels[src], entries)
            if any(n for _, n in entries):
                size = sum(a['symbol'] != 'H' for a in r['structure']['atoms'])
                band = '<=15' if size <= 15 else '16-25' if size <= 25 else '26-35' if size <= 35 else '36+'
                heavy_atoms[src][band] += 1
        for (src, group), entries in groups.items():
            accumulate(details[src][group], entries)
        for (src, _), entries in per_record.items():
            if {s['nucleus'] for s, n in entries if n == len(s['peaks'])} == set(NUCLEI):
                panels[src]['fully_assigned_hc_from_same_record'] += 1
        if i % 5000 == 0:
            print(f'Checked {i:,} candidate molecules', flush=True)
    return {'counts': panels, 'added_by_original_overlap': details, 'molecule_gains': gains,
            'peaks_by_assignment_origin_and_nucleus': assignments, 'heavy_atoms_among_added_labelled_molecules': heavy_atoms}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    snapshot = datetime.now(timezone(timedelta(hours=8))).isoformat(timespec='seconds')
    with closing(sqlite3.connect(core.DEFAULT_DB.resolve().as_uri()+'?mode=ro', uri=True)) as c:
        c.row_factory = sqlite3.Row
        c.execute('BEGIN')
        report = summarize(c, import_categories())
    report['counting_basis'] = 'Hard assigned/equivalent-group/overlap peaks in current quality-eligible spectra, active molecules'
    report['snapshot_as_of'] = snapshot
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
    print(args.output, flush=True)


if __name__ == '__main__':
    main()
