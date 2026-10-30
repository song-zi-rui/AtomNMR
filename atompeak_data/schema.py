"""One observable spectrum view with explicitly masked assignment supervision."""
from __future__ import annotations

import hashlib
import math
import re
from collections import defaultdict

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

SCHEMA_VERSION = 'atompeak-partial-v1'
BRANCHES = ('carbon', 'proton')
NUCLEI = {'carbon': '13C', 'proton': '1H'}
HARD_STATUSES = {'assigned', 'equivalent_group', 'overlap'}
HYBRIDS = {'UNSPECIFIED': 0, 'S': 1, 'SP': 2, 'SP2': 3, 'SP3': 4, 'SP3D': 5, 'SP3D2': 6, 'OTHER': 7}
CHIRALS = {'CHI_UNSPECIFIED': 0, 'CHI_TETRAHEDRAL_CW': 1, 'CHI_TETRAHEDRAL_CCW': 2, 'CHI_OTHER': 3}
BONDS = {'SINGLE': 1, 'DOUBLE': 2, 'TRIPLE': 3, 'AROMATIC': 4, 'DATIVE': 5}


def finite(value, default=None):
    if isinstance(value, dict):
        value = value.get('normalized_value')
    try:
        value = float(value)
        return value if math.isfinite(value) else default
    except (TypeError, ValueError):
        return default


def j_values(value):
    if value in (None, '', '_'):
        return []
    if isinstance(value, (int, float)):
        value = [value]
    if isinstance(value, str):
        value = re.findall(r'\d+(?:\.\d+)?', value)
    return sorted([x for v in value if (x := finite(v)) is not None and x >= 0], reverse=True)


def split_for_connectivity(key, seed=42, valid_ratio=.1, test_ratio=.1):
    if not key or not all(math.isfinite(v) for v in (valid_ratio, test_ratio)) or min(valid_ratio, test_ratio) < 0 or valid_ratio + test_ratio >= 1:
        raise ValueError('Invalid structure key or split ratios')
    draw = int(hashlib.sha256(f'{seed}:{key}'.encode()).hexdigest()[:16], 16) / 2**64
    return 'valid' if draw < valid_ratio else 'test' if draw < valid_ratio + test_ratio else 'train'


def graph_and_candidates(smiles, stored_atoms=None, require_connectivity=True):
    parser = Chem.SmilesParserParams()
    parser.removeHs = False
    mol = Chem.MolFromSmiles(smiles, parser)
    if mol is None:
        raise ValueError('Invalid stored SMILES')
    if stored_atoms is not None:
        if len(stored_atoms) != mol.GetNumAtoms() or any(
                a['index'] != i or a['symbol'] != mol.GetAtomWithIdx(i).GetSymbol()
                for i, a in enumerate(stored_atoms)):
            raise ValueError('Stored atom numbering does not match SMILES traversal')
    nodes, edges, groups = [], [], {b: [] for b in BRANCHES}
    for a in mol.GetAtoms():
        nodes.append({'atom_index': a.GetIdx(), 'element': a.GetSymbol(), 'atomic_number': a.GetAtomicNum(),
                      'formal_charge': a.GetFormalCharge(), 'degree': a.GetDegree(),
                      'total_h_count': a.GetTotalNumHs(includeNeighbors=True), 'hybridization': str(a.GetHybridization()),
                      'aromatic': a.GetIsAromatic(), 'in_ring': a.IsInRing(), 'chirality': str(a.GetChiralTag())})
        if a.GetSymbol() == 'H':
            continue
        for branch in BRANCHES:
            include = a.GetSymbol() == 'C' if branch == 'carbon' else a.GetTotalNumHs(includeNeighbors=True) > 0
            if include:
                groups[branch].append({'group_id': f'{branch}:{a.GetIdx()}', 'center_atom_indices': [a.GetIdx()],
                    'observable_h_count': a.GetTotalNumHs(includeNeighbors=True), 'transport_mass': 1.,
                    'grouping_policy': 'individual_carbon_or_attached_heavy_atom'})
    for bond in mol.GetBonds():
        edges.append({'source': bond.GetBeginAtomIdx(), 'target': bond.GetEndAtomIdx(), 'bond_type': str(bond.GetBondType())})
    inchikey = Chem.MolToInchiKey(mol)
    if not inchikey and require_connectivity:
        raise ValueError('Cannot determine connectivity InChIKey; not assigned a potentially leaking split')
    return mol, {'nodes': nodes, 'edges': edges}, groups, inchikey.split('-')[0]


def graph_arrays(graph):
    heavy = [n for n in graph['nodes'] if n['element'] != 'H']
    mapping = {n['atom_index']: i for i, n in enumerate(heavy)}
    features = [[n['atomic_number'], max(-5, min(5, n['formal_charge'])) + 5,
                 min(8, n['degree']), min(8, n['total_h_count']), HYBRIDS.get(n['hybridization'], 7),
                 int(n['aromatic']), int(n['in_ring']), CHIRALS.get(n['chirality'], 3)] for n in heavy]
    edges = [e for e in graph['edges'] if e['source'] in mapping and e['target'] in mapping]
    return {'node_features': np.asarray(features, dtype=np.int64).reshape(-1, 8),
            'node_valid': np.ones(len(heavy), dtype=np.int64),
            'edge_index': np.asarray([[mapping[e['source']], mapping[e['target']]] for e in edges], dtype=np.int64).reshape(-1, 2),
            'edge_types': np.asarray([BONDS.get(e['bond_type'], 6) for e in edges], dtype=np.int64)}, mapping


def conditions_for(spec):
    if not spec:
        return {}
    value = spec.get('conditions') or {}
    return {'solvent': spec.get('solvent') or value.get('solvent'),
            'temperature_k': finite(spec.get('temperature_k', value.get('temperature_k'))),
            'field_strength_mhz': finite(spec.get('frequency_mhz')) if value.get('frequency_scope') != 'source_reported_unspecified' else None,
            'ph': finite(value.get('ph')), 'analyte_concentration': finite(value.get('analyte_concentration'))}


def observed_peak(peak, nucleus):
    shift = finite(peak.get('shift_ppm'))
    if shift is None:
        raise ValueError('An observed peak must have a finite shift')
    value = {'peak_id': str(peak['peak_id']), 'shift_ppm': shift}
    if nucleus == '13C':
        value['delta (ppm)'] = shift
    else:
        value.update(centroid=shift, delta=shift,
                     rangeMax=finite(peak.get('range_max'), shift), rangeMin=finite(peak.get('range_min'), shift),
                     category=str(peak.get('multiplicity') or 'unk').strip().replace(' ', '_'),
                     nH=finite(peak.get('reported_nH')), j_values=j_values(peak.get('j_values')))
    width = finite(peak.get('line_width_hz'))
    if width is not None and width >= 0:
        value['line_width_hz'] = width
    # Raw assignment IDs, SDTQ inferred from targets, and label confidence are
    # deliberately absent from observable model inputs.
    return value


def select_spectra(record):
    spectra = {s['id']: s for s in record['spectra'] if s.get('manual_state', 'active') == 'active'
               and s.get('quality', {}).get('decision') not in ('quarantined', 'not_checked') and s.get('peaks')}
    pairs = [p for p in record.get('pairs', []) if p['c_id'] in spectra and p['h_id'] in spectra]
    pairs.sort(key=lambda p: (-sum(bool(spectra[s].get('version_info')) for s in (p['c_id'], p['h_id'])),
                             p['solvent_match'] != 'same_known', 'cross_database' in p['basis'], p['id']))
    pair = pairs[0] if pairs else None
    chosen = {}
    for branch in BRANCHES:
        options = sorted([s for s in spectra.values() if s['nucleus'] == NUCLEI[branch]],
                         key=lambda s: (not bool(s.get('version_info')), -(s.get('version_info') or {}).get('number', 0), s['id']))
        chosen[branch] = spectra[pair['c_id' if branch == 'carbon' else 'h_id']] if pair else next(iter(options), None)
    return chosen, pair


def branch_supervision(spec, groups, mol):
    raw = sorted(spec['peaks'], key=lambda p: -float(p['shift_ppm'])) if spec else []
    lookup = {g['center_atom_indices'][0]: i for i, g in enumerate(groups)}
    pids = {p['peak_id']: i for i, p in enumerate(raw)}
    shape = (len(raw), len(groups))
    target, known, allowed = (np.zeros(shape, dtype=np.float32) for _ in range(3))
    row_transport = np.zeros(len(raw), dtype=bool)
    weak_rows = np.zeros(len(raw), dtype=bool)
    label_status, atom_labels = [], []
    grouped = {pid for g in (spec or {}).get('assignment_groups', []) for pid in g['peak_ids']}
    for i, peak in enumerate(raw):
        label = peak.get('assignment') or {}
        status = label.get('status', 'unassigned')
        atoms = sorted(set(label.get('atom_indices') or []))
        valid = (label.get('source_identity_eligible', True) and bool(atoms)
                 and all(type(a) is int and a in lookup for a in atoms) and peak['peak_id'] not in grouped)
        label_status.append(status)
        atom_labels.append(atoms)
        if not valid:
            continue
        columns = [lookup[a] for a in atoms]
        if status in HARD_STATUSES:
            target[i, columns] = 1
            known[i, :] = 1
            row_transport[i] = len(columns) == 1
        elif status in ('source_group', 'ambiguous'):
            allowed[i, columns] = 1
            weak_rows[i] = True
    set_peaks, set_groups = [], []
    for group in (spec or {}).get('assignment_groups', []):
        if not group.get('peak_ids') or not group.get('atom_indices'):
            continue
        if not all(pid in pids for pid in group['peak_ids']) or not all(a in lookup for a in group['atom_indices']):
            continue
        pm, gm = np.zeros(len(raw)), np.zeros(len(groups))
        pm[[pids[p] for p in group['peak_ids']]] = 1
        gm[[lookup[a] for a in group['atom_indices']]] = 1
        set_peaks.append(pm.tolist()); set_groups.append(gm.tolist())
    shift_target = np.zeros(len(groups), dtype=np.float32)
    shift_valid = np.zeros(len(groups), dtype=bool)
    for j, group in enumerate(groups):
        rows = np.flatnonzero(target[:, j])
        if len(rows) != 1:
            continue
        i = int(rows[0]); atom = mol.GetAtomWithIdx(group['center_atom_indices'][0])
        nh = finite(raw[i].get('reported_nH'))
        # An isolated 1H assignment at CH2 need not be the environment mean.
        hydrogen_ok = atom.GetTotalNumHs(includeNeighbors=True) != 2 or nh == 2
        if spec['nucleus'] == '13C' or hydrogen_ok:
            shift_target[j] = raw[i]['shift_ppm']; shift_valid[j] = True
    label_complete = bool(raw) and bool(known.all()) and not weak_rows.any() and not set_peaks
    observable_complete = bool(spec and spec.get('completeness', {}).get('status') == 'complete')
    mass_valid = label_complete and observable_complete and bool(len(groups)) and bool((target.sum(0) > 0).all())
    return raw, {'target_matrix': target.tolist(), 'pair_known_mask': known.tolist(),
                 'allowed_target': allowed.tolist(), 'allowed_row_mask': weak_rows.tolist(),
                 'row_transport_mask': row_transport.tolist(), 'shift_target': shift_target.tolist(),
                 'shift_known_mask': shift_valid.tolist(), 'mass_supervision_valid': mass_valid,
                 'label_complete_by_branch': label_complete, 'spectrum_complete': observable_complete,
                 'assignment_status_by_peak': label_status, 'assigned_atom_indices_by_peak': atom_labels,
                 'assignment_origin_by_peak': [(p.get('assignment') or {}).get('origin') for p in raw],
                 'ai_provenance_by_peak': [(p.get('assignment') or {}).get('ai') for p in raw],
                 'set_peak_membership': set_peaks, 'set_group_membership': set_groups,
                 'assignment_groups': (spec or {}).get('assignment_groups', []),
                 'target_resolution': 'carbon_atom' if (spec or {}).get('nucleus') == '13C' else 'attached_heavy_atom'}


def build_record(record):
    if not record.get('training_eligible', True) or record.get('curation', {}).get('molecule_state', 'active') != 'active':
        raise ValueError('Molecule is excluded from training')
    mol, graph, groups, connectivity = graph_and_candidates(record['smiles'], record.get('structure', {}).get('atoms'))
    chosen, pair = select_spectra(record)
    if not any(chosen.values()):
        raise ValueError('No eligible nonempty experimental spectrum')
    output = {'schema_version': SCHEMA_VERSION, 'id': record['id'], 'origin': record.get('origin'),
              'smiles_canonical': record['smiles'], 'molecular_formula': rdMolDescriptors.CalcMolFormula(mol),
              'connectivity_key': connectivity, 'molecular_graph': graph,
              'candidate_environments': {'carbon_groups': groups['carbon'], 'proton_environments': groups['proton']},
              'spectra': {}, 'conditions': {'by_nucleus': {}}, 'alignment_supervision': {}, 'branch_present': {},
              'provenance': {'molecule_revision': record.get('revision'), 'selected_spectrum_ids': {},
                'pair_id': pair['id'] if pair else None,
                'pairing_basis': pair['basis'] if pair else 'unpaired_independent_experiments',
                'same_known_solvent': bool(pair and pair['solvent_match'] == 'same_known'),
                'atom_index_basis': 'stored_canonical_smiles_index', 'grouping_policy': 'individual_atoms_v1'}}
    for branch in BRANCHES:
        spec = chosen[branch]
        raw, labels = branch_supervision(spec, groups[branch], mol)
        peaks = [observed_peak(p, NUCLEI[branch]) for p in raw]
        output['spectra'][branch] = {'peaks': peaks, 'metadata': conditions_for(spec)}
        output['c_nmr_peaks' if branch == 'carbon' else 'h_nmr_peaks'] = peaks
        output['alignment_supervision'][branch] = labels
        output['conditions']['by_nucleus'][NUCLEI[branch]] = conditions_for(spec)
        output['branch_present'][branch] = bool(peaks)
        output['provenance']['selected_spectrum_ids'][branch] = spec['id'] if spec else None
    return output


def feature_arrays(record):
    if record.get('schema_version') != SCHEMA_VERSION:
        raise ValueError(f'Expected {SCHEMA_VERSION}; re-export using the AtomPeak training format')
    graph, mapping = graph_arrays(record['molecular_graph'])
    result = {'graph': graph}
    for branch in BRANCHES:
        peaks = record['spectra'][branch]['peaks']
        labels = record['alignment_supervision'][branch]
        groups = record['candidate_environments']['carbon_groups' if branch == 'carbon' else 'proton_environments']
        p, g = len(peaks), len(groups)
        member = np.zeros((g, len(mapping)), dtype=np.float32)
        structural_mass = np.ones(g, dtype=np.float32)
        for j, group in enumerate(groups):
            for atom in group['center_atom_indices']:
                member[j, mapping[atom]] = 1
            structural_mass[j] = max(1, group['observable_h_count']) if branch == 'proton' else 1
        def matrix(name):
            value = np.asarray(labels[name], dtype=np.float32)
            if not p:
                return np.zeros((0, g), dtype=np.float32)
            if value.shape != (p, g):
                raise ValueError(f'{branch}.{name}: matrix shape does not match peak/atom order')
            return value
        result[branch] = {'peak_shifts': np.asarray([v['shift_ppm'] for v in peaks], dtype=np.float32),
            'peak_valid': np.ones(p, dtype=bool), 'peak_mass': np.ones(p, dtype=np.float32),
            'peak_confidence': np.full(p, 1. if branch == 'carbon' else .7, dtype=np.float32),
            'group_membership': member, 'group_valid': np.ones(g, dtype=bool),
            'group_mass': np.ones(g, dtype=np.float32), 'group_structural_mass': structural_mass,
            'target': matrix('target_matrix'), 'pair_known_mask': matrix('pair_known_mask'),
            'allowed_target': matrix('allowed_target'),
            'allowed_row_mask': np.asarray(labels['allowed_row_mask'] if p else [], dtype=bool),
            'row_transport_mask': np.asarray(labels['row_transport_mask'] if p else [], dtype=bool),
            'shift_target': np.asarray(labels['shift_target'], dtype=np.float32),
            'shift_known_mask': np.asarray(labels['shift_known_mask'], dtype=bool),
            'mass_supervision_valid': np.asarray(labels['mass_supervision_valid'] and p > 0),
            'branch_present': np.asarray(p > 0),
            'set_peak_membership': np.asarray(labels['set_peak_membership'] if p else [], dtype=np.float32).reshape(-1, p) if p else np.zeros((0, 0), dtype=np.float32),
            'set_group_membership': np.asarray(labels['set_group_membership'] if p else [], dtype=np.float32).reshape(-1, g) if g else np.zeros((0, 0), dtype=np.float32)}
    return result
