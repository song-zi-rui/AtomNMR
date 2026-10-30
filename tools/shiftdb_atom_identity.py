"""NMRShiftDB2 *withsignals.sd* identity and assignment semantics.

Official SVN r2665: AdminTasks.getSpecfile2 writes getMdlNumberVisible();
DBAtom.getMdlNumberVisible is zero-based original visible-molfile order.
For 1H, invisible database hydrogen atoms are converted to their attached
visible atom; user-visible explicit H atoms remain H atom references.
This module does not interpret custom UI labels, NMReDATA or carbon ordinals.
"""
from __future__ import annotations
import json
import math
import re
from collections import Counter
from dataclasses import asdict, dataclass
from rdkit import Chem

POLICY = 'shiftdb-visible-atom-identity-v1'
SOURCE_TYPE = 'shiftdb_visible_molfile_index_0based'


@dataclass(frozen=True)
class AtomMapping:
    source_assignment_id: str
    source_atom_id: str
    molfile_position_1based: int
    original_rdkit_idx: int
    processed_rdkit_idx: int | None
    attached_heavy_processed_idx: int | None
    element: str
    isotope: int
    source_atom_map_number: int
    display_atom_label: str | None


@dataclass
class AtomIdentity:
    source: Chem.Mol
    processed: Chem.Mol
    smiles: str
    atoms: list[AtomMapping]
    source_key: str
    processed_index_basis: str = 'zero_based_canonical_smiles'
    assignment_id_basis: str = SOURCE_TYPE

    @property
    def legacy_mapping(self):
        return {a.original_rdkit_idx: (a.attached_heavy_processed_idx if a.element == 'H'
                else a.processed_rdkit_idx) for a in self.atoms
                if (a.attached_heavy_processed_idx if a.element == 'H' else a.processed_rdkit_idx) is not None}

    def to_dict(self):
        return {'policy': POLICY, 'source_assignment_type': self.assignment_id_basis,
                'source_key': self.source_key, 'processed_index_basis': self.processed_index_basis,
                'atoms': [asdict(a) for a in self.atoms]}


def canonical_identity(original, source_key='source', *, sanitize=False):
    """Tag before preprocessing; carry exact writer traversal across SMILES.

    No isomorphism is chosen, so symmetric atoms do not get arbitrary matches.
    Tags stay on a copy. The displayed SMILES and atom-map numbers are unchanged.
    """
    source = Chem.Mol(original)
    for a in source.GetAtoms(): a.SetIntProp('_OriginalSDFIndex', a.GetIdx())
    if sanitize: Chem.SanitizeMol(source)
    Chem.AssignStereochemistry(source, cleanIt=True, force=True)
    heavy = Chem.RemoveHs(source)
    # CIP properties copied across RemoveHs can be stale (especially rings).
    # Recompute on the actual processed graph before comparing with a SMILES
    # reparse; comparing cached R/S labels would falsely reject good structures.
    Chem.AssignStereochemistry(heavy, cleanIt=True, force=True)
    smi = Chem.MolToSmiles(heavy, canonical=True, isomericSmiles=True)
    order = json.loads(heavy.GetProp('_smilesAtomOutputOrder'))
    parser = Chem.SmilesParserParams()
    parser.removeHs = False  # Keep isotopic / user-visible explicit H identities.
    processed = Chem.MolFromSmiles(smi, parser)
    if processed is None or processed.GetNumAtoms() != len(order):
        raise ValueError('parser_error: canonical SMILES atom count changed')
    mapping = {heavy.GetAtomWithIdx(old).GetIntProp('_OriginalSDFIndex'): new
               for new, old in enumerate(order)}
    heavy_to_new = {old:new for new,old in enumerate(order)}
    for old,new in heavy_to_new.items():
        a,b=heavy.GetAtomWithIdx(old),processed.GetAtomWithIdx(new)
        if (a.GetAtomicNum(),a.GetIsotope(),a.GetFormalCharge(),a.GetAtomMapNum()) != (
                b.GetAtomicNum(),b.GetIsotope(),b.GetFormalCharge(),b.GetAtomMapNum()):
            raise ValueError('parser_error: canonical SMILES atom identity changed')
        if a.HasProp('_CIPCode') and (not b.HasProp('_CIPCode') or a.GetProp('_CIPCode')!=b.GetProp('_CIPCode')):
            raise ValueError('parser_error: canonical SMILES atom stereochemistry changed')
        b.SetIntProp('_OriginalSDFIndex', a.GetIntProp('_OriginalSDFIndex'))
    expected={(min(heavy_to_new[b.GetBeginAtomIdx()],heavy_to_new[b.GetEndAtomIdx()]),
               max(heavy_to_new[b.GetBeginAtomIdx()],heavy_to_new[b.GetEndAtomIdx()]),str(b.GetBondType()))
              for b in heavy.GetBonds()}
    actual={(min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),max(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),str(b.GetBondType()))
            for b in processed.GetBonds()}
    if actual!=expected: raise ValueError('parser_error: canonical SMILES connectivity changed')
    stereo_name={Chem.BondStereo.STEREOE:'E',Chem.BondStereo.STEREOTRANS:'E',
                 Chem.BondStereo.STEREOZ:'Z',Chem.BondStereo.STEREOCIS:'Z'}
    for bond in heavy.GetBonds():
        name=stereo_name.get(bond.GetStereo())
        rebuilt=processed.GetBondBetweenAtoms(heavy_to_new[bond.GetBeginAtomIdx()],heavy_to_new[bond.GetEndAtomIdx()])
        if name and stereo_name.get(rebuilt.GetStereo())!=name:
            raise ValueError('parser_error: canonical SMILES bond stereochemistry changed')
    identities=[]
    for a in source.GetAtoms():
        idx=a.GetIdx(); target=mapping.get(idx); parent=target
        if a.GetSymbol()=='H':
            neighbors=list(a.GetNeighbors())
            parent=mapping.get(neighbors[0].GetIdx()) if len(neighbors)==1 and neighbors[0].GetSymbol()!='H' else None
        identities.append(AtomMapping(str(idx),f'{source_key}:atom:{idx+1}',idx+1,idx,target,parent,
            a.GetSymbol(),a.GetIsotope(),a.GetAtomMapNum(),
            f'{processed.GetAtomWithIdx(parent).GetSymbol()}{parent}' if parent is not None else None))
    return AtomIdentity(source,processed,smi,identities,source_key)


def identity_for_current(mol, mapping=None, source_key='source'):
    """Compatibility for callers whose molecule is already in source order.

    An explicit old->new mapping must still use original full atom order.
    This fallback cannot validate removed explicit H; production uses identity.
    """
    mapping=mapping if mapping is not None else {a.GetIdx():a.GetIdx() for a in mol.GetAtoms()}
    atoms=[]
    for old,new in sorted(mapping.items()):
        a=mol.GetAtomWithIdx(new)
        parent=new
        if a.GetSymbol()=='H':
            ns=list(a.GetNeighbors())
            parent=ns[0].GetIdx() if len(ns)==1 and ns[0].GetSymbol()!='H' else None
        atoms.append(AtomMapping(str(old),f'{source_key}:atom:{old+1}',old+1,old,new,parent,
                                a.GetSymbol(),a.GetIsotope(),a.GetAtomMapNum(),
                                f'{mol.GetAtomWithIdx(parent).GetSymbol()}{parent}' if parent is not None else None))
    return AtomIdentity(mol,mol,'',atoms,source_key,'zero_based_current_molecule')


def parse_fragments(value, nucleus, identity):
    """Keep one fragment per source occurrence, with source and display identities."""
    lookup={a.source_assignment_id:a for a in identity.atoms}
    fragments=[]
    for ordinal,raw in enumerate(value.split('|')):
        if not raw.strip(): continue
        fields=raw.split(';'); issues=[]
        try:
            shift=float(fields[0])
            if not math.isfinite(shift): raise ValueError()
        except (ValueError,IndexError):
            fragments.append({'source_fragment_ordinal':ordinal,'raw_signal':raw,'shift_ppm':None,
                              'assignments':[],'issues':[{'code':'invalid_shift','severity':'error'}]})
            continue
        ids=re.split(r'[, _]+',fields[2].strip()) if len(fields)==3 and fields[2].strip() else []
        if len(fields)!=3 or not ids: issues.append({'code':'missing_or_malformed_assignment','severity':'error'})
        assignments=[]
        for token in ids:
            entry={'source_assignment_id':token,'source_assignment_type':identity.assignment_id_basis,
                   'original_atom_idx':None,'processed_atom_idx':None,'display_atom_label':None,
                   'transfer_status':'unresolved','confidence':'unresolved'}
            try:
                key=str(int(token))
                if int(token)<0: raise ValueError()
            except ValueError:
                issues.append({'code':'invalid_assignment_id','severity':'error','source_assignment_id':token})
                assignments.append(entry);continue
            a=lookup.get(key)
            if a is None:
                issues.append({'code':'source_atom_out_of_range','severity':'error','source_assignment_id':token})
                assignments.append(entry);continue
            entry.update(original_atom_idx=a.original_rdkit_idx,source_atom_id=a.source_atom_id,
                         molfile_position_1based=a.molfile_position_1based,source_element=a.element,
                         source_isotope=a.isotope,processed_atom_idx=a.processed_rdkit_idx,
                         attached_heavy_processed_idx=a.attached_heavy_processed_idx,
                         display_atom_label=a.display_atom_label,
                         reference_resolution='visible_explicit_h' if a.element=='H' else 'visible_atom')
            target=a.attached_heavy_processed_idx if nucleus=='1H' else a.processed_rdkit_idx
            error=None
            if nucleus=='13C' and a.element!='C': error='carbon_assignment_to_noncarbon'
            elif nucleus=='1H' and a.element=='H' and a.isotope not in (0,1): error='proton_assignment_to_nonprotium'
            elif target is None: error='unmapped_source_atom'
            elif nucleus=='1H' and identity.processed.GetAtomWithIdx(target).GetTotalNumHs(includeNeighbors=True)<=0:
                error='proton_assignment_without_h_capacity'
            if error:
                issues.append({'code':error,'severity':'error','source_assignment_id':token,
                               'original_atom_idx':a.original_rdkit_idx,'source_element':a.element})
                entry['transfer_status']='source_assignment_inconsistency'
            else:
                entry.update(target_atom_idx=target,transfer_status='mapped',
                    confidence='identity_verified' if identity.assignment_id_basis==SOURCE_TYPE else 'numbering_hypothesis_only')
                if nucleus=='13C' and len(fields)>1:
                    code=re.search(r'([SDTQ])$',fields[1].strip())
                    if code and identity.processed.GetAtomWithIdx(target).GetTotalNumHs(includeNeighbors=True)!={'S':0,'D':1,'T':2,'Q':3}[code[1]]:
                        issues.append({'code':'carbon_hcount_code_conflict','severity':'warning','source_assignment_id':token})
            assignments.append(entry)
        fragments.append({'source_fragment_ordinal':ordinal,'raw_signal':raw,'shift_ppm':shift,
                          'shift_text':fields[0].strip(),'source_signal_code':fields[1] if len(fields)>1 else None,
                          'assignments':assignments,'issues':issues})
    return fragments


def validate_fragments(fragments, nucleus):
    issues=[dict(issue,source_fragment_ordinal=f['source_fragment_ordinal']) for f in fragments for issue in f['issues']]
    # Carbon S/D/T/Q conventionally corresponds to H class, but a contributor
    # can enter the multiplicity. Disagreement is metadata warning, not alone
    # proof of bad atom identity (nor the 1H peak multiplicity).
    status='source_assignment_inconsistency' if any(i['severity']=='error' for i in issues) else 'needs_attention' if issues else 'identity_validated'
    return {'policy':POLICY,'status':status,'issues':issues,'fragment_count':len(fragments),
            'issue_counts':dict(Counter(i['code'] for i in issues)),
            'identity_validated_is_not_chemical_assignment_validation':True}


def source_peaks(value, nucleus, identity, *, record_conflict=False):
    """Source occurrence IDs are not discarded by equal ppm/code aggregation.

    A UI row is an exact-shift bucket, not a claim of a unique physical signal.
    Original proton occurrences, including identical ones, stay separately stored.
    """
    fragments=parse_fragments(value,nucleus,identity)
    validation=validate_fragments(fragments,nucleus)
    if record_conflict and validation['status']=='identity_validated':
        validation['status']='record_assignment_conflict'
    grouped={};peaks=[];labels={}
    for fragment in fragments:
        if fragment['shift_ppm'] is None: continue
        key=(fragment['shift_text'],fragment['source_signal_code'])
        if key not in grouped:
            pid=('C' if nucleus=='13C' else 'H')+str(fragment['source_fragment_ordinal']+1).zfill(3)
            p={'peak_id':pid,'shift_ppm':fragment['shift_ppm'],'raw_signal':'',
               'raw_source_atom_indices':[],'source_signal_code':fragment['source_signal_code'],
               'multiplicity':None,'j_values':None,'reported_nH':None,'source_assignment_fragments':[],
               'signal_grouping':'exact_shift_display_bucket_not_signal_identity',
               'source_assignment_validation':validation['status']}
            grouped[key]=p;peaks.append(p)
        p=grouped[key];p['raw_signal']+=('|' if p['raw_signal'] else '')+fragment['raw_signal']
        p['raw_source_atom_indices'].append(fragment['raw_signal'].split(';')[-1])
        p['source_assignment_fragments'].append(fragment)
    # If a spectrum has hard identity conflicts, do not transfer a valid-looking
    # subset. It could share the same lost source numbering as the invalid part.
    # Alternative numbering candidates must never create imported labels merely
    # because their structural checks pass. They need explicit human adoption.
    blocked=validation['status']!='identity_validated' or identity.assignment_id_basis!=SOURCE_TYPE
    if not blocked:
        for p in peaks:
            targets=sorted({a['target_atom_idx'] for f in p['source_assignment_fragments'] for a in f['assignments'] if a['transfer_status']=='mapped'})
            if targets: labels[p['peak_id']]={'status':'assigned' if len(targets)==1 else 'source_group',
                'atom_indices':targets,'note':'原始可见原子身份经显式映射；同位移展示分组不证明等价、重叠或积分。'}
    return peaks,labels,validation
