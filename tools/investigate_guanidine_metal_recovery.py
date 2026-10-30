"""Diagnostic-only recovery experiments; never write production assignments.

CXSMILES, original-order graph storage, and InChI canonical ordering are tested
as separate hypotheses. Passing a hypothesis does not authorize label transfer.
"""
from __future__ import annotations
import argparse
from collections import Counter
from dataclasses import replace
import json
from pathlib import Path
import re
import sys
if __package__ in (None,''):
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from rdkit import Chem,RDLogger
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import inchi
from tools.filter_nmrshiftdb2_experimental import iter_records,parse_record,classify,SPECTRUM,sha256
from tools.shiftdb_atom_identity import identity_for_current,source_peaks


def native_graph_signature(mol):
    """Compare tested atom/query/bond properties in stable SOURCE identity space.

    This is not a complete stereochemical equivalence test. No chemical
    normalisation, oxidation-state inference, or valence repair is performed.
    """
    ids={a.GetIdx():a.GetIntProp('_SourceDiagnosticIdx') for a in mol.GetAtoms()}
    atoms=sorted((ids[a.GetIdx()],a.GetAtomicNum(),a.GetIsotope(),a.GetFormalCharge(),
                  a.GetNumRadicalElectrons(),a.GetTotalNumHs(includeNeighbors=True),
                  a.DescribeQuery() if a.HasQuery() else '') for a in mol.GetAtoms())
    bonds=[]
    for b in mol.GetBonds():
        ends=(ids[b.GetBeginAtomIdx()],ids[b.GetEndAtomIdx()])
        if b.GetBondType()!=Chem.BondType.DATIVE:ends=tuple(sorted(ends))
        bonds.append((*ends,str(b.GetBondType()),b.DescribeQuery() if b.HasQuery() else ''))
    return atoms,sorted(bonds)


def native_storage_trial(mol):
    tagged=Chem.Mol(mol)
    for a in tagged.GetAtoms():a.SetIntProp('_SourceDiagnosticIdx',a.GetIdx())
    before=native_graph_signature(tagged)
    # Reorder the original graph, instead of recreating it from SMILES.
    reordered=Chem.RenumberAtoms(tagged,list(reversed(range(tagged.GetNumAtoms()))))
    block=Chem.MolToMolBlock(tagged)
    restored=Chem.MolFromMolBlock(block,sanitize=False,removeHs=False)
    if restored is not None:
        for a in restored.GetAtoms():a.SetIntProp('_SourceDiagnosticIdx',a.GetIdx())
        Chem.SanitizeMol(restored)
    return {'renumber_tested_properties_preserved':before==native_graph_signature(reordered),
            'molblock_tested_properties_preserved':restored is not None and before==native_graph_signature(restored),
            'query_bonds':[{'source_atoms':[b.GetBeginAtomIdx(),b.GetEndAtomIdx()],
                           'bond_type':str(b.GetBondType()),'query':b.DescribeQuery()}
                          for b in tagged.GetBonds() if b.HasQuery()],
            'chemical_readiness':'NOT established; retained valence and hydrogen counts may still be wrong',
            'scope':'atom element/isotope/charge/radicals/H capacity; bond type/dative direction/query predicate'}


def inchi_candidate_identity(source,order,index):
    """A candidate ID interpretation only; preserve actual molfile identities."""
    if len(set(order))!=len(order) or any(i<0 or i>=source.GetNumAtoms() for i in order):
        raise ValueError('Invalid AuxInfo original atom ordering')
    candidate=identity_for_current(source,dict(enumerate(order)),source_key=f'INCHI_HYPOTHESIS:{index}')
    candidate.atoms=[replace(a,original_rdkit_idx=a.processed_rdkit_idx,
                            molfile_position_1based=a.processed_rdkit_idx+1,
                            source_atom_id=f'shiftdb:{index}:atom:{a.processed_rdkit_idx+1}')
                     for a in candidate.atoms]
    candidate.assignment_id_basis='inchi_main_layer_canonical_index_0based_hypothesis'
    return candidate


def graph_roundtrip(heavy,kind):
    copy=Chem.Mol(heavy);copy.RemoveAllConformers()
    if kind=='cxsmiles':
        text=Chem.MolToCXSmiles(copy)
    else:text=Chem.MolToSmiles(copy)
    order=json.loads(copy.GetProp('_smilesAtomOutputOrder')) if copy.HasProp('_smilesAtomOutputOrder') else None
    params=Chem.SmilesParserParams();params.removeHs=False;params.parseName=False
    rebuilt=Chem.MolFromSmiles(text,params)
    if rebuilt is None:return {'parse_success':False,'text':text}
    if kind=='cxsmiles':
        original_to_heavy={a.GetProp('sourceAtomId'):a.GetIdx() for a in heavy.GetAtoms()}
        if any(not a.HasProp('sourceAtomId') for a in rebuilt.GetAtoms()):
            return {'parse_success':True,'text':text,'mapping_error':'CX sourceAtomId properties absent','graph_identity_preserved':False}
        old_to_new={original_to_heavy[a.GetProp('sourceAtomId')]:a.GetIdx() for a in rebuilt.GetAtoms()}
    else:old_to_new={old:new for new,old in enumerate(order)}
    differences=[]
    for old,new in old_to_new.items():
        a,b=heavy.GetAtomWithIdx(old),rebuilt.GetAtomWithIdx(new)
        for label,left,right in [('element',a.GetSymbol(),b.GetSymbol()),('isotope',a.GetIsotope(),b.GetIsotope()),
             ('charge',a.GetFormalCharge(),b.GetFormalCharge()),('radicals',a.GetNumRadicalElectrons(),b.GetNumRadicalElectrons()),
             ('atom_map',a.GetAtomMapNum(),b.GetAtomMapNum()),('hydrogen_capacity',a.GetTotalNumHs(includeNeighbors=True),b.GetTotalNumHs(includeNeighbors=True))]:
            if left!=right:differences.append({'kind':label,'source_atom_idx':a.GetIntProp('_SourceDiagnosticIdx'),'before':left,'after':right})
    for bond in heavy.GetBonds():
        a,b=old_to_new[bond.GetBeginAtomIdx()],old_to_new[bond.GetEndAtomIdx()]
        rb=rebuilt.GetBondBetweenAtoms(a,b)
        if rb is None or rb.GetBondType()!=bond.GetBondType():
            differences.append({'kind':'bond','source_atoms':[heavy.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetIntProp('_SourceDiagnosticIdx'),
                  heavy.GetAtomWithIdx(bond.GetEndAtomIdx()).GetIntProp('_SourceDiagnosticIdx')],
                'before':str(bond.GetBondType()),'after':str(rb.GetBondType()) if rb else None})
        elif bond.GetBondType()==Chem.BondType.DATIVE and (rb.GetBeginAtomIdx(),rb.GetEndAtomIdx())!=(a,b):
            differences.append({'kind':'dative_direction','source_atoms':[a,b]})
        if rb is not None and (bond.HasQuery()!=rb.HasQuery() or
              (bond.HasQuery() and bond.DescribeQuery()!=rb.DescribeQuery())):
            differences.append({'kind':'bond_query','source_atoms':[a,b],
                'before':bond.DescribeQuery() if bond.HasQuery() else '',
                'after':rb.DescribeQuery() if rb.HasQuery() else ''})
    properties_preserved=all(rebuilt.GetAtomWithIdx(new).HasProp('sourceAtomId') and
        rebuilt.GetAtomWithIdx(new).GetProp('sourceAtomId')==heavy.GetAtomWithIdx(old).GetProp('sourceAtomId') for old,new in old_to_new.items())
    return {'text':text,'parse_success':True,'n_atoms':rebuilt.GetNumAtoms(),'n_bonds':rebuilt.GetNumBonds(),
        'formula':rdMolDescriptors.CalcMolFormula(rebuilt),'difference_count':len(differences),
        'differences':differences,'source_atom_properties_preserved':properties_preserved,
        'graph_identity_preserved':not differences and rebuilt.GetNumAtoms()==heavy.GetNumAtoms() and rebuilt.GetNumBonds()==heavy.GetNumBonds(),
        'source_index_to_rebuilt_index':{str(heavy.GetAtomWithIdx(old).GetIntProp('_SourceDiagnosticIdx')):new for old,new in old_to_new.items()}}


def inspect_record(raw,index,category):
    block,props,_=parse_record(raw);original=Chem.MolFromMolBlock(block.decode(),sanitize=False,removeHs=False)
    for a in original.GetAtoms():a.SetIntProp('_SourceDiagnosticIdx',a.GetIdx());a.SetProp('sourceAtomId',str(a.GetIdx()))
    source=Chem.Mol(original);Chem.SanitizeMol(source);Chem.AssignStereochemistry(source,cleanIt=True,force=True)
    heavy=Chem.RemoveHs(source);Chem.AssignStereochemistry(heavy,cleanIt=True,force=True)
    result={'record_index':index,'category':category,'source_id':props.get('nmrshiftdb2 ID'),
        'properties':props,'source_formula_after_sanitize':rdMolDescriptors.CalcMolFormula(source),
        'original_atom_table':[{'source_idx':a.GetIdx(),'element':a.GetSymbol(),'formal_charge':a.GetFormalCharge(),
            'hydrogen_capacity':a.GetTotalNumHs(includeNeighbors=True),'neighbors':[n.GetIdx() for n in a.GetNeighbors()]} for a in source.GetAtoms()],
        'bond_types':dict(Counter(str(b.GetBondType()) for b in source.GetBonds())),
        'native_storage_trial':native_storage_trial(source),
        'source_order_assignments':{},'cxsmiles_trial':graph_roundtrip(heavy,'cxsmiles'),
        'plain_smiles_trial':graph_roundtrip(heavy,'plain')}
    source_identity=identity_for_current(source,source_key=f'shiftdb:{index}')
    # AuxInfo N gives canonical InChI ordering as 1-based ORIGINAL positions;
    # testing it is diagnostic, never a production interpretation of field 3.
    try:
        generated,aux=inchi.MolToInchiAndAuxInfo(source)
        m=re.search(r'/N:([^/]+)',aux)
        components=[[int(x)-1 for x in group.split(',')] for group in m[1].split(';')] if m else []
        order=[a for component in components for a in component]
        result['inchi_trial']={'generated_inchi':generated,'source_inchi':props.get('INChI'),
            'aux_info':aux,'canonical_to_original_source_index':order,'components':components,
            'inchi_exact_equal':generated==props.get('INChI'),'warning':'Hypothesis only; not a confirmed assignment numbering definition.'}
        lookup={new:old for new,old in enumerate(order)}
        alternate_identity=inchi_candidate_identity(source,order,index)
    except Exception as exc:
        result['inchi_trial']={'error':str(exc)};alternate_identity=None
    for field,value in props.items():
        m=SPECTRUM.fullmatch(field)
        if not m or m[1] not in ('13C','1H') or classify(props,m[2],m[1])[0]!='experimental':continue
        peaks,_,validation=source_peaks(value,m[1],source_identity)
        result['source_order_assignments'][field]={'raw':value,'validation':validation,'peaks':peaks}
        if alternate_identity is not None:
            ap,_,av=source_peaks(value,m[1],alternate_identity)
            result['source_order_assignments'][field]['inchi_hypothesis']={'validation':av,
                'references':[{'shift':p['shift_ppm'],'source_ids':[a['source_assignment_id'] for f in p['source_assignment_fragments'] for a in f['assignments']],
                    'original_source_indices':[a.get('original_atom_idx') for f in p['source_assignment_fragments'] for a in f['assignments']],
                    'candidate_elements':[a.get('source_element') for f in p['source_assignment_fragments'] for a in f['assignments']],
                    'actual_original_source_indices':[lookup.get(int(a['source_assignment_id'])) for f in p['source_assignment_fragments'] for a in f['assignments'] if a['source_assignment_id'].isdigit()]}
                     for p in ap]}
    return result


def investigate(source,affected,output):
    output.mkdir(parents=True,exist_ok=True);chosen={}
    for line in affected.read_text(encoding='utf-8').splitlines():
        r=json.loads(line);f=r['features']['flags']
        category='metal' if f['has_metal'] else 'guanidine' if f['guanidine'] else None
        if category:
            for s in r['spectra']:
                if s['record_index'] is not None:chosen[int(s['record_index'])]=category
    digest=sha256(source);results=[];RDLogger.DisableLog('rdApp.warning');RDLogger.DisableLog('rdApp.error')
    for index,_,raw in iter_records(source):
        if index not in chosen:continue
        try:
            record=inspect_record(raw,index,chosen[index])
        except Exception as exc:record={'record_index':index,'category':chosen[index],'error':str(exc)}
        results.append(record);(output/f'record_{index}.sdf').write_bytes(raw)
        (output/f'record_{index}.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
    report={'diagnostic_only':True,'source_sha256_unchanged':sha256(source)==digest,'source_sha256':digest,
        'records':results,'counts':dict(Counter(r['category'] for r in results)),
        'metal_cx_graph_preserved':sum(r['category']=='metal' and r.get('cxsmiles_trial',{}).get('graph_identity_preserved',False) for r in results)}
    report['metal_native_renumber_preserved']=sum(r['category']=='metal' and r.get('native_storage_trial',{}).get('renumber_tested_properties_preserved',False) for r in results)
    report['metal_native_molblock_preserved']=sum(r['category']=='metal' and r.get('native_storage_trial',{}).get('molblock_tested_properties_preserved',False) for r in results)
    (output/'recovery_trials.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='records'},ensure_ascii=False,indent=2))
    for r in results:
        print(r['record_index'],r['category'],r.get('source_id'),r.get('source_formula_after_sanitize'),
              'CX',r.get('cxsmiles_trial',{}).get('difference_count'), 'plain',r.get('plain_smiles_trial',{}).get('difference_count'),flush=True)
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sdf',type=Path,default=Path('data/nmrshiftdb2/nmrshiftdb2withsignals.sd'))
    p.add_argument('--affected',type=Path,default=Path('data/shiftdb_atom_identity_audit/exclusion_bias_20261003/affected_molecules.jsonl'))
    p.add_argument('--output-dir',type=Path,default=Path('data/shiftdb_atom_identity_audit/guanidine_metal_recovery_20261003'))
    args=p.parse_args();investigate(args.sdf,args.affected,args.output_dir)
