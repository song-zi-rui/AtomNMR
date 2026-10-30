"""Read-only source diagnostics, captured before any sanitization or reordering."""
from __future__ import annotations
import argparse
import io
import json
from pathlib import Path
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rdkit import Chem
from tools.filter_nmrshiftdb2_experimental import iter_records, parse_record, sha256


def atom_order(mol, *, original_order=False):
    return [{"original_molfile_position": (a.GetIntProp('_DiagnosticSourceIndex')+1 if a.HasProp('_DiagnosticSourceIndex') else a.GetIdx()+1 if original_order else None), "rdkit_idx": a.GetIdx(),
             "original_source_idx": a.GetIntProp('_DiagnosticSourceIndex') if a.HasProp('_DiagnosticSourceIndex') else a.GetIdx() if original_order else None,
             "element": a.GetSymbol(), "isotope": a.GetIsotope(),
             "atom_map_number": a.GetAtomMapNum(),
             "properties": a.GetPropsAsDict(includePrivate=True, includeComputed=False),
             "neighbors": [n.GetIdx() for n in a.GetNeighbors()]}
            for a in mol.GetAtoms()]


def diagnose(raw, index, line):
    block, props, _ = parse_record(raw)
    # No Sanitize, RemoveHs, AddHs, SMILES or drawing operation precedes this.
    original = next(Chem.ForwardSDMolSupplier(io.BytesIO(raw), sanitize=False, removeHs=False))
    if original is None:
        raise ValueError('RDKit cannot read this source record')
    result = {"record_index_1based": index, "source_start_line": line,
              "raw_molblock": block.decode(), "raw_properties": props,
              "rdkit_property_names": list(original.GetPropNames()),
              "immediate_unsanitized_atom_order": atom_order(original, original_order=True),
              "bonds": [{"begin_molfile_position": b.GetBeginAtomIdx()+1,
                         "end_molfile_position": b.GetEndAtomIdx()+1,
                         "order": str(b.GetBondType())} for b in original.GetBonds()]}
    work = Chem.Mol(original)
    for a in work.GetAtoms(): a.SetIntProp('_DiagnosticSourceIndex', a.GetIdx())
    result['copied_atom_order'] = atom_order(work)
    Chem.SanitizeMol(work)
    result['sanitized_atom_order'] = atom_order(work)
    expanded = Chem.AddHs(work)
    result['added_h_atom_order'] = atom_order(expanded)
    result['added_then_removed_h_atom_order'] = atom_order(Chem.RemoveHs(expanded))
    heavy = Chem.RemoveHs(work)
    result['removed_h_atom_order'] = atom_order(heavy)
    smi = Chem.MolToSmiles(heavy, canonical=True, isomericSmiles=True)
    result['canonical_smiles'] = smi
    result['smiles_output_order'] = json.loads(heavy.GetProp('_smilesAtomOutputOrder'))
    reconstructed = Chem.MolFromSmiles(smi)
    result['smiles_roundtrip_atom_order'] = atom_order(reconstructed)
    from annotation_tool_v3 import core
    _, _, mapping = core.canonicalize(work)
    result['existing_source_index_to_display_index'] = mapping
    for old,current in mapping.items():
        reconstructed.GetAtomWithIdx(current).SetIntProp('_DiagnosticSourceIndex',old)
    result['explicitly_remapped_atom_order'] = atom_order(reconstructed)
    from rdkit.Chem.Draw import rdMolDraw2D
    result['prepared_drawing_atom_order'] = atom_order(rdMolDraw2D.PrepareMolForDrawing(reconstructed))
    from tools.shiftdb_atom_identity import canonical_identity
    from tools.shiftdb_numbering import assess_numbering
    from tools.filter_nmrshiftdb2_experimental import SPECTRUM,classify
    identity=canonical_identity(original,f'shiftdb:{index}',sanitize=True)
    result['atom_identity']=identity.to_dict()
    result['spectrum_numbering']={field:assess_numbering(value,m[1],identity,source_inchi=props.get('INChI'))
        for field,value in props.items() if (m:=SPECTRUM.fullmatch(field)) and m[1] in ('13C','1H')
        and classify(props,m[2],m[1])[0]=='experimental'}
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--sdf',type=Path,default=Path('data/nmrshiftdb2/nmrshiftdb2withsignals.sd'))
    p.add_argument('--contains',default='185.084;0.0S;3')
    p.add_argument('--output-dir',type=Path,default=Path('data/shiftdb_atom_identity_audit'))
    p.add_argument('--quiet',action='store_true',help='Write full diagnostics but print only the locator')
    args=p.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    matches=[]
    for index,line,raw in iter_records(args.sdf):
        if args.contains.encode() not in raw: continue
        data=diagnose(raw,index,line)
        prefix=args.output_dir/f"diagnostic_record_{index}"
        prefix.with_suffix('.sdf').write_bytes(raw)
        prefix.with_suffix('.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
        print(json.dumps({'record':index,'source_id':data['raw_properties'].get('nmrshiftdb2 ID'),'file':str(prefix.with_suffix('.json'))} if args.quiet else data,ensure_ascii=False,indent=2),flush=True)
        matches.append(index)
    print(json.dumps({'matches':matches,'source_sha256':sha256(args.sdf)}),flush=True)


if __name__=='__main__':main()
