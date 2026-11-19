"""Map documented source atom indices to AtomNMR structure indices."""

from rdkit import Chem
from annotation_tool_v3 import core

MAX_MATCHES = 2048


def target_mappings(source, target, structure):
    if source is None:
        return [], {}
    # The application draws and labels this exact SMILES order. Check that the
    # stored atom table also describes it before accepting any computed targets.
    if structure['atoms'] != core.structure(target)['atoms']:
        return [], {}
    # Use only the declared stereochemistry for the graph match. A conformer
    # can add information about centers the source SMILES leaves unspecified.
    declared = Chem.Mol(source)
    declared.RemoveAllConformers()
    for atom in declared.GetAtoms():
        atom.SetIntProp('source_atom_index', atom.GetIdx())
    heavy = Chem.RemoveHs(declared)
    Chem.AssignStereochemistry(heavy, cleanIt=True, force=True)
    if Chem.MolToInchiKey(heavy) != Chem.MolToInchiKey(target):
        return [], {}
    matches = target.GetSubstructMatches(heavy, useChirality=True, uniquify=False, maxMatches=MAX_MATCHES)
    if not matches or len(matches) >= MAX_MATCHES or heavy.GetNumAtoms() != target.GetNumAtoms():
        return [], {}
    # RDKit substructure matching alone permits swapping aromatic nitrogens
    # with different H counts. The map must preserve the full atom identity.
    def signature(atom):
        return (atom.GetAtomicNum(), atom.GetIsotope(), atom.GetFormalCharge(),
                atom.GetIsAromatic(), atom.GetTotalNumHs())
    source_signatures = tuple(signature(a) for a in heavy.GetAtoms())
    target_signatures = tuple(signature(a) for a in target.GetAtoms())
    matches = [match for match in matches if tuple(target_signatures[i] for i in match) == source_signatures]
    if not matches:
        return [], {}
    source_indexes = [a.GetIntProp('source_atom_index') for a in heavy.GetAtoms()]
    mappings = [dict(zip(source_indexes, match)) for match in matches]
    attached = {}
    for atom in source.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 1 or neighbors[0].GetIdx() not in mappings[0]:
            return [], {}
        attached[atom.GetIdx()] = neighbors[0].GetIdx()
    return mappings, attached
