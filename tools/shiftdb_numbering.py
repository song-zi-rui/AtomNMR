"""Source-first numbering assessment; alternatives are never automatic labels.

The official visible-molfile basis remains the default. AuxInfo is a bounded
alternative hypothesis, not a redefinition of ShiftDB2's exported format.
"""
from __future__ import annotations

from dataclasses import replace
import hashlib
import json
import re

from rdkit.Chem import inchi
from tools.shiftdb_atom_identity import SOURCE_TYPE, source_peaks

POLICY = 'shiftdb-source-first-numbering-review-v1'
INCHI_BASIS = 'inchi_main_layer_canonical_index_0based_hypothesis'


def inchi_identity(identity, order):
    """Compose canonical InChI IDs with the already verified source/display map."""
    if not order or len(set(order))!=len(order):
        raise ValueError('AuxInfo order is empty or repeats a source atom')
    original={a.original_rdkit_idx:a for a in identity.atoms}
    if any(i not in original for i in order):
        raise ValueError('AuxInfo atom is absent from original molfile identity')
    return replace(identity,atoms=[replace(original[old],source_assignment_id=str(new))
                                  for new,old in enumerate(order)],assignment_id_basis=INCHI_BASIS)


def error_count(validation):
    return sum(i['severity']=='error' for i in validation['issues'])


def assess_numbering(value, nucleus, identity, *, source_inchi=None, record_conflict=False):
    peaks,_,validation=source_peaks(value,nucleus,identity,record_conflict=record_conflict)
    primary={'basis':SOURCE_TYPE,'validation':validation,'peaks':peaks,
             'mapping_has_hard_error':bool(error_count(validation)),
             'adoptable':validation['status']=='needs_attention' and error_count(validation)==0,
             'adoption_scope':'atom mapping only; raw multiplicity warnings remain unresolved'}
    primary['candidate_id']=hashlib.sha256(json.dumps({'policy':POLICY,'basis':SOURCE_TYPE,
        'raw':value,'nucleus':nucleus,'identity':identity.to_dict(),'record_conflict':record_conflict},
        sort_keys=True,separators=(',',':')).encode()).hexdigest()
    flow={'policy':POLICY,'default':primary,'fallback_attempted':False,'candidates':[],
          'automatic_adoption':False,
          'note':'Structural compatibility is not experimental assignment validation.'}
    if validation['status']=='identity_validated':return flow
    flow['fallback_attempted']=True
    try:
        generated,aux=inchi.MolToInchiAndAuxInfo(identity.source)
        found=re.search(r'/N:([^/]+)',aux)
        if not generated or found is None:raise ValueError('InChI/AuxInfo unavailable')
        components=[[int(n)-1 for n in part.split(',')] for part in found[1].split(';')]
        order=[n for component in components for n in component]
        candidate_identity=inchi_identity(identity,order)
        cp,_,cv=source_peaks(value,nucleus,candidate_identity)
        exact=bool(source_inchi and generated==source_inchi)
        changed=any(a.original_rdkit_idx!=int(a.source_assignment_id) for a in candidate_identity.atoms)
        reasons=[]
        if not exact:reasons.append('source_inchi_missing_or_not_exact_equal')
        if not changed:reasons.append('candidate_same_as_original_order')
        if cv['status']!='identity_validated':reasons.append('candidate_structural_checks_failed')
        if not cp or any(not f['assignments'] for p in cp for f in p['source_assignment_fragments']):
            reasons.append('candidate_has_empty_assignments')
        stamp=json.dumps({'policy':POLICY,'raw':value,'nucleus':nucleus,
            'generated_inchi':generated,'aux':aux,'mapping':candidate_identity.to_dict()},
            sort_keys=True,separators=(',',':'))
        flow['candidates']=[{'candidate_id':hashlib.sha256(stamp.encode()).hexdigest(),
            'basis':INCHI_BASIS,'generated_inchi':generated,'source_inchi_exact_equal':exact,
            'aux_info':aux,'canonical_to_original_source_index':order,
            'validation':cv,'peaks':cp,'identity':candidate_identity.to_dict(),
            'adoptable':not reasons,'blocked_reasons':reasons,
            'status':'awaiting_human_review' if not reasons else 'blocked',
            'note':'Passing checks does NOT confirm historical source numbering; independent source evidence required.'}]
    except Exception as exc:
        flow['fallback_error']=str(exc)
    return flow
