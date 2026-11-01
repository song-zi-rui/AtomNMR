"""Targeted AI reannotation of source group omissions; never expands by graph alone."""
import hashlib
import json
from collections import defaultdict

from rdkit import Chem
from . import core as db

POLICY='shiftdb-equivalence-repair-v1'
PROMPT='''本批复查原库可能遗漏的等价原子组成员。targets是需要独立重标的峰，旧的待修复归属没有提供给你。
equivalence_hints只是结构图对称筛查提示，不是NMR等价证据；必须结合整个结构、已指定立体化学及全部可用C/H谱判断。
确有多个等价原子共同对应一个峰时返回equivalent_group及完整集合；不能确认则返回ambiguous或unassigned，不能为补全而强行合并。
已有anchors不在本次改标范围内。只有一种核种时明确信息限制，不虚构另一种谱或二维数据。'''


def prepare(record):
    from . import ai_annotation as ai
    context,refs=ai.prepare(record,require_dual=False)
    mol=Chem.MolFromSmiles(record['smiles'])
    ranks=Chem.CanonicalRankAtoms(mol,breakTies=False,includeChirality=True)
    by_rank=defaultdict(list)
    for atom in mol.GetAtoms():by_rank[ranks[atom.GetIdx()]].append(atom.GetIdx())
    groups=[g for g in by_rank.values() if len(g)>1]
    wanted={};anchors=[]
    for s in record['spectra']:
        if s['id'] not in {v['spectrum_id'] for v in context['spectra']}:continue
        for p in s['peaks']:
            label=p.get('assignment') or {};atoms=set(label.get('atom_indices') or [])
            group=next((g for g in groups if atoms and atoms<set(g)
                and all(mol.GetAtomWithIdx(i).GetSymbol()=='C' if s['nucleus']=='13C'
                        else mol.GetAtomWithIdx(i).GetTotalNumHs(includeNeighbors=True)>0 for i in g)),None)
            key=(s['id'],p['peak_id'])
            if label.get('origin')=='nmrshiftdb2' and group and not p.get('assignment_group'):
                wanted[key]={'before':{k:label.get(k) for k in ('status','atom_indices','origin','updated_at','note')},'graph_group':group}
            elif label.get('status') in ai.HARD and atoms and label.get('origin'):
                anchors.append({'spectrum_id':s['id'],'peak_id':p['peak_id'],'status':label['status'],'atom_indices':sorted(atoms)})
    if not wanted:raise ValueError('当前没有仍需复查的原库等价组缺成员峰；人工及后续AI标签保留')
    context['requested_peaks']=[{'spectrum_id':s,'peak_id':p} for s,p in wanted]
    context['anchors']=anchors
    context['equivalence_hints']=[{'spectrum_id':s,'peak_id':p,'graph_atom_group':v['graph_group'],'meaning':'graph_symmetry_only_not_NMR_proof'} for (s,p),v in wanted.items()]
    for ref in refs:
        key=(ref['spectrum_id'],ref.get('peak_id'))
        if key in wanted:ref['repair']=wanted[key]
    if len(db.dumps(context).encode())>220000:raise ValueError('复查上下文过大，未截断谱表')
    return context,refs


def config(value):
    from .ai_annotation import SYSTEM
    prompt=value.get('system_prompt',SYSTEM)+'\n'+PROMPT
    return dict(value,repair_policy=POLICY,provider_parallelism_cap=4,
                system_prompt=prompt,prompt_sha256=hashlib.sha256(prompt.encode()).hexdigest())


def replacement_reason(peak,old,prediction,reference,spectrum,context):
    """None permits only a high-confidence source-group expansion."""
    repair=reference.get('repair') if reference else None
    if not repair or not old or old['origin']!='nmrshiftdb2':return 'existing_label_preserved'
    before=repair['before']
    if (old['status']!=before['status'] or json.loads(old['atoms_json'])!=before['atom_indices']
            or old['updated_at']!=before['updated_at']):return 'source_label_changed'
    if prediction['status']!='equivalent_group' or prediction['confidence']!='high':return 'repair_not_high_confidence_equivalent_group'
    original=set(before['atom_indices']);new=set(prediction['atom_indices'])
    if not original<new or not new<=set(repair['graph_group']):return 'repair_not_supported_group_expansion'
    # Check the actual final state, including existing labels not requested from AI.
    from .ai_consistency import review_issues
    projected=[]
    for p in spectrum['peaks']:
        label=prediction if p['peak_id']==peak['peak_id'] else p.get('assignment',{})
        status=label.get('status')
        if status=='source_group':status='overlap'  # Preserve its occupied atoms for conservative checks.
        projected.append({'spectrum_id':spectrum['id'],'peak_id':p['peak_id'],'status':status,'atom_indices':label.get('atom_indices',[])})
    if any(r['peak_id']==peak['peak_id'] for r in review_issues({'assignments':projected},context)):
        return 'repair_conflicts_with_existing_assignments'
    return None
