"""Conservative review flags for definite AI peak assignments.

These are consistency concerns, not proof that a source spectrum is wrong.
Conformers, slow exchange and other experimental circumstances can require
expert interpretation. This module never changes predictions or references.
"""
from collections import defaultdict


_HARD=frozenset(('assigned','equivalent_group','overlap'))
_GROUP=frozenset(('equivalent_group','overlap'))


def review_issues(output,context):
    """Return per-peak review issues; callers decide how to gate publication.

    Carbon reuse is checked within each 13C spectrum. Proton reuse is checked
    against the attached-H capacity of the corresponding heavy atom, so two
    separate proton resonances on a CH2 are allowed. Different spectra remain
    independent, and ambiguous candidate atoms never consume capacity.
    """
    atoms={a.get('index',i):a for i,a in enumerate(context.get('atoms',[]))
           if type(a.get('index',i)) is int}
    spectra={s['spectrum_id']:s for s in context.get('spectra',[]) if 'spectrum_id' in s}
    usage=defaultdict(lambda:defaultdict(set));issues={}

    def add(sid,pid,code,indices,reason):
        key=(sid,pid,code)
        if key not in issues:
            issues[key]={'spectrum_id':sid,'peak_id':pid,'code':code,
                         'atom_indices':[],'reason':reason}
        issues[key]['atom_indices']=sorted(set(issues[key]['atom_indices'])|set(indices))

    def hydrogens(index):
        number=atoms[index].get('hydrogens')
        return number if type(number) is int and number>=0 else None

    for row in output.get('assignments',[]):
        if not isinstance(row,dict) or row.get('status') not in _HARD:continue
        sid,pid=row.get('spectrum_id'),row.get('peak_id')
        if sid not in spectra or not isinstance(pid,str):continue
        nucleus=spectra[sid].get('nucleus')
        if nucleus not in ('13C','1H'):continue
        indices=row.get('atom_indices',[])
        if not isinstance(indices,list):continue
        # Invalid atom IDs have a separate validator; do not invent a capacity
        # for an atom which is absent from the given structure.
        indices=sorted({i for i in indices if type(i) is int and i in atoms})
        for index in indices:usage[sid][index].add(pid)
        if row['status'] in _GROUP:
            if nucleus=='13C':nuclei=len(indices)
            else:
                capacities=[hydrogens(i) for i in indices]
                nuclei=sum(capacities) if all(n is not None for n in capacities) else None
            if nuclei is not None and nuclei<2:
                add(sid,pid,'group_nucleus_count_review',indices,
                    '此确定的等价组/重叠组按结构仅涉及'+str(nuclei)+'个被观测核；请复核组状态含义与原子归属。')

    for sid,by_atom in usage.items():
        nucleus=spectra[sid]['nucleus']
        for index,peaks in by_atom.items():
            if nucleus=='13C' and len(peaks)>1:
                for pid in sorted(peaks):
                    add(sid,pid,'carbon_atom_multi_peak_review',[index],
                        '同一13C谱中这些原子被确定分配给多个不同峰；请复核归属，或核实是否存在多构象、慢交换等实验情况。')
            elif nucleus=='1H':
                capacity=hydrogens(index)
                if capacity is not None and len(peaks)>capacity:
                    for pid in sorted(peaks):
                        add(sid,pid,'proton_capacity_review',[index],
                            '同一1H谱中这些重原子的确定归属峰数超过连接氢数；请复核归属，或核实是否存在多构象、慢交换等实验情况。')
    return [issues[key] for key in sorted(issues)]
