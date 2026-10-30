"""Versioned, read-only export admission based on current effective spectra."""
from collections import Counter

DEFAULTS={'exclude_full':False,'exclude_partial':False,'exclude_unassigned':False,
          'equivalence':'none','require_dual':False,'require_same_solvent':False}
HARD={'assigned','equivalent_group','overlap'}

def validate(options=None):
    if options is None:options={}
    if not isinstance(options,dict) or set(options)-set(DEFAULTS):raise ValueError('导出筛选参数无效')
    value=DEFAULTS|options
    for k in DEFAULTS.keys()-{'equivalence'}:
        if type(value[k]) is not bool:raise ValueError('导出筛选开关应为布尔值')
    if value['equivalence'] not in ('none','carbon_missing','hc_missing','any_incomplete_group'):
        raise ValueError('等价原子筛选模式无效')
    if all(value[k] for k in ('exclude_full','exclude_partial','exclude_unassigned')):
        raise ValueError('至少保留一种归属完整度')
    return value

def current_spectra(record):
    return [s for s in record.get('spectra',[]) if s.get('manual_state','active')=='active'
            and s.get('quality',{}).get('decision')!='quarantined']

def coverage(spectra):
    peaks=[p for s in spectra for p in s['peaks']]
    assigned=sum(p.get('assignment',{}).get('status') in HARD
                 and bool(p['assignment'].get('atom_indices'))
                 and p['assignment'].get('source_identity_eligible',True) for p in peaks)
    return 'none' if not assigned else 'full' if assigned==len(peaks) else 'partial'

def equivalence_flags(record,spectra,mode):
    if mode=='none':return []
    groups=record.get('structure',{}).get('equivalent_groups',[])
    elements={a['index']:a['symbol'] for a in record.get('structure',{}).get('atoms',[])}
    flags=[]
    for s in spectra:
        if mode=='carbon_missing' and s['nucleus']!='13C':continue
        occupied={i for p in s['peaks'] for i in p.get('assignment',{}).get('atom_indices',[])}
        for p in s['peaks']:
            label=p.get('assignment',{})
            if label.get('origin')!='nmrshiftdb2':continue
            atoms=set(label.get('atom_indices',[]))
            if not atoms:continue
            for group in groups:
                g=set(group)
                if s['nucleus']=='13C' and any(elements.get(i)!='C' for i in g):continue
                if atoms<g and (g-occupied or mode=='any_incomplete_group'):
                    flags.append({'spectrum_id':s['id'],'peak_id':p['peak_id'],'nucleus':s['nucleus'],
                                  'source_atoms':sorted(atoms),'graph_group':sorted(g),'missing_atoms':sorted(g-occupied)})
    return flags

class Policy:
    def __init__(self,options=None):
        self.options=validate(options);self.scanned=0;self.kept=0;self.reasons=Counter();self.classes=Counter()
    def accept(self,record):
        self.scanned+=1;spectra=current_spectra(record);level=coverage(spectra);self.classes[level]+=1
        reasons=[];opt=self.options
        if not record.get('training_eligible',True) or not spectra:reasons.append('no_eligible_spectra')
        if opt[{'full':'exclude_full','partial':'exclude_partial','none':'exclude_unassigned'}[level]]:
            reasons.append('coverage_'+level)
        nuclei={s['nucleus'] for s in spectra if s['peaks']}
        if opt['require_dual'] and not {'13C','1H'}<=nuclei:reasons.append('missing_paired_nucleus')
        if opt['require_same_solvent']:
            ids={s['id'] for s in spectra}
            if not any(p['c_id'] in ids and p['h_id'] in ids and p.get('solvent_match')=='same_known' for p in record.get('pairs',[])):
                reasons.append('no_same_known_solvent_pair')
        if equivalence_flags(record,spectra,opt['equivalence']):reasons.append('suspected_source_equivalence')
        if reasons:self.reasons.update(reasons)
        else:self.kept+=1
        return not reasons
    def report(self):
        return {'policy':'export-admission-v1','options':self.options,'scanned':self.scanned,'kept':self.kept,
                'excluded':self.scanned-self.kept,'coverage_counts':dict(self.classes),'reasons':dict(self.reasons),
                'scope':'whole_molecule_all_current_eligible_spectra','reason_counts_may_overlap':True,
                'equivalence_is_suspected_not_confirmed':True}
