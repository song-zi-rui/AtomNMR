"""Versioned, read-only export admission based on current effective spectra."""
from collections import Counter
from . import computed_spectra

DEFAULTS={'exclude_full':False,'exclude_partial':False,'exclude_unassigned':False,
          'equivalence':'none','require_dual':False,'require_same_solvent':False,
          'spectrum_type':'experimental','spectrum_source':'','coverage':'all','nuclei':'all'}
HARD={'assigned','equivalent_group','overlap'}
CHOICES={'spectrum_type':('experimental','computed'),
         'coverage':('all','full','partial','none'),
         'nuclei':('all','dual','proton_only','carbon_only'),
         'equivalence':('none','carbon_missing','hc_missing','any_incomplete_group')}


def validate(options=None):
    if options is None:options={}
    if not isinstance(options,dict) or set(options)-set(DEFAULTS):raise ValueError('导出筛选参数无效')
    value=DEFAULTS|options
    for k in DEFAULTS.keys()-CHOICES.keys()-{'spectrum_source'}:
        if type(value[k]) is not bool:raise ValueError('导出筛选开关应为布尔值')
    if any(value[k] not in allowed for k, allowed in CHOICES.items()):
        raise ValueError('谱类型、核种、归属筛选或等价原子选项无效')
    if not isinstance(value['spectrum_source'],str) or len(value['spectrum_source'])>100:
        raise ValueError('谱来源无效')
    if value['require_same_solvent'] and (value['spectrum_type']=='computed' or
                                        value['nuclei'] in ('proton_only','carbon_only')):
        raise ValueError('已知同溶剂谱对条件仅适用于含双谱的实验数据')
    if all(value[k] for k in ('exclude_full','exclude_partial','exclude_unassigned')):
        raise ValueError('至少保留一种归属完整度')
    return value

def current_spectra(record):
    return [s for s in record.get('spectra',[]) if s.get('manual_state','active')=='active'
            and s.get('quality',{}).get('decision')!='quarantined']

def coverage(spectra):
    peaks=[p for s in spectra for p in s['peaks']]
    assigned=sum(hard_assignment(p) for p in peaks)
    return 'none' if not assigned else 'full' if assigned==len(peaks) else 'partial'


def hard_assignment(peak):
    label=peak.get('assignment') or {}
    return (label.get('status') in HARD and bool(label.get('atom_indices'))
            and label.get('source_identity_eligible',True))


def nuclear_class(spectra):
    nuclei={s['nucleus'] for s in spectra if s['peaks']}
    return {frozenset({'13C','1H'}):'dual',frozenset({'13C'}):'carbon_only',
            frozenset({'1H'}):'proton_only'}.get(frozenset(nuclei),'none')

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
        self.nuclear_classes=Counter()

    def accept(self,record):
        return self.select(record) is not None

    def select(self,record,*,single_view=False):
        self.scanned+=1
        opt=self.options
        kind,source=opt['spectrum_type'],opt['spectrum_source']
        spectra=[s for s in current_spectra(record) if s['peaks']
                 and s.get('spectrum_kind','experimental')==kind and (not source or s['source']==source)]
        if single_view:
            spectra=[s for s in spectra if s.get('provenance',{}).get('structure_match')!='tautomer_variant']
        ids={s['id'] for s in spectra}
        pairs=[p for p in record.get('pairs',[]) if p['c_id'] in ids and p['h_id'] in ids]
        if kind=='computed':
            pairs=computed_spectra.pairs(spectra)
        if opt['require_same_solvent']:
            pairs=[p for p in pairs if p.get('solvent_match')=='same_known']
        projected={**record,'spectra':spectra,'pairs':pairs}
        if single_view:
            from atompeak_data.schema import select_spectra
            chosen,_=select_spectra(projected)
            used=[s for s in chosen.values() if s]
        else:
            used=spectra
        level,nuclei=coverage(used),nuclear_class(used)
        self.classes[level]+=1
        self.nuclear_classes[nuclei]+=1
        reasons=[]
        active=record.get('curation',{}).get('molecule_state')=='active'
        allowed=record.get('training_eligible',True) or (kind=='computed' and active)
        if not allowed or not used:reasons.append('no_eligible_spectra')
        if ((opt['coverage']!='all' and opt['coverage']!=level)
            or opt[{'full':'exclude_full','partial':'exclude_partial','none':'exclude_unassigned'}[level]]):
            reasons.append('coverage_'+level)
        if opt['require_dual'] and nuclei!='dual':reasons.append('missing_paired_nucleus')
        if opt['nuclei']!='all' and opt['nuclei']!=nuclei:reasons.append('nuclei_'+nuclei)
        if opt['require_same_solvent'] and not pairs:reasons.append('no_same_known_solvent_pair')
        if equivalence_flags(record,used,opt['equivalence']):reasons.append('suspected_source_equivalence')
        if reasons:
            self.reasons.update(reasons)
            return None
        self.kept+=1
        projected.update(training_eligible=True, export_classification={'spectrum_type':kind,'coverage':level,'nuclei':nuclei})
        peaks=[p for s in used for p in s['peaks']]
        projected['annotation_coverage']={'status':level,'assigned_peaks':sum(map(hard_assignment,peaks)),
                                         'total_peaks':len(peaks),'basis':'exported_spectra_determined_assignments'}
        for nucleus,prefix in (('13C','c'),('1H','h')):
            branch=[p for s in used if s['nucleus']==nucleus for p in s['peaks']]
            projected[prefix+'_count']=len(branch)
            projected[prefix+'_assigned']=sum(map(hard_assignment,branch))
        projected['annotation_sources']=sorted({p['assignment']['origin'] for p in peaks if p['assignment'].get('origin')})
        return projected

    def report(self):
        return {'policy':'export-admission-v2','options':self.options,'scanned':self.scanned,'kept':self.kept,
                'excluded':self.scanned-self.kept,'coverage_counts':dict(self.classes),'reasons':dict(self.reasons),
                'nuclei_counts':dict(self.nuclear_classes),'scope':'current_spectra_from_selected_type_and_source',
                'reason_counts_may_overlap':True,
                'equivalence_is_suspected_not_confirmed':True}
