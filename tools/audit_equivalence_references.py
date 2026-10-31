"""Read-only legacy-run comparison and graph-symmetry screening (not automatic repair)."""
import csv
import json
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/assignment_database_v1/audits/equivalence_20261029'
RUN='533ae96877bc4d1da2631bf63eec25d1'

def connection(name):
    c=sqlite3.connect((ROOT/'data/assignment_database_v1/databases'/name).as_uri()+'?mode=ro',uri=True)
    c.row_factory=sqlite3.Row
    return c

def dump_csv(path,rows):
    if not rows:return
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader()
        for row in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in row.items()})

def graph_groups(mol):
    ranks=Chem.CanonicalRankAtoms(mol,breakTies=False,includeChirality=True)
    values=defaultdict(list)
    for a in mol.GetAtoms():values[ranks[a.GetIdx()]].append(a.GetIdx())
    return [g for g in values.values() if len(g)>1]

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    c=connection('ai_annotations.sqlite');run=dict(c.execute('SELECT * FROM runs WHERE id=?',(RUN,)).fetchone())
    snapshots=[];rows=[];summary=[]
    for item in c.execute('SELECT * FROM items WHERE run_id=? ORDER BY ordinal',(RUN,)):
        item=dict(item)
        context=json.loads(item['input_json']);refs=json.loads(item['reference_json']);output=json.loads(item['output_json'] or 'null')
        mol=Chem.MolFromSmiles(context['smiles']);groups=graph_groups(mol)
        predictions={(p['spectrum_id'],p['peak_id']):p for p in (output or {}).get('assignments',[])}
        shifts={(s['spectrum_id'],p['peak_id']):p.get('shift_ppm') for s in context['spectra'] for p in s['peaks']}
        refs_by_spec=defaultdict(set)
        for ref in refs:
            if ref['kind']=='peak':refs_by_spec[ref['spectrum_id']].update(ref['atom_indices'])
        counts=Counter()
        for ref in refs:
            if ref['kind']!='peak':continue
            p=predictions.get((ref['spectrum_id'],ref['peak_id']),{})
            truth=set(ref['atom_indices']);answer=set(p.get('atom_indices',[]))
            if not output:kind='no_output'
            elif p.get('status') not in ('assigned','equivalent_group','overlap'):kind='abstained'
            elif truth==answer:kind='exact_match'
            else:
                group=next((set(g) for g in groups if truth and truth<set(g) and answer==set(g)),None)
                kind='graph_group_expansion' if group else 'other_difference'
            counts[kind]+=1
            group=next((g for g in groups if truth and truth<set(g)),[])
            rows.append({'item_id':item['id'],'molecule_id':item['molecule_id'],'formula':item['formula'],
                'source':ref['source'],'nucleus':ref['nucleus'],'spectrum_id':ref['spectrum_id'],'peak_id':ref['peak_id'],
                'shift_ppm':shifts.get((ref['spectrum_id'],ref['peak_id'])),'reference_atoms':sorted(truth),
                'ai_status':p.get('status'),'ai_atoms':sorted(answer),'graph_group':group,
                'group_atoms_absent_from_all_snapshot_references':sorted(set(group)-refs_by_spec[ref['spectrum_id']]),'classification':kind})
        summary.append({'item_id':item['id'],'molecule_id':item['molecule_id'],'formula':item['formula'],'state':item['state'],
                        'sources':sorted({r['source'] for r in refs}),'counts':dict(counts)})
        snapshots.append({'item_id':item['id'],'input':context,'references':refs,'output':output,'error':item['error']})
    (OUT/'legacy_run_snapshots.json').write_text(json.dumps({'run_id':RUN,'created_at':run['created_at'],'items':snapshots},ensure_ascii=False,indent=2),encoding='utf-8')
    dump_csv(OUT/'legacy_peak_comparison.csv',rows);dump_csv(OUT/'legacy_molecules.csv',summary)
    main_db=connection('assignments.sqlite');stats=Counter();candidates=[]
    query='''SELECT m.id,m.formula,m.smiles,m.structure_json,s.id spectrum_id,s.nucleus,s.peaks_json,
      json_group_array(json_object('peak_id',a.peak_id,'status',a.status,'atoms',json(a.atoms_json),'origin',a.origin)) labels
      FROM annotations a CROSS JOIN spectra s ON s.id=a.spectrum_id CROSS JOIN molecules m ON m.id=s.molecule_id
      WHERE a.origin='nmrshiftdb2' AND m.origin='nmrshiftdb2' GROUP BY a.spectrum_id'''
    for r in main_db.execute(query):
        labels=json.loads(r['labels']);source=[a for a in labels if a['origin']=='nmrshiftdb2']
        if not source:continue
        stats['source_spectra_examined']+=1
        if stats['source_spectra_examined']%10000==0:print('Examined source spectra:',stats['source_spectra_examined'],flush=True)
        structure=json.loads(r['structure_json']);groups=structure.get('equivalent_groups',[])
        elements={a['index']:a['symbol'] for a in structure['atoms']}
        all_atoms={i for a in labels for i in a['atoms']}
        if groups:
            all_atoms={i for a in main_db.execute('SELECT atoms_json FROM annotations WHERE spectrum_id=?',(r['spectrum_id'],)) for i in json.loads(a[0])}
        peak_shifts={p['peak_id']:p.get('shift_ppm') for p in json.loads(r['peaks_json'])}
        for a in source:
            atoms=set(a['atoms'])
            for group in groups:
                g=set(group)
                if not atoms or not atoms<g:continue
                if r['nucleus']=='13C' and any(elements.get(i)!='C' for i in g):continue
                missing=sorted(g-all_atoms)
                flag='group_atoms_unassigned' if missing else 'group_distributed_across_labels'
                stats[r['nucleus']+':'+flag]+=1
                candidates.append({'molecule_id':r['id'],'formula':r['formula'],'spectrum_id':r['spectrum_id'],'nucleus':r['nucleus'],
                    'peak_id':a['peak_id'],'shift_ppm':peak_shifts.get(a['peak_id']),'source_atoms':sorted(atoms),'graph_group':group,
                    'unassigned_group_atoms':missing,'flag':flag,
                    'priority':'C_missing_member' if r['nucleus']=='13C' and missing else 'review_only',
                    'url':'http://127.0.0.1:8767/?molecule_id='+r['id']})
    dump_csv(OUT/'source_equivalence_candidates.csv',candidates)
    result={'legacy_run':RUN,'legacy_molecules':summary,'legacy_peak_classes':dict(Counter(r['classification'] for r in rows)),
            'source_screen':dict(stats),'candidate_rows':len(candidates),'candidate_molecules':len({r['molecule_id'] for r in candidates}),
            'c_missing_member_molecules':len({r['molecule_id'] for r in candidates if r['priority']=='C_missing_member'}),
            'warning':'Graph symmetry is a screening signal, not proof of NMR equivalence. No labels or scores were changed.'}
    (OUT/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False),flush=True)

if __name__=='__main__':main()
