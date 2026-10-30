"""Whole-snapshot assignment audit and reversible, source-label-only repair.

Run audit first. --apply reads that audit; it never invents a lost source map.
Raw SDF and human annotations are preserved. A rollback archive precedes each
mutation and repairs use short per-molecule transactions with revision checks.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import sqlite3
import time
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rdkit import Chem, RDLogger
from tools.filter_nmrshiftdb2_experimental import iter_records, parse_record, classify, SPECTRUM, sha256
from tools.shiftdb_atom_identity import canonical_identity, source_peaks, POLICY
from tools.shiftdb_numbering import assess_numbering

DEFAULT_OUT=Path('data/shiftdb_atom_identity_audit/full_snapshot')


def audit(source, folder):
    folder.mkdir(parents=True,exist_ok=True)
    target=folder/'audit.sqlite'
    if target.exists(): raise FileExistsError(f'Audit already exists: {target}; use another output folder')
    digest=sha256(source);counts=Counter();issue_counts=Counter();kinds=Counter(); started=time.monotonic()
    c=sqlite3.connect(target)
    c.executescript('''
      CREATE TABLE records(record_index INTEGER PRIMARY KEY,source_id TEXT,smiles TEXT,identity_json TEXT,error TEXT);
      CREATE TABLE spectra(record_index INTEGER,field TEXT,nucleus TEXT,kind TEXT,raw_value TEXT,
        validation_json TEXT,peaks_json TEXT,labels_json TEXT,PRIMARY KEY(record_index,field));
      CREATE TABLE meta(key TEXT PRIMARY KEY,value TEXT);
    ''')
    c.execute('INSERT INTO meta VALUES(?,?)',('source_sha256',digest))
    RDLogger.DisableLog('rdApp.warning');RDLogger.DisableLog('rdApp.error')
    with (folder/'issues.jsonl').open('w',encoding='utf-8') as issues_file:
        for index,line,raw in iter_records(source):
            counts['all_records']+=1
            props = {}
            try:
                block,props,_=parse_record(raw)
                # Stable identities assigned before sanitization/removal.
                original=Chem.MolFromMolBlock(block.decode(),sanitize=False,removeHs=False)
                if original is None: raise ValueError('invalid_structure')
                identity=canonical_identity(original,f'shiftdb:{index}',sanitize=True)
                fields=[]
                for field,value in props.items():
                    m=SPECTRUM.fullmatch(field)
                    if not m or m[1] not in ('13C','1H'):continue
                    kind=classify(props,m[2],m[1])[0]
                    peaks,labels,validation=source_peaks(value,m[1],identity)
                    fields.append([field,value,m[1],kind,peaks,labels,validation])
                conflict=any(nuc=='13C' and kind=='experimental' and v['status']=='source_assignment_inconsistency'
                             for _,_,nuc,kind,_,_,v in fields)
                counts['valid_records']+=1
                if conflict: counts['experimental_carbon_conflict_records']+=1
                for field,value,nuc,kind,peaks,labels,v in fields:
                    if nuc=='1H' and kind=='experimental' and conflict:
                        peaks,labels,v=source_peaks(value,nuc,identity,record_conflict=True)
                    if kind=='experimental' and v['status']!='identity_validated':
                        v['numbering_workflow']=assess_numbering(value,nuc,identity,source_inchi=props.get('INChI'),
                            record_conflict=nuc=='1H' and conflict)
                    kinds[f'{kind}:{nuc}:{v["status"]}']+=1
                    counts['all_ch_spectra']+=1
                    if v['status']!='identity_validated':
                        counts['flagged_spectra']+=1;issue_counts.update(v['issue_counts'])
                        issues_file.write(json.dumps({'record_index_1based':index,'source_id':props.get('nmrshiftdb2 ID'),
                            'field':field,'kind':kind,'raw_value':value,'validation':v},ensure_ascii=False)+'\n')
                    c.execute('INSERT INTO spectra VALUES(?,?,?,?,?,?,?,?)',(index,field,nuc,kind,value,
                        json.dumps(v,ensure_ascii=False),json.dumps(peaks,ensure_ascii=False) if kind=='experimental' else None,
                        json.dumps(labels,ensure_ascii=False)))
                c.execute('INSERT INTO records VALUES(?,?,?,?,?)',(index,props.get('nmrshiftdb2 ID'),identity.smiles,
                    json.dumps(identity.to_dict(),ensure_ascii=False),None))
            except Exception as exc:
                counts['invalid_or_mapping_error_records']+=1
                category = 'molfile_parse_failed' if str(exc)=='invalid_structure' else 'canonical_mapping_failed' if 'parser_error:' in str(exc) else 'sanitize_or_other_structure_failed'
                counts[category]+=1
                counts['unassessed_ch_spectra'] += sum(bool((m:=SPECTRUM.fullmatch(k)) and m[1] in ('13C','1H')) for k in props)
                issue_counts['parser_error' if 'parser_error:' in str(exc) else 'invalid_structure_or_record']+=1
                c.execute('INSERT INTO records VALUES(?,?,?,?,?)',(index,props.get('nmrshiftdb2 ID'),None,None,str(exc)))
                for field,value in props.items():
                    m=SPECTRUM.fullmatch(field)
                    if not m or m[1] not in ('13C','1H'):continue
                    validation={'policy':POLICY,'status':'atom_identity_mapping_failed','mapping_error':str(exc),
                                'issues':[],'issue_counts':{},'identity_validated_is_not_chemical_assignment_validation':True}
                    c.execute('INSERT INTO spectra VALUES(?,?,?,?,?,?,?,?)',(index,field,m[1],classify(props,m[2],m[1])[0],
                              value,json.dumps(validation),None,'{}'))
                issues_file.write(json.dumps({'record_index_1based':index,'source_start_line':line,'error':str(exc)},ensure_ascii=False)+'\n')
            if index%2000==0:
                c.commit();print(f'Audited {index:,} source records, {counts["flagged_spectra"]:,} flagged spectra, {time.monotonic()-started:.1f}s',flush=True)
    c.commit()
    end_digest=sha256(source)
    if end_digest!=digest: raise RuntimeError('Source changed during audit')
    report={'policy':POLICY,'source_sdf':str(source.resolve()),'source_sha256':digest,'source_sha256_unchanged':True,
            'completed_at':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':round(time.monotonic()-started,2),
            'counts':dict(counts),'spectrum_status_counts':dict(sorted(kinds.items())),'issue_counts':dict(issue_counts),
            'limitations':['Identity compatibility is not proof of chemical correctness.',
                'No inferred +1/-1, carbon-only, peak-order or chemical-shift repair is applied.',
                'An inconsistency in this snapshot may originate in original annotation or official export; internal DB cannot be inferred.',
                'Carbon S/D/T/Q hydrogen-class disagreement is a warning, not alone proof of wrong assignment.',
                'Proton occurrences remain separate; repeated values do not establish integral or equivalent nuclei.']}
    c.execute('INSERT INTO meta VALUES(?,?)',('completed_report',json.dumps(report,ensure_ascii=False)));c.commit();c.close()
    (folder/'report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    write_report(folder,report)
    print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)
    return report


def write_report(folder,report):
    lines=['# ShiftDB2 原子身份全量检查','',f'策略：`{POLICY}`','',
           '编号语义：withsignals.sd 第三个字段是从 0 开始的原始可见 molfile 原子位置；不是当前显示编号，也不是碳原子序号。',
           '1H 的不可见内部 H 转成相连可见原子；原文件中可见的显式 H 仍保留独立来源身份。','',
           '## 记录统计','', '| 项目 | 数量 |','|---|---:|']
    lines += [f'| {k} | {v:,} |' for k,v in report['counts'].items()]
    lines += ['','## 谱图统计','','| 来源类型 / 核种 / 状态 | 谱数 |','|---|---:|']
    lines += [f'| {k} | {v:,} |' for k,v in report['spectrum_status_counts'].items()]
    lines += ['','## 问题计数','','| 类型 | 发生次数 |','|---|---:|']
    lines += [f'| {k} | {v:,} |' for k,v in report['issue_counts'].items()]
    lines += ['','问题计数为问题发生次数；同一谱可有多个问题，不应相加当谱数。',
              '详细列表：`issues.jsonl`。完整谱级结果、原编号与显示编号映射：`audit.sqlite`。','',
              '## 解释边界','']+[f'- {x}' for x in report['limitations']]
    lines += ['',f'原始 SHA256：`{report["source_sha256"]}`；检查后不变。',
              '预测谱单独统计，不用于补实验标签。源冲突谱保留峰、原文与候选，但不自动输出可信归属。',
              '已有人工标注不被源库重新解析覆盖。已生成的旧 LMDB、环境数据库和标注教学缓存需重新生成或按问题清单排除。']
    (folder/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def apply_audit(path,folder,source_key_prefix=None):
    from annotation_tool_v3 import core as db
    from tools.shiftdb_atom_identity import POLICY
    audit_c=sqlite3.connect(f'file:{(folder/"audit.sqlite").resolve().as_posix()}?mode=ro',uri=True)
    audit_c.row_factory=sqlite3.Row
    if not audit_c.execute("SELECT 1 FROM meta WHERE key='completed_report'").fetchone():raise RuntimeError('Audit is not complete')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')
    backup=folder/f'repair_rollback_{stamp}.sqlite'
    archive=sqlite3.connect(backup)
    # Full durability, with fewer journal fsyncs than DELETE mode. Archives
    # stay readable while repair runs and checkpoint on clean connection close.
    archive.execute('PRAGMA journal_mode=WAL')
    archive.execute('PRAGMA synchronous=FULL')
    archive.executescript('''CREATE TABLE snapshots(spectrum_id TEXT PRIMARY KEY,molecule_id TEXT,spectrum_json TEXT,annotations_json TEXT,molecule_json TEXT);
        CREATE TABLE audit_before(spectrum_id TEXT PRIMARY KEY,row_json TEXT);
        CREATE TABLE applied(spectrum_id TEXT PRIMARY KEY,molecule_id TEXT,revision INTEGER,peaks_json TEXT,annotations_json TEXT,audit_json TEXT);
        CREATE TABLE meta(key TEXT PRIMARY KEY,value TEXT);''')
    archive.execute('INSERT INTO meta VALUES(?,?)',('production_db',str(path.resolve())));archive.commit()
    c=db.connect(path,timeout=2)
    c.execute('''CREATE TABLE IF NOT EXISTS spectrum_assignment_audit(spectrum_id TEXT PRIMARY KEY REFERENCES spectra(id),
      policy TEXT NOT NULL,validation_json TEXT NOT NULL,identity_json TEXT NOT NULL,updated_at TEXT NOT NULL)''');c.commit()
    counts=Counter();details=[]
    counts['audited_source_spectra_seen'] = 0
    print('Selecting source spectra using ' + ('source-key index' if source_key_prefix else 'source scan'),flush=True)
    if source_key_prefix:
        upper=source_key_prefix[:-1]+chr(ord(source_key_prefix[-1])+1)
        ids=c.execute("SELECT id,molecule_id FROM spectra WHERE source_key>=? AND source_key<? AND source='nmrshiftdb2'",
                      (source_key_prefix,upper)).fetchall()
        ids.sort(key=lambda r:(r['molecule_id'],r['id']))
    else:
        ids=c.execute("SELECT id,molecule_id FROM spectra WHERE source='nmrshiftdb2'").fetchall()
        ids.sort(key=lambda r:(r['molecule_id'],r['id']))
    print(f'Selected {len(ids):,} source spectra',flush=True)
    for item in ids:
        counts['audited_source_spectra_seen'] += 1
        old=c.execute('SELECT * FROM spectra WHERE id=?',(item['id'],)).fetchone()
        prov=json.loads(old['source_json']);ri=prov.get('source_record_index');field=prov.get('spectrum_field')
        try:ri=int(ri)
        except (ValueError,TypeError):counts['unresolved_provenance']+=1;continue
        a=audit_c.execute('SELECT s.*,r.identity_json,r.smiles FROM spectra s JOIN records r USING(record_index) WHERE s.record_index=? AND s.field=?',(ri,field)).fetchone()
        if not a or a['kind']!='experimental':counts['unresolved_provenance']+=1;continue
        m=c.execute('SELECT * FROM molecules WHERE id=?',(old['molecule_id'],)).fetchone()
        if prov.get('original_properties',{}).get(field)!=a['raw_value']:
            counts['source_or_structure_changed']+=1;continue
        old_peaks=json.loads(old['peaks_json'])
        validation=json.loads(a['validation_json'])
        validation_json=a['validation_json']
        identity=json.loads(a['identity_json']) if a['identity_json'] else {'policy':POLICY,'status':'unresolved','source_key':f'shiftdb:{ri}'}
        peaks=json.loads(a['peaks_json']) if a['peaks_json'] else old_peaks
        if a['peaks_json'] is None:
            for peak in peaks:peak['source_assignment_validation']='atom_identity_mapping_failed'
        structure_disagrees=a['smiles'] is not None and m['smiles']!=a['smiles']
        if structure_disagrees:
            # Do not choose an arbitrary isomorphism or change a display graph
            # on which human labels already depend. Withhold source labels.
            peaks=old_peaks
            validation={'policy':POLICY,'status':'stored_structure_identity_mismatch','issues':[],
                        'stored_smiles':m['smiles'],'source_canonical_smiles':a['smiles']}
            validation_json=db.dumps(validation)
            identity={'policy':POLICY,'status':'unresolved','source_key':f'shiftdb:{ri}'}
            for peak in peaks:peak['source_assignment_validation']=validation['status']
        # Retain stable peak IDs: match exact source text, never ppm heuristics.
        raw_to_old={p.get('raw_signal'):p['peak_id'] for p in old_peaks}
        labels={} if structure_disagrees else json.loads(a['labels_json'])
        if any(p['raw_signal'] not in raw_to_old for p in peaks):counts['peak_identity_changed']+=1;continue
        remapped_labels={raw_to_old[p['raw_signal']]:labels[p['peak_id']] for p in peaks if p['peak_id'] in labels}
        for p in peaks:p['peak_id']=raw_to_old[p['raw_signal']]
        annotations=[dict(x) for x in c.execute('SELECT * FROM annotations WHERE spectrum_id=?',(old['id'],))]
        previous_audit=c.execute('SELECT * FROM spectrum_assignment_audit WHERE spectrum_id=?',(old['id'],)).fetchone()
        if previous_audit and previous_audit['validation_json']==validation_json and old['peaks_json']==db.dumps(peaks):
            counts['already_updated']+=1;continue
        archive.execute('INSERT OR IGNORE INTO snapshots VALUES(?,?,?,?,?)',(old['id'],old['molecule_id'],db.dumps(dict(old)),db.dumps(annotations),db.dumps(dict(m))))
        archive.execute('INSERT OR IGNORE INTO audit_before VALUES(?,?)',(old['id'],db.dumps(dict(previous_audit)) if previous_audit else None));archive.commit()
        try:
            c.execute('BEGIN IMMEDIATE')
            current=c.execute('SELECT revision FROM molecules WHERE id=?',(m['id'],)).fetchone()[0]
            current_old=c.execute('SELECT peaks_json,source_json FROM spectra WHERE id=?',(old['id'],)).fetchone()
            if current!=m['revision'] or current_old[0]!=old['peaks_json'] or current_old[1]!=old['source_json']:
                c.rollback();counts['concurrent_edit_preserved']+=1;continue
            current_annotations=[dict(x) for x in c.execute('SELECT * FROM annotations WHERE spectrum_id=?',(old['id'],))]
            if current_annotations!=annotations:
                c.rollback();counts['concurrent_edit_preserved']+=1;continue
            for row in annotations:
                if row['origin']!='nmrshiftdb2':counts['human_or_other_labels_preserved']+=1;continue
                target=remapped_labels.get(row['peak_id'])
                if target:
                    if json.loads(row['atoms_json'])!=target['atom_indices'] or row['status']!=target['status']:
                        c.execute('UPDATE annotations SET atoms_json=?,status=?,note=?,updated_at=? WHERE spectrum_id=? AND peak_id=?',
                            (db.dumps(target['atom_indices']),target['status'],target['note'],db.now(),old['id'],row['peak_id']));counts['parser_mapping_labels_corrected']+=1
                elif row['status']!='unassigned' or json.loads(row['atoms_json']):
                    c.execute("UPDATE annotations SET status='unassigned',atoms_json='[]',note=?,updated_at=? WHERE spectrum_id=? AND peak_id=?",
                        ('源归属身份冲突或警告；原始候选在 source_assignment_fragments，旧值在修复归档。',db.now(),old['id'],row['peak_id']));counts['untrusted_source_labels_withdrawn']+=1
            existing_ids={row['peak_id'] for row in annotations}
            for pid,target in remapped_labels.items():
                if pid not in existing_ids:
                    c.execute('INSERT INTO annotations VALUES(?,?,?,?,?,?,?)',(old['id'],pid,target['status'],db.dumps(target['atom_indices']),target['note'],'nmrshiftdb2',db.now()));counts['verified_source_labels_restored']+=1
            c.execute('UPDATE spectra SET peaks_json=? WHERE id=?',(db.dumps(peaks),old['id']))
            c.execute('INSERT INTO spectrum_assignment_audit VALUES(?,?,?,?,?) ON CONFLICT(spectrum_id) DO UPDATE SET policy=excluded.policy,validation_json=excluded.validation_json,identity_json=excluded.identity_json,updated_at=excluded.updated_at',
                      (old['id'],POLICY,validation_json,db.dumps(identity),db.now()))
            db.refresh_molecule(c,m['id'])
            # Invalidate optimistic clients, but retain human completion decisions.
            c.execute('UPDATE molecules SET revision=revision+1 WHERE id=?',(m['id'],));c.commit()
            after_annotations=[dict(x) for x in c.execute('SELECT * FROM annotations WHERE spectrum_id=? ORDER BY peak_id',(old['id'],))]
            after_audit=dict(c.execute('SELECT * FROM spectrum_assignment_audit WHERE spectrum_id=?',(old['id'],)).fetchone())
            archive.execute('INSERT INTO applied VALUES(?,?,?,?,?,?)',(old['id'],m['id'],m['revision']+1,
                            db.dumps(peaks),db.dumps(after_annotations),db.dumps(after_audit)));archive.commit()
            counts['updated_spectra']+=1
            if validation['status']!='identity_validated':details.append({'molecule_id':m['id'],'spectrum_id':old['id'],'field':field,'validation':validation})
        except Exception:
            c.rollback();raise
        if counts['updated_spectra']%1000==0:print('Repair',dict(counts),flush=True)
    archive.close();c.close();audit_c.close()
    result={'policy':POLICY,'source_key_prefix':source_key_prefix,'counts':dict(counts),'rollback_archive':str(backup.resolve()),'flagged_production_spectra':details}
    (folder/f'production_repair_{stamp}.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='flagged_production_spectra'},ensure_ascii=False,indent=2),flush=True)
    return result


def rollback_repair(path, archive_path):
    """Restore a repair delta only when no subsequent molecule edits occurred.

    The database is not replaced, deleted, or rewound. Revisions remain monotone
    and concurrent/new human edits make the entire molecule ineligible to roll
    back. An incomplete applied-state archive is not guessed at.
    """
    from annotation_tool_v3 import core as db
    from contextlib import closing
    counts=Counter()
    with closing(sqlite3.connect(f'file:{archive_path.resolve().as_posix()}?mode=ro',uri=True)) as b, closing(db.connect(path,timeout=2)) as c:
        b.row_factory=sqlite3.Row
        recorded=b.execute("SELECT value FROM meta WHERE key='production_db'").fetchone()
        if not recorded or Path(recorded[0]).resolve()!=path.resolve():raise ValueError('Rollback archive belongs to another database')
        for row in b.execute('SELECT molecule_id,max(revision) revision FROM applied GROUP BY molecule_id'):
            saved=b.execute('SELECT a.*,s.spectrum_json,s.annotations_json before_annotations FROM applied a JOIN snapshots s USING(spectrum_id) WHERE a.molecule_id=?',(row['molecule_id'],)).fetchall()
            c.execute('BEGIN IMMEDIATE')
            try:
                current=c.execute('SELECT revision FROM molecules WHERE id=?',(row['molecule_id'],)).fetchone()
                intact=bool(current and current[0]==row['revision'])
                for s in saved:
                    p=c.execute('SELECT peaks_json FROM spectra WHERE id=?',(s['spectrum_id'],)).fetchone()
                    annotations=[dict(x) for x in c.execute('SELECT * FROM annotations WHERE spectrum_id=? ORDER BY peak_id',(s['spectrum_id'],))]
                    current_audit=c.execute('SELECT * FROM spectrum_assignment_audit WHERE spectrum_id=?',(s['spectrum_id'],)).fetchone()
                    intact &= bool(p and p[0]==s['peaks_json'] and db.dumps(annotations)==s['annotations_json']
                                   and current_audit and db.dumps(dict(current_audit))==s['audit_json'])
                if not intact:
                    c.rollback();counts['later_edits_preserved']+=1;continue
                for s in saved:
                    original=json.loads(s['spectrum_json'])
                    c.execute('UPDATE spectra SET peaks_json=? WHERE id=?',(original['peaks_json'],s['spectrum_id']))
                    c.execute('DELETE FROM annotations WHERE spectrum_id=?',(s['spectrum_id'],))
                    for a in json.loads(s['before_annotations']):
                        c.execute('INSERT INTO annotations VALUES(?,?,?,?,?,?,?)',tuple(a[k] for k in ('spectrum_id','peak_id','status','atoms_json','note','origin','updated_at')))
                    prior=b.execute('SELECT row_json FROM audit_before WHERE spectrum_id=?',(s['spectrum_id'],)).fetchone()[0]
                    c.execute('DELETE FROM spectrum_assignment_audit WHERE spectrum_id=?',(s['spectrum_id'],))
                    if prior:
                        a=json.loads(prior)
                        c.execute('INSERT INTO spectrum_assignment_audit VALUES(?,?,?,?,?)',tuple(a[k] for k in ('spectrum_id','policy','validation_json','identity_json','updated_at')))
                    counts['restored_spectra']+=1
                db.refresh_molecule(c,row['molecule_id'])
                c.execute('UPDATE molecules SET revision=revision+1 WHERE id=?',(row['molecule_id'],));c.commit()
                counts['restored_molecules']+=1
            except Exception:
                c.rollback();raise
    print(json.dumps(dict(counts),ensure_ascii=False),flush=True)
    return dict(counts)


def main():
    p=argparse.ArgumentParser();p.add_argument('--sdf',type=Path,default=Path('data/nmrshiftdb2/nmrshiftdb2withsignals.sd'))
    p.add_argument('--output-dir',type=Path,default=DEFAULT_OUT);p.add_argument('--apply',action='store_true')
    p.add_argument('--db',type=Path,default=Path('data/assignment_database_v1/assignments.sqlite'))
    p.add_argument('--rollback',type=Path,help='Explicitly restore a repair archive; refuses molecules edited since repair')
    p.add_argument('--source-key-prefix',help='Use indexed import source keys to avoid scanning unrelated spectra')
    args=p.parse_args()
    if args.rollback:rollback_repair(args.db,args.rollback)
    elif args.apply:apply_audit(args.db,args.output_dir,args.source_key_prefix)
    else:audit(args.sdf,args.output_dir)

if __name__=='__main__':main()
