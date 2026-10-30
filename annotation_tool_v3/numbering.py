"""Explicit, reversible numbering review. No automatic relabelling of source."""
from __future__ import annotations
import hashlib
import json

from rdkit import Chem
from . import core as db
from .curation import state_for
from tools.filter_nmrshiftdb2_experimental import parse_record,classify,SPECTRUM
from tools.shiftdb_atom_identity import canonical_identity,parse_fragments,validate_fragments
from tools.shiftdb_numbering import assess_numbering
from tools.shiftdb_source_index import DEFAULT_INDEX,load_record


def preview(c,mid,sid,*,source_index=DEFAULT_INDEX):
    row=c.execute('SELECT s.*,m.smiles,m.revision FROM spectra s JOIN molecules m ON m.id=s.molecule_id WHERE s.id=? AND s.molecule_id=?',(sid,mid)).fetchone()
    if row is None:raise ValueError('谱不属于当前分子。')
    if row['source']!='nmrshiftdb2':raise ValueError('此编号流程仅用于 ShiftDB2 源谱。')
    provenance=json.loads(row['source_json']);stored_props=provenance.get('original_properties',{})
    field=provenance.get('spectrum_field');value=stored_props.get(field)
    if not isinstance(value,str):raise ValueError('缺少原始谱字段，不能恢复编号。')
    origin={'raw_signal_sha256':hashlib.sha256(value.encode()).hexdigest()}
    if provenance.get('original_molblock'):
        block=provenance['original_molblock'];props=stored_props
        origin['molblock_sha256']=hashlib.sha256(block.encode()).hexdigest()
    else:
        raw,source_meta=load_record(source_index,source_id=provenance.get('nmrshiftdb2_id'))
        raw_block,props,_=parse_record(raw);block=raw_block.decode();origin.update(source_meta)
        if props.get(field)!=value:raise ValueError('原始 SDF 谱字段与当前存档不同，不能采用候选。')
    field_match=SPECTRUM.fullmatch(field)
    if not field_match or field_match[1]!=row['nucleus'] or classify(props,field_match[2],row['nucleus'])[0]!='experimental':
        raise ValueError('来源不是明确实验谱，或核种与谱字段不一致，不能采用编号结果。')
    identity=canonical_identity(Chem.MolFromMolBlock(block,sanitize=False,removeHs=False),
        f"shiftdb:{origin.get('record_index',provenance.get('source_record_index',row['source_key']))}",sanitize=True)
    if identity.smiles!=row['smiles']:raise ValueError('源结构与当前显示结构不同，必须先解决结构映射。')
    conflict=any(validate_fragments(parse_fragments(v,'13C',identity),'13C')['status']=='source_assignment_inconsistency'
        for k,v in props.items() if (m:=SPECTRUM.fullmatch(k)) and m[1]=='13C' and classify(props,m[2],m[1])[0]=='experimental')
    flow=assess_numbering(value,row['nucleus'],identity,source_inchi=props.get('INChI'),
                          record_conflict=conflict and row['nucleus']=='1H')
    stored_peaks=json.loads(row['peaks_json'])
    # Bind source buckets to persisted IDs by exact original fragments, never ppm.
    by_raw={p.get('raw_signal'):p['peak_id'] for p in stored_peaks}
    if len(by_raw)!=len(stored_peaks):raise ValueError('存档峰的原始片段不唯一，需要先完成源峰迁移。')
    for option in [flow['default'],*flow['candidates']]:
        for p in option['peaks']:
            if p['raw_signal'] not in by_raw:raise ValueError('源峰片段与存档不一致，不能按位移猜配峰。')
            p['peak_id']=by_raw[p['raw_signal']]
        if len(option['peaks'])!=len(stored_peaks):raise ValueError('候选峰数与源存档不同。')
    history=[dict(r) for r in c.execute('SELECT id,candidate_id,decision,note,created_at,target_review_id FROM spectrum_numbering_reviews WHERE spectrum_id=? ORDER BY id DESC',(sid,))]
    flow.update(molecule_id=mid,spectrum_id=sid,nucleus=row['nucleus'],revision=row['revision'],
                source=origin,history=history,atom_identity=identity.to_dict(),display_smiles=row['smiles'],
                stored_source_sha256=hashlib.sha256(row['source_json'].encode()).hexdigest(),
                stored_peaks_sha256=hashlib.sha256(row['peaks_json'].encode()).hexdigest())
    return flow


def review(c,mid,payload,*,source_index=DEFAULT_INDEX):
    # InChI generation/source retrieval occur outside the write lock.
    sid=payload.get('spectrum_id');decision=payload.get('decision')
    if decision not in ('confirm','reject','undo'):raise ValueError('编号审核操作无效。')
    note=str(payload.get('note','')).strip()
    if len(note)<8:raise ValueError('请填写至少 8 个字的原文依据或拒绝/撤销原因。')
    flow=preview(c,mid,sid,source_index=source_index)
    if payload.get('revision')!=flow['revision']:raise RuntimeError('记录已更新，请重新打开编号检查。')
    option=None
    if decision!='undo':
        option=next((o for o in [flow['default'],*flow['candidates']] if o['candidate_id']==payload.get('candidate_id')),None)
        if option is None:raise RuntimeError('候选已变化或不存在，请重新检查。')
        if decision=='confirm' and (not option['adoptable'] or payload.get('confirmed_source') is not True or payload.get('confirmed_mapping') is not True):
            raise ValueError('候选未通过检查，或尚未明确确认原始依据及逐峰对应。')
    try:
        c.execute('BEGIN IMMEDIATE')
        revision=c.execute('SELECT revision FROM molecules WHERE id=?',(mid,)).fetchone()
        if revision is None or revision[0]!=flow['revision']:raise RuntimeError('记录已更新，请重新检查。')
        current=c.execute('SELECT source_json,peaks_json FROM spectra WHERE id=?',(sid,)).fetchone()
        if (hashlib.sha256(current['source_json'].encode()).hexdigest()!=flow['stored_source_sha256'] or
                hashlib.sha256(current['peaks_json'].encode()).hexdigest()!=flow['stored_peaks_sha256']):
            raise RuntimeError('源谱已经变化，请重新检查。')
        active=c.execute('SELECT 1 FROM spectra WHERE id=? AND molecule_id=? AND'+db.ACTIVE,(sid,mid)).fetchone()
        if state_for(c,'molecule',mid)!='active' or not active:raise ValueError('分子或谱已停用，不能采用编号结果。')
        before=[];after=[];preserved=0;target_id=None
        if decision=='confirm':
            for p in option['peaks']:
                pid=p['peak_id'];old=c.execute('SELECT * FROM annotations WHERE spectrum_id=? AND peak_id=?',(sid,pid)).fetchone()
                if old is not None and old['origin']!='nmrshiftdb2':preserved+=1;continue
                atoms=sorted({a['target_atom_idx'] for f in p['source_assignment_fragments'] for a in f['assignments'] if a['transfer_status']=='mapped'})
                if not atoms:raise ValueError('存在无法解析的归属，不能采用。')
                before.append({'peak_id':pid,'annotation':dict(old) if old else None})
                c.execute('INSERT INTO annotations VALUES(?,?,?,?,?,?,?) ON CONFLICT(spectrum_id,peak_id) DO UPDATE SET status=excluded.status,atoms_json=excluded.atoms_json,note=excluded.note,origin=excluded.origin,updated_at=excluded.updated_at',
                    (sid,pid,'assigned' if len(atoms)==1 else 'source_group',db.dumps(atoms),
                     '编号审核确认（非自动推断）；'+note,'numbering_review',db.now()))
                after.append(dict(c.execute('SELECT * FROM annotations WHERE spectrum_id=? AND peak_id=?',(sid,pid)).fetchone()))
        elif decision=='undo':
            target_id=payload.get('review_id')
            previous=c.execute("SELECT * FROM spectrum_numbering_reviews WHERE id=? AND spectrum_id=? AND molecule_id=? AND decision='confirm'",(target_id,sid,mid)).fetchone()
            if previous is None:raise ValueError('没有对应的已采用编号记录。')
            if c.execute("SELECT 1 FROM spectrum_numbering_reviews WHERE target_review_id=? AND decision='undo'",(target_id,)).fetchone():raise ValueError('该次采用已经撤销。')
            old_by_pid={x['peak_id']:x['annotation'] for x in json.loads(previous['before_json'])}
            for expected in json.loads(previous['after_json']):
                pid=expected['peak_id'];current=c.execute('SELECT * FROM annotations WHERE spectrum_id=? AND peak_id=?',(sid,pid)).fetchone()
                if current is None or dict(current)!=expected:preserved+=1;continue
                before.append({'peak_id':pid,'annotation':dict(current)})
                restore=old_by_pid[pid]
                if restore is None:c.execute('DELETE FROM annotations WHERE spectrum_id=? AND peak_id=?',(sid,pid))
                else:c.execute('UPDATE annotations SET status=?,atoms_json=?,note=?,origin=?,updated_at=? WHERE spectrum_id=? AND peak_id=?',
                    (restore['status'],restore['atoms_json'],restore['note'],restore['origin'],restore['updated_at'],sid,pid))
                after.append({'peak_id':pid,'restored':restore})
        c.execute('INSERT INTO spectrum_numbering_reviews(molecule_id,spectrum_id,candidate_id,decision,note,flow_json,before_json,after_json,target_review_id,created_at) VALUES(?,?,?,?,?,?,?,?,?,?)',
            (mid,sid,option['candidate_id'] if option else previous['candidate_id'],decision,note,db.dumps(flow),db.dumps(before),db.dumps(after),target_id,db.now()))
        c.execute('UPDATE molecules SET revision=revision+1,task_done=0 WHERE id=?',(mid,))
        db.refresh_molecule(c,mid);db.touch_activity(c,mid,'numbering_'+decision)
        result=db.record(c,mid,include_archived=True)
        result['numbering_result']={'decision':decision,'changed_peaks':len(after),'preserved_manual_or_later_edits':preserved}
        c.commit();return result
    except Exception:c.rollback();raise
