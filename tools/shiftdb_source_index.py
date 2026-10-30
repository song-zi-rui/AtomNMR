"""Small seek index for the immutable SDF; never scan the full file in a save."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sqlite3
import sys
if __package__ in (None,''):
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from tools.filter_nmrshiftdb2_experimental import iter_records,parse_record

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_INDEX=ROOT/'data/shiftdb_atom_identity_audit/source_seek_index.sqlite'
DEFAULT_SOURCE=ROOT/'data/nmrshiftdb2/nmrshiftdb2withsignals.sd'


def build_index(source=DEFAULT_SOURCE,index=DEFAULT_INDEX):
    source,index=Path(source).resolve(),Path(index)
    if index.exists():raise FileExistsError('Source index exists; use a new output path')
    index.parent.mkdir(parents=True,exist_ok=True)
    start=source.stat();digest=hashlib.sha256();offset=0;count=0
    c=sqlite3.connect(index)
    try:
        c.executescript('CREATE TABLE meta(key TEXT PRIMARY KEY,value TEXT);'
            'CREATE TABLE records(record_index INTEGER PRIMARY KEY,source_id TEXT,byte_offset INTEGER,byte_length INTEGER,sha256 TEXT);'
            'CREATE INDEX records_source_id ON records(source_id);')
        for ordinal,_,raw in iter_records(source):
            digest.update(raw)
            try:_,props,_=parse_record(raw);source_id=props.get('nmrshiftdb2 ID')
            except Exception:source_id=None
            c.execute('INSERT INTO records VALUES(?,?,?,?,?)',(ordinal,source_id,offset,len(raw),hashlib.sha256(raw).hexdigest()))
            offset+=len(raw);count+=1
        finish=source.stat()
        if (start.st_size,start.st_mtime_ns)!=(finish.st_size,finish.st_mtime_ns) or offset!=start.st_size:
            raise RuntimeError('SDF changed during indexing or trailing bytes were not indexed')
        meta={'source_path':str(source),'size':start.st_size,'mtime_ns':start.st_mtime_ns,
              'sha256':digest.hexdigest(),'record_count':count,'complete':True}
        c.executemany('INSERT INTO meta VALUES(?,?)',[(k,json.dumps(v)) for k,v in meta.items()]);c.commit()
        return meta
    finally:c.close()


def load_record(index=DEFAULT_INDEX,*,source_id=None,record_index=None):
    index=Path(index)
    if not index.exists():raise ValueError('缺少原始 SDF 快速索引，请先运行 tools/shiftdb_source_index.py。')
    c=sqlite3.connect(f'file:{index.as_posix()}?mode=ro',uri=True)
    try:
        meta={k:json.loads(v) for k,v in c.execute('SELECT key,value FROM meta')}
        if not meta.get('complete'):raise ValueError('原始 SDF 索引未完成。')
        rows=c.execute('SELECT record_index,byte_offset,byte_length,sha256 FROM records WHERE '+
                       ('source_id=?' if source_id is not None else 'record_index=?'),
                       (str(source_id) if source_id is not None else record_index,)).fetchall()
        if len(rows)!=1:raise ValueError('原始 SDF 记录不存在或 source ID 不唯一，不能猜测来源。')
        ordinal,offset,length,digest=rows[0]
        source=Path(meta['source_path']);stat=source.stat()
        if (stat.st_size,stat.st_mtime_ns)!=(meta['size'],meta['mtime_ns']):
            raise ValueError('原始 SDF 已变化，必须重建快速索引。')
        with source.open('rb') as handle:handle.seek(offset);raw=handle.read(length)
        if hashlib.sha256(raw).hexdigest()!=digest:raise ValueError('原始 SDF 记录校验失败。')
        return raw,{'source_sha256':meta['sha256'],'record_sha256':digest,'record_index':ordinal}
    finally:c.close()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sdf',type=Path,default=DEFAULT_SOURCE)
    p.add_argument('--output',type=Path,default=DEFAULT_INDEX)
    a=p.parse_args();print(json.dumps(build_index(a.sdf,a.output),ensure_ascii=False,indent=2))
