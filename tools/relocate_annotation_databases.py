"""Offline relocation of three SQLite files with checkpoints and byte verification."""
from __future__ import annotations

import argparse
from contextlib import closing
import hashlib
import json
from pathlib import Path
import sqlite3
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from annotation_tool_v3.paths import DATABASE_NAMES,DATA_ROOT,DATABASE_DIR


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def describe(path,checkpoint=False):
    with closing(sqlite3.connect(path.as_uri()+('?mode=rw' if checkpoint else '?mode=ro'),uri=True,timeout=3)) as c:
        if checkpoint:
            result=c.execute('PRAGMA wal_checkpoint(TRUNCATE)').fetchone()
            if result[0]:raise ValueError('Database is still in use; stop all writers before moving: '+str(path))
        if c.execute('PRAGMA quick_check').fetchone()[0]!='ok':raise ValueError('SQLite quick_check failed: '+str(path))
        names={r[0] for r in c.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        wanted=('molecules','spectra','annotations','work_events','auto_jobs','auto_items','api_credentials','experiments','cases','attempts','profiles')
        counts={name:c.execute('SELECT count(*) FROM "'+name+'"').fetchone()[0] for name in wanted if name in names}
    return counts


def relocate(source,target,*,apply=False,root=ROOT,progress=lambda value:None):
    source,target,root=source.resolve(),target.resolve(),root.resolve()
    if not source.is_relative_to(root) or not target.is_relative_to(root):raise ValueError('Both folders must stay inside the explicitly selected workspace')
    if source==target or source.drive.lower()!=target.drive.lower():raise ValueError('Use a distinct folder on the same volume')
    if any(not (source/name).is_file() for name in DATABASE_NAMES):raise ValueError('All three source databases must exist')
    if any((target/(name+suffix)).exists() for name in DATABASE_NAMES for suffix in ('','-wal','-shm','-journal')):
        raise ValueError('Target database or companion already exists; refusing to overwrite')
    if not apply:return {'state':'planned','source':str(source),'target':str(target),'databases':list(DATABASE_NAMES)}
    before={}
    for name in DATABASE_NAMES:
        progress({'phase':'checking','database':name})
        counts=describe(source/name,checkpoint=True)
        progress({'phase':'hashing_source','database':name})
        before[name]={'counts':counts,'bytes':(source/name).stat().st_size,'sha256':digest(source/name)}
    target.mkdir(parents=True,exist_ok=True)
    moved=[]
    try:
        for name in DATABASE_NAMES:
            for suffix in ('','-wal','-shm','-journal'):
                src=source/(name+suffix);dest=target/(name+suffix)
                if src.exists():
                    if suffix in ('-wal','-journal') and src.stat().st_size:raise ValueError('Uncheckpointed journal remains: '+str(src))
                    if dest.exists():raise ValueError('Target appeared during migration: '+str(dest))
                    src.rename(dest);moved.append((src,dest))
        for name in DATABASE_NAMES:
            progress({'phase':'verifying_target','database':name})
            if digest(target/name)!=before[name]['sha256']:raise ValueError('Byte hash changed: '+name)
            if describe(target/name)!=before[name]['counts']:raise ValueError('Row counts changed: '+name)
    except Exception:
        for src,dest in reversed(moved):
            if dest.exists() and not src.exists():dest.rename(src)
        raise
    return {'state':'complete','source':str(source),'target':str(target),'databases':before,
            'method':'checkpoint_then_same_volume_rename','verified':'SHA256, quick_check and selected table counts'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-dir',type=Path,default=DATA_ROOT);p.add_argument('--target-dir',type=Path,default=DATABASE_DIR)
    p.add_argument('--apply',action='store_true');p.add_argument('--report',type=Path,required=True)
    args=p.parse_args()
    if args.report.exists():raise ValueError('Choose a new report path')
    def progress(value):print(json.dumps(value),flush=True)
    result=relocate(args.source_dir,args.target_dir,apply=args.apply,progress=progress)
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False),flush=True)


if __name__=='__main__':main()
