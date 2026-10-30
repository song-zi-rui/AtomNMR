"""Reassemble a downloaded backup into a new directory; never overwrite a DB."""
import argparse
import hashlib
import json
from pathlib import Path
import sqlite3
from contextlib import closing
from .backups import DATABASES


def restore(source, output):
    source=Path(source).resolve();output=Path(output).resolve()
    if output.exists():raise ValueError('恢复目标必须是尚不存在的新文件夹，避免覆盖当前数据库')
    manifest=json.loads((source/'manifest.json').read_text('utf-8'))
    if manifest.get('version')!=1:raise ValueError('不支持的备份格式')
    items=manifest.get('databases',[])
    if not items or len({x['file'] for x in items})!=len(items):raise ValueError('数据库清单无效')
    if 'assignments.sqlite' not in [x['file'] for x in items]:raise ValueError('备份中缺少主数据库')
    for item in items:
        if item['file'] not in DATABASES:raise ValueError('数据库文件名无效')
        for i,part in enumerate(item['parts']):
            if part['file']!=item['file']+'.part'+str(i).zfill(4):raise ValueError('分块文件名或顺序无效')
            if part['offset']!=sum(p['bytes'] for p in item['parts'][:i]):raise ValueError('分块偏移无效')
    output.mkdir(parents=True)
    try:
        for item in items:
            target=output/item['file'];total=hashlib.sha256();size=0
            with target.open('xb') as dst:
                for part in item['parts']:
                    h=hashlib.sha256();length=0
                    with (source/part['file']).open('rb') as f:
                        for block in iter(lambda:f.read(1024*1024),b''):
                            h.update(block);total.update(block);length+=len(block);dst.write(block)
                    if length!=part['bytes'] or h.hexdigest()!=part['sha256']:raise ValueError('分块校验失败：'+part['file'])
                    size+=length
            if size!=item['bytes'] or total.hexdigest()!=item['sha256']:raise ValueError('完整数据库校验失败')
            with closing(sqlite3.connect(target.resolve().as_uri()+'?mode=ro',uri=True)) as c:
                if c.execute('PRAGMA quick_check').fetchone()[0]!='ok':raise ValueError('SQLite完整性检查失败')
        (output/'RESTORE_OK.txt').write_text('所有分块、完整文件 SHA256 与 SQLite quick_check 通过。停止标注服务后再替换数据。',encoding='utf-8')
        return output
    except Exception:
        # Partial outputs are intentionally retained for diagnosis, without a
        # RESTORE_OK marker. They must not be installed into the running app.
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();print(restore(args.source,args.output))

if __name__=='__main__':main()
