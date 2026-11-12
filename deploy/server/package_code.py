"""Build an allowlisted code-only ZIP. Never package databases, keys or environments."""
import hashlib
import argparse
import json
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[2]
DEPLOY_FILES = {'package_code.py','snapshot_databases.py','api.example.json','atomnmr.env.example',
                'atomnmr.service','atomnmr-user.service','start-user.sh','requirements.txt','GITHUB_README.md',
                'sync-code.ps1','atomnmr-code-sync.service','atomnmr-code-sync.timer'}

def collect_files(root):
    files=[]
    for folder in ['annotation_tool_v3','atompeak_data','tools','dict','deploy/server']:
        for path in (root/folder).rglob('*'):
            if not path.is_file() or path.is_symlink() or not path.resolve().is_relative_to(root.resolve()):continue
            rel=path.relative_to(root).as_posix()
            if folder=='deploy/server':allowed=path.parent==root/folder and path.name in DEPLOY_FILES
            elif folder=='dict':allowed=path.suffix=='.txt'
            elif folder=='annotation_tool_v3':
                allowed=(path.suffix in ('.py','.html','.css','.js','.ps1','.cmd') or rel in
                         ('annotation_tool_v3/README.md','annotation_tool_v3/certificates/ai-cluster-ca.crt'))
            else:allowed=path.suffix=='.py'
            if allowed:files.append(path)
    files.append(root/'nmrpeak/utils/vocab_rules.py')
    doc=root/'docs/组内服务器部署与持续标注_20261030.md'
    if doc.exists():files.append(doc)
    user_doc=root/'docs/Ubuntu无sudo服务器部署_20261030.md'
    if user_doc.exists():files.append(user_doc)
    operations_doc=root/'docs/官方谷价收尾与数据整理_20261030.md'
    if operations_doc.exists():files.append(operations_doc)
    layout_doc=root/'docs/标注系统界面与数据库目录_20261031.md'
    if layout_doc.exists():files.append(layout_doc)
    repair_doc=root/'docs/AI等价组重标_暂停收尾_流式修复_20261031.md'
    if repair_doc.exists():files.append(repair_doc)
    current_doc=root/'docs/AtomNMR_组内服务器持续标注_20261112.md'
    if current_doc.exists():files.append(current_doc)
    return files

def main():
    root=Path(__file__).resolve().parents[2]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=root/'deploy/atomnmr-server-code-20261112.zip')
    target=parser.parse_args().output
    if target.exists():raise ValueError('Choose a new archive path; previous code packages are preserved')
    target.parent.mkdir(parents=True,exist_ok=True)
    files=collect_files(root)
    manifest={str(p.relative_to(root)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}
    with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(files):z.write(p,p.relative_to(root))
        z.writestr('server_code_manifest.json',json.dumps(manifest,ensure_ascii=False,indent=2))
    print(json.dumps({'archive':str(target),'files':len(files),'bytes':target.stat().st_size,'contains_database':False},ensure_ascii=False))

if __name__=='__main__':main()
