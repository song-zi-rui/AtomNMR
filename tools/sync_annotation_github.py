"""Publish only application code to the authorized GitHub repository. No DB access."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from deploy.server.package_code import collect_files


def sync(repo='https://github.com/song-zi-rui/AtomNMR.git'):
    if not re.fullmatch(r'(?:https://github\.com/|git@github\.com:)[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+(?:\.git)?',repo):raise ValueError('Invalid GitHub repository URL')
    mirror=ROOT/'.code_publish/AtomNMR';mirror.mkdir(parents=True,exist_ok=True)
    env=os.environ|{'GIT_TERMINAL_PROMPT':'0','GCM_INTERACTIVE':'Never','GIT_SSH_COMMAND':'ssh -o BatchMode=yes -o StrictHostKeyChecking=yes'}
    def git(*args,check=True):
        result=subprocess.run(['git',*args],cwd=mirror,env=env,stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180)
        if check and result.returncode:raise ValueError('Git failed: '+args[0]+'; check GitHub authorization, network and branch history (no forced push).')
        return result
    git('init','-b','atomnmr-backup')
    remote=git('remote','get-url','origin',check=False)
    if remote.returncode==0:
        current=remote.stdout.decode().strip()
        identity=lambda url:re.sub(r'^(?:https://github\.com/|git@github\.com:)','',url).removesuffix('.git').lower()
        if identity(current)!=identity(repo):raise ValueError('Existing mirror points to another repository')
        if current!=repo:git('remote','set-url','origin',repo)
    else:git('remote','add','origin',repo)
    mapping={path.relative_to(ROOT).as_posix():path for path in collect_files(ROOT)}
    template=ROOT/'deploy/server/GITHUB_README.md'
    if template.exists():mapping['README.md']=template
    if (ROOT/'LICENSE').is_file():mapping['LICENSE']=ROOT/'LICENSE'
    for relative,source in mapping.items():
        if not source.resolve().is_relative_to(ROOT) or source.is_symlink():raise ValueError('Unsafe source')
        if re.search(rb'\bsk-[A-Za-z0-9_-]{20,}',source.read_bytes()):raise ValueError('Possible API secret in '+relative)
        dest=mirror/relative
        if not dest.resolve().is_relative_to(mirror.resolve()) or dest.is_symlink():raise ValueError('Unsafe destination')
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,dest)
    listing=mirror/'code-files.json';old=json.loads(listing.read_text('utf-8')) if listing.exists() else []
    removed=sorted(set(old)-set(mapping))
    for relative in removed:
        dest=mirror/relative
        if dest.resolve().is_relative_to(mirror.resolve()) and '.git' not in dest.parts and dest.is_file():dest.unlink()
    listing.write_text(json.dumps(sorted(mapping),ensure_ascii=False,indent=2),encoding='utf-8')
    git('add','--',*sorted(mapping),'code-files.json')
    if removed:git('rm','--cached','--ignore-unmatch','--',*removed)
    staged=git('diff','--cached','--name-only','-z').stdout.decode().split('\x00')
    if any(p and p not in set(mapping)|set(old)|{'code-files.json'} for p in staged):raise ValueError('Unexpected staged files')
    changed=git('diff','--cached','--quiet',check=False)
    if changed.returncode==1:git('-c','user.name=AtomNMR Backup','-c','user.email=atomnmr-backup@localhost','commit','-m','Update AtomNMR annotation application')
    elif changed.returncode!=0:raise ValueError('Unable to inspect changes')
    git('push','-u','origin','HEAD:refs/heads/atomnmr-backup')
    commit=git('rev-parse','HEAD').stdout.decode().strip()
    remote_commit=git('ls-remote','origin','refs/heads/atomnmr-backup').stdout.decode().split()[0]
    if commit!=remote_commit:raise ValueError('Remote verification failed')
    return {'repository':repo,'branch':'atomnmr-backup','commit':commit,'files':len(mapping),'remote_verified':True,'database_uploaded':False}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',default='https://github.com/song-zi-rui/AtomNMR.git');args=p.parse_args()
    print(json.dumps(sync(args.repo),ensure_ascii=True))
