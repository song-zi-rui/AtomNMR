"""Versioned SQLite online snapshots and an authenticated WebDAV upload bridge.

OpenList can expose a logged-in Quark/Baidu drive through WebDAV. No drive
cookies or executable shell commands are accepted by this service.
"""
import base64
import hashlib
import http.client
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import threading
import time
import uuid
from urllib.parse import urlsplit, quote, unquote
from contextlib import closing
from . import local_secrets

DATABASES = ('assignments.sqlite', 'ai_annotations.sqlite', 'ai_benchmarks.sqlite')
CHUNK = 1024 ** 3


def stamp():
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    tmp = path.with_suffix('.tmp')
    with tmp.open('w', encoding='utf-8') as f:
        json.dump(value, f, ensure_ascii=False, indent=2)
        f.flush(); os.fsync(f.fileno())
    os.chmod(tmp, 0o600)
    tmp.replace(path)


def validate_url(value):
    u = urlsplit(value)
    if (u.scheme not in ('https', 'http') or not u.hostname or u.username or u.password
            or u.query or u.fragment or any(c in value for c in '\r\n')):
        raise ValueError('填写无密钥、无查询参数的 WebDAV 文件夹地址')
    if u.scheme == 'http' and u.hostname not in ('127.0.0.1', 'localhost', '::1'):
        raise ValueError('远程 WebDAV 必须使用 HTTPS；同机 OpenList 可使用 localhost HTTP')
    if '..' in unquote(u.path).split('/'):
        raise ValueError('WebDAV 路径不可包含上级目录')
    return value.rstrip('/')


class WebDAV:
    def __init__(self, url, username, password):
        self.url = validate_url(url)
        self.u = urlsplit(self.url)
        self.auth = 'Basic ' + base64.b64encode((username+':'+password).encode()).decode()

    def request(self, method, relative='', file=None, offset=0, length=None, body=None):
        cls = http.client.HTTPSConnection if self.u.scheme == 'https' else http.client.HTTPConnection
        conn = cls(self.u.hostname, self.u.port, timeout=120)
        path = quote(unquote(self.u.path).rstrip('/') + ('/'+relative if relative else ''), safe='/')
        size = length if file is not None else len(body or b'')
        try:
            conn.putrequest(method, path)
            conn.putheader('Authorization', self.auth)
            conn.putheader('Content-Length', str(size))
            conn.endheaders()
            if file is not None:
                with file.open('rb') as f:
                    f.seek(offset); remaining = size
                    while remaining:
                        block = f.read(min(1024*1024, remaining))
                        if not block: raise ValueError('备份文件在上传过程中发生变化')
                        conn.send(block); remaining -= len(block)
            elif body: conn.send(body)
            response = conn.getresponse()
            result = response.read(2*1024*1024)
            headers = dict((k.lower(), v) for k,v in response.getheaders())
            # Never follow redirects with credentials, nor expose response bodies.
            if response.status not in (200,201,204,207) and not (method=='MKCOL' and response.status==405):
                raise ValueError('WebDAV '+method+' 返回 HTTP '+str(response.status)+'；请检查挂载、写入权限和登录状态')
            return headers, result
        finally: conn.close()

    def mkdir(self, relative): self.request('MKCOL', relative)


def snapshot(source, target, progress=lambda _: None):
    """Each output is transactionally consistent, including committed WAL writes."""
    started=time.monotonic()
    with closing(sqlite3.connect(source.resolve().as_uri()+'?mode=ro', uri=True, timeout=5)) as src, closing(sqlite3.connect(target)) as dst:
        def update(status, remaining, total):
            if time.monotonic()-started > 3600: raise ValueError('单库快照超过一小时；请在低写入时段重试')
            progress({'pages_done':total-remaining,'pages_total':total})
        src.backup(dst, pages=512, progress=update, sleep=.05)
        if dst.execute('PRAGMA quick_check').fetchone()[0] != 'ok': raise ValueError('备份完整性检查未通过')


def digest_parts(path):
    parts=[]; total=hashlib.sha256()
    with path.open('rb') as f:
        index=0
        while True:
            h=hashlib.sha256(); size=0
            while size<CHUNK:
                block=f.read(min(1024*1024, CHUNK-size))
                if not block: break
                total.update(block); h.update(block); size+=len(block)
            if not size: break
            parts.append({'file':path.name+'.part'+str(index).zfill(4), 'offset':index*CHUNK,
                          'bytes':size,'sha256':h.hexdigest()}); index+=1
    return {'file':path.name,'bytes':path.stat().st_size,'sha256':total.hexdigest(),'parts':parts}


class Service:
    def __init__(self, database, autostart=True):
        self.database=Path(database).resolve()
        self.root=self.database.parent/'backup_control'; self.root.mkdir(exist_ok=True)
        os.chmod(self.root,0o700)
        self.config_file=self.root/'settings.json'; self.state_file=self.root/'history.json'
        self.lock=threading.RLock(); self.stop=threading.Event(); self.worker=None
        self.config=json.loads(self.config_file.read_text('utf-8')) if self.config_file.exists() else {
            'enabled':False,'interval_hours':24,'webdav_url':'','username':'','secret':None,
            'github_repo':'https://github.com/song-zi-rui/AtomNMR.git','github_enabled':False,'cleanup_uploaded':True}
        self.history=json.loads(self.state_file.read_text('utf-8')) if self.state_file.exists() else []
        for job in self.history:
            if job['state']=='running': job.update(state='interrupted',error='上次服务退出；可重新上传已完成快照')
        self.scheduler=threading.Thread(target=self._schedule,daemon=True,name='backup-schedule')
        if autostart:self.scheduler.start()

    def status(self):
        with self.lock:
            cfg={k:v for k,v in self.config.items() if k!='secret'}
            cfg['password_saved']=bool(self.config.get('secret'))
            jobs=[]
            for job in self.history[-30:]:
                folder=self.root/job['snapshot_id']
                jobs.append(dict(job,can_retry=job['state'] in ('error','interrupted') and not job['cloud_complete']
                    and (folder/'manifest.json').is_file() and (folder/'assignments.sqlite').is_file()))
            return {'settings':cfg,'jobs':list(reversed(jobs)),
                    'active':bool(self.worker and self.worker.is_alive()),
                    'note':'数据库使用在线快照；各库分别一致，三库不是同一瞬间的跨库事务。云端保留每次版本，不镜像删除。定时任务需服务运行。'}

    def configure(self, payload):
        if not isinstance(payload,dict): raise ValueError('备份设置应为对象')
        allowed={'enabled','interval_hours','webdav_url','username','password','github_repo','github_enabled','cleanup_uploaded'}
        if set(payload)-allowed: raise ValueError('备份设置字段无效')
        with self.lock:
            cfg=self.config|{k:v for k,v in payload.items() if k!='password'}
            if type(cfg['interval_hours']) is not int or not 1<=cfg['interval_hours']<=720: raise ValueError('备份间隔应为1–720小时')
            for k in ('enabled','github_enabled','cleanup_uploaded'):
                if type(cfg[k]) is not bool: raise ValueError('备份开关必须为布尔值')
            if cfg['webdav_url']: cfg['webdav_url']=validate_url(cfg['webdav_url'])
            if not isinstance(cfg['username'],str) or any(c in cfg['username'] for c in '\r\n:'): raise ValueError('WebDAV 用户名无效')
            if cfg['github_repo'] and not re.fullmatch(r'(?:git@github\.com:|https://github\.com/)[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+(?:\.git)?',cfg['github_repo']):
                raise ValueError('GitHub 地址应为 git@github.com:用户/仓库.git 或无密钥的 HTTPS 地址')
            password=payload.get('password')
            if password:
                if not isinstance(password,str) or len(password)>4096: raise ValueError('WebDAV 密码无效')
                cfg['secret']={'kind':'dpapi','value':local_secrets.protect(password)} if os.name=='nt' else {'kind':'private_file','value':password}
            if cfg['enabled'] and not(cfg['webdav_url'] and cfg.get('secret')): raise ValueError('先配置并保存 WebDAV 地址与密码，再启用自动备份')
            if cfg['github_enabled'] and not cfg['github_repo']: raise ValueError('先填写GitHub代码备份仓库')
            write_json(self.config_file,cfg); self.config=cfg
        return self.status()

    def _client(self,cfg):
        secret=cfg.get('secret')
        if not cfg.get('webdav_url') or not secret: raise ValueError('请先配置夸克的 WebDAV 地址和登录密码')
        password=local_secrets.unprotect(secret['value']) if secret['kind']=='dpapi' else secret['value']
        return WebDAV(cfg['webdav_url'],cfg['username'],password)

    def test(self):
        with self.lock: cfg=dict(self.config)
        client=self._client(cfg)
        # A small write + read probe, not just a successful login page.
        name='atomnmr-test-'+uuid.uuid4().hex+'.txt'; value=b'AtomNMR backup write/read probe'
        client.request('PUT',name,body=value)
        _,actual=client.request('GET',name)
        if actual!=value: raise ValueError('WebDAV读回内容不一致')
        client.request('DELETE',name)
        return {'ok':True,'message':'WebDAV写入、读回及删除测试文件成功；尚未备份数据库'}

    def start(self,payload=None):
        payload=payload or {}
        if not isinstance(payload,dict) or set(payload)-{'retry_id'}:raise ValueError('备份操作参数无效')
        with self.lock:
            if self.worker and self.worker.is_alive(): raise ValueError('已有备份任务在运行')
            cfg=dict(self.config)
            self._client(cfg) # Fail before creating a large snapshot.
            retry=payload.get('retry_id'); previous=next((j for j in self.history if j['id']==retry),None)
            if retry and not previous: raise ValueError('备份记录不存在')
            if retry:
                identifier=previous['snapshot_id']
                if previous['cloud_complete'] or not (self.root/identifier/'assignments.sqlite').is_file() or not (self.root/identifier/'manifest.json').is_file(): raise ValueError('该记录无可重传的本地快照，请创建新备份')
            else: identifier='snapshot-'+time.strftime('%Y%m%dT%H%M%SZ',time.gmtime())+'-'+uuid.uuid4().hex[:8]
            job={'id':uuid.uuid4().hex,'snapshot_id':identifier,'state':'running','phase':'queued','started_at':stamp(),'cloud_complete':False,'github_complete':False}
            self.history.append(job); write_json(self.state_file,self.history)
            self.worker=threading.Thread(target=self._run,args=(job,cfg,bool(retry)),daemon=True,name='database-backup');self.worker.start()
            return {'job_id':job['id']}

    def _update(self,job,**value):
        with self.lock:
            job.update(value); write_json(self.state_file,self.history)

    def start_code(self):
        with self.lock:
            if self.worker and self.worker.is_alive():raise ValueError('已有备份任务在运行')
            cfg=dict(self.config)
            if not cfg.get('github_repo'):raise ValueError('请先填写GitHub仓库地址')
            job={'id':uuid.uuid4().hex,'snapshot_id':'code-'+uuid.uuid4().hex[:8], 'kind':'code',
                 'state':'running','phase':'github','started_at':stamp(),'cloud_complete':False,'github_complete':False}
            self.history.append(job);write_json(self.state_file,self.history)
            def run():
                try:
                    self._github(cfg);self._update(job,state='done',phase='done',github_complete=True,completed_at=stamp())
                except Exception:
                    self._update(job,state='error',error='GitHub推送未完成，请配置Git授权并检查仓库写权限；数据库未受影响。',completed_at=stamp())
            self.worker=threading.Thread(target=run,daemon=True,name='code-backup');self.worker.start()
            return {'job_id':job['id']}

    def _github(self,cfg):
        """Push only the server packager's explicit code allowlist, never project data."""
        from deploy.server.package_code import collect_files
        root=Path(__file__).resolve().parents[1]; mirror=self.root/'github-code'
        mirror.mkdir(exist_ok=True)
        env=os.environ|{'GIT_TERMINAL_PROMPT':'0','GCM_INTERACTIVE':'Never','GIT_SSH_COMMAND':'ssh -o BatchMode=yes -o StrictHostKeyChecking=yes'}
        def git(*args,check=True):
            result=subprocess.run(['git',*args],cwd=mirror,env=env,stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180)
            if check and result.returncode: raise ValueError('GitHub代码备份失败；请在服务器配置SSH授权、可信主机及Git用户名邮箱。详细输出未写入网页以免泄露凭据。')
            return result
        git('init','-b','atomnmr-backup')
        remote=git('remote','get-url','origin',check=False)
        if remote.returncode==0:
            if remote.stdout.decode().strip()!=cfg['github_repo']: raise ValueError('已有代码备份仓库地址不同；请使用原地址，或另建备份目录后配置')
        else:git('remote','add','origin',cfg['github_repo'])
        files=collect_files(root)
        new=[]
        for source in files:
            if re.search(rb'\bsk-[A-Za-z0-9_-]{20,}',source.read_bytes()):
                raise ValueError('代码备份发现疑似API密钥，请先从源码中移除再上传')
            relative=source.relative_to(root); dest=mirror/relative
            if not dest.resolve().is_relative_to(mirror.resolve()) or dest.is_symlink():raise ValueError('代码备份路径不安全')
            dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,dest);new.append(relative.as_posix())
        # Only remove previously managed code paths, never user files or .git.
        listing=mirror/'code-files.json'; old=json.loads(listing.read_text('utf-8')) if listing.exists() else []
        for relative in set(old)-set(new):
            dest=mirror/relative
            if dest.resolve().is_relative_to(mirror.resolve()) and '.git' not in dest.parts and dest.is_file():dest.unlink()
        write_json(listing,sorted(new));git('add','--',*new,'code-files.json')
        removed=sorted(set(old)-set(new))
        if removed:git('rm','--cached','--ignore-unmatch','--',*removed)
        staged=git('diff','--cached','--name-only','-z').stdout.decode().split('\x00')
        if any(p and p not in set(new)|set(old)|{'code-files.json'} for p in staged):
            raise ValueError('代码备份区含非程序文件，拒绝提交')
        changed=git('diff','--cached','--quiet',check=False)
        if changed.returncode==1:git('-c','user.name=AtomNMR Backup','-c','user.email=atomnmr-backup@localhost','commit','-m','Update AtomNMR annotation service code')
        elif changed.returncode!=0:raise ValueError('无法检查代码变更')
        git('push','origin','HEAD:refs/heads/atomnmr-backup')

    def _run(self,job,cfg,retry):
        folder=self.root/job['snapshot_id']
        try:
            client=self._client(cfg)
            if not retry:
                folder.mkdir();os.chmod(folder,0o700)
                sources=[self.database if name==DATABASES[0] else self.database.parent/name for name in DATABASES]
                sources=[p for p in sources if p.is_file()]
                size=sum(p.stat().st_size+Path(str(p)+'-wal').stat().st_size if Path(str(p)+'-wal').exists() else p.stat().st_size for p in sources)
                if shutil.disk_usage(folder).free<size*1.15+512*1024**2:raise ValueError('临时快照磁盘空间不足；需要至少数据库总大小的1.15倍加512MB可用空间')
                manifest={'version':1,'snapshot_id':job['snapshot_id'],'started_at':stamp(),
                    'consistency':'Each SQLite file is individually consistent; databases are backed up sequentially.',
                    'restore':'Concatenate each file.parts in order, then verify bytes and SHA256. Stop the application before restore.', 'databases':[]}
                for source in sources:
                    target=folder/('assignments.sqlite' if source==self.database else source.name)
                    self._update(job,phase='snapshot',current_file=target.name)
                    last=[0.]
                    def progress(info):
                        if time.monotonic()-last[0]>2:
                            last[0]=time.monotonic();self._update(job,**info)
                    snapshot(source,target,progress);os.chmod(target,0o600)
                    self._update(job,phase='checksum')
                    manifest['databases'].append(digest_parts(target))
                manifest['completed_at']=stamp();write_json(folder/'manifest.json',manifest)
            else:manifest=json.loads((folder/'manifest.json').read_text('utf-8'))
            client.mkdir(job['snapshot_id']);self._update(job,phase='upload',uploaded_bytes=0)
            sent=0
            for item in manifest['databases']:
                path=folder/item['file']
                if not path.is_file() or digest_parts(path)!=item: raise ValueError('本地快照文件缺失或校验不一致；请创建新备份')
                for part in item['parts']:
                    self._update(job,current_file=part['file'])
                    remote=job['snapshot_id']+'/'+part['file']
                    client.request('PUT',remote,file=path,offset=part['offset'],length=part['bytes'])
                    headers,_=client.request('HEAD',remote)
                    if int(headers.get('content-length','-1'))!=part['bytes']: raise ValueError('云端分块大小校验失败')
                    sent+=part['bytes'];self._update(job,uploaded_bytes=sent)
            encoded=json.dumps(manifest,ensure_ascii=False,indent=2).encode()
            remote=job['snapshot_id']+'/manifest.json'
            client.request('PUT',remote,body=encoded);_,actual=client.request('GET',remote)
            if actual!=encoded:raise ValueError('云端备份清单读回校验失败')
            self._update(job,cloud_complete=True,cloud_completed_at=stamp())
            if cfg['cleanup_uploaded']:
                # Delete only our verified-upload snapshot outputs, not live DBs.
                for item in manifest['databases']:
                    path=folder/item['file']
                    if path.parent==folder and path.name in DATABASES:path.unlink()
            if cfg['github_enabled']:
                self._update(job,phase='github');self._github(cfg);self._update(job,github_complete=True)
            self._update(job,state='done',phase='done',completed_at=stamp())
        except Exception as exc:
            if folder.is_dir() and not (folder/'manifest.json').is_file():
                # Incomplete generated snapshots cannot be resumed. Free their
                # space so a failed attempt cannot fill the live DB disk.
                for name in DATABASES:
                    file=folder/name
                    if file.is_file():file.unlink()
            # Credential-bearing network exceptions/response bodies are never logged.
            message=str(exc) if isinstance(exc,ValueError) else type(exc).__name__+'：备份未完成，请检查磁盘、网盘连接及授权'
            secret=cfg.get('secret',{}).get('value','')
            if secret:message=message.replace(secret,'[REDACTED]')
            self._update(job,state='error',error=message[:800],completed_at=stamp())

    def _schedule(self):
        while not self.stop.wait(30):
            with self.lock:
                cfg=dict(self.config)
                last=next((j['started_at'] for j in reversed(self.history) if j.get('kind')!='code'),None)
            if not cfg['enabled']:continue
            from datetime import datetime,timezone
            elapsed=(datetime.now(timezone.utc)-datetime.fromisoformat(last)).total_seconds() if last else float('inf')
            if elapsed>=cfg['interval_hours']*3600:
                try:self.start()
                except Exception:pass # Existing active work/config error; next tick checks again.

    def shutdown(self):
        self.stop.set()
        if self.scheduler.is_alive():self.scheduler.join(2)
