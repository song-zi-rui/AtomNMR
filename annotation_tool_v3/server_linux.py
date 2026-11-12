"""Private localhost server with environment credentials; no Windows DPAPI dependency."""
import argparse
import json
import os
from pathlib import Path
import signal
import sqlite3
import threading
from contextlib import closing
from functools import partial

from .ai_annotation import DEFAULTS, Manager, normalize_base_url, request_api
from .server import App,Handler,ThreadingHTTPServer
from . import api_tls, core as db
from .paths import artifacts_root

def environment_key():
    value=os.environ.get('ATOMNMR_API_KEY','').strip()
    if not value or value=='replace_with_your_group_key':
        raise ValueError('服务器环境中缺少 ATOMNMR_API_KEY')
    if not 8<=len(value)<=512 or any(x.isspace() for x in value):
        raise ValueError('ATOMNMR_API_KEY 格式无效')
    return value


class EnvironmentManager(Manager):
    def __init__(self, path, *, base_url=None):
        self.bound_base_url=base_url
        super().__init__(path)

    def refresh_credential(self):
        self.secret=environment_key();self.credential=None;self.secret_error=None;self.secret_saved=False
    def settings(self):
        value=super().settings()
        return value|{'key_storage':'server_environment','credential_name':'服务器环境密钥','key_saved':False}
    def manage_credentials(self,payload):
        raise ValueError('服务器模式请修改私密环境文件并重启服务；不使用Windows网页密钥库')
    def key(self, config=None):
        if self.bound_base_url and normalize_base_url((config or self.settings())['base_url'])!=self.bound_base_url:
            raise ValueError('该任务的 API 地址与服务器密钥绑定地址不一致；未发送请求')
        return super().key(config)
    def configure(self,payload):
        if payload.get('api_key') or payload.get('clear_key'):
            raise ValueError('服务器模式请在私密环境文件中更换密钥')
        if self.bound_base_url and normalize_base_url(payload.get('base_url',self.settings()['base_url']))!=self.bound_base_url:
            raise ValueError('服务器 API 地址请修改 api.json 并重启')
        return super().configure(payload)

def resume_interrupted(app,job_id):
    """Only a configured interrupted job; never automatically replay failed items."""
    if not job_id:return '未配置自动恢复任务；可从网页开始标注'
    with closing(sqlite3.connect(app.ai.path)) as c:
        row=c.execute('SELECT state FROM auto_jobs WHERE id=? AND deleted_at IS NULL',(job_id,)).fetchone()
    if not row:raise ValueError('指定的恢复任务不存在，请检查 ATOMNMR_RESUME_JOB_ID')
    if row[0]!='interrupted':return '任务状态为 '+row[0]+'；保留原状态，不自动恢复'
    app.auto.resume(job_id,use_current_limits=True)
    return '已恢复指定任务的待处理项；失败及结果未知项不自动重试'


def active_recovery_job(app, configured_job):
    with app.ai.lock:
        if not configured_job or app.ai.active!=configured_job or app.ai.cancel.is_set():return None
        with closing(sqlite3.connect(app.ai.path)) as c:
            row=c.execute("SELECT state FROM auto_jobs WHERE id=? AND deleted_at IS NULL",(configured_job,)).fetchone()
        return configured_job if row and row[0]=='running' else None


def mark_stopped_job(app, job_id):
    # Run after the worker has recorded the terminal state of every issued request.
    if job_id:
        with closing(sqlite3.connect(app.ai.path)) as c:
            c.execute("UPDATE auto_jobs SET state='interrupted',error='服务器停止；已保存进度。',updated_at=? WHERE id=? AND state='paused'",
                      (db.now(),job_id))
            c.commit()


def check_api(config):
    settings=DEFAULTS | config
    base_url=normalize_base_url(config.get('base_url'))
    settings['base_url']=base_url
    settings.update(api_tls.validate_config(settings))
    response=request_api(environment_key(),'GET','/models',base_url=base_url,**api_tls.request_options(settings))
    models={row.get('id') for row in response.get('data',[]) if isinstance(row,dict) and isinstance(row.get('id'),str)}
    return {'models_query':'ok','selected_model':settings['model'],'selected_model_found':settings['model'] in models}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--db',type=Path,required=True);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--port',type=int,default=8767);p.add_argument('--check-api',action='store_true')
    args=p.parse_args()
    cfg=json.loads(args.config.read_text(encoding='utf-8'))
    if not isinstance(cfg,dict) or set(cfg)-{'api','execution','server'} or any(not isinstance(v,dict) for v in cfg.values()):
        p.error('配置应为 api、execution 和可选的 server 对象')
    if not cfg.get('api',{}).get('base_url'):p.error('请在 api.json 中填写 api.base_url')
    server=cfg.get('server',{})
    cooldown=server.get('transient_cooldown_s',0)
    if set(server)-{'transient_cooldown_s'} or type(cooldown) is not int or not 0<=cooldown<=1800:
        p.error('server.transient_cooldown_s 应为 0–1800 的整数')
    if args.check_api:
        result=check_api(cfg.get('api',{}))
        print(json.dumps(result,ensure_ascii=False),flush=True)
        return 0 if result['selected_model_found'] else 1
    if not args.db.is_file():p.error('主数据库不存在；拒绝创建空数据库，请先迁移真实快照')
    # Lock before App initialization, which recovers unfinished request states.
    lock=artifacts_root(args.db)/'annotation-server.lock';lock_file=lock.open('a+b')
    try:
        if os.name=='posix':
            import fcntl
            fcntl.flock(lock_file,fcntl.LOCK_EX|fcntl.LOCK_NB)
        else:
            import msvcrt
            if not lock.stat().st_size:lock_file.write(b'0');lock_file.flush()
            lock_file.seek(0);msvcrt.locking(lock_file.fileno(),msvcrt.LK_NBLCK,1)
    except (OSError,BlockingIOError):p.error('同一数据库已有服务器启动入口持有锁，拒绝重复运行')
    httpd=None;stopped_job=None
    try:
        httpd=ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
        base_url=normalize_base_url(cfg['api']['base_url'])
        app=App(args.db.resolve(),manager_factory=partial(EnvironmentManager,base_url=base_url));httpd.app=app
        app.ai.configure(cfg.get('api',{}));app.auto.execution(cfg.get('execution',{}))
        app.auto.unattended_cooldown=cooldown
        pinned_job=os.environ.get('ATOMNMR_RESUME_JOB_ID','').strip()
        stopping=threading.Event()
        def shutdown():
            nonlocal stopped_job
            try:
                stopped_job=active_recovery_job(app,pinned_job)
            finally:
                app.ai.cancel.set();app.ai.cancel_requests()
                httpd.shutdown()
        def stop(signum,frame):
            if not stopping.is_set():
                stopping.set()
                threading.Thread(target=shutdown,daemon=True).start()
        signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
        try:print(resume_interrupted(app,pinned_job),flush=True)
        except ValueError as exc:print('未自动恢复任务：'+str(exc),flush=True)
        print('AtomNMR private server: http://127.0.0.1:'+str(args.port),flush=True)
        httpd.serve_forever()
    finally:
        if httpd is not None:httpd.server_close()
        if hasattr(httpd,'app'):
            httpd.app.ai.shutdown();httpd.app.pool.shutdown(wait=True)
            mark_stopped_job(httpd.app,stopped_job)
        lock_file.close()

if __name__=='__main__':raise SystemExit(main())
