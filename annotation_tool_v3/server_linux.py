"""Private localhost server with environment credentials; no Windows DPAPI dependency."""
import argparse
import json
import os
from pathlib import Path
import signal
import sqlite3
import threading
from contextlib import closing

from .ai_annotation import Manager
from .server import App,Handler,ThreadingHTTPServer
from . import api_tls

class EnvironmentManager(Manager):
    def refresh_credential(self):
        value=os.environ.get('ATOMNMR_API_KEY','').strip()
        if not value:raise ValueError('服务器环境中缺少 ATOMNMR_API_KEY')
        if not 8<=len(value)<=512 or any(x.isspace() for x in value):raise ValueError('ATOMNMR_API_KEY 格式无效')
        self.secret=value;self.credential=None;self.secret_error=None;self.secret_saved=False
    def settings(self):
        value=super().settings()
        return value|{'key_storage':'server_environment','credential_name':'服务器环境密钥','key_saved':False}
    def manage_credentials(self,payload):
        raise ValueError('服务器模式请修改私密环境文件并重启服务；不使用Windows网页密钥库')

def resume_interrupted(app,job_id):
    """Only a configured interrupted job; never automatically replay failed items."""
    if not job_id:return '未配置自动恢复任务；可从网页开始标注'
    with closing(sqlite3.connect(app.ai.path)) as c:
        row=c.execute('SELECT state FROM auto_jobs WHERE id=?',(job_id,)).fetchone()
    if not row:raise ValueError('指定的恢复任务不存在，请检查 ATOMNMR_RESUME_JOB_ID')
    if row[0]!='interrupted':return '任务状态为 '+row[0]+'；保留原状态，不自动恢复'
    app.auto.resume(job_id,use_current_limits=True)
    return '已恢复指定任务的待处理项；失败及结果未知项不自动重试'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--db',type=Path,required=True);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--port',type=int,default=8767);p.add_argument('--check-api',action='store_true')
    args=p.parse_args()
    if not args.db.is_file():p.error('主数据库不存在；拒绝创建空数据库，请先迁移真实快照')
    cfg=json.loads(args.config.read_text(encoding='utf-8'))
    if not isinstance(cfg,dict) or set(cfg)-{'api','execution'}:p.error('配置应只包含api和execution')
    # Lock before App initialization, which recovers unfinished request states.
    lock=args.db.parent/'annotation-server.lock';lock_file=lock.open('a+b')
    try:
        if os.name=='posix':
            import fcntl
            fcntl.flock(lock_file,fcntl.LOCK_EX|fcntl.LOCK_NB)
        else:
            import msvcrt
            if not lock.stat().st_size:lock_file.write(b'0');lock_file.flush()
            lock_file.seek(0);msvcrt.locking(lock_file.fileno(),msvcrt.LK_NBLCK,1)
    except (OSError,BlockingIOError):p.error('同一数据库已有服务器启动入口持有锁，拒绝重复运行')
    httpd=ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    try:
        app=App(args.db.resolve(),manager_factory=EnvironmentManager);httpd.app=app
        app.ai.configure(cfg.get('api',{}));app.auto.execution(cfg.get('execution',{}))
        if args.check_api:
            from .ai_annotation import request_api
            s=app.ai.settings()
            value=request_api(app.ai.key(),'GET','/models',base_url=s['base_url'],**api_tls.request_options(s))
            models=[r['id'] for r in value.get('data',[]) if isinstance(r.get('id'),str)]
            print(json.dumps({'models_query':'ok','selected_model':s['model'],'selected_model_found':s['model'] in models},ensure_ascii=False),flush=True)
            return
        print(resume_interrupted(app,os.environ.get('ATOMNMR_RESUME_JOB_ID','').strip()),flush=True)
        def stop(signum,frame):
            app.ai.cancel.set();app.ai.cancel_requests()
            threading.Thread(target=httpd.shutdown,daemon=True).start()
        signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
        print('AtomNMR private server: http://127.0.0.1:'+str(args.port),flush=True)
        httpd.serve_forever()
    finally:
        httpd.server_close()
        if hasattr(httpd,'app'):
            httpd.app.backups.shutdown()
            httpd.app.ai.shutdown();httpd.app.pool.shutdown(wait=True)
        lock_file.close()

if __name__=='__main__':main()
