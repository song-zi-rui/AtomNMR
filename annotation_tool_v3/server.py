from __future__ import annotations

import argparse
import json
import mimetypes
import os
import sqlite3
import threading
import time
import uuid
import webbrowser
from contextlib import closing
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from datetime import datetime, timezone
from urllib.parse import parse_qs, unquote, urlparse

from . import core as db
from .exporters import export_database, preview_database
from .export_policy import validate as validate_export_options
from . import curation
from . import numbering
from . import workload
from .bmrb_reference import fetch_reference
from .ai_annotation import Manager as AIAnnotationManager
from .auto_annotation import Service as AutomaticAnnotationService
from .ai_benchmark import Service as BenchmarkService

STATIC = Path(__file__).parent / "static"


def migration_status(path):
    try:
        job = json.loads((Path(path).parent / "migration_job.json").read_text(encoding="utf-8"))
        if job.get("state") in ("running", "generating_report"):
            heartbeat = job.get("updated_at")
            if not heartbeat or (datetime.now(timezone.utc) - datetime.fromisoformat(heartbeat)).total_seconds() > 120:
                job = {**job, "state": "unconfirmed", "reported_state": job["state"]}
        return job
    except (OSError, ValueError, TypeError):
        return {"state": "unknown"}


class App:
    def __init__(self, path, manager_factory=None):
        self.path = Path(path)
        self.started_at = db.now()
        with closing(db.connect(path)) as c:
            db.initialize(c)
        self.jobs = {}
        self.pool = ThreadPoolExecutor(max_workers=1)
        self.lock = threading.Lock()
        self.stats_cache = None
        self.stats_at = 0.
        self.ai = (manager_factory or AIAnnotationManager)(self.path)
        self.auto = AutomaticAnnotationService(self.ai)
        self.benchmark = BenchmarkService(self.ai)

    def stats(self):
        # One count scan per interval, even when several pages request it.
        with self.lock:
            if self.stats_cache is None or time.monotonic() - self.stats_at >= 30:
                with closing(db.connect(self.path, timeout=.25)) as c:
                    c.execute("BEGIN")
                    self.stats_cache = db.stats(c)
                self.stats_at = time.monotonic()
            return {**self.stats_cache, "migration_job": migration_status(self.path)}

    def export(self, payload, preview=False):
        if not isinstance(payload,dict) or not isinstance(payload.get('filters') or {},dict):raise ValueError('导出参数及基础筛选应为对象')
        export_options=validate_export_options(payload.get('export_options'))
        if payload.get("format") not in ("jsonl", "parquet", "lmdb", "atompeak"):
            raise ValueError("导出格式无效")
        training_options = None
        if payload.get('format') == 'atompeak':
            options = payload.get('training_options') or {}
            if not isinstance(options,dict):
                raise ValueError('训练导出参数应为对象')
            allowed = {'seed','valid_ratio','test_ratio','max_seq_len','max_target_len','use_solvent','use_temperature','use_field_strength'}
            if set(options) - allowed:
                raise ValueError('训练导出参数无效')
            if any(k in options and not isinstance(options[k],bool) for k in ('use_solvent','use_temperature','use_field_strength')):
                raise ValueError('条件开关必须是布尔值')
            training_options = {'seed':int(options.get('seed',42)),
                'valid_ratio':float(options.get('valid_ratio',.1)),'test_ratio':float(options.get('test_ratio',.1)),
                'max_seq_len':int(options.get('max_seq_len',1024)),'max_target_len':int(options.get('max_target_len',512)),
                'use_solvent':bool(options.get('use_solvent',False)),
                'use_temperature':bool(options.get('use_temperature',False)),
                'use_field_strength':bool(options.get('use_field_strength',False))}
            if min(training_options['max_seq_len'],training_options['max_target_len'])<8 or max(training_options['max_seq_len'],training_options['max_target_len'])>8192:
                raise ValueError('训练序列长度应在 8–8192 之间')
            from atompeak_data.schema import split_for_connectivity
            split_for_connectivity('validate',training_options['seed'],training_options['valid_ratio'],training_options['test_ratio'])
        job = uuid.uuid4().hex
        self.jobs[job] = {"state": "running"}
        def work():
            try:
                def update(info):
                    self.jobs[job] = {'state':'running', **info}
                if preview:
                    result=preview_database(self.path,payload.get('filters'),payload.get('molecule_id'),export_options,update,require_quality=payload['format']=='atompeak')
                else:
                    result = export_database(self.path, self.path.parent / "exports", payload["format"], payload.get("filters"), payload.get("molecule_id"), training_options=training_options, progress=update,export_options=export_options)
                self.jobs[job] = {"state": "done", **result}
            except Exception as exc:
                self.jobs[job] = {"state": "error", "error": str(exc)}
        self.pool.submit(work)
        return {"job_id": job}


class Handler(BaseHTTPRequestHandler):
    server_version = "NMRDatabase/1"

    def disconnected(self):
        # A cancelled browser request is not an application/database failure.
        # Never send another response on this socket (including error replies).
        self.close_connection = True

    def handle_one_request(self):
        try:
            super().handle_one_request()
        except (ConnectionAbortedError, ConnectionResetError, BrokenPipeError):
            self.disconnected()

    def send(self, data, kind="application/json; charset=utf-8", status=200, attachment=None):
        body = data if isinstance(data, bytes) else json.dumps(data, ensure_ascii=False).encode()
        self.send_response(status)
        self.send_header("Content-Type", kind)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        if attachment:
            self.send_header("Content-Disposition", f'attachment; filename="{attachment}"')
        try:
            self.end_headers()
            self.wfile.write(body)
        except (ConnectionAbortedError, ConnectionResetError, BrokenPipeError):
            self.disconnected()

    def do_GET(self):
        parsed = urlparse(self.path)
        path = unquote(parsed.path)
        params = {k: v[-1] for k, v in parse_qs(parsed.query).items()}
        try:
            if path.startswith("/api/export/"):
                job = path.split("/")[3]
                value = self.server.app.jobs.get(job)
                if not value:
                    return self.send({"error": "没有该导出任务"}, status=404)
                if (path.endswith('/download') or path.endswith('/manifest')) and value['state']=='done':
                    filename=value.get('manifest') if path.endswith('/manifest') else value.get('file')
                    if not filename:return self.send({'error':'此任务只有预览结果，没有可下载文件'},status=400)
                    file = Path(filename)
                    self.send_response(200)
                    self.send_header("Content-Type", "application/octet-stream")
                    self.send_header("Content-Length", str(file.stat().st_size))
                    self.send_header("Content-Disposition", f'attachment; filename="{file.name}"')
                    self.end_headers()
                    with file.open("rb") as f:
                        for chunk in iter(lambda: f.read(1024*1024), b""):
                            self.wfile.write(chunk)
                    return
                return self.send(value)
            if path.startswith("/api/"):
                if path == '/api/ai/credentials':
                    return self.send(self.server.app.ai.vault.status())
                if path == '/api/health':return self.send({'status':'ok','service':'atomnmr','pid':os.getpid(),'started_at':self.server.app.started_at})
                if path == '/api/benchmark/profiles':return self.send(self.server.app.benchmark.profiles())
                if path == '/api/benchmark/sets':return self.send(self.server.app.benchmark.sets())
                if path.startswith('/api/benchmark/sets/'):return self.send(self.server.app.benchmark.sample_set(path.rsplit('/',1)[1]))
                if path == '/api/benchmark/runs':return self.send(self.server.app.benchmark.history())
                if path.startswith('/api/benchmark/runs/'):
                    identifier=path.rsplit('/',1)[1];return self.send(self.server.app.benchmark.detail(identifier))
                if path.startswith('/api/benchmark/cases/'):return self.send(self.server.app.benchmark.case(path.rsplit('/',1)[1]))
                if path == '/api/benchmark/compare':return self.send(self.server.app.benchmark.compare(params.get('ids','').split(',')))
                if path == '/api/benchmark/export':return self.send(self.server.app.benchmark.export(params.get('ids','').split(',')),attachment='benchmark_comparison.json')
                if path == '/api/auto/timing':return self.send(self.server.app.auto.timing(params.get('id')))
                if path == '/api/auto/execution':return self.send(self.server.app.auto.execution())
                if path == '/api/auto/jobs':
                    return self.send(self.server.app.auto.status(params.get('id'),params.get('result_state')))
                if path.startswith('/api/auto/items/'):
                    return self.send(self.server.app.auto.item(path.rsplit('/',1)[1]))
                if path == '/api/auto/results':
                    return self.send(self.server.app.auto.molecule_results(params.get('molecule_id')))
                if path.startswith('/api/auto/results/'):
                    return self.send(self.server.app.auto.result(path.removeprefix('/api/auto/results/')))
                if path.startswith('/api/ai/'):
                    if path == '/api/ai/settings':
                        return self.send(self.server.app.ai.settings())
                    if path == '/api/ai/runs':
                        return self.send(self.server.app.ai.history())
                    parts = path.strip('/').split('/')
                    if len(parts) >= 4 and parts[2] == 'runs':
                        if len(parts) == 6 and parts[4] == 'items':
                            return self.send(self.server.app.ai.item(parts[3],parts[5]))
                        if len(parts) == 5 and parts[4] == 'download':
                            result = self.server.app.ai.detail(parts[3])
                            result['results'] = [self.server.app.ai.item(parts[3],r['id']) for r in result['items']]
                            return self.send(result,attachment='ai_assignment_'+parts[3]+'.json')
                        if len(parts) == 4:
                            return self.send(self.server.app.ai.detail(parts[3]))
                    return self.send({'error':'AI 接口不存在'},status=404)
                if path == '/api/workload':
                    with closing(db.connect(self.server.app.path, timeout=.25)) as c:
                        c.execute('BEGIN')
                        return self.send(workload.summarize(c, params.get('date'), params.get('days', 14)))
                if path == "/api/stats":
                    return self.send(self.server.app.stats())
                with closing(db.connect(self.server.app.path, timeout=.25)) as c:
                    if path == "/api/molecules":
                        started = time.monotonic()
                        c.set_progress_handler(lambda: int(time.monotonic() - started > 10), 10000)
                        return self.send(db.listing(c, params))
                    if path.startswith("/api/molecule/"):
                        mid = path.split("/")[3]
                        if path.endswith('/numbering'):
                            c.execute('BEGIN')
                            return self.send(numbering.preview(c, mid, params.get('spectrum_id')))
                        if path.endswith('/management'):
                            c.execute('BEGIN')
                            return self.send(curation.details(c, mid))
                        if path.endswith("/svg"):
                            r = c.execute("SELECT smiles FROM molecules WHERE id=?", (mid,)).fetchone()
                            if not r:
                                raise KeyError(mid)
                            selected = [int(x) for x in params.get("atoms", "").split(",") if x]
                            mol = db.Chem.MolFromSmiles(r[0])
                            if any(x < 0 or x >= mol.GetNumAtoms() for x in selected):
                                raise ValueError("原子编号错误")
                            return self.send(db.molecule_svg(r[0], selected, params.get('view', 'auto')).encode(), "image/svg+xml")
                        return self.send(db.record(c, mid, include_archived=True))
                return self.send({"error": "接口不存在"}, status=404)
            if path == "/nmr_reference.js":
                file = db.ROOT / "annotation_tool_v2/static/nmr_reference.js"
            elif path.startswith("/figures/"):
                base = (db.ROOT / "docs/database_v1_figures").resolve()
                file = (base / path.removeprefix("/figures/")).resolve()
                if not file.is_relative_to(base):
                    raise ValueError("Invalid path")
            else:
                if path == '/ai.html' and not params.get('run_id') and not params.get('history'):
                    self.send_response(302);self.send_header('Location','/benchmark.html');self.send_header('Content-Length','0');self.end_headers();return
                file = (STATIC / ("index.html" if path == "/" else path.lstrip("/"))).resolve()
                if not file.is_relative_to(STATIC.resolve()):
                    raise ValueError("Invalid path")
            if not file.is_file():
                return self.send({"error": "文件不存在"}, status=404)
            return self.send(file.read_bytes(), mimetypes.guess_type(file.name)[0] or "application/octet-stream")
        except (ConnectionAbortedError, ConnectionResetError, BrokenPipeError):
            self.disconnected()
        except sqlite3.OperationalError as exc:
            code = getattr(exc, 'sqlite_errorcode', 0) & 255
            if code == sqlite3.SQLITE_INTERRUPT:
                return self.send({'error': '检索耗时过长，请缩小筛选范围或使用完整编号检索。'}, status=503)
            if db.is_busy_error(exc):
                return self.send({'error': '数据库正在短暂写入，请稍后重新读取。'}, status=503)
            self.send({'error': str(exc)}, status=500)
        except KeyError:
            self.send({"error": "记录不存在"}, status=404)
        except (ValueError, TypeError) as exc:
            self.send({"error": str(exc)}, status=400)
        except Exception as exc:
            self.send({"error": str(exc)}, status=500)

    def do_POST(self):
        try:
            origin = self.headers.get("Origin")
            expected = "http://" + self.headers.get("Host", "")
            if origin and origin != expected:
                return self.send({"error": "仅接受本页面提交"}, status=403)
            if not self.headers.get("Content-Type", "").startswith("application/json"):
                return self.send({"error": "只接受 JSON"}, status=415)
            size = int(self.headers.get("Content-Length", "0"))
            if not 0 < size <= 8*1024*1024:
                raise ValueError("请求体大小无效")
            payload = json.loads(self.rfile.read(size))
            if self.path.startswith('/api/benchmark/'):
                if not isinstance(payload,dict):raise ValueError('请求应为JSON对象')
                bench=self.server.app.benchmark
                if self.path == '/api/benchmark/profiles':return self.send(bench.save_profile(payload))
                if self.path == '/api/benchmark/sets':return self.send(bench.build_set(payload))
                if self.path == '/api/benchmark/runs':return self.send(bench.start(payload))
                parts=self.path.strip('/').split('/')
                if len(parts)==5:
                    if parts[2]=='profiles' and parts[4]=='test':return self.send(bench.test_profile(parts[3]))
                    if parts[2]=='runs' and parts[4]=='pause':return self.send(bench.pause(parts[3]))
                    if parts[2]=='runs' and parts[4]=='resume':return self.send(bench.resume(parts[3]))
                return self.send({'error':'测试接口不存在'},status=404)
            if self.path.startswith('/api/ai/'):
                if not isinstance(payload,dict):
                    raise ValueError('请求应为 JSON 对象')
                ai = self.server.app.ai
                if self.path == '/api/ai/credentials':
                    return self.send(ai.manage_credentials(payload))
                if self.path == '/api/ai/settings':
                    return self.send(ai.configure(payload))
                if self.path == '/api/ai/test':
                    return self.send(ai.test())
                if self.path == '/api/ai/preview':
                    return self.send(ai.preview(payload))
                if self.path == '/api/ai/runs':
                    return self.send(ai.start(payload))
                parts = self.path.strip('/').split('/')
                if len(parts) == 5 and parts[2] == 'runs':
                    if parts[4] == 'stop':
                        return self.send(ai.stop(parts[3]))
                    if parts[4] == 'resume':
                        return self.send(ai.resume(parts[3],payload.get('retry_errors') is True))
                    if parts[4] == 'retry-current':
                        return self.send(ai.retry_with_settings(parts[3]))
                return self.send({'error':'AI 接口不存在'},status=404)
            if self.path.startswith('/api/auto/'):
                if not isinstance(payload,dict):raise ValueError('请求应为 JSON 对象')
                auto=self.server.app.auto
                if self.path == '/api/auto/jobs':return self.send(auto.start(payload))
                if self.path == '/api/auto/execution':return self.send(auto.execution(payload))
                if self.path == '/api/auto/endpoint':return self.send(self.server.app.benchmark.adopt_for_production(payload.get('profile_id')))
                if self.path == '/api/auto/import-history':
                    result=auto.import_history();self.server.app.stats_at=0;return self.send(result)
                parts=self.path.strip('/').split('/')
                if len(parts)==5 and parts[2]=='jobs':
                    if parts[4]=='delete':return self.send(auto.delete(parts[3]))
                    if parts[4]=='pause':return self.send(auto.pause(parts[3]))
                    if parts[4]=='resume':return self.send(auto.resume(parts[3],payload.get('use_current_limits') is True))
                    if parts[4]=='retry':return self.send(auto.retry_failed(parts[3],payload.get('item_id'),payload.get('mode','failed')))
                    if parts[4]=='analyze':return self.send(auto.analyze(parts[3]))
                return self.send({'error':'自动标注接口不存在'},status=404)
            if self.path == '/api/export/preview':
                return self.send(self.server.app.export(payload,preview=True))
            if self.path == "/api/export":
                return self.send(self.server.app.export(payload))
            if self.path.startswith("/api/molecule/"):
                mid = self.path.split("/")[3]
                if self.path.endswith('/numbering'):
                    with closing(db.connect(self.server.app.path, timeout=.25)) as c:
                        value = numbering.review(c, mid, payload)
                    self.server.app.stats_at = 0
                    return self.send(value)
                if self.path.endswith('/management'):
                    reference = None
                    if payload.get('action') == 'import_bmrb':
                        with closing(db.connect(self.server.app.path, timeout=.25)) as c:
                            r = c.execute('SELECT smiles,revision FROM molecules WHERE id=?', (mid,)).fetchone()
                            if not r:
                                raise KeyError(mid)
                            if payload.get('revision') != r['revision']:
                                raise RuntimeError('记录已更新，请重新打开管理面板。')
                        # No database lock while waiting for the official source.
                        reference = fetch_reference(payload.get('entry_id'), r['smiles'])
                    with closing(db.connect(self.server.app.path, timeout=.25)) as c:
                        value = curation.action(c, mid, payload, reference)
                    self.server.app.stats_at = 0
                    return self.send(value)
                with closing(db.connect(self.server.app.path, timeout=.25)) as c:
                    value = db.save_annotations(c, mid, payload)
                # Invalidate without running global counts on the save path.
                self.server.app.stats_at = 0
                return self.send(value)
            return self.send({"error": "接口不存在"}, status=404)
        except (ConnectionAbortedError, ConnectionResetError, BrokenPipeError):
            self.disconnected()
        except sqlite3.OperationalError as exc:
            if db.is_busy_error(exc):
                return self.send({"error": "后台正在短暂写入，当前修改尚未提交，请稍后重试。", "retryable": True}, status=503)
            self.send({"error": str(exc)}, status=500)
        except RuntimeError as exc:
            self.send({"error": str(exc)}, status=409)
        except (ValueError, KeyError, TypeError) as exc:
            self.send({"error": str(exc)}, status=400)
        except Exception as exc:
            self.send({"error": str(exc)}, status=500)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--db", type=Path, default=db.DEFAULT_DB)
    p.add_argument("--port", type=int, default=8767)
    p.add_argument("--open", action="store_true")
    args = p.parse_args()
    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    server.app = App(args.db)
    print(f"NMR experimental database: http://127.0.0.1:{args.port}/", flush=True)
    if args.open:
        webbrowser.open(f"http://127.0.0.1:{args.port}/")
    try:
        server.serve_forever()
    finally:
        server.server_close()
        server.app.pool.shutdown(wait=True)
        server.app.ai.shutdown()


if __name__ == "__main__":
    main()
