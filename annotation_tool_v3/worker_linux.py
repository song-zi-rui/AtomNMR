"""Run the production annotation queue on Ubuntu without a web server."""
import argparse
import json
import os
import signal
import sqlite3
import threading
import time
from contextlib import closing
from pathlib import Path

from . import core as db, api_tls
from .auto_annotation import Service
from .paths import artifacts_root, companion_database
from .server_linux import EnvironmentManager, check_api


TASK_NAME = '组内 API · 服务器持续标注'
BASE_URL = 'https://10.8.130.88/v1'
MODEL = 'deepseek-v4.1-flash-fp8'
EXECUTION = {'parallelism': 3, 'parallel_strategy': 'fixed'}
PRODUCTION_SCOPES = ('pending', 'unlabeled', 'partial')
ATTENTION_EXIT = 75
CONFIG_EXIT = 78


def group_settings():
    return {
        'base_url': BASE_URL, 'model': MODEL, 'api_adapter': 'deepseek-local',
        'thinking': True, 'reasoning_effort': 'auto', 'max_tokens': 65536, 'timeout': 1400,
        'wire_format': 'compact-v3', 'assignment_mode': 'assisted', 'json_mode': True, 'stream': True,
        'tls_policy': 'legacy_ca',
        'ca_file': str(Path(__file__).resolve().parent / 'certificates' / 'ai-cluster-ca.crt'),
        'ca_sha256': 'b09dc7b6f3c4207a659466b6dce7d23775d7b4d24144f457346456dab7e2eaa7',
        'tls_origin': 'https://10.8.130.88:443',
    }


def job_row(conn, identifier=None):
    if identifier:
        row = conn.execute('SELECT * FROM auto_jobs WHERE id=? AND deleted_at IS NULL', (identifier,)).fetchone()
        if not row or row['scope'] not in PRODUCTION_SCOPES:
            raise ValueError('找不到指定的正式自动标注任务')
        return row
    return conn.execute("""SELECT * FROM auto_jobs WHERE name=? AND scope='pending' AND deleted_at IS NULL
        ORDER BY created_at DESC,id DESC LIMIT 1""", (TASK_NAME,)).fetchone()


def status(path, identifier=None):
    history = companion_database(path, 'ai_annotations.sqlite')
    if not history.is_file():
        return {'state': 'not_started'}
    with closing(sqlite3.connect(history.as_uri() + '?mode=ro', uri=True, timeout=5)) as conn:
        conn.row_factory = sqlite3.Row
        if not conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='auto_jobs'").fetchone():
            return {'state': 'not_started'}
        row = job_row(conn, identifier)
        if not row:
            return {'state': 'not_started'}
        config = json.loads(row['config_json'])
        counts = dict(conn.execute(
            'SELECT state,count(*) FROM auto_items WHERE job_id=? GROUP BY state', (row['id'],)))
        return {'id': row['id'], 'state': row['state'], 'model': config.get('model'),
                'parallelism': config.get('parallelism'),
                'done': counts.get('done', 0), 'failed': counts.get('error', 0), 'calling': counts.get('calling', 0),
                'cooldown_s': max(0, int(row['retry_cooldown_until'] - time.time())) if 'retry_cooldown_until' in row.keys() else 0,
                'updated_at': row['updated_at'], 'error': row['error']}


def recent_errors(path, identifier=None):
    history = companion_database(path, 'ai_annotations.sqlite')
    if not history.is_file():
        return {'state': 'not_started', 'errors': []}
    with closing(sqlite3.connect(history.as_uri() + '?mode=ro', uri=True, timeout=5)) as conn:
        conn.row_factory = sqlite3.Row
        if not conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='auto_jobs'").fetchone():
            return {'state': 'not_started', 'errors': []}
        job = job_row(conn, identifier)
        if not job:
            return {'state': 'not_started', 'errors': []}
        rows = conn.execute("""SELECT i.id,i.molecule_id,i.formula,i.updated_at AS at,i.elapsed,
                i.error,i.diagnostics_json,NULL AS archived_attempt
            FROM auto_items i WHERE i.job_id=? AND i.state='error'
            UNION ALL
            SELECT i.id,i.molecule_id,i.formula,a.archived_at,a.elapsed,a.error,a.diagnostics_json,a.attempt_no
            FROM auto_item_attempts a JOIN auto_items i ON i.id=a.item_id WHERE i.job_id=? AND a.state='error'
            ORDER BY at DESC,id DESC LIMIT 30""", (job['id'], job['id']))
        failures = []
        diagnostic_fields = ('http_status', 'provider_status', 'error_kind', 'failure_stage',
                             'request_id', 'request_outcome_uncertain', 'retry_after', 'automatic_retry')
        for row in rows:
            entry = dict(row)
            diagnostics = json.loads(entry.pop('diagnostics_json') or 'null') or {}
            provider = diagnostics.get('provider_error') or {}
            entry.update({key: diagnostics[key] for key in diagnostic_fields if diagnostics.get(key) is not None})
            if provider:
                entry['provider_error'] = provider
            failures.append(entry)
        return {'id': job['id'], 'state': job['state'], 'errors': failures}


def report(value):
    print(json.dumps(value, ensure_ascii=False), flush=True)


def run(path, identifier, stop, *, prepare_resume=False):
    with closing(db.connect(path)) as conn:
        db.initialize(conn)
    manager = EnvironmentManager(path, base_url=BASE_URL)
    try:
        service = Service(manager)
        if stop.is_set():
            return 0
        with closing(sqlite3.connect(manager.path, timeout=5)) as conn:
            conn.row_factory = sqlite3.Row
            row = job_row(conn, identifier)
        if prepare_resume:
            if not row:
                raise ValueError('找不到要恢复的正式自动标注任务')
            count = service.prepare_resume_after_fix(row['id'])
            report({'id': row['id'], 'state': 'paused', 'legacy_stream_errors_queued': count})
            return 0
        if row and row['state'] in ('completed', 'attention'):
            report(status(path, row['id']))
            return ATTENTION_EXIT if row['state'] == 'attention' else 0
        manager.configure(group_settings())
        service.execution(EXECUTION)
        if row:
            service.resume(row['id'], use_current_limits=True)
            job_id = row['id']
        else:
            job_id = service.start({'name': TASK_NAME, 'scope': 'pending'})['job_id']
        report(status(path, job_id))
        last_report = time.monotonic()
        while manager.active:
            if stop.is_set() and not manager.cancel.is_set():
                try:
                    service.pause(job_id)
                except ValueError:
                    if manager.active:
                        raise
                print('已停止派发，等待在途请求保存', flush=True)
            if time.monotonic() - last_report >= 60:
                report(status(path, job_id))
                last_report = time.monotonic()
            time.sleep(1)
        final = status(path, job_id)
        report(final)
        return ATTENTION_EXIT if final['state'] == 'attention' else 0
    finally:
        manager.shutdown()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', type=Path, required=True)
    actions = parser.add_mutually_exclusive_group()
    actions.add_argument('--check-api', action='store_true')
    actions.add_argument('--status', action='store_true')
    actions.add_argument('--errors', action='store_true')
    actions.add_argument('--resume-after-fix', action='store_true')
    args = parser.parse_args()
    identifier = os.environ.get('ATOMNMR_RESUME_JOB_ID') or None
    path = args.db.resolve()
    try:
        if args.check_api:
            info = check_api(group_settings())
            report(info)
            return 0 if info['selected_model_found'] else CONFIG_EXIT
        if args.status:
            report(status(path, identifier))
            return 0
        if args.errors:
            report(recent_errors(path, identifier))
            return 0
        if not path.is_file():
            raise ValueError('找不到标注主库 assignments.sqlite')
        api_tls.validate_config(group_settings())
        import fcntl
        stop = threading.Event()
        for signum in (signal.SIGTERM, signal.SIGINT):
            signal.signal(signum, lambda *_: stop.set())
        with (artifacts_root(path) / 'annotation-server.lock').open('a+b') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise ValueError('该数据库已有标注服务在运行') from exc
            return run(path, identifier, stop, prepare_resume=args.resume_after_fix)
    except ValueError as exc:
        report({'state': 'configuration_error', 'error': str(exc)})
        return CONFIG_EXIT


if __name__ == '__main__':
    raise SystemExit(main())
