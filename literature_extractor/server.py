"""Local browser interface for the literature extraction folder workflow."""

import argparse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from threading import Lock
from urllib.parse import quote, urlsplit
import json
import mimetypes
import os
import re
import shutil
import webbrowser

from rdkit import Chem

from .batch import BatchQueue, load
from .extract import diagram, digest, source_paths
from .grouping import grouping_plan
from .search import CATEGORIES, search

ROOT = Path(__file__).resolve().parents[1]
STATIC = Path(__file__).with_name('static')
OUTPUT = ROOT / 'outputs' / 'literature_extraction'
FILES = {'/': 'index.html', '/style.css': 'style.css', '/app.js': 'app.js'}
JOB = re.compile(r'/api/jobs/([0-9a-f]{12})(?:/(report|database|pdf/\d+))?')


def list_folders(folder):
    if not folder.is_absolute():
        raise ValueError('请输入完整的文件夹路径。')
    try:
        folder = folder.resolve()
        with os.scandir(folder) as entries:
            children = sorted((entry.name for entry in entries if entry.is_dir()), key=str.casefold)
    except OSError as exc:
        raise ValueError(f'无法读取该文件夹。{exc.strerror or exc}') from exc
    return {'folder': str(folder), 'parent': str(folder.parent) if folder.parent != folder else None,
            'children': [{'name': name, 'path': str(folder / name)} for name in children]}


def result_links(job_id, report, route_prefix=''):
    prefix = f'{route_prefix}/api/jobs/{job_id}'
    report['report_url'] = f'{prefix}/report'
    if 'database' in report:
        report['database_url'] = f'{prefix}/database'
    for table in report['tables']:
        table['pdf_url'] = f"{prefix}/pdf/{table['file_id']}#page={table['page']}"
        if 'supplement_file_id' in table:
            table['supplement_url'] = f"{prefix}/pdf/{table['supplement_file_id']}#page={table['supplement_page']}"
    return report


def report_pdf(root, file):
    path = root / file.get('relative_path', file['name'])
    if path.is_file() and path.resolve().is_relative_to(root):
        return path
    # Grouping may have moved a loose source from the location recorded in an earlier report.
    candidates = [p for p in root.iterdir() if p.is_dir() and not p.is_symlink() and p.resolve().parent == root]
    for child in candidates:
        same_name = child / path.name
        if not same_name.is_symlink() and same_name.is_file() and digest(same_name) == file['sha256']:
            return same_name
    for child in candidates:
        for candidate in source_paths(child):
            if candidate.suffix.lower() == '.pdf' and digest(candidate) == file['sha256']:
                return candidate
    raise FileNotFoundError(path)


class ExtractionService:
    def __init__(self, output, *, importer=None, prefix=''):
        self.output = output
        self.prefix = prefix
        self.queue = BatchQueue(output, importer)
        self.operations = Lock()

    def get(self, key, *, asset=False):
        if self.queue.state_file(key).is_file():
            state = self.queue.state(key)
            state['saved'] = len(state.pop('molecules'))
            state.pop('records')
            raw = state.pop('result_file', None)
            if not raw or (state['state'] == 'running' and not asset):
                return state
            path = Path(raw)
        else:
            path = next(self.output.glob(f'*_{key}/extraction_report.json'), None)
            if path is None:
                raise KeyError(key)
            state = {'state': 'complete'}
        report = load(path)
        report['output_folder'] = str(path.parent)
        if 'database' in report:
            report['database'] = str(path.with_name('assignments.sqlite'))
        if not asset:
            for table in report['tables']:
                if table.get('smiles') and table.get('mapping'):
                    mol = Chem.MolFromSmiles(table['smiles'])
                    if mol is not None:
                        mapping = {int(tag): idx for tag, idx in table['mapping'].items()}
                        display = {int(tag): label for tag, label in table.get('atom_display', {}).items()}
                        table['diagram'], table['points'] = diagram(mol, mapping, display)
        return {**state, 'folder': report['folder'], 'result': result_links(key, report, self.prefix)}


class Handler(BaseHTTPRequestHandler):
    def allowed_host(self):
        port = self.server.server_port
        return (self.client_address[0] == '127.0.0.1' and
                self.headers.get('Host') in (f'127.0.0.1:{port}', f'localhost:{port}'))

    def send_json(self, value, status=200):
        data = json.dumps(value, ensure_ascii=False).encode('utf-8')
        self.send_response(status)
        self.send_header('Content-Type', 'application/json; charset=utf-8')
        self.send_header('Content-Length', str(len(data)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(data)

    def send_file(self, path, download=False, static=False):
        with path.open('rb') as handle:
            self.send_response(200)
            self.send_header('Content-Type', mimetypes.guess_type(path.name)[0] or 'application/octet-stream')
            self.send_header('Content-Length', str(path.stat().st_size))
            self.send_header('X-Content-Type-Options', 'nosniff')
            if static:
                self.send_header('Cache-Control', 'no-store')
            if download:
                self.send_header('Content-Disposition', f"attachment; filename*=UTF-8''{quote(path.name)}")
            self.end_headers()
            shutil.copyfileobj(handle, self.wfile)

    def do_GET(self):
        self.get_literature(urlsplit(self.path).path, self.server.service)

    def get_literature(self, path, service):
        if not self.allowed_host():
            self.send_json({'error': '此界面仅供本机访问。'}, 403)
            return
        try:
            if path in FILES:
                self.send_file(STATIC / FILES[path], static=True)
                return
            route = JOB.fullmatch(path)
            if not route:
                raise KeyError(path)
            job_id, asset = route.groups()
            job = service.get(job_id, asset=bool(asset))
            if not asset:
                self.send_json(job)
                return
            result = job['result']
            if asset == 'report':
                self.send_file(Path(result['output_folder']) / 'extraction_report.json', download=True)
            elif asset == 'database':
                self.send_file(Path(result['database']), download=True)
            else:
                file = result['files'][int(asset.split('/')[1])]
                root = Path(job['folder']).resolve()
                self.send_file(report_pdf(root, file))
        except (KeyError, IndexError, FileNotFoundError):
            self.send_json({'error': '找不到该任务或文件。'}, 404)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def do_POST(self):
        self.post_literature(self.path, self.server.service)

    def post_literature(self, route, service):
        host = self.headers.get('Host')
        origin = self.headers.get('Origin')
        if not self.allowed_host() or origin not in (None, f'http://{host}'):
            self.send_json({'error': '此界面仅接受本机请求。'}, 403)
            return
        try:
            if self.headers.get_content_type() != 'application/json':
                raise ValueError('请求格式有误。')
            size = int(self.headers.get('Content-Length', '0'))
            if not 0 < size < 16000:
                raise ValueError('请求大小有误。')
            body = json.loads(self.rfile.read(size))
            raw_path = body.get('folder', '').strip().strip('"')
            folder = Path(raw_path).expanduser()
            if route == '/api/folders':
                self.send_json(list_folders(folder if raw_path else ROOT))
            elif route == '/api/search':
                name = body.get('name', '').strip()
                limit = body.get('limit', 8)
                categories = body.get('categories') or tuple(CATEGORIES)
                if len(name) > 150 or type(limit) is not int or not 1 <= limit <= 20 or any(c not in CATEGORIES for c in categories):
                    raise ValueError('请检查检索词、分类和每类数量。')
                self.send_json(search(name, limit, categories))
            elif route == '/api/pause':
                key = body.get('id', '')
                if not re.fullmatch(r'[0-9a-f]{12}', key):
                    raise ValueError('任务编号无效。')
                service.queue.pause(key)
                self.send_json({'ok': True})
            elif route in ('/api/extract', '/api/latest', '/api/group'):
                if not folder.is_absolute() or not folder.is_dir():
                    raise ValueError('请输入或选择一个有效的文件夹路径。')
                folder = folder.resolve()
                if route == '/api/latest':
                    self.send_json({'id': service.queue.latest(folder)})
                elif route == '/api/group':
                    with service.operations:
                        if service.queue.active:
                            raise ValueError('请先暂停提取，再对文件进行分组。')
                        self.send_json(grouping_plan(folder, apply=body.get('apply') is True))
                else:
                    limit = body.get('limit')
                    resume = body.get('resume')
                    if limit is not None and (type(limit) is not int or limit < 1):
                        raise ValueError('目标分子数应为正整数，留空表示不限。')
                    if resume and not re.fullmatch(r'[0-9a-f]{12}', resume):
                        raise ValueError('任务编号无效。')
                    with service.operations:
                        key = service.queue.start(folder, limit, resume)
                    self.send_json({'id': key}, 202)
            else:
                self.send_json({'error': '找不到该操作。'}, 404)
        except (BrokenPipeError, ConnectionResetError):
            pass
        except (ValueError, OSError, TypeError, AttributeError) as exc:
            self.send_json({'error': str(exc)}, 400)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port', type=int, default=8774)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--open', action='store_true')
    args = parser.parse_args()
    with ThreadingHTTPServer(('127.0.0.1', args.port), Handler) as server:
        server.service = ExtractionService(args.output.resolve())
        url = f'http://127.0.0.1:{server.server_port}/'
        print(f'AtomNMR 文献提取 · {url}', flush=True)
        if args.open:
            webbrowser.open(url)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
        finally:
            server.service.queue.shutdown()


if __name__ == '__main__':
    main()
