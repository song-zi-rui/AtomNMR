"""Durable folder queue for local literature extraction."""

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from pathlib import Path
from threading import Lock
import json
import uuid

from .dedup import source_bundle, spectrum_keys
from .extract import collect, db, digest, read_pdf, save_results, source_paths
from .importer import CHANGES

PARSER_VERSION = 1


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f'.{path.name}.{uuid.uuid4().hex}.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False), encoding='utf-8')
    temporary.replace(path)


def folders(root):
    children = [p for p in sorted(root.iterdir(), key=lambda p: p.name.casefold())
                if p.is_dir() and not p.is_symlink() and p.resolve().parent == root
                and any(f.suffix.lower() == '.pdf' for f in source_paths(p))]
    return children or ([root] if any(p.suffix.lower() == '.pdf' for p in source_paths(root)) else [])


def snapshot(folder):
    paths = source_paths(folder)
    stamp = file_stamp(paths)
    sources = {p.name: digest(p) for p in paths}
    key = sha256(json.dumps([PARSER_VERSION, str(folder), sources], sort_keys=True).encode()).hexdigest()
    return key, sources, stamp


def file_stamp(paths):
    stamps = {}
    for path in paths:
        stat = path.stat()
        stamps[path.name] = (stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino)
    return stamps


class BatchQueue:
    def __init__(self, output, importer=None):
        self.output = output
        self.importer = importer
        self.lock = Lock()
        self.pool = ThreadPoolExecutor(max_workers=1)
        self.active = None
        self.progress = {}

    def home(self, key):
        return self.output / 'batches' / key

    def state_file(self, key):
        return self.home(key) / 'batch.json'

    def state(self, key):
        with self.lock:
            state = load(self.state_file(key))
            if state['state'] == 'running' and self.active != key:
                state.update(state='paused', pause=False, reason='服务退出，任务可继续。')
            return {**state, **self.progress} if self.active == key else state

    def latest(self, folder):
        paths = sorted(self.output.glob('batches/*/batch.json'), key=lambda p: p.stat().st_mtime, reverse=True)
        return next((path.parent.name for path in paths if load(path)['folder'] == str(folder)), None)

    def start(self, folder, limit, resume=None):
        with self.lock:
            if self.active:
                raise ValueError('当前已有提取任务运行。')
            key = resume or uuid.uuid4().hex[:12]
            path = self.state_file(key)
            if resume:
                job = load(path)
                if job['state'] not in ('paused', 'running', 'error') or job['folder'] != str(folder):
                    raise ValueError('该任务无法继续，请开始新批次。')
                if limit is not None and limit <= len(job['molecules']):
                    raise ValueError('继续提取时，分子目标数应大于已经保存的数量。')
            else:
                pending = self.latest(folder)
                if pending and load(self.state_file(pending))['state'] in ('paused', 'running', 'error'):
                    raise ValueError('该目录有可继续的批次。请使用「查看上次任务」。')
                job = {'folder': str(folder), 'records': {}, 'molecules': [], 'errors': [], 'started_at': db.now()}
            job.update(state='running', limit=limit, pause=False, reason='查找论文文件夹')
            write(path, job)
            self.active = key
            self.progress = {}
            self.pool.submit(self.work, key)
            return key

    def pause(self, key):
        with self.lock:
            if not key or self.active != key:
                return
            path = self.state_file(key)
            state = load(path)
            state['pause'] = True
            write(path, state)

    def shutdown(self):
        self.pause(self.active)
        self.pool.shutdown()

    def statistics(self):
        result = self.importer.statistics()
        imports = result.pop('batches')
        batches, totals = [], defaultdict(int)
        paths = sorted(self.output.glob('batches/*/batch.json'), key=lambda p: p.stat().st_mtime, reverse=True)
        for path in paths:
            key = path.parent.name
            job = self.state(key)
            files = sum(r.get('pdf_count', 0) for r in job['records'].values())
            review_tables = sum(r.get('review_tables', 0) for r in job['records'].values())
            review_files = sum(r.get('review_files', 0) for r in job['records'].values())
            totals['pdf_count'] += files
            totals['review_tables'] += review_tables
            totals['review_files'] += review_files
            totals['skipped_folders'] += job.get('skipped', 0)
            batches.append({'id': key, 'state': job['state'], 'started_at': job.get('started_at'),
                            'pdf_count': files, 'review_tables': review_tables, 'review_files': review_files,
                            'errors': len(job['errors']), **imports.get(key, {})})
        return {**result, **totals, 'batch_count': len(batches), 'recent': batches[:15]}

    def checkpoint(self, path, job):
        with self.lock:
            job['pause'] = load(path)['pause']
            write(path, job)
            return job['pause']

    def cache_file(self, key):
        return self.output / 'cache' / f'{key}.json'

    def read_document(self, path, checksum, profiles):
        cache = self.output / 'documents' / f'{PARSER_VERSION}-{checksum}.json'
        if cache.is_file():
            return load(cache)
        document = read_pdf(path, checksum, profiles)
        write(cache, document)
        return document

    def history(self, key):
        claimed, closed = defaultdict(set), set()
        bundles, spectra, reports = set(), set(), []
        for path in self.output.glob('batches/*/batch.json'):
            if path.parent.name == key:
                continue
            job = load(path)
            result = job.get('result_file')
            if not result or not Path(result).is_file():
                continue
            if not self.importer and job.get('molecules') and not Path(result).with_name('assignments.sqlite').is_file():
                continue
            for fingerprint, record in job['records'].items():
                claimed[fingerprint].update(record['rows'])
                if record['closed']:
                    closed.add(fingerprint)
                    if record.get('bundle'):
                        bundles.add(record['bundle'])
            if 'exported_spectra' in job:
                spectra.update(job['exported_spectra'])
            else:
                reports.append(Path(result))
        for path in self.output.glob('*/extraction_report.json'):
            reports.append(path)
        for path in reports:
            report = load(path)
            if not report.get('saved') or not path.with_name('assignments.sqlite').is_file():
                continue
            for table in report['tables']:
                if table.get('molecule_id'):
                    spectra.update(spectrum_keys(table).values())
            # Old reports do not inventory SDFs. Only assume a PDF-only input
            # when these really were the original inputs, or reviewed profiles were used.
            files = {f['name']: f['sha256'] for f in report['files']}
            used_sdf = any(table.get('metadata', {}).get('mapping_file') for table in report['tables'])
            if (path.parent.parent == self.output and not used_sdf and files
                    and all(t.get('molecule_id') and not t['issues'] for t in report['tables'])
                    and all(f['status'] not in ('读取失败', '待复核版式', '需 OCR') for f in report['files'])):
                bundles.add(source_bundle(files))
        return claimed, closed, bundles, spectra

    def work(self, key):
        path = self.state_file(key)
        job = load(path)
        root = Path(job['folder'])
        selected = set(job['molecules'])
        job['errors'] = []
        try:
            claimed, closed, bundles, known_spectra = self.history(key)
            if not self.importer:
                for fingerprint, record in job['records'].items():
                    tables = load(self.cache_file(fingerprint))['tables']
                    for index in record['rows']:
                        keys = spectrum_keys(tables[index])
                        nuclei = record.get('spectra', {}).get(str(index), keys)
                        known_spectra.update(keys[nucleus] for nucleus in nuclei)
            queue = folders(root)
            if not queue:
                raise ValueError('根目录及其直接子文件夹都没有 PDF。')
            job['total'] = len(queue)
            job['skipped'] = 0
            reason = '当前目录的论文已全部检查。'
            for i, folder in enumerate(queue, 1):
                if self.checkpoint(path, job):
                    break
                if job['limit'] is not None and len(selected) >= job['limit']:
                    reason = '已达到本批分子目标。还有未保存的分子时，可开始新批次。'
                    break
                with self.lock:
                    self.progress = {'index': i, 'current': str(folder.relative_to(root)), 'filename': None}
                job['index'] = i
                try:
                    fingerprint, checksums, stamp = snapshot(folder)
                    current = job['records'].get(fingerprint)
                    cache = self.cache_file(fingerprint)
                    bundle = source_bundle(checksums)
                    if fingerprint in closed or bundle in bundles or (current and current['closed']):
                        job['skipped'] += 1
                        continue
                    if not cache.exists():
                        def progress(index, total, filename):
                            with self.lock:
                                self.progress['filename'] = filename

                        report = collect(folder, progress, checksums=checksums, read_document=self.read_document, draw=False)
                        if file_stamp(source_paths(folder)) != stamp:
                            raise ValueError('处理期间文件有变更，请重试。')
                        write(cache, report)
                    else:
                        report = load(cache)
                except (ValueError, OSError, RuntimeError) as exc:
                    job['errors'].append(f'{folder.name} · {exc}')
                    continue
                if self.checkpoint(path, job):
                    break
                record = job['records'].setdefault(fingerprint, {'rows': [], 'closed': False})
                record['bundle'] = bundle
                record['pdf_count'] = len(report['files'])
                record['review_files'] = sum(f['status'] in ('读取失败', '待复核版式', '需 OCR')
                                             for f in report['files'])
                used = claimed[fingerprint] | set(record['rows'])
                groups = defaultdict(list)
                for index, table in enumerate(report['tables']):
                    if index not in used:
                        identity = table.get('smiles') if not table['issues'] else None
                        groups[identity].append(index)
                for identity, indexes in groups.items():
                    for index in indexes:
                        table = report['tables'][index]
                        signatures = spectrum_keys(table)
                        fresh = signatures if self.importer else {
                            nucleus: code for nucleus, code in signatures.items() if code not in known_spectra}
                        if fresh and identity not in selected and job['limit'] is not None and len(selected) >= job['limit']:
                            continue
                        if self.importer and signatures:
                            outcomes = self.importer.ingest(key, table)
                            record.setdefault('imports', {})[str(index)] = outcomes
                            fresh = {n: signatures[n] for n, result in outcomes.items() if result['status'] in CHANGES}
                            record['duplicates'] = record.get('duplicates', 0) + sum(
                                item['status'] == 'duplicate' for item in outcomes.values())
                        elif self.importer:
                            record['review_tables'] = record.get('review_tables', 0) + 1
                        elif signatures:
                            record.setdefault('spectra', {})[str(index)] = list(fresh)
                            record['duplicates'] = record.get('duplicates', 0) + len(signatures) - len(fresh)
                        if fresh:
                            selected.add(identity)
                            if not self.importer:
                                known_spectra.update(fresh.values())
                        record['rows'].append(index)
                record['closed'] = len(set(record['rows']) | used) == len(report['tables'])
                if record['closed']:
                    bundles.add(bundle)
                job['molecules'] = sorted(selected)
                job['duplicate_spectra'] = sum(r.get('duplicates', 0) for r in job['records'].values())
                self.checkpoint(path, job)
                if job['limit'] is not None and len(selected) >= job['limit']:
                    reason = '已达到本批分子目标。还有未保存的分子时，可开始新批次。'
                    break
            self.publish(key, job, root)
            terminal = 'paused' if job['pause'] else 'complete'
            if terminal == 'paused':
                reason = '已暂停。已保存的分子保留在本批结果中。'
        except Exception as exc:
            terminal, reason = 'error', str(exc)
        finally:
            with self.lock:
                if terminal == 'complete' and load(path)['pause']:
                    terminal, reason = 'paused', '已暂停。已保存的分子保留在本批结果中。'
                job['molecules'] = sorted(selected)
                job['duplicate_spectra'] = sum(r.get('duplicates', 0) for r in job['records'].values())
                write(path, {**job, 'state': terminal, 'reason': reason, 'pause': False, 'updated_at': db.now()})
                self.progress = {}
                self.active = None

    def publish(self, key, job, root):
        files, tables, portable, warnings = [], [], [], list(job['errors'])
        for fingerprint, record in job['records'].items():
            report = load(self.cache_file(fingerprint))
            folder = Path(report['folder']).relative_to(root)
            offset = len(files)
            for file in report['files']:
                relative = folder / file['name']
                files.append({**file, 'id': file['id'] + offset,
                              'name': str(relative), 'relative_path': str(relative)})
            for index in sorted(record['rows']):
                table = report['tables'][index]
                nuclei = record.get('spectra', {}).get(str(index))
                if nuclei is not None:
                    if not nuclei:
                        continue
                    for nucleus in table['peaks']:
                        if nucleus not in nuclei:
                            table['peaks'][nucleus] = []
                if table.get('mapping'):
                    table['mapping'] = {int(tag): atom for tag, atom in table['mapping'].items()}
                for field in ('file_id', 'supplement_file_id'):
                    if field in table:
                        table[field] += offset
                if self.importer:
                    imported = record.get('imports', {}).get(str(index), {})
                    table['main_import'] = imported
                    retained = {n for n, item in imported.items() if item['status'] in CHANGES}
                    if retained:
                        portable.append({**table, 'peaks': {n: rows if n in retained else []
                                                            for n, rows in table['peaks'].items()}})
                tables.append(table)
            warnings.extend(f'{folder} · {warning}' for warning in report['warnings'])
        if not files and not warnings:
            return
        dest = self.home(key) / f'result-{uuid.uuid4().hex[:10]}'
        save_results(dest, tables, {'folder': str(root), 'files': files, 'warnings': warnings,
                                    'pdf_count': len(files), 'integrated': bool(self.importer),
                                    'duplicate_spectra': job.get('duplicate_spectra', 0),
                                    'peak_count': sum(len(peaks) for t in tables for peaks in t['peaks'].values())},
                     portable_tables=portable if self.importer else None)
        if not self.importer:
            job['exported_spectra'] = sorted({code for table in tables for code in spectrum_keys(table).values()})
        job['result_file'] = str(dest / 'extraction_report.json')
