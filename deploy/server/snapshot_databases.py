"""Stage AtomNMR SQLite files, without saved API keys, and a SHA256 list."""
import argparse
import hashlib
import shutil
import sqlite3
from contextlib import ExitStack, closing
from pathlib import Path
from tempfile import TemporaryDirectory


DATABASES = ('assignments.sqlite', 'ai_annotations.sqlite', 'ai_benchmarks.sqlite')
ANNOTATION_DATABASES = DATABASES[:2]


def sha256(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def strip_credentials(db):
    tables = {row[0] for row in db.execute("SELECT name FROM sqlite_schema WHERE type='table'")}
    if 'api_credentials' in tables:
        db.execute('DELETE FROM api_credentials')
    if 'settings' in tables:
        db.execute("DELETE FROM settings WHERE key IN ('active_credential','deepseek_api_key_dpapi')")
    db.commit()
    # Rebuild the file so deleted records cannot survive in unused pages.
    db.execute('VACUUM')


def snapshot(source, output, *, annotation_only=False):
    source, output = source.resolve(), output.resolve()
    names = ANNOTATION_DATABASES if annotation_only else DATABASES
    files = [source / name for name in names]
    if not all(path.is_file() for path in files):
        raise FileNotFoundError('Required AtomNMR databases: ' + ', '.join(names))
    if output.exists():
        raise FileExistsError('Choose a new destination directory; snapshots are never overwritten')
    output.parent.mkdir(parents=True, exist_ok=True)
    needed = sum(path.stat().st_size + sum(p.stat().st_size for p in path.parent.glob(path.name + '-wal')) for path in files)
    if shutil.disk_usage(output.parent).free < needed:
        raise OSError(f'Insufficient space for the databases and any WAL data ({needed:,} bytes)')
    output.mkdir()
    sums = []
    with ExitStack() as stack:
        private_dir = Path(stack.enter_context(TemporaryDirectory(prefix='atomnmr-snapshot-')))
        readers = [stack.enter_context(closing(sqlite3.connect(path.as_uri() + '?mode=ro', uri=True, timeout=20)))
                   for path in files]
        versions = [db.execute('PRAGMA data_version').fetchone()[0] for db in readers]
        for path, reader in zip(files, readers):
            target = output / path.name
            contains_keys = path.name == 'ai_annotations.sqlite'
            staging = private_dir / path.name if contains_keys else target
            print(f'Snapshotting {path.name} ({path.stat().st_size:,} bytes)', flush=True)
            with closing(sqlite3.connect(staging)) as writer:
                reader.backup(writer, pages=2048)
                # A transferable snapshot uses only its main file.
                writer.execute('PRAGMA journal_mode=DELETE')
                if contains_keys:
                    strip_credentials(writer)
                check = writer.execute('PRAGMA quick_check(1)').fetchone()[0]
                if check != 'ok':raise ValueError(f'{path.name}: SQLite check failed: {check}')
            if contains_keys:
                shutil.copyfile(staging, target)
            sums.append(f'{sha256(target)}  {path.name}\n')
        if versions != [db.execute('PRAGMA data_version').fetchone()[0] for db in readers]:
            raise RuntimeError('Source databases changed during backup. Stop the writers and use a new destination.')
    (output / 'SHA256SUMS').write_text(''.join(sums), encoding='ascii', newline='\n')
    print(f'Snapshot ready: {output}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, epilog='Stop all AtomNMR writers first to keep the three databases aligned.')
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--annotation-only', action='store_true', help='Save the main database and the automatic annotation history')
    args = parser.parse_args()
    snapshot(args.source, args.output, annotation_only=args.annotation_only)


if __name__ == '__main__':main()
