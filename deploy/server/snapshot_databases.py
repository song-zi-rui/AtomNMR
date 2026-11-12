"""Stage the three AtomNMR SQLite files and a SHA256 list for transfer or backup."""
import argparse
import hashlib
import shutil
import sqlite3
from contextlib import ExitStack, closing
from pathlib import Path


DATABASES = ('assignments.sqlite', 'ai_annotations.sqlite', 'ai_benchmarks.sqlite')


def sha256(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def snapshot(source, output):
    source, output = source.resolve(), output.resolve()
    files = [source / name for name in DATABASES]
    if not all(path.is_file() for path in files):
        raise FileNotFoundError('Source must contain all three AtomNMR databases')
    if output.exists():
        raise FileExistsError('Choose a new destination directory; snapshots are never overwritten')
    output.parent.mkdir(parents=True, exist_ok=True)
    needed = sum(path.stat().st_size + sum(p.stat().st_size for p in path.parent.glob(path.name + '-wal')) for path in files)
    if shutil.disk_usage(output.parent).free < needed:
        raise OSError(f'Insufficient space for the databases and any WAL data ({needed:,} bytes)')
    output.mkdir()
    sums = []
    with ExitStack() as stack:
        readers = [stack.enter_context(closing(sqlite3.connect(path.as_uri() + '?mode=ro', uri=True, timeout=20)))
                   for path in files]
        versions = [db.execute('PRAGMA data_version').fetchone()[0] for db in readers]
        for path, reader in zip(files, readers):
            target = output / path.name
            print(f'Snapshotting {path.name} ({path.stat().st_size:,} bytes)', flush=True)
            with closing(sqlite3.connect(target)) as writer:
                reader.backup(writer, pages=2048)
                # A transferable snapshot uses only its main file.
                writer.execute('PRAGMA journal_mode=DELETE')
                check = writer.execute('PRAGMA quick_check(1)').fetchone()[0]
                if check != 'ok':raise ValueError(f'{path.name}: SQLite check failed: {check}')
            sums.append(f'{sha256(target)}  {path.name}\n')
        if versions != [db.execute('PRAGMA data_version').fetchone()[0] for db in readers]:
            raise RuntimeError('Source databases changed during backup. Stop the writers and use a new destination.')
    (output / 'SHA256SUMS').write_text(''.join(sums), encoding='ascii', newline='\n')
    print(f'Snapshot ready: {output}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, epilog='Stop all AtomNMR writers first to keep the three databases aligned.')
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    snapshot(args.source, args.output)


if __name__ == '__main__':main()
