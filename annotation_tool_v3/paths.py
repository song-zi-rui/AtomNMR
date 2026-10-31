"""One location for the three live databases, separate from generated artifacts."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT / 'data/assignment_database_v1'
DATABASE_DIR = DATA_ROOT / 'databases'
DEFAULT_DB = DATABASE_DIR / 'assignments.sqlite'
LEGACY_DB = DATA_ROOT / 'assignments.sqlite'
DATABASE_NAMES = ('assignments.sqlite', 'ai_annotations.sqlite', 'ai_benchmarks.sqlite')


def main_database(path=DEFAULT_DB):
    """Keep old explicit local commands working, never create a second blank DB."""
    path = Path(path)
    if path.resolve() == LEGACY_DB.resolve() and DEFAULT_DB.is_file():
        if path.exists():
            raise ValueError('新旧数据库路径同时存在，请先核对；程序不会自行选择或覆盖。')
        return DEFAULT_DB
    return path


def artifacts_root(path):
    parent = main_database(path).parent
    return parent.parent if parent.name == 'databases' else parent


def companion_database(path, name):
    if name not in DATABASE_NAMES:
        raise ValueError('Unknown application database name')
    return main_database(path).parent / name
