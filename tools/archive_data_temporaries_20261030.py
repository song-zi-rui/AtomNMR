"""Reversible, explicit-list cleanup. Dry run by default; never delete data."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TARGET=ROOT/'trashbin/data_cleanup_20261030'


def plan():
    names={
        'data/blind_close_carbon_examples':'offline demonstration output; not used by annotation servers',
        'data/blind_complex_carbon_examples':'offline demonstration output; not used by annotation servers',
        'data/gpt_assignment_review':'old offline AI review examples',
        'data/assignment_database_v1/audits':'reproducible offline diagnostic reports; live filters use database records',
        'data/assignment_database_v1/diagnostics':'standalone diagnostic database and reports',
        'data/assignment_database_v1/backup_control':'removed OpenList state; Git mirror already moved outside data',
        'data/assignment_database_v1/benchmarks/production_timing_20261028.json':'old timing report; live timing reads SQLite',
        'data/assignment_database_v1/benchmarks/school_api_diagnostic_20261028.json':'old API diagnostic report',
        'data/assignment_database_v1/quality_batch_benchmark.json':'old performance benchmark',
        'data/assignment_database_v1/save_performance_fix_20260929.json':'old save-performance test report',
        'data/assignment_database_v1/quality_job.lock':'completed 2026-09-29 maintenance lock',
    }
    for name in ['nmrexp_override_smoke_1000','nmrexp_smoke_1000','nmrshiftdb2_smoke_100']:
        names['data/c13_decoupling_audit/'+name]='small smoke-test output; full audit results retained'
    for name in ['full_snapshot','full_snapshot_v2','full_snapshot_v3']:
        names['data/shiftdb_atom_identity_audit/'+name]='superseded audit iterations; final snapshot and runtime source index retained'
    for name in ['2026-09-29T050359_0000_44d9c3a9','2026-09-29T051209_0000_6c93bd8b','2026-09-29T051210_0000_49fa646d','2026-10-30T030123_0000_8cbafb7b']:
        names['data/assignment_database_v1/exports/'+name]='initial ethanol single-record format smoke tests or zero-record export check'
    for file in (ROOT/'data/assignment_database_v1').glob('*.log'):names[file.relative_to(ROOT).as_posix()]='closed historical service, migration or UI fixture log'
    logs=ROOT/'data/assignment_database_v1/service_logs'
    outputs=sorted(logs.glob('*.out.log'))
    active_prefix=outputs[-1].name.removesuffix('.out.log') if outputs else None
    for file in logs.glob('*.log'):
        if active_prefix and file.name.startswith(active_prefix+'.'):continue
        names[file.relative_to(ROOT).as_posix()]='closed historical service log; current log pair retained'
    result=[]
    for name,reason in names.items():
        source=(ROOT/name).resolve();dest=(TARGET/name).resolve()
        if not source.exists():continue
        if not source.is_relative_to(ROOT/'data') or not dest.is_relative_to(TARGET):raise ValueError('Unsafe resolved path')
        if source.name in ('assignments.sqlite','ai_annotations.sqlite','ai_benchmarks.sqlite','source_seek_index.sqlite'):raise ValueError('Protected runtime file')
        files=[p for p in source.rglob('*') if p.is_file()] if source.is_dir() else [source]
        if any(p.is_symlink() for p in files):raise ValueError('Symlink in archive source')
        result.append({'source':source.relative_to(ROOT).as_posix(),'destination':dest.relative_to(ROOT).as_posix(),
                       'files':len(files),'bytes':sum(p.stat().st_size for p in files),'reason':reason})
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--apply',action='store_true');args=p.parse_args()
    rows=plan();TARGET.mkdir(parents=True,exist_ok=True)
    manifest={'created_at':datetime.now(timezone.utc).isoformat(),'mode':'move, not delete','items':rows,
              'restore':'With annotation service stopped, move each destination back to its original source; do not overwrite existing files.'}
    for row in rows:
        row['state']='planned'
        if args.apply:
            source=(ROOT/row['source']).resolve();dest=(ROOT/row['destination']).resolve()
            if dest.exists():raise ValueError('Archive target already exists: '+row['destination'])
            dest.parent.mkdir(parents=True,exist_ok=True)
            try:source.rename(dest);row['state']='moved'
            except OSError as exc:row.update(state='not_moved',error=type(exc).__name__)
    manifest_path=TARGET/('manifest.json' if args.apply else 'plan.json')
    manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'manifest':str(manifest_path),'items':len(rows),'files':sum(r['files'] for r in rows),'bytes':sum(r['bytes'] for r in rows),
                      'moved':sum(r['state']=='moved' for r in rows),'not_moved':[r['source'] for r in rows if r['state']=='not_moved']}))


if __name__=='__main__':main()
