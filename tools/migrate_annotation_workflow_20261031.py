"""Apply the explicitly requested historical completeness policy online."""
import argparse
import json
from pathlib import Path
import sys
import time

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from annotation_tool_v3 import core as db,workflow_policy


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--db',type=Path,default=db.DEFAULT_DB)
    p.add_argument('--report',type=Path,required=True)
    args=p.parse_args()
    if args.report.exists():raise ValueError('Choose a new report file')
    started=time.monotonic()
    with db.connect(args.db,timeout=10) as c:
        c.execute('BEGIN IMMEDIATE')
        result=workflow_policy.backfill(c,catch_up=True)
        c.commit()
    result.update(database=str(db.main_database(args.db).resolve()),elapsed_seconds=round(time.monotonic()-started,2),
        before_values='workflow_completeness_history',note='Only current usable spectra of finished or historically human/AI-annotated molecules; no annotation labels changed.')
    args.report.parent.mkdir(parents=True,exist_ok=True)
    with args.report.open('x',encoding='utf-8') as handle:json.dump(result,handle,ensure_ascii=False,indent=2)
    print(json.dumps(result,ensure_ascii=False))


if __name__=='__main__':main()
