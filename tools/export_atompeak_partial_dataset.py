"""Export current experimental annotations into a new, split AtomPeak dataset."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from atompeak_data.export import export_training_dataset
from annotation_tool_v3.paths import DEFAULT_DB


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--db',type=Path,default=DEFAULT_DB)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--valid-ratio',type=float,default=.1);p.add_argument('--test-ratio',type=float,default=.1)
    p.add_argument('--max-seq-len',type=int,default=1024);p.add_argument('--max-target-len',type=int,default=512)
    p.add_argument('--use-solvent',action='store_true');p.add_argument('--use-temperature',action='store_true');p.add_argument('--use-field-strength',action='store_true')
    p.add_argument('--annotated-only',action='store_true',help='Optional: omit molecules without any saved definite peak assignments')
    p.add_argument('--molecule-ids',help='Optional newline-separated molecule IDs, e.g. a small export validation selection')
    p.add_argument('--archive',action='store_true')
    p.add_argument('--estimate-only',action='store_true',help='Check/serialize candidates in memory only; do not write training data')
    p.add_argument('--estimate-report',type=Path,help='Optional JSON report for --estimate-only')
    args=p.parse_args()
    ids=Path(args.molecule_ids).read_text(encoding='utf-8-sig').split() if args.molecule_ids else None
    filters={}
    selection_started=time.monotonic()
    if args.annotated_only:
        import sqlite3
        c=sqlite3.connect(f'file:{args.db.resolve().as_posix()}?mode=ro',uri=True)
        # These are the only active queues that can contain saved assignments.
        # Use their existing index instead of scanning ~600k large molecule rows.
        marked={r[0] for r in c.execute("SELECT id FROM molecules WHERE queue IN ('complete','complete_remaining','waiting_spectrum') AND c_assigned+h_assigned>0")};c.close()
        ids=sorted(marked&set(ids) if ids is not None else marked)
    selection_seconds=time.monotonic()-selection_started
    if args.estimate_only:
        from atompeak_data.estimate import estimate_training_dataset
        result=estimate_training_dataset(args.db,molecule_ids=ids,seed=args.seed,
            valid_ratio=args.valid_ratio,test_ratio=args.test_ratio,max_seq_len=args.max_seq_len,max_target_len=args.max_target_len,
            use_solvent=args.use_solvent,use_temperature=args.use_temperature,use_field_strength=args.use_field_strength,
            require_supervision=args.annotated_only,progress=lambda v:print(json.dumps(v,ensure_ascii=False),flush=True))
        result['candidate_selection_seconds']=round(selection_seconds,2)
        result['estimated_export_seconds_range']=[round(s+selection_seconds) for s in result['estimated_export_seconds_range']]
        if args.estimate_report:
            args.estimate_report.parent.mkdir(parents=True,exist_ok=True)
            args.estimate_report.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
        print(json.dumps(result,ensure_ascii=False,indent=2));return
    result=export_training_dataset(args.db,args.output_dir,molecule_ids=ids,filters=filters,seed=args.seed,
        valid_ratio=args.valid_ratio,test_ratio=args.test_ratio,max_seq_len=args.max_seq_len,max_target_len=args.max_target_len,
        archive=args.archive,use_solvent=args.use_solvent,use_temperature=args.use_temperature,use_field_strength=args.use_field_strength,
        require_supervision=args.annotated_only,
        progress=lambda v:print(json.dumps(v,ensure_ascii=False),flush=True))
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
