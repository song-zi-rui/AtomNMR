"""Read-only production verification and a searchable affected-spectrum list."""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sqlite3
import sys
import urllib.request

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from annotation_tool_v3 import core
from annotation_tool_v3.exporters import training_view
from tools.shiftdb_atom_identity import POLICY


def verify(database, output, port=None):
    output.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(f"file:{database.resolve().as_posix()}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    conn.execute("BEGIN")
    counts, statuses = Counter(), Counter()
    query = """SELECT s.id,s.molecule_id,s.source_json,a.policy,a.validation_json
        FROM spectra s LEFT JOIN spectrum_assignment_audit a ON a.spectrum_id=s.id
        WHERE s.source_key>='shiftdb:' AND s.source_key<'shiftdb;' AND s.source='nmrshiftdb2'"""
    with (output / "production_affected_spectra.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["molecule_id", "spectrum_id", "shiftdb_id", "source_record", "field", "status", "issue_counts"])
        for row in conn.execute(query):
            counts['source_spectra'] += 1
            if row['policy'] != POLICY:
                counts['missing_current_audit'] += 1
                continue
            validation = json.loads(row['validation_json'])
            statuses[validation['status']] += 1
            if validation['status'] != 'identity_validated':
                counts['flagged_spectra'] += 1
                src = json.loads(row['source_json'])
                writer.writerow([row['molecule_id'],row['id'],src.get('nmrshiftdb2_id',src.get('source_id','')),
                    src.get('source_record_index'),src.get('spectrum_field',src.get('field','')),validation['status'],
                    json.dumps(validation.get('issue_counts',{}),ensure_ascii=False)])
                for annotation in conn.execute("SELECT origin,atoms_json FROM annotations WHERE spectrum_id=?", (row['id'],)):
                    if annotation['origin'] == 'nmrshiftdb2' and json.loads(annotation['atoms_json']):
                        counts['flagged_source_labels_still_assigned'] += 1
    mid = 'Mb21e911341dfcb8ffdbd7df1'
    record = core.record(conn, mid)
    carbon = next(s for s in record['spectra'] if s['id']=='S541e6fe7a968beef5264c6b1')
    proton = next(s for s in record['spectra'] if s['id']=='S16b143fb224172e8beb1901e')
    assert [p['assignment']['atom_indices'] for p in carbon['peaks']]==[[7],[1],[4],[6],[5]]
    assert all(p['assignment']['origin']=='manual' for p in carbon['peaks'])
    assert all(not p['assignment']['atom_indices'] for p in proton['peaks'])
    assert sum(len(p.get('source_assignment_fragments',[])) for p in proton['peaks'])==6
    view = training_view(record)
    assert not any(view['alignment_supervision']['proton']['assignment_present_by_peak'])
    assert all(view['alignment_supervision']['carbon']['assignment_valid_by_peak'])
    forbidden = {'assignment','source_assignment_fragments','raw_text','atom_indices','source_atom_id'}
    assert all(not forbidden.intersection(p) for name in ('h_nmr_peaks','c_nmr_peaks') for p in view[name])
    (output/'specific_record_after_repair.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
    counts['missing_current_audit'] += 0
    counts['flagged_source_labels_still_assigned'] += 0
    assert counts['missing_current_audit']==0
    assert counts['flagged_source_labels_still_assigned']==0
    result = {'policy':POLICY,'database':str(database.resolve()),'counts':dict(counts),
              'source_spectrum_status_counts':dict(statuses),'specific_manual_carbon_labels_preserved':5,
              'specific_hydrogen_occurrences_preserved':6,'model_peak_input_assignment_leak_check':'passed'}
    conn.close()
    if port:
        with urllib.request.urlopen(f'http://127.0.0.1:{port}/api/molecule/{mid}',timeout=15) as response:
            live = json.load(response)
        assert live['revision']==record['revision']
        assert live['spectra']==record['spectra']
        result['live_server_port_verified']=port
    (output/'production_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db',type=Path,default=core.DEFAULT_DB)
    parser.add_argument('--output-dir',type=Path,default=Path('data/shiftdb_atom_identity_audit/full_snapshot_final'))
    parser.add_argument('--port',type=int)
    args=parser.parse_args()
    verify(args.db,args.output_dir,args.port)
