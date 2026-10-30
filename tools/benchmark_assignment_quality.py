"""Compare identical QC policies on disposable subsets; never modify production."""
import json
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from annotation_tool_v3 import core as db, migrate
from annotation_tool_v3.quality_job import recover_batched


def main():
    source = db.connect()
    rows = source.execute("SELECT * FROM candidates WHERE state='awaiting_source' AND molecule_id IS NULL ORDER BY smiles LIMIT 500").fetchall()
    wanted = sorted({s for r in rows for s in (r["smiles"],r["raw_smiles"]) if s})
    sources = source.execute("SELECT * FROM source_rows WHERE raw_smiles IN (" + ",".join("?" for _ in wanted) + ")",wanted).fetchall()
    source.close()
    timings, outputs = {}, []
    with tempfile.TemporaryDirectory(prefix="nmr_quality_benchmark_") as directory:
        for name, runner in (("per_candidate_transactions",migrate.recover),("prefetch_small_batches",recover_batched)):
            c = db.connect(Path(directory)/(name+".sqlite"))
            db.initialize(c)
            c.executemany("INSERT INTO candidates VALUES("+",".join("?" for _ in rows[0])+")",[tuple(r) for r in rows])
            c.executemany("INSERT INTO source_rows VALUES(?,?,?,?,?,?,?,?,?)",[tuple(r) for r in sources])
            c.commit()
            start=time.perf_counter()
            runner(c,limit=len(rows))
            timings[name]=round(time.perf_counter()-start,3)
            outputs.append([tuple(r) for r in c.execute("SELECT smiles,state,reason,molecule_id,source_pair_json FROM candidates ORDER BY smiles")])
            c.close()
    result={"candidates":len(rows),"source_rows":len(sources),"seconds":timings,
            "same_candidate_decisions":outputs[0]==outputs[1],
            "speedup":round(timings["per_candidate_transactions"]/timings["prefetch_small_batches"],2),
            "caveat":"Disposable small subset, both use the same quality checks; not a full-database speedup guarantee."}
    if not result["same_candidate_decisions"]:
        raise RuntimeError("Batching changed matching results")
    output=db.DEFAULT_DB.parent/"quality_batch_benchmark.json"
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
