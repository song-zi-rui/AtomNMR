"""Finish the full migration, then publish final checks and descriptive reports."""
import sys
import os
import json
import hashlib
import threading
from . import core as db
from . import migrate
from .maintenance import report
from .generate_figures import main as figures
from .paths import artifacts_root


def main():
    # The checkpointed CLI owns all data writes; source files remain read-only.
    state={"pid":os.getpid(),"started_at":db.now(),"state":"running"}
    state_file=artifacts_root(db.DEFAULT_DB)/"migration_job.json"
    state_lock=threading.Lock()
    stopped=threading.Event()
    def save_state():
        with state_lock:
            state["updated_at"]=db.now()
            temporary=state_file.with_suffix(".json.tmp")
            temporary.write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding="utf-8")
            temporary.replace(state_file)
    def heartbeat():
        while not stopped.wait(5):
            try:
                save_state()
            except OSError:
                # The UI marks a stale heartbeat as unconfirmed, not running.
                pass
    worker=threading.Thread(target=heartbeat,daemon=True)
    try:
        save_state()
        worker.start()
        c=db.connect()
        db.initialize(c)
        db.set_meta(c,"implementation_sha256",{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (db.ROOT/"annotation_tool_v3").glob("*.py")})
        migrate.migrate_legacy(c)
        c.commit()
        c.close()
        sys.argv=[sys.argv[0],"all"]
        migrate.main()
        state["state"]="generating_report"
        save_state()
        report(db.DEFAULT_DB)
        figures()
        state.update(state="finished",finished_at=db.now())
        save_state()
        print("FINISHED: migration, reports and integrity checks written.",flush=True)
    except Exception as exc:
        state.update(state="failed",error=str(exc),failed_at=db.now())
        save_state()
        raise
    finally:
        stopped.set()
        if worker.is_alive():
            worker.join(timeout=2)


if __name__=="__main__":
    main()
