#!/usr/bin/env python
"""Repair structure/scaffold split leakage in an existing assignment bank.

This standalone compatibility entry point also works with the first version of
``nmrpeak.shift_assignment.database`` that predates the ``repair-splits`` CLI
subcommand.  It avoids reparsing the source SDF.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
os.environ["NMRPEAK_LIGHT_IMPORT"] = "1"

from nmrpeak.shift_assignment.chemistry import stable_split  # noqa: E402
from nmrpeak.shift_assignment.database import (  # noqa: E402
    connect_database,
    read_metadata,
    rebuild_environment_stats,
)


def repair(database_path: Path) -> dict:
    conn = connect_database(database_path)
    try:
        metadata = read_metadata(conn)
        seed = int(metadata.get("seed", 20260828))
        train_ratio = float(metadata.get("train_ratio", 0.8))
        valid_ratio = float(metadata.get("valid_ratio", 0.1))
        rows = conn.execute(
            "SELECT id,structure_key,split_group,split FROM source_spectra ORDER BY id"
        ).fetchall()
        parent: dict[str, str] = {}

        def find(value: str) -> str:
            parent.setdefault(value, value)
            while parent[value] != value:
                parent[value] = parent[parent[value]]
                value = parent[value]
            return value

        def union(left: str, right: str) -> None:
            left_root, right_root = find(left), find(right)
            if left_root == right_root:
                return
            if left_root < right_root:
                parent[right_root] = left_root
            else:
                parent[left_root] = right_root

        for row in rows:
            union("structure:" + row["structure_key"], "group:" + row["split_group"])
        component_members: dict[str, list[str]] = {}
        for node in parent:
            component_members.setdefault(find(node), []).append(node)
        component_split = {
            root: stable_split(min(members), seed, train_ratio, valid_ratio)
            for root, members in component_members.items()
        }
        updates = []
        moved = 0
        for row in rows:
            target = component_split[find("structure:" + row["structure_key"])]
            moved += int(target != row["split"])
            updates.append((target, int(row["id"])))
        conn.executemany("UPDATE source_spectra SET split=? WHERE id=?", updates)
        conn.execute(
            """UPDATE observations
               SET split=(SELECT s.split FROM source_spectra s
                          WHERE s.id=observations.source_spectrum_id)"""
        )
        conn.commit()
        print("[repair-splits] rebuilding train-only environment statistics", flush=True)
        rebuild_environment_stats(conn)
        structure_cross = int(
            conn.execute(
                """SELECT COUNT(*) FROM (
                       SELECT structure_key FROM source_spectra
                       GROUP BY structure_key HAVING COUNT(DISTINCT split)>1)"""
            ).fetchone()[0]
        )
        group_cross = int(
            conn.execute(
                """SELECT COUNT(*) FROM (
                       SELECT split_group FROM source_spectra
                       GROUP BY split_group HAVING COUNT(DISTINCT split)>1)"""
            ).fetchone()[0]
        )
        report = {
            "schema_version": "nmr-assignment-split-repair-v1",
            "database": str(database_path.resolve()),
            "repaired_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "component_count": len(component_members),
            "spectra_moved": moved,
            "structure_cross_split_after": structure_cross,
            "split_group_cross_split_after": group_cross,
            "spectra_by_split": {
                row["split"]: int(row["n"])
                for row in conn.execute("SELECT split,COUNT(*) n FROM source_spectra GROUP BY split")
            },
            "observations_by_split": {
                row["split"]: int(row["n"])
                for row in conn.execute("SELECT split,COUNT(*) n FROM observations GROUP BY split")
            },
        }
        conn.commit()
        return report
    finally:
        conn.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--database",
        type=Path,
        default=ROOT / "data/shift_assignment_assistant/nmrshift_environment.sqlite",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "data/shift_assignment_assistant/split_repair_report.json",
    )
    args = parser.parse_args()
    report = repair(args.database.resolve())
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if report["structure_cross_split_after"] or report["split_group_cross_split_after"]:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

