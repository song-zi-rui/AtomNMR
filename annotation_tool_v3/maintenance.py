"""Consistent backup and descriptive reports; no edits to annotations or sources."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sqlite3
import platform
from rdkit import rdBase
from . import core as db
from .paths import artifacts_root, main_database


def read_connection(path):
    path=main_database(path)
    c=sqlite3.connect(f"file:{Path(path).as_posix()}?mode=ro",uri=True,timeout=60)
    c.row_factory=sqlite3.Row
    return c


def backup(path):
    path=main_database(path)
    folder=artifacts_root(path)/"backups"
    folder.mkdir(exist_ok=True)
    dest=folder/("assignments_"+db.now().replace(":","").replace("+","_")+".sqlite")
    if dest.exists():
        raise ValueError("Backup already exists; never overwrite")
    source=read_connection(path)
    target=sqlite3.connect(dest)
    try:
        source.backup(target,pages=2048)
    finally:
        target.close()
        source.close()
    print(dest)


def report(path):
    c=read_connection(path)
    c.execute("PRAGMA cache_size=-131072")
    c.execute("BEGIN")
    print("Report: counts and input metadata",flush=True)
    result=db.stats(c)
    result["created_at"]=db.now()
    result["schema"]=db.SCHEMA
    result["runtime"]={"python":platform.python_version(),"rdkit":rdBase.rdkitVersion}
    result["input_metadata"]={r[0]:json.loads(r[1]) for r in c.execute("SELECT key,value FROM meta")}
    def grouped(sql):
        return [dict(r) for r in c.execute(sql)]
    result["spectra_by_source_nucleus"]=grouped("SELECT source,nucleus,count(*) count,sum(json_array_length(peaks_json)) signals,sum(solvent IS NOT NULL) solvent_known,sum(temperature_k IS NOT NULL) temperature_reported,sum(frequency_mhz IS NOT NULL) frequency_reported FROM spectra GROUP BY source,nucleus")
    result["molecules_by_source"]=grouped("SELECT origin,count(*) count FROM molecules GROUP BY origin")
    result["pairs"]=grouped("SELECT basis,solvent_match,count(*) count FROM pairs GROUP BY basis,solvent_match")
    result["annotations"]=grouped("SELECT origin,status,count(*) count FROM annotations GROUP BY origin,status")
    print("Report: source recovery and labels",flush=True)
    result["source_recovery_reasons"]=grouped("SELECT state,reason,count(*) count FROM candidates NOT INDEXED GROUP BY state,reason")
    result["legacy_files"]=grouped("SELECT state,count(*) count FROM legacy_annotations GROUP BY state")
    result["import_issue_reasons"]=grouped("SELECT stage,reason,count(*) count FROM import_issues GROUP BY stage,reason")
    result["mw_distribution"]=grouped("SELECT origin,CASE WHEN mw<150 THEN '<150' WHEN mw<250 THEN '150-249' WHEN mw<350 THEN '250-349' WHEN mw<500 THEN '350-499' ELSE '>=500' END bin,count(*) count FROM molecules GROUP BY origin,bin")
    result["ring_distribution"]=grouped("SELECT origin,json_extract(structure_json,'$.ring_count') rings,count(*) count FROM molecules GROUP BY origin,rings")
    result["solvent_distribution"]=grouped("SELECT source,nucleus,coalesce(solvent,'unknown') solvent,count(*) count FROM spectra GROUP BY source,nucleus,solvent ORDER BY count DESC")
    print("Report: distributions",flush=True)
    element_counts=Counter()
    for row in c.execute("SELECT origin,structure_json FROM molecules"):
        for element in {a["symbol"] for a in json.loads(row[1])["atoms"]}:
            element_counts[(row[0],element)]+=1
    result["element_distribution"]=[{"origin":k[0],"element":k[1],"molecules":v} for k,v in sorted(element_counts.items())]
    print("Report: foreign key check",flush=True)
    result["foreign_key_errors"]=[tuple(r) for r in c.execute("PRAGMA foreign_key_check")]
    result["invalid_empty_assignments"]=c.execute("SELECT count(*) FROM annotations WHERE status IN ('assigned','equivalent_group','overlap','source_group') AND json_array_length(atoms_json)=0").fetchone()[0]
    c.close()
    c=read_connection(path)
    print("Report: SQLite integrity scan",flush=True)
    result["sqlite_quick_check"]=[r[0] for r in c.execute("PRAGMA quick_check")]
    c.close()
    print("Report: source hash",flush=True)
    source=Path(result["input_metadata"]["shiftdb_path"])
    digest=hashlib.sha256()
    with source.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            digest.update(block)
    result["filtered_sdf_hash_unchanged"]=digest.hexdigest()==result["input_metadata"]["shiftdb_sha256"]
    folder=artifacts_root(path)
    (folder/"database_report.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    labels={"complete_remaining":"补充部分归属","both_unassigned":"优先：H/C 均未归属","waiting_spectrum":"缺少一种实验谱","complete":"已存谱均有归属"}
    lines=["# 实验 NMR 归属数据库迁移报告","",f"生成时间：{result['created_at']}。这是主库快照，不是模型性能或人工标签正确率评估。","",f"主库 **{result['molecules']:,}** 个分子，**{result['spectra']:,}** 条实验谱。","","## 工作队列","","| 队列 | 分子数 |","|---|---:|"]
    lines += [f"| {labels.get(k,k)} | {v:,} |" for k,v in result["queues"].items()]
    lines += ["","## 来源与字段覆盖","","| 来源 | 核种 | 实验谱 | 信号数 | 已报溶剂 | 已报温度 | 已报频率 |","|---|---|---:|---:|---:|---:|---:|"]
    lines += [f"| {r['source']} | {r['nucleus']} | {r['count']:,} | {r['signals']:,} | {r['solvent_known']:,} | {r['temperature_reported']:,} | {r['frequency_reported']:,} |" for r in result["spectra_by_source_nucleus"]]
    lines += ["","信号数不是原子数；同一分子的多次实验各自计数。ShiftDB 源信号组不是自动确证的等价组。频率保留来源口径，未知不补默认值。","","## NMRPeak 候选目录的去向","","| 状态 | 数量 |","|---|---:|"]
    lines += [f"| {k} | {v:,} |" for k,v in result["candidates"].items()]
    lines += ["","`imported` 是新增结构；`complemented` 是已有结构补实验谱；`duplicate_complete` 已由实验 ShiftDB 覆盖；`overlap_not_complementable` 未找到需要且溶剂相容的补谱位置；`source_unresolved` 未恢复唯一原表对，未捏造进入主库。","","## 旧标注文件","","| 迁移状态 | 文件数 |","|---|---:|"]
    lines += [f"| {r['state']} | {r['count']} |" for r in result["legacy_files"]]
    lines += ["","`reference_only_not_human_labels` 是原文件明确带 GPT / 伪标签标记；保留副本但不纳入人工归属。旧标注文件从未改写。","","## 分布概况","","| 来源 | 分子量区间 | 分子数 |","|---|---|---:|"]
    lines += [f"| {r['origin']} | {r['bin']} | {r['count']:,} |" for r in result["mw_distribution"]]
    lines += ["","环数、元素和完整溶剂分布见同目录 JSON。分布仅描述本地来源，不能据此宣称代表真实世界。回源要求双谱和已报告溶剂会引入额外选择偏差。","","## 完整性检查","",f"- 外键错误：{len(result['foreign_key_errors'])}",f"- SQLite quick_check：{', '.join(result['sqlite_quick_check'])}",f"- 导入前后实验 SDF SHA256 相同：{result['filtered_sdf_hash_unchanged']}",f"- 导入异常记录：{result['issues']}","- 原始 SDF、NMRexp LMDB / Parquet 和旧标注只读；所有新标注、迁移目录与导出存入新路径。","- 全量恢复完毕与否请同时查看 JSON 中 progress，以及 candidates 中是否仍有 indexed / awaiting_source。","- 取消复核和预测模块，不意味着自动证明归属正确；歧义和源组通过标签与掩码明确区分。",""]
    (folder/"database_report.md").write_text("\n".join(lines),encoding="utf-8")
    print(json.dumps({k:result[k] for k in ("molecules","spectra","queues","candidates","legacy_files","foreign_key_errors","sqlite_quick_check","filtered_sdf_hash_unchanged")},ensure_ascii=False,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("action",choices=["report","backup"])
    p.add_argument("--db",type=Path,default=db.DEFAULT_DB)
    a=p.parse_args()
    (report if a.action=="report" else backup)(a.db)


if __name__=="__main__":
    main()
