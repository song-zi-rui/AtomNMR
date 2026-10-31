# 实验 NMR 归属数据库 v1

新的 SQLite 标注系统，默认端口 8767；不改动旧 v2 系统或源文件。

## 开始使用

在项目根目录：

```powershell
python -m annotation_tool_v3.server --port 8767 --open
```

打开 http://127.0.0.1:8767/ 。也可运行本目录 `start.ps1`。

- `data/assignment_database_v1/databases/assignments.sqlite`：日常主库，同目录另有 `ai_annotations.sqlite` 和 `ai_benchmarks.sqlite`。
- [AtomNMR 标注系统完整上手教程](../docs/AtomNMR_标注系统完整上手教程.md)：首次使用、原子与整组选取、状态判断、保存快捷键、数据管理、导出和常见问题。
- `docs/实验NMR数据库建设与使用_v2_20260929.md`：建设规则、字段字典、操作说明与限制。
- `/overview.html`：三张 SVG 与真实分子示例。
- `/auto.html`：持续自动标注、暂停/恢复、正式密钥管理及直接入库；[使用指南](../docs/持续自动标注与密钥管理_20261027.md)。AI 标签记录具体模型，只补空缺；原人工与原库标签保留。`/ai.html?history=1` 查看旧批次。
- 界面支持 JSONL / Parquet / LMDB 导出，范围为当前分子、筛选结果或全部已迁入主库。

## 重新迁移 / 中断后继续

```powershell
python -u -m annotation_tool_v3.quality_job full
```

当前库已完成源文件和候选索引。`quality_job full` 会跳过已检查的候选，继续质量检查、重新补谱和剩余回源，最终执行完整性检查。迁移到全新空库时才需要先运行原 `migrate` 的 `shiftdb`、`candidates`、`sources` 索引步骤。

旧 `finish_migration` 是此前的迁移入口；本次修复后使用上面的 `quality_job full`。不要同时开两个迁移进程。新版结果见 `quality_full_report.json` 和 `quality_final_database_report.json`，只有通过新检查的 NMRexp 谱才能导出。

后台任务状态在 `data/assignment_database_v1/migration_job.json`，包含进程号、阶段、阶段 ETA、心跳和运行状态。本次活动日志是 `quality_full_20260929.log` 与 `quality_full_20260929_error.log`；旧日志保留，不覆盖。电脑关机后重新运行续跑命令，不会自动向系统注册开机任务。

输入位置固定于 `migrate.py`，目前使用实验筛选后的 ShiftDB SDF、已有 NMRexp LMDB、已有结构索引、四个原始 Parquet 分片。依赖 `rdkit`、`numpy`、`scipy`、`pyarrow`、`lmdb`；本机使用 `C:\Python314\python.exe`，未另装服务端数据库。

## 备份与自检

```powershell
python -m annotation_tool_v3.maintenance backup
python -m annotation_tool_v3.maintenance report
python -m unittest tests.test_assignment_database_v1 -v
```

备份使用 SQLite backup API。运行中请不要单独复制主文件而遗漏 WAL。每日操作只有归属、备注、保存、结束任务，没有复核/审批/预测/DFT 步骤。

## 2026-09-29 保存性能修复

- 启动时自动补建 `pairs(molecule_id,id)` 和列表排序索引，不改动归属或原始谱。
- 保存只读取一次完整记录，并在同一事务内生成确认结果；返回错误时不会出现“实际已提交却报告保存失败”的数据库读取路径。
- 保存遇到写锁时快速返回可重试状态，页面显示等待原因并有限次重试；网络断开、版本冲突不会盲目重放。失败时请保留页面，不要刷新。
- 自动保存、手动保存、结束任务共用顺序保存流程；请求期间的新修改不会被旧响应覆盖，切换分子和导出会等待新修改保存。
- NMRexp 回源的峰匹配在写事务之外执行，每个候选独立提交；锁冲突回滚后重试，失败候选不留半条数据。
- 后台状态根据任务状态文件及心跳显示。历史 `failed` 不会因“还有候选”而显示成运行中。继续时运行上面的新版 `quality_job full` 命令。

回归测试：

```powershell
python -m unittest tests.test_assignment_database_v1 tests.test_assignment_save_performance -v
node --test tests/test_assignment_save_ui.cjs
```

本次完整备份位于 `data/assignment_database_v1/backups/assignments_2026-09-29T062141_0000.sqlite`，性能及标签一致性验证见 `data/assignment_database_v1/save_performance_fix_20260929.json`。

## 数据含义与限制

1. 一个分子拥有多条实验谱，归属绑定谱和信号。
2. H 标签指向带 H 的重原子，不声称单个 H 的立体专属归属。
3. 源多原子标签保持 `source_group`，默认不当作无歧义训练标签。
4. 溶剂相同不等于所有条件相同。原表定位不唯一的候选不伪装为已回源记录。
5. 旧 GPT 伪标签只留原文件 / 迁移副本，不改名为人工标签。
6. LMDB 提供通用归属出口及训练掩码，尚不是所有 NMRPeak 训练入口的直接替换；不生成坐标，不自动切分数据。
## 2026-09-29 补谱质量检查更新

新增 `nmrexp-link-v2` 检查、无损隔离和跨库同核种比较。修复前生成的导出不自动更新，应重新导出。

推荐全量检查与续跑：`C:\Python314\python.exe -u -m annotation_tool_v3.quality_job full`。不要与旧迁移命令同时运行。`migration_job.json` 显示阶段、心跳和阶段 ETA；限量测试不代表全量完成。

详细结果、阈值、备份和性能测量见 [补谱质量修复与迁移加速](../docs/补谱质量修复与迁移加速_20260929.md)。自动通过不代表已核对原论文，温度与样品身份未报告时仍为未知。

