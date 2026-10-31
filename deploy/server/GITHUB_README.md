# AtomNMR

实验核磁原子归属数据库与标注程序，提供 C/H 双谱标注、AI 自动归属、独立模型测试与数据导出。

本仓库只包含程序代码、字典、部署模板及文档。数据库、实验数据、真实密钥、测试输出、日志和运行环境不随代码上传。

## 启动

在 Python 3.11/3.12 环境安装 `deploy/server/requirements.txt`，将自己的三个数据库放入 `data/assignment_database_v1/databases/`。日志、导出与历史备份位于外层目录：

- `assignments.sqlite`：分子、谱峰、归属及来源。
- `ai_annotations.sqlite`：正式 AI 任务、设置、响应记录。
- `ai_benchmarks.sqlite`：独立测试数据与结果。

桌面使用：`python -m annotation_tool_v3.server --port 8767`，打开 `http://127.0.0.1:8767/`。

Ubuntu 无 sudo 部署请看 [部署说明](docs/Ubuntu无sudo服务器部署_20261030.md)。服务仅监听 localhost，通过 SSH 隧道访问；API Key 通过私密环境文件配置。

## 运行与维护

- 官方 API 支持固定并行与仅谷价运行。仅谷价策略根据历史耗时，在谷价结束前提前停止派发，已发出的请求继续完成。
- 数据库备份由用户自己的工具管理，程序不集成 OpenList 或网盘登录。
- `tools/sync_annotation_github.py` 只按程序文件清单同步代码，不上传 data、outputs、trashbin、真实 .env 或密钥。
- 当前发布分支：`atomnmr-backup`。推送不强制覆盖远端历史。

详细运行策略与数据目录说明见 [维护说明](docs/官方谷价收尾与数据整理_20261030.md)。
