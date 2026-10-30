# Ubuntu服务器，无sudo权限：按此安装

以SSH普通用户、目录 `~/atomnmr/` 为准。代码包和数据库准备方法见原部署手册。下面不需要创建系统账号，不改系统Python，也不安装系统服务。

## 一、上传

本机先暂停标注、保存修改，准备三个一致性快照；不要直接复制仍有未合并WAL的活跃数据库。在本机项目目录PowerShell执行，替换用户名和服务器地址：

```powershell
ssh 用户名@服务器地址 "mkdir -p ~/atomnmr_upload"
scp .\deploy\atomnmr-server-code-20261030.zip 用户名@服务器地址:~/atomnmr_upload/
scp E:\atomnmr_transfer\assignments.sqlite E:\atomnmr_transfer\ai_annotations.sqlite E:\atomnmr_transfer\ai_benchmarks.sqlite 用户名@服务器地址:~/atomnmr_upload/
```

`E:\atomnmr_transfer`只是示例快照目录，按实际位置替换。服务器端以下命令按Bash执行。

## 二、用户目录与Python环境

```bash
mkdir -p ~/atomnmr/{app,data,secrets,logs}
unzip ~/atomnmr_upload/atomnmr-server-code-20261030.zip -d ~/atomnmr/app
cp -n ~/atomnmr_upload/assignments.sqlite ~/atomnmr/data/
cp -n ~/atomnmr_upload/ai_annotations.sqlite ~/atomnmr/data/
cp -n ~/atomnmr_upload/ai_benchmarks.sqlite ~/atomnmr/data/
chmod 700 ~/atomnmr/secrets
chmod 600 ~/atomnmr/data/*.sqlite
python3 --version
python3 -m venv ~/atomnmr/venv
~/atomnmr/venv/bin/python -m pip install --upgrade pip
~/atomnmr/venv/bin/python -m pip install -r ~/atomnmr/app/deploy/server/requirements.txt
```

推荐Python3.11/3.12。`cp -n`不覆盖已有库；这套命令针对首次安装，已有运行库时不要当成升级数据库的流程。

若没有unzip，用已有Python解压：

```bash
python3 -m zipfile -e ~/atomnmr_upload/atomnmr-server-code-20261030.zip ~/atomnmr/app
```

若venv提示缺少ensurepip或python3-venv，无sudo时不要执行apt安装。先检查组里已有的conda、mamba或其他Python3.12路径。如果已有conda，可在指定路径建立环境：

```bash
conda create -p ~/atomnmr/venv python=3.12 pip -y
~/atomnmr/venv/bin/python -m pip install -r ~/atomnmr/app/deploy/server/requirements.txt
```

两种环境创建方式选其一，不叠加。若创建失败留下了半成品环境，换一个空目录建立并同步修改启动脚本中的解释器路径；不直接删除或覆盖系统环境。

## 三、生成用户路径配置

```bash
cp -n ~/atomnmr/app/deploy/server/atomnmr.env.example ~/atomnmr/secrets/atomnmr.env
chmod 600 ~/atomnmr/secrets/atomnmr.env
python3 - <<'PY'
import json
from pathlib import Path
root=Path.home()/'atomnmr'
target=root/'api.json'
if target.exists():
    raise SystemExit('api.json already exists; edit it instead of overwriting')
cfg=json.loads((root/'app/deploy/server/api.example.json').read_text())
cfg['api']['ca_file']=str(root/'app/annotation_tool_v3/certificates/ai-cluster-ca.crt')
target.write_text(json.dumps(cfg,indent=2))
PY
nano ~/atomnmr/secrets/atomnmr.env
```

将ATOMNMR_API_KEY改成真实组内密钥，例如 `ATOMNMR_API_KEY='实际密钥'`，保留私密文件权限，不把密钥放在命令行参数中。首次ATOMNMR_RESUME_JOB_ID留空。api.json默认固定2并行、组内本地适配、V3、关闭思考。

## 四、先前台启动验证

```bash
bash ~/atomnmr/app/deploy/server/start-user.sh
```

另开一个SSH窗口检查：

```bash
curl --fail http://127.0.0.1:8767/api/health
```

确认成功后，在前台窗口Ctrl+C停止，再选择下一节的一种后台方式。不要前台、nohup和systemd同时启动。同库新入口有文件锁，但不要同时用旧server入口绕过锁。

## 五、后台方式二选一

### 优先：用户级systemd（可用时）

先检查，不需要sudo：

```bash
systemctl --user status
loginctl show-user "$(id -un)" -p Linger
```

如果Linger=no，可先尝试为自己的用户启用（是否允许取决于服务器策略，不需要在命令前加sudo）：

```bash
loginctl enable-linger "$(id -un)"
loginctl show-user "$(id -un)" -p Linger
```

如提示权限不足，不要当成已启用，改用下方备选方式并保留常驻限制。用户服务可用且Linger=yes时：

```bash
mkdir -p ~/.config/systemd/user
cp ~/atomnmr/app/deploy/server/atomnmr-user.service ~/.config/systemd/user/atomnmr.service
systemctl --user daemon-reload
systemctl --user enable --now atomnmr
systemctl --user status atomnmr --no-pager
journalctl --user -u atomnmr -f
```

这比nohup更适合常驻：进程失败时重启，并在用户管理器存在时自动启动。**无sudo不代表用户级systemd一定可用；Linger=no时，不能保证退出全部会话或服务器重启后继续运行。** 如果站点策略不允许启用linger，不要把它当成已配置的开机启动。

### 通用备选：nohup

```bash
nohup bash ~/atomnmr/app/deploy/server/start-user.sh >> ~/atomnmr/logs/server.log 2>&1 < /dev/null &
echo $! > ~/atomnmr/server.pid
curl --fail http://127.0.0.1:8767/api/health
```

通常可在断开SSH后继续运行，但不会提供进程故障重启、服务器重启恢复；站点若会在退出时清理用户进程，也可能被终止。不要把nohup等同于可靠的开机常驻。

查看日志：`tail -n 100 ~/atomnmr/logs/server.log`。长期运行要轮转日志，不要任其无限增长。若需结束，先在网页暂停任务，再核对`/api/health`返回的PID与启动命令，确认是自己的标注服务后再结束对应进程。

## 六、本机查看、开始任务

本机PowerShell：

```powershell
ssh -L 18767:127.0.0.1:8767 用户名@服务器地址
```

打开 `http://127.0.0.1:18767/auto.html`，沿用当前配置，检查组内模型后开始持续标注。服务器只监听本机，不需要开放防火墙端口。

开始后记录持续任务ID，将其写入 `~/atomnmr/secrets/atomnmr.env` 的ATOMNMR_RESUME_JOB_ID。用户服务故障重启时，仅自动恢复此任务的interrupted状态，失败/结果未知项不自动重试；paused和attention仍需人工处理。正常停服后paused也需网页恢复。

## 七、验收

先跑少量任务、核对数据库确实增加结果，再断开SSH后重连查看健康接口和任务进度。若选择用户systemd，还需确认Linger及站点的进程生命周期政策，才能宣称具备注销后常驻或开机启动能力。此项目的任务异常保护不会因为使用systemd而消失。

本机与服务器不同时写两份分叉库；每日将服务器三库的一致性快照备份到另一存储位置。全部操作都应针对自己的用户目录；如果这是Slurm等计算集群，长期服务还需遵守登录节点/计算节点运行规则。
