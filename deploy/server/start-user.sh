#!/usr/bin/env bash
set -euo pipefail
# Installed at ~/atomnmr/app/deploy/server/start-user.sh.
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
task_root="${ATOMNMR_ROOT:-$(cd -- "$script_dir/../../.." && pwd)}"
cd -- "$task_root/app"
test -f "$task_root/secrets/atomnmr.env" || { echo 'Missing secrets/atomnmr.env' >&2; exit 1; }
set -a
source "$task_root/secrets/atomnmr.env"
set +a
exec "$task_root/venv/bin/python" -u -m annotation_tool_v3.server_linux \
  --db "$task_root/data/assignments.sqlite" --config "$task_root/api.json" --port 8767
