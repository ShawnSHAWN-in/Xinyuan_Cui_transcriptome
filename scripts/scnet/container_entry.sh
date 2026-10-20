#!/usr/bin/env bash
set -euo pipefail
umask 077

# 日志持久化
mkdir -p "$SCNET_WORK_ROOT/runs/$SCNET_RUN_KEY"
run_dir="$SCNET_WORK_ROOT/runs/$SCNET_RUN_KEY"
project_dir="$SCNET_WORK_ROOT/project"
exec >> "$run_dir/run.log" 2>&1
finish_run() {
    exit_code=$?
    trap - EXIT
    printf '%s\n' "$exit_code" > "$run_dir/exit_code.tmp"
    mv "$run_dir/exit_code.tmp" "$run_dir/exit_code"
    date -u +%FT%TZ > "$run_dir/finished_at"
    exit "$exit_code"
}
trap finish_run EXIT
date -u +%FT%TZ > "$run_dir/started_at"
printf '%s\n' "$SCNET_COMMIT" > "$run_dir/commit"
printf 'Commit: %s\n' "$SCNET_COMMIT"

# 串行更新项目
command -v git
command -v flock
exec 9> "$SCNET_WORK_ROOT/.project.lock"
flock -x 9
if [[ ! -e "$project_dir" ]]; then
    git clone "$SCNET_REPO_URL" "$project_dir"
fi
cd "$project_dir"
[[ "$(git rev-parse --show-toplevel)" == "$(pwd -P)" ]]
if [[ -n "$(git status --porcelain --untracked-files=no)" ]]; then
    printf 'Tracked files have local changes; commit or move them before running.\n' >&2
    exit 1
fi

# 固定触发版本
git fetch --no-tags "$SCNET_REPO_URL" "$SCNET_COMMIT"
git checkout --detach "$SCNET_COMMIT"
[[ "$(git rev-parse HEAD)" == "$SCNET_COMMIT" ]]
printf 'Project: %s\nCommand: %s\n' "$project_dir" "$SCNET_RUN_COMMAND"
bash -euo pipefail -c "$SCNET_RUN_COMMAND"
