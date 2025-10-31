#!/usr/bin/env bash
# 一键在本地创建仓库并推送（依赖 gh 已成功登录）
set -euo pipefail

OWNER="ShawnSHAWN-in"
REPO="Xinyuan_Cui_transcriptome"
MSG="Add transcriptome pipeline scripts and README"

# 检查 gh 是否登录
if ! command -v gh >/dev/null 2>&1; then
    echo "请先安装 GitHub CLI (gh): https://cli.github.com/" >&2
    exit 1
fi

if ! gh auth status >/dev/null 2>&1; then
    echo "请先运行: gh auth login" >&2
    exit 1
fi

# 创建仓库（如果已存在，gh 会返回非0；这里忽略错误以继续）
gh repo create "${OWNER}/${REPO}" --public --confirm || true

# 初始化 git 并推送
git init
git add .
git commit -m "${MSG}" || true
git branch -M main || true

# 确保 remote 指向
if ! git remote get-url origin >/dev/null 2>&1; then
    git remote add origin "git@github.com:${OWNER}/${REPO}.git" || git remote add origin "https://github.com/${OWNER}/${REPO}.git"
fi

git push -u origin main
echo "Pushed to github.com/${OWNER}/${REPO}"