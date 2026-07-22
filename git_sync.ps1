#!/usr/bin/env pwsh
# ============================================================
# git_sync.ps1 - 一键提交并推送到 GitHub
# ============================================================
# 使用方法：
#   1. 修改完代码后，在仓库根目录运行：  .\git_sync.ps1
#   2. 或带提交信息：                  .\git_sync.ps1 -Message "fix: xxx"
#
# 功能：
#   - 清理 __pycache__ 等临时缓存（保持仓库干净）
#   - 校验是否有真实变更（无变更则跳过）
#   - 自动 git add + commit（无变更时静默退出）
#   - 自动 git push（如未配置 remote 则提示）
#
# 安全提示：
#   - 切勿把 Personal Access Token 直接写进此脚本或 commit
#   - 推荐首次 push 后配置凭据存储（git credential manager）
#   - 如需更换 token，编辑 .git/config 中的 remote URL
# ============================================================

[CmdletBinding()]
param(
    [string]$Message = "",
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

# 切到脚本所在目录（即仓库根）
$RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location -LiteralPath $RepoRoot

Write-Host "[git_sync] repo root: $RepoRoot" -ForegroundColor Cyan

# 检查是否在 git 仓库中
$gitDir = Join-Path $RepoRoot ".git"
if (-not (Test-Path -LiteralPath $gitDir)) {
    Write-Host "[ERROR] 当前目录不是 git 仓库（缺少 .git）。" -ForegroundColor Red
    exit 1
}

# 清理临时缓存（不删除 .gitignore 已声明的源文件）
$pycacheDirs = Get-ChildItem -Path $RepoRoot -Recurse -Directory -Filter "__pycache__" -ErrorAction SilentlyContinue
foreach ($d in $pycacheDirs) {
    Remove-Item -LiteralPath $d.FullName -Recurse -Force -ErrorAction SilentlyContinue
    Write-Host "[clean] removed $($d.FullName)" -ForegroundColor DarkGray
}

# 查看是否有变更
$status = git status --porcelain
if ([string]::IsNullOrWhiteSpace($status)) {
    Write-Host "[git_sync] 没有变更，无需提交。" -ForegroundColor Yellow
    exit 0
}

Write-Host "[git_sync] 待提交变更：" -ForegroundColor Cyan
git status --short

# 默认提交信息
if ([string]::IsNullOrWhiteSpace($Message)) {
    $ts = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    $Message = "update: $ts"
}

if ($DryRun) {
    Write-Host "[dry-run] 将使用提交信息：$Message" -ForegroundColor Yellow
    exit 0
}

# add + commit
git add .
git commit -m $Message
if ($LASTEXITCODE -ne 0) {
    Write-Host "[ERROR] git commit 失败。" -ForegroundColor Red
    exit $LASTEXITCODE
}

# push
$remoteUrl = git remote get-url origin 2>$null
if ([string]::IsNullOrWhiteSpace($remoteUrl)) {
    Write-Host "[WARN] 未配置 origin remote，跳过 push。" -ForegroundColor Yellow
    Write-Host "       本次 commit 仍保留在本地。" -ForegroundColor Yellow
    exit 0
}

Write-Host "[git_sync] pushing to $remoteUrl ..." -ForegroundColor Cyan
git push origin main
if ($LASTEXITCODE -ne 0) {
    Write-Host "[ERROR] git push 失败。请检查网络和凭据。" -ForegroundColor Red
    exit $LASTEXITCODE
}

Write-Host "[git_sync] done." -ForegroundColor Green