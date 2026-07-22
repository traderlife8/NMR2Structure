#!/usr/bin/env pwsh
# ============================================================
# setup.ps1 - 一次性安装脚本（每个新 clone 后跑一次）
# ============================================================
# 用途：
#   - 把 git hooks 路径指向 .githooks/（仓库内置、可分享）
#   - 验证 git 配置 + 远程可达性
#
# 使用：
#   cd E:\NMR2Structure\NMR2Structure
#   .\setup.ps1
# ============================================================

$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location -LiteralPath $RepoRoot

Write-Host "[setup] repo root: $RepoRoot" -ForegroundColor Cyan

# 1. 检查 .git 存在
if (-not (Test-Path -LiteralPath (Join-Path $RepoRoot ".git"))) {
    Write-Host "[ERROR] 当前目录不是 git 仓库。" -ForegroundColor Red
    exit 1
}

# 2. 设置 hooks 路径
git config core.hooksPath .githooks
Write-Host "[setup] core.hooksPath = .githooks" -ForegroundColor Green

# 3. 验证 hook 文件
$postCommit = Join-Path $RepoRoot ".githooks/post-commit"
if (-not (Test-Path -LiteralPath $postCommit)) {
    Write-Host "[ERROR] .githooks/post-commit 缺失。" -ForegroundColor Red
    exit 1
}

# 4. 检查 remote
$remoteUrl = git remote get-url origin 2>$null
if ([string]::IsNullOrWhiteSpace($remoteUrl)) {
    Write-Host "[WARN] 未配置 origin remote，跳过推送验证。" -ForegroundColor Yellow
} else {
    Write-Host "[setup] origin = $remoteUrl" -ForegroundColor Cyan
    # 把 token 隐藏显示
    $sanitized = $remoteUrl -replace '://[^@]+@', '://***@'
    Write-Host "[setup]        = $sanitized" -ForegroundColor DarkGray
}

# 5. 列出可用脚本
Write-Host ""
Write-Host "[setup] 可用脚本：" -ForegroundColor Cyan
Write-Host "  .\git_sync.ps1                # 手动：clean + commit + push" -ForegroundColor White
Write-Host "  .\setup.ps1                   # 一次性安装（已完成）" -ForegroundColor White
Write-Host ""
Write-Host "[setup] 从此 git commit 后会自动 push 到 origin（无需手动 git push）。" -ForegroundColor Green
Write-Host "[setup] done." -ForegroundColor Green