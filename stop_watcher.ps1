#!/usr/bin/env pwsh
# ============================================================
# stop_watcher.ps1 - 停止 auto_commit.py 守护进程
# ============================================================
# 终止策略：
#   1. 读 .auto_commit.pid，结束该进程
#   2. 兜底：扫描所有含 "auto_commit.py" 的 python.exe 进程
# ============================================================

$ErrorActionPreference = "Continue"

$RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$pidFile = Join-Path $RepoRoot ".auto_commit.pid"

$killed = $false

# 1. PID 文件
if (Test-Path -LiteralPath $pidFile) {
    $pid = Get-Content -LiteralPath $pidFile -ErrorAction SilentlyContinue
    if ($pid -and (Get-Process -Id $pid -ErrorAction SilentlyContinue)) {
        Write-Host "[stop_watcher] killing PID $pid..." -ForegroundColor Cyan
        Stop-Process -Id $pid -Force
        $killed = $true
    }
    Remove-Item -LiteralPath $pidFile -Force
}

# 2. 兜底：扫所有 python 进程
$pyProcs = Get-CimInstance Win32_Process -Filter "Name='python.exe'" |
    Where-Object { $_.CommandLine -like "*auto_commit.py*" }
foreach ($p in $pyProcs) {
    Write-Host "[stop_watcher] killing python PID $($p.ProcessId)..." -ForegroundColor Cyan
    Stop-Process -Id $p.ProcessId -Force
    $killed = $true
}

if ($killed) {
    Write-Host "[stop_watcher] done." -ForegroundColor Green
} else {
    Write-Host "[stop_watcher] no watcher process found." -ForegroundColor Yellow
}