#!/usr/bin/env pwsh
# ============================================================
# start_watcher.ps1 - 启动 auto_commit.py 文件监视器
# ============================================================
# 行为：
#   - 在新 PowerShell 窗口中启动 watcher（不阻塞当前窗口）
#   - 防抖默认 5 秒，可通过 -Debounce 调整
#   - 写入 PID 文件 .auto_commit.pid，方便 stop_watcher.ps1 终止
#
# 使用：
#   .\start_watcher.ps1                       # 默认防抖 5 秒
#   .\start_watcher.ps1 -Debounce 10          # 防抖 10 秒
#   .\start_watcher.ps1 -DryRun               # 只观察，不提交
# ============================================================

[CmdletBinding()]
param(
    [double]$Debounce = 5.0,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location -LiteralPath $RepoRoot

$pidFile = Join-Path $RepoRoot ".auto_commit.pid"

# 如果已经在运行，先提示
if (Test-Path -LiteralPath $pidFile) {
    $existing = Get-Content -LiteralPath $pidFile -ErrorAction SilentlyContinue
    if ($existing -and (Get-Process -Id $existing -ErrorAction SilentlyContinue)) {
        Write-Host "[start_watcher] already running (PID=$existing). Use stop_watcher.ps1 first." -ForegroundColor Yellow
        exit 1
    }
    Remove-Item -LiteralPath $pidFile -Force
}

# 检查依赖
try {
    python -X utf8 -c "import watchdog" 2>$null
    if ($LASTEXITCODE -ne 0) {
        throw "watchdog not installed"
    }
} catch {
    Write-Host "[start_watcher] watchdog not installed. Installing..." -ForegroundColor Yellow
    python -X utf8 -m pip install watchdog --index-url https://pypi.tuna.tsinghua.edu.cn/simple
    if ($LASTEXITCODE -ne 0) {
        Write-Host "[ERROR] failed to install watchdog." -ForegroundColor Red
        exit 1
    }
}

# 启动参数
$argsList = @("-X", "utf8", "auto_commit.py", "--debounce", "$Debounce")
if ($DryRun) { $argsList += "--dry-run" }

# 在新窗口中启动
$newTitle = "auto_commit watcher (debounce=${Debounce}s)"
$startArgs = @{
    FilePath     = "powershell.exe"
    ArgumentList = @("-NoProfile", "-Command", "Set-Location -LiteralPath '$RepoRoot'; Write-Host '[watcher] starting...'; python $($argsList -join ' '); Read-Host 'press Enter to close'")
    WindowTitle  = $newTitle
    PassThru     = $true
}

Write-Host "[start_watcher] launching watcher in new window..." -ForegroundColor Cyan
$proc = Start-Process @startArgs

# 等一下让进程起来，再尝试拿 PID
Start-Sleep -Seconds 2
# 注：Start-Process 返回的是 powershell wrapper 的 PID，
#     真正跑 python 的子进程 PID 我们用 tasklist 找
$pythonPid = (Get-CimInstance Win32_Process -Filter "Name='python.exe'" |
    Where-Object { $_.CommandLine -like "*auto_commit.py*" } |
    Select-Object -First 1).ProcessId

if ($pythonPid) {
    Set-Content -LiteralPath $pidFile -Value $pythonPid
    Write-Host "[start_watcher] watcher started (PID=$pythonPid)" -ForegroundColor Green
} else {
    Set-Content -LiteralPath $pidFile -Value $proc.Id
    Write-Host "[start_watcher] wrapper started (PID=$($proc.Id)); python PID not yet visible." -ForegroundColor Yellow
}

Write-Host ""
Write-Host "Stop with:  .\stop_watcher.ps1" -ForegroundColor White
Write-Host "View log:   Get-Content auto_commit.log -Wait" -ForegroundColor White