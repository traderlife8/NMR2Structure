"""auto_commit.py - 文件系统监视器，文件变更后自动 git add + commit。

设计原则：
  - 触发即 commit，push 由 .githooks/post-commit 接管
  - 防抖（debounce）5 秒：5 秒内无新变更才触发，避免编辑器自动保存产生噪音 commit
  - 排除：.git/、__pycache__/、临时文件、隐藏目录、IDE 锁文件
  - 仅跟踪当前仓库内文件

使用方法：
  python -X utf8 auto_commit.py [--debounce 5] [--dry-run]

停止：Ctrl+C（在启动它的 PowerShell 窗口里）
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

from watchdog.events import FileSystemEvent, FileSystemEventHandler
from watchdog.observers import Observer

REPO_ROOT = Path(__file__).resolve().parent

# 应排除的路径片段（任一匹配则忽略事件）
EXCLUDED_PATH_PARTS = frozenset({
    ".git",
    "__pycache__",
    ".pytest_cache",
    ".venv",
    "venv",
    "node_modules",
    ".idea",
    ".vscode",
    "output",
    "logs",
})

# 应排除的文件后缀（编辑器临时文件、编译产物等）
EXCLUDED_SUFFIXES = frozenset({
    ".pyc", ".pyo", ".pyd",
    ".swp", ".swo",  # vim
    ".tmp", ".bak",
    ".lock",  # IDE 锁文件
})

# 应排除的文件名前缀
EXCLUDED_PREFIXES = frozenset({
    ".#",  # emacs
    "~$",  # Office 临时
})

# 应排除的特定文件名（watcher 自身的运行时产物，防止自触发）
EXCLUDED_NAMES = frozenset({
    "auto_commit.log",      # watcher 自身的日志
    ".auto_commit.pid",     # watcher 自身的 PID 文件
})


def is_excluded(path: Path) -> bool:
    """检查路径是否应被忽略。"""
    rel = path.resolve()
    try:
        parts = set(rel.relative_to(REPO_ROOT).parts)
    except ValueError:
        # 不在仓库内
        return True
    if parts & EXCLUDED_PATH_PARTS:
        return True
    name = path.name
    if name in EXCLUDED_NAMES:
        return True
    if any(name.startswith(p) for p in EXCLUDED_PREFIXES):
        return True
    if any(name.endswith(s) for s in EXCLUDED_SUFFIXES):
        return True
    return False


class DebouncedHandler(FileSystemEventHandler):
    """防抖事件处理器：N 秒内无新事件才触发一次 commit。"""

    def __init__(self, debounce_sec: float, dry_run: bool, log_path: Path):
        super().__init__()
        self.debounce_sec = debounce_sec
        self.dry_run = dry_run
        self.log_path = log_path
        self._timer: threading.Timer | None = None
        self._lock = threading.Lock()
        self._pending: set[str] = set()

    def _log(self, msg: str) -> None:
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {msg}"
        print(line, flush=True)
        try:
            with self.log_path.open("a", encoding="utf-8") as f:
                f.write(line + "\n")
        except OSError:
            pass  # 日志写入失败不影响主流程

    def _on_event(self, event_path: str) -> None:
        p = Path(event_path)
        if is_excluded(p):
            return
        rel = p.relative_to(REPO_ROOT) if p.is_absolute() else p
        with self._lock:
            self._pending.add(str(rel))
            if self._timer is not None:
                self._timer.cancel()
            self._timer = threading.Timer(self.debounce_sec, self._commit)
            self._timer.daemon = True
            self._timer.start()
        self._log(f"queued: {rel} (debounce {self.debounce_sec}s)")

    def _commit(self) -> None:
        with self._lock:
            pending = sorted(self._pending)
            self._pending.clear()
            self._timer = None

        if not pending:
            return

        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        n = len(pending)
        msg = f"auto: update {n} file{'s' if n != 1 else ''} at {ts}"
        # 限制 commit message 长度（git 推荐 < 72 字符首行）
        body = "\n".join(f"  - {p}" for p in pending[:20])
        if len(pending) > 20:
            body += f"\n  ... and {len(pending) - 20} more"
        full_msg = f"{msg}\n\n{body}"

        self._log(f"committing {n} change(s)...")
        if self.dry_run:
            self._log(f"[dry-run] would commit with message:\n{full_msg}")
            return

        try:
            subprocess.run(
                ["git", "add", "-A"],
                cwd=REPO_ROOT, check=True,
                capture_output=True, text=True,
            )
            # -- 防止 commit message 含 shell 特殊字符
            result = subprocess.run(
                ["git", "commit", "-m", msg, "-m", body],
                cwd=REPO_ROOT, check=False,
                capture_output=True, text=True,
            )
            if result.returncode == 0:
                self._log(f"committed: {result.stdout.strip().splitlines()[0]}")
                # post-commit hook 会自动 push
            elif "nothing to commit" in (result.stdout + result.stderr):
                self._log("nothing to commit (changes already staged?)")
            else:
                self._log(f"commit failed: {result.stderr.strip()[:200]}")
        except subprocess.CalledProcessError as e:
            self._log(f"git command failed: {e}")

    # watchdog 事件回调（统一走 _on_event）
    def on_modified(self, event: FileSystemEvent) -> None:
        if not event.is_directory:
            self._on_event(event.src_path)

    def on_created(self, event: FileSystemEvent) -> None:
        if not event.is_directory:
            self._on_event(event.src_path)

    def on_deleted(self, event: FileSystemEvent) -> None:
        if not event.is_directory:
            self._on_event(event.src_path)

    def on_moved(self, event: FileSystemEvent) -> None:
        if not event.is_directory:
            self._on_event(event.dest_path)

    def shutdown(self) -> None:
        with self._lock:
            if self._timer is not None:
                self._timer.cancel()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--debounce", type=float, default=5.0,
        help="debounce window in seconds (default: 5.0)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="detect changes but don't actually commit",
    )
    args = parser.parse_args()

    if not (REPO_ROOT / ".git").exists():
        print(f"ERROR: {REPO_ROOT} is not a git repo", file=sys.stderr)
        return 1

    log_path = REPO_ROOT / "auto_commit.log"
    handler = DebouncedHandler(args.debounce, args.dry_run, log_path)
    observer = Observer()
    observer.schedule(handler, str(REPO_ROOT), recursive=True)
    observer.start()

    print(f"[auto_commit] watching: {REPO_ROOT}")
    print(f"[auto_commit] debounce: {args.debounce}s")
    print(f"[auto_commit] log file: {log_path}")
    print(f"[auto_commit] press Ctrl+C to stop")
    print(f"[auto_commit] post-commit hook handles push")

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("\n[auto_commit] shutting down...")
        handler.shutdown()
        observer.stop()
    observer.join()
    print("[auto_commit] stopped.")
    return 0


if __name__ == "__main__":
    sys.exit(main())