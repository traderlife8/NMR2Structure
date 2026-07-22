"""项目统一日志初始化.

设计要点:
  1. 每次调用 setup_logging() 在指定目录生成一个独立日志文件
     (默认 E:\\NMR2Structure\\logs\\NMR2Structure_YYYY-MM-DD_HH-MM-SS.log),
     避免多次运行相互覆盖.
  2. 通过 ``_Tee`` 包装 sys.stdout / sys.stderr,
     把所有 print() 输出和所有 logging 输出都同时写到控制台 + 日志文件,
     无需修改现有 25+ 处 print 调用.
  3. 包装 sys.excepthook, 把未捕获异常 (含 traceback) 写入日志.
  4. 项目 logger 名: ``nmr_elucidator``, 不向 root 冒泡,
     避免污染 matplotlib / rdkit 等第三方库.

使用:
    from nmr_elucidator.log_setup import setup_logging
    log_file = setup_logging()                # 用默认路径
    log_file = setup_logging(log_dir='D:/mylogs')
    log_file = setup_logging(also_tee=False)   # 只写文件, 不重定向屏幕
"""
from __future__ import annotations

import datetime as _dt
import io
import logging
import os
import sys
import traceback
from pathlib import Path
from typing import IO, Optional, Tuple


# ---------- 默认值 ----------

DEFAULT_LOG_DIR = r'E:\NMR2Structure\logs'
DEFAULT_NAME_PREFIX = 'NMR2Structure'
LOGGER_NAME = 'nmr_elucidator'

# 已安装的 Tee 实例 (防止 main() 嵌套调用时重复包装 stdout)
_TEE_INSTALLED: Optional[Tuple[int, int]] = None

# 当前活动 log file 路径 (全局读取, 给子进程或调试用)
_CURRENT_LOG_FILE: Optional[Path] = None


# ---------- Tee 包装 ----------

class _TeeTextIO(io.TextIOBase):
    """把文本流同时写到原 stream + 日志文件.

    兼容 sys.stdout / sys.stderr, 同时作为 file= 参数传给 print() 使用.
    不修改 isatty/fileno/encoding 等元属性, 让 IPython / pytest 等工具
    仍能识别它为真实终端.
    """

    def __init__(self, original: IO[str], log_path: Path) -> None:
        self._orig = original
        self._log_path = log_path
        self._fh = open(log_path, 'a', encoding='utf-8', buffering=1)

    @property
    def log_path(self) -> Path:
        return self._log_path

    def write(self, s: str) -> int:
        # 终端写出 (容忍 UnicodeEncodeError, 如 Windows GBK cmd)
        try:
            self._orig.write(s)
        except (UnicodeEncodeError, OSError):
            # 退化为 errors='replace'
            try:
                enc = getattr(self._orig, 'encoding', 'utf-8') or 'utf-8'
                self._orig.write(s.encode(enc, errors='replace').decode(enc))
            except Exception:
                pass
        # 日志写出 (永远 utf-8)
        try:
            self._fh.write(s)
        except Exception:
            pass
        return len(s)

    def flush(self) -> None:
        try:
            self._orig.flush()
        except Exception:
            pass
        try:
            self._fh.flush()
        except Exception:
            pass

    def __del__(self) -> None:
        # 程序退出时确保文件句柄关闭
        try:
            self._fh.close()
        except Exception:
            pass

    # 透传元属性 (保持 IPython / pytest / IDE 调试器兼容性)
    def isatty(self) -> bool:
        return getattr(self._orig, 'isatty', lambda: False)()

    def fileno(self) -> int:
        return self._orig.fileno()

    @property
    def encoding(self) -> str:
        return getattr(self._orig, 'encoding', 'utf-8') or 'utf-8'

    def writable(self) -> bool:
        return True

    def __getattr__(self, name: str):
        return getattr(self._orig, name)


# ---------- 日志初始化 ----------

def setup_logging(
    log_dir: Optional[str] = None,
    name_prefix: Optional[str] = None,
    also_tee: bool = True,
    console_level: int = logging.INFO,
    file_level: int = logging.DEBUG,
) -> Path:
    """初始化项目日志. 每次调用生成新的日志文件.

    Args:
        log_dir: 日志目录 (默认 E:\\NMR2Structure\\logs).
        name_prefix: 文件名前缀 (默认 NMR2Structure).
        also_tee: 是否重定向 stdout / stderr 到日志文件
                 (True = 屏幕上的所有内容都同步写入日志, 包括 print).
        console_level: 控制台输出级别 (默认 INFO).
        file_level: 日志文件记录级别 (默认 DEBUG).

    Returns:
        生成的日志文件绝对路径.
    """
    global _CURRENT_LOG_FILE, _TEE_INSTALLED

    log_dir_str = log_dir or os.environ.get('NMR_ELUCIDATOR_LOG_DIR') or DEFAULT_LOG_DIR
    prefix = name_prefix or os.environ.get('NMR_ELUCIDATOR_LOG_NAME') or DEFAULT_NAME_PREFIX

    log_dir_path = Path(log_dir_str)
    try:
        log_dir_path.mkdir(parents=True, exist_ok=True)
    except OSError as ex:
        # 用户目录无写权限时, 退回到当前目录
        fallback = Path.cwd() / 'logs'
        fallback.mkdir(parents=True, exist_ok=True)
        print(f'[log_setup] WARN: 无法创建 {log_dir_path} ({ex}); 退回到 {fallback}', file=sys.stderr)
        log_dir_path = fallback

    # 文件名: NMR2Structure_2026-07-22_14-30-15.log
    ts = _dt.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    log_file = log_dir_path / f'{prefix}_{ts}.log'

    # 确保文件存在后再包装 (避免两次写入导致路径不一致)
    log_file.touch()

    # 1) Tee 输出 (捕获所有 print 原文)
    if also_tee and _TEE_INSTALLED is None:
        try:
            sys.stdout = _TeeTextIO(sys.stdout, log_file)
            sys.stderr = _TeeTextIO(sys.stderr, log_file)
            _TEE_INSTALLED = (id(sys.stdout), id(sys.stderr))
        except Exception as ex:
            print(f'[log_setup] WARN: stdout/stderr 包装失败: {ex}', file=sys.stderr)
    elif also_tee and _TEE_INSTALLED is not None:
        # 已包装过 (嵌套调用), 跳过
        pass

    # 2) 配置项目 logger
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(logging.DEBUG)
    logger.propagate = False  # 不向 root 冒泡

    # 清掉可能存在的旧 handler (避免重复)
    for h in list(logger.handlers):
        logger.removeHandler(h)

    # File handler: 详细格式 + DEBUG
    try:
        fh = logging.FileHandler(str(log_file), encoding='utf-8')
        fh.setLevel(file_level)
        fh.setFormatter(logging.Formatter(
            fmt='%(asctime)s [%(levelname)-7s] %(name)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
        ))
        logger.addHandler(fh)
    except Exception as ex:
        print(f'[log_setup] WARN: file handler 注册失败: {ex}', file=sys.stderr)

    # Stream handler: 输出到 stderr (统一用 stderr 与 print 视觉一致)
    # 关键: 必须走 sys.__stderr__ (原 stream), 避免被 _TeeTextIO 包装造成双写
    orig_stderr = sys.__stderr__ if also_tee and _TEE_INSTALLED is None else sys.stderr
    sh = logging.StreamHandler(stream=orig_stderr)
    sh.setLevel(console_level)
    sh.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(sh)

    # 3) 包装未捕获异常 (traceback 写入日志)
    _install_excepthook(logger)

    # 4) 写入文件头
    logger.info('=' * 78)
    logger.info(f'日志文件: {log_file}')
    logger.info(f'启动时间: {_dt.datetime.now().isoformat(timespec="seconds")}')
    logger.info(f'命令行参数: {" ".join(sys.argv)!r}')
    logger.info(f'Python: {sys.version.split()[0]}  平台: {sys.platform}')
    logger.info(f'工作目录: {os.getcwd()}')
    logger.info('=' * 78)

    _CURRENT_LOG_FILE = log_file
    return log_file


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """获取项目 logger.

    Args:
        name: 子模块名 (如 'evidence'); 不传则返回根项目 logger.
    """
    if name is None:
        return logging.getLogger(LOGGER_NAME)
    return logging.getLogger(f'{LOGGER_NAME}.{name}')


def get_current_log_file() -> Optional[Path]:
    """返回当前活动日志文件路径; 若未初始化则返回 None."""
    return _CURRENT_LOG_FILE


def _install_excepthook(logger: logging.Logger) -> None:
    """把未捕获异常写入日志, 再调用系统默认 handler."""
    default_hook = sys.excepthook

    def _hook(exc_type, exc_value, exc_tb):
        try:
            tb_lines = ''.join(traceback.format_exception(exc_type, exc_value, exc_tb))
            logger.critical('未捕获异常:\n%s', tb_lines)
            logger.handlers[0].flush() if logger.handlers else None  # 强制 flush 文件
        except Exception:
            pass
        default_hook(exc_type, exc_value, exc_tb)

    sys.excepthook = _hook


# ---------- 便捷自动初始化 (供子模块 import 时使用) ----------

_AUTO_INITIALIZED = False


def auto_setup(*args, **kwargs) -> Path:
    """仅在尚未初始化时执行 setup_logging, 用于子模块顶部调用.

    多次调用是安全的 (幂等). 仅首次生效.
    """
    global _AUTO_INITIALIZED
    if _AUTO_INITIALIZED:
        return _CURRENT_LOG_FILE or Path('.')
    _AUTO_INITIALIZED = True
    return setup_logging(*args, **kwargs)


if __name__ == '__main__':
    # 单独运行此脚本用于调试 / 自检
    p = setup_logging(console_level=logging.DEBUG)
    log = get_logger('demo')
    log.debug('debug message (仅文件)')
    log.info('info message (屏幕+文件)')
    log.warning('warning message')
    print('这条 print() 也应同时出现在屏幕上和日志中')
    print(f'日志文件路径: {p}')
