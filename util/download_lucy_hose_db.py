#!/usr/bin/env python
"""Download lucy-ng HOSE database (lucy-ng-derep.db, ~2.8 GB).

该脚本提供 4 种下载策略, 按推荐度排序:

1. (--strategy colab)   输出 Google Colab Notebook JSON, 在浏览器一键运行
2. (--strategy aria2)   生成 aria2 多线程命令行 (国内最快, 需先安装 aria2)
3. (--strategy curl)    直接 curl/wget (适合有代理/VPN 的环境)
4. (--strategy proxy)   通过 HTTP/SOCKS5 代理下载

figshare 没有官方国内镜像. 对于无法访问外网的环境, 推荐使用 --strategy colab,
在 Google Colab 免费 GPU 环境中下载到 Google Drive, 再从 Google Drive 网页下载到本地.

Usage:
    # 输出 Colab notebook JSON + 打印操作步骤
    python util/download_lucy_hose_db.py --strategy colab

    # 生成 aria2 命令行 (复制到 PowerShell 直接运行)
    python util/download_lucy_hose_db.py --strategy aria2 --output aria2_cmd.txt

    # 用代理下载
    python util/download_lucy_hose_db.py --strategy proxy \\
        --proxy http://127.0.0.1:7890 \\
        --output E:\\NMR2Structure\\external\\lucy-ng\\data\\reference\\
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from textwrap import dedent

# 数据库元数据 (来自 lucy-ng 源码 DATABASE_URL)
DATABASE_URL = 'https://ndownloader.figshare.com/files/61078393'
DATABASE_DOI = '10.6084/m9.figshare.31073554'
DATABASE_SIZE_MB = 830
DATABASE_UNZIPPED_GB = 2.8
DEFAULT_OUTPUT_DIR = r'E:\NMR2Structure\external\lucy-ng\data\reference'


# ---------- 策略 1: Google Colab Notebook ----------

COLAB_NOTEBOOK = {
    'cells': [
        {
            'cell_type': 'markdown',
            'metadata': {},
            'source': [
                '# lucy-ng HOSE 数据库下载 (Google Colab 曲线下载)\n',
                '\n',
                '**适用人群**: 国内用户, 无法直接访问 figshare, 但能访问 colab.research.google.com\n',
                '\n',
                '**目标**: 把 ~830 MB 的 `lucy-ng-derep.db.zip` 下载到你的 Google Drive, 再从网页下载到本地.\n',
                '\n',
                '**步骤**:\n',
                '1. 点击下方 ▶ 运行第一个 Cell (挂载 Google Drive)\n',
                '2. 按提示登录 Google 账号授权\n',
                '3. 点击下方 ▶ 运行第二个 Cell (开始下载, 约 5-10 分钟)\n',
                '4. 打开 https://drive.google.com/drive/u/0/my-drive 找到 `lucy-ng-derep.db.zip`\n',
                '5. 右键 → 下载 (Google Drive 国内下载速度约 1-3 MB/s, ~5 分钟)\n',
                '\n',
                '**直链来源**: lucy-ng 源码 `src/lucy_ng/cli/database.py:20`\n',
                '**DOI**: 10.6084/m9.figshare.31073554\n',
            ],
        },
        {
            'cell_type': 'code',
            'execution_count': None,
            'metadata': {},
            'outputs': [],
            'source': [
                'from google.colab import drive\n',
                'drive.mount("/content/drive")\n',
            ],
        },
        {
            'cell_type': 'code',
            'execution_count': None,
            'metadata': {},
            'outputs': [],
            'source': [
                '# 下载到 Google Drive 根目录\n',
                '!mkdir -p /content/drive/MyDrive/lucy-ng\n',
                f'!wget -c -O /content/drive/MyDrive/lucy-ng/lucy-ng-derep.db.zip {DATABASE_URL}\n',
                '!ls -lh /content/drive/MyDrive/lucy-ng/\n',
            ],
        },
        {
            'cell_type': 'code',
            'execution_count': None,
            'metadata': {},
            'outputs': [],
            'source': [
                '# 验证 zip 完整性 (可选)\n',
                '!unzip -t /content/drive/MyDrive/lucy-ng/lucy-ng-derep.db.zip 2>&1 | tail -20\n',
            ],
        },
    ],
    'metadata': {
        'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
        'language_info': {'name': 'python', 'version': '3.10'},
        'colab': {
            'provenance': [],
            'gpu_visible': False,
        },
    },
    'nbformat': 4,
    'nbformat_minor': 0,
}


def strategy_colab(output_dir: str) -> int:
    """生成 Google Colab notebook JSON 文件."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    nb_path = out_dir / 'download_lucy_hose_db.ipynb'

    nb_path.write_text(
        json.dumps(COLAB_NOTEBOOK, indent=1, ensure_ascii=False),
        encoding='utf-8',
    )

    print(f'[OK] 已生成 Colab Notebook: {nb_path}', file=sys.stderr)
    print('', file=sys.stderr)
    print('=' * 70, file=sys.stderr)
    print('使用步骤:', file=sys.stderr)
    print('=' * 70, file=sys.stderr)
    print('1. 打开 https://colab.research.google.com/', file=sys.stderr)
    print('2. 菜单 File → Upload Notebook → 选择上述 .ipynb 文件', file=sys.stderr)
    print('3. 依次运行三个 Cell (首个 Cell 会要求 Google 账号授权)', file=sys.stderr)
    print('4. 下载完成后到 https://drive.google.com/drive/u/0/my-drive 找文件', file=sys.stderr)
    print('5. 右键 `lucy-ng-derep.db.zip` → 下载到本地 Windows', file=sys.stderr)
    print('6. 用 7-Zip 解压到: ' + DEFAULT_OUTPUT_DIR, file=sys.stderr)
    print('7. 解压后路径: ' + DEFAULT_OUTPUT_DIR + '\\lucy-ng-derep.db', file=sys.stderr)
    print('', file=sys.stderr)
    print('=' * 70, file=sys.stderr)
    print('或者一键 Google Drive 直传链接:', file=sys.stderr)
    print('=' * 70, file=sys.stderr)
    print(dedent('''
        如果你愿意, 我可以提供 .ipynb 的内容, 你直接粘到 Colab 单元格中运行.

        也可以使用 Google Drive 直传 (无需 Colab, 需 Google 账号):
        - 把 figshare 下载链接发到自己的 Gmail
        - 在 Gmail 里点击附件, "添加到 Drive"
        - 从 Drive 下载到本地

        最简单方式: 把 figshare URL 发给任何海外的朋友/同学/老师,
        让他们 wget 下来后通过网盘传给你.''').strip(), file=sys.stderr)
    return 0


# ---------- 策略 2: aria2 多线程 ----------

def strategy_aria2(output_dir: str, output_file: str | None) -> int:
    """生成 aria2 命令行脚本."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / 'lucy-ng-derep.db.zip'

    script_lines = [
        '# aria2 多线程下载脚本 (Windows PowerShell)',
        '# 用法: 先安装 aria2 (choco install aria2 或 scoop install aria2)',
        '# 然后运行: aria2c --input=aria2_cmd.txt',
        '',
        f'{DATABASE_URL}',
        f'  --out=lucy-ng-derep.db.zip',
        f'  --dir={out_dir}',
        '  --max-connection-per-server=16',
        '  --split=16',
        '  --min-split-size=1M',
        '  --continue=true',
        '  --user-agent=Mozilla/5.0',
        '',
        '# 解压:',
        f'# Expand-Archive -Path "{target}" -DestinationPath "{out_dir}"',
    ]

    out_file = Path(output_file) if output_file else out_dir / 'aria2_cmd.txt'
    out_file.write_text('\n'.join(script_lines), encoding='utf-8')

    print(f'[OK] 已生成 aria2 命令脚本: {out_file}', file=sys.stderr)
    print('', file=sys.stderr)
    print('=' * 70, file=sys.stderr)
    print('aria2 安装:', file=sys.stderr)
    print('=' * 70, file=sys.stderr)
    print(dedent('''
        # 方法 1: scoop (推荐, 已自带)
        scoop install aria2

        # 方法 2: choco
        choco install aria2

        # 方法 3: 直接下载
        # https://github.com/aria2/aria2/releases
    ''').strip(), file=sys.stderr)
    print('=' * 70, file=sys.stderr)
    print('运行:', file=sys.stderr)
    print('=' * 70, file=sys.stderr)
    print(f'  aria2c --input="{out_file}"', file=sys.stderr)
    return 0


# ---------- 策略 3: curl/wget 直接下载 ----------

def strategy_curl(output_dir: str) -> int:
    """打印 curl/wget 命令行, 用户自行执行."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / 'lucy-ng-derep.db.zip'

    print('curl/wget 直接下载命令 (需要能访问 figshare 的网络):', file=sys.stderr)
    print('', file=sys.stderr)
    print(f'  curl -L -C - -o "{target}" {DATABASE_URL}', file=sys.stderr)
    print('', file=sys.stderr)
    print(f'  # 或 wget', file=sys.stderr)
    print(f'  wget -c -O "{target}" {DATABASE_URL}', file=sys.stderr)
    print('', file=sys.stderr)
    print(f'# 解压: Expand-Archive -Path "{target}" -DestinationPath "{out_dir}"', file=sys.stderr)
    return 0


# ---------- 策略 4: 代理下载 ----------

def strategy_proxy(output_dir: str, proxy_url: str) -> int:
    """通过 HTTP/SOCKS5 代理下载."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / 'lucy-ng-derep.db.zip'

    print(f'[INFO] 使用代理: {proxy_url}', file=sys.stderr)
    print('', file=sys.stderr)
    print('curl 命令 (设置代理环境变量):', file=sys.stderr)
    print('', file=sys.stderr)
    print(f'  $env:HTTP_PROXY = "{proxy_url}"', file=sys.stderr)
    print(f'  $env:HTTPS_PROXY = "{proxy_url}"', file=sys.stderr)
    print(f'  curl -L -C - -o "{target}" {DATABASE_URL}', file=sys.stderr)
    print('', file=sys.stderr)
    print('# 或在 Python requests 中设置 proxies 参数 (见 lucy-ng CLI 源码)', file=sys.stderr)
    return 0


# ---------- 主入口 ----------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description='下载 lucy-ng HOSE 数据库 (~830 MB / ~2.8 GB 解压)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--strategy', default='colab',
        choices=['colab', 'aria2', 'curl', 'proxy'],
        help='下载策略 (默认 colab, 国内最稳)',
    )
    parser.add_argument(
        '--output-dir', default=DEFAULT_OUTPUT_DIR,
        help=f'输出目录 (默认 {DEFAULT_OUTPUT_DIR})',
    )
    parser.add_argument(
        '--output-file', default=None,
        help='额外输出文件 (如 aria2 脚本路径)',
    )
    parser.add_argument(
        '--proxy', default=None,
        help='代理 URL (仅 --strategy proxy 生效, 如 http://127.0.0.1:7890)',
    )

    args = parser.parse_args(argv)

    print(f'数据库 URL: {DATABASE_URL}', file=sys.stderr)
    print(f'数据库 DOI: {DATABASE_DOI}', file=sys.stderr)
    print(f'大小: ~{DATABASE_SIZE_MB} MB (压缩) / ~{DATABASE_UNZIPPED_GB} GB (解压)', file=sys.stderr)
    print(f'输出目录: {args.output_dir}', file=sys.stderr)
    print('', file=sys.stderr)

    if args.strategy == 'colab':
        return strategy_colab(args.output_dir)
    elif args.strategy == 'aria2':
        return strategy_aria2(args.output_dir, args.output_file)
    elif args.strategy == 'curl':
        return strategy_curl(args.output_dir)
    elif args.strategy == 'proxy':
        if not args.proxy:
            print('ERROR: --strategy proxy 需要 --proxy 参数', file=sys.stderr)
            return 1
        return strategy_proxy(args.output_dir, args.proxy)
    return 1


if __name__ == '__main__':
    raise SystemExit(main())