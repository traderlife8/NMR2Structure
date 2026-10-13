"""MestReNova .mnova 文档 JSON 导出物读取与主流程桥接.

数据流:
  MestReNova.exe --nogui -w <file.mnova> --sf mnovaBatchExportJson
      → %USERPROFILE%/mnova_export/<doc>/<doc>__pageN__itemK__<nucleus[_2D]>.json
      → 本模块读取 JSON → BrukerSpectrum (1D) / Spectrum2D (2D)
      → 复用主流程 _build_problem_from_bruker / bridge_nmr2struct / correlation_2d

设计约束:
  - 只读: 不改写 .mnova 与导出 JSON
  - 错误显式暴露: MnovaJSONError 带 stage 字段 (遵循项目错误约定)
  - DEPT 判定复用 io_bruker._detect_pulse_role 的语义 (nucleus='13C-DEPT')
  - 2D 实验类型由核素对 + title 判定: 13C/1H → HSQC 或 HMBC (title/pulse 辅助), 1H/1H → COSY 或 ROESY
    MNova JSON 的 title 通常含原始数据文件名/脉冲程序, 优先用 title 关键词.
"""
from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np

from .io_bruker import BrukerSpectrum, IntegralRegion, Peak

# 2D 数据结构直接复用主流程定义, 避免重复实现
from .io_bruker_2d import Peak2D, Spectrum2D


class MnovaJSONError(Exception):
    """mnova JSON 读取/桥接错误."""

    def __init__(self, msg: str, path: str = '', stage: str = ''):
        super().__init__(f'[{stage}] {msg} (path={path})')
        self.path = path
        self.stage = stage


# MNova 导出 JSON 的 2D 核素对 → 实验类型判定.
# 优先用 title/pulse 关键词, 核素对仅作兜底 (HSQC 与 HMBC 核素对相同).
_2D_TYPE_KEYWORDS = (
    ('hmbc', 'HMBC'),
    ('hsqc', 'HSQC'),
    ('hmqc', 'HSQC'),
    ('cosy', 'COSY'),
    ('tocsy', 'TOCSY'),
    ('noesy', 'NOESY'),
    ('roesy', 'ROESY'),
)


def _detect_2d_type(doc: dict, nuc1: str, nuc2: str) -> str:
    """从 MNova experimentType (优先) / title 关键词 / 核素对 兜底 推断 2D 实验类型."""
    et = (doc.get('experimentType') or '').upper()
    # MNova 命名: '2D-HSQC' / '2D-HMBC' / '2D-COSY' / '2D-ROESY' ...
    for known in ('HMBC', 'HSQC', 'HMQC', 'COSY', 'TOCSY', 'NOESY', 'ROESY'):
        if known in et:
            return 'HSQC' if known == 'HMQC' else known
    t = (doc.get('title') or '').lower()
    for kw, etype in _2D_TYPE_KEYWORDS:
        if kw in t:
            return etype
    # 兜底: 核素对
    if nuc1 == '13C' and nuc2 == '1H':
        # 13C-1H 异核相关: 无法区分 HSQC/HMBC, 默认 HSQC (更常见)
        return 'HSQC'
    if nuc1 == '1H' and nuc2 == '1H':
        return 'COSY'
    return ''


def _parse_nucleus(raw: str) -> str:
    """MNova 核素标签 (如 '1H'/'13C') 直接可用; 防御性清洗."""
    n = (raw or '').strip()
    return n


def _pulse_role_from_nucleus_and_title(nucleus: str, doc: dict) -> str:
    """仿 io_bruker._detect_pulse_role: 返回 '1H'/'13C'/'13C-DEPT'."""
    nuc = nucleus
    # MNova experimentType 优先: '1D-DEPT-90'/'1D-DEPT-135' 明确
    et = (doc.get('experimentType') or '').upper()
    if 'DEPT' in et:
        return '13C-DEPT'
    t = (doc.get('title') or '').lower()
    if nuc == '13C':
        if 'dept' in t:
            return '13C-DEPT'
        return '13C'
    if nuc == '1H':
        return '1H'
    # 未知核素: 由 title 判断
    if 'dept' in t:
        return '13C-DEPT'
    if '13c' in t:
        return '13C'
    return nuc or '1H'


def _build_1d_spectrum(doc: dict, json_path: str) -> BrukerSpectrum:
    """MNova 1D JSON → BrukerSpectrum."""
    nuc_raw = _parse_nucleus(doc.get('nucleus', ''))
    role = _pulse_role_from_nucleus_and_title(nuc_raw, doc)
    nucleus = '13C-DEPT' if role == '13C-DEPT' else (nuc_raw or role)

    peaks: List[Peak] = []
    for pk in doc.get('peaks', []):
        ppm = pk.get('ppm1')
        if ppm is None:
            continue
        peaks.append(Peak(
            shift=float(ppm),
            intensity=float(pk.get('intensity', 0.0) or 0.0),
        ))
    if not peaks:
        raise MnovaJSONError(
            f'1D JSON 无峰 (nucleus={nuc_raw})', json_path, stage='peaks-empty')

    # 溶剂归一化: MNova 一般给 'CDCl3'/'Pyr', 与项目溶剂表一致, 只去空白
    solvent = (doc.get('solvent') or '').strip()

    # firstXppm/lastXppm → 积分区不可得, 留空 (主流程对积分缺失有兜底)
    integrals: List[IntegralRegion] = []

    return BrukerSpectrum(
        exp_no=str(os.path.basename(json_path)),  # 用文件名充当实验编号
        name=os.path.splitext(os.path.basename(json_path))[0],
        solvent=solvent,
        nucleus=nucleus,
        field_mhz=float(doc.get('frequency_MHz', 0.0) or 0.0),
        temperature_k=None,
        pulse_program='',  # MNova JSON 不导出脉冲程序; 2D 类型另由 title 判定
        peaks=peaks,
        integrals=integrals,
        source_path=os.path.abspath(json_path),
        parm_path=None,
        peaklist_path=None,
    )


def _build_2d_spectrum(doc: dict, json_path: str) -> Spectrum2D:
    """MNova 2D JSON → Spectrum2D.

    MNova 导出只含峰表 (无矩阵), data 置空数组; 主流程 correlation_2d
    只消费 peaks, 不依赖矩阵本身. f1/f2 ppm 轴用线性近似 (两端点 + 点数).
    """
    nuc1 = _parse_nucleus(doc.get('nucleus', ''))
    nuc2 = _parse_nucleus(doc.get('nucleus2', ''))
    title = doc.get('title', '')
    etype = _detect_2d_type(doc, nuc1, nuc2)
    if not etype:
        raise MnovaJSONError(
            f'无法判定 2D 实验类型 (nuc1={nuc1} nuc2={nuc2} title={title!r})',
            json_path, stage='2d-type')

    f1_scale = np.linspace(doc['firstXppm'], doc['lastXppm'], doc['pointCount'])

    peaks = [
        Peak2D(
            f1_ppm=float(pk['ppm1']),
            f2_ppm=float(pk.get('ppm2', 0.0)),
            intensity=float(pk.get('intensity', 0.0) or 0.0),
        )
        for pk in doc.get('peaks', [])
        if pk.get('ppm1') is not None
    ]
    if not peaks:
        raise MnovaJSONError('2D JSON 无峰', json_path, stage='peaks-empty')

    return Spectrum2D(
        exp_no=str(os.path.basename(json_path)),
        experiment_type=etype,
        f1_nucleus=nuc1,
        f2_nucleus=nuc2,
        f1_ppm_scale=f1_scale,
        f2_ppm_scale=np.linspace(0.0, 1.0, 2),  # 占位: 主流程仅用 peaks
        data=np.zeros((0, 0)),
        peaks=peaks,
        solvent=(doc.get('solvent') or '').strip(),
        field_mhz=float(doc.get('frequency_MHz', 0.0) or 0.0),
        pulse_program='',
        source_path=os.path.abspath(json_path),
        name=os.path.splitext(os.path.basename(json_path))[0],
    )


def read_mnova_json_folder(json_dir: str) -> tuple[List[BrukerSpectrum], Dict[str, Spectrum2D]]:
    """读取 mnovaBatchExportJson 导出目录, 返回 (specs_1d, spectra_2d).

    json_dir: 含 '<doc>__pageN__itemK__<nucleus[_2D]>.json' 的目录
    返回:
      specs_1d: List[BrukerSpectrum] (1H / 13C / 13C-DEPT)
      spectra_2d: {etype: Spectrum2D} 同类型取峰数最多的一条
    抛出: MnovaJSONError (目录不存在 / 无有效 JSON)
    """
    json_dir = os.path.abspath(json_dir)
    if not os.path.isdir(json_dir):
        raise MnovaJSONError(f'导出目录不存在', json_dir, stage='folder')

    specs: List[BrukerSpectrum] = []
    candidates_2d: Dict[str, List[Spectrum2D]] = {}
    n_json = 0
    for fn in sorted(os.listdir(json_dir)):
        if not fn.endswith('.json') or fn.startswith('_'):
            continue
        path = os.path.join(json_dir, fn)
        try:
            with open(path, 'r', encoding='utf-8') as f:
                doc = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            raise MnovaJSONError(f'JSON 解析失败: {e}', path, stage='json-parse')
        n_json += 1
        dim = doc.get('dim', 1)
        if dim == 1:
            specs.append(_build_1d_spectrum(doc, path))
        elif dim == 2:
            spec2d = _build_2d_spectrum(doc, path)
            candidates_2d.setdefault(spec2d.experiment_type, []).append(spec2d)
        # 其他维度跳过 (当前不支持 3D)
    if n_json == 0:
        raise MnovaJSONError('目录中无任何 .json 导出文件', json_dir, stage='folder-empty')

    # 同类型 2D 取峰数最多的一条 (与 read_bruker_2d_folder "同类型取第一个" 语义一致但更稳)
    spectra_2d: Dict[str, Spectrum2D] = {}
    for etype, lst in candidates_2d.items():
        spectra_2d[etype] = max(lst, key=lambda s: len(s.peaks))

    return specs, spectra_2d


def export_mnova_json(mnova_path: str, quiet: bool = False) -> str:
    """调用 MestReNova.exe 无界面导出 .mnova → JSON 目录.

    返回导出目录路径 (C:/Users/<user>/mnova_export/<doc_name>).
    抛出: MnovaJSONError (MestReNova.exe 不存在 / 导出失败)
    """
    import subprocess

    exe = r'D:\Program Files\Mestrelab Research S.L\MestReNova\MestReNova.exe'
    if not os.path.isfile(exe):
        raise MnovaJSONError('未找到 MestReNova.exe', exe, stage='mnova-exe')
    if not os.path.isfile(mnova_path):
        raise MnovaJSONError('.mnova 文件不存在', mnova_path, stage='mnova-file')

    # 脚本须已部署到 Mnova scripts 目录; 不存在时从项目 util 复制
    qs_deployed = r'D:\Program Files\Mestrelab Research S.L\MestReNova\scripts\mnovaBatchExportJson.qs'
    qs_source = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             'util', 'mnovaBatchExportJson.qs')
    if not os.path.isfile(qs_deployed):
        if not os.path.isfile(qs_source):
            raise MnovaJSONError(
                f'导出脚本缺失: {qs_deployed} 与 {qs_source} 均不存在',
                qs_deployed, stage='qs-missing')
        try:
            os.makedirs(os.path.dirname(qs_deployed), exist_ok=True)
            with open(qs_source, 'rb') as src, open(qs_deployed, 'wb') as dst:
                dst.write(src.read())
        except OSError as e:
            raise MnovaJSONError(
                f'导出脚本部署失败 (需管理员权限): {e}', qs_deployed, stage='qs-deploy')

    home = os.path.expanduser('~')
    doc_name = os.path.splitext(os.path.basename(mnova_path))[0]
    out_dir = os.path.join(home, 'mnova_export', doc_name)
    if os.path.isdir(out_dir):
        # 清理旧导出, 防止混入旧谱
        import shutil
        shutil.rmtree(out_dir, ignore_errors=True)

    cmd = [exe, '--nogui', '-w', os.path.abspath(mnova_path), '--sf', 'mnovaBatchExportJson']
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if proc.returncode != 0:
        raise MnovaJSONError(
            f'MestReNova 导出失败 exit={proc.returncode} stderr={proc.stderr[:500]}',
            mnova_path, stage='mnova-export')
    if not os.path.isdir(out_dir) or not any(
            fn.endswith('.json') for fn in os.listdir(out_dir)):
        raise MnovaJSONError('导出完成但未发现 JSON 产物', out_dir, stage='mnova-export-empty')
    if not quiet:
        n = sum(1 for fn in os.listdir(out_dir) if fn.endswith('.json') and not fn.startswith('_'))
        print(f'[mnova] 导出 {n} 条谱 JSON → {out_dir}', file=sys.stderr)
    return out_dir


def read_mnova_file(mnova_path: str, quiet: bool = False) -> tuple[List[BrukerSpectrum], Dict[str, Spectrum2D]]:
    """端到端: .mnova → (specs_1d, spectra_2d)."""
    json_dir = export_mnova_json(mnova_path, quiet=quiet)
    return read_mnova_json_folder(json_dir)
