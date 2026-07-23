"""NMR 反推结构核心包.

设计原则：
- 数据原始保留：所有读取只读，禁止修改原始 Bruker/MNova 文件
- 确定性逻辑：峰匹配、打分、排序均由可复核代码控制
- 证据链：每个候选输出 matched/unmatched/conflict 证据
- 失败显式化：路径不存在、编码错、字段缺失、值域非法分别抛明确异常
- 多引擎支持：候选库 (nmrshiftdb2) + ML (NMR2Struct) 通过 bridge 串接
"""
from __future__ import annotations

from .bridge_nmr2struct import convert_bruker_folder_to_nmr2struct
from .constraints import Constraints, parse_constraints_from_acqus
from .drawing import (
    draw_single_mol,
    draw_topn_grid,
    draw_unmatched_overlay,
    write_unmatched_summary_csv,
)
from .evidence import (
    CandidateEvidence,
    ScoredCandidate,
    build_evidence,
    score_candidates,
)
from .io_bruker import (
    BrukerSpectrum,
    read_bruker_experiment,
    read_bruker_folder,
)
from .io_bruker_2d import (
    Peak2D,
    Spectrum2D,
    read_bruker_2d_experiment,
    read_bruker_2d_folder,
)
from .correlation_2d import (
    AssignmentRow,
    Correlation2DScore,
    FullAssignmentTable,
    build_full_assignment_table,
    score_candidate_2d,
    write_assignment_table_csv,
)
from .interpret_2d import (
    CHPair,
    HMBCConnection,
    Interpretation2D,
    OverlapRegion,
    SpinSystem,
    interpret_2d_sequential,
    write_step_files,
)
from .io_sdf import (
    NmrShiftDb2Entry,
    read_nmrshiftdb2_sdf,
    stream_nmrshiftdb2_sdf,
)
from .parse_nmr2struct_output import (
    parse_predictions_h5,
    write_candidates_markdown,
)
from .report import write_report_json, write_report_markdown
from .retrieve import RetrievalHit, retrieve_candidates
from .run_nmr2struct import build_inference_config, run_nmr2struct_inference

__all__ = [
    "BrukerSpectrum",
    "CandidateEvidence",
    "Constraints",
    "Correlation2DScore",
    "NmrShiftDb2Entry",
    "Peak2D",
    "RetrievalHit",
    "ScoredCandidate",
    "Spectrum2D",
    "build_evidence",
    "build_inference_config",
    "convert_bruker_folder_to_nmr2struct",
    "draw_single_mol",
    "draw_topn_grid",
    "draw_unmatched_overlay",
    "parse_constraints_from_acqus",
    "parse_predictions_h5",
    "read_bruker_2d_experiment",
    "read_bruker_2d_folder",
    "read_bruker_experiment",
    "read_bruker_folder",
    "read_nmrshiftdb2_sdf",
    "retrieve_candidates",
    "run_nmr2struct_inference",
    "score_candidate_2d",
    "score_candidates",
    "stream_nmrshiftdb2_sdf",
    "write_candidates_markdown",
    "write_report_json",
    "write_report_markdown",
    "write_unmatched_summary_csv",
]