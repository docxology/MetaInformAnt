"""Alternative splicing analysis submodule for RNA-seq data.

This subpackage provides tools for detecting and quantifying alternative
splicing events from RNA-seq data, including:
- Splice junction detection and classification
- Percent Spliced In (PSI) computation
- Differential splicing analysis
- Novel junction identification
- Splice site strength scoring
- Isoform quantification via EM algorithm
- Splice graph construction and isoform enumeration
- Isoform diversity and usage comparison
"""

from __future__ import annotations

from .detection import (
    classify_splicing_events,
    compute_psi,
    compute_splice_site_strength,
    detect_splice_junctions,
    differential_splicing,
    find_novel_junctions,
)
from .isoforms import (
    build_isoform_graph,
    compare_isoform_usage,
    compute_isoform_diversity,
    enumerate_isoforms,
    quantify_isoforms,
)

__all__ = [
    # Detection
    "detect_splice_junctions",
    "classify_splicing_events",
    "compute_psi",
    "differential_splicing",
    "find_novel_junctions",
    "compute_splice_site_strength",
    # Isoforms
    "quantify_isoforms",
    "build_isoform_graph",
    "enumerate_isoforms",
    "compute_isoform_diversity",
    "compare_isoform_usage",
]
