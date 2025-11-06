"""Epigenome domain functionality.

Public API:
- read_bedgraph: load bedGraph tracks into a DataFrame
- load_cpg_table / compute_beta_values / summarize_beta_by_chromosome: basic methylation utilities
- differential_methylation: identify differentially methylated sites
- mqtl_analysis: methylation quantitative trait locus analysis
- ChIP-seq: peak calling and analysis
- ATAC-seq: chromatin accessibility analysis
"""

from .atac import (
    calculate_accessibility_scores,
    compare_accessibility,
    identify_accessible_regions,
)
from .chipseq import (
    analyze_peak_overlap,
    call_peaks_simple,
    peak_enrichment_analysis,
)
from .methylation import (
    compute_beta_values,
    differential_methylation,
    load_cpg_table,
    mqtl_analysis,
    summarize_beta_by_chromosome,
)
from .tracks import read_bedgraph

try:
    from .workflow import (
        run_atacseq_workflow,
        run_chipseq_workflow,
        run_methylation_workflow,
    )
    _workflow_available = True
except ImportError:
    _workflow_available = False

__all__ = [
    # Track I/O
    "read_bedgraph",
    # Methylation
    "load_cpg_table",
    "compute_beta_values",
    "summarize_beta_by_chromosome",
    "differential_methylation",
    "mqtl_analysis",
    # ChIP-seq
    "call_peaks_simple",
    "analyze_peak_overlap",
    "peak_enrichment_analysis",
    # ATAC-seq
    "identify_accessible_regions",
    "calculate_accessibility_scores",
    "compare_accessibility",
]

if _workflow_available:
    __all__.extend([
        "run_methylation_workflow",
        "run_chipseq_workflow",
        "run_atacseq_workflow",
    ])
