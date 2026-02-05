"""Epigenome analysis module for METAINFORMANT.

This module provides tools for analyzing epigenetic data, including
DNA methylation, histone modifications (ChIP-seq), and chromatin
accessibility (ATAC-seq), along with peak calling, chromatin state
learning, and genome browser track generation.
"""

from __future__ import annotations

# Import subpackages
from . import analysis, assays, chromatin_state, peak_calling, visualization, workflow
from .analysis import (
    tracks,
)
from .analysis.tracks import read_bedgraph

# Import modules from subpackages for backward compatibility
from .assays import (
    atacseq,
    chipseq,
)
from .assays.methylation import (
    compute_beta_values,
    load_cpg_table,
    summarize_beta_by_chromosome,
)

# Chromatin state public API
from .chromatin_state.state_learning import (
    assign_states,
    compare_chromatin_states,
    compute_state_enrichment,
    interpret_states,
    learn_chromatin_states,
    segment_genome,
)

# Peak calling public API
from .peak_calling.peak_detection import (
    call_peaks_broad,
    call_peaks_simple,
    compute_frip,
    differential_peaks,
    filter_peaks,
    merge_peaks,
)
from .visualization import visualization as visualization_module
from .workflow import workflow as workflow_module

# Optional imports with graceful fallbacks
try:
    from .assays import atacseq
except ImportError:
    atacseq = None

try:
    from .assays import chipseq
except ImportError:
    chipseq = None

try:
    from .assays import methylation
except ImportError:
    methylation = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "analysis",
    "assays",
    "chromatin_state",
    "peak_calling",
    "visualization",
    "workflow",
    # Epigenetic assays
    "atacseq",
    "chipseq",
    "methylation",
    "compute_beta_values",
    "load_cpg_table",
    "summarize_beta_by_chromosome",
    "read_bedgraph",
    # Peak calling
    "call_peaks_simple",
    "call_peaks_broad",
    "merge_peaks",
    "filter_peaks",
    "compute_frip",
    "differential_peaks",
    # Chromatin state learning
    "learn_chromatin_states",
    "assign_states",
    "interpret_states",
    "compute_state_enrichment",
    "segment_genome",
    "compare_chromatin_states",
    # Data visualization and tracks
    "tracks",
    "visualization_module",
    # Workflow
    "workflow_module",
]
