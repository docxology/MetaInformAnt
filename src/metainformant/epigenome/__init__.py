"""Epigenomic analysis module for METAINFORMANT.

This module provides comprehensive tools for epigenomic data analysis,
including DNA methylation, ChIP-seq, ATAC-seq, chromatin accessibility,
and epigenetic track visualization.
"""

from __future__ import annotations

# Import all epigenome submodules
from . import (
    chipseq,
    methylation,
    tracks,
    visualization,
    workflow,
)

# Direct imports of commonly used functions
from .methylation import compute_beta_values, load_cpg_table, summarize_beta_by_chromosome
from .tracks import load_bedgraph_track as read_bedgraph

# Optional imports with graceful fallbacks
try:
    from . import workflow
except ImportError:
    workflow = None

try:
    from . import tracks
except ImportError:
    tracks = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core epigenomic assays
    "methylation",
    "chipseq",

    # Data visualization and workflow
    "tracks",
    "visualization",
    "workflow",
]



