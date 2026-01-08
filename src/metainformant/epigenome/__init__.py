"""Epigenome analysis module for METAINFORMANT.

This module provides tools for analyzing epigenetic data, including
DNA methylation, histone modifications (ChIP-seq), and chromatin
accessibility (ATAC-seq), along with genome browser track generation.
"""

from __future__ import annotations

# Import subpackages
from . import analysis
from . import assays
from . import visualization
from . import workflow

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
from .analysis.tracks import read_bedgraph
from .analysis import (
    tracks,
)
from .workflow import workflow as workflow_module
from .visualization import visualization as visualization_module

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

    # Data visualization and tracks
    "tracks",
    "visualization_module",
    
    # Workflow
    "workflow_module",
]
