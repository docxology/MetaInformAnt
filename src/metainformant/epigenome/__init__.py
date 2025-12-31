"""Epigenomic analysis module for METAINFORMANT.

This module provides comprehensive tools for epigenomic data analysis,
including DNA methylation, ChIP-seq, ATAC-seq, chromatin accessibility,
and epigenetic track visualization.
"""

from __future__ import annotations

# Import all epigenome submodules
from . import (
    atac,
    chipseq,
    methylation,
    tracks,
    workflow,
)

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
    "atac",

    # Data visualization and workflow
    "tracks",
    "workflow",
]


