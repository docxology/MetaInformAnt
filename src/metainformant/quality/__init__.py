"""Data quality control and assessment module for METAINFORMANT.

This module provides comprehensive quality control tools for biological data,
including FASTQ quality assessment, contamination detection, quality metrics,
and data validation across all supported data types.
"""

from __future__ import annotations

# Import all quality control submodules
from . import (
    contamination,
    fastq,
    metrics,
)

# Optional imports with graceful fallbacks
try:
    from . import contamination
except ImportError:
    contamination = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Quality assessment tools
    "fastq",
    "metrics",

    # Specialized analysis
    "contamination",
]


