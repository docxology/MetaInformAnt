"""Quality control analysis module for METAINFORMANT.

This module provides tools for assessing data quality, including
contamination detection, sequencing metrics, and FastQ file handling.
"""

from __future__ import annotations

# Import subpackages
from . import analysis
from . import io

# Import modules from subpackages for backward compatibility
from .analysis import (
    contamination,
    metrics,
)
from .io import (
    fastq,
)

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "analysis",
    "io",
    # Analysis
    "contamination",
    "metrics",
    # IO
    "fastq",
]
