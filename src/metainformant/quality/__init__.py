"""Quality control analysis module for METAINFORMANT.

This module provides tools for assessing data quality, including
contamination detection, sequencing metrics, and FastQ file handling.
"""

from __future__ import annotations

# Import subpackages
from . import analysis
from . import io
from . import reporting

# Import modules from subpackages for backward compatibility
from .analysis import (
    contamination,
    metrics,
)
from .io import (
    fastq,
)
from .reporting import (
    generate_qc_report,
    aggregate_sample_qc,
    check_qc_thresholds,
    default_qc_thresholds,
    qc_trend_analysis,
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
    # Reporting
    "reporting",
    "generate_qc_report",
    "aggregate_sample_qc",
    "check_qc_thresholds",
    "default_qc_thresholds",
    "qc_trend_analysis",
]
