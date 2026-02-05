"""QC reporting subpackage.

Provides comprehensive quality control report generation, sample
aggregation, threshold checking, and trend analysis.
"""

from __future__ import annotations

from .multiqc_integration import (
    generate_qc_report,
    aggregate_sample_qc,
    check_qc_thresholds,
    default_qc_thresholds,
    qc_trend_analysis,
)

__all__ = [
    "generate_qc_report",
    "aggregate_sample_qc",
    "check_qc_thresholds",
    "default_qc_thresholds",
    "qc_trend_analysis",
]
