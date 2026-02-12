"""QC reporting subpackage.

Provides comprehensive quality control report generation, sample
aggregation, threshold checking, and trend analysis."""
from __future__ import annotations

from . import multiqc_integration

__all__ = ['multiqc_integration']
