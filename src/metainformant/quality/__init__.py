"""Quality control analysis module for METAINFORMANT.

This module provides tools for assessing data quality, including
contamination detection, sequencing metrics, and FastQ file handling.
"""

from __future__ import annotations

from . import analysis, batch, io, reporting

__all__ = ["analysis", "batch", "io", "reporting"]
