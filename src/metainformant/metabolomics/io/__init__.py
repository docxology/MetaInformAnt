"""Metabolomics I/O sub-package.

File format readers and writers for metabolomics data (intensity-matrix CSV,
MGF mass spectra).
"""

from __future__ import annotations

from . import formats

__all__ = ["formats"]
