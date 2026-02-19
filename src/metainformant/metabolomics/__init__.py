"""Metabolomics analysis module for METAINFORMANT.

This module provides tools for metabolite identification, mass spectrometry
data processing, pathway mapping, and metabolite-gene integration analysis.
"""

from __future__ import annotations

from . import analysis, io, pathways, visualization

__all__ = ["analysis", "io", "pathways", "visualization"]
