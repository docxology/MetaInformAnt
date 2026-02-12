"""Single-cell analysis module for METAINFORMANT.

This module provides tools for single-cell data analysis, including
dimensionality reduction, clustering, trajectory inference, multi-sample
integration, cell type annotation, differential expression, and RNA
velocity estimation.
"""

from __future__ import annotations

from . import analysis, celltyping, data, differential, velocity, visualization

__all__ = ["analysis", "celltyping", "data", "differential", "velocity", "visualization"]
