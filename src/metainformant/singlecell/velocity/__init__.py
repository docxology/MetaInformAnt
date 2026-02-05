"""RNA velocity analysis for single-cell data.

Provides RNA velocity estimation from spliced/unspliced count matrices,
velocity embedding projection, velocity-based pseudotime, confidence
metrics, and dynamical model fitting.
"""

from __future__ import annotations

from . import rna_velocity

__all__ = [
    "rna_velocity",
]
