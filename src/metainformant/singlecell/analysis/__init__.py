"""Single-cell analysis sub-package (clustering, dimensionality, trajectory)."""
from __future__ import annotations

from . import clustering, dimensionality, nonlinear_methods, pca_methods, trajectory

__all__ = [
    "clustering",
    "dimensionality",
    "nonlinear_methods",
    "pca_methods",
    "trajectory",
]
