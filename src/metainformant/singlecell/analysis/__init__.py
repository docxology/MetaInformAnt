"""Single-cell analysis sub-package (clustering, dimensionality reduction, trajectory)."""

from __future__ import annotations

from . import clustering, nonlinear_methods, pca_methods, trajectory

__all__ = [
    "clustering",
    "nonlinear_methods",
    "pca_methods",
    "trajectory",
]
