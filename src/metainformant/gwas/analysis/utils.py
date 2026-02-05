"""Shared utility functions for GWAS analysis.

Provides common computations used across multiple analysis and visualization modules,
extracted to eliminate code duplication.
"""

from __future__ import annotations

import math
from typing import List

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore


def compute_r_squared(geno_a: List[int], geno_b: List[int]) -> float:
    """Compute R-squared (squared Pearson correlation) between two genotype vectors.

    Missing values (< 0) are excluded from the computation. When fewer than 3
    valid (non-missing) pairs remain, returns 0.0. When either vector has zero
    variance (constant), returns 0.0.

    Args:
        geno_a: Genotype values for variant A (0, 1, 2 or -1 for missing).
        geno_b: Genotype values for variant B (0, 1, 2 or -1 for missing).

    Returns:
        R-squared value in [0, 1].
    """
    # Filter to non-missing pairs
    pairs_a: List[float] = []
    pairs_b: List[float] = []
    for a, b in zip(geno_a, geno_b):
        if a >= 0 and b >= 0:
            pairs_a.append(float(a))
            pairs_b.append(float(b))

    n = len(pairs_a)
    if n < 3:
        return 0.0

    if HAS_NUMPY:
        arr_a = np.array(pairs_a)
        arr_b = np.array(pairs_b)
        if np.std(arr_a) < 1e-12 or np.std(arr_b) < 1e-12:
            return 0.0
        corr_matrix = np.corrcoef(arr_a, arr_b)
        r = float(corr_matrix[0, 1])
        if math.isnan(r):
            return 0.0
        return r * r

    # Pure Python fallback
    mean_a = sum(pairs_a) / n
    mean_b = sum(pairs_b) / n

    cov = sum((a - mean_a) * (b - mean_b) for a, b in zip(pairs_a, pairs_b))
    var_a = sum((a - mean_a) ** 2 for a in pairs_a)
    var_b = sum((b - mean_b) ** 2 for b in pairs_b)

    denom = math.sqrt(var_a * var_b)
    if denom < 1e-12:
        return 0.0

    r = cov / denom
    return r * r
