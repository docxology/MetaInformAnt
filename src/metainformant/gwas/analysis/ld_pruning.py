"""Linkage disequilibrium (LD) pruning for GWAS.

This module provides functions for LD-based variant pruning using a sliding
window approach, commonly used before PCA computation to remove correlated SNPs.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Import numpy with graceful fallback
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore


def ld_prune(
    genotype_matrix: List[List[int]],
    variant_positions: Optional[List[int]] = None,
    variant_chroms: Optional[List[int]] = None,
    window_size: int = 50,
    step_size: int = 5,
    r2_threshold: float = 0.2,
) -> List[int]:
    """Prune variants in LD using a sliding window approach.

    Within each window, computes pairwise R-squared between variants and greedily
    removes one variant from each correlated pair (preferring to remove the variant
    with higher missingness).

    Args:
        genotype_matrix: Genotype matrix (variants x samples), values 0/1/2/-1
        variant_positions: Optional list of variant positions (for window-based pruning)
        variant_chroms: Optional list of chromosome assignments per variant
        window_size: Number of variants in each window
        step_size: Number of variants to advance the window
        r2_threshold: R-squared threshold above which variants are pruned

    Returns:
        List of kept variant indices (0-based)
    """
    if not genotype_matrix:
        return []

    n_variants = len(genotype_matrix)
    n_samples = len(genotype_matrix[0]) if genotype_matrix else 0

    if n_variants == 0 or n_samples == 0:
        return []

    logger.info(f"LD pruning: {n_variants} variants, window={window_size}, " f"step={step_size}, r2={r2_threshold}")

    # Track which variants are removed
    removed = set()

    # Precompute missingness for each variant
    missingness = []
    for variant_gts in genotype_matrix:
        n_missing = sum(1 for g in variant_gts if g < 0)
        missingness.append(n_missing / n_samples if n_samples > 0 else 0.0)

    # Slide window across variants
    start = 0
    while start < n_variants:
        end = min(start + window_size, n_variants)

        # Get indices of non-removed variants in this window
        window_indices = [i for i in range(start, end) if i not in removed]

        # Check all pairs in the window
        for idx_a in range(len(window_indices)):
            i = window_indices[idx_a]
            if i in removed:
                continue

            for idx_b in range(idx_a + 1, len(window_indices)):
                j = window_indices[idx_b]
                if j in removed:
                    continue

                # Skip pairs on different chromosomes if chromosome info provided
                if variant_chroms is not None and variant_chroms[i] != variant_chroms[j]:
                    continue

                r2 = _compute_r_squared_pair(genotype_matrix[i], genotype_matrix[j])

                if r2 >= r2_threshold:
                    # Remove the variant with higher missingness
                    if missingness[i] >= missingness[j]:
                        removed.add(i)
                        break  # i is removed, move to next i
                    else:
                        removed.add(j)

        start += step_size

    kept = sorted(i for i in range(n_variants) if i not in removed)
    logger.info(f"LD pruning: kept {len(kept)}/{n_variants} variants")
    return kept


def _compute_r_squared_pair(geno_a: List[int], geno_b: List[int]) -> float:
    """Compute R-squared (squared Pearson correlation) between two genotype vectors.

    Delegates to the shared implementation in analysis.utils to avoid duplication.

    Args:
        geno_a: Genotype values for variant A (0, 1, 2 or -1 for missing)
        geno_b: Genotype values for variant B

    Returns:
        R-squared value in [0, 1]
    """
    from metainformant.gwas.analysis.utils import compute_r_squared

    return compute_r_squared(geno_a, geno_b)
