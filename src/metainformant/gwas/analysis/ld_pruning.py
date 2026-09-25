"""Linkage disequilibrium (LD) pruning for GWAS.

This module provides functions for LD-based variant pruning using a sliding
window approach, commonly used before PCA computation to remove correlated SNPs.
"""

from __future__ import annotations

from typing import List, Optional

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Import numpy with graceful fallback
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None


def ld_prune(
    genotype_matrix: List[List[int]],
    variant_positions: Optional[List[int]] = None,
    variant_chroms: Optional[List[int]] = None,
    window_size: int = 50,
    step_size: int = 5,
    r2_threshold: float = 0.2,
) -> List[int]:
    """Prune variants in LD using a sliding window approach.

    Two window modes:

    * Variant-count mode (default, ``variant_positions=None``): a sliding
      window of ``window_size`` consecutive variants advanced by
      ``step_size`` variants; pairs are only formed inside a window.
    * Physical-distance mode (``variant_positions`` supplied; PLINK
      ``--indep-pairwise`` convention): ``window_size`` is interpreted as
      kilobases and positions are base pairs. A pair is only tested when
      both variants sit within ``window_size`` kb of each other (and on the
      same chromosome, when ``variant_chroms`` is given). ``step_size`` is
      not used in this mode.

    Within each window, computes pairwise R-squared between variants and greedily
    removes one variant from each correlated pair (preferring to remove the variant
    with higher missingness).

    Args:
        genotype_matrix: Genotype matrix (variants x samples), values 0/1/2/-1
        variant_positions: Optional list of variant positions in base pairs;
            switches on physical-distance mode
        variant_chroms: Optional list of chromosome assignments per variant
        window_size: Number of variants per window (variant-count mode) or
            window span in kilobases (physical-distance mode)
        step_size: Number of variants to advance the window (variant-count
            mode; ignored in physical-distance mode)
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

    physical_mode = (
        variant_positions is not None and len(variant_positions) == n_variants
    )

    logger.info(
        f"LD pruning: {n_variants} variants, window={window_size}, "
        f"step={step_size}, r2={r2_threshold}"
    )

    # Track which variants are removed
    removed = set()

    # Precompute missingness for each variant
    missingness = []
    for variant_gts in genotype_matrix:
        n_missing = sum(1 for g in variant_gts if g < 0)
        missingness.append(n_missing / n_samples if n_samples > 0 else 0.0)

    if physical_mode:
        positions = [float(pos) for pos in variant_positions]
        window_span_bp = float(window_size) * 1000.0
        # Visit variants in ascending positional order so the inner scan can
        # stop as soon as the physical window is exceeded.
        order = sorted(range(n_variants), key=lambda i: positions[i])
        for a in range(n_variants):
            i = order[a]
            if i in removed:
                continue
            for b in range(a + 1, n_variants):
                j = order[b]
                if j in removed:
                    continue
                if abs(positions[j] - positions[i]) > window_span_bp:
                    break  # sorted positions; all further pairs are farther apart
                if (
                    variant_chroms is not None
                    and variant_chroms[i] != variant_chroms[j]
                ):
                    continue
                r2 = _compute_r_squared_pair(genotype_matrix[i], genotype_matrix[j])
                if r2 >= r2_threshold:
                    # Remove the variant with higher missingness
                    if missingness[i] >= missingness[j]:
                        removed.add(i)
                        break  # i is removed, move to next i
                    else:
                        removed.add(j)
    else:
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
                    if (
                        variant_chroms is not None
                        and variant_chroms[i] != variant_chroms[j]
                    ):
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
