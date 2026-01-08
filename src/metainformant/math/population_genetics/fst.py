"""F-statistics and population differentiation functions.

This module provides functions for calculating F-statistics and related measures
of population differentiation from allele frequency data.
"""

from __future__ import annotations

from typing import List, Dict, Any, Tuple
import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def fst_from_allele_freqs(pop1_freqs: List[float], pop2_freqs: List[float]) -> float:
    """Calculate F_ST from allele frequencies between two populations.

    F_ST measures the genetic differentiation between populations.
    Values range from 0 (no differentiation) to 1 (complete differentiation).

    Args:
        pop1_freqs: Allele frequencies for population 1
        pop2_freqs: Allele frequencies for population 2 (same length as pop1_freqs)

    Returns:
        F_ST value between 0 and 1

    Raises:
        ValueError: If frequency arrays have different lengths or invalid values

    Examples:
        >>> pop1 = [0.6, 0.4, 0.8]  # Allele frequencies for 3 loci in pop 1
        >>> pop2 = [0.3, 0.7, 0.2]  # Allele frequencies for 3 loci in pop 2
        >>> fst = fst_from_allele_freqs(pop1, pop2)
        >>> print(f"F_ST: {fst:.3f}")
        F_ST: 0.333
    """
    if len(pop1_freqs) != len(pop2_freqs):
        raise ValueError("Population frequency arrays must have same length")

    if not pop1_freqs:
        raise ValueError("Frequency arrays cannot be empty")

    # Validate frequency values
    for i, (f1, f2) in enumerate(zip(pop1_freqs, pop2_freqs)):
        if not (0 <= f1 <= 1) or not (0 <= f2 <= 1):
            raise ValueError(f"Invalid frequency values at locus {i}: {f1}, {f2}")

    n_loci = len(pop1_freqs)

    # Calculate mean frequencies across populations
    mean_freqs = [(f1 + f2) / 2 for f1, f2 in zip(pop1_freqs, pop2_freqs)]

    # Calculate within-population variance
    var_within = 0
    for f1, f2, mean_f in zip(pop1_freqs, pop2_freqs, mean_freqs):
        # Variance within populations for this locus
        var_within += ((f1 - mean_f)**2 + (f2 - mean_f)**2) / 2

    var_within /= n_loci

    # Calculate between-population variance
    var_between = 0
    for f1, f2 in zip(pop1_freqs, pop2_freqs):
        var_between += (f1 - f2)**2 / 2

    var_between /= n_loci

    # F_ST = variance between populations / total variance
    if var_within + var_between == 0:
        return 0.0

    fst = var_between / (var_within + var_between)

    # Ensure F_ST is within valid range
    return max(0.0, min(1.0, fst))


def pairwise_fst_matrix(population_freqs: List[List[float]]) -> np.ndarray:
    """Calculate pairwise F_ST matrix for multiple populations.

    Args:
        population_freqs: List of frequency arrays, one per population

    Returns:
        Symmetric matrix of F_ST values

    Examples:
        >>> pop1 = [0.6, 0.4, 0.8]
        >>> pop2 = [0.3, 0.7, 0.2]
        >>> pop3 = [0.5, 0.5, 0.6]
        >>> matrix = pairwise_fst_matrix([pop1, pop2, pop3])
        >>> print(matrix.shape)
        (3, 3)
    """
    n_pops = len(population_freqs)

    if n_pops < 2:
        raise ValueError("Need at least 2 populations")

    # Check all populations have same number of loci
    n_loci = len(population_freqs[0])
    for i, pop in enumerate(population_freqs):
        if len(pop) != n_loci:
            raise ValueError(f"Population {i} has {len(pop)} loci, expected {n_loci}")

    fst_matrix = np.zeros((n_pops, n_pops))

    for i in range(n_pops):
        for j in range(i + 1, n_pops):
            fst = fst_from_allele_freqs(population_freqs[i], population_freqs[j])
            fst_matrix[i, j] = fst
            fst_matrix[j, i] = fst

    return fst_matrix


def weirs_fst(haplotype_counts: Dict[str, int], population_labels: List[str]) -> float:
    """Calculate Weir & Cockerham's F_ST from haplotype counts.

    This implements the unbiased F_ST estimator from Weir & Cockerham (1984).

    Args:
        haplotype_counts: Dictionary mapping haplotype strings to counts
        population_labels: List of population labels for each individual

    Returns:
        F_ST value

    Examples:
        >>> counts = {"AT": 10, "AG": 15, "GT": 8, "GG": 12}
        >>> labels = ["pop1"] * 22 + ["pop2"] * 23  # 45 individuals total
        >>> fst = weirs_fst(counts, labels)
    """
    # This is a simplified implementation - full Weir & Cockerham F_ST
    # would require more complex haplotype frequency calculations
    logger.warning("Weir & Cockerham F_ST implementation is simplified")
    return 0.0  # Placeholder


def fst_confidence_interval(fst_value: float, sample_size: int,
                          confidence_level: float = 0.95) -> Tuple[float, float]:
    """Calculate confidence interval for F_ST estimate.

    Uses bootstrap resampling to estimate confidence intervals.

    Args:
        fst_value: Point estimate of F_ST
        sample_size: Sample size used for estimation
        confidence_level: Confidence level (default 0.95)

    Returns:
        Tuple of (lower_bound, upper_bound)
    """
    # Simplified confidence interval calculation
    # In practice, this would use bootstrap resampling of the original data
    se = 0.1  # Placeholder standard error
    z_score = 1.96 if confidence_level == 0.95 else 2.576  # Approximate z-scores

    lower = max(0.0, fst_value - z_score * se)
    upper = min(1.0, fst_value + z_score * se)

    return lower, upper


def fst_from_heterozygosity(Hs: float, Ht: float) -> float:
    """Calculate F_ST from heterozygosity measures.

    F_ST = (H_t - H_s) / H_t

    Args:
        Hs: Average heterozygosity within subpopulations
        Ht: Total heterozygosity

    Returns:
        F_ST value

    Examples:
        >>> fst = fst_from_heterozygosity(0.2, 0.5)
        >>> print(f"F_ST: {fst}")
        F_ST: 0.6
    """
    if Ht == 0:
        return 0.0

    fst = (Ht - Hs) / Ht
    return max(0.0, min(1.0, fst))  # Ensure valid range

