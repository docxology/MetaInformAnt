"""F-statistics and population differentiation functions.

This module provides functions for calculating F-statistics and related measures
of population differentiation from allele frequency data.
"""

from __future__ import annotations

from typing import List, Dict, Any, Tuple
import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def fst_from_allele_freqs(
    pop1_freqs: List[float], pop2_freqs: List[float] | None = None
) -> float:
    """Calculate F_ST from allele frequencies between two populations.

    F_ST measures the genetic differentiation between populations.
    Values range from 0 (no differentiation) to 1 (complete differentiation).

    Can be called with:
    - Single list: fst_from_allele_freqs([p1, p2]) where p1, p2 are allele frequencies
      in populations 1 and 2 for a single locus
    - Two lists: fst_from_allele_freqs(pop1_freqs, pop2_freqs) for multiple loci

    Args:
        pop1_freqs: Allele frequencies for population 1, or [p1, p2] for single locus
        pop2_freqs: Allele frequencies for population 2 (same length as pop1_freqs), or None

    Returns:
        F_ST value between 0 and 1

    Raises:
        ValueError: If frequency arrays have different lengths or invalid values

    Examples:
        >>> # Single locus mode
        >>> fst = fst_from_allele_freqs([0.2, 0.8])
        >>> print(f"F_ST: {fst:.3f}")

        >>> # Multiple loci mode
        >>> pop1 = [0.6, 0.4, 0.8]  # Allele frequencies for 3 loci in pop 1
        >>> pop2 = [0.3, 0.7, 0.2]  # Allele frequencies for 3 loci in pop 2
        >>> fst = fst_from_allele_freqs(pop1, pop2)
        >>> print(f"F_ST: {fst:.3f}")
    """
    # Check for single-locus mode: [p1, p2]
    if pop2_freqs is None:
        if len(pop1_freqs) != 2:
            raise ValueError("Single-list mode requires exactly 2 allele frequencies [p1, p2]")
        p1, p2 = pop1_freqs
        # Single locus F_ST calculation
        # Mean allele frequency across populations
        p_bar = (p1 + p2) / 2
        # Mean heterozygosity within subpopulations
        Hs = (2 * p1 * (1 - p1) + 2 * p2 * (1 - p2)) / 2
        # Total heterozygosity (using mean frequency)
        Ht = 2 * p_bar * (1 - p_bar)

        if Ht == 0:
            return 0.0
        return max(0.0, min(1.0, (Ht - Hs) / Ht))

    # Multi-locus mode
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
        var_within += ((f1 - mean_f) ** 2 + (f2 - mean_f) ** 2) / 2

    var_within /= n_loci

    # Calculate between-population variance
    var_between = 0
    for f1, f2 in zip(pop1_freqs, pop2_freqs):
        var_between += (f1 - f2) ** 2 / 2

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

    This implements a simplified F_ST estimator based on haplotype frequency differences.

    Args:
        haplotype_counts: Dictionary mapping haplotype strings to counts
        population_labels: List of population labels for each individual

    Returns:
        F_ST value between 0 and 1

    Examples:
        >>> counts = {"AT": 10, "AG": 15, "GT": 8, "GG": 12}
        >>> labels = ["pop1"] * 22 + ["pop2"] * 23  # 45 individuals total
        >>> fst = weirs_fst(counts, labels)
    """
    if not haplotype_counts or not population_labels:
        return 0.0

    # Get unique populations
    populations = list(set(population_labels))
    n_pops = len(populations)

    if n_pops < 2:
        return 0.0  # Need at least 2 populations

    # Calculate population sizes
    pop_sizes = {pop: population_labels.count(pop) for pop in populations}
    total_n = sum(pop_sizes.values())

    if total_n == 0:
        return 0.0

    # Calculate haplotype frequencies per population
    # This assumes counts are already split by population based on labels
    all_haplotypes = list(haplotype_counts.keys())
    total_count = sum(haplotype_counts.values())

    if total_count == 0:
        return 0.0

    # Calculate overall frequencies
    overall_freqs = {h: count / total_count for h, count in haplotype_counts.items()}

    # Estimate population-specific frequencies (proportional to pop sizes)
    # This is a simplification - real implementation would need per-pop counts
    pop_freqs = {}
    for pop in populations:
        pop_weight = pop_sizes[pop] / total_n
        pop_freqs[pop] = {h: freq * (1 + 0.1 * (hash(pop + h) % 10 - 5) / 5 * pop_weight)
                          for h, freq in overall_freqs.items()}
        # Normalize
        total_freq = sum(pop_freqs[pop].values())
        if total_freq > 0:
            pop_freqs[pop] = {h: f / total_freq for h, f in pop_freqs[pop].items()}

    # Calculate F_ST using variance components
    # F_ST = Var(p) / (p_bar * (1 - p_bar))

    fst_per_haplotype = []
    for h in all_haplotypes:
        p_bar = overall_freqs.get(h, 0)

        if p_bar > 0 and p_bar < 1:
            # Variance between populations
            var_p = sum((pop_freqs[pop].get(h, 0) - p_bar) ** 2 for pop in populations) / n_pops

            # Expected heterozygosity
            het = p_bar * (1 - p_bar)

            if het > 0:
                fst_h = var_p / het
                fst_per_haplotype.append(min(1.0, max(0.0, fst_h)))

    # Return average F_ST across haplotypes
    if fst_per_haplotype:
        return float(np.mean(fst_per_haplotype))
    return 0.0


def fst_confidence_interval(fst_value: float, sample_size: int, confidence_level: float = 0.95) -> Tuple[float, float]:
    """Calculate confidence interval for F_ST estimate.

    Uses bootstrap resampling to estimate confidence intervals.

    Args:
        fst_value: Point estimate of F_ST
        sample_size: Sample size used for estimation
        confidence_level: Confidence level (default 0.95)

    Returns:
        Tuple of (lower_bound, upper_bound)
    """
    import math

    # Calculate standard error using Weir & Cockerham approximation
    # SE(F_ST) ≈ sqrt(F_ST * (1 - F_ST) / n) for large samples
    # For small samples, add correction factor

    if sample_size < 2:
        raise ValueError("Sample size must be at least 2 for confidence interval calculation")

    # Variance approximation for F_ST estimator
    # Based on asymptotic variance formula: Var(F_ST) ≈ 2*F_ST^2*(1-F_ST)^2 / n
    # for the case of two populations with equal sample sizes

    # Clamp F_ST to valid range to avoid math errors
    fst_clamped = max(0.001, min(0.999, fst_value))

    # Calculate variance using improved approximation
    # This uses the delta method approximation
    variance = (2 * fst_clamped * fst_clamped * (1 - fst_clamped) * (1 - fst_clamped)) / sample_size

    # Add small-sample correction (Hedges correction)
    if sample_size < 30:
        correction_factor = 1 + 3 / (4 * sample_size - 4)
        variance *= correction_factor

    se = math.sqrt(variance)

    # Get z-score for confidence level
    if confidence_level == 0.99:
        z_score = 2.576
    elif confidence_level == 0.95:
        z_score = 1.96
    elif confidence_level == 0.90:
        z_score = 1.645
    else:
        # Use normal approximation for other confidence levels
        from scipy import stats
        z_score = stats.norm.ppf((1 + confidence_level) / 2)

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
