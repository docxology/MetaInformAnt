"""F-statistics and population differentiation functions.

This module provides functions for calculating F-statistics and related measures
of population differentiation from allele frequency data.
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np

from metainformant.core.utils import logging

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

    The estimator pools the between-population variance of allele
    frequencies across loci:

        F_ST = sum_j var_p_j / sum_j (var_p_j + mean_i p_ij (1 - p_ij))

    For two populations this is algebraically identical to the standard
    heterozygosity form (Ht - Hs) / Ht, so the single-locus and
    multi-locus code paths agree for one locus.

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
            raise ValueError(
                "Single-list mode requires exactly 2 allele frequencies [p1, p2]"
            )
        for p in pop1_freqs:
            if not 0 <= p <= 1:
                raise ValueError(f"Invalid frequency values: {list(pop1_freqs)}")
        # One locus, two populations: the same estimator as the multi-locus path.
        return fst_from_allele_freq_matrix([[pop1_freqs[0]], [pop1_freqs[1]]])

    # Multi-locus mode
    if len(pop1_freqs) != len(pop2_freqs):
        raise ValueError("Population frequency arrays must have same length")

    if not pop1_freqs:
        raise ValueError("Frequency arrays cannot be empty")

    # Validate frequency values
    for i, (f1, f2) in enumerate(zip(pop1_freqs, pop2_freqs)):
        if not (0 <= f1 <= 1) or not (0 <= f2 <= 1):
            raise ValueError(f"Invalid frequency values at locus {i}: {f1}, {f2}")

    return fst_from_allele_freq_matrix([list(pop1_freqs), list(pop2_freqs)])


def fst_from_allele_freq_matrix(pop_freqs: List[List[float]]) -> float:
    """Calculate the multi-locus moment F_ST across any number of populations.

    Standard moment estimator pooled over loci:

        F_ST = sum_j var_p_j / sum_j (var_p_j + mean_i p_ij (1 - p_ij))

    where var_p_j is the variance of allele frequencies across populations
    at locus j. For more than two populations var_p uses the unbiased
    k/(k-1) corrected sample variance; for exactly two populations it is
    the population variance, which reproduces the classic
    (Ht - Hs) / Ht estimator and keeps the single-locus and multi-locus
    code paths consistent.

    Args:
        pop_freqs: One per-locus frequency list per population; all
            populations must have the same number of loci.

    Returns:
        F_ST value between 0 and 1; 0.0 when no variance can be estimated.

    Raises:
        ValueError: If fewer than 2 populations are supplied, the locus
            lists are empty or mismatched, or frequencies are invalid.

    Examples:
        >>> pops = [[0.6, 0.4], [0.3, 0.7], [0.9, 0.1]]
        >>> fst = fst_from_allele_freq_matrix(pops)
    """
    n_pops = len(pop_freqs)
    if n_pops < 2:
        raise ValueError("Need at least 2 populations")

    n_loci = len(pop_freqs[0])
    if n_loci == 0:
        raise ValueError("Frequency arrays cannot be empty")

    for i, pop in enumerate(pop_freqs):
        if len(pop) != n_loci:
            raise ValueError(f"Population {i} has {len(pop)} loci, expected {n_loci}")
        for j, p in enumerate(pop):
            if not 0 <= p <= 1:
                raise ValueError(
                    f"Invalid frequency value at population {i}, locus {j}: {p}"
                )

    variance_sum = 0.0
    total_sum = 0.0
    for j in range(n_loci):
        locus_freqs = [pop[j] for pop in pop_freqs]
        p_bar = sum(locus_freqs) / n_pops
        if n_pops > 2:
            var_p = sum((p - p_bar) ** 2 for p in locus_freqs) / (n_pops - 1)
        else:
            var_p = sum((p - p_bar) ** 2 for p in locus_freqs) / n_pops
        within = sum(p * (1.0 - p) for p in locus_freqs) / n_pops
        variance_sum += var_p
        total_sum += var_p + within

    if total_sum <= 0.0:
        return 0.0

    # Ensure F_ST is within valid range
    return max(0.0, min(1.0, variance_sum / total_sum))


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


def weirs_fst(population_counts: Dict[str, Dict[str, int]]) -> float:
    """Calculate Weir & Cockerham's (1984) F_ST from per-population haplotype counts.

    Implements the single-locus Weir & Cockerham (1984) variance-component
    estimator. For every distinct haplotype ``u`` the components are

        a_u: among-population variance component
        b_u: within-population sampling component
        c_u: within-individual component

    computed from per-population gene-copy counts with the observed
    heterozygosity term set to zero: haplotype counts carry no diploid
    genotype information, so ``c`` is 0 and ``b`` captures the binomial
    sampling variance of allele frequencies within populations. The pooled
    estimator is

        theta = sum_u a_u / sum_u (a_u + b_u + c_u)

    clamped to [0, 1] (the unbiased components can be slightly negative
    for homogeneous samples).

    Args:
        population_counts: Mapping of population label to haplotype-count
            mapping. Counts are gene copies (haplotypes), not diploid individuals.

    Returns:
        F_ST value between 0 and 1; 0.0 when it cannot be estimated.

    Raises:
        ValueError: If any haplotype count is negative.

    Examples:
        >>> counts = {"pop1": {"AT": 10, "AG": 5}, "pop2": {"AT": 4, "AG": 11}}
        >>> fst = weirs_fst(counts)
    """
    if not population_counts:
        return 0.0

    # Sorted iteration keeps the estimate independent of dict ordering.
    populations = sorted(population_counts)
    if any(
        count < 0 for pop in populations for count in population_counts[pop].values()
    ):
        raise ValueError("Haplotype counts must be non-negative")
    n_i = {pop: sum(population_counts[pop].values()) for pop in populations}
    # Populations with no sampled gene copies cannot inform the estimator.
    populations = [pop for pop in populations if n_i[pop] > 0]
    n_i = {pop: n_i[pop] for pop in populations}
    r = len(populations)

    if r < 2:
        return 0.0  # Need at least 2 populations

    total_n = sum(n_i.values())
    n_bar = total_n / r
    if n_bar <= 1.0:
        return 0.0  # Need at least 2 gene copies per population on average

    # Weir & Cockerham (1984) sample-size terms.
    sum_n_squared = sum(n * n for n in n_i.values())
    n_c = (total_n - sum_n_squared / total_n) / (r - 1)
    if n_c <= 0.0:
        return 0.0

    alleles = sorted(
        {haplotype for pop in populations for haplotype in population_counts[pop]}
    )

    a_total = 0.0
    b_total = 0.0
    c_total = 0.0
    for allele in alleles:
        p_hat = {
            pop: population_counts[pop].get(allele, 0) / n_i[pop] for pop in populations
        }
        p_bar = sum(n_i[pop] * p_hat[pop] for pop in populations) / total_n
        s_squared = sum(n_i[pop] * (p_hat[pop] - p_bar) ** 2 for pop in populations) / (
            (r - 1) * n_bar
        )
        het = p_bar * (1.0 - p_bar)

        a = (n_bar / n_c) * (
            s_squared - (1.0 / (n_bar - 1.0)) * (het - ((r - 1.0) / r) * s_squared)
        )
        b = (n_bar / (n_bar - 1.0)) * (het - ((r - 1.0) / r) * s_squared)
        c = 0.0

        a_total += a
        b_total += b
        c_total += c

    denominator = a_total + b_total + c_total
    if denominator <= 0.0:
        return 0.0

    return max(0.0, min(1.0, a_total / denominator))


def fst_confidence_interval(
    fst_value: float, sample_size: int, confidence_level: float = 0.95
) -> Tuple[float, float]:
    """Calculate confidence interval for F_ST estimate.

    Uses bootstrap resampling to estimate confidence intervals.

    Args:
        fst_value: Point estimate of F_ST
        sample_size: Sample size used for estimation
        confidence_level: Confidence level (default 0.95)

    Returns:
        Tuple of (lower_bound, upper_bound)
    """
    # Calculate standard error using Weir & Cockerham approximation
    # SE(F_ST) ≈ sqrt(F_ST * (1 - F_ST) / n) for large samples
    # For small samples, add correction factor

    if sample_size < 2:
        raise ValueError(
            "Sample size must be at least 2 for confidence interval calculation"
        )

    # Variance approximation for F_ST estimator
    # Based on asymptotic variance formula: Var(F_ST) ≈ 2*F_ST^2*(1-F_ST)^2 / n
    # for the case of two populations with equal sample sizes

    # Clamp F_ST to valid range to avoid math errors
    fst_clamped = max(0.001, min(0.999, fst_value))

    # Calculate variance using improved approximation
    # This uses the delta method approximation
    variance = (
        2 * fst_clamped * fst_clamped * (1 - fst_clamped) * (1 - fst_clamped)
    ) / sample_size

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
