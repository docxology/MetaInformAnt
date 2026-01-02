"""Population genetics statistical functions.

This module provides mathematical functions for population genetics statistics,
coalescent theory, and evolutionary computations.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def expected_pairwise_diversity(sample_size: int, theta: float) -> float:
    """Calculate expected pairwise nucleotide diversity (π).

    Args:
        sample_size: Number of sequences (n)
        theta: Population mutation parameter (4Neμ)

    Returns:
        Expected pairwise diversity
    """
    if sample_size < 2:
        return 0.0

    # E[π] = θ * (n/(n-1)) * sum(1/i for i in 1 to n-1)
    harmonic_sum = sum(1.0 / i for i in range(1, sample_size))
    return theta * (sample_size / (sample_size - 1)) * harmonic_sum


def expected_segregating_sites(sample_size: int, theta: float) -> float:
    """Calculate expected number of segregating sites.

    Args:
        sample_size: Number of sequences (n)
        theta: Population mutation parameter (4Neμ)

    Returns:
        Expected number of segregating sites
    """
    if sample_size < 2:
        return 0.0

    # E[S] = θ * sum(1/i for i in 1 to n-1)
    harmonic_sum = sum(1.0 / i for i in range(1, sample_size))
    return theta * harmonic_sum


def expected_coalescent_waiting_times(sample_size: int, Ne: float) -> List[float]:
    """Calculate expected waiting times for coalescent events.

    Args:
        sample_size: Starting number of lineages
        Ne: Effective population size

    Returns:
        List of expected waiting times for each coalescent interval
    """
    waiting_times = []
    current_lineages = sample_size

    while current_lineages > 1:
        # Expected waiting time for next coalescence
        # Time is in units of 4Ne generations
        rate = current_lineages * (current_lineages - 1) / 2
        expected_time = 1.0 / rate
        waiting_times.append(expected_time)
        current_lineages -= 1

    return waiting_times


def expected_r2_from_Ne_c(recombination_rate: float, Ne: float,
                         distance_bp: float) -> float:
    """Calculate expected LD (r²) from effective population size and recombination.

    Args:
        recombination_rate: Recombination rate per bp per generation
        Ne: Effective population size
        distance_bp: Distance between sites in base pairs

    Returns:
        Expected r² value
    """
    # r² decays as 1/(1 + 4Ne*c*d) where c is recombination rate, d is distance
    c = recombination_rate
    d = distance_bp

    denominator = 1 + 4 * Ne * c * d
    return 1.0 / denominator


def equilibrium_heterozygosity_infinite_alleles(theta: float) -> float:
    """Calculate equilibrium heterozygosity under infinite alleles model.

    Args:
        theta: Population mutation parameter (4Neμ)

    Returns:
        Equilibrium heterozygosity (H)
    """
    # H = θ / (1 + θ) for infinite alleles model
    return theta / (1 + theta)


def fixation_probability(selection_coefficient: float, population_size: int) -> float:
    """Calculate fixation probability under selection.

    Args:
        selection_coefficient: Selection coefficient (s)
        population_size: Population size (N)

    Returns:
        Probability of fixation
    """
    if selection_coefficient == 0:
        return 1.0 / (2 * population_size)  # Neutral case

    s = selection_coefficient
    N = population_size

    if s > 0:  # Beneficial mutation
        return (1 - math.exp(-2 * s)) / (1 - math.exp(-4 * N * s))
    else:  # Deleterious mutation
        return (1 - math.exp(-2 * s)) / (1 - math.exp(-4 * N * s))


def bottleneck_effective_size(initial_size: int, bottleneck_size: int,
                             bottleneck_duration: int, final_size: int) -> float:
    """Calculate effective population size after bottleneck.

    Args:
        initial_size: Initial population size
        bottleneck_size: Population size during bottleneck
        bottleneck_duration: Duration of bottleneck in generations
        final_size: Final population size

    Returns:
        Effective population size
    """
    # Simplified approximation
    # Ne = 1 / (1/N_initial + 1/N_bottleneck * duration + 1/N_final)
    Ne = 1.0 / (1.0/initial_size + bottleneck_duration * 1.0/bottleneck_size + 1.0/final_size)
    return Ne


def effective_size_from_family_size_variance(family_sizes: List[int]) -> float:
    """Estimate effective population size from variance in family sizes.

    Args:
        family_sizes: List of offspring counts per parent

    Returns:
        Estimated effective population size
    """
    if not family_sizes:
        return 0.0

    n = len(family_sizes)
    mean_family_size = sum(family_sizes) / n
    variance = sum((x - mean_family_size) ** 2 for x in family_sizes) / n

    # Ne = k / (V_k - 1) where k is mean family size, V_k is variance
    if variance > 1:
        return mean_family_size / (variance - 1)
    else:
        return float('inf')  # No variance means infinite Ne


def bootstrap_confidence_interval(values: List[float], confidence_level: float = 0.95,
                                n_bootstraps: int = 1000) -> Tuple[float, float]:
    """Calculate bootstrap confidence interval.

    Args:
        values: Original sample values
        confidence_level: Confidence level (0.95 for 95% CI)
        n_bootstraps: Number of bootstrap resamples

    Returns:
        Tuple of (lower_bound, upper_bound)
    """
    if not values:
        return (0.0, 0.0)

    np.random.seed(42)  # For reproducibility
    bootstrap_means = []

    for _ in range(n_bootstraps):
        # Resample with replacement
        sample = np.random.choice(values, size=len(values), replace=True)
        bootstrap_means.append(np.mean(sample))

    # Calculate confidence interval
    lower_percentile = (1 - confidence_level) / 2 * 100
    upper_percentile = (1 + confidence_level) / 2 * 100

    lower_bound = np.percentile(bootstrap_means, lower_percentile)
    upper_bound = np.percentile(bootstrap_means, upper_percentile)

    return (lower_bound, upper_bound)
