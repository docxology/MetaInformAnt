"""Effective population size calculations.

This module provides functions for calculating effective population size
under various demographic scenarios.
"""

from __future__ import annotations

from typing import List

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def effective_size_sex_ratio(male_size: float, female_size: float, sex_ratio: float = 0.5) -> float:
    """Calculate effective population size for unequal numbers of breeding males and females.

    Wright's (1931) formula Ne = 4 * Nm * Nf / (Nm + Nf). The sex ratio of the breeders is already
    fixed by the two counts (Nm / (Nm + Nf)), so this one expression covers every sex ratio: it
    equals the census size Nm + Nf when Nm == Nf and falls below it as the sexes diverge, and it is
    never larger than Nm + Nf.

    Args:
        male_size: Number of breeding males (Nm)
        female_size: Number of breeding females (Nf)
        sex_ratio: Proportion of the population that is male, strictly between 0 and 1 (default
            0.5). Accepted and validated for compatibility, but it does not enter the formula:
            the breeding sex ratio is already given by ``male_size`` and ``female_size``, so every
            value returns the same Ne

    Returns:
        Effective population size

    Raises:
        ValueError: If a population size is not positive or ``sex_ratio`` is not strictly between 0 and 1
    """
    if male_size <= 0 or female_size <= 0:
        raise ValueError("Population sizes must be positive")
    if not (0 < sex_ratio < 1):
        raise ValueError("Sex ratio must be between 0 and 1")

    return 4 * male_size * female_size / (male_size + female_size)


def harmonic_mean_effective_size(population_sizes: List[float]) -> float:
    """Calculate harmonic mean effective population size over multiple generations.

    Args:
        population_sizes: List of population sizes for each generation

    Returns:
        Harmonic mean effective population size

    Examples:
        >>> sizes = [100, 50, 200]
        >>> Ne = harmonic_mean_effective_size(sizes)
        >>> print(f"Harmonic mean Ne: {Ne:.1f}")
    """
    if not population_sizes:
        raise ValueError("Population sizes list cannot be empty")

    if any(size <= 0 for size in population_sizes):
        raise ValueError("All population sizes must be positive")

    n_generations = len(population_sizes)

    # Harmonic mean: n / sum(1/x_i)
    harmonic_sum = sum(1.0 / size for size in population_sizes)

    if harmonic_sum == 0:
        return float("inf")

    return n_generations / harmonic_sum


def effective_population_size_from_heterozygosity(heterozygosity: float, mutation_rate: float = 1e-8) -> float:
    """Estimate effective population size from heterozygosity.

    Uses the relationship H = 4*Ne*mu / (1 + 4*Ne*mu) for diploid organisms,
    rearranged to solve for Ne.

    Args:
        heterozygosity: Observed heterozygosity (0 < H < 1)
        mutation_rate: Mutation rate per site per generation (default 1e-8)

    Returns:
        Estimated effective population size

    Raises:
        ValueError: If heterozygosity is not in valid range (0, 1) or mutation rate <= 0
    """
    if heterozygosity <= 0 or heterozygosity >= 1:
        raise ValueError(f"Heterozygosity must be between 0 and 1 (exclusive), got {heterozygosity}")

    if mutation_rate <= 0:
        raise ValueError(f"Mutation rate must be positive, got {mutation_rate}")

    # From H = 4*Ne*mu / (1 + 4*Ne*mu), solve for Ne:
    # H * (1 + 4*Ne*mu) = 4*Ne*mu
    # H + 4*H*Ne*mu = 4*Ne*mu
    # H = 4*Ne*mu - 4*H*Ne*mu
    # H = 4*Ne*mu * (1 - H)
    # Ne = H / (4*mu * (1 - H))

    Ne = heterozygosity / (4 * mutation_rate * (1 - heterozygosity))
    return Ne
