"""Effective population size calculations.

This module provides functions for calculating effective population size
under various demographic scenarios.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)


def effective_size_sex_ratio(male_size: float, female_size: float, sex_ratio: float = 0.5) -> float:
    """Calculate effective population size accounting for unequal sex ratios.

    Args:
        male_size: Number of males
        female_size: Number of females
        sex_ratio: Proportion of population that is male (default 0.5)

    Returns:
        Effective population size
    """
    if male_size <= 0 or female_size <= 0:
        raise ValueError("Population sizes must be positive")
    if not (0 < sex_ratio < 1):
        raise ValueError("Sex ratio must be between 0 and 1")

    total_size = male_size + female_size

    # Wright's formula for unequal sex ratios
    # Ne = 4 * Nm * Nf / (Nm + Nf)
    # Adjusted for sex ratio

    # If sex ratio is 0.5, this reduces to standard formula
    if abs(sex_ratio - 0.5) < 1e-10:
        return 4 * male_size * female_size / (male_size + female_size)

    # General case with unequal sex contribution
    # This is a simplified approximation
    harmonic_mean = 2 * male_size * female_size / (male_size + female_size)
    sex_ratio_adjustment = 1 / (sex_ratio * (1 - sex_ratio))

    return harmonic_mean * sex_ratio_adjustment


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
