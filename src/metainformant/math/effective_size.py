"""Effective population size calculations.

This module provides functions for calculating effective population size
under various demographic scenarios.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)


def effective_size_sex_ratio(male_size: float, female_size: float,
                           sex_ratio: float = 0.5) -> float:
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
        return float('inf')

    return n_generations / harmonic_sum
