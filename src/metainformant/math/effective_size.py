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
