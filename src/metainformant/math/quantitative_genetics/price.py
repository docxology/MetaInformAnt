"""Price equation and quantitative genetics functions.

This module provides functions for Price equation analysis and quantitative genetics.
"""

from __future__ import annotations

import statistics
from typing import List

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def price_equation(
    fitness: List[float], parent: List[float], offspring: List[float] | None = None
) -> tuple[float, float, float]:
    """Decompose evolutionary change using the Price equation.

    Args:
        fitness: Fitness values for each individual
        parent: Trait values in parent generation
        offspring: Trait values in offspring generation

    Returns:
        Tuple of (covariance_term, transmission_term, total_change)
    """
    if len(fitness) != len(parent) or (offspring is not None and len(parent) != len(offspring)):
        raise ValueError("All input lists must have the same length")

    # Handle empty inputs gracefully
    if not fitness or not parent:
        return (0.0, 0.0, 0.0)

    fitness_mean = statistics.mean(fitness)
    parent_mean = statistics.mean(parent)
    covariance_term = sum((f - fitness_mean) * (p - parent_mean) for f, p in zip(fitness, parent)) / len(fitness)
    covariance_term = covariance_term / fitness_mean if fitness_mean else 0.0

    if offspring is None:
        return covariance_term, 0.0, covariance_term

    offspring_mean = statistics.mean(offspring)
    total_change = offspring_mean - parent_mean
    transmission_term = total_change - covariance_term
    return covariance_term, transmission_term, total_change


def delta_mean_trait(
    fitness_or_trait: List[float], parent_or_fitness: List[float], offspring: List[float] | None = None
) -> float:
    """Calculate change in mean trait value due to selection.

    Args:
        trait_values: List of trait values
        fitness_values: Corresponding fitness values

    Returns:
        Expected change in mean trait value
    """
    if offspring is not None:
        return price_equation(fitness_or_trait, parent_or_fitness, offspring)[2]

    trait_values = fitness_or_trait
    fitness_values = parent_or_fitness
    if len(trait_values) != len(fitness_values):
        raise ValueError("Trait and fitness lists must have same length")

    if not trait_values:
        raise ValueError("Lists cannot be empty")

    # Δz = Cov(trait, fitness) / mean(fitness)
    trait_mean = statistics.mean(trait_values)
    fitness_mean = statistics.mean(fitness_values)

    covariance = sum((t - trait_mean) * (f - fitness_mean) for t, f in zip(trait_values, fitness_values))
    covariance /= len(trait_values)

    if fitness_mean == 0:
        return 0.0

    return covariance / fitness_mean


def expectation(values: List[float], weights: List[float] | None = None) -> float:
    """Calculate expectation (weighted mean).

    Args:
        values: List of values
        weights: Optional weights (defaults to uniform)

    Returns:
        Expected value
    """
    if not values:
        return 0.0

    if weights is None:
        return statistics.mean(values)

    if len(values) != len(weights):
        raise ValueError("Values and weights must have same length")

    total_weight = sum(weights)
    if total_weight == 0:
        return 0.0

    return sum(v * w for v, w in zip(values, weights)) / total_weight


def variance(values: List[float]) -> float:
    """Calculate variance.

    Args:
        values: List of values

    Returns:
        Variance
    """
    if len(values) < 2:
        return 0.0

    return statistics.variance(values)


def covariance(x: List[float], y: List[float]) -> float:
    """Calculate covariance between two lists.

    Args:
        x: First list of values
        y: Second list of values

    Returns:
        Covariance
    """
    if len(x) != len(y):
        raise ValueError("Lists must have same length")

    if len(x) < 2:
        return 0.0

    x_mean = statistics.mean(x)
    y_mean = statistics.mean(y)

    return sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(x, y)) / len(x)


def relative_fitness(fitness_values: List[float]) -> List[float]:
    """Convert absolute fitness to relative fitness.

    Args:
        fitness_values: List of absolute fitness values

    Returns:
        List of relative fitness values
    """
    if not fitness_values:
        return []

    fitness_mean = statistics.mean(fitness_values)
    if fitness_mean == 0:
        return [0.0] * len(fitness_values)

    return [f / fitness_mean for f in fitness_values]


def selection_intensity(fitness_values: List[float]) -> float:
    """Calculate selection intensity (i) from fitness values.

    Selection intensity is the standardized selection differential.

    Args:
        fitness_values: List of relative fitness values

    Returns:
        Selection intensity (i)
    """
    if not fitness_values:
        return 0.0

    fitness_array = [f for f in fitness_values if f > 0]  # Remove zero fitness
    if not fitness_array:
        return 0.0

    # Selection intensity = SD / mean, where SD is selection differential
    # For relative fitness, this is approximately the coefficient of variation
    fitness_mean = statistics.mean(fitness_array)
    if fitness_mean == 0:
        return 0.0

    fitness_std = statistics.stdev(fitness_array) if len(fitness_array) > 1 else 0.0

    return fitness_std / fitness_mean


def standard_deviation(values: List[float]) -> float:
    """Calculate standard deviation.

    Args:
        values: List of values

    Returns:
        Standard deviation
    """
    if len(values) < 2:
        return 0.0

    return statistics.stdev(values)


def weighted_variance(values: List[float], weights: List[float]) -> float:
    """Calculate weighted variance.

    Args:
        values: List of values
        weights: Corresponding weights

    Returns:
        Weighted variance
    """
    if len(values) != len(weights):
        raise ValueError("Values and weights must have same length")

    if not values:
        return 0.0

    total_weight = sum(weights)
    if total_weight == 0:
        return 0.0

    # Weighted mean
    weighted_mean = sum(v * w for v, w in zip(values, weights)) / total_weight

    # Weighted variance
    weighted_var = sum(w * (v - weighted_mean) ** 2 for v, w in zip(values, weights)) / total_weight

    return weighted_var


def weighted_covariance(x: List[float], y: List[float], weights: List[float]) -> float:
    """Calculate weighted covariance.

    Args:
        x: First list of values
        y: Second list of values
        weights: Corresponding weights

    Returns:
        Weighted covariance
    """
    if len(x) != len(y) or len(x) != len(weights):
        raise ValueError("All input lists must have same length")

    if not x:
        return 0.0

    total_weight = sum(weights)
    if total_weight == 0:
        return 0.0

    # Weighted means
    x_mean = sum(xi * w for xi, w in zip(x, weights)) / total_weight
    y_mean = sum(yi * w for yi, w in zip(y, weights)) / total_weight

    # Weighted covariance
    weighted_cov = sum(w * (xi - x_mean) * (yi - y_mean) for xi, yi, w in zip(x, y, weights)) / total_weight

    return weighted_cov


def weighted_correlation(x: List[float], y: List[float], weights: List[float]) -> float:
    """Calculate weighted correlation coefficient.

    Args:
        x: First list of values
        y: Second list of values
        weights: Corresponding weights

    Returns:
        Weighted correlation coefficient
    """
    weighted_cov_xy = weighted_covariance(x, y, weights)
    weighted_var_x = weighted_variance(x, weights)
    weighted_var_y = weighted_variance(y, weights)

    if weighted_var_x <= 0 or weighted_var_y <= 0:
        return 0.0

    return weighted_cov_xy / (weighted_var_x**0.5 * weighted_var_y**0.5)
