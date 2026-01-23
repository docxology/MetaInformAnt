"""Selection theory functions.

This module provides mathematical functions for evolutionary selection analysis,
including quantitative genetics and kin selection theory.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def breeders_equation_response(heritability: float, selection_differential: float) -> float:
    """Calculate response to selection using breeder's equation.

    Args:
        heritability: Narrow-sense heritability (h²)
        selection_differential: Selection differential (S)

    Returns:
        Response to selection (R)
    """
    # R = h² * S
    return heritability * selection_differential


def kin_selection_response(individual_fitness: float, kin_fitness: float, relatedness: float) -> float:
    """Calculate kin selection response.

    Args:
        individual_fitness: Fitness cost to individual
        kin_fitness: Fitness benefit to kin
        relatedness: Coefficient of relatedness (r)

    Returns:
        Net selection response
    """
    # Hamilton's rule: rb - c > 0
    # Net response = benefit to kin - cost to individual
    return relatedness * kin_fitness - individual_fitness


def mutation_update(allele_frequency: float, mutation_rate_forward: float, mutation_rate_backward: float) -> float:
    """Update allele frequency due to mutation.

    Args:
        allele_frequency: Current allele frequency (p)
        mutation_rate_forward: Forward mutation rate (u: A→a)
        mutation_rate_backward: Backward mutation rate (v: a→A)

    Returns:
        New allele frequency after mutation
    """
    if not (0 <= allele_frequency <= 1):
        raise ValueError("Allele frequency must be between 0 and 1")

    # Mutation balance: Δp = -u*p + v*(1-p)
    delta_p = -mutation_rate_forward * allele_frequency + mutation_rate_backward * (1 - allele_frequency)

    new_frequency = allele_frequency + delta_p

    # Ensure bounds
    return max(0, min(1, new_frequency))


def selection_update(
    allele_frequency: float, selection_coefficient: float, dominance_coefficient: float = 0.5
) -> float:
    """Update allele frequency due to selection.

    Args:
        allele_frequency: Current allele frequency (p)
        selection_coefficient: Selection coefficient (s)
        dominance_coefficient: Dominance coefficient (h)

    Returns:
        New allele frequency after selection
    """
    if not (0 <= allele_frequency <= 1):
        raise ValueError("Allele frequency must be between 0 and 1")

    # Selection: Δp = s*p*(1-p)*(p*h + (1-p))/mean_fitness
    # Simplified for haploid or additive case
    if dominance_coefficient == 0.5:  # Additive
        # Δp ≈ s*p*(1-p)/2 for weak selection
        delta_p = selection_coefficient * allele_frequency * (1 - allele_frequency) / 2
    else:
        # General case: Δp = s*p*(1-p)*[p*h + (1-p)] / w_bar
        # Approximate for weak selection
        fitness_effect = allele_frequency * dominance_coefficient + (1 - allele_frequency)
        delta_p = selection_coefficient * allele_frequency * (1 - allele_frequency) * fitness_effect

    new_frequency = allele_frequency + delta_p

    # Ensure bounds
    return max(0, min(1, new_frequency))


def selection_differential(trait_values: List[float], fitness_values: List[float]) -> float:
    """Calculate selection differential.

    Args:
        trait_values: List of trait values
        fitness_values: Corresponding fitness values

    Returns:
        Selection differential (S)
    """
    if len(trait_values) != len(fitness_values):
        raise ValueError("Trait and fitness lists must have same length")

    if not trait_values:
        raise ValueError("Lists cannot be empty")

    # S = Cov(trait, fitness) / mean(fitness)
    trait_mean = sum(trait_values) / len(trait_values)
    fitness_mean = sum(fitness_values) / len(fitness_values)

    covariance = sum((t - trait_mean) * (f - fitness_mean) for t, f in zip(trait_values, fitness_values))
    covariance /= len(trait_values)

    if fitness_mean == 0:
        return 0.0

    return covariance / fitness_mean


def selection_gradient(trait_values: List[float], fitness_values: List[float]) -> float:
    """Calculate selection gradient.

    Args:
        trait_values: List of trait values
        fitness_values: Corresponding fitness values

    Returns:
        Selection gradient (β)
    """
    if len(trait_values) != len(fitness_values):
        raise ValueError("Trait and fitness lists must have same length")

    if not trait_values:
        raise ValueError("Lists cannot be empty")

    # β = Cov(trait, relative_fitness) / Var(trait)
    relative_fitness = [f / sum(fitness_values) * len(fitness_values) for f in fitness_values]

    trait_mean = sum(trait_values) / len(trait_values)
    fitness_rel_mean = sum(relative_fitness) / len(relative_fitness)

    covariance = sum((t - trait_mean) * (f - fitness_rel_mean) for t, f in zip(trait_values, relative_fitness))
    covariance /= len(trait_values)

    trait_variance = sum((t - trait_mean) ** 2 for t in trait_values) / len(trait_values)

    if trait_variance == 0:
        return 0.0

    return covariance / trait_variance


def relative_fitness(fitness_values: List[float]) -> List[float]:
    """Calculate relative fitness values.

    Args:
        fitness_values: List of absolute fitness values

    Returns:
        List of relative fitness values (w / mean(w))
    """
    if not fitness_values:
        return []

    mean_fitness = sum(fitness_values) / len(fitness_values)
    if mean_fitness == 0:
        return [0.0] * len(fitness_values)

    return [w / mean_fitness for w in fitness_values]


def selection_intensity(trait_values: List[float], fitness_values: List[float]) -> float:
    """Calculate selection intensity.

    Args:
        trait_values: List of trait values
        fitness_values: Corresponding fitness values

    Returns:
        Selection intensity (i)
    """
    # i = S / sigma_trait
    from .statistics import standard_deviation

    S = selection_differential(trait_values, fitness_values)
    sigma = standard_deviation(trait_values)

    if sigma == 0:
        return 0.0

    if sigma == 0:
        return 0.0

    return S / sigma


def mutation_selection_balance_recessive(mutation_rate: float, selection_coefficient: float) -> float:
    """Calculate equilibrium frequency for recessive deleterious allele.

    Args:
        mutation_rate: Mutation rate to deleterious allele (mu)
        selection_coefficient: Selection coefficient against homozygotes (s)

    Returns:
        Equilibrium frequency (q_eq)
    """
    if selection_coefficient <= 0:
        raise ValueError("Selection coefficient must be positive for deleterious allele")
    if mutation_rate < 0:
        raise ValueError("Mutation rate cannot be negative")

    # q ~ sqrt(mu/s)
    return (mutation_rate / selection_coefficient) ** 0.5


def mutation_selection_balance_dominant(mutation_rate: float, selection_coefficient: float) -> float:
    """Calculate equilibrium frequency for dominant deleterious allele.

    Args:
        mutation_rate: Mutation rate to deleterious allele (mu)
        selection_coefficient: Selection coefficient against heterozygotes (s)

    Returns:
        Equilibrium frequency (q_eq)
    """
    if selection_coefficient <= 0:
        raise ValueError("Selection coefficient must be positive for deleterious allele")
    if mutation_rate < 0:
        raise ValueError("Mutation rate cannot be negative")

    # q ~ mu/s (assuming h=1 for simplicity in this context, or generalized s)
    # Ideally q = mu / (h*s), if fully dominant h=1 -> mu/s
    return mutation_rate / selection_coefficient


def multilevel_selection_decomposition(
    group_trait_mean: float, individual_trait_variance: float, group_trait_variance: float, group_size: int
) -> Dict[str, float]:
    """Decompose selection into individual and group-level components.

    Args:
        group_trait_mean: Mean trait value across groups
        individual_trait_variance: Variance within groups
        group_trait_variance: Variance between groups
        group_size: Average group size

    Returns:
        Dictionary with multilevel selection components
    """
    if group_size <= 0:
        raise ValueError("Group size must be positive")

    # Total phenotypic variance
    total_variance = individual_trait_variance + group_trait_variance * group_size

    # Individual-level selection
    individual_selection = individual_trait_variance / total_variance if total_variance > 0 else 0

    # Group-level selection
    group_selection = (group_trait_variance * group_size) / total_variance if total_variance > 0 else 0

    return {
        "total_variance": total_variance,
        "individual_level_selection": individual_selection,
        "group_level_selection": group_selection,
        "multilevel_selection_ratio": group_selection / individual_selection if individual_selection > 0 else 0,
    }
