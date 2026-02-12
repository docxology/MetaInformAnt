"""Selection theory functions.

This module provides mathematical functions for evolutionary selection analysis,
including quantitative genetics and kin selection theory.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core.utils import logging

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


def kin_selection_response(
    relatedness: float | None = None,
    benefit: float | None = None,
    cost: float | None = None,
    *,
    r: float | None = None,
    b: float | None = None,
    c: float | None = None,
    individual_fitness: float | None = None,
    kin_fitness: float | None = None,
) -> float:
    """Calculate kin selection response using Hamilton's rule.

    Supports multiple calling conventions:
    - kin_selection_response(r, b, c) - relatedness, benefit, cost
    - kin_selection_response(relatedness=r, benefit=b, cost=c)
    - kin_selection_response(individual_fitness, kin_fitness, relatedness) - legacy

    Hamilton's rule: r*b - c > 0 for altruism to evolve

    Args:
        relatedness: Coefficient of relatedness (r)
        benefit: Fitness benefit to kin (b)
        cost: Fitness cost to individual (c)
        r, b, c: Aliases for relatedness, benefit, cost
        individual_fitness: Legacy alias for cost
        kin_fitness: Legacy alias for benefit

    Returns:
        Net selection response (r*b - c)
    """
    # Handle parameter aliases
    if r is not None:
        relatedness = r
    if b is not None:
        benefit = b
    if c is not None:
        cost = c
    if individual_fitness is not None and cost is None:
        cost = individual_fitness
    if kin_fitness is not None and benefit is None:
        benefit = kin_fitness

    if relatedness is None or benefit is None or cost is None:
        raise ValueError("Must provide relatedness, benefit, and cost")

    # Hamilton's rule: rb - c > 0
    return relatedness * benefit - cost


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
    allele_frequency: float,
    selection_coefficient: float | None = None,
    dominance_coefficient: float = 0.5,
    *,
    fitness_AA: float | None = None,
    fitness_Aa: float | None = None,
    fitness_aa: float | None = None,
) -> float:
    """Update allele frequency due to selection.

    Supports two parameter styles:
    1. selection_coefficient + dominance_coefficient (traditional)
    2. fitness_AA, fitness_Aa, fitness_aa (explicit fitness values)

    Args:
        allele_frequency: Current allele frequency (p) for allele A
        selection_coefficient: Selection coefficient (s) - used if fitness values not provided
        dominance_coefficient: Dominance coefficient (h) - default 0.5 for additive
        fitness_AA: Fitness of AA homozygote
        fitness_Aa: Fitness of Aa heterozygote
        fitness_aa: Fitness of aa homozygote

    Returns:
        New allele frequency after selection
    """
    if not (0 <= allele_frequency <= 1):
        raise ValueError("Allele frequency must be between 0 and 1")

    p = allele_frequency
    q = 1 - p

    # If fitness values are provided, use them directly
    if fitness_AA is not None and fitness_Aa is not None and fitness_aa is not None:
        # Calculate genotype frequencies under HWE
        freq_AA = p * p
        freq_Aa = 2 * p * q
        freq_aa = q * q

        # Mean fitness
        w_bar = freq_AA * fitness_AA + freq_Aa * fitness_Aa + freq_aa * fitness_aa

        if w_bar == 0:
            return allele_frequency

        # New allele frequency: p' = (p² * w_AA + p*q * w_Aa) / w_bar
        p_new = (freq_AA * fitness_AA + 0.5 * freq_Aa * fitness_Aa) / w_bar

        return max(0, min(1, p_new))

    # Otherwise use selection coefficient approach
    if selection_coefficient is None:
        selection_coefficient = 0.0

    # Selection: Δp = s*p*(1-p)*(p*h + (1-p))/mean_fitness
    # Simplified for haploid or additive case
    if dominance_coefficient == 0.5:  # Additive
        # Δp ≈ s*p*(1-p)/2 for weak selection
        delta_p = selection_coefficient * p * q / 2
    else:
        # General case: Δp = s*p*(1-p)*[p*h + (1-p)] / w_bar
        # Approximate for weak selection
        fitness_effect = p * dominance_coefficient + q
        delta_p = selection_coefficient * p * q * fitness_effect

    new_frequency = p + delta_p

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
