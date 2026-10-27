"""Population demography and growth models.

This module provides mathematical models for population growth,
age structure dynamics, and demographic processes.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Sequence

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def exponential_growth_model(initial_size: float, growth_rate: float, generations: int) -> List[float]:
    """Simulate exponential population growth.

    Args:
        initial_size: Initial population size
        growth_rate: Growth rate per generation
        generations: Number of generations to simulate

    Returns:
        List of population sizes for each generation

    Raises:
        ValueError: If parameters are invalid
    """
    if initial_size <= 0:
        raise ValueError("Initial population size must be positive")
    if generations < 0:
        raise ValueError("Number of generations cannot be negative")

    population_sizes = [initial_size]

    for _ in range(generations):
        current_size = population_sizes[-1]
        new_size = current_size * (1 + growth_rate)
        population_sizes.append(new_size)

    return population_sizes


def logistic_growth_model(
    initial_size: float, carrying_capacity: float, growth_rate: float, generations: int
) -> List[float]:
    """Simulate logistic population growth.

    Args:
        initial_size: Initial population size
        carrying_capacity: Maximum population size (K)
        growth_rate: Intrinsic growth rate (r)
        generations: Number of generations to simulate

    Returns:
        List of population sizes for each generation

    Raises:
        ValueError: If parameters are invalid
    """
    if initial_size <= 0:
        raise ValueError("Initial population size must be positive")
    if carrying_capacity <= 0:
        raise ValueError("Carrying capacity must be positive")
    if generations < 0:
        raise ValueError("Number of generations cannot be negative")

    population_sizes = [initial_size]

    for _ in range(generations):
        current_size = population_sizes[-1]
        # Logistic growth equation: dN/dt = rN(1 - N/K)
        # Discrete form: N(t+1) = N(t) + r*N(t)*(1 - N(t)/K)
        growth = growth_rate * current_size * (1 - current_size / carrying_capacity)
        new_size = current_size + growth
        new_size = max(0, new_size)  # Prevent negative populations
        population_sizes.append(new_size)

    return population_sizes


def age_structure_model(
    fertility_rates: Sequence[float],
    survival_rates: Sequence[float],
    initial_age_structure: Sequence[float],
    generations: int,
) -> Dict[str, Any]:
    """Simulate age-structured population dynamics.

    Args:
        fertility_rates: Fertility rates for each age class
        survival_rates: Survival rates between age classes
        initial_age_structure: Initial population distribution across age classes
        generations: Number of generations to simulate

    Returns:
        Dictionary containing population dynamics results

    Raises:
        ValueError: If parameters are invalid
    """
    if len(fertility_rates) != len(survival_rates) or len(survival_rates) != len(initial_age_structure):
        raise ValueError("All parameter arrays must have the same length")
    if generations < 0:
        raise ValueError("Number of generations cannot be negative")

    n_ages = len(fertility_rates)
    population_history = [list(initial_age_structure)]

    for gen in range(generations):
        current_pop = population_history[-1]

        # Calculate new population for each age class
        new_pop = [0.0] * n_ages

        # Age 0: newborns from all reproductive age classes
        newborns = 0.0
        for age in range(n_ages):
            if age < len(fertility_rates):  # Reproductive ages
                newborns += current_pop[age] * fertility_rates[age]

        new_pop[0] = newborns

        # Ages 1+: survivors from previous age classes
        for age in range(1, n_ages):
            new_pop[age] = current_pop[age - 1] * survival_rates[age - 1]

        population_history.append(new_pop)

    # Calculate summary statistics
    total_populations = [sum(gen_pop) for gen_pop in population_history]
    age_distributions = population_history

    results = {
        "population_history": population_history,
        "total_populations": total_populations,
        "age_distributions": age_distributions,
        "final_population": total_populations[-1] if total_populations else 0,
        "growth_rate": (
            (total_populations[-1] / total_populations[0]) ** (1 / generations) - 1
            if generations > 0 and total_populations[0] > 0
            else 0
        ),
        "stable_age_distribution": age_distributions[-1] if age_distributions else [],
    }

    return results


def _harmonic_mean_size(sizes: Sequence[float]) -> float:
    """Wright's (1938) effective size of a population whose census size varies by generation."""
    return len(sizes) / math.fsum(1.0 / size for size in sizes)


def bottleneck_effective_size(
    pre_bottleneck_size: float, bottleneck_size: float, post_bottleneck_size: float, recovery_generations: int = 0
) -> float:
    """Calculate the effective population size across a one-generation bottleneck.

    Ne is the harmonic mean of the census sizes (Wright 1938) over one generation
    at each of the three sizes plus ``recovery_generations`` further generations
    at the post-bottleneck size:
    Ne = (3 + k) / (1/N_pre + 1/N_bottleneck + (1 + k)/N_post).

    Args:
        pre_bottleneck_size: Population size before the bottleneck
        bottleneck_size: Population size during the bottleneck generation
        post_bottleneck_size: Population size after the bottleneck
        recovery_generations: Additional generations spent at the post-bottleneck size

    Returns:
        Effective population size

    Raises:
        ValueError: If a size is not positive or ``recovery_generations`` is negative
    """
    if min(pre_bottleneck_size, bottleneck_size, post_bottleneck_size) <= 0:
        raise ValueError("Population sizes must be positive")
    if recovery_generations < 0:
        raise ValueError("Recovery generations cannot be negative")

    return _harmonic_mean_size(
        [pre_bottleneck_size, bottleneck_size] + [post_bottleneck_size] * (1 + recovery_generations)
    )


def exponential_growth_effective_size(initial_size: float, growth_rate: float, generations: int) -> float:
    """Calculate the effective size of an exponentially growing or declining population.

    The census size follows :func:`exponential_growth_model`, N_i = N_0 (1 + r)^i,
    and Ne over generations 0..t-1 is the harmonic mean of those sizes
    (Wright 1938). Growth therefore gives an Ne above N_0 but weighted toward the
    early, small generations; decline gives an Ne below N_0.

    Args:
        initial_size: Population size at generation 0 (N_0)
        growth_rate: Growth rate per generation (r), greater than -1
        generations: Number of generations (t)

    Returns:
        Effective population size; ``initial_size`` when ``generations`` is 0

    Raises:
        ValueError: If ``initial_size`` is not positive, ``generations`` is
            negative or ``growth_rate`` is -1 or less
    """
    if initial_size <= 0:
        raise ValueError("Initial population size must be positive")
    if generations < 0:
        raise ValueError("Generations cannot be negative")
    if growth_rate <= -1:
        raise ValueError("Growth rate must exceed -1 so population sizes stay positive")

    if generations == 0 or growth_rate == 0:
        return float(initial_size)
    return _harmonic_mean_size([initial_size * (1 + growth_rate) ** i for i in range(generations)])


def two_epoch_effective_size(ancient_size: float, recent_size: float, split_time: int) -> float:
    """Calculate the coalescent effective size of a two-epoch population.

    The population had ``ancient_size`` until ``split_time`` generations ago and
    ``recent_size`` since. Ne is the constant size with the same expected pairwise
    coalescence time, and hence the same expected nucleotide diversity. For a
    diploid population,
    Ne = N_recent + (N_ancient - N_recent) exp(-T / (2 N_recent)),
    so diversity still reflects the ancient size shortly after a change and the
    recent size long after it.

    Args:
        ancient_size: Population size before the change
        recent_size: Population size since the change
        split_time: Generations since the size changed (T)

    Returns:
        Effective population size; ``ancient_size`` when ``split_time`` is 0

    Raises:
        ValueError: If a size is not positive or ``split_time`` is negative
    """
    if ancient_size <= 0 or recent_size <= 0:
        raise ValueError("Population sizes must be positive")
    if split_time < 0:
        raise ValueError("Split time cannot be negative")

    return recent_size + (ancient_size - recent_size) * math.exp(-split_time / (2.0 * recent_size))


def island_model_update(allele_frequency: float, migration_rate: float, migrant_allele_frequency: float) -> float:
    """Update allele frequency under the island model of migration.

    Args:
        allele_frequency: Current allele frequency (p)
        migration_rate: Proportion of migrants (m)
        migrant_allele_frequency: Allele frequency in migrant pool (pm)

    Returns:
        New allele frequency
    """
    if not (0 <= allele_frequency <= 1):
        raise ValueError("Allele frequency must be between 0 and 1")
    if not (0 <= migrant_allele_frequency <= 1):
        raise ValueError("Migrant allele frequency must be between 0 and 1")
    if not (0 <= migration_rate <= 1):
        raise ValueError("Migration rate must be between 0 and 1")

    # p' = (1-m)p + m*pm = p - mp + m*pm = p + m(pm - p)
    delta_p = migration_rate * (migrant_allele_frequency - allele_frequency)
    return allele_frequency + delta_p
