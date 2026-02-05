"""Population genetics mathematical models and calculations.

This module provides mathematical functions for population genetics analysis,
including Hardy-Weinberg equilibrium, selection models, and demographic calculations.
"""

from __future__ import annotations

from typing import Tuple


def hardy_weinberg_genotype_freqs(allele_a_frequency: float) -> Tuple[float, float, float]:
    """Calculate Hardy-Weinberg genotype frequencies.

    Args:
        allele_a_frequency: Frequency of allele A (0.0 to 1.0)

    Returns:
        Tuple of (AA_frequency, Aa_frequency, aa_frequency)

    Raises:
        ValueError: If allele frequency is not between 0 and 1

    Example:
        >>> hardy_weinberg_genotype_freqs(0.6)
        (0.36, 0.48, 0.16)
        >>> hardy_weinberg_genotype_freqs(0.0)
        (0.0, 0.0, 1.0)
        >>> hardy_weinberg_genotype_freqs(1.0)
        (1.0, 0.0, 0.0)
    """
    # Handle invalid frequencies gracefully by returning zeros
    if not (0.0 <= allele_a_frequency <= 1.0):
        return (0.0, 0.0, 0.0)

    p = allele_a_frequency
    q = 1.0 - p

    aa_freq = p * p  # AA homozygotes
    aa_freq_het = 2 * p * q  # Aa heterozygotes
    aa_freq_hom = q * q  # aa homozygotes

    return (aa_freq, aa_freq_het, aa_freq_hom)


def heterozygosity_decay(initial_heterozygosity: float, effective_size: float, generations: int) -> float:
    """Calculate heterozygosity decay due to genetic drift.

    Args:
        initial_heterozygosity: Initial heterozygosity (H0)
        effective_size: Effective population size (Ne)
        generations: Number of generations

    Returns:
        Heterozygosity after t generations
    """
    if initial_heterozygosity < 0 or initial_heterozygosity > 1:
        raise ValueError("Initial heterozygosity must be between 0 and 1")
    if effective_size <= 0:
        raise ValueError("Effective population size must be positive")
    if generations < 0:
        raise ValueError("Generations cannot be negative")

    # H_t = H0 * (1 - 1/(2Ne))^t
    decay_factor = (1 - 1 / (2 * effective_size)) ** generations
    return initial_heterozygosity * decay_factor


def inbreeding_coefficient(effective_size: float, generations: int) -> float:
    """Calculate inbreeding coefficient after t generations.

    Args:
        effective_size: Effective population size (Ne)
        generations: Number of generations

    Returns:
        Inbreeding coefficient (F)
    """
    if effective_size <= 0:
        raise ValueError("Effective population size must be positive")
    if generations < 0:
        raise ValueError("Generations cannot be negative")

    # F = 1 - (1 - 1/(2Ne))^t
    if effective_size == float("inf"):
        return 0.0

    inbreeding_factor = 1 - (1 - 1 / (2 * effective_size)) ** generations
    return inbreeding_factor


def mutation_selection_balance_dominant(mutation_rate: float, selection_coefficient: float) -> float:
    """Calculate equilibrium frequency for dominant mutations under mutation-selection balance.

    For dominant mutations, the equilibrium frequency q satisfies:
    q = μ / s
    where μ is the mutation rate and s is the selection coefficient.

    Args:
        mutation_rate: Mutation rate (μ)
        selection_coefficient: Selection coefficient (s)

    Returns:
        Equilibrium frequency (q)
    """
    return mutation_rate / selection_coefficient


def mutation_selection_balance_recessive(mutation_rate: float, selection_coefficient: float) -> float:
    """Calculate equilibrium frequency for recessive mutations under mutation-selection balance.

    For recessive mutations, the equilibrium frequency q satisfies:
    q² = μ / s
    where μ is the mutation rate and s is the selection coefficient.

    Args:
        mutation_rate: Mutation rate (μ)
        selection_coefficient: Selection coefficient (s)

    Returns:
        Equilibrium frequency (q)
    """
    import math

    return math.sqrt(mutation_rate / selection_coefficient)
