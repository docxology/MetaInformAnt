"""Linkage disequilibrium functions.

This module provides mathematical functions for analyzing linkage disequilibrium (LD).
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def ld_coefficients(genotypes: List[List[int]]) -> Dict[str, float]:
    """Calculate linkage disequilibrium coefficients.

    Args:
        genotypes: 2D list of genotypes [sample][locus]

    Returns:
        Dictionary with LD coefficients (D, D', r²)
    """
    if len(genotypes) < 2 or len(genotypes[0]) != 2:
        raise ValueError("Need at least 2 samples with 2 loci each")

    # Extract alleles for each locus
    locus1 = [row[0] for row in genotypes]
    locus2 = [row[1] for row in genotypes]

    # Calculate allele frequencies
    p1 = sum(locus1) / len(locus1)  # Frequency of allele 1 at locus 1
    q1 = 1 - p1  # Frequency of allele 0 at locus 1
    p2 = sum(locus2) / len(locus2)  # Frequency of allele 1 at locus 2
    q2 = 1 - p2  # Frequency of allele 0 at locus 2

    # Calculate haplotype frequencies (simplified)
    # Assuming 0/1 coding
    haplotype_freqs = {}
    for g1, g2 in zip(locus1, locus2):
        key = f"{g1}{g2}"
        haplotype_freqs[key] = haplotype_freqs.get(key, 0) + 1

    for key in haplotype_freqs:
        haplotype_freqs[key] /= len(genotypes)

    # D = P(AB) - P(A)P(B)
    D = haplotype_freqs.get('11', 0) - p1 * p2

    # D' = D / D_max
    D_max = min(p1*q2, q1*p2) if D < 0 else min(p1*p2, q1*q2)
    D_prime = D / D_max if D_max > 0 else 0

    # r² = D² / (p1*q1*p2*q2)
    r_squared = (D ** 2) / (p1 * q1 * p2 * q2) if (p1 * q1 * p2 * q2) > 0 else 0

    return {
        'D': D,
        'D_prime': D_prime,
        'r_squared': r_squared
    }


def ld_decay_r2(distances: List[float], r_squared_values: List[float],
               max_distance: Optional[float] = None) -> Dict[str, float]:
    """Analyze LD decay with distance.

    Args:
        distances: Physical distances between SNP pairs
        r_squared_values: Corresponding r² values
        max_distance: Maximum distance to consider

    Returns:
        Dictionary with LD decay statistics
    """
    if len(distances) != len(r_squared_values):
        raise ValueError("Distances and r² values must have same length")

    # Filter by max distance if specified
    if max_distance is not None:
        filtered = [(d, r) for d, r in zip(distances, r_squared_values) if d <= max_distance]
        distances, r_squared_values = zip(*filtered) if filtered else ([], [])

    if not distances:
        return {'decay_rate': 0.0, 'half_decay_distance': float('inf')}

    # Fit exponential decay: r² = r²₀ * exp(-d/d₀)
    # Use simple binning approach
    bins = np.logspace(0, np.log10(max(distances)), 20)
    bin_means = []

    for i in range(len(bins) - 1):
        mask = (np.array(distances) >= bins[i]) & (np.array(distances) < bins[i+1])
        if np.any(mask):
            mean_r2 = np.mean(np.array(r_squared_values)[mask])
            bin_means.append((bins[i], mean_r2))

    if len(bin_means) < 3:
        return {'decay_rate': 0.0, 'half_decay_distance': float('inf')}

    # Estimate half-decay distance (where r² drops to 0.5 of initial)
    initial_r2 = bin_means[0][1]
    half_decay_value = initial_r2 * 0.5

    half_decay_distance = float('inf')
    for dist, r2 in bin_means:
        if r2 <= half_decay_value:
            half_decay_distance = dist
            break

    # Estimate decay rate (rough approximation)
    if len(bin_means) >= 2:
        decay_rate = (bin_means[0][1] - bin_means[-1][1]) / (bin_means[-1][0] - bin_means[0][0])
    else:
        decay_rate = 0.0

    return {
        'decay_rate': decay_rate,
        'half_decay_distance': half_decay_distance,
        'initial_r2': initial_r2
    }


def haldane_c_to_d(recombination_fraction: float) -> float:
    """Convert recombination fraction to genetic distance using Haldane's mapping function.

    Args:
        recombination_fraction: Recombination fraction (c) between 0 and 0.5

    Returns:
        Genetic distance in Morgans (d)
    """
    if not (0 <= recombination_fraction <= 0.5):
        raise ValueError("Recombination fraction must be between 0 and 0.5")

    if recombination_fraction == 0:
        return 0.0
    elif recombination_fraction == 0.5:
        return float('inf')

    # Haldane's mapping function: d = -0.5 * ln(1 - 2c)
    return -0.5 * np.log(1 - 2 * recombination_fraction)


def haldane_d_to_c(genetic_distance: float) -> float:
    """Convert genetic distance (Morgans) to recombination fraction using Haldane's mapping function.

    Haldane's mapping function: c = 0.5(1 - exp(-2d))
    where d is genetic distance in Morgans and c is recombination fraction.

    Args:
        genetic_distance: Genetic distance in Morgans

    Returns:
        Recombination fraction (0 <= c <= 0.5)
    """
    import math
    return 0.5 * (1 - math.exp(-2 * genetic_distance))

def kosambi_c_to_d(recombination_fraction: float) -> float:
    """Convert recombination fraction to genetic distance using Kosambi mapping function.

    The Kosambi mapping function: c = 0.5 * tanh(2d)
    where c is recombination fraction and d is genetic distance in Morgans.

    Args:
        recombination_fraction: Recombination fraction (0 <= c <= 0.5)

    Returns:
        Genetic distance in Morgans
    """
    import math

    if not (0 <= recombination_fraction <= 0.5):
        raise ValueError("Recombination fraction must be between 0 and 0.5")

    if recombination_fraction == 0.5:
        return float('inf')  # Infinite distance
    elif recombination_fraction == 0:
        return 0.0

    return 0.25 * math.log((1 + 2 * recombination_fraction) / (1 - 2 * recombination_fraction))


def kosambi_d_to_c(genetic_distance: float) -> float:
    """Convert genetic distance to recombination fraction using Kosambi mapping function.

    The Kosambi mapping function: c = 0.5 * tanh(2d)
    where d is genetic distance in Morgans and c is recombination fraction.

    Args:
        genetic_distance: Genetic distance in Morgans

    Returns:
        Recombination fraction (0 <= c <= 0.5)
    """
    import math

    if genetic_distance < 0:
        raise ValueError("Genetic distance cannot be negative")

    if genetic_distance == 0:
        return 0.0

    return 0.5 * math.tanh(2 * genetic_distance)
