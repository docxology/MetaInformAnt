"""Coalescent theory mathematical models.

This module provides functions for coalescent theory calculations,
including expected times to most recent common ancestor and
related demographic parameters.
"""

from __future__ import annotations

import math
from typing import Any, Dict

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def expected_time_to_mrca(n_samples: int, effective_size: float) -> float:
    """Calculate expected time to most recent common ancestor.

    Args:
        n_samples: Number of sampled lineages
        effective_size: Effective population size (diploid)

    Returns:
        Expected time in generations (4Ne for diploid populations)

    Raises:
        ValueError: If n_samples < 2
    """
    if n_samples < 2:
        raise ValueError("Need at least 2 samples")
    # Handle zero/invalid effective size gracefully
    if effective_size <= 0:
        return 0.0

    # For diploid populations, time is scaled by 4Ne
    # Expected time for n lineages to coalesce to 1 is sum_{k=2 to n} 4Ne / (k(k-1))
    expected_time = 0.0
    for k in range(2, n_samples + 1):
        expected_time += 4 * effective_size / (k * (k - 1))

    return expected_time


def watterson_theta(
    S: int | None = None,
    n: int | None = None,
    sequence_length: int = 1,
    *,
    n_sites: int | None = None,
    n_segregating_sites: int | None = None,
    n_sequences: int | None = None,
) -> float:
    """Calculate Watterson's θ estimator.

    Supports multiple calling conventions:
    - watterson_theta(S, n) - S = segregating sites, n = number of sequences
    - watterson_theta(n_sites=N, n_segregating_sites=S) - explicit names

    Args:
        S: Number of segregating sites (positional)
        n: Number of sequences sampled (positional)
        sequence_length: Length of sequence (default 1 for per-site)
        n_sites: Alias for n (number of sequences)
        n_segregating_sites: Alias for S
        n_sequences: Alias for n

    Returns:
        Watterson's θ estimate

    Raises:
        ValueError: If parameters are invalid
    """
    # Handle parameter aliases
    # When called as watterson_theta(S, n), S goes to first positional, n to second
    segregating = S
    num_seqs = n

    # Check for keyword aliases
    if n_segregating_sites is not None:
        segregating = n_segregating_sites
    if n_sites is not None:
        num_seqs = n_sites
    if n_sequences is not None:
        num_seqs = n_sequences

    if segregating is None:
        raise ValueError("Number of segregating sites must be provided")
    if num_seqs is None:
        raise ValueError("Number of sequences must be provided")

    if num_seqs <= 1:
        raise ValueError("Number of sequences must be > 1")
    if segregating < 0:
        raise ValueError("Number of segregating sites cannot be negative")

    # Harmonic number a_1 = sum(1/i for i in 1..n-1)
    harmonic_sum = 0.0
    try:
        from scipy import special

        # psi is digamma function. H_n = psi(n+1) + gamma
        harmonic_sum = special.digamma(num_seqs) + np.euler_gamma
    except ImportError:
        # Fallback: exact sum or approximation
        if num_seqs > 100:
            harmonic_sum = math.log(num_seqs - 1) + 0.5772156649015329  # γ
        else:
            harmonic_sum = sum(1.0 / i for i in range(1, num_seqs))

    if harmonic_sum == 0:
        return 0.0

    return segregating / (harmonic_sum * sequence_length)


def expected_coalescent_waiting_times(n_samples: int, effective_size: float) -> List[float]:
    """Calculate expected waiting times between coalescent events.

    Args:
        n_samples: Number of sampled lineages
        effective_size: Effective population size

    Returns:
        List of expected waiting times for each coalescent interval

    Raises:
        ValueError: If parameters are invalid
    """
    if n_samples < 2:
        raise ValueError("Need at least 2 samples")
    if effective_size <= 0:
        raise ValueError("Effective size must be positive")

    waiting_times = []
    for k in range(n_samples, 1, -1):
        # Expected time for k lineages to coalesce to k-1
        expected_time = 4 * effective_size / (k * (k - 1))
        waiting_times.append(expected_time)

    return waiting_times


def expected_pairwise_diversity(effective_size: float, mutation_rate: float) -> float:
    """Calculate expected pairwise nucleotide diversity.

    Args:
        effective_size: Effective population size (Ne)
        mutation_rate: Per-site mutation rate (μ)

    Returns:
        Expected pairwise diversity π = 4Neμ
    """
    if effective_size <= 0:
        raise ValueError("Effective size must be positive")
    if mutation_rate < 0:
        raise ValueError("Mutation rate cannot be negative")

    return 4 * effective_size * mutation_rate


def expected_pairwise_diversity_from_theta(theta: float) -> float:
    """Alias for expected_pairwise_diversity with fixed n_samples."""
    return theta


def expected_segregating_sites(n_samples: int, theta: float, sequence_length: int = 1) -> float:
    """Calculate expected number of segregating sites.

    Args:
        n_samples: Number of sequences
        theta: Population mutation parameter
        sequence_length: Length of sequence (default 1 for per-site)

    Returns:
        Expected number of segregating sites
    """
    if n_samples < 2:
        return 0.0

    # For neutral infinite sites model, E[S] ≈ θ * H_{n-1} * L
    # where H is the harmonic number and L is sequence length
    harmonic_sum = sum(1.0 / i for i in range(1, n_samples))
    return theta * harmonic_sum * sequence_length


def expected_sfs_counts(n_samples: int, theta: float) -> List[float]:
    """Calculate expected site frequency spectrum counts.

    Args:
        n_samples: Number of sequences (haploid)
        theta: Population mutation parameter

    Returns:
        List of expected counts for each frequency class
    """
    if n_samples < 2:
        return []

    # Neutral SFS: E[ξ_i] = θ/i for i=1 to n-1
    sfs = []
    for i in range(1, n_samples):
        expected_count = theta / i
        sfs.append(expected_count)

    return sfs


def expected_total_branch_length(n_samples: int, effective_size: float) -> float:
    """Calculate expected total branch length in coalescent tree.

    Args:
        n_samples: Number of sampled lineages
        effective_size: Effective population size

    Returns:
        Expected total branch length
    """
    if n_samples < 2:
        return 0.0

    # Calculate harmonic sum H_{n-1}
    if n_samples - 1 > 100:
        harmonic_sum = math.log(n_samples - 1) + 0.5772156649015329  # γ
    else:
        harmonic_sum = sum(1.0 / i for i in range(1, n_samples))

    return 4 * effective_size * harmonic_sum


def tajima_constants(n: int) -> Dict[str, float]:
    """Calculate Tajima's D constants.

    Args:
        n: Sample size

    Returns:
        Dictionary with Tajima's constants a1, a2, b1, b2, c1, c2, e1, e2
    """
    if n < 2:
        return {"a1": 0.0, "a2": 0.0, "b1": 0.0, "b2": 0.0, "c1": 0.0, "c2": 0.0, "e1": 0.0, "e2": 0.0}

    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i * i) for i in range(1, n))

    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))

    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)

    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)

    return {"a1": a1, "a2": a2, "b1": b1, "b2": b2, "c1": c1, "c2": c2, "e1": e1, "e2": e2}


def tajimas_D(pi: float, s: int, n: int, tajima_constants: Dict[str, float] | None = None) -> float:
    """Calculate Tajima's D statistic.

    Args:
        pi: Nucleotide diversity
        s: Number of segregating sites
        n: Sample size
        tajima_constants: Pre-computed Tajima constants (optional)

    Returns:
        Tajima's D value
    """
    if n < 2 or s == 0:
        return 0.0

    if tajima_constants is None:
        tajima_constants = tajima_constants(n)

    a1 = tajima_constants["a1"]

    # Difference between diversity and segregating sites
    diff = pi - s / a1

    # Variance of the difference
    b1 = tajima_constants["b1"]
    b2 = tajima_constants["b2"]
    c1 = tajima_constants["c1"]
    c2 = tajima_constants["c2"]
    e1 = tajima_constants["e1"]
    e2 = tajima_constants["e2"]

    var = e1 * s + e2 * s * (s - 1)

    if var <= 0:
        return 0.0

    return diff / math.sqrt(var)


def simulate_coalescent(n_samples: int, effective_size: float = 10000, mutation_rate: float = 1e-8) -> Dict[str, Any]:
    """Simulate coalescent process.

    Args:
        n_samples: Number of sampled lineages
        effective_size: Effective population size
        mutation_rate: Per-site mutation rate

    Returns:
        Dictionary with coalescent simulation results

    Raises:
        ValueError: If parameters are invalid
    """
    if n_samples < 2:
        raise ValueError("Need at least 2 samples")
    if effective_size <= 0:
        raise ValueError("Effective size must be positive")
    if mutation_rate < 0:
        raise ValueError("Mutation rate cannot be negative")

    # Simple coalescent simulation
    lineages = list(range(n_samples))
    coalescence_times = []
    current_time = 0.0

    while len(lineages) > 1:
        # Coalescence rate for k lineages: k(k-1)/(4Ne)
        k = len(lineages)
        rate = k * (k - 1) / (4 * effective_size)

        if rate <= 0:
            break

        # Time to next coalescence (exponential distribution)
        time_increment = -math.log(np.random.random()) / rate
        current_time += time_increment

        coalescence_times.append(current_time)

        # Randomly merge two lineages
        i, j = np.random.choice(len(lineages), 2, replace=False)
        # Remove the second one (keep the first)
        lineages.pop(j)

    # Calculate summary statistics
    results = {
        "n_samples": n_samples,
        "effective_size": effective_size,
        "mutation_rate": mutation_rate,
        "coalescence_times": coalescence_times,
        "total_time": current_time if coalescence_times else 0.0,
        "theta_watterson": watterson_theta(1000, len(coalescence_times)) if coalescence_times else 0.0,
    }

    return results


def coalescent_time_to_mrca(n_samples: int, effective_size: float) -> float:
    """Alias for expected_time_to_mrca for backward compatibility."""
    return expected_time_to_mrca(n_samples, effective_size)
