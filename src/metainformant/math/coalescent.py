"""Coalescent theory mathematical models.

This module provides functions for coalescent theory calculations,
including expected times to most recent common ancestor and
related demographic parameters.
"""

from __future__ import annotations

from typing import Dict, Any
import math

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
        ValueError: If parameters are invalid
    """
    if n_samples < 2:
        raise ValueError("Need at least 2 samples")
    if effective_size <= 0:
        raise ValueError("Effective size must be positive")

    # For diploid populations, time is scaled by 4Ne
    # Expected time for n lineages to coalesce to 1 is sum_{k=2 to n} 4Ne / (k(k-1))
    expected_time = 0.0
    for k in range(2, n_samples + 1):
        expected_time += 4 * effective_size / (k * (k - 1))

    return expected_time


def watterson_theta(n_sites: int, n_segregating_sites: int, sequence_length: int = 1) -> float:
    """Calculate Watterson's θ estimator.

    Args:
        n_sites: Number of sequences sampled
        n_segregating_sites: Number of segregating sites
        sequence_length: Length of sequence (default 1 for per-site)

    Returns:
        Watterson's θ estimate

    Raises:
        ValueError: If parameters are invalid
    """
    if n_sites <= 1:
        raise ValueError("Number of sites must be > 1")
    if n_segregating_sites < 0:
        raise ValueError("Number of segregating sites cannot be negative")
    if n_segregating_sites > sequence_length:
        raise ValueError("Segregating sites cannot exceed sequence length")

    # Harmonic number H_{n-1} for n sequences
    # H_{n-1} ≈ ln(n-1) + γ where γ is Euler-Mascheroni constant
    n_sequences = n_sites
    if n_sequences > 100:
        harmonic_sum = math.log(n_sequences - 1) + 0.5772156649015329  # γ
    else:
        harmonic_sum = sum(1.0 / i for i in range(1, n_sequences))

    if harmonic_sum == 0:
        return 0.0

    return n_segregating_sites / (harmonic_sum * sequence_length)


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
        return {'a1': 0.0, 'a2': 0.0, 'b1': 0.0, 'b2': 0.0,
                'c1': 0.0, 'c2': 0.0, 'e1': 0.0, 'e2': 0.0}

    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i * i) for i in range(1, n))

    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))

    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)

    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)

    return {
        'a1': a1, 'a2': a2, 'b1': b1, 'b2': b2,
        'c1': c1, 'c2': c2, 'e1': e1, 'e2': e2
    }


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

    a1 = tajima_constants['a1']

    # Difference between diversity and segregating sites
    diff = pi - s / a1

    # Variance of the difference
    b1 = tajima_constants['b1']
    b2 = tajima_constants['b2']
    c1 = tajima_constants['c1']
    c2 = tajima_constants['c2']
    e1 = tajima_constants['e1']
    e2 = tajima_constants['e2']

    var = (e1 * s + e2 * s * (s - 1))

    if var <= 0:
        return 0.0

    return diff / math.sqrt(var)


def simulate_coalescent(n_samples: int, effective_size: float = 10000,
                       mutation_rate: float = 1e-8) -> Dict[str, Any]:
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
        'n_samples': n_samples,
        'effective_size': effective_size,
        'mutation_rate': mutation_rate,
        'coalescence_times': coalescence_times,
        'total_time': current_time if coalescence_times else 0.0,
        'theta_watterson': watterson_theta(1000, len(coalescence_times)) if coalescence_times else 0.0
    }

    return results


def coalescent_time_to_mrca(n_samples: int, effective_size: float) -> float:
    """Alias for expected_time_to_mrca for backward compatibility."""
    return expected_time_to_mrca(n_samples, effective_size)
