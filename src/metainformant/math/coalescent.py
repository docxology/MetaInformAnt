from __future__ import annotations

import math


def expected_time_to_mrca(sample_size: int, effective_population_size: float) -> float:
    r"""Expected total time to MRCA for a sample under Kingman's coalescent.

    T_MRCA = 4N \sum_{k=2}^n 1/(k(k-1)). Uses diploid scaling (4N).
    Returns 0 for invalid inputs.
    """
    n = int(sample_size)
    N = float(effective_population_size)
    if n < 2 or N <= 0:
        return 0.0
    harmonic_like = sum(1.0 / (k * (k - 1)) for k in range(2, n + 1))
    return 4.0 * N * harmonic_like


def expected_total_branch_length(sample_size: int, effective_population_size: float) -> float:
    r"""Expected total tree length under Kingman's coalescent (diploid scaling).

    E[L] = 4N \sum_{k=2}^n 1/(k-1) = 4N H_{n-1}.
    """
    n = int(sample_size)
    N = float(effective_population_size)
    if n < 2 or N <= 0:
        return 0.0
    H = sum(1.0 / k for k in range(1, n))
    return 4.0 * N * H



def expected_pairwise_diversity(effective_population_size: float, mutation_rate: float) -> float:
    """Expected pairwise nucleotide diversity π under neutral equilibrium: π ≈ 4Neμ (diploid)."""
    Ne = max(0.0, float(effective_population_size))
    mu = max(0.0, float(mutation_rate))
    return 4.0 * Ne * mu


def tajima_constants(sample_size: int) -> dict[str, float]:
    """Return standard constants used in Tajima's D for a given sample size n.

    Includes a1, a2, b1, b2, c1, c2, e1, e2 following Tajima (1989).
    """
    n = int(sample_size)
    if n < 2:
        return {k: 0.0 for k in ("a1","a2","b1","b2","c1","c2","e1","e2")}
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i * i) for i in range(1, n))
    b1 = (n + 1) / (3.0 * (n - 1))
    b2 = (2.0 * (n * n + n + 3)) / (9.0 * n * (n - 1))
    c1 = b1 - (1.0 / a1)
    c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1))
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    return {"a1": a1, "a2": a2, "b1": b1, "b2": b2, "c1": c1, "c2": c2, "e1": e1, "e2": e2}


def tajimas_D(num_segregating_sites: int, pairwise_diversity: float, sample_size: int) -> float:
    """Compute Tajima's D given S, π, and n using standard normalizing constants.

    If the variance term is zero or inputs invalid, returns 0.0.
    """
    S = max(0, int(num_segregating_sites))
    pi = max(0.0, float(pairwise_diversity))
    const = tajima_constants(sample_size)
    a1 = const["a1"]
    if a1 <= 0:
        return 0.0
    D_num = pi - (S / a1)
    var = const["e1"] * S + const["e2"] * S * (S - 1)
    if var <= 0:
        return 0.0
    return D_num / (var ** 0.5)


