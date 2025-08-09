from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass


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


def expected_pairwise_diversity_from_theta(theta: float) -> float:
    """Expected pairwise diversity π from population-scaled mutation rate θ.

    Under the standard neutral coalescent with infinite sites, E[π] = θ (per-site).
    """
    th = max(0.0, float(theta))
    return th


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


def watterson_theta(
    num_segregating_sites: int,
    sample_size: int,
    *,
    sequence_length: float | None = None,
) -> float:
    """Watterson's estimator of θ (per-site if sequence_length provided).

    θ_W = S / a1 (total) or θ_W = S / (a1 * L) (per-site) with L=sequence_length.
    Returns 0.0 for invalid inputs.
    """
    S = max(0, int(num_segregating_sites))
    n = int(sample_size)
    if n < 2 or S <= 0:
        return 0.0
    const = tajima_constants(n)
    a1 = const["a1"]
    if a1 <= 0:
        return 0.0
    if sequence_length is not None and sequence_length > 0:
        return float(S) / (a1 * float(sequence_length))
    return float(S) / a1


def expected_segregating_sites(
    theta: float,
    sample_size: int,
    *,
    sequence_length: float | None = None,
) -> float:
    """E[S] under the standard neutral coalescent.

    E[S] = a1 * θ (total) or E[S] = a1 * θ * L (if per-site θ and sequence length L provided).
    Returns 0.0 for invalid inputs.
    """
    th = max(0.0, float(theta))
    n = int(sample_size)
    if n < 2 or th <= 0:
        return 0.0
    a1 = tajima_constants(n)["a1"]
    if sequence_length is not None and sequence_length > 0:
        return a1 * th * float(sequence_length)
    return a1 * th


def expected_sfs_counts(
    sample_size: int,
    theta: float,
    *,
    sequence_length: float | None = None,
) -> list[float]:
    """Expected site frequency spectrum counts under neutrality.

    For derived allele counts i=1..n-1, E[X_i] = θ / i (per-site) or θ L / i (with length L).
    Returns a list of length n-1.
    """
    n = int(sample_size)
    th = max(0.0, float(theta))
    if n < 2 or th <= 0:
        return [0.0] * max(0, n - 1)
    factor = th
    if sequence_length is not None and sequence_length > 0:
        factor *= float(sequence_length)
    return [factor / i for i in range(1, n)]


def expected_coalescent_waiting_times(
    sample_size: int, effective_population_size: float
) -> list[float]:
    """Expected waiting times between coalescent events T_k for k=n..2 (diploid scaling).

    E[T_k] = 4N / (k(k-1)). Returns a list of length n-1 in descending k.
    """
    n = int(sample_size)
    N = float(effective_population_size)
    if n < 2 or N <= 0:
        return []
    return [4.0 * N / (k * (k - 1)) for k in range(n, 1, -1)]

@dataclass(slots=True)
class CoalescentSummary:
    sample_size: int
    effective_population_size: float
    mutation_rate: float

    def theta(self) -> float:
        return expected_pairwise_diversity(self.effective_population_size, self.mutation_rate)

    def total_branch_length(self) -> float:
        return expected_total_branch_length(self.sample_size, self.effective_population_size)

    def time_to_mrca(self) -> float:
        return expected_time_to_mrca(self.sample_size, self.effective_population_size)

    def tajimas_D(self, num_segregating_sites: int) -> float:
        return tajimas_D(num_segregating_sites, self.theta(), self.sample_size)


def wattersons_theta(num_segregating_sites: int, sample_size: int) -> float:
    """Watterson's estimator: theta_W = S / a1.

    Returns 0.0 if invalid inputs.
    """
    S = max(0, int(num_segregating_sites))
    n = int(sample_size)
    if n < 2:
        return 0.0
    a1 = sum(1.0 / i for i in range(1, n))
    if a1 <= 0:
        return 0.0
    return S / a1


def site_frequency_spectrum_counts(derived_counts: Iterable[int], sample_size: int) -> list[int]:
    """Compute folded site-frequency spectrum (SFS) counts for i=1..⌊n/2⌋.

    derived_counts: iterable of derived allele counts per segregating site (0..n)
    Returns a list of length floor(n/2), with minor allele counts aggregated.
    """
    n = int(sample_size)
    if n < 2:
        return []
    half = n // 2
    sfs = [0] * half
    for c in derived_counts:
        k = int(c)
        if 0 < k < n:  # polymorphic
            minor = min(k, n - k)
            if 1 <= minor <= half:
                sfs[minor - 1] += 1
    return sfs

 


