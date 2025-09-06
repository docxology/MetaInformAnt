from __future__ import annotations

from typing import Iterable


def harmonic_mean_effective_size(census_sizes: Iterable[float]) -> float:
    """Effective size across generations via harmonic mean of census sizes."""
    Ns = [max(0.0, float(N)) for N in census_sizes]
    if not Ns or any(N == 0 for N in Ns):
        return 0.0
    return len(Ns) / sum(1.0 / N for N in Ns)


def effective_size_sex_ratio(num_males: float, num_females: float) -> float:
    """Ne under unequal sex ratio: Ne = 4 Nm Nf / (Nm + Nf)."""
    Nm = max(0.0, num_males)
    Nf = max(0.0, num_females)
    if Nm + Nf <= 0:
        return 0.0
    return (4.0 * Nm * Nf) / (Nm + Nf)


def effective_size_from_family_size_variance(census_size: float, variance_offspring_number: float) -> float:
    """Crow and Denniston approximation: Ne ≈ (4N - 2) / (Vk + 2).

    N is census size (diploid), Vk is variance in family size. Guards for non-negative values.
    """
    N = max(0.0, float(census_size))
    Vk = max(0.0, float(variance_offspring_number))
    denom = Vk + 2.0
    if denom <= 0:
        return 0.0
    return (4.0 * N - 2.0) / denom
