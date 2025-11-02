from __future__ import annotations

from typing import Iterable


def harmonic_mean_effective_size(census_sizes: Iterable[float]) -> float:
    """Calculate effective population size across multiple generations.
    
    Uses harmonic mean of census sizes across generations. The harmonic mean
    gives more weight to smaller population sizes, reflecting their greater
    impact on genetic drift.
    
    Args:
        census_sizes: Iterable of census population sizes across generations
        
    Returns:
        Harmonic mean effective size. Returns 0.0 if any size is zero or
        input is empty.
        Formula: Ne = n / sum(1/N_i) where n is number of generations.
        
    Examples:
        >>> harmonic_mean_effective_size([1000, 500, 2000, 800])
        800.0...
        >>> harmonic_mean_effective_size([100, 100, 100])
        100.0
        
    References:
        Crow, J. F., & Kimura, M. (1970). An introduction to population
        genetics theory. Harper & Row.
    """
    Ns = [max(0.0, float(N)) for N in census_sizes]
    if not Ns or any(N == 0 for N in Ns):
        return 0.0
    return len(Ns) / sum(1.0 / N for N in Ns)


def effective_size_sex_ratio(num_males: float, num_females: float) -> float:
    """Calculate effective population size with unequal sex ratio.
    
    When the number of breeding males and females differ, effective size
    is reduced compared to census size. This accounts for variance in
    reproductive success due to sex ratio imbalance.
    
    Args:
        num_males: Number of breeding males (Nm)
        num_females: Number of breeding females (Nf)
        
    Returns:
        Effective population size. Returns 0.0 if both are zero.
        Formula: Ne = 4 × Nm × Nf / (Nm + Nf)
        
    Examples:
        >>> effective_size_sex_ratio(num_males=100, num_females=100)
        200.0
        >>> effective_size_sex_ratio(num_males=10, num_females=90)
        36.0
        
    References:
        Wright, S. (1931). Evolution in Mendelian populations. Genetics,
        16(2), 97-159.
    """
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
