from __future__ import annotations

import math


def hardy_weinberg_genotype_freqs(allele_a_frequency: float) -> tuple[float, float, float]:
    """Return Hardy–Weinberg genotype frequencies (AA, Aa, aa) given allele A frequency.

    Returns (p2, 2pq, q2). If input is invalid, returns zeros.
    """
    p = allele_a_frequency
    if not (0.0 <= p <= 1.0):
        return 0.0, 0.0, 0.0
    q = 1.0 * (1.0 - p)
    return p * p, 2.0 * p * q, q * q


def selection_update(
    allele_a_frequency: float,
    fitness_AA: float,
    fitness_Aa: float,
    fitness_aa: float,
) -> float:
    """One-step deterministic selection update for allele A frequency.

    Uses standard viability selection in a diploid Wright–Fisher model.
    """
    p = allele_a_frequency
    if not (0.0 <= p <= 1.0):
        return p
    q = 1.0 - p
    mean_fitness = (p * p) * fitness_AA + 2.0 * p * q * fitness_Aa + (q * q) * fitness_aa
    if mean_fitness <= 0.0:
        return p
    next_p_numerator = (p * p) * fitness_AA + p * q * fitness_Aa
    return next_p_numerator / mean_fitness


def mutation_update(allele_a_frequency: float, mu: float, nu: float) -> float:
    """One-step mutation update with forward rate mu (A→a) and back rate nu (a→A)."""
    p = allele_a_frequency
    if not (0.0 <= p <= 1.0):
        return p
    mu = max(0.0, mu)
    nu = max(0.0, nu)
    # p' = p(1-mu) + (1-p)nu
    return p * (1.0 - mu) + (1.0 - p) * nu


def fixation_probability(
    initial_frequency: float, effective_population_size: int, selection_coefficient: float = 0.0
) -> float:
    """Approximate fixation probability under selection (Kimura formula).

    For s = 0, returns the initial frequency. For s != 0, uses
    u(p) = (1 - exp(-2 N s p)) / (1 - exp(-2 N s)). This is a common approximation;
    values are clamped to [0, 1].
    """
    p = initial_frequency
    N = max(1, int(effective_population_size))
    s = selection_coefficient
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    if abs(s) < 1e-12:
        return max(0.0, min(1.0, p))
    try:
        numerator = 1.0 - math.exp(-2.0 * N * s * p)
        denominator = 1.0 - math.exp(-2.0 * N * s)
        if abs(denominator) < 1e-18:
            return 1.0 if s > 0 else 0.0
        u = numerator / denominator
    except OverflowError:
        return 1.0 if s > 0 else 0.0
    return max(0.0, min(1.0, u))


def watterson_theta(num_segregating_sites: int, sample_size: int) -> float:
    """Watterson's theta estimate: theta_W = S / a1, a1 = sum_{i=1}^{n-1} 1/i."""
    S = max(0, int(num_segregating_sites))
    n = int(sample_size)
    if n < 2:
        return 0.0
    a1 = sum(1.0 / i for i in range(1, n))
    if a1 <= 0:
        return 0.0
    return S / a1


def heterozygosity_decay(initial_heterozygosity: float, effective_population_size: float, generations: int) -> float:
    """Expected heterozygosity decay under drift: H_t = H_0 (1 - 1/(2Ne))^t.

    Values are clamped to [0, 1].
    """
    H0 = max(0.0, min(1.0, initial_heterozygosity))
    Ne = float(effective_population_size)
    t = max(0, int(generations))
    if Ne <= 0:
        return 0.0
    factor = 1.0 - 1.0 / (2.0 * Ne)
    return max(0.0, min(1.0, H0 * (factor**t)))


def inbreeding_coefficient(effective_population_size: float, generations: int) -> float:
    """Inbreeding coefficient under drift: F_t = 1 - (1 - 1/(2Ne))^t."""
    Ne = float(effective_population_size)
    t = max(0, int(generations))
    if Ne <= 0:
        return 0.0
    factor = 1.0 - 1.0 / (2.0 * Ne)
    Ft = 1.0 - (factor**t)
    return max(0.0, min(1.0, Ft))


def equilibrium_heterozygosity_infinite_alleles(effective_population_size: float, mutation_rate: float) -> float:
    """Equilibrium heterozygosity under infinite-alleles model: Heq = 4Neμ / (1 + 4Neμ)."""
    Ne = max(0.0, float(effective_population_size))
    mu = max(0.0, float(mutation_rate))
    x = 4.0 * Ne * mu
    if x <= 0:
        return 0.0
    return x / (1.0 + x)


def island_model_update(local_frequency: float, migration_rate: float, migrant_pool_frequency: float) -> float:
    """One generation of Wright's island model: p' = (1-m) p + m p_m.

    All inputs are clamped to valid ranges [0, 1] where applicable.
    """
    p = max(0.0, min(1.0, local_frequency))
    m = max(0.0, min(1.0, migration_rate))
    pm = max(0.0, min(1.0, migrant_pool_frequency))
    return (1.0 - m) * p + m * pm


def mutation_selection_balance_recessive(mutation_rate: float, selection_coefficient: float) -> float:
    """Mutation-selection balance for fully recessive deleterious allele: q ≈ sqrt(μ / s)."""
    mu = max(0.0, float(mutation_rate))
    s = max(0.0, float(selection_coefficient))
    if s <= 0:
        return 0.0
    return (mu / s) ** 0.5


def mutation_selection_balance_dominant(mutation_rate: float, selection_coefficient: float) -> float:
    """Mutation-selection balance for fully dominant deleterious allele: q ≈ μ / s."""
    mu = max(0.0, float(mutation_rate))
    s = max(0.0, float(selection_coefficient))
    if s <= 0:
        return 0.0
    return mu / s


def effective_population_size_from_heterozygosity(observed_heterozygosity: float, theta: float = 4) -> float:
    """Estimate effective population size from observed heterozygosity.

    Args:
        observed_heterozygosity: Observed heterozygosity (H_o)
        theta: Theta parameter (default 4 for diploids)

    Returns:
        Estimated effective population size
    """
    if observed_heterozygosity <= 0 or observed_heterozygosity >= 1:
        raise ValueError("Heterozygosity must be between 0 and 1")

    # H_e = theta / (theta + 1), so theta = H_e / (1 - H_e)
    # Ne = theta / (4 * mu) for diploids, but here we estimate from H_o
    # This is an approximation
    expected_heterozygosity = theta / (theta + 1)
    return observed_heterozygosity / (4 * (1 - observed_heterozygosity)) if observed_heterozygosity < 1 else float('inf')


def inbreeding_coefficient_from_fst(fst: float, subpopulations: int = 2) -> float:
    """Estimate inbreeding coefficient from F_ST.

    Args:
        fst: F_ST value
        subpopulations: Number of subpopulations

    Returns:
        Estimated inbreeding coefficient
    """
    if not 0 <= fst <= 1:
        raise ValueError("F_ST must be between 0 and 1")

    # F_IT = F_ST + F_IS * (1 - F_ST)
    # For Wright's F-statistics, this gives an estimate
    return fst / (1 - fst) if fst < 1 else float('inf')


def linkage_disequilibrium_decay_distance(r_squared: float, recombination_rate: float) -> float:
    """Estimate physical distance over which LD decays to a given r².

    Args:
        r_squared: Target r² value
        recombination_rate: Recombination rate per base pair per generation

    Returns:
        Physical distance in base pairs
    """
    if r_squared <= 0 or r_squared >= 1:
        raise ValueError("r² must be between 0 and 1")

    # LD decay: r² = e^(-2 * recombination_rate * distance)
    # So distance = -ln(r²) / (2 * recombination_rate)
    import math
    return -math.log(r_squared) / (2 * recombination_rate) if recombination_rate > 0 else float('inf')


def coalescent_time_to_mrca(sample_size: int, effective_size: float) -> float:
    """Calculate expected time to most recent common ancestor.

    Args:
        sample_size: Number of samples
        effective_size: Effective population size

    Returns:
        Expected TMRCA in generations
    """
    if sample_size < 2:
        return 0.0

    # Expected TMRCA for neutral coalescent
    # E[T_MRCA] = 2 * Ne * sum_{k=2 to n} 1/k
    import math
    harmonic_sum = sum(1/k for k in range(2, sample_size + 1))
    return 2 * effective_size * harmonic_sum
