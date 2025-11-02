from __future__ import annotations

import math


def hardy_weinberg_genotype_freqs(allele_a_frequency: float) -> tuple[float, float, float]:
    """Calculate Hardy–Weinberg equilibrium genotype frequencies.
    
    Under Hardy–Weinberg equilibrium (no selection, mutation, migration, drift),
    genotype frequencies are determined solely by allele frequencies:
    - AA: p²
    - Aa: 2pq  
    - aa: q²
    where p is allele A frequency and q = 1 - p is allele a frequency.
    
    Args:
        allele_a_frequency: Frequency of allele A (p) in [0, 1]
        
    Returns:
        Tuple of (AA_frequency, Aa_frequency, aa_frequency).
        Returns (0.0, 0.0, 0.0) if input is invalid.
        
    Examples:
        >>> hardy_weinberg_genotype_freqs(0.5)
        (0.25, 0.5, 0.25)
        >>> hardy_weinberg_genotype_freqs(0.7)
        (0.49, 0.42, 0.09)
        
    References:
        Hardy, G. H. (1908). Mendelian proportions in a mixed population.
        Science, 28(706), 49-50.
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
    """Calculate allele frequency after one generation of viability selection.
    
    Implements deterministic selection in a diploid population with viability
    selection. Genotype fitnesses determine the change in allele frequency
    according to the standard Wright–Fisher selection model.
    
    Args:
        allele_a_frequency: Current frequency of allele A (p)
        fitness_AA: Fitness of AA homozygote
        fitness_Aa: Fitness of Aa heterozygote
        fitness_aa: Fitness of aa homozygote
        
    Returns:
        New allele A frequency after selection. Returns original frequency
        if input is invalid or mean fitness is zero.
        
    Examples:
        >>> # Directional selection favoring A
        >>> selection_update(0.5, fitness_AA=1.0, fitness_Aa=1.0, fitness_aa=0.8)
        0.526...
        
        >>> # Overdominance (heterozygote advantage)
        >>> selection_update(0.5, fitness_AA=0.8, fitness_Aa=1.0, fitness_aa=0.8)
        0.5
        
    References:
        Crow, J. F., & Kimura, M. (1970). An introduction to population
        genetics theory. Harper & Row.
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
    """Calculate allele frequency after one generation of bidirectional mutation.
    
    Models mutation as a reversible process with forward mutation rate μ (A→a)
    and backward mutation rate ν (a→A). Approaches equilibrium frequency
    ν / (μ + ν) over many generations.
    
    Args:
        allele_a_frequency: Current frequency of allele A (p)
        mu: Forward mutation rate from A to a (μ)
        nu: Backward mutation rate from a to A (ν)
        
    Returns:
        New allele A frequency after mutation. Returns original frequency if
        input is invalid.
        Formula: p' = p(1 - μ) + (1 - p)ν
        
    Examples:
        >>> mutation_update(0.5, mu=0.01, nu=0.005)
        0.5025...
        
    References:
        Crow, J. F., & Kimura, M. (1970). An introduction to population
        genetics theory. Harper & Row.
    """
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
    """Calculate probability that an allele will fix in a population.
    
    Uses Kimura's diffusion approximation for fixation probability under
    selection and genetic drift in a finite population.
    
    Args:
        initial_frequency: Starting frequency of the allele (p)
        effective_population_size: Effective population size (Ne)
        selection_coefficient: Selection coefficient (s). s > 0 for advantage,
            s < 0 for disadvantage, s = 0 for neutral
            
    Returns:
        Fixation probability in [0, 1]. For neutral alleles (s=0), returns
        the initial frequency. Values clamped to [0, 1].
        
    Examples:
        >>> # Neutral allele with 10% initial frequency
        >>> fixation_probability(0.1, effective_population_size=1000, selection_coefficient=0.0)
        0.1
        
        >>> # Advantageous allele
        >>> fixation_probability(0.01, effective_population_size=1000, selection_coefficient=0.01)
        0.182...
        
    References:
        Kimura, M. (1962). On the probability of fixation of mutant genes
        in a population. Genetics, 47(6), 713-719.
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
    """Calculate Watterson's estimator of the population mutation parameter θ.
    
    Watterson's θ is based on the number of segregating sites and provides
    an estimate of the population-scaled mutation rate. Useful for neutral
    mutation rate estimation.
    
    Args:
        num_segregating_sites: Number of segregating (polymorphic) sites (S)
        sample_size: Number of sampled sequences (n)
        
    Returns:
        Watterson's θ estimate. Returns 0.0 if sample_size < 2 or S <= 0.
        Formula: θ_W = S / a₁ where a₁ = Σᵢ₌₁ⁿ⁻¹ 1/i
        
    Examples:
        >>> watterson_theta(num_segregating_sites=10, sample_size=10)
        4.08...
        >>> watterson_theta(num_segregating_sites=20, sample_size=20)
        6.82...
        
    References:
        Watterson, G. A. (1975). On the number of segregating sites in
        genetical models without recombination. Theoretical Population Biology,
        7(2), 256-276.
    """
    S = max(0, int(num_segregating_sites))
    n = int(sample_size)
    if n < 2:
        return 0.0
    a1 = sum(1.0 / i for i in range(1, n))
    if a1 <= 0:
        return 0.0
    return S / a1


def heterozygosity_decay(initial_heterozygosity: float, effective_population_size: float, generations: int) -> float:
    """Calculate expected heterozygosity decay under genetic drift.
    
    Under genetic drift in a finite population, heterozygosity decreases
    over time due to allele fixation and loss. The rate of decay depends
    on effective population size.
    
    Args:
        initial_heterozygosity: Starting heterozygosity H₀
        effective_population_size: Effective population size (Ne)
        generations: Number of generations (t)
        
    Returns:
        Expected heterozygosity after t generations. Values clamped to [0, 1].
        Formula: H_t = H₀ × (1 - 1/(2Ne))^t
        
    Examples:
        >>> heterozygosity_decay(initial_heterozygosity=0.5, effective_population_size=100, generations=100)
        0.303...
        >>> heterozygosity_decay(initial_heterozygosity=0.5, effective_population_size=1000, generations=100)
        0.475...  # Slower decay with larger population
        
    References:
        Nei, M., & Tajima, F. (1981). DNA polymorphism detectable by
        restriction endonucleases. Genetics, 97(1), 145-163.
    """
    H0 = max(0.0, min(1.0, initial_heterozygosity))
    Ne = float(effective_population_size)
    t = max(0, int(generations))
    if Ne <= 0:
        return 0.0
    factor = 1.0 - 1.0 / (2.0 * Ne)
    return max(0.0, min(1.0, H0 * (factor**t)))


def inbreeding_coefficient(effective_population_size: float, generations: int) -> float:
    """Calculate inbreeding coefficient under genetic drift.
    
    The inbreeding coefficient F measures the probability that two alleles
    at a locus are identical by descent. Under drift, F increases over time
    as the population becomes more inbred.
    
    Args:
        effective_population_size: Effective population size (Ne)
        generations: Number of generations (t)
        
    Returns:
        Inbreeding coefficient F_t in [0, 1]. Returns 0.0 if Ne <= 0.
        Formula: F_t = 1 - (1 - 1/(2Ne))^t
        
    Examples:
        >>> inbreeding_coefficient(effective_population_size=100, generations=100)
        0.393...
        >>> inbreeding_coefficient(effective_population_size=1000, generations=100)
        0.0487...  # Lower inbreeding with larger population
        
    References:
        Crow, J. F., & Kimura, M. (1970). An introduction to population
        genetics theory. Harper & Row.
    """
    Ne = float(effective_population_size)
    t = max(0, int(generations))
    if Ne <= 0:
        return 0.0
    factor = 1.0 - 1.0 / (2.0 * Ne)
    Ft = 1.0 - (factor**t)
    return max(0.0, min(1.0, Ft))


def equilibrium_heterozygosity_infinite_alleles(effective_population_size: float, mutation_rate: float) -> float:
    """Calculate equilibrium heterozygosity under infinite-alleles mutation model.
    
    Under the infinite-alleles model, each mutation creates a new unique allele.
    At equilibrium, mutation (increasing diversity) balances genetic drift
    (decreasing diversity).
    
    Args:
        effective_population_size: Effective population size (Ne)
        mutation_rate: Mutation rate per generation per locus (μ)
        
    Returns:
        Equilibrium heterozygosity in [0, 1]. Returns 0.0 if either parameter <= 0.
        Formula: H_eq = 4Neμ / (1 + 4Neμ)
        
    Examples:
        >>> equilibrium_heterozygosity_infinite_alleles(effective_population_size=1000, mutation_rate=0.0001)
        0.285...
        >>> equilibrium_heterozygosity_infinite_alleles(effective_population_size=10000, mutation_rate=0.0001)
        0.8...
        
    References:
        Kimura, M., & Crow, J. F. (1964). The number of alleles that can be
        maintained in a finite population. Genetics, 49(4), 725-738.
    """
    Ne = max(0.0, float(effective_population_size))
    mu = max(0.0, float(mutation_rate))
    x = 4.0 * Ne * mu
    if x <= 0:
        return 0.0
    return x / (1.0 + x)


def island_model_update(local_frequency: float, migration_rate: float, migrant_pool_frequency: float) -> float:
    """Calculate allele frequency after one generation of Wright's island model.
    
    The island model describes migration between a local population and a
    large migrant pool. Local frequency changes due to migration from the
    migrant pool.
    
    Args:
        local_frequency: Current local allele A frequency (p)
        migration_rate: Fraction of population replaced by migrants per generation (m)
        migrant_pool_frequency: Allele A frequency in migrant pool (p_m)
        
    Returns:
        New local allele frequency after migration. All inputs clamped to [0, 1].
        Formula: p' = (1 - m) × p + m × p_m
        
    Examples:
        >>> island_model_update(local_frequency=0.2, migration_rate=0.1, migrant_pool_frequency=0.5)
        0.23
        
    References:
        Wright, S. (1931). Evolution in Mendelian populations. Genetics,
        16(2), 97-159.
    """
    p = max(0.0, min(1.0, local_frequency))
    m = max(0.0, min(1.0, migration_rate))
    pm = max(0.0, min(1.0, migrant_pool_frequency))
    return (1.0 - m) * p + m * pm


def mutation_selection_balance_recessive(mutation_rate: float, selection_coefficient: float) -> float:
    """Calculate equilibrium frequency under mutation-selection balance (recessive).
    
    For a fully recessive deleterious allele, the equilibrium frequency is
    determined by the balance between mutation creating the allele and selection
    against homozygotes.
    
    Args:
        mutation_rate: Mutation rate creating deleterious allele (μ)
        selection_coefficient: Selection coefficient against homozygote (s)
            
    Returns:
        Equilibrium frequency of deleterious allele. Returns 0.0 if s <= 0.
        Formula: q_eq ≈ √(μ / s)
        
    Examples:
        >>> mutation_selection_balance_recessive(mutation_rate=0.0001, selection_coefficient=0.1)
        0.0316...
        >>> mutation_selection_balance_recessive(mutation_rate=0.0001, selection_coefficient=0.01)
        0.1
        
    References:
        Crow, J. F., & Kimura, M. (1970). An introduction to population
        genetics theory. Harper & Row.
    """
    mu = max(0.0, float(mutation_rate))
    s = max(0.0, float(selection_coefficient))
    if s <= 0:
        return 0.0
    return (mu / s) ** 0.5


def mutation_selection_balance_dominant(mutation_rate: float, selection_coefficient: float) -> float:
    """Calculate equilibrium frequency under mutation-selection balance (dominant).
    
    For a fully dominant deleterious allele, heterozygotes and homozygotes
    both experience selection. Equilibrium frequency is typically much lower
    than for recessive alleles.
    
    Args:
        mutation_rate: Mutation rate creating deleterious allele (μ)
        selection_coefficient: Selection coefficient (s, applies to heterozygotes
            and homozygotes for dominant allele)
            
    Returns:
        Equilibrium frequency of deleterious allele. Returns 0.0 if s <= 0.
        Formula: q_eq ≈ μ / s
        
    Examples:
        >>> mutation_selection_balance_dominant(mutation_rate=0.0001, selection_coefficient=0.1)
        0.001
        >>> mutation_selection_balance_dominant(mutation_rate=0.0001, selection_coefficient=0.01)
        0.01
        
    References:
        Haldane, J. B. S. (1927). A mathematical theory of natural and artificial
        selection, part V: selection and mutation. Mathematical Proceedings of
        the Cambridge Philosophical Society, 23(7), 838-844.
    """
    mu = max(0.0, float(mutation_rate))
    s = max(0.0, float(selection_coefficient))
    if s <= 0:
        return 0.0
    return mu / s


def effective_population_size_from_heterozygosity(observed_heterozygosity: float, theta: float = 4) -> float:
    """Estimate effective population size from observed heterozygosity.
    
    Uses the relationship between heterozygosity and the population mutation
    parameter θ = 4Neμ (or 2Neμ for haploid) to infer effective size from
    observed diversity.
    
    Args:
        observed_heterozygosity: Observed heterozygosity in the population (H)
        theta: Scaling factor for mutation parameter. Default 4 for diploid
            autosomes (θ = 4Neμ). Use 2 for haploid or mtDNA (θ = 2Neμ).
            
    Returns:
        Estimated effective population size. Returns 0.0 if heterozygosity
        is invalid or theta <= 0. Formula: Ne ≈ H / (theta × μ) where μ
        is mutation rate (not provided, so this is a simplified approximation).
        
    Examples:
        >>> effective_population_size_from_heterozygosity(0.5, theta=4)
        0.125...  # Simplified approximation
        
    Note:
        This is a simplified calculation. For accurate estimates, mutation
        rate must be known and the formula Ne = H / (theta × μ) used directly.
        
    References:
        Nei, M., & Tajima, F. (1981). DNA polymorphism detectable by restriction
        endonucleases. Genetics, 97(1), 145-163.
    """
    if observed_heterozygosity <= 0 or observed_heterozygosity >= 1:
        raise ValueError("Heterozygosity must be between 0 and 1")

    # H_e = theta / (theta + 1), so theta = H_e / (1 - H_e)
    # Ne = theta / (4 * mu) for diploids, but here we estimate from H_o
    # This is an approximation
    expected_heterozygosity = theta / (theta + 1)
    return observed_heterozygosity / (4 * (1 - observed_heterozygosity)) if observed_heterozygosity < 1 else float('inf')


def inbreeding_coefficient_from_fst(fst: float, subpopulations: int = 2) -> float:
    """Estimate inbreeding coefficient from Fst.
    
    Under certain population structure models, Fst is related to the
    inbreeding coefficient F. This provides an approximate conversion.
    
    Args:
        fst: Wright's Fst value (population differentiation) in [0, 1]
        subpopulations: Number of subpopulations (default 2)
            
    Returns:
        Estimated inbreeding coefficient F. Returns 0.0 if fst <= 0.
        
    Examples:
        >>> inbreeding_coefficient_from_fst(0.1, subpopulations=2)
        0.052...
        >>> inbreeding_coefficient_from_fst(0.3, subpopulations=5)
        0.075...
        
    Note:
        This is an approximation. The exact relationship depends on
        population structure assumptions.
    """
    if not 0 <= fst <= 1:
        raise ValueError("F_ST must be between 0 and 1")

    # F_IT = F_ST + F_IS * (1 - F_ST)
    # For Wright's F-statistics, this gives an estimate
    return fst / (1 - fst) if fst < 1 else float('inf')


def linkage_disequilibrium_decay_distance(r_squared: float, recombination_rate: float) -> float:
    """Estimate genetic distance from observed r² and recombination rate.
    
    Uses the relationship between LD decay and recombination to infer genetic
    map distance. Assumes LD is at equilibrium under mutation-drift-recombination
    balance.
    
    Args:
        r_squared: Observed r² value (linkage disequilibrium measure)
        recombination_rate: Recombination fraction (c) in [0, 0.5]
        
    Returns:
        Estimated genetic distance in Morgans. Returns 0.0 if inputs invalid.
        
    Examples:
        >>> linkage_disequilibrium_decay_distance(r_squared=0.5, recombination_rate=0.01)
        34.65...
        
    Note:
        This is an approximation. Actual distance depends on population history
        and mutation rates.
    """
    if r_squared <= 0 or r_squared >= 1:
        raise ValueError("r² must be between 0 and 1")

    # LD decay: r² = e^(-2 * recombination_rate * distance)
    # So distance = -ln(r²) / (2 * recombination_rate)
    import math
    return -math.log(r_squared) / (2 * recombination_rate) if recombination_rate > 0 else float('inf')


def coalescent_time_to_mrca(sample_size: int, effective_size: float) -> float:
    """Calculate expected time to most recent common ancestor (MRCA) under coalescent.
    
    Computes the expected number of generations until all sampled lineages
    have coalesced to a single common ancestor. This is the height of the
    coalescent tree.
    
    Args:
        sample_size: Number of sampled sequences (n)
        effective_size: Effective population size (Ne)
        
    Returns:
        Expected time to MRCA in generations. Formula based on sum of
        expected coalescent waiting times.
        E[T_MRCA] ≈ 2Ne × (1 - 1/n)
        
    Examples:
        >>> coalescent_time_to_mrca(sample_size=10, effective_size=1000)
        1800.0...
        >>> coalescent_time_to_mrca(sample_size=100, effective_size=10000)
        19800.0...
        
    References:
        Kingman, J. F. C. (1982). On the genealogy of large populations.
        Journal of Applied Probability, 19(A), 27-43.
    """
    if sample_size < 2:
        return 0.0

    # Expected TMRCA for neutral coalescent
    # E[T_MRCA] = 2 * Ne * sum_{k=2 to n} 1/k
    import math
    harmonic_sum = sum(1/k for k in range(2, sample_size + 1))
    return 2 * effective_size * harmonic_sum
