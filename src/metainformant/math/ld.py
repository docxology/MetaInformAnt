from __future__ import annotations

from typing import Tuple


def ld_coefficients(pA: float, pa: float, pB: float, pb: float, haplotype_pAB: float) -> Tuple[float, float]:
    """Compute linkage disequilibrium D and normalized D' given haplotype frequencies.
    
    Calculates the linkage disequilibrium (LD) coefficient D, which measures
    the deviation from random association between alleles at two loci. Also
    computes D', the normalized LD coefficient adjusted for allele frequencies.
    
    Args:
        pA: Frequency of allele A at first locus
        pa: Frequency of allele a at first locus (should satisfy pA + pa = 1)
        pB: Frequency of allele B at second locus
        pb: Frequency of allele b at second locus (should satisfy pB + pb = 1)
        haplotype_pAB: Observed frequency of AB haplotype
        
    Returns:
        Tuple of (D, D_prime):
        - D: Linkage disequilibrium coefficient. D = pAB - pA × pB.
          Positive D indicates AB and ab haplotypes are more common than expected.
        - D_prime: Normalized LD coefficient in [-1, 1]. D' = D / D_max where
          D_max is the maximum possible |D| given allele frequencies.
        
        Returns (0.0, 0.0) if inputs violate constraints or are invalid.
        
    Examples:
        >>> ld_coefficients(pA=0.6, pa=0.4, pB=0.7, pb=0.3, haplotype_pAB=0.5)
        (0.08, 0.333...)  # Positive LD
        >>> ld_coefficients(pA=0.5, pa=0.5, pB=0.5, pb=0.5, haplotype_pAB=0.25)
        (0.0, 0.0)  # Linkage equilibrium
        
    References:
        Lewontin, R. C. (1964). The interaction of selection and linkage.
        I. General considerations; heterotic models. Genetics, 49(1), 49-67.
    """
    if not (0 <= pA <= 1 and 0 <= pa <= 1 and abs(pA + pa - 1.0) < 1e-8):
        return 0.0, 0.0
    if not (0 <= pB <= 1 and 0 <= pb <= 1 and abs(pB + pb - 1.0) < 1e-8):
        return 0.0, 0.0
    if not (0.0 <= haplotype_pAB <= 1.0):
        return 0.0, 0.0
    D = haplotype_pAB - pA * pB
    Dmax = min(pA * pb, pa * pB) if D >= 0 else min(pA * pB, pa * pb)
    if Dmax <= 0:
        return D, 0.0
    return D, D / Dmax


def r_squared(pA: float, pa: float, pB: float, pb: float, haplotype_pAB: float) -> float:
    """Calculate r², the squared correlation coefficient between two biallelic loci.
    
    r² is a widely used measure of linkage disequilibrium that ranges from 0
    (no LD, alleles are independent) to 1 (complete LD, alleles are perfectly
    correlated). It is the square of the correlation between allele presence
    at the two loci.
    
    Args:
        pA: Frequency of allele A at first locus
        pa: Frequency of allele a at first locus
        pB: Frequency of allele B at second locus
        pb: Frequency of allele b at second locus
        haplotype_pAB: Observed frequency of AB haplotype
        
    Returns:
        r² value in [0, 1]. Returns 0.0 if denominator is zero or inputs invalid.
        Formula: r² = D² / (pA × pa × pB × pb)
        
    Examples:
        >>> r_squared(pA=0.6, pa=0.4, pB=0.7, pb=0.3, haplotype_pAB=0.5)
        0.095...
        >>> r_squared(pA=0.5, pa=0.5, pB=0.5, pb=0.5, haplotype_pAB=0.25)
        0.0  # No LD
        
    References:
        Hill, W. G., & Robertson, A. (1968). Linkage disequilibrium in finite
        populations. Theoretical and Applied Genetics, 38(6), 226-231.
    """
    D, _ = ld_coefficients(pA, pa, pB, pb, haplotype_pAB)
    denom = pA * pa * pB * pb
    if denom <= 0:
        return 0.0
    return (D * D) / denom


def ld_decay_r2(r2_initial: float, recombination_rate: float, generations: int) -> float:
    """Calculate expected r² decay over generations due to recombination.
    
    Models the decay of linkage disequilibrium as recombination breaks down
    haplotype associations. Under random mating, LD decreases exponentially
    with time.
    
    Args:
        r2_initial: Initial r² value (LD at generation 0)
        recombination_rate: Recombination fraction per generation (c) in [0, 0.5]
        generations: Number of generations (t)
        
    Returns:
        Expected r² after t generations. Values clamped to [0, 1].
        Formula: r²_t ≈ r²_0 × (1 - c)^(2t)
        
    Examples:
        >>> ld_decay_r2(r2_initial=0.8, recombination_rate=0.01, generations=50)
        0.296...
        >>> ld_decay_r2(r2_initial=0.8, recombination_rate=0.1, generations=10)
        0.010...  # Faster decay with higher recombination
        
    References:
        Sved, J. A. (1971). Linkage disequilibrium and homozygosity of
        chromosome segments in finite populations. Theoretical Population Biology,
        2(2), 125-141.
    """
    r2 = max(0.0, min(1.0, r2_initial))
    c = max(0.0, min(1.0, recombination_rate))
    t = max(0, generations)
    return r2 * ((1.0 - c) ** (2 * t))


def haldane_d_to_c(map_distance_morgans: float) -> float:
    """Convert genetic map distance to recombination fraction using Haldane's mapping function.
    
    Haldane's function assumes no interference between crossovers. Suitable for
    chromosomes where crossovers occur independently.
    
    Args:
        map_distance_morgans: Genetic map distance in Morgans (d)
        
    Returns:
        Recombination fraction in [0, 0.5]. Formula: c = 0.5 × (1 - exp(-2d))
        
    Examples:
        >>> haldane_d_to_c(0.1)  # 10 cM
        0.090...
        >>> haldane_d_to_c(0.5)  # 50 cM (approaches 0.5)
        0.316...
        
    References:
        Haldane, J. B. S. (1919). The combination of linkage values and the
        calculation of distances between the loci of linked factors.
        Journal of Genetics, 8(4), 299-309.
    """
    d = max(0.0, float(map_distance_morgans))
    from math import exp

    c = 0.5 * (1.0 - exp(-2.0 * d))
    return max(0.0, min(0.5, c))


def haldane_c_to_d(recombination_fraction: float) -> float:
    """Convert recombination fraction to map distance using inverse Haldane function.
    
    Inverse of Haldane's mapping function. Converts observed recombination
    rates to genetic map distances assuming no interference.
    
    Args:
        recombination_fraction: Recombination fraction (c) in [0, 0.5]
        
    Returns:
        Genetic map distance in Morgans. Returns inf if c >= 0.5 (unlinked).
        Formula: d = -0.5 × ln(1 - 2c)
        
    Examples:
        >>> haldane_c_to_d(0.1)
        0.111...
        >>> haldane_c_to_d(0.25)
        0.346...
        >>> haldane_c_to_d(0.5)
        inf  # Unlinked loci
    """
    c = max(0.0, min(0.5, float(recombination_fraction)))
    from math import log

    if c >= 0.5:
        return float("inf")
    return -0.5 * log(1.0 - 2.0 * c)


def kosambi_d_to_c(map_distance_morgans: float) -> float:
    """Convert genetic map distance to recombination fraction using Kosambi's mapping function.
    
    Kosambi's function accounts for crossover interference, where one crossover
    reduces the probability of nearby crossovers. More realistic for many
    chromosomes than Haldane's function.
    
    Args:
        map_distance_morgans: Genetic map distance in Morgans (d)
        
    Returns:
        Recombination fraction in [0, 0.5]. Formula: c = 0.5 × tanh(2d)
        
    Examples:
        >>> kosambi_d_to_c(0.1)  # 10 cM
        0.099...
        >>> kosambi_d_to_c(0.5)  # 50 cM
        0.462...
        
    References:
        Kosambi, D. D. (1944). The estimation of map distances from
        recombination values. Annals of Eugenics, 12(1), 172-175.
    """
    d = max(0.0, float(map_distance_morgans))
    from math import tanh

    c = 0.5 * tanh(2.0 * d)
    return max(0.0, min(0.5, c))


def kosambi_c_to_d(recombination_fraction: float) -> float:
    """Convert recombination fraction to map distance using inverse Kosambi function.
    
    Inverse of Kosambi's mapping function. Accounts for crossover interference
    when converting recombination rates to genetic map distances.
    
    Args:
        recombination_fraction: Recombination fraction (c) in [0, 0.5]
        
    Returns:
        Genetic map distance in Morgans. Returns inf if c >= 0.5 (unlinked).
        Formula: d = 0.25 × ln((1 + 2c) / (1 - 2c))
        
    Examples:
        >>> kosambi_c_to_d(0.1)
        0.100...
        >>> kosambi_c_to_d(0.25)
        0.549...
    """
    c = max(0.0, min(0.5, float(recombination_fraction)))
    from math import log

    if c >= 0.5:
        return float("inf")
    return 0.25 * log((1.0 + 2.0 * c) / (1.0 - 2.0 * c))


def expected_r2_from_Ne_c(effective_population_size: float, recombination_fraction: float) -> float:
    """Calculate expected r² under mutation-drift-recombination balance.
    
    Under neutral evolution, the expected LD between unlinked loci depends
    on effective population size and recombination rate. Higher Ne or higher
    recombination leads to lower expected LD.
    
    Args:
        effective_population_size: Effective population size (Ne)
        recombination_fraction: Recombination fraction (c) in [0, 0.5]
        
    Returns:
        Expected r² value in [0, 1]. Returns 0.0 if Ne <= 0.
        Formula: E[r²] ≈ 1 / (1 + 4Ne × c)
        
    Examples:
        >>> expected_r2_from_Ne_c(effective_population_size=1000, recombination_fraction=0.01)
        0.961...
        >>> expected_r2_from_Ne_c(effective_population_size=10000, recombination_fraction=0.1)
        0.024...  # Much lower LD with larger Ne and higher recombination
        
    References:
        Hill, W. G., & Robertson, A. (1968). Linkage disequilibrium in finite
        populations. Theoretical and Applied Genetics, 38(6), 226-231.
    """
    Ne = max(0.0, float(effective_population_size))
    c = max(0.0, min(0.5, float(recombination_fraction)))
    denom = 1.0 + 4.0 * Ne * c
    if denom <= 0:
        return 0.0
    return 1.0 / denom
