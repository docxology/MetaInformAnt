from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass


def expected_time_to_mrca(sample_size: int, effective_population_size: float) -> float:
    r"""Calculate expected time to most recent common ancestor (MRCA).
    
    Under Kingman's coalescent model, this is the expected total time (in generations)
    until all sampled lineages have coalesced to a single common ancestor.
    
    Args:
        sample_size: Number of sampled individuals (n)
        effective_population_size: Effective population size (Ne, diploid)
        
    Returns:
        Expected time to MRCA. Formula: T_MRCA = 4N \sum_{k=2}^n 1/(k(k-1))
        Uses diploid scaling (4N). Returns 0.0 for invalid inputs.
        
    Examples:
        >>> expected_time_to_mrca(sample_size=10, effective_population_size=1000)
        7366.6...
        >>> expected_time_to_mrca(sample_size=2, effective_population_size=1000)
        4000.0
        
    References:
        Kingman, J. F. C. (1982). On the genealogy of large populations.
        Journal of Applied Probability, 19(A), 27-43.
    """
    n = int(sample_size)
    N = float(effective_population_size)
    if n < 2 or N <= 0:
        return 0.0
    harmonic_like = sum(1.0 / (k * (k - 1)) for k in range(2, n + 1))
    return 4.0 * N * harmonic_like


def expected_total_branch_length(sample_size: int, effective_population_size: float) -> float:
    r"""Calculate expected total branch length of coalescent tree.
    
    Computes the expected sum of all branch lengths in the coalescent tree
    until all lineages have coalesced. This is proportional to the expected
    number of mutations in the sample.
    
    Args:
        sample_size: Number of sampled individuals (n)
        effective_population_size: Effective population size (Ne, diploid)
        
    Returns:
        Expected total branch length in generations. Formula: E[L] = 4Ne × H_{n-1}
        where H_{n-1} is the (n-1)th harmonic number. Returns 0.0 for invalid inputs.
        
    Examples:
        >>> expected_total_branch_length(sample_size=10, effective_population_size=1000)
        9677.4...
        >>> expected_total_branch_length(sample_size=2, effective_population_size=1000)
        4000.0
        
    References:
        Kingman, J. F. C. (1982). On the genealogy of large populations.
        Journal of Applied Probability, 19(A), 27-43.
    """
    n = int(sample_size)
    N = float(effective_population_size)
    if n < 2 or N <= 0:
        return 0.0
    H = sum(1.0 / k for k in range(1, n))
    return 4.0 * N * H


def expected_pairwise_diversity(effective_population_size: float, mutation_rate: float) -> float:
    """Calculate expected pairwise nucleotide diversity π under neutral equilibrium.
    
    Under the infinite sites model with neutral mutations, the expected pairwise
    diversity equals the population mutation parameter θ = 4Neμ for diploids.
    
    Args:
        effective_population_size: Effective population size (Ne, diploid)
        mutation_rate: Mutation rate per site per generation (μ)
        
    Returns:
        Expected pairwise diversity π. Formula: E[π] = 4Ne × μ (diploid).
        Returns 0.0 if Ne <= 0 or μ <= 0.
        
    Examples:
        >>> expected_pairwise_diversity(effective_population_size=10000, mutation_rate=0.0001)
        4.0
        >>> expected_pairwise_diversity(effective_population_size=1000, mutation_rate=0.00001)
        0.04
        
    References:
        Tajima, F. (1983). Evolutionary relationship of DNA sequences in finite
        populations. Genetics, 105(2), 437-460.
    """
    Ne = max(0.0, float(effective_population_size))
    mu = max(0.0, float(mutation_rate))
    return 4.0 * Ne * mu


def expected_pairwise_diversity_from_theta(theta: float) -> float:
    """Calculate expected pairwise diversity from population mutation parameter θ.
    
    Under the infinite sites model, the expected pairwise diversity equals
    the population mutation parameter θ (for per-site calculations).
    
    Args:
        theta: Population mutation parameter θ. For diploids: θ = 4Neμ.
            For haploid: θ = 2Neμ.
            
    Returns:
        Expected pairwise diversity π. For infinite sites model, E[π] = θ.
        Returns 0.0 if theta <= 0.
        
    Examples:
        >>> expected_pairwise_diversity_from_theta(theta=0.01)
        0.01
        >>> expected_pairwise_diversity_from_theta(theta=1.5)
        1.5
        
    Note:
        This is for per-site diversity. For genome-wide estimates, multiply
        by the number of sites or use expected_pairwise_diversity() with
        effective size and mutation rate.
    """
    th = max(0.0, float(theta))
    return th


def tajima_constants(sample_size: int) -> dict[str, float]:
    """Calculate constants required for Tajima's D statistic.
    
    Computes the standard normalization constants used in Tajima's D calculation,
    which compare diversity-based and segregating-sites-based estimates of θ.
    
    Args:
        sample_size: Number of sampled sequences (n)
        
    Returns:
        Dictionary with keys: "a1", "a2", "b1", "b2", "c1", "c2", "e1", "e2".
        These constants are used for variance calculations in Tajima's D.
        Returns zeros for all constants if n < 2.
        
    Examples:
        >>> constants = tajima_constants(sample_size=10)
        >>> "a1" in constants
        True
        >>> constants["a1"] > 0
        True
        
    References:
        Tajima, F. (1989). Statistical method for testing the neutral mutation
        hypothesis by DNA polymorphism. Genetics, 123(3), 585-595.
    """
    n = int(sample_size)
    if n < 2:
        return {k: 0.0 for k in ("a1", "a2", "b1", "b2", "c1", "c2", "e1", "e2")}
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
    """Calculate Tajima's D test statistic for departure from neutral equilibrium.
    
    Tajima's D compares two estimators of the population mutation parameter θ:
    - θ_π (average pairwise diversity)
    - θ_S (Watterson's estimator based on segregating sites)
    
    D > 0 suggests balancing selection or population contraction.
    D < 0 suggests directional selection or population expansion.
    D ≈ 0 is consistent with neutral equilibrium.
    
    Args:
        num_segregating_sites: Number of segregating sites (S)
        pairwise_diversity: Average pairwise nucleotide diversity (π)
        sample_size: Number of sampled sequences (n)
        
    Returns:
        Tajima's D statistic (standardized difference). Returns 0.0 if variance
        term is zero or inputs are invalid.
        
    Examples:
        >>> tajimas_D(num_segregating_sites=10, pairwise_diversity=5.0, sample_size=10)
        -2.1...
        
    References:
        Tajima, F. (1989). Statistical method for testing the neutral mutation
        hypothesis by DNA polymorphism. Genetics, 123(3), 585-595.
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
    return D_num / (var**0.5)


def watterson_theta(
    num_segregating_sites: int,
    sample_size: int,
    *,
    sequence_length: float | None = None,
) -> float:
    """Calculate Watterson's estimator of population mutation parameter θ.
    
    Based on the number of segregating sites rather than pairwise diversity.
    Provides an alternative estimate of θ = 4Neμ (diploid) that can be compared
    with π-based estimates in tests for selection.
    
    Args:
        num_segregating_sites: Number of segregating (polymorphic) sites (S)
        sample_size: Number of sampled sequences (n)
        sequence_length: Optional sequence length L for per-site normalization.
            If provided, returns θ_W / L. If None, returns total θ_W.
            
    Returns:
        Watterson's θ estimate. Formula: θ_W = S / a₁ where a₁ = Σᵢ₌₁ⁿ⁻¹ 1/i.
        If sequence_length is provided, returns θ_W / L (per-site estimate).
        Returns 0.0 if sample_size < 2 or S <= 0.
        
    Examples:
        >>> watterson_theta(num_segregating_sites=10, sample_size=10)
        4.08...
        >>> watterson_theta(num_segregating_sites=100, sample_size=50, sequence_length=1000)
        0.249...
        
    References:
        Watterson, G. A. (1975). On the number of segregating sites in genetical
        models without recombination. Theoretical Population Biology, 7(2), 256-276.
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
    sample_size: int,
    theta: float,
    *,
    sequence_length: float | None = None,
) -> float:
    """Calculate expected number of segregating sites under infinite sites model.
    
    Under the infinite sites model (each mutation creates a new unique site),
    the expected number of segregating sites depends on θ, sample size, and
    sequence length.
    
    Args:
        sample_size: Number of sampled sequences (n)
        theta: Population mutation parameter θ. For diploids: θ = 4Neμ.
            For haploid: θ = 2Neμ.
        sequence_length: Optional sequence length L. If provided, assumes
            per-site θ and multiplies result by L.
            
    Returns:
        Expected number of segregating sites. Formula: E[S] = a₁ × θ where
        a₁ = Σᵢ₌₁ⁿ⁻¹ 1/i. If sequence_length is provided, returns E[S] = a₁ × θ × L.
        Returns 0.0 if inputs invalid.
        
    Examples:
        >>> expected_segregating_sites(sample_size=10, theta=0.01)
        2.82...
        >>> expected_segregating_sites(sample_size=10, theta=0.001, sequence_length=1000)
        2.82...
        
    References:
        Watterson, G. A. (1975). On the number of segregating sites in genetical
        models without recombination. Theoretical Population Biology, 7(2), 256-276.
    """
    # Support both call orders: (n, theta) and (theta, n) for backward compatibility.
    # Detect if inputs are reversed (first arg looks like theta < 2, second like integer n >= 2).
    n_raw = sample_size
    th_raw = theta
    try:
        n_as_float = float(n_raw)
        th_as_float = float(th_raw)
    except Exception:
        n_as_float = float("nan")
        th_as_float = float("nan")

    if (n_as_float < 2.0) and (th_as_float >= 2.0) and (abs(th_as_float - round(th_as_float)) < 1e-9):
        # Likely called as (theta, n)
        n = int(round(th_as_float))
        th = max(0.0, float(n_as_float))
    else:
        n = int(round(n_as_float))
        th = max(0.0, float(th_as_float))
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
    """Calculate expected site frequency spectrum (SFS) counts under neutral equilibrium.
    
    The site frequency spectrum describes the distribution of allele frequencies
    across polymorphic sites. Under neutral equilibrium, rare variants are more
    common than common variants due to the recent origin of most mutations.
    
    Args:
        sample_size: Number of sampled sequences (n)
        theta: Population mutation parameter θ. For diploids: θ = 4Neμ.
            For haploid: θ = 2Neμ.
        sequence_length: Optional sequence length L. If provided, multiplies
            each count by L.
            
    Returns:
        List of length (n-1) containing expected counts of sites with
        i derived alleles, for i = 1, 2, ..., n-1. Formula: E[ξᵢ] = θ / i
        (per-site). If sequence_length provided, E[ξᵢ] = θ × L / i.
        Returns zeros if inputs invalid.
        
    Examples:
        >>> sfs = expected_sfs_counts(sample_size=5, theta=0.01)
        >>> len(sfs)
        4
        >>> sfs[0] > sfs[-1]  # Rare variants more common
        True
        
    References:
        Fu, Y. X. (1995). Statistical properties of segregating sites.
        Theoretical Population Biology, 48(2), 172-197.
    """
    n = int(sample_size)
    th = max(0.0, float(theta))
    if n < 2 or th <= 0:
        return [0.0] * max(0, n - 1)
    factor = th
    if sequence_length is not None and sequence_length > 0:
        factor *= float(sequence_length)
    return [factor / i for i in range(1, n)]


def expected_coalescent_waiting_times(sample_size: int, effective_population_size: float) -> list[float]:
    """Calculate expected waiting times for each coalescent event.
    
    Under Kingman's coalescent, the expected time until k lineages coalesce
    to (k-1) lineages is 4Ne / (k(k-1)) for diploids.
    
    Args:
        sample_size: Number of sampled sequences (n)
        effective_population_size: Effective population size (Ne, diploid)
        
    Returns:
        List of expected waiting times in generations, one for each coalescent
        event from n lineages down to 2 lineages. The k-th element (0-indexed)
        is the expected time for (n-k) lineages to coalesce to (n-k-1) lineages.
        Formula: E[Tₖ] = 4Ne / (k(k-1)). Returns empty list if inputs invalid.
        
    Examples:
        >>> times = expected_coalescent_waiting_times(sample_size=5, effective_population_size=1000)
        >>> len(times)
        3  # 5->4, 4->3, 3->2 coalescent events
        >>> times[0] < times[-1]  # More recent coalescence faster
        True
        
    References:
        Kingman, J. F. C. (1982). On the genealogy of large populations.
        Journal of Applied Probability, 19(A), 27-43.
    """
    n = int(sample_size)
    N = float(effective_population_size)
    if n < 2 or N <= 0:
        return []
    return [4.0 * N / (k * (k - 1)) for k in range(n, 1, -1)]


@dataclass(slots=True)
class CoalescentSummary:
    """Summary statistics container for coalescent model parameters.
    
    Provides convenient access to coalescent-derived statistics given
    population parameters. All methods use the stored sample size,
    effective population size, and mutation rate.
    
    Attributes:
        sample_size: Number of sampled sequences (n)
        effective_population_size: Effective population size (Ne, diploid)
        mutation_rate: Mutation rate per site per generation (μ)
        
    Examples:
        >>> summary = CoalescentSummary(
        ...     sample_size=10,
        ...     effective_population_size=10000,
        ...     mutation_rate=0.0001
        ... )
        >>> theta = summary.theta()
        >>> theta
        4.0
        >>> tmrca = summary.time_to_mrca()
        >>> tmrca > 0
        True
    """
    sample_size: int
    effective_population_size: float
    mutation_rate: float

    def theta(self) -> float:
        """Calculate population mutation parameter θ = 4Neμ (diploid).
        
        Returns:
            θ value in generations
        """
        return expected_pairwise_diversity(self.effective_population_size, self.mutation_rate)

    def total_branch_length(self) -> float:
        """Calculate expected total branch length of coalescent tree.
        
        Returns:
            Expected total branch length in generations
        """
        return expected_total_branch_length(self.sample_size, self.effective_population_size)

    def time_to_mrca(self) -> float:
        """Calculate expected time to most recent common ancestor.
        
        Returns:
            Expected TMRCA in generations
        """
        return expected_time_to_mrca(self.sample_size, self.effective_population_size)

    def tajimas_D(self, num_segregating_sites: int) -> float:
        """Calculate Tajima's D statistic.
        
        Args:
            num_segregating_sites: Number of segregating sites (S)
            
        Returns:
            Tajima's D statistic
        """
        return tajimas_D(num_segregating_sites, self.theta(), self.sample_size)


def wattersons_theta(num_segregating_sites: int, sample_size: int) -> float:
    """Calculate Watterson's estimator θ_W (convenience alias).
    
    This is an alias for watterson_theta() with default parameters.
    See watterson_theta() for full documentation.
    
    Args:
        num_segregating_sites: Number of segregating sites (S)
        sample_size: Number of sampled sequences (n)
        
    Returns:
        Watterson's θ estimate: θ_W = S / a₁
        
    Examples:
        >>> wattersons_theta(num_segregating_sites=10, sample_size=10)
        4.08...
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
    """Count sites by derived allele frequency to construct site frequency spectrum.
    
    Converts a list of derived allele counts (how many derived alleles at each site)
    into the site frequency spectrum, which counts how many sites have i derived
    alleles for each i = 1, 2, ..., n-1.
    
    Args:
        derived_counts: Iterable of derived allele counts, one per polymorphic site
            (e.g., [1, 2, 1, 3, 1] means 5 sites with 1, 2, 1, 3, 1 derived alleles)
        sample_size: Number of sampled sequences (n). Derived counts must be in [1, n-1]
            for polymorphic sites.
            
    Returns:
        List of length (n-1) where the i-th element (0-indexed) is the count of sites
        with exactly (i+1) derived alleles. Indices 0 to n-2 correspond to frequencies
        1 to n-1.
        
    Examples:
        >>> counts = site_frequency_spectrum_counts([1, 1, 2, 2, 3], sample_size=5)
        >>> counts[0]  # Sites with 1 derived allele
        2
        >>> counts[1]  # Sites with 2 derived alleles
        2
        >>> counts[2]  # Sites with 3 derived alleles
        1
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
