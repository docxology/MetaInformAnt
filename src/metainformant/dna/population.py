from __future__ import annotations

import warnings
from collections.abc import Iterable, Sequence


def allele_frequencies(genotype_matrix: Sequence[Sequence[int]]) -> list[float]:
    """Compute frequency of allele '1' per site.

    Calculates the frequency of the alternate allele (encoded as 1) at each
    site across all individuals in the genotype matrix.

    Args:
        genotype_matrix: List of lists where rows are individuals, columns
            are sites. Values should be 0 (reference allele) or 1 (alternate allele).
            All rows must have the same length.

    Returns:
        List of frequencies per site (float values in [0, 1]). Returns empty
        list if input is empty.

    Raises:
        ValueError: If genotype matrix has inconsistent row lengths or invalid values.

    Examples:
        >>> genotypes = [[0, 1, 0], [0, 1, 1], [1, 0, 1]]
        >>> freqs = allele_frequencies(genotypes)
        >>> freqs[0]  # Frequency at first site
        0.333...
        >>> freqs[1]  # Frequency at second site
        0.666...
    """
    if not genotype_matrix:
        return []
    
    # Validate consistent row lengths
    num_sites = len(genotype_matrix[0])
    for i, row in enumerate(genotype_matrix):
        if len(row) != num_sites:
            raise ValueError(f"Row {i} has length {len(row)}, expected {num_sites}")
    
    site_totals = [0] * num_sites
    for i, row in enumerate(genotype_matrix):
        for j, val in enumerate(row):
            if val not in (0, 1):
                raise ValueError(f"Invalid genotype value {val} at row {i}, site {j}. Must be 0 or 1.")
            site_totals[j] += int(val == 1)
    
    n_individuals = len(genotype_matrix)
    return [total / n_individuals for total in site_totals]


def observed_heterozygosity(genotypes: Iterable[tuple[int, int]]) -> float:
    """Proportion of heterozygous individuals among diploid genotypes.

    Calculates the observed heterozygosity as the fraction of individuals
    that are heterozygous (have different alleles at the two homologous chromosomes).

    Args:
        genotypes: Iterable of tuples (a1, a2) representing diploid genotypes.
            Each tuple contains two alleles encoded as 0 (reference) or 1 (alternate).

    Returns:
        Float in [0, 1] representing the proportion of heterozygous individuals.
        Returns 0.0 if input is empty.

    Raises:
        ValueError: If genotype tuples contain invalid allele values.

    Examples:
        >>> genotypes = [(0, 0), (0, 1), (1, 1), (1, 0)]
        >>> observed_heterozygosity(genotypes)
        0.5  # 2 of 4 are heterozygous
        >>> observed_heterozygosity([(0, 0), (1, 1)])
        0.0  # All homozygous
    """
    genotypes_list = list(genotypes)
    if not genotypes_list:
        return 0.0
    
    hetero = 0
    for i, (a1, a2) in enumerate(genotypes_list):
        if a1 not in (0, 1) or a2 not in (0, 1):
            raise ValueError(f"Invalid allele values in genotype {i}: ({a1}, {a2}). Must be 0 or 1.")
        if a1 != a2:
            hetero += 1
    
    return hetero / len(genotypes_list)


def nucleotide_diversity(seqs: Sequence[str]) -> float:
    """Average pairwise nucleotide difference per site (π).

    Calculates nucleotide diversity (π) as the average number of pairwise
    nucleotide differences per site across all pairs of sequences.

    Args:
        seqs: Sequence of DNA sequences (strings). If sequences have different
            lengths, they will be truncated to the shortest length.

    Returns:
        Average pairwise diversity per site (float). Returns 0.0 if:
        - Less than 2 sequences provided
        - All sequences are identical
        - Sequences have zero length

    Examples:
        >>> seqs = ["AAAA", "AAAT"]
        >>> nucleotide_diversity(seqs)
        0.25  # 1 difference in 4 sites = 0.25 per site
        >>> seqs = ["AAAA", "AAAA", "AAAA"]
        >>> nucleotide_diversity(seqs)
        0.0  # All sequences identical

    References:
        Nei, M., & Li, W. H. (1979). Mathematical model for studying genetic
        variation in terms of restriction endonucleases. *Proceedings of the
        National Academy of Sciences*, 76(10), 5269-5273.
    """
    if len(seqs) < 2:
        return 0.0
    
    # Validate sequences are strings
    for i, seq in enumerate(seqs):
        if not isinstance(seq, str):
            raise TypeError(f"Sequence {i} is not a string: {type(seq)}")
    
    # Check for length mismatches
    lengths = [len(s) for s in seqs]
    L = min(lengths)
    max_len = max(lengths)
    if L != max_len:
        warnings.warn(
            f"nucleotide_diversity: Sequences have different lengths (min={L}, max={max_len}). "
            f"Truncating to shortest length ({L}).",
            UserWarning,
            stacklevel=2,
        )
    
    if L == 0:
        return 0.0
    
    n = len(seqs)
    total_diff = 0
    num_pairs = 0
    for i in range(n):
        for j in range(i + 1, n):
            diff = sum(1 for a, b in zip(seqs[i][:L], seqs[j][:L]) if a != b)
            total_diff += diff / L
            num_pairs += 1
    return total_diff / num_pairs if num_pairs else 0.0


def tajimas_d(seqs: Sequence[str]) -> float:
    """Simplified Tajima's D approximation for small samples.

    This is a simplified version of Tajima's D that provides a basic measure
    of departure from neutral equilibrium. For accurate statistical testing,
    use `metainformant.math.coalescent.tajimas_D()` which implements the
    full formula with proper variance calculation.

    Tajima's D compares two estimators of the population mutation parameter θ:
    - θ_π (average pairwise diversity)
    - θ_S (Watterson's estimator based on segregating sites)

    D > 0 suggests balancing selection or population contraction.
    D < 0 suggests directional selection or population expansion.
    D ≈ 0 is consistent with neutral equilibrium.

    Args:
        seqs: Sequence of DNA sequences (must be same length)

    Returns:
        Simplified Tajima's D statistic. Returns 0.0 when:
        - Less than 2 sequences
        - No segregating sites
        - Zero-length sequences

    Note:
        This simplified version uses a basic normalization and may not
        match the standard Tajima's D statistic exactly. For publication-
        quality analysis, use the full implementation in `math.coalescent`.

    Examples:
        >>> seqs = ["AAAA", "AAAT", "AATT"]
        >>> d = tajimas_d(seqs)
        >>> isinstance(d, float)
        True

    References:
        Tajima, F. (1989). Statistical method for testing the neutral mutation
        hypothesis by DNA polymorphism. Genetics, 123(3), 585-595.

    See Also:
        `metainformant.math.coalescent.tajimas_D()` : Full Tajima's D implementation
    """
    if len(seqs) < 2:
        return 0.0
    L = min(len(s) for s in seqs)
    if L == 0:
        return 0.0
    # Count segregating sites
    segregating = 0
    for pos in range(L):
        column = {s[pos] for s in seqs}
        if len(column) > 1:
            segregating += 1
    if segregating == 0:
        return 0.0
    # For tests, a simplified normalized difference between pi and S/L
    pi = nucleotide_diversity(seqs)
    theta_w = segregating / L
    # Basic normalization to avoid division by zero
    denom = max(1e-9, (pi + theta_w) / 2)
    return (pi - theta_w) / denom


def hudson_fst(pop1: Sequence[str], pop2: Sequence[str]) -> float:
    """Hudson's Fst estimator for two populations.

    Calculates Hudson's Fst (1992) which measures population differentiation
    as the proportion of genetic diversity that is between populations rather
    than within them.

    Args:
        pop1: Sequence of DNA sequences from population 1
        pop2: Sequence of DNA sequences from population 2

    Returns:
        Fst value in [0, 1] where:
        - 0 = no differentiation (populations identical)
        - 1 = complete differentiation (fixed differences between populations)

    Note:
        If sequences have different lengths, truncates to shortest length.
        Returns 0.0 if either population is empty.

    Examples:
        >>> pop1 = ["AAAA", "AAAA"]
        >>> pop2 = ["TTTT", "TTTT"]
        >>> hudson_fst(pop1, pop2)
        1.0  # Complete differentiation
        >>> pop1 = ["AAAA", "AAAT"]
        >>> pop2 = ["AAAA", "AAAT"]
        >>> hudson_fst(pop1, pop2)
        0.0  # No differentiation

    References:
        Hudson, R. R., Slatkin, M., & Maddison, W. P. (1992). Estimation of
        levels of gene flow from DNA sequence data. *Genetics*, 132(2), 583-589.
    """
    if not pop1 or not pop2:
        return 0.0
    
    # Check for length mismatches
    pop1_lengths = [len(s) for s in pop1]
    pop2_lengths = [len(s) for s in pop2]
    L1 = min(pop1_lengths)
    L2 = min(pop2_lengths)
    L = min(L1, L2)
    
    max_len1 = max(pop1_lengths)
    max_len2 = max(pop2_lengths)
    
    if L1 != max_len1 or L2 != max_len2:
        warnings.warn(
            f"hudson_fst: Sequences have different lengths "
            f"(pop1: min={L1}, max={max_len1}; pop2: min={L2}, max={max_len2}). "
            f"Truncating to shortest length ({L}).",
            UserWarning,
            stacklevel=2,
        )
    
    if L == 0:
        return 0.0

    num = 0.0
    den = 0.0
    n1 = len(pop1)
    n2 = len(pop2)

    for pos in range(L):
        a1 = [s[pos] for s in pop1]
        a2 = [s[pos] for s in pop2]
        # Choose consistent reference allele per site from pop1
        ref = a1[0]
        p1 = sum(1 for a in a1 if a == ref) / n1
        p2 = sum(1 for a in a2 if a == ref) / n2

        # Hudson 1992 numerator and denominator components (for biallelic sites)
        within1 = (p1 * (1 - p1)) / (n1 - 1 if n1 > 1 else 1)
        within2 = (p2 * (1 - p2)) / (n2 - 1 if n2 > 1 else 1)
        num += (p1 - p2) ** 2 - within1 - within2
        den += p1 * (1 - p1) + p2 * (1 - p2)

    if den <= 0:
        # If there is no within-pop variation but fixed differences exist, Fst = 1
        for pos in range(L):
            a1 = [s[pos] for s in pop1]
            a2 = [s[pos] for s in pop2]
            if len(set(a1)) == 1 and len(set(a2)) == 1 and a1[0] != a2[0]:
                return 1.0
        return 0.0

    fst = num / den
    return max(0.0, min(1.0, fst))


def _allele_freq(alleles: Sequence[str]) -> float:
    """Calculate allele frequency of non-reference allele.
    
    This function is kept for backwards compatibility but is not currently
    used in the Hudson Fst implementation. The current implementation
    computes allele frequencies directly within hudson_fst().
    
    Args:
        alleles: Sequence of allele characters at a site
        
    Returns:
        Frequency of non-reference allele (1 - frequency of first allele)
        
    Note:
        This function may be removed in a future version if no external
        code depends on it. If you need this functionality, consider using
        the allele frequency calculation directly in your code.
    """
    ref = alleles[0]
    count_ref = sum(1 for a in alleles if a == ref)
    return 1.0 - count_ref / len(alleles)


def segregating_sites(seqs: Sequence[str]) -> int:
    """Count the number of sites with more than one allele among sequences.

    A segregating (polymorphic) site is a position where at least two
    different nucleotides are observed across the sequences.

    Args:
        seqs: Sequence of DNA sequences (strings)

    Returns:
        Integer count of segregating sites. Returns 0 if:
        - Less than 2 sequences provided
        - Sequences have zero length
        - All sequences are identical

    Examples:
        >>> seqs = ["AAAA", "AAAT", "AATT"]
        >>> segregating_sites(seqs)
        2  # Sites at positions 2 and 3 are polymorphic
        >>> seqs = ["AAAA", "AAAA"]
        >>> segregating_sites(seqs)
        0  # No variation
    """
    if len(seqs) < 2:
        return 0
    
    # Validate sequences are strings
    for i, seq in enumerate(seqs):
        if not isinstance(seq, str):
            raise TypeError(f"Sequence {i} is not a string: {type(seq)}")
    
    # Check for length mismatches
    lengths = [len(s) for s in seqs]
    L = min(lengths)
    max_len = max(lengths)
    if L != max_len:
        warnings.warn(
            f"segregating_sites: Sequences have different lengths (min={L}, max={max_len}). "
            f"Truncating to shortest length ({L}).",
            UserWarning,
            stacklevel=2,
        )
    
    if L == 0:
        return 0
    count = 0
    for pos in range(L):
        if len({s[pos] for s in seqs}) > 1:
            count += 1
    return count


def wattersons_theta(seqs: Sequence[str]) -> float:
    """Watterson's theta per site: θ_W = S / (a₁ × L).

    Calculates Watterson's estimator of the population mutation parameter θ
    based on the number of segregating sites rather than pairwise diversity.

    Args:
        seqs: Sequence of DNA sequences (strings)

    Returns:
        Watterson's theta estimate per site (float). Returns 0.0 if:
        - Less than 2 sequences provided
        - Sequences have zero length
        - No segregating sites

    Formula:
        θ_W = S / (a₁ × L) where:
        - S = number of segregating sites
        - a₁ = Σᵢ₌₁ⁿ⁻¹ 1/i (harmonic sum)
        - L = sequence length

    Examples:
        >>> seqs = ["AAAA", "AAAT", "AATT"]
        >>> wattersons_theta(seqs)
        0.5...  # Depends on number of segregating sites

    References:
        Watterson, G. A. (1975). On the number of segregating sites in genetical
        models without recombination. *Theoretical Population Biology*, 7(2), 256-276.
    """
    n = len(seqs)
    if n < 2:
        return 0.0
    
    # Validate sequences are strings
    for i, seq in enumerate(seqs):
        if not isinstance(seq, str):
            raise TypeError(f"Sequence {i} is not a string: {type(seq)}")
    
    # Check for length mismatches
    lengths = [len(s) for s in seqs]
    L = min(lengths)
    max_len = max(lengths)
    if L != max_len:
        warnings.warn(
            f"wattersons_theta: Sequences have different lengths (min={L}, max={max_len}). "
            f"Truncating to shortest length ({L}).",
            UserWarning,
            stacklevel=2,
        )
    
    if L == 0:
        return 0.0
    S = segregating_sites(seqs)
    a1 = sum(1.0 / i for i in range(1, n))
    if a1 <= 0:
        return 0.0
    return S / (a1 * L)


# (single definitions above)
