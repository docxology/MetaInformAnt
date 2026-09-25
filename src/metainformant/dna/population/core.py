"""Population genetics analysis utilities.

This module provides functions for calculating population genetics statistics,
including nucleotide diversity, neutrality tests, F-statistics, and other
metrics used in population genomics.
"""

from __future__ import annotations

import math
from typing import Dict, Iterable, List, Sequence, Tuple, cast

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def allele_frequencies(
    genotype_matrix: Sequence[Sequence[int]] | Sequence[str],
) -> List[float] | List[Dict[str, float]]:
    """Calculate allele frequencies from genotype counts or aligned sequences.

    Args:
        genotype_matrix: Either a genotype matrix where each row is a locus
            and each column is an individual (0/1/2 allele counts), or a
            sequence of aligned DNA strings.

    Returns:
        Genotype input returns one alternate-allele frequency per locus.
        Sequence input returns one base-frequency mapping per alignment site.
    """
    if len(genotype_matrix) == 0:
        return []

    first_item = genotype_matrix[0]
    if isinstance(first_item, str):
        seqs = [str(seq).upper() for seq in genotype_matrix]
        if not _check_alignment(seqs):
            raise ValueError("Sequences must be aligned (same length)")

        site_frequencies: List[Dict[str, float]] = []
        for site_idx in range(len(seqs[0])):
            counts: Dict[str, int] = {}
            for seq in seqs:
                base = seq[site_idx]
                if base in "ATCG":
                    counts[base] = counts.get(base, 0) + 1
            total = sum(counts.values())
            site_frequencies.append(
                {base: count / total for base, count in counts.items()} if total else {}
            )
        return site_frequencies

    genotype_counts = cast("Sequence[Sequence[int]]", genotype_matrix)

    frequencies: List[float] = []

    for locus in genotype_counts:
        total_alleles = sum(locus)
        total_possible = len(locus) * 2  # diploid

        if total_possible > 0:
            freq = total_alleles / total_possible
        else:
            freq = 0.0

        frequencies.append(freq)

    return frequencies


def observed_heterozygosity(genotypes: Iterable[Tuple[int, int]]) -> float:
    """Calculate observed heterozygosity from genotype data.

    Args:
        genotypes: Iterable of (allele1, allele2) tuples

    Returns:
        Observed heterozygosity (0.0 to 1.0)
    """
    hetero_count = 0
    total_count = 0

    for allele1, allele2 in genotypes:
        if allele1 != allele2:  # heterozygous
            hetero_count += 1
        total_count += 1

    return hetero_count / total_count if total_count > 0 else 0.0


def nucleotide_diversity(seqs: Sequence[str]) -> float:
    """Calculate nucleotide diversity (π) from sequence alignment.

    π = average number of nucleotide differences per site

    Args:
        seqs: Sequence of aligned DNA sequences

    Returns:
        Nucleotide diversity (π)
    """
    if len(seqs) < 2:
        return 0.0

    if not _check_alignment(seqs):
        raise ValueError("Sequences must be aligned (same length)")

    total_differences = 0.0
    total_comparisons = 0
    seq_length = len(seqs[0])

    # Compare all pairs of sequences
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            seq1 = seqs[i].upper()
            seq2 = seqs[j].upper()

            differences = 0
            valid_sites = 0

            for pos in range(seq_length):
                base1 = seq1[pos]
                base2 = seq2[pos]

                # Skip gaps and ambiguous bases
                if base1 in "ATCG" and base2 in "ATCG":
                    if base1 != base2:
                        differences += 1
                    valid_sites += 1

            if valid_sites > 0:
                total_differences += differences / valid_sites
                total_comparisons += 1

    return total_differences / total_comparisons if total_comparisons > 0 else 0.0


def tajimas_d(seqs: Sequence[str]) -> float:
    """Calculate Tajima's D statistic.

    Tajima's D compares nucleotide diversity (π) with the number of
    segregating sites to detect deviations from neutral evolution.

    Args:
        seqs: Sequence of aligned DNA sequences

    Returns:
        Tajima's D value

    Raises:
        ValueError: If insufficient sequences or data
    """
    if len(seqs) < 4:
        raise ValueError("Tajima's D requires at least 4 sequences")

    if not _check_alignment(seqs):
        raise ValueError("Sequences must be aligned")

    # Calculate π (nucleotide diversity)
    pi = nucleotide_diversity(seqs)

    # Calculate S (segregating sites)
    s = segregating_sites(seqs)

    # Calculate θ (Watterson's estimator)
    n = len(seqs)
    theta = wattersons_theta(seqs)

    if theta == 0:
        return 0.0

    # Tajima's D = (π - θ) / sqrt(Var(π - θ))
    # Simplified calculation (approximation)
    d = (pi - theta) / math.sqrt(_variance_pi_theta(n, s))

    return d


def wattersons_theta(seqs: Sequence[str]) -> float:
    """Calculate Watterson's θ (theta) estimator.

    θ = S / a_n where a_n is the sum of 1/i for i=1 to n-1

    Args:
        seqs: Sequence of aligned DNA sequences

    Returns:
        Watterson's θ
    """
    if len(seqs) < 2:
        return 0.0

    if not _check_alignment(seqs):
        raise ValueError("Sequences must be aligned")

    s = segregating_sites(seqs)
    n = len(seqs)

    # Calculate a_n = sum(1/i for i in 1 to n-1)
    a_n = sum(1.0 / i for i in range(1, n))

    return s / a_n if a_n > 0 else 0.0


def segregating_sites(seqs: Sequence[str]) -> int:
    """Count the number of segregating sites in aligned sequences.

    A segregating site is a position where at least two different nucleotides
    are observed (excluding gaps and ambiguous bases).

    Args:
        seqs: Sequence of aligned DNA sequences

    Returns:
        Number of segregating sites
    """
    if len(seqs) < 2:
        return 0

    if not _check_alignment(seqs):
        raise ValueError("Sequences must be aligned")

    seq_length = len(seqs[0])
    segregating_count = 0

    for pos in range(seq_length):
        bases_at_pos = set()

        for seq in seqs:
            base = seq[pos].upper()
            if base in "ATCG":
                bases_at_pos.add(base)

        if len(bases_at_pos) > 1:
            segregating_count += 1

    return segregating_count


def hudson_fst(pop1: Sequence[str], pop2: Sequence[str]) -> float:
    """Calculate Hudson's F_ST between two populations.

    Implements Hudson's (1992) moment estimator as popularised by Bhatia et
    al. (2013). At each site the reference allele is the most common pooled
    allele (ties broken alphabetically); with ``p1``/``p2`` its frequency in
    each population and ``n1``/``n2`` the number of valid gene copies
    sampled per population,

        num = (p1 - p2)^2 - p1*(1 - p1)/(n1 - 1) - p2*(1 - p2)/(n2 - 1)
        den = p1*(1 - p2) + p2*(1 - p1)

    and ``F_ST = sum(num) / sum(den)`` across usable sites. Sites with
    fewer than two valid copies in either population are skipped. When no
    site contributes a positive denominator, the result is 1.0 if the
    populations showed a frequency difference (fixed for different
    alleles) and 0.0 otherwise (identically fixed).

    Note: the estimate can be slightly negative when the populations are
    effectively undifferentiated -- a known property of the unbiased
    moment estimator. This is not the heterozygosity-based per-site
    estimator in
    :func:`metainformant.dna.population.analysis.calculate_fst`.

    Args:
        pop1: Sequences from population 1
        pop2: Sequences from population 2

    Returns:
        F_ST value

    Raises:
        ValueError: If populations have different sequence lengths
    """
    if not pop1 or not pop2:
        raise ValueError("Both populations must contain sequences")

    if not _check_alignment(pop1) or not _check_alignment(pop2):
        raise ValueError("Sequences within each population must be aligned")

    if len(pop1[0]) != len(pop2[0]):
        raise ValueError("Populations must have same sequence length")

    seq_length = len(pop1[0])
    numerator_sum = 0.0
    denominator_sum = 0.0
    saw_valid_pair = False
    saw_shared_allele = False

    for pos in range(seq_length):
        # Keep only unambiguous bases (ATCG, case-insensitive)
        alleles_pop1 = [seq[pos].upper() for seq in pop1 if seq[pos].upper() in "ATCG"]
        alleles_pop2 = [seq[pos].upper() for seq in pop2 if seq[pos].upper() in "ATCG"]

        # Track whether the populations are comparable and ever fixed for
        # different alleles, for the degenerate no-usable-site fallback.
        if alleles_pop1 and alleles_pop2:
            saw_valid_pair = True
            if set(alleles_pop1) & set(alleles_pop2):
                saw_shared_allele = True

        n1 = len(alleles_pop1)
        n2 = len(alleles_pop2)
        if n1 < 2 or n2 < 2:
            continue  # Cannot estimate the sampling correction

        reference = _most_common_allele(alleles_pop1 + alleles_pop2)
        p1 = alleles_pop1.count(reference) / n1
        p2 = alleles_pop2.count(reference) / n2

        numerator_sum += (
            (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
        )
        denominator_sum += p1 * (1 - p2) + p2 * (1 - p1)

    if denominator_sum > 0:
        return numerator_sum / denominator_sum
    if saw_valid_pair and not saw_shared_allele:
        # No site was usable for the corrected estimator, but every
        # comparable site is fixed for different alleles: complete
        # differentiation.
        return 1.0
    return 0.0


def fu_and_li_d_star_from_sequences(seqs: Sequence[str]) -> float:
    """Calculate Fu and Li's D* statistic from sequences.

    D* compares the number of singletons with nucleotide diversity.

    Args:
        seqs: Sequence of aligned DNA sequences

    Returns:
        Fu and Li's D* value
    """
    if len(seqs) < 4:
        raise ValueError("Fu and Li's D* requires at least 4 sequences")

    # Count singletons (mutations that appear only once)
    singletons = _count_singletons(seqs)

    # Calculate π
    pi = nucleotide_diversity(seqs)

    # Calculate D* = (n/(n-1)) * singletons - π
    n = len(seqs)
    if n <= 1:
        return 0.0

    d_star = (n / (n - 1)) * singletons - pi

    return d_star


def fu_and_li_f_star_from_sequences(seqs: Sequence[str]) -> float:
    """Calculate Fu and Li's F* statistic from sequences.

    F* compares the number of singletons with the number of segregating sites.

    Args:
        seqs: Sequence of aligned DNA sequences

    Returns:
        Fu and Li's F* value
    """
    if len(seqs) < 4:
        raise ValueError("Fu and Li's F* requires at least 4 sequences")

    # Count singletons and segregating sites
    singletons = _count_singletons(seqs)
    s = segregating_sites(seqs)

    if s == 0:
        return 0.0

    # Calculate F* = (n/(n-1)) * singletons - 1 + 1/(n-1)
    n = len(seqs)
    f_star = (n / (n - 1)) * singletons - 1 + 1 / (n - 1)

    return f_star


def fay_wu_h_from_sequences(seqs: Sequence[str], outgroup: str | None = None) -> float:
    """Calculate Fay and Wu's H statistic from sequences.

    H = π - θ_H, where θ_H weights each SNP by the square of its derived
    allele frequency. This statistic detects positive selection.

    θ_H is normalized per site (divided by the alignment length) so that H
    is length-independent and directly comparable to the per-site π.

    When no outgroup is provided, the most frequent allele at each site
    is assumed to be ancestral (parsimony assumption).

    Args:
        seqs: Sequence of aligned DNA sequences
        outgroup: Optional outgroup sequence for ancestral state inference

    Returns:
        Fay and Wu's H value (negative values suggest positive selection)
    """
    if len(seqs) < 4:
        raise ValueError("Fay and Wu's H requires at least 4 sequences")

    n = len(seqs)
    if not seqs[0]:
        raise ValueError("Sequences cannot be empty")

    seq_len = len(seqs[0])

    # Check all sequences have same length
    if not all(len(s) == seq_len for s in seqs):
        raise ValueError("All sequences must have the same length")

    # Calculate π (nucleotide diversity)
    pi = nucleotide_diversity(seqs)

    # Calculate θ_H (Fay and Wu's theta)
    # θ_H = Σ 2 * i^2 * S_i / (n * (n-1))
    # where S_i is the number of sites where derived allele is at frequency i/n

    theta_h = 0.0
    valid_sites = 0

    for pos in range(seq_len):
        # Get nucleotides at this position
        nucs = [s[pos].upper() for s in seqs if s[pos].upper() in "ACGT"]

        if len(nucs) < 2:
            continue

        # Count allele frequencies
        from collections import Counter

        allele_counts = Counter(nucs)

        if len(allele_counts) < 2:
            # Not a polymorphic site
            continue

        valid_sites += 1

        # Determine ancestral allele
        if outgroup and pos < len(outgroup):
            ancestral = outgroup[pos].upper()
            if ancestral not in "ACGT":
                # If outgroup has ambiguous base, use most frequent
                ancestral = allele_counts.most_common(1)[0][0]
        else:
            # Assume most frequent allele is ancestral (parsimony)
            ancestral = allele_counts.most_common(1)[0][0]

        # Calculate contribution to theta_H
        # Sum over derived alleles
        for allele, count in allele_counts.items():
            if allele != ancestral:
                # This is a derived allele
                i = count  # Frequency count of derived allele
                theta_h += (2 * i * i) / (n * (n - 1))

    # Normalize θ_H per site so H is length-independent and comparable to
    # the per-site π (raw accumulation above is a sum over sites).
    theta_h /= seq_len

    # Fay and Wu's H = π - θ_H
    h = pi - theta_h

    return h


def expected_heterozygosity(genotype_matrix: Sequence[Sequence[int]]) -> float:
    """Calculate expected heterozygosity (gene diversity) from genotype matrix.

    Args:
        genotype_matrix: Matrix where each row is a locus and each column is an individual

    Returns:
        Expected heterozygosity (0.0 to 1.0)
    """
    if not genotype_matrix:
        return 0.0

    total_he = 0.0

    for locus in genotype_matrix:
        freqs = cast(
            "List[float]", allele_frequencies([locus])
        )  # Wrap in list for single locus
        p = freqs[0]  # Allele frequency
        q = 1 - p  # Other allele frequency

        # H_E = 2pq for diploid organisms
        he = 2 * p * q
        total_he += he

    return total_he / len(genotype_matrix) if genotype_matrix else 0.0


def fixation_index(genotypes: Iterable[Tuple[int, int]], expected_freq: float) -> float:
    """Calculate fixation index (F) from genotype data.

    F = 1 - (observed heterozygosity / expected heterozygosity)

    Args:
        genotypes: Iterable of (allele1, allele2) tuples
        expected_freq: Expected heterozygosity under Hardy-Weinberg

    Returns:
        Fixation index F
    """
    if expected_freq <= 0:
        return 0.0

    observed_het = observed_heterozygosity(genotypes)
    f = 1 - (observed_het / expected_freq)

    return f


def hardy_weinberg_allele_freqs(p: float, q: float) -> Tuple[float, float, float]:
    """Calculate Hardy-Weinberg genotype frequencies.

    Args:
        p: Frequency of allele A
        q: Frequency of allele a (q = 1 - p)

    Returns:
        Tuple of (AA_freq, Aa_freq, aa_freq)
    """
    if not (0 <= p <= 1) or not (0 <= q <= 1) or abs(p + q - 1) > 1e-6:
        raise ValueError("Allele frequencies must sum to 1")

    homo_dom = p * p  # AA frequency
    het_freq = 2 * p * q  # Aa frequency
    homo_rec = q * q  # aa frequency

    return (homo_dom, het_freq, homo_rec)


def linkage_disequilibrium(seqs: Sequence[str], pos1: int, pos2: int) -> float:
    """Calculate linkage disequilibrium (D) between two positions.

    Args:
        seqs: Sequence of aligned DNA sequences
        pos1: First position
        pos2: Second position

    Returns:
        Linkage disequilibrium coefficient
        ``D = f_AB - p_A * p_B``, where the reference alleles A and B are
        the most common allele at each site (ties broken alphabetically) and
        f_AB counts haplotypes carrying both reference alleles. Positive D
        means the reference alleles are coupled (in phase); negative D means
        they are in repulsion (e.g. haplotypes Ab and aB only, which yields
        D = -0.25 at equal frequencies).
    """
    if len(seqs) < 2:
        return 0.0

    if not _check_alignment(seqs):
        raise ValueError("Sequences must be aligned")

    if pos1 >= len(seqs[0]) or pos2 >= len(seqs[0]):
        raise ValueError("Position out of sequence bounds")

    # Keep a sequence only if it has valid (non-ambiguous, non-gap) alleles
    # at BOTH positions, so haplotypes are always paired from the same
    # sequence.
    alleles1: List[str] = []
    alleles2: List[str] = []
    for seq in seqs:
        a1 = seq[pos1].upper()
        a2 = seq[pos2].upper()
        if a1 in "ATCG" and a2 in "ATCG":
            alleles1.append(a1)
            alleles2.append(a2)

    n = len(alleles1)

    if n < 2:
        return 0.0

    # Calculate haplotype frequencies
    haplotypes: Dict[Tuple[str, str], int] = {}
    for a1, a2 in zip(alleles1, alleles2):
        hap = (a1, a2)
        haplotypes[hap] = haplotypes.get(hap, 0) + 1

    # D is only defined for two segregating sites
    if len(set(alleles1)) < 2 or len(set(alleles2)) < 2:
        return 0.0  # No variation

    # Reference alleles ("A" and "B" below): the most common allele at each
    # site, with ties broken alphabetically. The choice must be a pure
    # function of the input -- iterating over a `set` (the previous
    # implementation) made the sign of D depend on Python's per-process
    # string hash order.
    ref_a = _most_common_allele(alleles1)
    ref_b = _most_common_allele(alleles2)

    # D = f_AB - p_A * p_B, where f_AB is the frequency of haplotypes
    # carrying the reference allele at BOTH sites (phase taken from the
    # same sequence). Coupled reference alleles give D > 0, repulsion
    # (e.g. haplotypes Ab / aB in equal frequency) gives D < 0.
    p_a = alleles1.count(ref_a) / n
    p_b = alleles2.count(ref_b) / n
    p_ab = haplotypes.get((ref_a, ref_b), 0) / n

    return p_ab - p_a * p_b


def _most_common_allele(alleles: List[str]) -> str:
    """Return the most common allele, breaking ties alphabetically.

    Deterministic on the input order: unlike ``max(set(alleles), key=...)``
    this never depends on set iteration order.
    """
    counts: Dict[str, int] = {}
    for allele in alleles:
        counts[allele] = counts.get(allele, 0) + 1
    return min(counts, key=lambda allele: (-counts[allele], allele))


def _check_alignment(seqs: Sequence[str]) -> bool:
    """Check if sequences are properly aligned (same length)."""
    if not seqs:
        return True

    length = len(seqs[0])
    return all(len(seq) == length for seq in seqs)


def _count_singletons(seqs: Sequence[str]) -> int:
    """Count singleton mutations (alleles that appear only once)."""
    if len(seqs) < 2:
        return 0

    seq_length = len(seqs[0])
    singletons = 0

    for pos in range(seq_length):
        alleles = [seq[pos].upper() for seq in seqs if seq[pos].upper() in "ATCG"]

        if len(alleles) < 2:
            continue

        # Count frequency of each allele
        from collections import Counter

        counts = Counter(alleles)

        # Check if any allele appears exactly once
        if 1 in counts.values():
            singletons += 1

    return singletons


def _variance_pi_theta(n: int, s: int) -> float:
    """Calculate variance of π - θ for Tajima's D."""
    # Simplified variance calculation
    # Full calculation involves complex formulas
    if n < 2 or s < 1:
        return 1.0

    # Approximation
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i * i) for i in range(1, n))

    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))

    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)

    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)

    variance = e1 * s + e2 * s * (s - 1)

    return max(variance, 0.01)  # Avoid division by zero
