"""Population genetics analysis utilities.

This module provides functions for calculating population genetics statistics,
including nucleotide diversity, neutrality tests, F-statistics, and other
metrics used in population genomics.
"""

from __future__ import annotations

import math
from typing import Iterable, List, Sequence, Tuple

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def allele_frequencies(genotype_matrix: Sequence[Sequence[int]]) -> List[float]:
    """Calculate allele frequencies from genotype matrix.

    Args:
        genotype_matrix: Matrix where each row is a locus and each column is an individual.
                        Values should be 0, 1, or 2 representing allele counts.

    Returns:
        List of allele frequencies (one per locus)
    """
    frequencies = []

    for locus in genotype_matrix:
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

    total_differences = 0
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

    F_ST measures the genetic differentiation between populations.

    Args:
        pop1: Sequences from population 1
        pop2: Sequences from population 2

    Returns:
        F_ST value (0.0 to 1.0)

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
    total_fst = 0.0
    valid_sites = 0

    for pos in range(seq_length):
        # Get alleles at this position for both populations
        alleles_pop1 = [seq[pos].upper() for seq in pop1 if seq[pos].upper() in "ATCG"]
        alleles_pop2 = [seq[pos].upper() for seq in pop2 if seq[pos].upper() in "ATCG"]

        if not alleles_pop1 or not alleles_pop2:
            continue

        # Calculate allele frequencies
        freq_pop1 = _allele_frequencies_from_list(alleles_pop1)
        freq_pop2 = _allele_frequencies_from_list(alleles_pop2)

        # Calculate F_ST for this locus
        fst_locus = _fst_single_locus(freq_pop1, freq_pop2)

        if not math.isnan(fst_locus):
            total_fst += fst_locus
            valid_sites += 1

    return total_fst / valid_sites if valid_sites > 0 else 0.0


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
                # Contribution: 2 * i^2 / (n * (n-1))
                theta_h += (2 * i * i) / (n * (n - 1))

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
        freqs = allele_frequencies([locus])  # Wrap in list for single locus
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
        Linkage disequilibrium coefficient D
    """
    if len(seqs) < 2:
        return 0.0

    if not _check_alignment(seqs):
        raise ValueError("Sequences must be aligned")

    if pos1 >= len(seqs[0]) or pos2 >= len(seqs[0]):
        raise ValueError("Position out of sequence bounds")

    # Get alleles at both positions
    alleles1 = [seq[pos1].upper() for seq in seqs if seq[pos1].upper() in "ATCG"]
    alleles2 = [seq[pos2].upper() for seq in seqs if seq[pos2].upper() in "ATCG"]

    if len(alleles1) != len(alleles2):
        raise ValueError("Positions have different numbers of valid alleles")

    # Calculate haplotype frequencies
    haplotypes = {}
    for a1, a2 in zip(alleles1, alleles2):
        hap = (a1, a2)
        haplotypes[hap] = haplotypes.get(hap, 0) + 1

    n = len(alleles1)

    # Calculate D = p_AB - p_A * p_B
    # Where A and B are derived alleles (assuming first allele is ancestral)
    alleles_a = set(alleles1)
    alleles_b = set(alleles2)

    if len(alleles_a) < 2 or len(alleles_b) < 2:
        return 0.0  # No variation

    # Use most common allele as ancestral
    ancestral_a = max(set(alleles1), key=alleles1.count)
    ancestral_b = max(set(alleles2), key=alleles2.count)

    # Find derived alleles (non-ancestral)
    derived_a = next((a for a in alleles_a if a != ancestral_a), None)
    derived_b = next((b for b in alleles_b if b != ancestral_b), None)

    if not derived_a or not derived_b:
        return 0.0

    # Calculate frequencies
    p_a = alleles1.count(ancestral_a) / n
    p_b = alleles2.count(ancestral_b) / n
    p_ab = haplotypes.get((ancestral_a, ancestral_b), 0) / n

    d = p_ab - p_a * p_b

    return d


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


def _allele_frequencies_from_list(alleles: List[str]) -> dict[str, float]:
    """Calculate allele frequencies from list of alleles."""
    if not alleles:
        return {}

    from collections import Counter

    counts = Counter(alleles)
    total = len(alleles)

    return {allele: count / total for allele, count in counts.items()}


def _fst_single_locus(freq1: dict[str, float], freq2: dict[str, float]) -> float:
    """Calculate F_ST for a single locus."""
    # Get all alleles
    all_alleles = set(freq1.keys()) | set(freq2.keys())

    # Calculate heterozygosity within populations
    h1 = 1 - sum(freq1.get(allele, 0) ** 2 for allele in all_alleles)
    h2 = 1 - sum(freq2.get(allele, 0) ** 2 for allele in all_alleles)

    # Average heterozygosity within populations
    h_s = (h1 + h2) / 2

    if h_s == 0:
        return 0.0

    # Calculate heterozygosity in total population
    total_freq = {}
    for allele in all_alleles:
        f1 = freq1.get(allele, 0)
        f2 = freq2.get(allele, 0)
        total_freq[allele] = (f1 + f2) / 2

    h_t = 1 - sum(freq**2 for freq in total_freq.values())

    # F_ST = (H_T - H_S) / H_T
    fst = (h_t - h_s) / h_t if h_t > 0 else 0.0

    return fst


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
