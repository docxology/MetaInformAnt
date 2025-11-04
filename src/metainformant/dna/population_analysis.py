"""Population genetics analysis orchestrators and workflows."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from .population import (
    allele_frequencies,
    hudson_fst,
    nucleotide_diversity,
    observed_heterozygosity,
    segregating_sites,
    tajimas_d,
    wattersons_theta,
)


def calculate_summary_statistics(
    sequences: Sequence[str] | None = None,
    genotype_matrix: Sequence[Sequence[int]] | None = None,
) -> dict[str, Any]:
    """Calculate comprehensive population genetics summary statistics.
    
    This function computes multiple population genetics statistics in a single
    call, providing a convenient way to get a complete overview of genetic
    diversity and structure.
    
    Args:
        sequences: Optional sequence of DNA sequences (strings). If provided,
            calculates sequence-based statistics (π, S, θ_W, Tajima's D).
        genotype_matrix: Optional genotype matrix (individuals × sites) with
            values 0/1. If provided, calculates allele frequency and heterozygosity
            statistics.
    
    Returns:
        Dictionary containing:
        - nucleotide_diversity: Average pairwise nucleotide differences per site (π)
        - segregating_sites: Number of polymorphic sites
        - wattersons_theta: Watterson's theta per site
        - tajimas_d: Tajima's D statistic (simplified)
        - allele_frequencies: List of allele frequencies per site (if genotype_matrix provided)
        - observed_heterozygosity: Proportion of heterozygous individuals (if genotype_matrix provided)
        - sample_size: Number of sequences/individuals
        - sequence_length: Length of sequences (if sequences provided)
    
    Examples:
        >>> seqs = ["AAAA", "AAAT", "AATT", "ATTT"]
        >>> stats = calculate_summary_statistics(sequences=seqs)
        >>> stats["nucleotide_diversity"]
        0.5...
        >>> stats["segregating_sites"]
        3
        
        >>> genotypes = [[0, 1, 0], [0, 1, 1], [1, 0, 1]]
        >>> stats = calculate_summary_statistics(genotype_matrix=genotypes)
        >>> stats["allele_frequencies"]
        [0.333..., 0.666..., 0.666...]
    """
    result: dict[str, Any] = {}
    
    if sequences:
        result["nucleotide_diversity"] = nucleotide_diversity(sequences)
        result["segregating_sites"] = segregating_sites(sequences)
        result["wattersons_theta"] = wattersons_theta(sequences)
        result["tajimas_d"] = tajimas_d(sequences)
        result["sample_size"] = len(sequences)
        if sequences:
            result["sequence_length"] = min(len(s) for s in sequences)
    
    if genotype_matrix:
        result["allele_frequencies"] = allele_frequencies(genotype_matrix)
        # Convert genotype matrix to diploid format for heterozygosity
        # Assuming each row is diploid, convert 0/1/2 encoding to (a1, a2) tuples
        diploid_genotypes = []
        for row in genotype_matrix:
            for site_val in row:
                if site_val == 0:
                    diploid_genotypes.append((0, 0))
                elif site_val == 1:
                    diploid_genotypes.append((0, 1))
                elif site_val == 2:
                    diploid_genotypes.append((1, 1))
                else:
                    diploid_genotypes.append((0, 0))  # Default for invalid
        
        if diploid_genotypes:
            result["observed_heterozygosity"] = observed_heterozygosity(diploid_genotypes)
        result["sample_size"] = len(genotype_matrix)
        if genotype_matrix:
            result["num_sites"] = len(genotype_matrix[0])
    
    return result


def compare_populations(
    pop1_sequences: Sequence[str] | None = None,
    pop2_sequences: Sequence[str] | None = None,
    pop1_genotypes: Sequence[Sequence[int]] | None = None,
    pop2_genotypes: Sequence[Sequence[int]] | None = None,
) -> dict[str, Any]:
    """Compare two populations using multiple statistics.
    
    Calculates within-population diversity and between-population differentiation
    for two populations. Can use either sequence data or genotype matrices.
    
    Args:
        pop1_sequences: Sequences from population 1 (strings)
        pop2_sequences: Sequences from population 2 (strings)
        pop1_genotypes: Genotype matrix for population 1 (individuals × sites)
        pop2_genotypes: Genotype matrix for population 2 (individuals × sites)
    
    Returns:
        Dictionary containing:
        - pop1_stats: Summary statistics for population 1
        - pop2_stats: Summary statistics for population 2
        - fst: Hudson's Fst between populations
        - differentiation: Qualitative assessment ("none", "low", "moderate", "high")
    
    Examples:
        >>> pop1 = ["AAAA", "AAAA", "AAAT"]
        >>> pop2 = ["TTTT", "TTTT", "TTTA"]
        >>> comparison = compare_populations(pop1_sequences=pop1, pop2_sequences=pop2)
        >>> comparison["fst"]
        0.9...
        >>> comparison["differentiation"]
        'high'
    """
    result: dict[str, Any] = {}
    
    if pop1_sequences and pop2_sequences:
        result["pop1_stats"] = calculate_summary_statistics(sequences=pop1_sequences)
        result["pop2_stats"] = calculate_summary_statistics(sequences=pop2_sequences)
        result["fst"] = hudson_fst(pop1_sequences, pop2_sequences)
    elif pop1_genotypes and pop2_genotypes:
        result["pop1_stats"] = calculate_summary_statistics(genotype_matrix=pop1_genotypes)
        result["pop2_stats"] = calculate_summary_statistics(genotype_matrix=pop2_genotypes)
        # For genotype matrices, we'd need to convert to sequences or use a different Fst method
        # For now, we'll note that Fst requires sequence data
        result["fst"] = None
        result["note"] = "Fst calculation requires sequence data, not genotype matrices"
    else:
        raise ValueError(
            "Must provide either (pop1_sequences, pop2_sequences) or "
            "(pop1_genotypes, pop2_genotypes)"
        )
    
    # Qualitative assessment of differentiation
    fst_val = result.get("fst", 0.0)
    if fst_val is None:
        result["differentiation"] = "unknown"
    elif fst_val < 0.05:
        result["differentiation"] = "none"
    elif fst_val < 0.15:
        result["differentiation"] = "low"
    elif fst_val < 0.25:
        result["differentiation"] = "moderate"
    else:
        result["differentiation"] = "high"
    
    return result


def neutrality_test_suite(
    sequences: Sequence[str],
) -> dict[str, Any]:
    """Run a suite of neutrality tests on sequence data.
    
    Computes multiple neutrality test statistics to assess whether observed
    variation is consistent with neutral evolution. Combines Tajima's D,
    nucleotide diversity, and Watterson's theta.
    
    Args:
        sequences: Sequence of DNA sequences (strings)
    
    Returns:
        Dictionary containing:
        - tajimas_d: Tajima's D statistic
        - nucleotide_diversity: π (average pairwise differences)
        - wattersons_theta: θ_W (Watterson's estimator)
        - segregating_sites: Number of polymorphic sites
        - pi_theta_ratio: Ratio of π to θ_W (should be ~1 under neutrality)
        - interpretation: Qualitative assessment of neutrality
    
    Examples:
        >>> seqs = ["AAAA", "AAAT", "AATT", "ATTT"]
        >>> results = neutrality_test_suite(seqs)
        >>> results["tajimas_d"]
        -0.5...
        >>> results["interpretation"]
        'negative_d' or 'neutral' or 'positive_d'
    """
    pi = nucleotide_diversity(sequences)
    theta_w = wattersons_theta(sequences)
    seg_sites = segregating_sites(sequences)
    tajima_d = tajimas_d(sequences)
    
    # Calculate ratio
    pi_theta_ratio = pi / theta_w if theta_w > 0 else 0.0
    
    # Interpretation
    if tajima_d < -2.0:
        interpretation = "strong_negative_d"
    elif tajima_d < -1.0:
        interpretation = "negative_d"
    elif tajima_d > 2.0:
        interpretation = "strong_positive_d"
    elif tajima_d > 1.0:
        interpretation = "positive_d"
    else:
        interpretation = "neutral"
    
    return {
        "tajimas_d": tajima_d,
        "nucleotide_diversity": pi,
        "wattersons_theta": theta_w,
        "segregating_sites": seg_sites,
        "pi_theta_ratio": pi_theta_ratio,
        "interpretation": interpretation,
        "sample_size": len(sequences),
    }

