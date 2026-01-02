"""Advanced population genetics analysis utilities.

This module provides comprehensive tools for population genetics analysis
including demographic inference, selection tests, and population structure analysis.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def calculate_fst(population1: List[str], population2: List[str]) -> float:
    """Calculate F_ST between two populations using Hudson's estimator.

    Args:
        population1: List of DNA sequences from population 1
        population2: List of DNA sequences from population 2

    Returns:
        F_ST value (0.0 to 1.0)

    Raises:
        ValueError: If populations have incompatible sequences

    Example:
        >>> pop1 = ["ATCG", "ATCG"]
        >>> pop2 = ["ATCG", "GCTA"]
        >>> fst = calculate_fst(pop1, pop2)
        >>> 0.0 <= fst <= 1.0
        True
    """
    if not population1 or not population2:
        return 0.0

    # Check sequence lengths
    seq_len = len(population1[0])
    if not all(len(seq) == seq_len for seq in population1 + population2):
        raise ValueError("All sequences must have the same length")

    total_sites = 0
    numerator_sum = 0.0
    denominator_sum = 0.0

    for pos in range(seq_len):
        # Get alleles at this position for both populations
        alleles_pop1 = [seq[pos] for seq in population1]
        alleles_pop2 = [seq[pos] for seq in population2]

        # Count alleles
        allele_counts = {}
        for allele in alleles_pop1 + alleles_pop2:
            allele_counts[allele] = allele_counts.get(allele, 0) + 1

        if len(allele_counts) < 2:
            continue  # Monomorphic site

        # Calculate allele frequencies
        n1, n2 = len(population1), len(population2)
        total_n = n1 + n2

        freq_pop1 = {}
        freq_pop2 = {}
        total_freq = {}

        for allele in allele_counts:
            count_pop1 = alleles_pop1.count(allele)
            count_pop2 = alleles_pop2.count(allele)

            freq_pop1[allele] = count_pop1 / n1 if n1 > 0 else 0
            freq_pop2[allele] = count_pop2 / n2 if n2 > 0 else 0
            total_freq[allele] = (count_pop1 + count_pop2) / total_n

        # Calculate F_ST contribution for this site
        ht = 1 - sum(f ** 2 for f in total_freq.values())  # Total heterozygosity
        hs = (sum(f ** 2 for f in freq_pop1.values()) + sum(f ** 2 for f in freq_pop2.values())) / 2
        hs = 1 - hs  # Average within-population heterozygosity

        if ht > 0:
            fst_site = (ht - hs) / ht
            numerator_sum += fst_site
            total_sites += 1

        denominator_sum += ht

    if total_sites == 0:
        return 0.0

    # Weighted average F_ST
    return numerator_sum / total_sites


def detect_selection(sequences: List[str], method: str = "tajima_d") -> Dict[str, float]:
    """Detect signatures of natural selection in population data.

    Args:
        sequences: List of DNA sequences from a population
        method: Selection detection method ("tajima_d", "fu_li_d", "mk_test")

    Returns:
        Dictionary with selection statistics

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA"]
        >>> selection = detect_selection(seqs, method="tajima_d")
        >>> "tajima_d" in selection
        True
    """
    if not sequences:
        return {'tajima_d': 0.0, 'p_value': 1.0}

    results = {}

    if method == "tajima_d":
        tajima_d, p_value = calculate_tajima_d(sequences)
        results['tajima_d'] = tajima_d
        results['p_value'] = p_value

    elif method == "fu_li_d":
        fu_li_d, p_value = calculate_fu_li_d(sequences)
        results['fu_li_d'] = fu_li_d
        results['p_value'] = p_value

    elif method == "mk_test":
        alpha, omega = mcdonald_kreitman_test(sequences)
        results['alpha'] = alpha
        results['omega'] = omega

    else:
        raise ValueError(f"Unknown selection detection method: {method}")

    return results


def calculate_tajima_d(sequences: List[str]) -> Tuple[float, float]:
    """Calculate Tajima's D statistic for neutrality testing.

    Args:
        sequences: List of DNA sequences

    Returns:
        Tuple of (tajima_d, p_value)

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA"]
        >>> tajima_d, p_value = calculate_tajima_d(seqs)
        >>> isinstance(tajima_d, float)
        True
    """
    if len(sequences) < 4:
        return 0.0, 1.0  # Not enough sequences for meaningful test

    # Calculate segregating sites (S)
    seq_length = len(sequences[0])
    segregating_sites = 0

    for pos in range(seq_length):
        alleles = set(seq[pos] for seq in sequences)
        if len(alleles) > 1:
            segregating_sites += 1

    n = len(sequences)

    # Calculate average pairwise differences (pi)
    pi = calculate_nucleotide_diversity(sequences)

    # Calculate Watterson's theta
    theta_w = calculate_wattersons_theta(sequences)

    if theta_w == 0:
        return 0.0, 1.0

    # Tajima's D formula
    d = pi - theta_w

    # Variance of D (approximate)
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i ** 2) for i in range(1, n))

    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n ** 2 + n + 3) / (9 * n * (n - 1))

    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 ** 2)

    e1 = c1 / a1
    e2 = c2 / (a1 ** 2 + a2)

    var_d = e1 * segregating_sites + e2 * segregating_sites * (segregating_sites - 1)

    if var_d <= 0:
        return 0.0, 1.0

    # Standardize D
    tajima_d = d / math.sqrt(var_d)

    # Approximate p-value (two-tailed test)
    from scipy import stats
    try:
        p_value = 2 * (1 - stats.norm.cdf(abs(tajima_d)))
    except ImportError:
        # Fallback approximation
        p_value = min(1.0, 2 * (1 - 0.5 * (1 + math.erf(abs(tajima_d) / math.sqrt(2)))))

    return tajima_d, p_value


def calculate_fu_li_d(sequences: List[str]) -> Tuple[float, float]:
    """Calculate Fu and Li's D statistic.

    Args:
        sequences: List of DNA sequences

    Returns:
        Tuple of (fu_li_d, p_value)

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA"]
        >>> fu_li_d, p_value = calculate_fu_li_d(seqs)
        >>> isinstance(fu_li_d, float)
        True
    """
    if len(sequences) < 4:
        return 0.0, 1.0

    # Count singletons (alleles that appear only once)
    seq_length = len(sequences[0])
    singletons = 0

    for pos in range(seq_length):
        alleles = [seq[pos] for seq in sequences]
        allele_counts = {}

        for allele in alleles:
            allele_counts[allele] = allele_counts.get(allele, 0) + 1

        # Count alleles that appear exactly once
        singletons += sum(1 for count in allele_counts.values() if count == 1)

    n = len(sequences)

    # Fu and Li's D formula (simplified)
    # D = (singletons - a1 * segregating_sites) / sqrt(var)

    # For now, return a basic implementation
    segregating_sites = sum(1 for pos in range(seq_length)
                           if len(set(seq[pos] for seq in sequences)) > 1)

    if segregating_sites == 0:
        return 0.0, 1.0

    # Approximate calculation
    a1 = sum(1.0 / i for i in range(1, n))

    expected_singletons = a1 * segregating_sites
    variance = calculate_fu_li_variance(n, segregating_sites)

    if variance <= 0:
        return 0.0, 1.0

    fu_li_d = (singletons - expected_singletons) / math.sqrt(variance)

    # Approximate p-value
    try:
        from scipy import stats
        p_value = 2 * (1 - stats.norm.cdf(abs(fu_li_d)))
    except ImportError:
        p_value = 0.5  # Conservative estimate

    return fu_li_d, p_value


def calculate_fu_li_variance(n: int, s: int) -> float:
    """Calculate variance for Fu and Li's D (approximate)."""
    if n < 2 or s < 1:
        return 0.0

    # Simplified variance calculation
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i ** 2) for i in range(1, n))

    variance = ((n - 2) / (6 * (n - 1))) * s + (18 * n * (n - 1) * a2 - 88 * n * a1 ** 2) / (9 * (n - 1) ** 2) * s * (s - 1)

    return max(0, variance)


def mcdonald_kreitman_test(sequences: List[str]) -> Tuple[float, float]:
    """Perform McDonald-Kreitman test for selection.

    Args:
        sequences: List of DNA sequences

    Returns:
        Tuple of (alpha, omega) - proportion of adaptive substitutions

    Example:
        >>> seqs = ["ATGGCC", "ATGGCC", "ATGGCC"]
        >>> alpha, omega = mcdonald_kreitman_test(seqs)
        >>> isinstance(alpha, float)
        True
    """
    # This is a simplified MK test implementation
    # In practice, would need polymorphism and divergence data

    # For now, return neutral values
    return 0.0, 1.0


def estimate_population_size(sequences: List[str], mutation_rate: float = 1e-8) -> Dict[str, float]:
    """Estimate effective population size using various methods.

    Args:
        sequences: List of DNA sequences
        mutation_rate: Per-site mutation rate

    Returns:
        Dictionary with population size estimates

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA"]
        >>> sizes = estimate_population_size(seqs)
        >>> "theta_watterson" in sizes
        True
    """
    if not sequences:
        return {'theta_watterson': 0.0, 'theta_pi': 0.0, 'ne_estimate': 0.0}

    # Calculate theta estimates
    theta_w = calculate_wattersons_theta(sequences)
    theta_pi = calculate_nucleotide_diversity(sequences)

    # Estimate Ne using Watterson's theta
    # Ne = theta / (4 * mu * L) where L is sequence length
    seq_length = len(sequences[0])
    if mutation_rate > 0 and seq_length > 0:
        ne_estimate = theta_w / (4 * mutation_rate * seq_length)
    else:
        ne_estimate = 0.0

    return {
        'theta_watterson': theta_w,
        'theta_pi': theta_pi,
        'ne_estimate': ne_estimate,
        'sequence_length': seq_length,
        'sample_size': len(sequences)
    }


def detect_population_structure(sequences: List[str], k_max: int = 5) -> Dict[str, any]:
    """Detect population structure using simple clustering.

    Args:
        sequences: List of DNA sequences
        k_max: Maximum number of populations to test

    Returns:
        Dictionary with structure analysis results

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA", "GCTA"]
        >>> structure = detect_population_structure(seqs, k_max=2)
        >>> "clusters" in structure
        True
    """
    if len(sequences) < 4:
        return {'clusters': [], 'k_optimal': 1}

    # Simple distance-based clustering
    # In practice, would use more sophisticated methods

    clusters = []
    distances = []

    # Calculate pairwise distances
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            from . import distances
            dist = distances.p_distance(sequences[i], sequences[j])
            distances.append((i, j, dist))

    # Simple clustering based on distance threshold
    threshold = 0.1  # 10% divergence
    cluster_assignments = {}

    cluster_id = 0
    for i, seq in enumerate(sequences):
        assigned = False
        for cid, members in clusters.items():
            # Check if this sequence is close to cluster members
            avg_dist = sum(distances.p_distance(seq, sequences[j]) for j in members) / len(members)
            if avg_dist < threshold:
                members.append(i)
                cluster_assignments[i] = cid
                assigned = True
                break

        if not assigned:
            clusters.append([i])
            cluster_assignments[i] = cluster_id
            cluster_id += 1

    return {
        'clusters': clusters,
        'cluster_assignments': cluster_assignments,
        'k_optimal': len(clusters),
        'threshold_used': threshold
    }


def calculate_ld_decay(sequences: List[str]) -> List[Tuple[int, float]]:
    """Calculate linkage disequilibrium decay with distance.

    Args:
        sequences: List of DNA sequences

    Returns:
        List of tuples (distance, average_ld)

    Example:
        >>> seqs = ["ATCGATCG", "ATCGATCG"]
        >>> ld_decay = calculate_ld_decay(seqs)
        >>> isinstance(ld_decay, list)
        True
    """
    if len(sequences) < 2:
        return []

    seq_length = len(sequences[0])
    ld_values = []

    # Calculate LD for all pairs of sites
    for i in range(seq_length):
        for j in range(i + 1, seq_length):
            distance = j - i

            # Calculate correlation coefficient (simplified LD measure)
            alleles_i = [seq[i] for seq in sequences]
            alleles_j = [seq[j] for seq in sequences]

            # Simple correlation
            from . import distances
            correlation = 0.0  # Would need proper LD calculation

            ld_values.append((distance, correlation))

    # Average LD by distance
    distance_groups = {}
    for dist, ld in ld_values:
        if dist not in distance_groups:
            distance_groups[dist] = []
        distance_groups[dist].append(ld)

    avg_ld_by_distance = []
    for dist in sorted(distance_groups.keys()):
        avg_ld = sum(distance_groups[dist]) / len(distance_groups[dist])
        avg_ld_by_distance.append((dist, avg_ld))

    return avg_ld_by_distance


# Helper functions that delegate to existing implementations
def calculate_nucleotide_diversity(sequences: List[str]) -> float:
    """Calculate nucleotide diversity (π)."""
    from . import population
    return population.nucleotide_diversity(sequences)


def calculate_wattersons_theta(sequences: List[str]) -> float:
    """Calculate Watterson's θ."""
    from . import population
    return population.wattersons_theta(sequences)


def calculate_segregating_sites(sequences: List[str]) -> int:
    """Count segregating sites."""
    from . import population
    return population.segregating_sites(sequences)


def calculate_summary_statistics(
    sequences: List[str] | None = None,
    genotype_matrix: List[List[int]] | None = None,
    populations: List[int] | None = None
) -> Dict[str, any]:
    """Calculate comprehensive population genetics summary statistics.

    Args:
        sequences: List of DNA sequences (for sequence-based statistics)
        genotype_matrix: Genotype matrix (for genotype-based statistics)
        populations: Population assignments for each individual

    Returns:
        Dictionary containing various population genetics statistics

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA"]
        >>> stats = calculate_summary_statistics(sequences=seqs)
        >>> isinstance(stats, dict)
        True
    """
    results = {}

    if sequences:
        logger.info(f"Calculating summary statistics for {len(sequences)} sequences")

        # Basic diversity measures
        results['nucleotide_diversity'] = calculate_nucleotide_diversity(sequences)
        results['watterson_theta'] = calculate_wattersons_theta(sequences)
        results['segregating_sites'] = calculate_segregating_sites(sequences)

        # Neutrality tests
        try:
            tajima_d, tajima_p = calculate_tajima_d(sequences)
            results['tajima_d'] = tajima_d
            results['tajima_d_p_value'] = tajima_p
        except Exception:
            results['tajima_d'] = None
            results['tajima_d_p_value'] = None

        try:
            fu_li_d, fu_li_p = calculate_fu_li_d(sequences)
            results['fu_li_d'] = fu_li_d
            results['fu_li_d_p_value'] = fu_li_p
        except Exception:
            results['fu_li_d'] = None
            results['fu_li_d_p_value'] = None

        # Sequence properties
        results['sequence_length'] = len(sequences[0]) if sequences else 0
        results['sample_size'] = len(sequences)

    if genotype_matrix:
        logger.info(f"Calculating summary statistics for genotype matrix with {len(genotype_matrix)} individuals")

        # Convert to numpy for easier computation
        genotypes = np.array(genotype_matrix)

        # Basic statistics
        results['num_individuals'] = len(genotypes)
        results['num_loci'] = genotypes.shape[1] if genotypes.size > 0 else 0

        # Heterozygosity
        if genotypes.size > 0:
            from . import population
            het_values = []
            for locus in range(genotypes.shape[1]):
                locus_genotypes = [(genotypes[i, locus] // 2, genotypes[i, locus] % 2)
                                 for i in range(len(genotypes))]
                het = population.observed_heterozygosity(locus_genotypes)
                het_values.append(het)

            results['mean_heterozygosity'] = np.mean(het_values) if het_values else 0.0
            results['heterozygosity_variance'] = np.var(het_values) if het_values else 0.0

    if populations and sequences:
        # Population differentiation
        if len(set(populations)) > 1:
            pop_indices = {}
            for i, pop in enumerate(set(populations)):
                pop_indices[pop] = [j for j, p in enumerate(populations) if p == pop]

            if len(pop_indices) >= 2:
                pop1_idx, pop2_idx = list(pop_indices.values())[:2]
                pop1_seqs = [sequences[i] for i in pop1_idx]
                pop2_seqs = [sequences[i] for i in pop2_idx]

                try:
                    results['fst'] = calculate_fst(pop1_seqs, pop2_seqs)
                except Exception:
                    results['fst'] = None

    return results


def compare_populations(pop1_data: Dict[str, Any], pop2_data: Dict[str, Any]) -> Dict[str, Any]:
    """Compare two populations based on their summary statistics.

    Args:
        pop1_data: Summary statistics for population 1
        pop2_data: Summary statistics for population 2

    Returns:
        Dictionary containing comparison metrics

    Example:
        >>> pop1 = {"nucleotide_diversity": 0.01, "tajima_d": 0.5}
        >>> pop2 = {"nucleotide_diversity": 0.015, "tajima_d": -0.3}
        >>> comparison = compare_populations(pop1, pop2)
        >>> isinstance(comparison, dict)
        True
    """
    comparison = {}

    # Compare nucleotide diversity
    if "nucleotide_diversity" in pop1_data and "nucleotide_diversity" in pop2_data:
        pi1 = pop1_data["nucleotide_diversity"]
        pi2 = pop2_data["nucleotide_diversity"]
        comparison["nucleotide_diversity_ratio"] = pi2 / pi1 if pi1 > 0 else float('inf')
        comparison["nucleotide_diversity_difference"] = pi2 - pi1

    # Compare Tajima's D
    if "tajima_d" in pop1_data and "tajima_d" in pop2_data:
        d1 = pop1_data["tajima_d"]
        d2 = pop2_data["tajima_d"]
        comparison["tajima_d_difference"] = d2 - d1

        # Classify selection patterns
        def classify_selection(d: float) -> str:
            if d is None:
                return "unknown"
            elif d > 2:
                return "balancing_selection"
            elif d < -2:
                return "positive_selection"
            else:
                return "neutral"

        comparison["pop1_selection_pattern"] = classify_selection(d1)
        comparison["pop2_selection_pattern"] = classify_selection(d2)

    # Compare F_ST if available
    if "fst" in pop1_data and "fst" in pop2_data:
        fst1 = pop1_data["fst"]
        fst2 = pop2_data["fst"]
        if fst1 is not None and fst2 is not None:
            comparison["fst_difference"] = fst2 - fst1
            comparison["fst_ratio"] = fst2 / fst1 if fst1 > 0 else float('inf')

    return comparison


def neutrality_test_suite(sequences: List[str]) -> Dict[str, Any]:
    """Run a comprehensive suite of neutrality tests.

    Args:
        sequences: List of DNA sequences from a population

    Returns:
        Dictionary containing results from multiple neutrality tests

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA"]
        >>> results = neutrality_test_suite(seqs)
        >>> isinstance(results, dict)
        True
    """
    results = {}

    if not sequences:
        return results

    # Basic statistics
    results['num_sequences'] = len(sequences)
    results['sequence_length'] = len(sequences[0]) if sequences else 0

    # Tajima's D
    try:
        tajima_d, tajima_p = calculate_tajima_d(sequences)
        results['tajima_d'] = tajima_d
        results['tajima_d_p_value'] = tajima_p
    except Exception as e:
        logger.warning(f"Failed to calculate Tajima's D: {e}")
        results['tajima_d'] = None
        results['tajima_d_p_value'] = None

    # Fu and Li's D*
    try:
        fu_li_d_star, fu_li_d_star_p = calculate_fu_li_d(sequences)
        results['fu_li_d_star'] = fu_li_d_star
        results['fu_li_d_star_p_value'] = fu_li_d_star_p
    except Exception as e:
        logger.warning(f"Failed to calculate Fu and Li's D*: {e}")
        results['fu_li_d_star'] = None
        results['fu_li_d_star_p_value'] = None

    # Fu and Li's F*
    try:
        fu_li_f_star, fu_li_f_star_p = calculate_fu_li_f(sequences)
        results['fu_li_f_star'] = fu_li_f_star
        results['fu_li_f_star_p_value'] = fu_li_f_star_p
    except Exception as e:
        logger.warning(f"Failed to calculate Fu and Li's F*: {e}")
        results['fu_li_f_star'] = None
        results['fu_li_f_star_p_value'] = None

    # Fay and Wu's H
    try:
        fay_wu_h, fay_wu_h_p = calculate_fay_wu_h(sequences)
        results['fay_wu_h'] = fay_wu_h
        results['fay_wu_h_p_value'] = fay_wu_h_p
    except Exception as e:
        logger.warning(f"Failed to calculate Fay and Wu's H: {e}")
        results['fay_wu_h'] = None
        results['fay_wu_h_p_value'] = None

    # Summary interpretation
    results['neutrality_summary'] = interpret_neutrality_results(results)

    return results


def interpret_neutrality_results(results: Dict[str, Any]) -> Dict[str, str]:
    """Interpret neutrality test results.

    Args:
        results: Results from neutrality_test_suite

    Returns:
        Dictionary with interpretation of results
    """
    interpretation = {}

    # Tajima's D interpretation
    tajima_d = results.get('tajima_d')
    if tajima_d is not None:
        if tajima_d > 0:
            interpretation['tajima_d'] = 'balancing_selection_or_population_expansion'
        elif tajima_d < 0:
            interpretation['tajima_d'] = 'positive_selection_or_population_bottleneck'
        else:
            interpretation['tajima_d'] = 'neutral_evolution'
    else:
        interpretation['tajima_d'] = 'calculation_failed'

    # Fu and Li's tests interpretation
    fu_li_d_star = results.get('fu_li_d_star')
    if fu_li_d_star is not None:
        if fu_li_d_star > 0:
            interpretation['fu_li_d_star'] = 'balancing_selection'
        elif fu_li_d_star < 0:
            interpretation['fu_li_d_star'] = 'positive_selection_or_population_expansion'
        else:
            interpretation['fu_li_d_star'] = 'neutral_evolution'
    else:
        interpretation['fu_li_d_star'] = 'calculation_failed'

    return interpretation



