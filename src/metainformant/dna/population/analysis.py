"""Advanced population genetics analysis utilities.

This module provides comprehensive tools for population genetics analysis
including demographic inference, selection tests, and population structure analysis.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core.utils import logging
from metainformant.dna.annotation.gene_finding import GENETIC_CODES

logger = logging.get_logger(__name__)


def calculate_fst(population1: List[str], population2: List[str]) -> float:
    """Calculate F_ST between two populations from aligned sequences.

    Uses the heterozygosity-based per-site estimator: at each polymorphic
    site ``F_ST = (H_T - H_S) / H_T``, with ``H_T`` the heterozygosity of
    the pooled populations and ``H_S`` the mean within-population
    heterozygosity, averaged over polymorphic sites (Nei's G_ST style
    estimator). This is not Hudson's estimator; for Hudson's F_ST based on
    between- versus within-population pairwise diversity use
    :func:`metainformant.dna.population.core.hudson_fst`.

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

    for pos in range(seq_len):
        # Keep only unambiguous bases (ATCG, case-insensitive) so gaps and
        # ambiguous characters are not treated as alleles
        alleles_pop1 = [seq[pos].upper() for seq in population1 if seq[pos].upper() in "ATCG"]
        alleles_pop2 = [seq[pos].upper() for seq in population2 if seq[pos].upper() in "ATCG"]

        if not alleles_pop1 or not alleles_pop2:
            continue  # No usable data at this site

        # Count alleles
        allele_counts: Dict[str, int] = {}
        for allele in alleles_pop1 + alleles_pop2:
            allele_counts[allele] = allele_counts.get(allele, 0) + 1

        if len(allele_counts) < 2:
            continue  # Monomorphic site

        # Calculate allele frequencies
        n1, n2 = len(alleles_pop1), len(alleles_pop2)
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
        ht = 1 - sum(f**2 for f in total_freq.values())  # Total heterozygosity
        hs = (sum(f**2 for f in freq_pop1.values()) + sum(f**2 for f in freq_pop2.values())) / 2
        hs = 1 - hs  # Average within-population heterozygosity

        if ht > 0:
            fst_site = (ht - hs) / ht
            numerator_sum += fst_site
            total_sites += 1

    if total_sites == 0:
        return 0.0

    # Mean F_ST across polymorphic sites with usable data
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
        return {"tajima_d": 0.0, "p_value": 1.0}

    results = {}

    if method == "tajima_d":
        tajima_d, p_value = calculate_tajima_d(sequences)
        results["tajima_d"] = tajima_d
        results["p_value"] = p_value

    elif method == "fu_li_d":
        fu_li_d, p_value = calculate_fu_li_d(sequences)
        results["fu_li_d"] = fu_li_d
        results["p_value"] = p_value

    elif method == "mk_test":
        alpha, omega = mcdonald_kreitman_test(sequences)
        results["alpha"] = alpha
        results["omega"] = omega

    else:
        raise ValueError(f"Unknown selection detection method: {method}")

    return results


def _two_sided_normal_p(statistic: float) -> float:
    """Two-sided p-value of a standardized statistic under the normal approximation."""
    return math.erfc(abs(statistic) / math.sqrt(2.0))


def calculate_tajima_d(sequences: List[str]) -> Tuple[float, float]:
    """Calculate Tajima's D and its two-sided p-value.

    See :func:`metainformant.dna.population.core.tajimas_d`. The p-value uses
    the normal approximation, which is rough for small samples.

    Args:
        sequences: List of aligned DNA sequences

    Returns:
        Tuple of (tajima_d, p_value); ``(0.0, 1.0)`` for fewer than 4 sequences

    Example:
        >>> seqs = ["AAAA", "AAAA", "AACT", "AAAT"]
        >>> tajima_d, p_value = calculate_tajima_d(seqs)
        >>> round(tajima_d, 4)
        0.5916
    """
    if len(sequences) < 4:
        return 0.0, 1.0  # Not enough sequences for meaningful test

    from . import core as population

    tajima_d = population.tajimas_d(sequences)
    return tajima_d, _two_sided_normal_p(tajima_d)


def calculate_fu_li_d(sequences: List[str]) -> Tuple[float, float]:
    """Calculate Fu and Li's D* (no outgroup) and its two-sided p-value.

    See :func:`metainformant.dna.population.core.fu_and_li_d_star_from_sequences`.
    The p-value uses the normal approximation.

    Args:
        sequences: List of aligned DNA sequences

    Returns:
        Tuple of (fu_li_d_star, p_value); ``(0.0, 1.0)`` for fewer than 4 sequences

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA", "GCTT"]
        >>> fu_li_d, p_value = calculate_fu_li_d(seqs)
        >>> isinstance(fu_li_d, float)
        True
    """
    if len(sequences) < 4:
        return 0.0, 1.0

    from . import core as population

    d_star = population.fu_and_li_d_star_from_sequences(sequences)
    return d_star, _two_sided_normal_p(d_star)


def _fisher_exact_two_sided_2x2(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p-value for the 2x2 table ``[[a, b], [c, d]]``.

    Uses :func:`scipy.stats.fisher_exact` when SciPy is installed; otherwise
    computes the exact hypergeometric two-sided p-value using SciPy's
    definition: the sum of the probabilities of all tables with the same
    margins whose probability does not exceed the observed table's.
    """
    if min(a, b, c, d) < 0:
        raise ValueError("Counts must be non-negative")
    if a + b == 0 or c + d == 0 or a + c == 0 or b + d == 0:
        return 1.0  # A degenerate margin leaves nothing to test

    scipy_stats: Any  # None when SciPy is unavailable
    try:
        from scipy import stats as scipy_stats
    except ImportError:
        scipy_stats = None

    if scipy_stats is not None:
        return float(scipy_stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")[1])

    n = a + b + c + d
    row1 = a + b
    col1 = a + c

    def _prob(k: int) -> float:
        return (math.comb(row1, k) * math.comb(n - row1, col1 - k)) / math.comb(n, col1)

    lo = max(0, col1 - (n - row1))
    hi = min(row1, col1)
    observed = _prob(a)
    tolerance = observed * 1e-9
    return min(
        1.0,
        sum(_prob(k) for k in range(lo, hi + 1) if _prob(k) <= observed + tolerance),
    )


def mcdonald_kreitman_contingency(sequences: List[str]) -> Dict[str, Any]:
    """Build the McDonald-Kreitman 2x2 contingency table from coding sequences.

    The final sequence is the outgroup (divergence reference); all preceding
    sequences form the ingroup population sample (polymorphism). Sequences
    must be aligned coding sequences of equal length; bases beyond the last
    complete codon are ignored, and codon columns containing ambiguous bases
    are skipped. A site polymorphic within the ingroup is counted only as
    polymorphism (standard site-based MK counting), never additionally as
    divergence.

    Returns:
        Dict with the synonymous/nonsynonymous polymorphism (``Ps``, ``Pn``)
        and divergence (``Ds``, ``Dn``) counts, the neutrality index
        ``NI = (Pn / Ps) / (Dn / Ds)`` (``None`` when the ratio is
        undefined), the two-sided Fisher exact p-value for the table
        ``[[Pn, Ps], [Dn, Ds]]``, and the ``alpha`` / ``omega`` summary
        returned by :func:`mcdonald_kreitman_test`.
    """
    if len(sequences) < 3:
        # Fewer than two ingroup sequences plus one outgroup carries no
        # joint polymorphism/divergence information: report an empty table.
        return {
            "Ps": 0,
            "Pn": 0,
            "Ds": 0,
            "Dn": 0,
            "neutral_ratio": None,
            "alpha": 0.0,
            "omega": 0.0,
            "fisher_p": 1.0,
        }

    seq_len = len(sequences[0])
    if any(len(seq) != seq_len for seq in sequences):
        raise ValueError("Sequences must be aligned (equal length)")

    ingroup = [seq.upper() for seq in sequences[:-1]]
    outgroup = sequences[-1].upper()
    code = GENETIC_CODES.get(1, GENETIC_CODES[1])
    bases = frozenset("ATCG")

    ps = pn = ds = dn = 0
    for i in range(0, seq_len - 2, 3):
        in_codons = [seq[i : i + 3] for seq in ingroup]
        out_codon = outgroup[i : i + 3]
        if any(set(codon) - bases for codon in in_codons) or set(out_codon) - bases:
            continue  # Ambiguous or gapped codon column is unusable

        distinct = set(in_codons)
        if len(distinct) > 1:
            # Polymorphic within the ingroup
            amino_acids = {code[codon] for codon in distinct}
            if len(amino_acids) > 1:
                pn += 1
            else:
                ps += 1
        elif in_codons[0] != out_codon:
            # Monomorphic in the ingroup but divergent from the outgroup
            if code[in_codons[0]] != code[out_codon]:
                dn += 1
            else:
                ds += 1

    fisher_p = _fisher_exact_two_sided_2x2(pn, ps, dn, ds)

    if ps > 0 and dn > 0:
        neutral_ratio = (pn / ps) * (ds / dn)
        alpha = 1.0 - neutral_ratio
    else:
        neutral_ratio = None
        alpha = 0.0

    omega = dn / ds if ds > 0 else 0.0

    return {
        "Ps": ps,
        "Pn": pn,
        "Ds": ds,
        "Dn": dn,
        "neutral_ratio": neutral_ratio,
        "alpha": alpha,
        "omega": omega,
        "fisher_p": fisher_p,
    }


def mcdonald_kreitman_test(sequences: List[str]) -> Tuple[float, float]:
    """Perform the McDonald-Kreitman test for selection on coding sequences.

    The last sequence is treated as the outgroup and the preceding aligned
    coding sequences as the ingroup population sample (see
    :func:`mcdonald_kreitman_contingency` for the exact counting rules).

    Returns:
        Tuple of ``(alpha, omega)``:

        - ``alpha = 1 - NI`` estimates the proportion of adaptive
          substitutions, where ``NI = (Pn/Ps) / (Dn/Ds)`` is the neutrality
          index; 0.0 when the ratio is undefined (no adaptive signal
          estimable).
        - ``omega = Dn / Ds`` is the nonsynonymous-to-synonymous divergence
          ratio; 0.0 when Ds == 0.

    Example:
        >>> seqs = ["ATGGCC", "ATGGCC", "ATGGCC"]
        >>> alpha, omega = mcdonald_kreitman_test(seqs)
        >>> isinstance(alpha, float) and isinstance(omega, float)
        True
    """
    results = mcdonald_kreitman_contingency(sequences)
    return results["alpha"], results["omega"]


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
        return {"theta_watterson": 0.0, "theta_pi": 0.0, "ne_estimate": 0.0}

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
        "theta_watterson": theta_w,
        "theta_pi": theta_pi,
        "ne_estimate": ne_estimate,
        "sequence_length": seq_length,
        "sample_size": len(sequences),
    }


def detect_population_structure(sequences: List[str], k_max: int = 5) -> Dict[str, Any]:
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
        return {"clusters": [], "k_optimal": 1}

    # Simple distance-based clustering
    # In practice, would use more sophisticated methods

    def _p_distance_simple(seq_a: str, seq_b: str) -> float:
        """Calculate p-distance between two sequences."""
        if len(seq_a) != len(seq_b):
            min_len = min(len(seq_a), len(seq_b))
            seq_a = seq_a[:min_len]
            seq_b = seq_b[:min_len]
        if len(seq_a) == 0:
            return 0.0
        mismatches = sum(1 for a, b in zip(seq_a, seq_b) if a != b)
        return mismatches / len(seq_a)

    # Calculate pairwise distances
    pairwise_dists: List[Tuple[int, int, float]] = []
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            dist = _p_distance_simple(sequences[i], sequences[j])
            pairwise_dists.append((i, j, dist))

    # Simple clustering based on distance threshold
    threshold = 0.1  # 10% divergence
    cluster_map: Dict[int, List[int]] = {}
    cluster_assignments: Dict[int, int] = {}

    cluster_id = 0
    for i in range(len(sequences)):
        assigned = False
        for cid, members in cluster_map.items():
            # Check if this sequence is close to cluster members
            avg_dist = sum(_p_distance_simple(sequences[i], sequences[j]) for j in members) / len(members)
            if avg_dist < threshold:
                members.append(i)
                cluster_assignments[i] = cid
                assigned = True
                break

        if not assigned:
            cluster_map[cluster_id] = [i]
            cluster_assignments[i] = cluster_id
            cluster_id += 1

    clusters = list(cluster_map.values())

    return {
        "clusters": clusters,
        "cluster_assignments": cluster_assignments,
        "k_optimal": len(clusters),
        "threshold_used": threshold,
    }


def calculate_ld_decay(sequences: List[str], max_distance: int = 0) -> List[Tuple[int, float]]:
    """Calculate linkage disequilibrium (r²) decay with physical distance.

    Computes pairwise r² between biallelic polymorphic sites and averages by distance.

    Args:
        sequences: List of aligned DNA sequences
        max_distance: Maximum distance to consider (0 = all)

    Returns:
        List of tuples (distance, average_r_squared)

    Example:
        >>> seqs = ["ATCGATCG", "GCTAGCTA", "ATCGATCG", "GCTAGCTA"]
        >>> ld_decay = calculate_ld_decay(seqs)
        >>> isinstance(ld_decay, list)
        True
    """
    from collections import Counter

    if len(sequences) < 2:
        return []

    seq_length = len(sequences[0])

    # Find biallelic polymorphic sites
    polymorphic = []
    for pos in range(seq_length):
        alleles = set(seq[pos].upper() for seq in sequences if seq[pos].upper() in "ACGT")
        if len(alleles) == 2:
            polymorphic.append(pos)

    if len(polymorphic) < 2:
        return []

    ld_by_distance: Dict[int, List[float]] = {}

    for idx_i in range(len(polymorphic)):
        for idx_j in range(idx_i + 1, len(polymorphic)):
            pos_i = polymorphic[idx_i]
            pos_j = polymorphic[idx_j]
            distance = pos_j - pos_i

            if max_distance > 0 and distance > max_distance:
                continue

            alleles_i = [seq[pos_i].upper() for seq in sequences if seq[pos_i].upper() in "ACGT"]
            alleles_j = [seq[pos_j].upper() for seq in sequences if seq[pos_j].upper() in "ACGT"]

            if len(alleles_i) != len(alleles_j) or len(alleles_i) < 2:
                continue

            # Encode as 0/1 (major=0, minor=1)
            count_i = Counter(alleles_i)
            count_j = Counter(alleles_j)
            major_i = count_i.most_common(1)[0][0]
            major_j = count_j.most_common(1)[0][0]

            x = [0 if a == major_i else 1 for a in alleles_i]
            y = [0 if a == major_j else 1 for a in alleles_j]

            n_hap = len(x)
            p_a = sum(x) / n_hap
            p_b = sum(y) / n_hap

            if p_a == 0 or p_a == 1 or p_b == 0 or p_b == 1:
                continue

            p_ab = sum(1 for xi, yi in zip(x, y) if xi == 1 and yi == 1) / n_hap
            d = p_ab - p_a * p_b
            r_squared = (d * d) / (p_a * (1 - p_a) * p_b * (1 - p_b))

            if distance not in ld_by_distance:
                ld_by_distance[distance] = []
            ld_by_distance[distance].append(r_squared)

    avg_ld_by_distance = []
    for dist in sorted(ld_by_distance.keys()):
        avg_r2 = sum(ld_by_distance[dist]) / len(ld_by_distance[dist])
        avg_ld_by_distance.append((dist, avg_r2))

    return avg_ld_by_distance


# Helper functions that delegate to existing implementations
def calculate_nucleotide_diversity(sequences: List[str]) -> float:
    """Calculate nucleotide diversity (π)."""
    from . import core as population

    return population.nucleotide_diversity(sequences)


def calculate_wattersons_theta(sequences: List[str]) -> float:
    """Calculate Watterson's θ."""
    from . import core as population

    return population.wattersons_theta(sequences)


def calculate_segregating_sites(sequences: List[str]) -> int:
    """Count segregating sites."""
    from . import core as population

    return population.segregating_sites(sequences)


def calculate_summary_statistics(
    sequences: List[str] | None = None,
    genotype_matrix: List[List[int]] | None = None,
    populations: List[int] | None = None,
) -> Dict[str, Any]:
    """Calculate comprehensive population genetics summary statistics.

    Args:
        sequences: List of DNA sequences (for sequence-based statistics)
        genotype_matrix: Individuals x loci matrix of diploid 0/1/2 alternate-allele
            counts (for genotype-based statistics); other values count as missing
        populations: Population assignments for each individual

    Returns:
        Dictionary containing various population genetics statistics; ``fu_li_d``
        holds Fu and Li's D*, which needs no outgroup

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA"]
        >>> stats = calculate_summary_statistics(sequences=seqs)
        >>> isinstance(stats, dict)
        True
    """
    results: Dict[str, Any] = {}

    if sequences:
        logger.info(f"Calculating summary statistics for {len(sequences)} sequences")

        # Basic diversity measures
        results["nucleotide_diversity"] = calculate_nucleotide_diversity(sequences)
        results["watterson_theta"] = calculate_wattersons_theta(sequences)
        results["segregating_sites"] = calculate_segregating_sites(sequences)

        # Neutrality tests
        try:
            tajima_d, tajima_p = calculate_tajima_d(sequences)
            results["tajima_d"] = tajima_d
            results["tajima_d_p_value"] = tajima_p
        except Exception:
            results["tajima_d"] = None
            results["tajima_d_p_value"] = None

        try:
            fu_li_d, fu_li_p = calculate_fu_li_d(sequences)
            results["fu_li_d"] = fu_li_d
            results["fu_li_d_p_value"] = fu_li_p
        except Exception:
            results["fu_li_d"] = None
            results["fu_li_d_p_value"] = None

        # Sequence properties
        results["sequence_length"] = len(sequences[0]) if sequences else 0
        results["sample_size"] = len(sequences)

    if genotype_matrix:
        logger.info(f"Calculating summary statistics for genotype matrix with {len(genotype_matrix)} individuals")

        # Convert to numpy for easier computation
        genotypes = np.array(genotype_matrix)

        # Basic statistics
        results["num_individuals"] = len(genotypes)
        results["num_loci"] = genotypes.shape[1] if genotypes.size > 0 else 0

        # Observed heterozygosity: entries are diploid 0/1/2 alternate-allele
        # counts, so only 1 is heterozygous; other values count as missing.
        if genotypes.size > 0:
            called = np.isin(genotypes, (0, 1, 2)).sum(axis=0)
            heterozygous = (genotypes == 1).sum(axis=0)
            het_values = np.divide(heterozygous, called, out=np.zeros(genotypes.shape[1]), where=called > 0)

            results["mean_heterozygosity"] = float(np.mean(het_values))
            results["heterozygosity_variance"] = float(np.var(het_values))

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
                    results["fst"] = calculate_fst(pop1_seqs, pop2_seqs)
                except Exception:
                    results["fst"] = None

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
        comparison["nucleotide_diversity_ratio"] = pi2 / pi1 if pi1 > 0 else float("inf")
        comparison["nucleotide_diversity_difference"] = pi2 - pi1

    # Compare Tajima's D
    if "tajima_d" in pop1_data and "tajima_d" in pop2_data:
        d1 = pop1_data["tajima_d"]
        d2 = pop2_data["tajima_d"]
        comparison["tajima_d_difference"] = d2 - d1

        # Classify selection patterns
        def classify_selection(d: float | None) -> str:
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
            comparison["fst_ratio"] = fst2 / fst1 if fst1 > 0 else float("inf")

    return comparison


def calculate_fay_wu_h(sequences: List[str], outgroup: Optional[str] = None) -> Tuple[float, Optional[float]]:
    """Calculate Fay and Wu's H and, with an outgroup, its two-sided p-value.

    ``h`` is :func:`metainformant.dna.population.core.fay_wu_h_from_sequences`
    (per-site theta_pi - theta_H). Its null distribution needs derived-allele
    frequencies, so the p-value comes from Zeng et al.'s (2006) normalized H
    (:func:`metainformant.dna.population.core.normalized_fay_wu_h`) under the
    normal approximation and is None when no outgroup polarizes the sites of
    a polymorphic sample. A monomorphic sample has H = 0 and p-value 1.0.

    Args:
        sequences: List of aligned DNA sequences
        outgroup: Optional outgroup sequence aligned to ``sequences``

    Returns:
        Tuple of (fay_wu_h, p_value); ``(0.0, 1.0)`` for fewer than 4 sequences

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA", "GCTA"]
        >>> h, p = calculate_fay_wu_h(seqs)
        >>> isinstance(h, float), p is None
        (True, True)
    """
    if len(sequences) < 4:
        return 0.0, 1.0

    from . import core as population

    h = population.fay_wu_h_from_sequences(sequences, outgroup)
    if outgroup is not None:
        return h, _two_sided_normal_p(population.normalized_fay_wu_h(sequences, outgroup))
    if population.segregating_sites(sequences) == 0:
        return h, 1.0
    return h, None


def calculate_fu_li_f(sequences: List[str]) -> Tuple[float, float]:
    """Calculate Fu and Li's F* (no outgroup) and its two-sided p-value.

    See :func:`metainformant.dna.population.core.fu_and_li_f_star_from_sequences`.
    The p-value uses the normal approximation.

    Args:
        sequences: List of aligned DNA sequences

    Returns:
        Tuple of (fu_li_f_star, p_value); ``(0.0, 1.0)`` for fewer than 4 sequences

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA", "GCTA"]
        >>> f_star, p = calculate_fu_li_f(seqs)
        >>> isinstance(f_star, float)
        True
    """
    if len(sequences) < 4:
        return 0.0, 1.0

    from . import core as population

    f_star = population.fu_and_li_f_star_from_sequences(sequences)
    return f_star, _two_sided_normal_p(f_star)


def neutrality_test_suite(sequences: List[str], outgroup: Optional[str] = None) -> Dict[str, Any]:
    """Run a comprehensive suite of neutrality tests.

    Args:
        sequences: List of DNA sequences from a population
        outgroup: Optional outgroup sequence; without it ``fay_wu_h_p_value`` is None

    Returns:
        Dictionary containing results from multiple neutrality tests

    Example:
        >>> seqs = ["ATCG", "ATCG", "GCTA"]
        >>> results = neutrality_test_suite(seqs)
        >>> isinstance(results, dict)
        True
    """
    results: Dict[str, Any] = {}

    if not sequences:
        return results

    # Basic statistics
    results["num_sequences"] = len(sequences)
    results["sequence_length"] = len(sequences[0]) if sequences else 0

    # Tajima's D
    try:
        tajima_d, tajima_p = calculate_tajima_d(sequences)
        results["tajima_d"] = tajima_d
        results["tajima_d_p_value"] = tajima_p
    except Exception as e:
        logger.warning(f"Failed to calculate Tajima's D: {e}")
        results["tajima_d"] = None
        results["tajima_d_p_value"] = None

    # Fu and Li's D*
    try:
        fu_li_d_star, fu_li_d_star_p = calculate_fu_li_d(sequences)
        results["fu_li_d_star"] = fu_li_d_star
        results["fu_li_d_star_p_value"] = fu_li_d_star_p
    except Exception as e:
        logger.warning(f"Failed to calculate Fu and Li's D*: {e}")
        results["fu_li_d_star"] = None
        results["fu_li_d_star_p_value"] = None

    # Fu and Li's F*
    try:
        fu_li_f_star, fu_li_f_star_p = calculate_fu_li_f(sequences)
        results["fu_li_f_star"] = fu_li_f_star
        results["fu_li_f_star_p_value"] = fu_li_f_star_p
    except Exception as e:
        logger.warning(f"Failed to calculate Fu and Li's F*: {e}")
        results["fu_li_f_star"] = None
        results["fu_li_f_star_p_value"] = None

    # Fay and Wu's H
    try:
        fay_wu_h, fay_wu_h_p = calculate_fay_wu_h(sequences, outgroup)
        results["fay_wu_h"] = fay_wu_h
        results["fay_wu_h_p_value"] = fay_wu_h_p
    except Exception as e:
        logger.warning(f"Failed to calculate Fay and Wu's H: {e}")
        results["fay_wu_h"] = None
        results["fay_wu_h_p_value"] = None

    # Summary interpretation
    results["neutrality_summary"] = interpret_neutrality_results(results)

    return results


_NEUTRALITY_LABELS = {
    "tajima_d": ("balancing_selection_or_population_bottleneck", "positive_selection_or_population_expansion"),
    "fu_li_d_star": ("balancing_selection", "positive_selection_or_population_expansion"),
}


def interpret_neutrality_results(results: Dict[str, Any], alpha: float = 0.05) -> Dict[str, str]:
    """Interpret neutrality test results.

    A statistic gets a directional label only when it departs from neutrality:
    its ``<name>_p_value`` is below ``alpha``, or the results carry no p-value
    and the statistic is nonzero. Otherwise it is ``"neutral_evolution"``,
    meaning consistent with neutrality.

    Args:
        results: Results from neutrality_test_suite
        alpha: Significance level for the p-values

    Returns:
        Dictionary with interpretation of results
    """
    interpretation = {}
    for name, (positive, negative) in _NEUTRALITY_LABELS.items():
        value = results.get(name)
        p_value = results.get(f"{name}_p_value")
        if value is None:
            interpretation[name] = "calculation_failed"
        elif (p_value is not None and p_value >= alpha) or value == 0:
            interpretation[name] = "neutral_evolution"
        else:
            interpretation[name] = positive if value > 0 else negative
    return interpretation
