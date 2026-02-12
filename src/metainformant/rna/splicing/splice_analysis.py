"""Splicing event classification, quantification, and differential analysis.

This module provides tools for classifying alternative splicing events,
computing Percent Spliced In (PSI) values, performing differential splicing
tests between conditions, and discovering novel splice junctions.

All implementations are pure Python with optional numpy/scipy acceleration.
"""

from __future__ import annotations

import math
import random
from collections import defaultdict
from typing import Any, Dict, List, Literal, Optional, Sequence, Tuple

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependency handling
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]

try:
    from scipy import stats as scipy_stats
    from scipy.special import beta as beta_func
    from scipy.special import betainc

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None  # type: ignore[assignment]


# =============================================================================
# Type Definitions
# =============================================================================

SplicingEventType = Literal[
    "exon_skipping",
    "intron_retention",
    "alt_5ss",
    "alt_3ss",
    "mutually_exclusive",
    "unknown",
]


# =============================================================================
# Splicing Event Classification
# =============================================================================


def classify_splicing_events(
    junctions: list[dict],
    gene_model: dict | None = None,
) -> list[dict]:
    """Classify splice junctions into alternative splicing event types.

    Uses junction coordinates and optionally a gene model to classify each
    junction as one of the standard alternative splicing event types.

    Without a gene model, classification uses coordinate-based heuristics
    comparing junctions that share donor or acceptor sites. With a gene
    model, exon annotations enable precise event typing.

    Args:
        junctions: List of junction dicts from detect_splice_junctions(),
            each with keys: chrom, start, end, strand, read_count.
        gene_model: Optional gene annotation dict with structure:
            {
                "genes": {
                    "gene_id": {
                        "chrom": str,
                        "strand": str,
                        "exons": [{"start": int, "end": int, "exon_id": str}, ...],
                        "introns": [{"start": int, "end": int}, ...]
                    }
                }
            }
            If None, classification uses coordinate-based inference only.

    Returns:
        List of event dicts, each containing:
            - event_type (str): One of "exon_skipping", "intron_retention",
              "alt_5ss", "alt_3ss", "mutually_exclusive", "unknown"
            - junctions (list[dict]): The junction(s) involved in this event
            - chrom (str): Chromosome
            - region_start (int): Start of the affected region
            - region_end (int): End of the affected region
            - confidence (float): Classification confidence score (0-1)

    Example:
        >>> junctions = [
        ...     {"chrom": "chr1", "start": 100, "end": 300, "strand": "+", "read_count": 10},
        ...     {"chrom": "chr1", "start": 100, "end": 200, "strand": "+", "read_count": 5},
        ...     {"chrom": "chr1", "start": 200, "end": 300, "strand": "+", "read_count": 5},
        ... ]
        >>> events = classify_splicing_events(junctions)
    """
    if not junctions:
        return []

    events: list[dict] = []

    # Group junctions by chromosome and strand
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for junc in junctions:
        key = (junc["chrom"], junc.get("strand", "+"))
        grouped[key].append(junc)

    for (chrom, strand), chrom_junctions in grouped.items():
        # Sort by start position
        sorted_juncs = sorted(chrom_junctions, key=lambda x: (x["start"], x["end"]))

        # If gene model available, use annotation-based classification
        if gene_model is not None:
            annotation_events = _classify_with_gene_model(sorted_juncs, chrom, strand, gene_model)
            events.extend(annotation_events)
            continue

        # Coordinate-based classification
        coord_events = _classify_by_coordinates(sorted_juncs, chrom, strand)
        events.extend(coord_events)

    logger.info(f"Classified {len(events)} splicing events from {len(junctions)} junctions")

    return events


def _classify_with_gene_model(
    junctions: list[dict],
    chrom: str,
    strand: str,
    gene_model: dict,
) -> list[dict]:
    """Classify junctions using gene model annotations.

    Args:
        junctions: Sorted junctions for one chromosome/strand.
        chrom: Chromosome name.
        strand: Strand.
        gene_model: Gene annotation dictionary.

    Returns:
        List of classified event dicts.
    """
    events: list[dict] = []
    genes = gene_model.get("genes", {})

    # Find relevant genes on this chromosome/strand
    relevant_genes = {
        gid: ginfo for gid, ginfo in genes.items() if ginfo.get("chrom") == chrom and ginfo.get("strand", "+") == strand
    }

    for junc in junctions:
        junc_start = junc["start"]
        junc_end = junc["end"]
        best_event_type: SplicingEventType = "unknown"
        best_confidence = 0.0

        for gid, ginfo in relevant_genes.items():
            exons = sorted(ginfo.get("exons", []), key=lambda e: e["start"])
            introns = ginfo.get("introns", [])

            # Check for intron retention: junction spans exactly one annotated intron
            for intron in introns:
                if abs(junc_start - intron["start"]) <= 5 and abs(junc_end - intron["end"]) <= 5:
                    best_event_type = "intron_retention"
                    best_confidence = 0.9
                    break

            # Check for exon skipping: junction skips one or more internal exons
            if best_event_type == "unknown":
                skipped_exons = [e for e in exons if e["start"] > junc_start and e["end"] < junc_end]
                if skipped_exons:
                    best_event_type = "exon_skipping"
                    best_confidence = 0.85

            # Check for alt 5' or 3' splice sites
            if best_event_type == "unknown":
                for other_junc in junctions:
                    if other_junc is junc:
                        continue
                    # Same acceptor, different donor -> alt 5'ss
                    if abs(other_junc["end"] - junc_end) <= 5 and abs(other_junc["start"] - junc_start) > 5:
                        best_event_type = "alt_5ss"
                        best_confidence = 0.7
                        break
                    # Same donor, different acceptor -> alt 3'ss
                    if abs(other_junc["start"] - junc_start) <= 5 and abs(other_junc["end"] - junc_end) > 5:
                        best_event_type = "alt_3ss"
                        best_confidence = 0.7
                        break

        events.append(
            {
                "event_type": best_event_type,
                "junctions": [junc],
                "chrom": chrom,
                "region_start": junc_start,
                "region_end": junc_end,
                "confidence": best_confidence,
            }
        )

    return events


def _classify_by_coordinates(
    junctions: list[dict],
    chrom: str,
    strand: str,
) -> list[dict]:
    """Classify junctions using coordinate-based heuristics.

    When no gene model is available, infers event types from the spatial
    relationships between junctions sharing donor or acceptor sites.

    Args:
        junctions: Sorted junctions for one chromosome/strand.
        chrom: Chromosome name.
        strand: Strand.

    Returns:
        List of classified event dicts.
    """
    events: list[dict] = []

    # Index junctions by donor (start) and acceptor (end) sites
    by_donor: dict[int, list[dict]] = defaultdict(list)
    by_acceptor: dict[int, list[dict]] = defaultdict(list)

    for junc in junctions:
        by_donor[junc["start"]].append(junc)
        by_acceptor[junc["end"]].append(junc)

    classified_juncs: set[int] = set()

    # Detect exon skipping: a long junction that spans the region of two
    # shorter junctions sharing an internal acceptor/donor
    for i, junc in enumerate(junctions):
        if id(junc) in classified_juncs:
            continue

        junc_start = junc["start"]
        junc_end = junc["end"]

        # Look for shorter junctions that together span the same region
        for j, other in enumerate(junctions):
            if i == j or id(other) in classified_juncs:
                continue

            # Check if other junction is contained within this one
            if other["start"] >= junc_start and other["end"] <= junc_end:
                # Find a complementary junction
                for k, third in enumerate(junctions):
                    if k in (i, j) or id(third) in classified_juncs:
                        continue
                    # Complementary: third starts near other's end, ends near junc's end
                    if abs(third["start"] - other["end"]) <= 10 and abs(third["end"] - junc_end) <= 10:
                        events.append(
                            {
                                "event_type": "exon_skipping",
                                "junctions": [junc, other, third],
                                "chrom": chrom,
                                "region_start": junc_start,
                                "region_end": junc_end,
                                "confidence": 0.7,
                            }
                        )
                        classified_juncs.update([id(junc), id(other), id(third)])
                        break

    # Detect alternative 5' splice sites: same acceptor, different donors
    for acceptor, juncs in by_acceptor.items():
        unclassified = [j for j in juncs if id(j) not in classified_juncs]
        if len(unclassified) >= 2:
            events.append(
                {
                    "event_type": "alt_5ss",
                    "junctions": unclassified,
                    "chrom": chrom,
                    "region_start": min(j["start"] for j in unclassified),
                    "region_end": acceptor,
                    "confidence": 0.6,
                }
            )
            for j in unclassified:
                classified_juncs.add(id(j))

    # Detect alternative 3' splice sites: same donor, different acceptors
    for donor, juncs in by_donor.items():
        unclassified = [j for j in juncs if id(j) not in classified_juncs]
        if len(unclassified) >= 2:
            events.append(
                {
                    "event_type": "alt_3ss",
                    "junctions": unclassified,
                    "chrom": chrom,
                    "region_start": donor,
                    "region_end": max(j["end"] for j in unclassified),
                    "confidence": 0.6,
                }
            )
            for j in unclassified:
                classified_juncs.add(id(j))

    # Remaining junctions are unclassified
    for junc in junctions:
        if id(junc) not in classified_juncs:
            events.append(
                {
                    "event_type": "unknown",
                    "junctions": [junc],
                    "chrom": chrom,
                    "region_start": junc["start"],
                    "region_end": junc["end"],
                    "confidence": 0.0,
                }
            )

    return events


# =============================================================================
# PSI Computation
# =============================================================================


def compute_psi(
    inclusion_reads: int,
    exclusion_reads: int,
    effective_length_ratio: float = 1.0,
) -> dict:
    """Compute Percent Spliced In (PSI) with confidence interval.

    PSI measures the fraction of transcripts that include a particular
    exon or splice junction. Uses a Beta distribution model to estimate
    the PSI value and its confidence interval.

    The formula is:
        PSI = (inclusion / effective_length_ratio) /
              (inclusion / effective_length_ratio + exclusion)

    Args:
        inclusion_reads: Number of reads supporting the inclusion isoform
            (e.g., reads spanning the included exon).
        exclusion_reads: Number of reads supporting the exclusion isoform
            (e.g., reads spanning the skipping junction).
        effective_length_ratio: Ratio of effective lengths between inclusion
            and exclusion isoforms. Accounts for length bias in read
            assignment. Default 1.0 assumes equal effective lengths.

    Returns:
        Dictionary with keys:
            - psi (float): Point estimate of PSI (0.0 to 1.0)
            - ci_low (float): Lower bound of 95% confidence interval
            - ci_high (float): Upper bound of 95% confidence interval
            - total_reads (int): Total supporting reads
            - inclusion_reads (int): Input inclusion count
            - exclusion_reads (int): Input exclusion count
            - method (str): "beta_distribution"

    Raises:
        ValueError: If read counts are negative or effective_length_ratio <= 0.

    Example:
        >>> result = compute_psi(inclusion_reads=30, exclusion_reads=10)
        >>> 0.7 < result["psi"] < 0.8
        True
        >>> result["ci_low"] < result["psi"] < result["ci_high"]
        True
    """
    if inclusion_reads < 0 or exclusion_reads < 0:
        raise ValueError(
            f"Read counts must be non-negative: inclusion={inclusion_reads}, " f"exclusion={exclusion_reads}"
        )
    if effective_length_ratio <= 0:
        raise ValueError(f"effective_length_ratio must be positive, got {effective_length_ratio}")

    total_reads = inclusion_reads + exclusion_reads

    if total_reads == 0:
        logger.warning("Zero total reads, PSI undefined")
        return {
            "psi": 0.0,
            "ci_low": 0.0,
            "ci_high": 1.0,
            "total_reads": 0,
            "inclusion_reads": 0,
            "exclusion_reads": 0,
            "method": "beta_distribution",
        }

    # Adjust inclusion reads by effective length ratio
    adjusted_inclusion = inclusion_reads / effective_length_ratio
    psi = adjusted_inclusion / (adjusted_inclusion + exclusion_reads)

    # Beta distribution parameters (Bayesian with uniform prior)
    # Using adjusted counts as pseudo-observations
    alpha = adjusted_inclusion + 1  # +1 for uniform prior
    beta_param = exclusion_reads + 1

    # Compute 95% confidence interval using Beta distribution
    ci_low, ci_high = _beta_confidence_interval(alpha, beta_param, confidence=0.95)

    return {
        "psi": round(psi, 6),
        "ci_low": round(ci_low, 6),
        "ci_high": round(ci_high, 6),
        "total_reads": total_reads,
        "inclusion_reads": inclusion_reads,
        "exclusion_reads": exclusion_reads,
        "method": "beta_distribution",
    }


def _beta_confidence_interval(
    alpha: float,
    beta_param: float,
    confidence: float = 0.95,
) -> tuple[float, float]:
    """Compute confidence interval from Beta distribution.

    Args:
        alpha: Alpha parameter of Beta distribution.
        beta_param: Beta parameter of Beta distribution.
        confidence: Confidence level (e.g., 0.95 for 95% CI).

    Returns:
        Tuple of (lower_bound, upper_bound).
    """
    tail = (1 - confidence) / 2

    if HAS_SCIPY:
        ci_low = float(scipy_stats.beta.ppf(tail, alpha, beta_param))
        ci_high = float(scipy_stats.beta.ppf(1 - tail, alpha, beta_param))
    else:
        # Pure Python approximation using normal approximation to Beta
        mean = alpha / (alpha + beta_param)
        variance = (alpha * beta_param) / ((alpha + beta_param) ** 2 * (alpha + beta_param + 1))
        std = math.sqrt(variance)

        # z-score for 95% CI
        z = 1.96 if confidence == 0.95 else _norm_ppf(1 - tail)
        ci_low = max(0.0, mean - z * std)
        ci_high = min(1.0, mean + z * std)

    return ci_low, ci_high


def _norm_ppf(p: float) -> float:
    """Approximate normal distribution percent point function (inverse CDF).

    Uses the Abramowitz and Stegun rational approximation.

    Args:
        p: Probability (0 < p < 1).

    Returns:
        z-score corresponding to the given probability.
    """
    if p <= 0 or p >= 1:
        return 0.0

    # Rational approximation constants
    a = [
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577518672690e02,
        -3.066479806614716e01,
        2.506628277459239e00,
    ]
    b = [
        -5.447609879822406e01,
        1.615858368580409e02,
        -1.556989798598866e02,
        6.680131188771972e01,
        -1.328068155288572e01,
    ]

    if p < 0.5:
        q = math.sqrt(-2 * math.log(p))
        num = ((((a[0] * q + a[1]) * q + a[2]) * q + a[3]) * q + a[4]) * q + a[5]
        den = ((((b[0] * q + b[1]) * q + b[2]) * q + b[3]) * q + b[4]) * q + 1
        return num / den
    else:
        q = math.sqrt(-2 * math.log(1 - p))
        num = ((((a[0] * q + a[1]) * q + a[2]) * q + a[3]) * q + a[4]) * q + a[5]
        den = ((((b[0] * q + b[1]) * q + b[2]) * q + b[3]) * q + b[4]) * q + 1
        return -(num / den)


# =============================================================================
# Differential Splicing
# =============================================================================


def differential_splicing(
    psi_group1: list[float],
    psi_group2: list[float],
    method: str = "empirical_bayes",
) -> dict:
    """Test for differential splicing between two conditions.

    Compares PSI distributions between two groups using delta-PSI
    (difference in mean PSI) and a statistical test to assess significance.

    Supports empirical Bayes shrinkage to improve estimates with small
    sample sizes, as well as standard non-parametric tests.

    Args:
        psi_group1: PSI values for condition 1 (e.g., control).
            Each value should be between 0 and 1.
        psi_group2: PSI values for condition 2 (e.g., treatment).
            Each value should be between 0 and 1.
        method: Statistical method for testing:
            - "empirical_bayes": Empirical Bayes shrinkage with moderated
              t-test (recommended for small samples)
            - "wilcoxon": Wilcoxon rank-sum (Mann-Whitney U) test
            - "permutation": Permutation test on delta-PSI

    Returns:
        Dictionary with keys:
            - delta_psi (float): Mean PSI difference (group2 - group1)
            - p_value (float): Statistical significance
            - mean_psi_group1 (float): Mean PSI in group 1
            - mean_psi_group2 (float): Mean PSI in group 2
            - sd_group1 (float): Standard deviation in group 1
            - sd_group2 (float): Standard deviation in group 2
            - method (str): Method used for testing
            - significant (bool): Whether delta_psi is significant
              (|delta_psi| >= 0.1 and p_value < 0.05)
            - n_group1 (int): Sample size group 1
            - n_group2 (int): Sample size group 2

    Raises:
        ValueError: If either group is empty or contains values outside [0, 1].

    Example:
        >>> control_psi = [0.8, 0.75, 0.82, 0.78]
        >>> treatment_psi = [0.4, 0.45, 0.38, 0.42]
        >>> result = differential_splicing(control_psi, treatment_psi)
        >>> result["delta_psi"] < -0.3
        True
    """
    if not psi_group1 or not psi_group2:
        raise ValueError("Both groups must have at least one PSI value")

    # Validate PSI ranges
    for i, val in enumerate(psi_group1):
        if not (0.0 <= val <= 1.0):
            raise ValueError(f"PSI values must be in [0,1], got {val} in group1[{i}]")
    for i, val in enumerate(psi_group2):
        if not (0.0 <= val <= 1.0):
            raise ValueError(f"PSI values must be in [0,1], got {val} in group2[{i}]")

    n1 = len(psi_group1)
    n2 = len(psi_group2)

    mean1 = sum(psi_group1) / n1
    mean2 = sum(psi_group2) / n2
    delta_psi = mean2 - mean1

    # Standard deviations
    if n1 > 1:
        var1 = sum((x - mean1) ** 2 for x in psi_group1) / (n1 - 1)
        sd1 = math.sqrt(var1)
    else:
        var1 = 0.0
        sd1 = 0.0

    if n2 > 1:
        var2 = sum((x - mean2) ** 2 for x in psi_group2) / (n2 - 1)
        sd2 = math.sqrt(var2)
    else:
        var2 = 0.0
        sd2 = 0.0

    # Statistical test
    if method == "empirical_bayes":
        p_value = _empirical_bayes_test(psi_group1, psi_group2, mean1, mean2, var1, var2)
    elif method == "wilcoxon":
        p_value = _wilcoxon_test(psi_group1, psi_group2)
    elif method == "permutation":
        p_value = _permutation_test(psi_group1, psi_group2, n_permutations=10000)
    else:
        raise ValueError(f"Unknown method: {method}. Valid: empirical_bayes, wilcoxon, permutation")

    # Significance threshold: |delta_PSI| >= 0.1 and p < 0.05
    significant = abs(delta_psi) >= 0.1 and p_value < 0.05

    return {
        "delta_psi": round(delta_psi, 6),
        "p_value": round(p_value, 8),
        "mean_psi_group1": round(mean1, 6),
        "mean_psi_group2": round(mean2, 6),
        "sd_group1": round(sd1, 6),
        "sd_group2": round(sd2, 6),
        "method": method,
        "significant": significant,
        "n_group1": n1,
        "n_group2": n2,
    }


def _empirical_bayes_test(
    group1: list[float],
    group2: list[float],
    mean1: float,
    mean2: float,
    var1: float,
    var2: float,
) -> float:
    """Moderated t-test with empirical Bayes variance shrinkage.

    Shrinks per-group variances toward a pooled prior, improving
    statistical power for small sample sizes (limma-style approach).

    Args:
        group1: PSI values group 1.
        group2: PSI values group 2.
        mean1: Mean of group 1.
        mean2: Mean of group 2.
        var1: Variance of group 1.
        var2: Variance of group 2.

    Returns:
        p-value from the moderated t-test.
    """
    n1 = len(group1)
    n2 = len(group2)

    if n1 < 2 or n2 < 2:
        # Cannot perform t-test with n=1; use simple threshold
        return 1.0 if abs(mean2 - mean1) < 0.1 else 0.05

    # Pool variance with prior (shrinkage toward common variance)
    pooled_var = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)

    # Empirical Bayes shrinkage: blend sample variance with prior
    prior_df = 3.0  # Prior degrees of freedom (tunable)
    prior_var = pooled_var  # Use pooled as prior

    shrunk_var1 = (prior_df * prior_var + (n1 - 1) * var1) / (prior_df + n1 - 1)
    shrunk_var2 = (prior_df * prior_var + (n2 - 1) * var2) / (prior_df + n2 - 1)

    # Moderated t-statistic
    se = math.sqrt(shrunk_var1 / n1 + shrunk_var2 / n2)

    if se < 1e-10:
        return 1.0 if abs(mean2 - mean1) < 1e-10 else 0.0

    t_stat = (mean2 - mean1) / se

    # Degrees of freedom (Welch-Satterthwaite with inflated df from prior)
    df = n1 + n2 - 2 + 2 * prior_df

    # p-value from t-distribution
    if HAS_SCIPY:
        p_value = float(2 * scipy_stats.t.sf(abs(t_stat), df))
    else:
        # Approximate using normal for large df
        p_value = 2 * (1 - _normal_cdf(abs(t_stat)))

    return min(p_value, 1.0)


def _wilcoxon_test(group1: list[float], group2: list[float]) -> float:
    """Wilcoxon rank-sum (Mann-Whitney U) test for PSI comparison.

    Args:
        group1: PSI values group 1.
        group2: PSI values group 2.

    Returns:
        p-value from the test.
    """
    if HAS_SCIPY:
        try:
            _, p_value = scipy_stats.mannwhitneyu(group1, group2, alternative="two-sided")
            return float(p_value)
        except ValueError:
            return 1.0

    # Pure Python Mann-Whitney U implementation
    n1 = len(group1)
    n2 = len(group2)
    combined = [(val, 0) for val in group1] + [(val, 1) for val in group2]
    combined.sort(key=lambda x: x[0])

    # Assign ranks (handle ties by average rank)
    ranks: list[float] = []
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        avg_rank = (i + j + 1) / 2  # 1-based average rank
        for _ in range(i, j):
            ranks.append(avg_rank)
        i = j

    # U statistic
    r1 = sum(ranks[k] for k in range(len(combined)) if combined[k][1] == 0)
    u1 = r1 - n1 * (n1 + 1) / 2

    # Normal approximation for p-value
    mu = n1 * n2 / 2
    sigma = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)

    if sigma < 1e-10:
        return 1.0

    z = (u1 - mu) / sigma
    p_value = 2 * (1 - _normal_cdf(abs(z)))
    return min(p_value, 1.0)


def _permutation_test(
    group1: list[float],
    group2: list[float],
    n_permutations: int = 10000,
) -> float:
    """Permutation test for differential splicing.

    Args:
        group1: PSI values group 1.
        group2: PSI values group 2.
        n_permutations: Number of permutations.

    Returns:
        Permutation p-value.
    """
    observed_diff = abs(sum(group2) / len(group2) - sum(group1) / len(group1))
    combined = group1 + group2
    n1 = len(group1)
    count_extreme = 0

    for _ in range(n_permutations):
        random.shuffle(combined)
        perm_g1 = combined[:n1]
        perm_g2 = combined[n1:]
        perm_diff = abs(sum(perm_g2) / len(perm_g2) - sum(perm_g1) / len(perm_g1))
        if perm_diff >= observed_diff:
            count_extreme += 1

    p_value = (count_extreme + 1) / (n_permutations + 1)
    return p_value


def _normal_cdf(x: float) -> float:
    """Approximate standard normal CDF using error function approximation.

    Args:
        x: z-score value.

    Returns:
        Cumulative probability P(Z <= x).
    """
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


# =============================================================================
# Novel Junction Discovery
# =============================================================================


def find_novel_junctions(
    junctions: list[dict],
    known_junctions: list[dict],
    max_distance: int = 10,
) -> list[dict]:
    """Identify unannotated splice junctions not present in a reference set.

    Compares detected junctions against a set of known/annotated junctions
    and reports those that are not within max_distance of any known junction.

    Args:
        junctions: Detected junctions from detect_splice_junctions().
        known_junctions: Reference set of annotated junctions, each with
            keys: chrom, start, end. Strand is optional.
        max_distance: Maximum coordinate distance (in bp) for a detected
            junction to be considered matching a known junction. Both start
            and end must be within this distance.

    Returns:
        List of novel junction dicts. Each dict is a copy of the input
        junction with an additional key:
            - nearest_known_distance (int): Distance to the nearest known
              junction (sum of start and end offsets), or -1 if no known
              junctions exist on the same chromosome.

    Example:
        >>> detected = [
        ...     {"chrom": "chr1", "start": 1000, "end": 2000, "strand": "+", "read_count": 5},
        ...     {"chrom": "chr1", "start": 3000, "end": 4000, "strand": "+", "read_count": 8},
        ... ]
        >>> known = [
        ...     {"chrom": "chr1", "start": 1002, "end": 1998},
        ... ]
        >>> novel = find_novel_junctions(detected, known, max_distance=10)
        >>> len(novel)  # Only the chr1:3000-4000 junction is novel
        1
    """
    if not junctions:
        return []

    if not known_junctions:
        logger.info("No known junctions provided, all junctions are novel")
        results = []
        for junc in junctions:
            novel = dict(junc)
            novel["nearest_known_distance"] = -1
            results.append(novel)
        return results

    # Index known junctions by chromosome for fast lookup
    known_by_chrom: dict[str, list[dict]] = defaultdict(list)
    for kj in known_junctions:
        known_by_chrom[kj["chrom"]].append(kj)

    novel: list[dict] = []

    for junc in junctions:
        chrom = junc["chrom"]
        junc_start = junc["start"]
        junc_end = junc["end"]

        chrom_known = known_by_chrom.get(chrom, [])

        if not chrom_known:
            result = dict(junc)
            result["nearest_known_distance"] = -1
            novel.append(result)
            continue

        # Find nearest known junction
        min_dist = float("inf")
        is_known = False

        for kj in chrom_known:
            start_dist = abs(junc_start - kj["start"])
            end_dist = abs(junc_end - kj["end"])

            if start_dist <= max_distance and end_dist <= max_distance:
                is_known = True
                break

            total_dist = start_dist + end_dist
            if total_dist < min_dist:
                min_dist = total_dist

        if not is_known:
            result = dict(junc)
            result["nearest_known_distance"] = int(min_dist) if min_dist != float("inf") else -1
            novel.append(result)

    logger.info(f"Found {len(novel)} novel junctions out of {len(junctions)} total " f"(max_distance={max_distance})")

    return novel
