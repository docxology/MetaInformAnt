"""Methylation calling from long-read nanopore signal data.

Provides signal-level methylation calling using threshold or likelihood-based
models, per-site aggregation across reads, differentially methylated region
detection, single-read methylation pattern analysis, and summary statistics.

The calling pipeline:
1. Per-read methylation calls from signal data (log-likelihood ratios)
2. Aggregation to per-site methylation frequencies
3. DMR detection via segmentation and statistical testing
4. Epiallele analysis for single-read patterns

Optional dependencies:
    - numpy: Numerical computation
    - scipy: Statistical tests
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependencies
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]

try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None  # type: ignore[assignment]


def call_methylation_from_signal(
    signal_data: list[dict[str, Any]],
    model: str = "threshold",
    threshold: float = 0.5,
) -> list[dict[str, Any]]:
    """Call methylation status from Nanopore signal-level data.

    Each input entry contains signal-derived features including a
    log-likelihood ratio (LLR) for methylation vs unmethylation. The LLR
    is converted to a probability using the logistic function and
    compared against the threshold.

    Supported models:
        - ``"threshold"``: Simple LLR threshold classification.
        - ``"bayesian"``: Bayesian posterior with a Beta(1,1) prior on
          the methylation rate, integrating the LLR evidence.

    Args:
        signal_data: List of per-read signal dictionaries, each containing:
            - ``read_id`` (str): Read identifier.
            - ``chrom`` (str): Chromosome name.
            - ``position`` (int): 0-based genomic position.
            - ``strand`` (str): ``"+"`` or ``"-"``.
            - ``log_likelihood_ratio`` (float): Log-likelihood ratio for
              methylation (positive = more likely methylated).
        model: Classification model. ``"threshold"`` or ``"bayesian"``.
        threshold: Probability threshold for calling a site as methylated.

    Returns:
        List of call dictionaries, each containing all original fields plus:
            - ``methylated`` (bool): Whether the site is called methylated.
            - ``confidence`` (float): Classification confidence (0-1).

    Raises:
        ValueError: If model is unrecognized.
    """
    if model not in ("threshold", "bayesian"):
        raise ValueError(f"Unknown model: {model}. Use 'threshold' or 'bayesian'.")

    calls: list[dict[str, Any]] = []

    for entry in signal_data:
        llr = float(entry.get("log_likelihood_ratio", 0.0))

        # Convert LLR to probability
        prob = _llr_to_probability(llr)

        if model == "bayesian":
            # Bayesian posterior with Beta(1,1) prior
            # posterior = Beta(1 + prob, 1 + (1-prob))
            # posterior mean = (1 + prob) / (3)
            # simplified as weighted prior
            prob = (1.0 + prob) / 3.0 if prob > 0.5 else prob

        methylated = prob >= threshold
        confidence = abs(prob - 0.5) * 2.0  # 0 at boundary, 1 at extremes

        call = dict(entry)
        call["methylated"] = methylated
        call["confidence"] = confidence
        calls.append(call)

    n_methylated = sum(1 for c in calls if c["methylated"])
    logger.info(
        "Called methylation for %d sites: %d methylated, %d unmethylated " "(model=%s, threshold=%.2f)",
        len(calls),
        n_methylated,
        len(calls) - n_methylated,
        model,
        threshold,
    )

    return calls


def aggregate_methylation(
    calls: list[dict[str, Any]],
    min_coverage: int = 5,
) -> dict[str, Any]:
    """Aggregate per-read methylation calls to per-site methylation frequencies.

    Groups calls by genomic position (chrom, position, strand) and computes
    the methylation frequency (fraction of reads calling the site methylated)
    at each position with sufficient coverage.

    Args:
        calls: List of per-read call dictionaries as returned by
            ``call_methylation_from_signal``. Each must contain ``chrom``,
            ``position``, ``strand``, ``methylated`` (bool).
        min_coverage: Minimum number of reads covering a site for it to
            be included in the output.

    Returns:
        Dictionary with keys:
            - ``sites``: List of site dicts, each containing ``chrom``,
              ``position``, ``strand``, ``methylation_freq`` (float, 0-1),
              ``coverage`` (int), ``n_methylated`` (int),
              ``n_unmethylated`` (int).
            - ``n_sites``: Total number of sites passing coverage filter.
            - ``mean_coverage``: Mean coverage across reported sites.

    Raises:
        ValueError: If min_coverage < 1.
    """
    if min_coverage < 1:
        raise ValueError("min_coverage must be at least 1")

    # Group calls by (chrom, position, strand)
    site_data: dict[tuple[str, int, str], list[bool]] = defaultdict(list)

    for call in calls:
        chrom = call.get("chrom", "")
        position = int(call.get("position", 0))
        strand = call.get("strand", "+")
        methylated = bool(call.get("methylated", False))
        site_data[(chrom, position, strand)].append(methylated)

    sites: list[dict[str, Any]] = []
    total_coverage = 0

    for (chrom, position, strand), meth_calls in sorted(site_data.items()):
        coverage = len(meth_calls)
        if coverage < min_coverage:
            continue

        n_methylated = sum(meth_calls)
        n_unmethylated = coverage - n_methylated
        methylation_freq = n_methylated / coverage

        sites.append(
            {
                "chrom": chrom,
                "position": position,
                "strand": strand,
                "methylation_freq": methylation_freq,
                "coverage": coverage,
                "n_methylated": n_methylated,
                "n_unmethylated": n_unmethylated,
            }
        )
        total_coverage += coverage

    n_sites = len(sites)
    mean_cov = total_coverage / n_sites if n_sites > 0 else 0.0

    logger.info(
        "Aggregated %d calls to %d sites (min_coverage=%d, mean_coverage=%.1f)",
        len(calls),
        n_sites,
        min_coverage,
        mean_cov,
    )

    return {
        "sites": sites,
        "n_sites": n_sites,
        "mean_coverage": mean_cov,
    }


def detect_dmrs(
    methylation_a: dict[str, Any],
    methylation_b: dict[str, Any],
    min_cpgs: int = 3,
    min_diff: float = 0.2,
    max_distance: int = 1000,
) -> list[dict[str, Any]]:
    """Detect differentially methylated regions between two samples.

    Segments CpG sites into candidate regions based on genomic proximity
    (``max_distance``), then tests each region for significant methylation
    differences between the two samples using a Wilcoxon rank-sum test
    (or Fisher's exact test as fallback).

    Args:
        methylation_a: Aggregated methylation from sample A as returned by
            ``aggregate_methylation``. Must contain ``sites`` list.
        methylation_b: Aggregated methylation from sample B.
        min_cpgs: Minimum number of CpG sites in a region to test.
        min_diff: Minimum absolute mean methylation difference to report.
        max_distance: Maximum distance between consecutive CpGs to belong
            to the same candidate region.

    Returns:
        List of DMR dictionaries, each containing:
            - ``chrom`` (str): Chromosome.
            - ``start`` (int): Region start position.
            - ``end`` (int): Region end position.
            - ``n_cpgs`` (int): Number of CpGs in the region.
            - ``mean_diff`` (float): Mean methylation difference (B - A).
            - ``p_value`` (float): Statistical p-value.
            - ``direction`` (str): ``"hyper"`` if B > A, ``"hypo"`` otherwise.

    Raises:
        ValueError: If input samples have no sites.
    """
    sites_a = methylation_a.get("sites", [])
    sites_b = methylation_b.get("sites", [])

    if not sites_a or not sites_b:
        logger.warning("One or both samples have no sites for DMR detection")
        return []

    # Index sites by (chrom, position)
    idx_a: dict[tuple[str, int], dict[str, Any]] = {}
    for site in sites_a:
        idx_a[(site["chrom"], site["position"])] = site

    idx_b: dict[tuple[str, int], dict[str, Any]] = {}
    for site in sites_b:
        idx_b[(site["chrom"], site["position"])] = site

    # Find shared positions
    shared_keys = sorted(set(idx_a.keys()) & set(idx_b.keys()))
    if not shared_keys:
        logger.warning("No shared CpG positions between samples")
        return []

    # Segment shared CpGs into candidate regions
    regions = _segment_cpgs(shared_keys, max_distance)

    dmrs: list[dict[str, Any]] = []

    for region_keys in regions:
        if len(region_keys) < min_cpgs:
            continue

        chrom = region_keys[0][0]
        start = region_keys[0][1]
        end = region_keys[-1][1]

        # Collect methylation frequencies
        freqs_a = [idx_a[k]["methylation_freq"] for k in region_keys]
        freqs_b = [idx_b[k]["methylation_freq"] for k in region_keys]

        mean_a = sum(freqs_a) / len(freqs_a)
        mean_b = sum(freqs_b) / len(freqs_b)
        mean_diff = mean_b - mean_a

        if abs(mean_diff) < min_diff:
            continue

        # Statistical test
        p_value = _test_methylation_difference(freqs_a, freqs_b)

        direction = "hyper" if mean_diff > 0 else "hypo"

        dmrs.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "n_cpgs": len(region_keys),
                "mean_diff": mean_diff,
                "p_value": p_value,
                "direction": direction,
            }
        )

    logger.info(
        "Detected %d DMRs from %d shared CpG positions (%d candidate regions)",
        len(dmrs),
        len(shared_keys),
        len(regions),
    )

    return dmrs


def methylation_pattern_analysis(
    calls: list[dict[str, Any]],
    region: dict[str, Any],
) -> dict[str, Any]:
    """Analyze methylation patterns at single-read level within a region.

    Examines per-read methylation patterns (epialleles) within the specified
    genomic region. Computes co-methylation between CpG sites, identifies
    distinct epialleles, and measures pattern entropy to quantify
    epigenetic heterogeneity.

    Args:
        calls: Per-read methylation calls, each containing ``read_id``,
            ``chrom``, ``position``, ``methylated`` (bool).
        region: Region dictionary with ``chrom`` (str), ``start`` (int),
            ``end`` (int).

    Returns:
        Dictionary with keys:
            - ``patterns``: List of unique epiallele strings (e.g.,
              ``"MMUUM"`` for methylated/unmethylated pattern).
            - ``entropy``: Shannon entropy of the epiallele distribution.
            - ``n_reads``: Number of reads covering the region.
            - ``epialleles``: Dictionary mapping epiallele pattern string
              to its count.

    Raises:
        ValueError: If region is missing required keys.
    """
    chrom = region.get("chrom", "")
    start = int(region.get("start", 0))
    end = int(region.get("end", 0))

    # Filter calls to region
    region_calls: list[dict[str, Any]] = []
    for call in calls:
        if call.get("chrom", "") == chrom and start <= call.get("position", -1) < end:
            region_calls.append(call)

    if not region_calls:
        return {
            "patterns": [],
            "entropy": 0.0,
            "n_reads": 0,
            "epialleles": {},
        }

    # Group by read_id
    reads_data: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for call in region_calls:
        read_id = call.get("read_id", "")
        if read_id:
            reads_data[read_id].append(call)

    # Get sorted CpG positions in region
    all_positions = sorted(set(call["position"] for call in region_calls if "position" in call))
    pos_to_idx = {pos: i for i, pos in enumerate(all_positions)}
    n_cpgs = len(all_positions)

    if n_cpgs == 0:
        return {
            "patterns": [],
            "entropy": 0.0,
            "n_reads": len(reads_data),
            "epialleles": {},
        }

    # Build epiallele patterns
    epiallele_counts: dict[str, int] = defaultdict(int)
    n_reads = 0

    for read_id, read_calls in reads_data.items():
        # Build pattern string for this read
        pattern = ["?"] * n_cpgs
        for call in read_calls:
            pos = call.get("position", -1)
            if pos in pos_to_idx:
                idx = pos_to_idx[pos]
                pattern[idx] = "M" if call.get("methylated", False) else "U"

        # Only count reads covering at least half the CpGs
        covered = sum(1 for p in pattern if p != "?")
        if covered >= max(1, n_cpgs // 2):
            pattern_str = "".join(pattern)
            epiallele_counts[pattern_str] += 1
            n_reads += 1

    # Compute Shannon entropy of epiallele distribution
    entropy = 0.0
    total = sum(epiallele_counts.values())
    if total > 0:
        for count in epiallele_counts.values():
            freq = count / total
            if freq > 0:
                entropy -= freq * math.log2(freq)

    patterns = sorted(epiallele_counts.keys(), key=lambda k: -epiallele_counts[k])

    logger.info(
        "Methylation pattern analysis: %d reads, %d unique epialleles, " "entropy=%.3f in %s:%d-%d",
        n_reads,
        len(patterns),
        entropy,
        chrom,
        start,
        end,
    )

    return {
        "patterns": patterns,
        "entropy": entropy,
        "n_reads": n_reads,
        "epialleles": dict(epiallele_counts),
    }


def compute_methylation_stats(
    sites: list[dict[str, Any]],
) -> dict[str, Any]:
    """Compute summary statistics for aggregated methylation data.

    Calculates global methylation level, distribution of methylation
    frequencies, coverage statistics, and optional context-stratified
    proportions (CpG, CHG, CHH).

    Args:
        sites: List of aggregated site dictionaries as returned by the
            ``sites`` field of ``aggregate_methylation``. Each should
            contain ``methylation_freq``, ``coverage``, and optionally
            ``context``.

    Returns:
        Dictionary with keys:
            - ``global_methylation``: Mean methylation across all sites.
            - ``n_sites``: Total number of sites.
            - ``mean_coverage``: Mean read coverage.
            - ``methylation_distribution``: Dictionary with ``min``, ``max``,
              ``median``, ``std``, ``q25``, ``q75``.
            - ``context_proportions``: Dictionary mapping context string
              (e.g., ``"CpG"``) to its fraction of total sites.
            - ``highly_methylated_fraction``: Fraction of sites with
              methylation >= 0.8.
            - ``lowly_methylated_fraction``: Fraction of sites with
              methylation <= 0.2.

    Returns:
        Empty summary dict if no sites provided.
    """
    if not sites:
        return {
            "global_methylation": 0.0,
            "n_sites": 0,
            "mean_coverage": 0.0,
            "methylation_distribution": {},
            "context_proportions": {},
            "highly_methylated_fraction": 0.0,
            "lowly_methylated_fraction": 0.0,
        }

    freqs = [s["methylation_freq"] for s in sites]
    coverages = [s.get("coverage", 0) for s in sites]
    n = len(freqs)

    global_meth = sum(freqs) / n
    mean_cov = sum(coverages) / n if coverages else 0.0

    # Sort for percentiles
    sorted_freqs = sorted(freqs)
    q25_idx = max(0, int(n * 0.25) - 1)
    q75_idx = min(n - 1, int(n * 0.75))
    median_idx = n // 2

    if n % 2 == 0 and n > 1:
        median_val = (sorted_freqs[median_idx - 1] + sorted_freqs[median_idx]) / 2.0
    else:
        median_val = sorted_freqs[median_idx]

    variance = sum((f - global_meth) ** 2 for f in freqs) / max(n - 1, 1)
    std_val = math.sqrt(variance)

    distribution = {
        "min": sorted_freqs[0],
        "max": sorted_freqs[-1],
        "median": median_val,
        "std": std_val,
        "q25": sorted_freqs[q25_idx],
        "q75": sorted_freqs[q75_idx],
    }

    # Context proportions
    context_counts: dict[str, int] = defaultdict(int)
    for site in sites:
        ctx = site.get("context", "unknown")
        context_counts[ctx] += 1

    context_proportions = {ctx: count / n for ctx, count in context_counts.items()}

    # Methylation level fractions
    highly_meth = sum(1 for f in freqs if f >= 0.8) / n
    lowly_meth = sum(1 for f in freqs if f <= 0.2) / n

    logger.info(
        "Methylation stats: %d sites, global=%.3f, mean_coverage=%.1f",
        n,
        global_meth,
        mean_cov,
    )

    return {
        "global_methylation": global_meth,
        "n_sites": n,
        "mean_coverage": mean_cov,
        "methylation_distribution": distribution,
        "context_proportions": context_proportions,
        "highly_methylated_fraction": highly_meth,
        "lowly_methylated_fraction": lowly_meth,
    }


# ---------------------------------------------------------------------------
# Internal helper functions
# ---------------------------------------------------------------------------


def _llr_to_probability(llr: float) -> float:
    """Convert log-likelihood ratio to probability.

    Uses the logistic function: ``p = 1 / (1 + exp(-llr))``.

    Args:
        llr: Log-likelihood ratio (positive favors methylation).

    Returns:
        Probability between 0 and 1.
    """
    if llr > 500:
        return 1.0
    if llr < -500:
        return 0.0
    return 1.0 / (1.0 + math.exp(-llr))


def _segment_cpgs(
    positions: list[tuple[str, int]],
    max_distance: int,
) -> list[list[tuple[str, int]]]:
    """Segment sorted CpG positions into contiguous regions.

    Consecutive CpGs on the same chromosome within ``max_distance`` are
    grouped into the same region.

    Args:
        positions: Sorted list of (chrom, position) tuples.
        max_distance: Maximum gap between consecutive CpGs in the same region.

    Returns:
        List of regions, each a list of (chrom, position) tuples.
    """
    if not positions:
        return []

    regions: list[list[tuple[str, int]]] = [[positions[0]]]

    for i in range(1, len(positions)):
        prev_chrom, prev_pos = positions[i - 1]
        curr_chrom, curr_pos = positions[i]

        if curr_chrom == prev_chrom and (curr_pos - prev_pos) <= max_distance:
            regions[-1].append(positions[i])
        else:
            regions.append([positions[i]])

    return regions


def _test_methylation_difference(
    freqs_a: list[float],
    freqs_b: list[float],
) -> float:
    """Test for methylation difference between two groups.

    Uses Wilcoxon rank-sum test if scipy is available, otherwise falls back
    to a simple t-test approximation.

    Args:
        freqs_a: Methylation frequencies from sample A.
        freqs_b: Methylation frequencies from sample B.

    Returns:
        Two-sided p-value.
    """
    if HAS_SCIPY and scipy_stats is not None:
        try:
            stat, p_val = scipy_stats.mannwhitneyu(freqs_a, freqs_b, alternative="two-sided")
            return float(p_val)
        except Exception:
            pass

    # Fallback: Welch's t-test
    n_a = len(freqs_a)
    n_b = len(freqs_b)

    if n_a < 2 or n_b < 2:
        return 1.0

    mean_a = sum(freqs_a) / n_a
    mean_b = sum(freqs_b) / n_b

    var_a = sum((x - mean_a) ** 2 for x in freqs_a) / (n_a - 1)
    var_b = sum((x - mean_b) ** 2 for x in freqs_b) / (n_b - 1)

    se = math.sqrt(var_a / n_a + var_b / n_b)
    if se == 0:
        return 1.0

    t_stat = abs(mean_a - mean_b) / se

    # Approximate p-value from t-distribution using normal approx
    # for large samples
    z = t_stat
    p_value = 2.0 * (1.0 - _normal_cdf(z))
    return min(1.0, max(0.0, p_value))


def _normal_cdf(x: float) -> float:
    """Approximate standard normal CDF using Abramowitz and Stegun formula."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))
