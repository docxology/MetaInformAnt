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
np: Any  # None when numpy is unavailable (optional dependency)
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None

scipy_stats: Any  # None when scipy is unavailable (optional dependency)
try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None

# Numerical integration settings of the Bayesian model (see _site_posteriors)
_RATE_LOGIT_LIMIT = 60.0  # the site rate is integrated over logit(theta) in [-60, 60]
_LLR_LIMIT = 700.0  # larger |LLR| is certain evidence; clipping keeps the arithmetic finite
_MIN_PRIOR_PARAMETER = 0.25  # smaller Beta parameters put prior mass beyond the integration range
_WINDOW_NATS = 30.0  # nodes whose log density is further below the peak carry less than 1e-13 of the mass
_COARSE_NODES = 121
_FINE_NODES = 401
_MAX_REFINEMENTS = 4
_READ_CHUNK = 2048  # reads per block, bounding the size of the nodes x reads arrays


def call_methylation_from_signal(
    signal_data: list[dict[str, Any]],
    model: str = "threshold",
    threshold: float = 0.5,
    prior_alpha: float = 1.0,
    prior_beta: float = 1.0,
) -> list[dict[str, Any]]:
    """Call methylation status from Nanopore signal-level data.

    Each input entry contains signal-derived features including a
    log-likelihood ratio (LLR) for methylation vs unmethylation, that is
    ``log P(signal | methylated) - log P(signal | unmethylated)``. A call is
    ``methylated`` when the probability of methylation reaches ``threshold``.

    Supported models:
        - ``"threshold"``: the probability is the logistic function of the LLR,
          ``1 / (1 + exp(-LLR))``, which is the posterior under a flat 50 % prior.
        - ``"bayesian"``: a hierarchical Bayesian call. Reads at the same
          ``(chrom, position, strand)`` share an unknown methylation rate
          ``theta ~ Beta(prior_alpha, prior_beta)``; each read is methylated with
          probability ``theta`` and its LLR is the likelihood ratio of its
          signal. The probability of a read is its posterior
          ``P(methylated | LLRs of all reads at the site)``, the expectation of
          ``logistic(LLR + logit(theta))`` over the posterior of ``theta``,
          computed by numerical integration (accurate to better than 1e-6). Reads
          that disagree with the rest of a site are pulled toward it, a site
          with one read has probability ``logistic(LLR + logit(prior mean))``
          (identical to the threshold model for the default Beta(1, 1) prior),
          and results do not depend on the order of the reads. Entries without
          ``chrom`` and ``position`` are not pooled. Requires numpy.

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
        prior_alpha: Alpha of the Beta prior on the site methylation rate
            (``"bayesian"`` only). At least 0.25.
        prior_beta: Beta of the Beta prior on the site methylation rate
            (``"bayesian"`` only). At least 0.25. The prior mean is
            ``prior_alpha / (prior_alpha + prior_beta)``.

    Returns:
        List of call dictionaries, each containing all original fields plus:
            - ``methylated`` (bool): Whether the site is called methylated.
            - ``confidence`` (float): Classification confidence (0-1): twice
              the distance of the probability from 0.5. An entry whose LLR is
              NaN is called unmethylated with a NaN confidence and does not
              contribute to its site.

    Raises:
        ValueError: If model is unrecognized, an entry has no
            ``log_likelihood_ratio`` (a missing LLR used to count as 0.0, a probability
            of 0.5, which the default threshold calls methylated), or (``"bayesian"`` only)
            a prior parameter is below 0.25.
        ImportError: If model is ``"bayesian"`` and numpy is not installed.
    """
    if model not in ("threshold", "bayesian"):
        raise ValueError(f"Unknown model: {model}. Use 'threshold' or 'bayesian'.")

    llrs = [_required_llr(entry, index) for index, entry in enumerate(signal_data)]

    if model == "bayesian":
        if not (prior_alpha >= _MIN_PRIOR_PARAMETER and prior_beta >= _MIN_PRIOR_PARAMETER):
            raise ValueError(
                f"prior_alpha and prior_beta must be at least {_MIN_PRIOR_PARAMETER}, "
                f"got {prior_alpha} and {prior_beta}"
            )
        probabilities = _bayesian_probabilities(signal_data, llrs, prior_alpha, prior_beta)
    else:
        probabilities = [_llr_to_probability(llr) for llr in llrs]

    calls: list[dict[str, Any]] = []

    for entry, prob in zip(signal_data, probabilities):
        methylated = bool(prob >= threshold)
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
            ``call_methylation_from_signal``. Each must contain ``position`` and
            ``methylated`` (bool); ``chrom`` (default ``""``) and ``strand``
            (default ``"+"``) name the site when present.
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
        ValueError: If min_coverage < 1, or a call has no ``position`` or no ``methylated`` (they
            used to default to position 0 and unmethylated, which put the call at a site it was
            never made at, or counted a missing call as unmethylated).
    """
    if min_coverage < 1:
        raise ValueError("min_coverage must be at least 1")

    # Group calls by (chrom, position, strand)
    site_data: dict[tuple[str, int, str], list[bool]] = defaultdict(list)

    for index, call in enumerate(calls):
        for required in ("position", "methylated"):
            if call.get(required) is None:
                raise ValueError(f"calls[{index}] has no {required!r}")
        chrom = call.get("chrom", "")
        position = int(call["position"])
        strand = call.get("strand", "+")
        methylated = bool(call["methylated"])
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
    differences between the two samples using the Mann-Whitney U (Wilcoxon
    rank-sum) test, or a normal approximation of Welch's t statistic when scipy
    is unavailable or the exact test fails (a warning says so). The p-values are
    per region and are not corrected for the number of regions tested.

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

        An empty list, with a warning, when either sample has no sites or the two share no position.
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
            - ``n_reads``: Number of reads counted: those that call at least half
              of the CpGs seen in the region.
            - ``epialleles``: Dictionary mapping epiallele pattern string
              to its count.

    Raises:
        ValueError: If region is missing ``chrom``, ``start`` or ``end``.
    """
    missing = [key for key in ("chrom", "start", "end") if region.get(key) is None]
    if missing:
        raise ValueError(f"region is missing {', '.join(missing)}")
    chrom = region["chrom"]
    start = int(region["start"])
    end = int(region["end"])

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

    # Sort for percentiles (linear interpolation between order statistics, like numpy.percentile)
    sorted_freqs = sorted(freqs)
    median_val = _percentile(sorted_freqs, 50.0)

    variance = sum((f - global_meth) ** 2 for f in freqs) / max(n - 1, 1)
    std_val = math.sqrt(variance)

    distribution = {
        "min": sorted_freqs[0],
        "max": sorted_freqs[-1],
        "median": median_val,
        "std": std_val,
        "q25": _percentile(sorted_freqs, 25.0),
        "q75": _percentile(sorted_freqs, 75.0),
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


def _percentile(sorted_values: list[float], percent: float) -> float:
    """``percent``-th percentile of ascending ``sorted_values`` by linear interpolation (numpy's default)."""
    rank = (len(sorted_values) - 1) * percent / 100.0
    lower = math.floor(rank)
    upper = min(lower + 1, len(sorted_values) - 1)
    return sorted_values[lower] + (sorted_values[upper] - sorted_values[lower]) * (rank - lower)


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


def _required_llr(entry: dict[str, Any], index: int) -> float:
    """The ``log_likelihood_ratio`` of one signal entry; an entry without one cannot be called."""
    value = entry.get("log_likelihood_ratio")
    if value is None:
        raise ValueError(f"signal_data[{index}] has no 'log_likelihood_ratio'; the site cannot be called without it")
    return float(value)


def _bayesian_probabilities(
    signal_data: list[dict[str, Any]],
    llrs: list[float],
    prior_alpha: float,
    prior_beta: float,
) -> list[float]:
    """Posterior probability of methylation for every entry, pooling the reads of each site.

    Entries are grouped by ``(chrom, position, strand)``; an entry without ``chrom``
    and ``position`` forms its own group. An entry whose LLR is NaN gets a NaN
    probability and is left out of its site's evidence.

    Raises:
        ImportError: If numpy is not installed.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required for model='bayesian': uv pip install numpy")

    groups: dict[Any, list[int]] = {}
    for index, entry in enumerate(signal_data):
        if "chrom" in entry and "position" in entry:
            key: Any = (entry["chrom"], int(entry["position"]), entry.get("strand", "+"))
        else:
            key = index
        groups.setdefault(key, []).append(index)

    probabilities = [math.nan] * len(signal_data)
    for indices in groups.values():
        usable = [i for i in indices if not math.isnan(llrs[i])]
        if not usable:
            continue
        for i, probability in zip(usable, _site_posteriors([llrs[i] for i in usable], prior_alpha, prior_beta)):
            probabilities[i] = probability

    return probabilities


def _softplus(x: Any) -> Any:
    """``log(1 + exp(x))`` without overflow."""
    return np.logaddexp(0.0, x)


def _sigmoid(x: Any) -> Any:
    """``1 / (1 + exp(-x))`` without overflow."""
    return np.exp(-np.logaddexp(0.0, -x))


def _log_rate_density(u: Any, llr: Any, alpha: float, beta: float) -> Any:
    """Unnormalised log posterior density of ``u = logit(theta)`` for a site's methylation rate.

    The density of ``theta`` is ``theta^(alpha-1) (1-theta)^(beta-1) prod_i (1 - theta + theta e^LLR_i)``:
    the Beta prior times, for every read, its likelihood ratio averaged over whether the read is
    methylated. The change of variable to ``u`` adds the Jacobian ``theta (1-theta)``, which gives
    ``-alpha softplus(-u) - beta softplus(u) + sum_i [softplus(u + LLR_i) - softplus(u)]``.
    """
    total = -alpha * _softplus(-u) - beta * _softplus(u)
    baseline = _softplus(u)
    for start in range(0, len(llr), _READ_CHUNK):
        chunk = llr[start : start + _READ_CHUNK]
        total = total + (_softplus(u[:, None] + chunk[None, :]) - baseline[:, None]).sum(axis=1)
    return total


def _occupied_window(nodes: Any, log_density: Any) -> tuple[float, float]:
    """Interval of the grid that holds every node within 30 nats of the peak, padded by one node."""
    occupied = np.nonzero(log_density >= log_density.max() - _WINDOW_NATS)[0]
    return (
        float(nodes[max(int(occupied[0]) - 1, 0)]),
        float(nodes[min(int(occupied[-1]) + 1, len(nodes) - 1)]),
    )


def _site_posteriors(llrs: list[float], alpha: float, beta: float) -> list[float]:
    """Posterior probability that each read at one site is methylated, given every read's LLR.

    Model: ``theta ~ Beta(alpha, beta)``; read ``i`` is methylated with probability ``theta``;
    its LLR is the log likelihood ratio of its signal. The posterior of read ``i`` is the
    expectation of ``sigmoid(LLR_i + logit(theta))`` over the posterior of ``theta`` given all
    reads, its own included (read ``i`` is conditionally independent of the others given
    ``theta``, so this is exact). One read has the closed form
    ``sigmoid(LLR + log(alpha / beta))``.

    The posterior of ``u = logit(theta)`` is integrated on a uniform grid. A coarse grid over the
    whole range finds the window that holds its mass (every node within 30 nats of the peak, which
    is all but 1e-13 of it); a fine grid over that window integrates it, and while the window found
    on the fine grid is less than half as wide as the grid itself (many reads) the integration is
    repeated on it. Checked against exact enumeration of all methylation states, the probabilities
    agree to better than 1e-6 for every admissible prior (about 1e-13 for parameters of at least 0.5).
    """
    if len(llrs) == 1:
        return [_llr_to_probability(llrs[0] + math.log(alpha / beta))]

    llr = np.clip(np.asarray(llrs, dtype=float), -_LLR_LIMIT, _LLR_LIMIT)

    nodes = np.linspace(-_RATE_LOGIT_LIMIT, _RATE_LOGIT_LIMIT, _COARSE_NODES)
    log_density = _log_rate_density(nodes, llr, alpha, beta)
    low, high = _occupied_window(nodes, log_density)

    for _ in range(_MAX_REFINEMENTS):
        nodes = np.linspace(low, high, _FINE_NODES)
        log_density = _log_rate_density(nodes, llr, alpha, beta)
        narrower_low, narrower_high = _occupied_window(nodes, log_density)
        if narrower_high - narrower_low > 0.5 * (high - low):
            break
        low, high = narrower_low, narrower_high

    weights = np.exp(log_density - log_density.max())
    weights /= weights.sum()

    posteriors = np.empty(len(llr))
    for start in range(0, len(llr), _READ_CHUNK):
        chunk = llr[start : start + _READ_CHUNK]
        posteriors[start : start + len(chunk)] = weights @ _sigmoid(nodes[:, None] + chunk[None, :])

    return [float(min(1.0, max(0.0, p))) for p in posteriors]


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

    Uses the Mann-Whitney U (Wilcoxon rank-sum) test if scipy is available. Without scipy, or when scipy
    rejects the data (a warning says so), it falls back to a normal approximation of Welch's t statistic.

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
        except ValueError as exc:
            logger.warning("Mann-Whitney U test failed (%s); using the normal approximation of Welch's t instead", exc)

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
