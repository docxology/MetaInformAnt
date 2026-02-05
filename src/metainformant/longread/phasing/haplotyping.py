"""Haplotype phasing for long-read sequencing data.

Provides read-based phasing using heterozygous variants, phase block
construction, switch error computation, read haplotagging, and allele-specific
analysis from phased long reads.

The phasing algorithm assigns long reads to one of two haplotypes based on
their observed alleles at heterozygous variant sites, using a greedy
consistency-maximizing approach.

Optional dependencies:
    - numpy: Numerical computation (for allele-specific analysis)
    - scipy: Statistical tests (for allelic imbalance testing)
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


def phase_reads(
    variants: list[dict[str, Any]],
    reads: list[dict[str, Any]],
    method: str = "greedy",
) -> dict[str, Any]:
    """Phase long reads into haplotypes using heterozygous variant information.

    For each read, examines which alleles it carries at heterozygous variant
    sites. Reads are then assigned to one of two haplotypes using a greedy
    algorithm that maximizes consistency: reads sharing the same alleles at
    co-observed variants are placed on the same haplotype.

    Args:
        variants: List of heterozygous variant dicts, each containing:
            - ``chrom`` (str): Chromosome name.
            - ``position`` (int): 0-based position.
            - ``ref`` (str): Reference allele.
            - ``alt`` (str): Alternate allele.
        reads: List of read dicts, each containing:
            - ``read_id`` (str): Unique read identifier.
            - ``chrom`` (str): Mapped chromosome.
            - ``start`` (int): Alignment start position.
            - ``end`` (int): Alignment end position.
            - ``alleles`` (list[dict]): Observed alleles at variant sites,
              each with ``position`` (int) and ``allele`` (int, 0=ref, 1=alt).
        method: Phasing algorithm. Currently supports ``"greedy"``.

    Returns:
        Dictionary with keys:
            - ``haplotype_assignments``: Dict mapping read_id to haplotype
              (1 or 2).
            - ``phase_blocks``: List of phase block dicts with ``chrom``,
              ``start``, ``end``, ``n_variants``, ``n_reads``.
            - ``switch_error_rate``: Estimated switch error rate (0-1).
            - ``n_phased_reads``: Number of successfully phased reads.

    Raises:
        ValueError: If method is unrecognized.
    """
    if method != "greedy":
        raise ValueError(f"Unknown phasing method: {method}. Use 'greedy'.")

    if not variants or not reads:
        return {
            "haplotype_assignments": {},
            "phase_blocks": [],
            "switch_error_rate": 0.0,
            "n_phased_reads": 0,
        }

    # Index variants by (chrom, position)
    var_index: dict[tuple[str, int], int] = {}
    for i, v in enumerate(variants):
        var_index[(v["chrom"], v["position"])] = i

    n_variants = len(variants)

    # Extract read-variant observations
    read_observations: dict[str, list[tuple[int, int]]] = {}
    for read in reads:
        read_id = read.get("read_id", "")
        if not read_id:
            continue

        alleles = read.get("alleles", [])
        obs: list[tuple[int, int]] = []
        for a in alleles:
            key = (read.get("chrom", ""), a.get("position", -1))
            if key in var_index:
                obs.append((var_index[key], int(a.get("allele", 0))))

        if obs:
            read_observations[read_id] = obs

    if not read_observations:
        logger.warning("No read-variant observations found for phasing")
        return {
            "haplotype_assignments": {},
            "phase_blocks": [],
            "switch_error_rate": 0.0,
            "n_phased_reads": 0,
        }

    # Build variant co-observation graph
    edge_weights: dict[tuple[int, int], float] = defaultdict(float)

    for read_id, obs_list in read_observations.items():
        for i in range(len(obs_list)):
            for j in range(i + 1, len(obs_list)):
                vi, ai = obs_list[i]
                vj, aj = obs_list[j]
                if vi > vj:
                    vi, vj = vj, vi
                    ai, aj = aj, ai

                # Same allele = in-phase, different = anti-phase
                weight = 1.0 if ai == aj else -1.0
                edge_weights[(vi, vj)] += weight

    # Greedy max-cut for haplotype assignment
    variant_haplotype = _greedy_phase(n_variants, edge_weights)

    # Assign reads to haplotypes via majority vote
    haplotype_assignments: dict[str, int] = {}

    for read_id, obs_list in read_observations.items():
        h1_votes = 0
        h2_votes = 0
        for vi, ai in obs_list:
            if variant_haplotype[vi] == ai:
                h1_votes += 1
            else:
                h2_votes += 1

        if h1_votes > h2_votes:
            haplotype_assignments[read_id] = 1
        elif h2_votes > h1_votes:
            haplotype_assignments[read_id] = 2
        # Ties left unassigned

    # Build phase blocks from connected components
    phase_blocks = _build_blocks_from_graph(variants, edge_weights, variant_haplotype, read_observations)

    # Estimate switch error rate from edge consistency
    consistent = 0
    total = 0
    for (vi, vj), weight in edge_weights.items():
        same_hap = variant_haplotype[vi] == variant_haplotype[vj]
        total += 1
        if (weight > 0 and same_hap) or (weight < 0 and not same_hap):
            consistent += 1

    switch_error_rate = 1.0 - (consistent / total) if total > 0 else 0.0

    logger.info(
        "Phased %d reads into 2 haplotypes across %d phase blocks " "(switch_error_rate=%.4f)",
        len(haplotype_assignments),
        len(phase_blocks),
        switch_error_rate,
    )

    return {
        "haplotype_assignments": haplotype_assignments,
        "phase_blocks": phase_blocks,
        "switch_error_rate": switch_error_rate,
        "n_phased_reads": len(haplotype_assignments),
    }


def build_phase_blocks(
    phased_variants: list[dict[str, Any]],
    max_gap: int = 1000000,
) -> list[dict[str, Any]]:
    """Group phased variants into contiguous phase blocks.

    Variants on the same chromosome within ``max_gap`` bases of each other
    are grouped into the same phase block.

    Args:
        phased_variants: List of phased variant dicts, each containing:
            - ``chrom`` (str): Chromosome.
            - ``position`` (int): Genomic position.
            - ``haplotype`` (int): Assigned haplotype (1 or 2).
            - ``quality`` (float, optional): Phasing quality score.
        max_gap: Maximum gap between consecutive variants in the same block.

    Returns:
        List of phase block dicts, each containing:
            - ``chrom`` (str): Chromosome.
            - ``start`` (int): Block start position.
            - ``end`` (int): Block end position.
            - ``n_variants`` (int): Number of variants in the block.
            - ``n_reads`` (int): Set to 0 (not available without read data).
            - ``phase_quality`` (float): Mean quality of variants in block.
    """
    if not phased_variants:
        return []

    # Group by chromosome
    by_chrom: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for v in phased_variants:
        by_chrom[v["chrom"]].append(v)

    blocks: list[dict[str, Any]] = []

    for chrom in sorted(by_chrom.keys()):
        sorted_vars = sorted(by_chrom[chrom], key=lambda v: v["position"])

        current_block: list[dict[str, Any]] = [sorted_vars[0]]

        for i in range(1, len(sorted_vars)):
            gap = sorted_vars[i]["position"] - sorted_vars[i - 1]["position"]
            if gap > max_gap:
                if len(current_block) >= 2:
                    blocks.append(_summarize_block(chrom, current_block))
                current_block = [sorted_vars[i]]
            else:
                current_block.append(sorted_vars[i])

        if len(current_block) >= 2:
            blocks.append(_summarize_block(chrom, current_block))

    logger.info("Built %d phase blocks from %d variants", len(blocks), len(phased_variants))
    return blocks


def compute_switch_errors(
    truth_haplotypes: list[int],
    predicted_haplotypes: list[int],
) -> dict[str, Any]:
    """Compute switch error rate and phasing quality metrics.

    A switch error occurs at position i if the phase relationship between
    positions i-1 and i differs between truth and prediction. Formally,
    ``switch_error[i] = (truth[i] == truth[i-1]) != (pred[i] == pred[i-1])``.

    Args:
        truth_haplotypes: Ground truth haplotype assignments (list of 1 or 2).
        predicted_haplotypes: Predicted haplotype assignments (list of 1 or 2).

    Returns:
        Dictionary with keys:
            - ``switch_error_rate``: Fraction of switch positions with errors.
            - ``n_switches``: Total number of switch positions evaluated.
            - ``n_errors``: Number of switch errors.
            - ``longest_correct_run``: Longest consecutive correct switches.
            - ``hamming_distance``: Fraction of positions with different
              haplotype assignment (accounting for global flip).

    Raises:
        ValueError: If lists have different lengths.
    """
    if len(truth_haplotypes) != len(predicted_haplotypes):
        raise ValueError(
            f"Length mismatch: truth has {len(truth_haplotypes)}, " f"predicted has {len(predicted_haplotypes)}"
        )

    n = len(truth_haplotypes)
    if n < 2:
        return {
            "switch_error_rate": 0.0,
            "n_switches": 0,
            "n_errors": 0,
            "longest_correct_run": 0,
            "hamming_distance": 0.0,
        }

    # Count switches
    n_errors = 0
    n_switches = n - 1
    current_run = 0
    longest_run = 0

    for i in range(1, n):
        truth_same = truth_haplotypes[i] == truth_haplotypes[i - 1]
        pred_same = predicted_haplotypes[i] == predicted_haplotypes[i - 1]

        if truth_same != pred_same:
            n_errors += 1
            current_run = 0
        else:
            current_run += 1
            longest_run = max(longest_run, current_run)

    switch_error_rate = n_errors / n_switches if n_switches > 0 else 0.0

    # Hamming distance (accounting for possible global flip)
    direct_match = sum(1 for t, p in zip(truth_haplotypes, predicted_haplotypes) if t == p)
    flipped_match = sum(1 for t, p in zip(truth_haplotypes, predicted_haplotypes) if t != p)

    best_match = max(direct_match, flipped_match)
    hamming = 1.0 - (best_match / n)

    return {
        "switch_error_rate": switch_error_rate,
        "n_switches": n_switches,
        "n_errors": n_errors,
        "longest_correct_run": longest_run,
        "hamming_distance": hamming,
    }


def haplotag_reads(
    reads: list[dict[str, Any]],
    phase_blocks: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Tag reads with haplotype assignment based on phase blocks.

    For each read, finds the overlapping phase block and assigns a
    haplotype tag (HP) and phase set tag (PS). Reads that do not overlap
    any phase block receive HP=0, PS=0.

    Args:
        reads: List of read dicts with ``read_id``, ``chrom``, ``start``,
            ``end``, and optionally ``alleles`` (list of dicts with
            ``position`` and ``allele``).
        phase_blocks: List of phase block dicts with ``chrom``, ``start``,
            ``end``, and optionally ``variant_haplotypes`` (dict mapping
            position to haplotype assignment for that block).

    Returns:
        List of read dicts with added ``HP`` (int, 0/1/2) and ``PS``
        (int, phase set index) fields.
    """
    # Index phase blocks by chromosome
    block_index: dict[str, list[tuple[int, dict[str, Any]]]] = defaultdict(list)
    for idx, block in enumerate(phase_blocks):
        block_index[block["chrom"]].append((idx, block))

    # Sort blocks by start
    for chrom in block_index:
        block_index[chrom].sort(key=lambda x: x[1]["start"])

    tagged: list[dict[str, Any]] = []

    for read in reads:
        read_dict = dict(read)
        chrom = read_dict.get("chrom", "")
        read_start = read_dict.get("start", 0)
        read_end = read_dict.get("end", read_start)

        hp = 0
        ps = 0

        if chrom in block_index:
            for block_idx, block in block_index[chrom]:
                if read_start < block["end"] and read_end > block["start"]:
                    # Read overlaps block - assign haplotype from alleles
                    var_haps = block.get("variant_haplotypes", {})
                    alleles = read_dict.get("alleles", [])

                    h1_votes = 0
                    h2_votes = 0
                    for a in alleles:
                        pos = a.get("position", -1)
                        obs_allele = a.get("allele", -1)
                        if pos in var_haps:
                            if var_haps[pos] == obs_allele:
                                h1_votes += 1
                            else:
                                h2_votes += 1

                    if h1_votes > h2_votes:
                        hp = 1
                        ps = block_idx
                    elif h2_votes > h1_votes:
                        hp = 2
                        ps = block_idx
                    break  # Use first overlapping block

        read_dict["HP"] = hp
        read_dict["PS"] = ps
        tagged.append(read_dict)

    n_tagged = sum(1 for r in tagged if r["HP"] > 0)
    logger.info(
        "Haplotagged %d/%d reads with haplotype assignments",
        n_tagged,
        len(tagged),
    )

    return tagged


def allele_specific_analysis(
    reads: list[dict[str, Any]],
    haplotype_assignments: dict[str, int],
    region: dict[str, Any],
) -> dict[str, Any]:
    """Analyze allele-specific expression or methylation using phased reads.

    Partitions reads by haplotype and compares their signal (expression or
    methylation level) within the specified genomic region to detect allelic
    imbalance.

    Args:
        reads: List of read dicts with ``read_id``, ``chrom``, ``start``,
            ``end``, and ``signal`` (float, e.g., expression level or
            methylation probability).
        haplotype_assignments: Dict mapping read_id to haplotype (1 or 2).
        region: Region dict with ``chrom``, ``start``, ``end``.

    Returns:
        Dictionary with keys:
            - ``haplotype_1_signal``: Mean signal on haplotype 1.
            - ``haplotype_2_signal``: Mean signal on haplotype 2.
            - ``allelic_imbalance``: Absolute difference in signal between
              haplotypes, normalized by total.
            - ``p_value``: Statistical p-value from Mann-Whitney U test
              (or Welch's t-test fallback).

    Raises:
        ImportError: If numpy is not available.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required: uv pip install numpy")

    chrom = region.get("chrom", "")
    start = int(region.get("start", 0))
    end = int(region.get("end", 0))

    h1_signals: list[float] = []
    h2_signals: list[float] = []

    for read in reads:
        read_id = read.get("read_id", "")
        if read_id not in haplotype_assignments:
            continue

        read_chrom = read.get("chrom", "")
        read_start = read.get("start", 0)
        read_end = read.get("end", read_start)

        # Check region overlap
        if read_chrom != chrom or read_start >= end or read_end <= start:
            continue

        signal = float(read.get("signal", 0.0))
        hp = haplotype_assignments[read_id]

        if hp == 1:
            h1_signals.append(signal)
        elif hp == 2:
            h2_signals.append(signal)

    if not h1_signals or not h2_signals:
        return {
            "haplotype_1_signal": 0.0,
            "haplotype_2_signal": 0.0,
            "allelic_imbalance": 0.0,
            "p_value": 1.0,
        }

    h1_mean = float(np.mean(h1_signals))
    h2_mean = float(np.mean(h2_signals))

    # Allelic imbalance: |h1 - h2| / (h1 + h2)
    total = h1_mean + h2_mean
    imbalance = abs(h1_mean - h2_mean) / total if total > 0 else 0.0

    # Statistical test
    p_value = _test_allelic_difference(h1_signals, h2_signals)

    logger.info(
        "Allele-specific analysis in %s:%d-%d: H1=%.3f (n=%d), H2=%.3f (n=%d), " "imbalance=%.3f, p=%.4f",
        chrom,
        start,
        end,
        h1_mean,
        len(h1_signals),
        h2_mean,
        len(h2_signals),
        imbalance,
        p_value,
    )

    return {
        "haplotype_1_signal": h1_mean,
        "haplotype_2_signal": h2_mean,
        "allelic_imbalance": imbalance,
        "p_value": p_value,
    }


# ---------------------------------------------------------------------------
# Internal helper functions
# ---------------------------------------------------------------------------


def _greedy_phase(
    n_variants: int,
    edge_weights: dict[tuple[int, int], float],
) -> list[int]:
    """Greedy max-cut for bipartite variant phasing.

    Args:
        n_variants: Number of variants.
        edge_weights: Dict mapping (vi, vj) to weight.

    Returns:
        List of haplotype assignments (0 or 1) per variant.
    """
    assignment = [0] * n_variants

    adj: dict[int, list[tuple[int, float]]] = defaultdict(list)
    for (vi, vj), w in edge_weights.items():
        adj[vi].append((vj, w))
        adj[vj].append((vi, w))

    for node in range(1, n_variants):
        score_0 = 0.0
        score_1 = 0.0

        for neighbor, weight in adj[node]:
            if neighbor >= node:
                continue

            if weight > 0:
                # In-phase: same partition preferred
                if assignment[neighbor] == 0:
                    score_0 += weight
                else:
                    score_1 += weight
            else:
                # Anti-phase: different partition preferred
                if assignment[neighbor] == 0:
                    score_1 += abs(weight)
                else:
                    score_0 += abs(weight)

        assignment[node] = 0 if score_0 >= score_1 else 1

    return assignment


def _build_blocks_from_graph(
    variants: list[dict[str, Any]],
    edge_weights: dict[tuple[int, int], float],
    variant_haplotype: list[int],
    read_observations: dict[str, list[tuple[int, int]]],
) -> list[dict[str, Any]]:
    """Build phase blocks from the variant connection graph.

    Connected components of the variant graph become phase blocks.

    Args:
        variants: List of variant dicts.
        edge_weights: Edge weight dict.
        variant_haplotype: Per-variant haplotype assignment.
        read_observations: Per-read variant observations.

    Returns:
        List of phase block dicts.
    """
    n_variants = len(variants)

    # Find connected components
    adj: dict[int, set[int]] = defaultdict(set)
    for vi, vj in edge_weights.keys():
        adj[vi].add(vj)
        adj[vj].add(vi)

    visited: set[int] = set()
    components: list[list[int]] = []

    for node in range(n_variants):
        if node in visited or node not in adj:
            continue

        component: list[int] = []
        stack = [node]
        visited.add(node)

        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in adj[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)

        if len(component) >= 2:
            components.append(sorted(component))

    # Count reads per block
    variant_to_component: dict[int, int] = {}
    for comp_idx, component in enumerate(components):
        for vi in component:
            variant_to_component[vi] = comp_idx

    reads_per_component: dict[int, set[str]] = defaultdict(set)
    for read_id, obs_list in read_observations.items():
        for vi, _ in obs_list:
            if vi in variant_to_component:
                reads_per_component[variant_to_component[vi]].add(read_id)

    blocks: list[dict[str, Any]] = []
    for comp_idx, component in enumerate(components):
        comp_variants = [variants[i] for i in component]
        chrom = comp_variants[0].get("chrom", "")
        positions = [v["position"] for v in comp_variants]

        blocks.append(
            {
                "chrom": chrom,
                "start": min(positions),
                "end": max(positions),
                "n_variants": len(component),
                "n_reads": len(reads_per_component.get(comp_idx, set())),
            }
        )

    return blocks


def _summarize_block(
    chrom: str,
    block_variants: list[dict[str, Any]],
) -> dict[str, Any]:
    """Summarize a phase block from its constituent variants.

    Args:
        chrom: Chromosome name.
        block_variants: List of variant dicts in this block.

    Returns:
        Phase block summary dict.
    """
    positions = [v["position"] for v in block_variants]
    qualities = [v.get("quality", 1.0) for v in block_variants]
    mean_quality = sum(qualities) / len(qualities) if qualities else 0.0

    return {
        "chrom": chrom,
        "start": min(positions),
        "end": max(positions),
        "n_variants": len(block_variants),
        "n_reads": 0,
        "phase_quality": mean_quality,
    }


def _test_allelic_difference(
    signals_a: list[float],
    signals_b: list[float],
) -> float:
    """Test for significant difference between two signal distributions.

    Uses Mann-Whitney U test if scipy is available, otherwise Welch's t-test.

    Args:
        signals_a: Signal values from haplotype 1.
        signals_b: Signal values from haplotype 2.

    Returns:
        Two-sided p-value.
    """
    if HAS_SCIPY and scipy_stats is not None:
        try:
            _, p_val = scipy_stats.mannwhitneyu(signals_a, signals_b, alternative="two-sided")
            return float(p_val)
        except Exception:
            pass

    # Welch's t-test fallback
    n_a = len(signals_a)
    n_b = len(signals_b)

    if n_a < 2 or n_b < 2:
        return 1.0

    mean_a = sum(signals_a) / n_a
    mean_b = sum(signals_b) / n_b

    var_a = sum((x - mean_a) ** 2 for x in signals_a) / (n_a - 1)
    var_b = sum((x - mean_b) ** 2 for x in signals_b) / (n_b - 1)

    se = math.sqrt(var_a / n_a + var_b / n_b)
    if se == 0:
        return 1.0

    t_stat = abs(mean_a - mean_b) / se
    p_value = 2.0 * (1.0 - _normal_cdf(t_stat))
    return min(1.0, max(0.0, p_value))


def _normal_cdf(x: float) -> float:
    """Approximate standard normal CDF."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))
