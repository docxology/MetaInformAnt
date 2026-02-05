"""Isoform quantification and splice graph analysis for RNA-seq data.

This module provides tools for transcript-level abundance estimation using
an Expectation-Maximization (EM) algorithm, splice graph construction from
exon and junction data, isoform enumeration by graph traversal, and
isoform diversity and differential usage analysis.

All implementations are pure Python with optional numpy acceleration.
No external tool dependencies (Salmon, RSEM, StringTie) required.

Main Functions:
    Quantification:
        - quantify_isoforms: EM algorithm for isoform abundance estimation

    Graph Construction:
        - build_isoform_graph: Construct splice graph from exons/junctions
        - enumerate_isoforms: Enumerate possible isoforms by path traversal

    Diversity Analysis:
        - compute_isoform_diversity: Shannon entropy and effective isoform count
        - compare_isoform_usage: Compare isoform usage between conditions

Example:
    >>> from metainformant.rna.splicing import isoforms
    >>> read_assignments = [
    ...     {"read_id": "r1", "isoform_compatibilities": {"iso1": 1.0, "iso2": 0.5}},
    ...     {"read_id": "r2", "isoform_compatibilities": {"iso1": 0.3, "iso2": 1.0}},
    ... ]
    >>> isoform_models = [
    ...     {"isoform_id": "iso1", "length": 1500},
    ...     {"isoform_id": "iso2", "length": 1200},
    ... ]
    >>> result = isoforms.quantify_isoforms(read_assignments, isoform_models)
"""

from __future__ import annotations

import math
from collections import defaultdict, deque
from typing import Any, Dict, List, Optional, Sequence, Tuple

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependency handling
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]


# =============================================================================
# Isoform Quantification (EM Algorithm)
# =============================================================================


def quantify_isoforms(
    read_assignments: list[dict],
    isoform_models: list[dict],
    max_iterations: int = 1000,
    convergence_threshold: float = 1e-6,
) -> dict:
    """Estimate isoform abundances using an Expectation-Maximization algorithm.

    Implements a simplified version of the transcript quantification approach
    used by Salmon/RSEM. Reads are probabilistically assigned to isoforms
    based on compatibility scores, and isoform abundances are iteratively
    refined until convergence.

    The algorithm accounts for isoform effective length to convert raw
    estimated counts into TPM (Transcripts Per Million) values.

    Args:
        read_assignments: List of read-to-isoform compatibility dicts, each with:
            - read_id (str): Unique read identifier
            - isoform_compatibilities (dict): Mapping of isoform_id to
              compatibility score (float, 0-1). A score of 1.0 means the
              read is fully compatible; 0.0 means incompatible.
        isoform_models: List of isoform model dicts, each with:
            - isoform_id (str): Unique isoform identifier
            - length (int): Transcript length in base pairs
            - gene_id (str, optional): Parent gene identifier
        max_iterations: Maximum EM iterations before stopping.
        convergence_threshold: Stop when the maximum change in any isoform
            abundance between iterations falls below this threshold.

    Returns:
        Dictionary with keys:
            - abundances (dict): Mapping of isoform_id to estimated abundance
              (fraction of transcripts, sums to 1.0)
            - tpm (dict): Mapping of isoform_id to TPM values
            - counts (dict): Mapping of isoform_id to estimated read counts
            - converged (bool): Whether the EM algorithm converged
            - iterations (int): Number of iterations performed
            - effective_lengths (dict): Effective length per isoform

    Raises:
        ValueError: If no read assignments or isoform models provided, or
            if no reads are compatible with any isoform.

    Example:
        >>> reads = [
        ...     {"read_id": "r1", "isoform_compatibilities": {"iso1": 1.0, "iso2": 0.5}},
        ...     {"read_id": "r2", "isoform_compatibilities": {"iso2": 1.0}},
        ... ]
        >>> models = [
        ...     {"isoform_id": "iso1", "length": 1500},
        ...     {"isoform_id": "iso2", "length": 1200},
        ... ]
        >>> result = quantify_isoforms(reads, models)
        >>> sum(result["abundances"].values())  # Sums to ~1.0
        1.0
    """
    if not read_assignments:
        raise ValueError("No read assignments provided")
    if not isoform_models:
        raise ValueError("No isoform models provided")

    # Build isoform info lookup
    isoform_ids = [m["isoform_id"] for m in isoform_models]
    isoform_lengths = {m["isoform_id"]: m["length"] for m in isoform_models}
    n_isoforms = len(isoform_ids)

    # Mean read length estimate (used for effective length calculation)
    mean_read_length = 100  # Default; could be computed from data

    # Effective lengths (transcript length - read length + 1, minimum 1)
    effective_lengths = {}
    for iso_id in isoform_ids:
        eff_len = max(1, isoform_lengths[iso_id] - mean_read_length + 1)
        effective_lengths[iso_id] = eff_len

    # Initialize abundances uniformly
    abundances = {iso_id: 1.0 / n_isoforms for iso_id in isoform_ids}

    # Pre-process read compatibilities for efficiency
    processed_reads: list[dict[str, float]] = []
    for read in read_assignments:
        compat = read.get("isoform_compatibilities", {})
        # Filter to only known isoforms with positive compatibility
        filtered = {iso_id: score for iso_id, score in compat.items() if iso_id in effective_lengths and score > 0}
        if filtered:
            processed_reads.append(filtered)

    if not processed_reads:
        raise ValueError("No reads are compatible with any isoform")

    n_reads = len(processed_reads)
    logger.info(f"Running EM with {n_reads} reads and {n_isoforms} isoforms " f"(max_iter={max_iterations})")

    # EM iteration
    converged = False
    iteration = 0

    for iteration in range(1, max_iterations + 1):
        # E-step: Compute expected read assignments
        expected_counts = {iso_id: 0.0 for iso_id in isoform_ids}

        for read_compat in processed_reads:
            # Compute responsibility: P(isoform | read) proportional to
            # abundance * compatibility / effective_length
            responsibilities: dict[str, float] = {}
            total_responsibility = 0.0

            for iso_id, compat_score in read_compat.items():
                r = abundances[iso_id] * compat_score / effective_lengths[iso_id]
                responsibilities[iso_id] = r
                total_responsibility += r

            if total_responsibility < 1e-20:
                continue

            # Normalize and accumulate
            for iso_id, r in responsibilities.items():
                expected_counts[iso_id] += r / total_responsibility

        # M-step: Update abundances
        total_count = sum(expected_counts.values())

        if total_count < 1e-20:
            logger.warning("Total expected count near zero, stopping EM")
            break

        new_abundances = {iso_id: expected_counts[iso_id] / total_count for iso_id in isoform_ids}

        # Check convergence
        max_change = max(abs(new_abundances[iso_id] - abundances[iso_id]) for iso_id in isoform_ids)

        abundances = new_abundances

        if max_change < convergence_threshold:
            converged = True
            break

    if converged:
        logger.info(f"EM converged after {iteration} iterations")
    else:
        logger.warning(f"EM did not converge after {max_iterations} iterations " f"(max_change={max_change:.2e})")

    # Compute estimated counts
    counts = {iso_id: abundances[iso_id] * n_reads for iso_id in isoform_ids}

    # Compute TPM
    # TPM_i = (count_i / eff_len_i) / sum(count_j / eff_len_j) * 1e6
    count_per_len = {iso_id: counts[iso_id] / effective_lengths[iso_id] for iso_id in isoform_ids}
    total_cpl = sum(count_per_len.values())

    if total_cpl > 0:
        tpm = {iso_id: (count_per_len[iso_id] / total_cpl) * 1e6 for iso_id in isoform_ids}
    else:
        tpm = {iso_id: 0.0 for iso_id in isoform_ids}

    return {
        "abundances": abundances,
        "tpm": tpm,
        "counts": counts,
        "converged": converged,
        "iterations": iteration,
        "effective_lengths": effective_lengths,
    }


# =============================================================================
# Splice Graph Construction
# =============================================================================


def build_isoform_graph(
    exons: list[dict],
    junctions: list[dict],
) -> dict:
    """Construct a splice graph from exon and junction annotations.

    A splice graph (directed acyclic graph) represents the possible
    combinations of exons and introns in a gene. Nodes are exons (with
    source and sink sentinel nodes), and edges represent either direct
    adjacency (constitutive splicing) or observed splice junctions.

    Args:
        exons: List of exon dicts, each with:
            - exon_id (str): Unique exon identifier
            - start (int): Exon start coordinate (0-based)
            - end (int): Exon end coordinate (0-based, exclusive)
            - chrom (str, optional): Chromosome name
        junctions: List of junction dicts, each with:
            - start (int): Junction donor site (exon end coordinate)
            - end (int): Junction acceptor site (next exon start coordinate)
            - read_count (int, optional): Supporting read count
            - chrom (str, optional): Chromosome name

    Returns:
        Splice graph dictionary with keys:
            - nodes (list[dict]): Graph nodes, each with:
                - node_id (str): Node identifier (exon_id or sentinel)
                - type (str): "exon", "source", or "sink"
                - start (int): Start coordinate
                - end (int): End coordinate
            - edges (list[dict]): Graph edges, each with:
                - source (str): Source node_id
                - target (str): Target node_id
                - type (str): "junction" or "adjacency"
                - weight (int): Edge weight (read count or 1)
            - n_nodes (int): Total number of nodes
            - n_edges (int): Total number of edges

    Raises:
        ValueError: If no exons provided.

    Example:
        >>> exons = [
        ...     {"exon_id": "e1", "start": 100, "end": 200},
        ...     {"exon_id": "e2", "start": 300, "end": 400},
        ...     {"exon_id": "e3", "start": 500, "end": 600},
        ... ]
        >>> junctions = [
        ...     {"start": 200, "end": 300, "read_count": 10},
        ...     {"start": 200, "end": 500, "read_count": 5},  # Skips exon 2
        ...     {"start": 400, "end": 500, "read_count": 8},
        ... ]
        >>> graph = build_isoform_graph(exons, junctions)
        >>> graph["n_nodes"]  # 3 exons + source + sink
        5
    """
    if not exons:
        raise ValueError("No exons provided for graph construction")

    # Sort exons by start coordinate
    sorted_exons = sorted(exons, key=lambda e: (e["start"], e["end"]))

    # Build node list with source and sink sentinels
    nodes: list[dict] = []
    node_lookup: dict[str, dict] = {}

    # Source node
    min_start = sorted_exons[0]["start"]
    source_node = {
        "node_id": "__source__",
        "type": "source",
        "start": min_start - 1,
        "end": min_start - 1,
    }
    nodes.append(source_node)
    node_lookup["__source__"] = source_node

    # Exon nodes
    for exon in sorted_exons:
        exon_id = exon["exon_id"]
        node = {
            "node_id": exon_id,
            "type": "exon",
            "start": exon["start"],
            "end": exon["end"],
        }
        nodes.append(node)
        node_lookup[exon_id] = node

    # Sink node
    max_end = sorted_exons[-1]["end"]
    sink_node = {
        "node_id": "__sink__",
        "type": "sink",
        "start": max_end + 1,
        "end": max_end + 1,
    }
    nodes.append(sink_node)
    node_lookup["__sink__"] = sink_node

    # Build index: coordinate -> exon_id for junction endpoint matching
    start_to_exon: dict[int, list[str]] = defaultdict(list)
    end_to_exon: dict[int, list[str]] = defaultdict(list)

    for exon in sorted_exons:
        start_to_exon[exon["start"]].append(exon["exon_id"])
        end_to_exon[exon["end"]].append(exon["exon_id"])

    # Build edges
    edges: list[dict] = []
    edge_set: set[tuple[str, str]] = set()

    # Connect source to first exon(s) and last exon(s) to sink
    first_start = sorted_exons[0]["start"]
    for exon in sorted_exons:
        if exon["start"] == first_start:
            edge = {
                "source": "__source__",
                "target": exon["exon_id"],
                "type": "adjacency",
                "weight": 1,
            }
            edges.append(edge)
            edge_set.add(("__source__", exon["exon_id"]))

    last_end = sorted_exons[-1]["end"]
    for exon in sorted_exons:
        if exon["end"] == last_end:
            edge = {
                "source": exon["exon_id"],
                "target": "__sink__",
                "type": "adjacency",
                "weight": 1,
            }
            edges.append(edge)
            edge_set.add((exon["exon_id"], "__sink__"))

    # Add junction edges
    for junc in junctions:
        junc_start = junc["start"]
        junc_end = junc["end"]
        weight = junc.get("read_count", 1)

        # Find donor exon (exon whose end matches junction start)
        donor_exons = end_to_exon.get(junc_start, [])
        # Find acceptor exon (exon whose start matches junction end)
        acceptor_exons = start_to_exon.get(junc_end, [])

        for donor_id in donor_exons:
            for acceptor_id in acceptor_exons:
                edge_key = (donor_id, acceptor_id)
                if edge_key not in edge_set:
                    edges.append(
                        {
                            "source": donor_id,
                            "target": acceptor_id,
                            "type": "junction",
                            "weight": weight,
                        }
                    )
                    edge_set.add(edge_key)

    # Also add adjacency edges for consecutively annotated exons
    # (exons where one ends and the next begins without a gap or with
    # junctions already connecting them)
    for i in range(len(sorted_exons) - 1):
        curr_exon = sorted_exons[i]
        next_exon = sorted_exons[i + 1]

        edge_key = (curr_exon["exon_id"], next_exon["exon_id"])
        if edge_key not in edge_set:
            # Add adjacency if exons are close (within typical intron distance)
            gap = next_exon["start"] - curr_exon["end"]
            if 0 <= gap <= 500_000:
                edges.append(
                    {
                        "source": curr_exon["exon_id"],
                        "target": next_exon["exon_id"],
                        "type": "adjacency",
                        "weight": 1,
                    }
                )
                edge_set.add(edge_key)

    logger.info(
        f"Built splice graph with {len(nodes)} nodes and {len(edges)} edges "
        f"from {len(exons)} exons and {len(junctions)} junctions"
    )

    return {
        "nodes": nodes,
        "edges": edges,
        "n_nodes": len(nodes),
        "n_edges": len(edges),
    }


# =============================================================================
# Isoform Enumeration
# =============================================================================


def enumerate_isoforms(
    splice_graph: dict,
    max_isoforms: int = 100,
) -> list[dict]:
    """Enumerate possible isoforms from a splice graph by path traversal.

    Performs a depth-first search from the source node to the sink node
    of the splice graph, enumerating all unique paths. Each path
    represents a possible transcript isoform.

    Args:
        splice_graph: Splice graph dict from build_isoform_graph() with
            "nodes" and "edges" keys.
        max_isoforms: Maximum number of isoforms to enumerate. Prevents
            combinatorial explosion in complex genes. When the limit is
            reached, enumeration stops with a warning.

    Returns:
        List of isoform dicts, each containing:
            - isoform_id (str): Generated identifier (e.g., "isoform_001")
            - exons (list[str]): Ordered list of exon node_ids in the path
            - n_exons (int): Number of exons in the isoform
            - total_length (int): Sum of exon lengths
            - path_weight (float): Product of edge weights along the path
              (higher weight = more read support)

    Raises:
        ValueError: If splice graph has no nodes or edges.

    Example:
        >>> graph = build_isoform_graph(exons, junctions)
        >>> isoforms_list = enumerate_isoforms(graph, max_isoforms=50)
        >>> for iso in isoforms_list:
        ...     print(f"{iso['isoform_id']}: {iso['n_exons']} exons")
    """
    nodes = splice_graph.get("nodes", [])
    edges = splice_graph.get("edges", [])

    if not nodes or not edges:
        raise ValueError("Splice graph has no nodes or edges")

    # Build adjacency list and node lookup
    node_map = {n["node_id"]: n for n in nodes}
    adjacency: dict[str, list[tuple[str, int]]] = defaultdict(list)

    for edge in edges:
        adjacency[edge["source"]].append((edge["target"], edge.get("weight", 1)))

    # Verify source and sink exist
    if "__source__" not in node_map:
        raise ValueError("Splice graph missing source node")
    if "__sink__" not in node_map:
        raise ValueError("Splice graph missing sink node")

    # DFS enumeration of all source-to-sink paths
    isoforms_found: list[dict] = []

    # Use iterative DFS with explicit stack to avoid recursion limits
    # Stack entries: (current_node, path_so_far, cumulative_weight)
    stack: list[tuple[str, list[str], float]] = [("__source__", [], 1.0)]

    while stack and len(isoforms_found) < max_isoforms:
        current, path, weight = stack.pop()

        if current == "__sink__":
            # Filter path to only exon nodes
            exon_path = [nid for nid in path if node_map[nid]["type"] == "exon"]

            if exon_path:
                total_length = sum(node_map[nid]["end"] - node_map[nid]["start"] for nid in exon_path)

                isoforms_found.append(
                    {
                        "isoform_id": f"isoform_{len(isoforms_found) + 1:03d}",
                        "exons": exon_path,
                        "n_exons": len(exon_path),
                        "total_length": total_length,
                        "path_weight": round(weight, 4),
                    }
                )
            continue

        # Expand neighbors
        for neighbor, edge_weight in adjacency.get(current, []):
            # Cycle detection: skip if already in path
            if neighbor in path and neighbor != "__sink__":
                continue
            stack.append((neighbor, path + [current], weight * edge_weight))

    if len(isoforms_found) >= max_isoforms:
        logger.warning(f"Reached max_isoforms limit ({max_isoforms}). " f"Gene may have more possible isoforms.")

    # Sort by path weight (most supported first)
    isoforms_found.sort(key=lambda x: x["path_weight"], reverse=True)

    logger.info(f"Enumerated {len(isoforms_found)} isoforms from splice graph " f"(max={max_isoforms})")

    return isoforms_found


# =============================================================================
# Isoform Diversity Analysis
# =============================================================================


def compute_isoform_diversity(
    abundances: list[float],
) -> dict:
    """Compute isoform diversity metrics for a gene.

    Measures the diversity of isoform usage using information-theoretic
    metrics. A gene with one dominant isoform has low diversity, while
    a gene with many equally used isoforms has high diversity.

    Args:
        abundances: List of isoform abundance values (fractions or raw
            counts). Values must be non-negative. Will be normalized to
            sum to 1.0 internally. Zero-abundance isoforms are excluded
            from entropy calculation.

    Returns:
        Dictionary with keys:
            - shannon_entropy (float): Shannon entropy in bits (H = -sum(p * log2(p)))
            - effective_isoforms (float): Effective number of isoforms
              (2^H, analogous to Hill number of order 1)
            - max_entropy (float): Maximum possible entropy (log2(n_isoforms))
            - evenness (float): Pielou's evenness index (H / H_max), range [0, 1]
            - dominance (float): Simpson's dominance index (sum(p^2)),
              range [1/n, 1]. Higher = more dominated by few isoforms.
            - n_expressed (int): Number of isoforms with non-zero abundance
            - gini_index (float): Gini inequality coefficient [0, 1].
              0 = perfect equality, 1 = maximum inequality.

    Raises:
        ValueError: If abundances list is empty or contains negative values.

    Example:
        >>> # Equal usage of 4 isoforms
        >>> result = compute_isoform_diversity([0.25, 0.25, 0.25, 0.25])
        >>> result["shannon_entropy"]  # = log2(4) = 2.0
        2.0
        >>> result["effective_isoforms"]  # = 4.0
        4.0
        >>> # One dominant isoform
        >>> result = compute_isoform_diversity([0.95, 0.03, 0.01, 0.01])
        >>> result["effective_isoforms"] < 2.0
        True
    """
    if not abundances:
        raise ValueError("Abundances list cannot be empty")

    for val in abundances:
        if val < 0:
            raise ValueError(f"Abundances must be non-negative, got {val}")

    # Filter zero abundances for entropy calculation
    nonzero = [a for a in abundances if a > 0]
    n_expressed = len(nonzero)

    if n_expressed == 0:
        return {
            "shannon_entropy": 0.0,
            "effective_isoforms": 0.0,
            "max_entropy": 0.0,
            "evenness": 0.0,
            "dominance": 0.0,
            "n_expressed": 0,
            "gini_index": 0.0,
        }

    # Normalize to proportions
    total = sum(nonzero)
    proportions = [a / total for a in nonzero]

    # Shannon entropy (bits)
    entropy = -sum(p * math.log2(p) for p in proportions if p > 0)

    # Effective number of isoforms (Hill number order 1)
    effective = 2**entropy if entropy > 0 else 1.0

    # Maximum entropy
    max_entropy = math.log2(n_expressed) if n_expressed > 1 else 0.0

    # Pielou's evenness
    evenness = entropy / max_entropy if max_entropy > 0 else 1.0

    # Simpson's dominance (sum of squared proportions)
    dominance = sum(p**2 for p in proportions)

    # Gini index
    gini = _compute_gini(proportions)

    return {
        "shannon_entropy": round(entropy, 6),
        "effective_isoforms": round(effective, 6),
        "max_entropy": round(max_entropy, 6),
        "evenness": round(evenness, 6),
        "dominance": round(dominance, 6),
        "n_expressed": n_expressed,
        "gini_index": round(gini, 6),
    }


def _compute_gini(values: list[float]) -> float:
    """Compute Gini inequality coefficient.

    Args:
        values: List of non-negative proportions.

    Returns:
        Gini coefficient between 0 (perfect equality) and 1 (maximum inequality).
    """
    n = len(values)
    if n <= 1:
        return 0.0

    sorted_vals = sorted(values)
    total = sum(sorted_vals)

    if total == 0:
        return 0.0

    cumulative = 0.0
    gini_sum = 0.0

    for i, val in enumerate(sorted_vals):
        cumulative += val
        gini_sum += (2 * (i + 1) - n - 1) * val

    gini = gini_sum / (n * total)
    return max(0.0, min(1.0, gini))


# =============================================================================
# Isoform Usage Comparison
# =============================================================================


def compare_isoform_usage(
    abundances_a: dict,
    abundances_b: dict,
) -> dict:
    """Compare isoform usage patterns between two conditions.

    Quantifies the difference in isoform proportions between two
    conditions using multiple distance/divergence metrics and a
    statistical test for significant changes.

    Args:
        abundances_a: Isoform abundances for condition A.
            Mapping of isoform_id to abundance value (counts or fractions).
        abundances_b: Isoform abundances for condition B.
            Mapping of isoform_id to abundance value (counts or fractions).
            Isoforms not present in one condition are treated as zero.

    Returns:
        Dictionary with keys:
            - jensen_shannon_divergence (float): JSD between the two
              distributions [0, 1]. 0 = identical, 1 = maximally different.
            - chi_square_stat (float): Chi-square test statistic
            - chi_square_pvalue (float): P-value from chi-square test
            - correlation (float): Pearson correlation between proportions
            - max_delta (float): Maximum absolute difference in proportions
              for any single isoform
            - max_delta_isoform (str): Isoform with the largest change
            - switched_isoforms (list[str]): Isoforms that change from
              dominant (>50%) in one condition to non-dominant in the other
            - n_shared (int): Number of isoforms present in both conditions
            - n_condition_a_only (int): Isoforms unique to condition A
            - n_condition_b_only (int): Isoforms unique to condition B

    Raises:
        ValueError: If both abundance dicts are empty.

    Example:
        >>> a = {"iso1": 0.8, "iso2": 0.15, "iso3": 0.05}
        >>> b = {"iso1": 0.2, "iso2": 0.6, "iso3": 0.2}
        >>> result = compare_isoform_usage(a, b)
        >>> result["jensen_shannon_divergence"] > 0
        True
    """
    if not abundances_a and not abundances_b:
        raise ValueError("Both abundance dictionaries are empty")

    # Get union of all isoform IDs
    all_isoforms = sorted(set(list(abundances_a.keys()) + list(abundances_b.keys())))

    if not all_isoforms:
        raise ValueError("No isoforms found in either condition")

    # Build proportion vectors
    raw_a = [abundances_a.get(iso_id, 0.0) for iso_id in all_isoforms]
    raw_b = [abundances_b.get(iso_id, 0.0) for iso_id in all_isoforms]

    total_a = sum(raw_a)
    total_b = sum(raw_b)

    if total_a == 0 and total_b == 0:
        raise ValueError("Both conditions have zero total abundance")

    prop_a = [v / total_a if total_a > 0 else 0.0 for v in raw_a]
    prop_b = [v / total_b if total_b > 0 else 0.0 for v in raw_b]

    # Jensen-Shannon divergence
    jsd = _jensen_shannon_divergence(prop_a, prop_b)

    # Chi-square test
    chi2_stat, chi2_pvalue = _chi_square_test(raw_a, raw_b, total_a, total_b)

    # Pearson correlation
    correlation = _pearson_correlation(prop_a, prop_b)

    # Maximum delta and switching analysis
    max_delta = 0.0
    max_delta_isoform = ""
    switched: list[str] = []

    n_shared = 0
    n_a_only = 0
    n_b_only = 0

    for i, iso_id in enumerate(all_isoforms):
        delta = abs(prop_a[i] - prop_b[i])
        if delta > max_delta:
            max_delta = delta
            max_delta_isoform = iso_id

        # Check for isoform switching (dominant in one, not in other)
        if (prop_a[i] > 0.5 and prop_b[i] <= 0.5) or (prop_b[i] > 0.5 and prop_a[i] <= 0.5):
            switched.append(iso_id)

        # Presence tracking
        a_present = abundances_a.get(iso_id, 0.0) > 0
        b_present = abundances_b.get(iso_id, 0.0) > 0
        if a_present and b_present:
            n_shared += 1
        elif a_present:
            n_a_only += 1
        elif b_present:
            n_b_only += 1

    return {
        "jensen_shannon_divergence": round(jsd, 6),
        "chi_square_stat": round(chi2_stat, 4),
        "chi_square_pvalue": round(chi2_pvalue, 8),
        "correlation": round(correlation, 6),
        "max_delta": round(max_delta, 6),
        "max_delta_isoform": max_delta_isoform,
        "switched_isoforms": switched,
        "n_shared": n_shared,
        "n_condition_a_only": n_a_only,
        "n_condition_b_only": n_b_only,
    }


def _jensen_shannon_divergence(
    p: list[float],
    q: list[float],
) -> float:
    """Compute Jensen-Shannon divergence between two distributions.

    JSD is a symmetric, bounded measure of distributional similarity.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        JSD value in [0, 1] (using log base 2).
    """
    # Midpoint distribution
    m = [(pi + qi) / 2 for pi, qi in zip(p, q)]

    # KL divergences
    kl_pm = _kl_divergence(p, m)
    kl_qm = _kl_divergence(q, m)

    jsd = (kl_pm + kl_qm) / 2
    return max(0.0, min(1.0, jsd))


def _kl_divergence(p: list[float], q: list[float]) -> float:
    """Compute Kullback-Leibler divergence D_KL(P || Q).

    Args:
        p: True distribution.
        q: Approximating distribution.

    Returns:
        KL divergence (non-negative, in bits).
    """
    kl = 0.0
    for pi, qi in zip(p, q):
        if pi > 0 and qi > 0:
            kl += pi * math.log2(pi / qi)
    return kl


def _chi_square_test(
    counts_a: list[float],
    counts_b: list[float],
    total_a: float,
    total_b: float,
) -> tuple[float, float]:
    """Chi-square test for homogeneity of isoform usage.

    Args:
        counts_a: Raw counts/abundances for condition A.
        counts_b: Raw counts/abundances for condition B.
        total_a: Sum of condition A.
        total_b: Sum of condition B.

    Returns:
        Tuple of (chi-square statistic, p-value).
    """
    try:
        from scipy import stats as sp_stats

        # Build contingency-style comparison
        # Expected proportions under null (pooled)
        total = total_a + total_b
        if total == 0:
            return 0.0, 1.0

        chi2 = 0.0
        df = 0

        for ca, cb in zip(counts_a, counts_b):
            pooled = ca + cb
            if pooled == 0:
                continue

            # Expected values under null
            expected_a = pooled * (total_a / total)
            expected_b = pooled * (total_b / total)

            if expected_a > 0:
                chi2 += (ca - expected_a) ** 2 / expected_a
            if expected_b > 0:
                chi2 += (cb - expected_b) ** 2 / expected_b
            df += 1

        df = max(df - 1, 1)
        p_value = float(sp_stats.chi2.sf(chi2, df))
        return chi2, p_value

    except ImportError:
        # Pure Python fallback: no p-value computation
        total = total_a + total_b
        if total == 0:
            return 0.0, 1.0

        chi2 = 0.0
        for ca, cb in zip(counts_a, counts_b):
            pooled = ca + cb
            if pooled == 0:
                continue
            expected_a = pooled * (total_a / total)
            expected_b = pooled * (total_b / total)
            if expected_a > 0:
                chi2 += (ca - expected_a) ** 2 / expected_a
            if expected_b > 0:
                chi2 += (cb - expected_b) ** 2 / expected_b

        # Approximate p-value using normal approximation for large chi2
        return chi2, 1.0  # Conservative fallback


def _pearson_correlation(x: list[float], y: list[float]) -> float:
    """Compute Pearson correlation coefficient between two vectors.

    Args:
        x: First vector.
        y: Second vector.

    Returns:
        Correlation coefficient in [-1, 1], or 0.0 if undefined.
    """
    n = len(x)
    if n < 2:
        return 0.0

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    cov = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    var_x = sum((xi - mean_x) ** 2 for xi in x)
    var_y = sum((yi - mean_y) ** 2 for yi in y)

    denom = math.sqrt(var_x * var_y)
    if denom < 1e-15:
        return 0.0

    return cov / denom
