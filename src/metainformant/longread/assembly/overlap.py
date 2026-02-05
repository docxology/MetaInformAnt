"""Overlap computation for long-read assembly using minimizer sketching.

Implements minimizer-based sequence sketching and all-vs-all overlap
detection for de novo assembly. The minimizer approach (as used in minimap2)
provides efficient seed-based overlap finding.

Algorithm overview:
1. Compute minimizer sketches for all reads
2. Build a hash index of minimizer positions
3. Find read pairs sharing minimizers
4. Score overlaps by chaining shared minimizers
5. Build overlap graph and filter contained reads

Optional dependencies:
    - numpy: For efficient numerical computation
"""

from __future__ import annotations

import hashlib
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

try:
    import numpy as np  # type: ignore[import-untyped]
except ImportError:
    np = None  # type: ignore[assignment]


@dataclass
class Minimizer:
    """A minimizer (k-mer hash, position, strand) from a sequence.

    Attributes:
        hash_value: Integer hash of the k-mer.
        position: Position in the sequence (0-based).
        is_reverse: Whether the minimizer is from the reverse complement.
        kmer: The actual k-mer string (optional, for debugging).
    """

    hash_value: int
    position: int
    is_reverse: bool = False
    kmer: str = ""


@dataclass
class Overlap:
    """An overlap between two reads.

    Attributes:
        query_name: Name of the query read.
        query_length: Length of the query read.
        query_start: Start of the overlap on the query (0-based).
        query_end: End of the overlap on the query.
        target_name: Name of the target read.
        target_length: Length of the target read.
        target_start: Start of the overlap on the target.
        target_end: End of the overlap on the target.
        strand: Relative strand ('+' = same, '-' = reverse complement).
        num_matches: Number of matching minimizers.
        overlap_length: Length of the overlap region.
        identity: Estimated sequence identity.
        is_contained: Whether the query is fully contained in the target.
    """

    query_name: str
    query_length: int
    query_start: int
    query_end: int
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    strand: str = "+"
    num_matches: int = 0
    overlap_length: int = 0
    identity: float = 0.0
    is_contained: bool = False


@dataclass
class OverlapGraph:
    """Graph representation of read overlaps for assembly.

    Attributes:
        nodes: List of read names.
        edges: List of (query, target, overlap) tuples.
        adjacency: Adjacency list mapping read_name -> list of (neighbor, overlap).
        num_nodes: Number of nodes.
        num_edges: Number of edges.
    """

    nodes: list[str] = field(default_factory=list)
    edges: list[tuple[str, str, Overlap]] = field(default_factory=list)
    adjacency: dict[str, list[tuple[str, Overlap]]] = field(default_factory=dict)
    num_nodes: int = 0
    num_edges: int = 0


def minimizer_sketch(
    sequence: str,
    k: int = 15,
    w: int = 10,
) -> list[Minimizer]:
    """Compute the minimizer sketch of a sequence.

    Minimizers are the smallest k-mer hashes within each window of w
    consecutive k-mers. This provides a compressed representation that
    maintains locality-sensitive properties for overlap detection.

    The algorithm slides a window of w k-mers across the sequence and
    selects the k-mer with the minimum hash value in each window. When
    the minimum changes, a new minimizer is emitted.

    Uses canonical k-mers (minimum of forward and reverse complement hash)
    for strand-agnostic matching.

    Args:
        sequence: DNA sequence string.
        k: K-mer size (default 15).
        w: Window size (number of k-mers per window, default 10).

    Returns:
        List of Minimizer objects representing the sketch.
    """
    if not sequence or len(sequence) < k:
        return []

    seq_upper = sequence.upper()
    minimizers: list[Minimizer] = []

    # Compute hash for each k-mer position
    kmer_hashes: list[tuple[int, int, bool, str]] = []  # (hash, position, is_rc, kmer)

    for i in range(len(seq_upper) - k + 1):
        kmer = seq_upper[i : i + k]

        # Skip k-mers with N
        if "N" in kmer:
            kmer_hashes.append((2**64, i, False, kmer))
            continue

        # Canonical k-mer: min of forward and reverse complement hash
        fwd_hash = _hash_kmer(kmer)
        rc_kmer = _reverse_complement(kmer)
        rc_hash = _hash_kmer(rc_kmer)

        if fwd_hash <= rc_hash:
            kmer_hashes.append((fwd_hash, i, False, kmer))
        else:
            kmer_hashes.append((rc_hash, i, True, rc_kmer))

    if len(kmer_hashes) < w:
        # Sequence too short for windowed minimizers, return all
        for h, pos, is_rc, kmer in kmer_hashes:
            if h < 2**64:
                minimizers.append(Minimizer(hash_value=h, position=pos, is_reverse=is_rc, kmer=kmer))
        return minimizers

    # Slide window and extract minimizers
    prev_min_hash = -1
    prev_min_pos = -1

    for window_start in range(len(kmer_hashes) - w + 1):
        window = kmer_hashes[window_start : window_start + w]

        # Find minimum hash in window
        min_hash = 2**64
        min_pos = -1
        min_rc = False
        min_kmer = ""

        for h, pos, is_rc, kmer in window:
            if h < min_hash:
                min_hash = h
                min_pos = pos
                min_rc = is_rc
                min_kmer = kmer

        # Only emit if this is a new minimizer
        if min_pos != prev_min_pos:
            minimizers.append(Minimizer(
                hash_value=min_hash,
                position=min_pos,
                is_reverse=min_rc,
                kmer=min_kmer,
            ))
            prev_min_hash = min_hash
            prev_min_pos = min_pos

    return minimizers


def find_overlaps(
    reads: Sequence[dict[str, Any] | str],
    min_overlap: int = 2000,
    k: int = 15,
    w: int = 10,
    min_minimizer_matches: int = 3,
    max_overhang: int = 1000,
) -> list[Overlap]:
    """Find all pairwise overlaps between reads using minimizer sketching.

    Algorithm:
    1. Compute minimizer sketches for all reads
    2. Build an inverted index (minimizer hash -> list of (read, position))
    3. For each read pair sharing minimizers, compute overlap coordinates
    4. Filter by minimum overlap length and maximum overhang

    Args:
        reads: Sequence of read dictionaries (with 'sequence' and 'read_id' keys)
            or plain sequence strings.
        min_overlap: Minimum overlap length in base pairs.
        k: K-mer size for minimizer computation.
        w: Window size for minimizer computation.
        min_minimizer_matches: Minimum shared minimizers for a candidate pair.
        max_overhang: Maximum unaligned overhang at overlap ends.

    Returns:
        List of Overlap objects for all detected overlaps.
    """
    # Normalize reads
    read_data: list[tuple[str, str]] = []  # (read_id, sequence)
    for i, read in enumerate(reads):
        if isinstance(read, str):
            read_data.append((f"read_{i}", read))
        elif isinstance(read, dict):
            rid = read.get("read_id", read.get("read_name", f"read_{i}"))
            seq = read.get("sequence", "")
            if seq:
                read_data.append((str(rid), seq))
        elif hasattr(read, "sequence") and read.sequence:
            rid = getattr(read, "read_id", getattr(read, "read_name", f"read_{i}"))
            read_data.append((str(rid), read.sequence))

    if len(read_data) < 2:
        return []

    logger.info("Computing minimizer sketches for %d reads...", len(read_data))

    # Step 1: Compute minimizer sketches
    sketches: dict[str, list[Minimizer]] = {}
    read_lengths: dict[str, int] = {}

    for rid, seq in read_data:
        sketches[rid] = minimizer_sketch(seq, k=k, w=w)
        read_lengths[rid] = len(seq)

    # Step 2: Build inverted index
    # hash_value -> list of (read_id, position)
    index: dict[int, list[tuple[str, int]]] = defaultdict(list)

    for rid, mins in sketches.items():
        for m in mins:
            index[m.hash_value].append((rid, m.position))

    # Step 3: Find candidate pairs
    # Count shared minimizers per read pair
    pair_matches: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)

    for hash_val, entries in index.items():
        # Skip overly common minimizers (repetitive sequence)
        if len(entries) > max(50, len(read_data) // 2):
            continue

        for i in range(len(entries)):
            for j in range(i + 1, len(entries)):
                rid_i, pos_i = entries[i]
                rid_j, pos_j = entries[j]
                if rid_i == rid_j:
                    continue
                key = (min(rid_i, rid_j), max(rid_i, rid_j))
                if key[0] == rid_i:
                    pair_matches[key].append((pos_i, pos_j))
                else:
                    pair_matches[key].append((pos_j, pos_i))

    # Step 4: Score and filter overlaps
    overlaps: list[Overlap] = []

    for (rid_a, rid_b), matches in pair_matches.items():
        if len(matches) < min_minimizer_matches:
            continue

        len_a = read_lengths[rid_a]
        len_b = read_lengths[rid_b]

        # Compute overlap coordinates by chaining matches
        overlap = _chain_minimizer_matches(
            rid_a, len_a, rid_b, len_b, matches, min_overlap, max_overhang,
        )

        if overlap is not None:
            overlaps.append(overlap)

    logger.info("Found %d overlaps from %d candidate pairs", len(overlaps), len(pair_matches))
    return overlaps


def _chain_minimizer_matches(
    query_name: str,
    query_length: int,
    target_name: str,
    target_length: int,
    matches: list[tuple[int, int]],
    min_overlap: int,
    max_overhang: int,
) -> Overlap | None:
    """Chain shared minimizer positions to compute overlap coordinates.

    Uses the diagonal consistency of minimizer matches to determine
    overlap extent. Matches on the same diagonal (pos_q - pos_t = constant)
    indicate colinear overlap.

    Args:
        query_name: Query read name.
        query_length: Query read length.
        target_name: Target read name.
        target_length: Target read length.
        matches: List of (query_pos, target_pos) minimizer matches.
        min_overlap: Minimum overlap length.
        max_overhang: Maximum overhang.

    Returns:
        Overlap object if valid overlap detected, None otherwise.
    """
    if not matches:
        return None

    # Compute diagonals (query_pos - target_pos)
    diagonals: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for qpos, tpos in matches:
        # Bin diagonals to allow for indels
        diag_bin = (qpos - tpos) // 100 * 100
        diagonals[diag_bin].append((qpos, tpos))

    # Find the best diagonal (most matches)
    best_diag = max(diagonals.keys(), key=lambda d: len(diagonals[d]))
    best_matches = sorted(diagonals[best_diag], key=lambda m: m[0])

    if len(best_matches) < 2:
        return None

    # Compute overlap from chained matches on best diagonal
    q_start = best_matches[0][0]
    q_end = best_matches[-1][0]
    t_start = best_matches[0][1]
    t_end = best_matches[-1][1]

    overlap_length = max(q_end - q_start, t_end - t_start)

    if overlap_length < min_overlap:
        return None

    # Check overhangs
    q_overhang_start = q_start
    q_overhang_end = query_length - q_end
    t_overhang_start = t_start
    t_overhang_end = target_length - t_end

    # A valid overlap should have small overhang on one side for each read
    min_q_overhang = min(q_overhang_start, q_overhang_end)
    min_t_overhang = min(t_overhang_start, t_overhang_end)

    if min_q_overhang > max_overhang and min_t_overhang > max_overhang:
        return None

    # Estimate identity from match density
    expected_minimizers = overlap_length // 50  # Rough estimate
    identity = min(1.0, len(best_matches) / max(expected_minimizers, 1))

    # Check for containment
    is_contained = (q_start <= max_overhang and
                    query_length - q_end <= max_overhang and
                    overlap_length >= query_length * 0.8)

    return Overlap(
        query_name=query_name,
        query_length=query_length,
        query_start=q_start,
        query_end=q_end,
        target_name=target_name,
        target_length=target_length,
        target_start=t_start,
        target_end=t_end,
        strand="+",
        num_matches=len(best_matches),
        overlap_length=overlap_length,
        identity=identity,
        is_contained=is_contained,
    )


def compute_overlap_graph(overlaps: Sequence[Overlap]) -> OverlapGraph:
    """Build a string graph from detected overlaps.

    Creates a directed graph where nodes are reads and edges represent
    dovetail overlaps. Contained reads are identified but not removed.

    Args:
        overlaps: Sequence of Overlap objects.

    Returns:
        OverlapGraph with nodes, edges, and adjacency list.
    """
    nodes_set: set[str] = set()
    edges: list[tuple[str, str, Overlap]] = []
    adjacency: dict[str, list[tuple[str, Overlap]]] = defaultdict(list)

    for ovl in overlaps:
        nodes_set.add(ovl.query_name)
        nodes_set.add(ovl.target_name)

        edges.append((ovl.query_name, ovl.target_name, ovl))
        adjacency[ovl.query_name].append((ovl.target_name, ovl))
        adjacency[ovl.target_name].append((ovl.query_name, ovl))

    nodes = sorted(nodes_set)

    graph = OverlapGraph(
        nodes=nodes,
        edges=edges,
        adjacency=dict(adjacency),
        num_nodes=len(nodes),
        num_edges=len(edges),
    )

    logger.info("Built overlap graph: %d nodes, %d edges", graph.num_nodes, graph.num_edges)
    return graph


def filter_contained_reads(overlaps: Sequence[Overlap]) -> list[Overlap]:
    """Remove overlaps involving contained reads.

    A read is contained if it is fully encompassed by another longer read
    (both ends of the shorter read overlap within the longer read).
    Contained reads are redundant for assembly and should be removed.

    Args:
        overlaps: Sequence of Overlap objects.

    Returns:
        List of overlaps with contained reads removed.
    """
    # Identify contained reads
    contained: set[str] = set()

    for ovl in overlaps:
        if ovl.is_contained:
            contained.add(ovl.query_name)

    if contained:
        logger.info("Identified %d contained reads for removal", len(contained))

    # Filter out overlaps involving contained reads
    filtered = [
        ovl for ovl in overlaps
        if ovl.query_name not in contained and ovl.target_name not in contained
    ]

    logger.info(
        "Filtered overlaps: %d -> %d (removed %d involving contained reads)",
        len(overlaps), len(filtered), len(overlaps) - len(filtered),
    )
    return filtered


# --- Internal helper functions ---


def _hash_kmer(kmer: str) -> int:
    """Compute a consistent integer hash for a k-mer string.

    Uses a custom polynomial rolling hash for speed, with FNV-1a
    as the base algorithm for good distribution.
    """
    # FNV-1a inspired hash
    h = 0x811C9DC5
    for c in kmer:
        h ^= ord(c)
        h = (h * 0x01000193) & 0xFFFFFFFF
    return h


def _reverse_complement(seq: str) -> str:
    """Compute the reverse complement of a DNA sequence."""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(c, "N") for c in reversed(seq))
