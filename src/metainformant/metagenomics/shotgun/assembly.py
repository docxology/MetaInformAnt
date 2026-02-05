"""Metagenome assembly from shotgun sequencing data.

Implements a simplified de Bruijn graph assembler that supports multiple k-mer
sizes, contig scaffolding using paired-end information, and comprehensive
assembly quality statistics (N50, L50, GC content, etc.).

The de Bruijn graph approach:
1. Extract k-mers from reads to build a directed graph where nodes are
   (k-1)-mers and edges are k-mers.
2. Simplify the graph by collapsing non-branching paths into unitigs.
3. Traverse the graph to produce contigs.
4. Optionally merge results from multiple k-mer sizes.
"""

from __future__ import annotations

import os
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "META_"


@dataclass
class Contig:
    """Assembled contig.

    Attributes:
        contig_id: Unique identifier for this contig.
        sequence: Nucleotide sequence of the contig.
        length: Length of the contig in base pairs.
        coverage: Average read coverage depth.
        gc_content: GC content as a fraction (0.0 to 1.0).
    """

    contig_id: str
    sequence: str
    length: int = 0
    coverage: float = 0.0
    gc_content: float = 0.0

    def __post_init__(self) -> None:
        if self.length == 0 and self.sequence:
            self.length = len(self.sequence)
        if self.gc_content == 0.0 and self.sequence:
            self.gc_content = _calculate_gc(self.sequence)


@dataclass
class AssemblyStats:
    """Assembly quality statistics.

    Attributes:
        total_contigs: Number of contigs.
        total_length: Sum of all contig lengths.
        largest_contig: Length of the largest contig.
        smallest_contig: Length of the smallest contig.
        n50: N50 statistic (50% of total assembly is in contigs >= this length).
        l50: L50 (minimum number of contigs whose length sum >= 50% of total).
        n90: N90 statistic.
        l90: L90 statistic.
        gc_content: Overall GC content as fraction.
        mean_coverage: Mean coverage across all contigs.
        mean_length: Average contig length.
        median_length: Median contig length.
    """

    total_contigs: int = 0
    total_length: int = 0
    largest_contig: int = 0
    smallest_contig: int = 0
    n50: int = 0
    l50: int = 0
    n90: int = 0
    l90: int = 0
    gc_content: float = 0.0
    mean_coverage: float = 0.0
    mean_length: float = 0.0
    median_length: float = 0.0


@dataclass
class Scaffold:
    """Scaffolded assembly.

    Attributes:
        scaffold_id: Unique identifier.
        contigs: Ordered list of contigs in this scaffold.
        gaps: Gap sizes between consecutive contigs (in bp, estimated from paired-end data).
        sequence: Full scaffold sequence (with N's for gaps).
        total_length: Total scaffold length including gaps.
    """

    scaffold_id: str
    contigs: list[Contig] = field(default_factory=list)
    gaps: list[int] = field(default_factory=list)
    sequence: str = ""
    total_length: int = 0


def _calculate_gc(sequence: str) -> float:
    """Calculate GC content of a nucleotide sequence.

    Args:
        sequence: DNA sequence string.

    Returns:
        GC fraction (0.0 to 1.0).
    """
    if not sequence:
        return 0.0
    seq_upper = sequence.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    valid_bases = sum(1 for c in seq_upper if c in "ACGT")
    return gc_count / valid_bases if valid_bases > 0 else 0.0


def _build_de_bruijn_graph(reads: list[str], k: int) -> dict[str, list[str]]:
    """Build a de Bruijn graph from reads using k-mers.

    Each k-mer is an edge from its (k-1)-prefix to its (k-1)-suffix.

    Args:
        reads: List of sequencing read strings.
        k: K-mer size.

    Returns:
        Adjacency list: prefix -> [list of suffixes].
    """
    graph: dict[str, list[str]] = defaultdict(list)
    for read in reads:
        seq = read.upper()
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            if any(c not in "ACGT" for c in kmer):
                continue
            prefix = kmer[:-1]
            suffix = kmer[1:]
            graph[prefix].append(suffix)
    return dict(graph)


def _compute_node_degrees(graph: dict[str, list[str]]) -> tuple[dict[str, int], dict[str, int]]:
    """Compute in-degree and out-degree for each node.

    Args:
        graph: Adjacency list representation.

    Returns:
        Tuple of (in_degree dict, out_degree dict).
    """
    in_degree: dict[str, int] = defaultdict(int)
    out_degree: dict[str, int] = defaultdict(int)

    for node, neighbors in graph.items():
        out_degree[node] = len(neighbors)
        for neighbor in neighbors:
            in_degree[neighbor] += 1
            if neighbor not in out_degree:
                out_degree[neighbor] = 0
        if node not in in_degree:
            in_degree[node] = 0

    return dict(in_degree), dict(out_degree)


def _collapse_non_branching_paths(graph: dict[str, list[str]]) -> list[str]:
    """Collapse non-branching paths in the de Bruijn graph into contigs.

    A non-branching path is a sequence of nodes where each has exactly
    one incoming and one outgoing edge. These are collapsed into a single
    contig string.

    Args:
        graph: De Bruijn graph as adjacency list.

    Returns:
        List of contig sequences.
    """
    if not graph:
        return []

    in_deg, out_deg = _compute_node_degrees(graph)
    all_nodes = set(in_deg.keys()) | set(out_deg.keys())

    # Find start nodes: nodes that are not internal to a non-branching path
    # A node is a start node if: in_degree != 1 or out_degree != 1, or it starts a cycle
    start_nodes: set[str] = set()
    for node in all_nodes:
        in_d = in_deg.get(node, 0)
        out_d = out_deg.get(node, 0)
        if in_d != 1 or out_d != 1:
            start_nodes.add(node)

    # If no start nodes found (everything is in cycles), pick arbitrary starts
    if not start_nodes:
        start_nodes = {next(iter(all_nodes))} if all_nodes else set()

    visited_edges: set[tuple[str, str]] = set()
    contigs: list[str] = []

    for start in start_nodes:
        if start not in graph:
            continue
        for neighbor in graph[start]:
            if (start, neighbor) in visited_edges:
                continue

            # Build path
            path = [start, neighbor]
            visited_edges.add((start, neighbor))
            current = neighbor

            while True:
                out_d = out_deg.get(current, 0)
                in_d = in_deg.get(current, 0)
                if out_d != 1 or in_d != 1:
                    break
                if current not in graph:
                    break
                next_node = graph[current][0]
                if (current, next_node) in visited_edges:
                    break
                visited_edges.add((current, next_node))
                path.append(next_node)
                current = next_node

            # Build contig sequence from path of (k-1)-mers
            if len(path) > 1:
                contig_seq = path[0]
                for node in path[1:]:
                    contig_seq += node[-1]
                contigs.append(contig_seq)

    return contigs


def assemble_contigs(
    reads: dict[str, str] | list[str],
    k_range: list[int] | None = None,
    min_contig_length: int = 200,
    min_kmer_coverage: int = 2,
) -> list[Contig]:
    """Assemble reads into contigs using de Bruijn graph approach.

    Builds de Bruijn graphs at multiple k-mer sizes, extracts contigs by
    collapsing non-branching paths, and merges results. Larger k-mers
    produce more specific (but potentially shorter) contigs, while smaller
    k-mers produce longer (but potentially more chimeric) contigs.

    Args:
        reads: Sequencing reads, either as dict (id->sequence) or list of sequences.
        k_range: List of k-mer sizes to use. Default [21, 33, 55].
            Must be odd numbers. Multiple sizes improve assembly by capturing
            different repeat structures.
        min_contig_length: Minimum contig length to report (default 200).
        min_kmer_coverage: Minimum k-mer coverage to include in graph (default 2).
            Helps filter sequencing errors.

    Returns:
        List of Contig objects sorted by length (descending).

    Raises:
        ValueError: If reads is empty or k values are invalid.

    Examples:
        >>> reads = ["ATCGATCGATCGATCGATCG", "CGATCGATCGATCGATCGTT",
        ...          "ATCGATCGATCGATCGATCG"]
        >>> contigs = assemble_contigs(reads, k_range=[11])
        >>> all(c.length >= 11 for c in contigs)
        True
    """
    if not reads:
        raise ValueError("Input reads must not be empty")

    # Normalize input
    if isinstance(reads, dict):
        read_list = list(reads.values())
    else:
        read_list = list(reads)

    if k_range is None:
        # Default k-mer sizes; adjust based on read length
        avg_len = sum(len(r) for r in read_list) / len(read_list)
        if avg_len < 100:
            k_range = [21]
        elif avg_len < 200:
            k_range = [21, 33]
        else:
            k_range = [21, 33, 55]

    for k in k_range:
        if k < 11:
            raise ValueError(f"K-mer size must be >= 11, got {k}")
        if k % 2 == 0:
            logger.warning(f"K-mer size {k} is even; odd values recommended for avoiding palindromes")

    logger.info(f"Assembling {len(read_list)} reads with k-mer sizes {k_range}")

    all_contigs: list[Contig] = []
    contig_counter = 0

    for k in k_range:
        logger.info(f"Building de Bruijn graph with k={k}")

        # Filter k-mers by coverage
        kmer_counts: Counter[str] = Counter()
        for read in read_list:
            seq = read.upper()
            for i in range(len(seq) - k + 1):
                kmer = seq[i : i + k]
                if all(c in "ACGT" for c in kmer):
                    kmer_counts[kmer] += 1

        # Build filtered graph
        graph: dict[str, list[str]] = defaultdict(list)
        edge_coverage: dict[tuple[str, str], int] = {}
        for kmer, count in kmer_counts.items():
            if count >= min_kmer_coverage:
                prefix = kmer[:-1]
                suffix = kmer[1:]
                graph[prefix].append(suffix)
                edge_coverage[(prefix, suffix)] = count

        if not graph:
            logger.warning(f"No k-mers passed coverage filter at k={k}")
            continue

        # Collapse non-branching paths
        contig_seqs = _collapse_non_branching_paths(dict(graph))

        for seq in contig_seqs:
            if len(seq) < min_contig_length:
                continue

            # Estimate coverage from k-mer counts
            kmer_coverages: list[int] = []
            for i in range(len(seq) - k + 1):
                kmer = seq[i : i + k]
                kmer_coverages.append(kmer_counts.get(kmer, 0))
            avg_cov = sum(kmer_coverages) / len(kmer_coverages) if kmer_coverages else 0.0

            contig = Contig(
                contig_id=f"contig_{contig_counter:06d}_k{k}",
                sequence=seq,
                coverage=avg_cov,
            )
            all_contigs.append(contig)
            contig_counter += 1

    # Remove redundant contigs from different k-mer sizes
    # Simple approach: remove contigs that are substrings of longer contigs
    all_contigs.sort(key=lambda c: -c.length)
    filtered_contigs: list[Contig] = []
    seen_sequences: set[str] = set()

    for contig in all_contigs:
        # Check if this contig is a substring of any already-accepted contig
        is_redundant = False
        for accepted_seq in seen_sequences:
            if contig.sequence in accepted_seq:
                is_redundant = True
                break
        if not is_redundant:
            filtered_contigs.append(contig)
            seen_sequences.add(contig.sequence)

    # Re-number contigs
    for idx, contig in enumerate(filtered_contigs):
        contig.contig_id = f"contig_{idx:06d}"

    logger.info(f"Assembly complete: {len(filtered_contigs)} contigs (from {len(all_contigs)} raw)")
    return filtered_contigs


def scaffold_contigs(
    contigs: list[Contig],
    paired_reads: dict[str, tuple[str, str]] | None = None,
    insert_size: int = 500,
    insert_std: int = 100,
    min_links: int = 3,
) -> list[Scaffold]:
    """Scaffold contigs using paired-end read linkage information.

    For each pair of contigs, counts the number of paired-end reads where
    one read maps to one contig and the mate maps to the other. Contig
    pairs with sufficient linking evidence (>= min_links) are joined
    into scaffolds with estimated gap sizes based on insert size.

    Args:
        contigs: List of assembled Contig objects.
        paired_reads: Dictionary mapping pair IDs to (forward, reverse) tuples.
            If None, each contig becomes its own scaffold.
        insert_size: Expected insert size in bp (default 500).
        insert_std: Standard deviation of insert size (default 100).
        min_links: Minimum paired-end links to join contigs (default 3).

    Returns:
        List of Scaffold objects.

    Examples:
        >>> c1 = Contig("c1", "ATCGATCGATCG" * 20)
        >>> c2 = Contig("c2", "TTTTGGGGCCCC" * 20)
        >>> scaffolds = scaffold_contigs([c1, c2])
        >>> len(scaffolds) >= 1
        True
    """
    if not contigs:
        return []

    logger.info(f"Scaffolding {len(contigs)} contigs")

    if paired_reads is None or not paired_reads:
        # No paired-end data: each contig is its own scaffold
        scaffolds = []
        for idx, contig in enumerate(contigs):
            scaffold = Scaffold(
                scaffold_id=f"scaffold_{idx:06d}",
                contigs=[contig],
                gaps=[],
                sequence=contig.sequence,
                total_length=contig.length,
            )
            scaffolds.append(scaffold)
        return scaffolds

    # Map reads to contigs using k-mer anchoring
    k = 31
    contig_kmer_index: dict[str, list[int]] = {}  # kmer -> [contig indices]

    for idx, contig in enumerate(contigs):
        seq = contig.sequence.upper()
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            if all(c in "ACGT" for c in kmer):
                if kmer not in contig_kmer_index:
                    contig_kmer_index[kmer] = []
                contig_kmer_index[kmer].append(idx)

    # Count paired-end links between contigs
    link_counts: dict[tuple[int, int], int] = defaultdict(int)

    for pair_id, (fwd, rev) in paired_reads.items():
        fwd_upper = fwd.upper()
        rev_upper = rev.upper()

        # Find which contigs the forward and reverse reads map to
        fwd_contigs: set[int] = set()
        for i in range(len(fwd_upper) - k + 1):
            kmer = fwd_upper[i : i + k]
            if kmer in contig_kmer_index:
                fwd_contigs.update(contig_kmer_index[kmer])

        rev_contigs: set[int] = set()
        for i in range(len(rev_upper) - k + 1):
            kmer = rev_upper[i : i + k]
            if kmer in contig_kmer_index:
                rev_contigs.update(contig_kmer_index[kmer])

        # Count links between different contigs
        for fc in fwd_contigs:
            for rc in rev_contigs:
                if fc != rc:
                    pair_key = (min(fc, rc), max(fc, rc))
                    link_counts[pair_key] += 1

    # Build scaffold graph
    scaffold_edges: list[tuple[int, int, int]] = []  # (contig_a, contig_b, link_count)
    for (a, b), count in link_counts.items():
        if count >= min_links:
            scaffold_edges.append((a, b, count))

    # Sort edges by link count (greedy scaffolding)
    scaffold_edges.sort(key=lambda x: -x[2])

    # Build scaffolds using union-find
    parent: dict[int, int] = {i: i for i in range(len(contigs))}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    scaffold_members: dict[int, list[int]] = {i: [i] for i in range(len(contigs))}

    for a, b, count in scaffold_edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            union(ra, rb)
            new_root = find(ra)
            # Merge member lists
            if new_root == find(ra):
                members_a = scaffold_members.pop(ra, [ra])
                members_b = scaffold_members.pop(rb, [rb])
                scaffold_members[new_root] = members_a + members_b
            else:
                members_a = scaffold_members.pop(ra, [ra])
                members_b = scaffold_members.pop(rb, [rb])
                scaffold_members[new_root] = members_b + members_a

    # Rebuild scaffold member sets based on final parent
    final_groups: dict[int, list[int]] = defaultdict(list)
    for i in range(len(contigs)):
        final_groups[find(i)].append(i)

    # Create scaffold objects
    scaffolds: list[Scaffold] = []
    for scaffold_idx, (root, members) in enumerate(sorted(final_groups.items())):
        # Sort members by contig length (largest first)
        members.sort(key=lambda i: -contigs[i].length)
        scaffold_contigs_list = [contigs[i] for i in members]

        # Estimate gaps between consecutive contigs
        gaps: list[int] = []
        for i in range(len(scaffold_contigs_list) - 1):
            # Estimate gap from insert size and contig lengths
            estimated_gap = max(10, insert_size - min(scaffold_contigs_list[i].length, insert_size))
            gaps.append(estimated_gap)

        # Build scaffold sequence with N-gaps
        scaffold_seq_parts: list[str] = []
        for i, contig in enumerate(scaffold_contigs_list):
            scaffold_seq_parts.append(contig.sequence)
            if i < len(gaps):
                scaffold_seq_parts.append("N" * gaps[i])
        scaffold_seq = "".join(scaffold_seq_parts)

        scaffold = Scaffold(
            scaffold_id=f"scaffold_{scaffold_idx:06d}",
            contigs=scaffold_contigs_list,
            gaps=gaps,
            sequence=scaffold_seq,
            total_length=len(scaffold_seq),
        )
        scaffolds.append(scaffold)

    scaffolds.sort(key=lambda s: -s.total_length)
    logger.info(f"Scaffolding complete: {len(scaffolds)} scaffolds from {len(contigs)} contigs")
    return scaffolds


def calculate_assembly_stats(contigs: list[Contig]) -> AssemblyStats:
    """Calculate comprehensive assembly quality statistics.

    Computes standard assembly metrics including N50, L50, N90, L90,
    total length, GC content, and coverage statistics.

    N50 is defined as: the minimum contig length such that 50% of the
    total assembly length is contained in contigs of this length or longer.

    Args:
        contigs: List of assembled Contig objects.

    Returns:
        AssemblyStats with all computed metrics.

    Raises:
        ValueError: If contigs list is empty.

    Examples:
        >>> contigs = [Contig("c1", "A" * 1000), Contig("c2", "T" * 500)]
        >>> stats = calculate_assembly_stats(contigs)
        >>> stats.total_contigs
        2
        >>> stats.n50
        1000
    """
    if not contigs:
        raise ValueError("Contigs list must not be empty")

    lengths = sorted([c.length for c in contigs], reverse=True)
    total_length = sum(lengths)
    coverages = [c.coverage for c in contigs if c.coverage > 0]

    # N50 and L50
    n50, l50 = _compute_nx(lengths, total_length, 0.50)
    # N90 and L90
    n90, l90 = _compute_nx(lengths, total_length, 0.90)

    # GC content (weighted by contig length)
    total_gc = 0
    total_valid = 0
    for contig in contigs:
        seq_upper = contig.sequence.upper()
        total_gc += seq_upper.count("G") + seq_upper.count("C")
        total_valid += sum(1 for c in seq_upper if c in "ACGT")
    overall_gc = total_gc / total_valid if total_valid > 0 else 0.0

    # Median length
    n = len(lengths)
    if n % 2 == 0:
        median_length = (lengths[n // 2 - 1] + lengths[n // 2]) / 2.0
    else:
        median_length = float(lengths[n // 2])

    stats = AssemblyStats(
        total_contigs=len(contigs),
        total_length=total_length,
        largest_contig=lengths[0],
        smallest_contig=lengths[-1],
        n50=n50,
        l50=l50,
        n90=n90,
        l90=l90,
        gc_content=overall_gc,
        mean_coverage=sum(coverages) / len(coverages) if coverages else 0.0,
        mean_length=total_length / len(contigs),
        median_length=median_length,
    )

    logger.info(
        f"Assembly stats: {stats.total_contigs} contigs, {stats.total_length} bp, "
        f"N50={stats.n50}, GC={stats.gc_content:.2%}"
    )
    return stats


def _compute_nx(lengths: list[int], total_length: int, fraction: float) -> tuple[int, int]:
    """Compute Nx and Lx statistics.

    Args:
        lengths: Sorted (descending) list of contig lengths.
        total_length: Total assembly length.
        fraction: Fraction threshold (e.g., 0.5 for N50).

    Returns:
        Tuple of (Nx value, Lx count).
    """
    target = total_length * fraction
    cumulative = 0
    for idx, length in enumerate(lengths):
        cumulative += length
        if cumulative >= target:
            return length, idx + 1
    return lengths[-1] if lengths else 0, len(lengths)


__all__ = [
    "AssemblyStats",
    "Contig",
    "Scaffold",
    "assemble_contigs",
    "calculate_assembly_stats",
    "scaffold_contigs",
]
