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

import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from metainformant.core.sequence import COMPLEMENT_TABLE_UPPER, reverse_complement
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "META_"

# K-mer size used to anchor paired reads on contigs
_ANCHOR_K = 31

# Smallest gap written into a scaffold. It is also the gap of a junction whose pairs give no
# usable position (see ``scaffold_contigs``), where the size is genuinely unknown.
_MIN_GAP = 10


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
        sequence: Full scaffold sequence (with N's for gaps); a contig whose orientation
            is "-" appears reverse-complemented.
        total_length: Total scaffold length including gaps.
        orientations: One entry per contig: "+" if the contig is used as given, "-" if its
            reverse complement is used.
        link_counts: One entry per junction: the number of read pairs that link the two
            contigs in the way the scaffold joins them.
        orientation_support: One entry per junction: how many of those pairs have both reads
            placed on a definite strand of their contig. 0 means no pair could resolve the
            relative orientation (for example contigs that are their own reverse complement),
            so the junction keeps the orientation the contigs were given in.
    """

    scaffold_id: str
    contigs: list[Contig] = field(default_factory=list)
    gaps: list[int] = field(default_factory=list)
    sequence: str = ""
    total_length: int = 0
    orientations: list[str] = field(default_factory=list)
    link_counts: list[int] = field(default_factory=list)
    orientation_support: list[int] = field(default_factory=list)


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

    # If no start nodes found (everything is in cycles), pick a start (the smallest node,
    # so the result does not depend on set iteration order)
    if not start_nodes:
        start_nodes = {min(all_nodes)} if all_nodes else set()

    visited_edges: set[tuple[str, str]] = set()
    contigs: list[str] = []

    # Sorted so that contig order (and hence contig IDs) is reproducible: iterating a set of
    # strings follows the per-process hash seed.
    for start in sorted(start_nodes):
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


@dataclass(frozen=True)
class _ReadPlacement:
    """Where a read lies on a contig.

    Attributes:
        contig: Index of the contig.
        strand: "+" when the read matches the contig as given, "-" when it matches the
            reverse complement, None when both strands match equally well.
        start: Index of the read's first base on the matched strand (the "+" strand when
            ``strand`` is None); negative or beyond the contig when the read overhangs an end.
    """

    contig: int
    strand: str | None
    start: int


def _build_anchor_index(contigs: list[Contig], k: int) -> dict[str, list[tuple[int, int, int]]]:
    """Index every k-mer of both strands of every contig.

    Both strands are indexed because a read (typically the reverse mate) can come from either
    strand. For each k-mer only its first position on each strand of each contig is kept.

    Args:
        contigs: Contigs to index.
        k: K-mer size.

    Returns:
        Mapping of k-mer to ``(contig index, strand, position)`` where strand is 0 for the
        contig as given and 1 for its reverse complement.
    """
    index: dict[str, list[tuple[int, int, int]]] = {}
    for idx, contig in enumerate(contigs):
        seq = contig.sequence.upper()
        for strand, strand_seq in enumerate((seq, seq.translate(COMPLEMENT_TABLE_UPPER)[::-1])):
            for pos in range(len(strand_seq) - k + 1):
                kmer = strand_seq[pos : pos + k]
                if not all(c in "ACGT" for c in kmer):
                    continue
                entries = index.setdefault(kmer, [])
                if not any(e[0] == idx and e[1] == strand for e in entries):
                    entries.append((idx, strand, pos))
    return index


def _place_read(read: str, index: dict[str, list[tuple[int, int, int]]], k: int) -> _ReadPlacement | None:
    """Place a read on the contig most of its k-mers hit, and on the strand they favour.

    Args:
        read: Upper-case read sequence.
        index: Output of ``_build_anchor_index``.
        k: K-mer size of the index.

    Returns:
        The placement, or None if no k-mer hits a contig or two contigs hit equally well (a
        repeat shared between contigs anchors the read nowhere in particular).
    """
    votes: dict[tuple[int, int], int] = defaultdict(int)
    starts: dict[tuple[int, int], list[int]] = defaultdict(list)
    for i in range(len(read) - k + 1):
        for contig, strand, pos in index.get(read[i : i + k], ()):
            votes[(contig, strand)] += 1
            starts[(contig, strand)].append(pos - i)
    if not votes:
        return None

    best_per_contig: dict[int, int] = defaultdict(int)
    for (contig, _), count in votes.items():
        best_per_contig[contig] = max(best_per_contig[contig], count)
    top = max(best_per_contig.values())
    winners = [contig for contig, count in best_per_contig.items() if count == top]
    if len(winners) != 1:
        return None

    contig = winners[0]
    forward, reverse = votes.get((contig, 0), 0), votes.get((contig, 1), 0)
    if forward > reverse:
        return _ReadPlacement(contig, "+", int(statistics.median(starts[(contig, 0)])))
    if reverse > forward:
        return _ReadPlacement(contig, "-", int(statistics.median(starts[(contig, 1)])))
    return _ReadPlacement(contig, None, int(statistics.median(starts[(contig, 0)])))


def scaffold_contigs(
    contigs: list[Contig],
    paired_reads: dict[str, tuple[str, str]] | None = None,
    insert_size: int = 500,
    insert_std: int = 100,
    min_links: int = 3,
) -> list[Scaffold]:
    """Scaffold contigs using paired-end read linkage information.

    Each read is anchored on a contig by its 31-mers, on either strand of the contig (a mate
    sequenced from the opposite strand matches the reverse complement). A read is used only if
    one contig attracts clearly more of its k-mers than any other.

    The library is assumed to be a standard paired-end one: the forward read lies on one
    strand of the fragment and the reverse read on the other, facing each other (FR, as
    Illumina paired-end libraries are). The strand each mate matches then fixes which end of
    its contig faces the other contig, so a pair says how two contigs lie relative to each
    other and in which orientation (a mate-pair or long-insert library, where the reads face
    away from each other, is not supported). The relative orientation of every join is
    inferred from these strands and the scaffold sequence contains the reverse complement of
    contigs that have to be flipped.

    Joins are chosen greedily: pairs vote for a join of two contig ends, joins with at least
    ``min_links`` pairs are taken in order of the support for their orientation, each contig
    end is joined at most once and no join closes a cycle, so every scaffold is a linear
    chain of oriented contigs. A contig linked to several others keeps only its best join
    at each end. The position of each mate gives the distance from its read to the contig
    end, so a junction's gap is ``insert_size`` minus those two distances, taken as the median
    over the pairs and raised to at least 10 bp. That estimate is biased low as the gap
    approaches the insert size, because only the longer fragments of the insert distribution can
    span a wide gap. A pair whose gap would be shorter than ``-3 * insert_std`` is inconsistent
    with the two contigs being neighbours and is ignored.

    Limitations: a read whose two strands match equally well (a contig that is its own reverse
    complement) cannot resolve the orientation. Such pairs still link the contigs, but they
    support every orientation equally, so the contigs are joined in the orientation given and
    ``orientation_support`` is 0 for that junction; its gap is the 10 bp minimum because no
    position is available. Repeats are not modelled beyond ignoring reads that two contigs
    claim equally.

    Args:
        contigs: List of assembled Contig objects.
        paired_reads: Dictionary mapping pair IDs to (forward, reverse) tuples.
            If None, each contig becomes its own scaffold.
        insert_size: Expected fragment length in bp, from the outer end of one read to the
            outer end of its mate (default 500).
        insert_std: Standard deviation of insert size (default 100).
        min_links: Minimum paired-end links to join contigs (default 3).

    Returns:
        List of Scaffold objects, longest first.

    Examples:
        >>> c1 = Contig("c1", "ATCGATCGATCG" * 20)
        >>> c2 = Contig("c2", "TTTTGGGGCCCC" * 20)
        >>> example_scaffolds = scaffold_contigs([c1, c2])
        >>> len(example_scaffolds) >= 1
        True
    """
    if not contigs:
        return []

    logger.info(f"Scaffolding {len(contigs)} contigs")

    if paired_reads is None or not paired_reads:
        # No paired-end data: each contig is its own scaffold
        early_scaffolds: list[Scaffold] = []
        for idx, contig in enumerate(contigs):
            scaffold = Scaffold(
                scaffold_id=f"scaffold_{idx:06d}",
                contigs=[contig],
                gaps=[],
                sequence=contig.sequence,
                total_length=contig.length,
                orientations=["+"],
            )
            early_scaffolds.append(scaffold)
        return early_scaffolds

    k = _ANCHOR_K
    index = _build_anchor_index(contigs, k)

    # Vote for joins of two contig ends. An end is (contig index, "L" or "R"): the left end is
    # the start of the contig as given, the right end its last base. A read matching the
    # contig as given faces its right end, one matching the reverse complement its left end.
    End = tuple[int, str]
    Edge = tuple[End, End]
    edge_votes: dict[Edge, int] = defaultdict(int)
    edge_resolved: dict[Edge, int] = defaultdict(int)
    edge_gaps: dict[Edge, list[int]] = defaultdict(list)

    def facing_ends(placement: _ReadPlacement) -> list[End]:
        strands = ("+", "-") if placement.strand is None else (placement.strand,)
        return [(placement.contig, "R" if strand == "+" else "L") for strand in strands]

    for pair_id, (fwd, rev) in paired_reads.items():
        first = _place_read(fwd.upper(), index, k)
        second = _place_read(rev.upper(), index, k)
        if first is None or second is None or first.contig == second.contig:
            continue

        resolved = first.strand is not None and second.strand is not None
        gap_estimate = 0
        if resolved:
            # bases of each contig covered by the fragment, from its read to the facing end
            length_first = len(contigs[first.contig].sequence)
            length_second = len(contigs[second.contig].sequence)
            covered_first = min(max(length_first - first.start, 0), length_first)
            covered_second = min(max(length_second - second.start, 0), length_second)
            gap_estimate = insert_size - covered_first - covered_second
            if gap_estimate < -3 * insert_std:
                continue

        for end_a in facing_ends(first):
            for end_b in facing_ends(second):
                edge: Edge = (end_a, end_b) if end_a < end_b else (end_b, end_a)
                edge_votes[edge] += 1
                if resolved:
                    edge_resolved[edge] += 1
                    edge_gaps[edge].append(gap_estimate)

    candidates = [edge for edge, votes in edge_votes.items() if votes >= min_links]
    # Best-supported orientation first; among equals prefer "first contig then second, both as given"
    candidates.sort(key=lambda e: (-edge_resolved[e], -edge_votes[e], 0 if (e[0][1], e[1][1]) == ("R", "L") else 1, e))

    parent: dict[int, int] = {i: i for i in range(len(contigs))}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    joined: dict[End, tuple[End, Edge]] = {}
    for edge in candidates:
        end_a, end_b = edge
        if end_a in joined or end_b in joined:
            continue
        root_a, root_b = find(end_a[0]), find(end_b[0])
        if root_a == root_b:
            continue
        parent[root_a] = root_b
        joined[end_a] = (end_b, edge)
        joined[end_b] = (end_a, edge)

    components: dict[int, list[int]] = defaultdict(list)
    for i in range(len(contigs)):
        components[find(i)].append(i)

    scaffolds: list[Scaffold] = []
    for scaffold_idx, members in enumerate(sorted(components.values(), key=min)):
        # A chain has two contigs with a free end (one contig with two free ends when alone)
        degree = {m: sum(1 for end in ("L", "R") if (m, end) in joined) for m in members}
        first_contig = min(m for m in members if degree[m] < 2)
        # Read the chain from a free end: orient the first contig so its free end is on the left
        orientation = "-" if (first_contig, "L") in joined else "+"

        chain: list[tuple[int, str]] = []
        junctions: list[Edge] = []
        current = first_contig
        while True:
            chain.append((current, orientation))
            exit_end: End = (current, "R" if orientation == "+" else "L")
            if exit_end not in joined:
                break
            (next_contig, entry_side), edge = joined[exit_end]
            junctions.append(edge)
            orientation = "+" if entry_side == "L" else "-"
            current = next_contig

        gaps: list[int] = []
        for edge in junctions:
            estimates = edge_gaps.get(edge)
            gaps.append(max(_MIN_GAP, round(statistics.median(estimates))) if estimates else _MIN_GAP)

        parts: list[str] = []
        for position, (member, member_orientation) in enumerate(chain):
            sequence = contigs[member].sequence
            parts.append(sequence if member_orientation == "+" else reverse_complement(sequence))
            if position < len(gaps):
                parts.append("N" * gaps[position])
        scaffold_seq = "".join(parts)

        scaffolds.append(
            Scaffold(
                scaffold_id=f"scaffold_{scaffold_idx:06d}",
                contigs=[contigs[member] for member, _ in chain],
                gaps=gaps,
                sequence=scaffold_seq,
                total_length=len(scaffold_seq),
                orientations=[member_orientation for _, member_orientation in chain],
                link_counts=[edge_votes[edge] for edge in junctions],
                orientation_support=[edge_resolved[edge] for edge in junctions],
            )
        )

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
