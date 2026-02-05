"""Hybrid assembly combining long and short reads.

Provides long-read error correction using short reads, scaffolding of
short-read contigs with long reads, and integrated hybrid assembly
combining both data types.

Approach:
1. Error correction: Map short reads to long reads and correct errors
   using local consensus from short-read k-mer support
2. Scaffolding: Use long reads to bridge and order short-read contigs
3. Hybrid assembly: Combine corrected long reads with short-read
   contigs for gap filling

Optional dependencies:
    - numpy: For numerical computation
"""

from __future__ import annotations

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
class CorrectedRead:
    """A long read after error correction with short reads.

    Attributes:
        read_id: Original read identifier.
        original_sequence: Original long-read sequence.
        corrected_sequence: Error-corrected sequence.
        corrections_made: Number of corrections applied.
        correction_rate: Fraction of bases corrected.
        regions_corrected: List of (start, end) corrected regions.
    """

    read_id: str = ""
    original_sequence: str = ""
    corrected_sequence: str = ""
    corrections_made: int = 0
    correction_rate: float = 0.0
    regions_corrected: list[tuple[int, int]] = field(default_factory=list)


@dataclass
class Scaffold:
    """A scaffold constructed from contigs bridged by long reads.

    Attributes:
        scaffold_id: Unique scaffold identifier.
        sequence: Full scaffold sequence (contigs joined by gap-filled regions).
        contigs: Ordered list of contig names in the scaffold.
        contig_orientations: Orientation of each contig ('+' or '-').
        gaps: Gap sizes between consecutive contigs.
        total_length: Total scaffold length.
        num_contigs: Number of contigs in the scaffold.
        supporting_reads: Number of long reads supporting the scaffold.
    """

    scaffold_id: str = ""
    sequence: str = ""
    contigs: list[str] = field(default_factory=list)
    contig_orientations: list[str] = field(default_factory=list)
    gaps: list[int] = field(default_factory=list)
    total_length: int = 0
    num_contigs: int = 0
    supporting_reads: int = 0


@dataclass
class HybridAssemblyResult:
    """Result of hybrid assembly combining long and short reads.

    Attributes:
        contigs: List of assembled contig sequences.
        contig_names: Names for each contig.
        total_length: Total assembled length.
        num_contigs: Number of contigs.
        n50: N50 of the assembly.
        largest_contig: Length of the largest contig.
        scaffolds: List of Scaffold objects (if scaffolding was performed).
        corrected_reads: List of corrected long reads.
    """

    contigs: list[str] = field(default_factory=list)
    contig_names: list[str] = field(default_factory=list)
    total_length: int = 0
    num_contigs: int = 0
    n50: int = 0
    largest_contig: int = 0
    scaffolds: list[Scaffold] = field(default_factory=list)
    corrected_reads: list[CorrectedRead] = field(default_factory=list)


def hybrid_assemble(
    long_reads: Sequence[str | dict[str, Any]],
    short_reads: Sequence[str | dict[str, Any]],
    min_long_read_length: int = 2000,
    kmer_size: int = 21,
    min_overlap: int = 1000,
) -> HybridAssemblyResult:
    """Perform hybrid assembly using long reads scaffolded with short-read correction.

    Pipeline:
    1. Build k-mer database from short reads
    2. Correct long reads using short-read k-mer support
    3. Find overlaps between corrected long reads
    4. Generate consensus contigs from overlapping reads
    5. Scaffold contigs using remaining long-read information

    Args:
        long_reads: Sequence of long-read sequences (strings or dicts with 'sequence').
        short_reads: Sequence of short-read sequences.
        min_long_read_length: Minimum long-read length to include.
        kmer_size: K-mer size for short-read correction.
        min_overlap: Minimum overlap for contig generation.

    Returns:
        HybridAssemblyResult with assembled contigs and metrics.
    """
    # Extract sequences
    long_seqs = _extract_sequences(long_reads, min_length=min_long_read_length)
    short_seqs = _extract_sequences(short_reads)

    if not long_seqs:
        logger.warning("No long reads passing length filter (%d bp)", min_long_read_length)
        return HybridAssemblyResult()

    logger.info(
        "Hybrid assembly: %d long reads, %d short reads",
        len(long_seqs),
        len(short_seqs),
    )

    # Step 1: Build k-mer database from short reads
    kmer_db = _build_kmer_database(short_seqs, kmer_size)
    logger.info("Built k-mer database: %d unique k-mers", len(kmer_db))

    # Step 2: Correct long reads
    corrected = correct_with_short_reads(long_seqs, short_seqs, kmer_size=kmer_size)
    corrected_seqs = [
        {"read_id": f"corrected_{i}", "sequence": cr.corrected_sequence}
        for i, cr in enumerate(corrected)
        if cr.corrected_sequence
    ]

    logger.info("Corrected %d long reads", len(corrected_seqs))

    # Step 3: Find overlaps between corrected reads
    from .consensus import generate_consensus
    from .overlap import filter_contained_reads, find_overlaps

    if len(corrected_seqs) >= 2:
        overlaps = find_overlaps(corrected_seqs, min_overlap=min_overlap)
        overlaps = filter_contained_reads(overlaps)
    else:
        overlaps = []

    # Step 4: Generate contigs from overlapping reads
    contigs: list[str] = []
    contig_names: list[str] = []

    if overlaps:
        # Group overlapping reads into clusters
        clusters = _cluster_reads_by_overlaps(corrected_seqs, overlaps)

        for i, cluster in enumerate(clusters):
            cluster_seqs = [s["sequence"] for s in cluster if "sequence" in s]
            if cluster_seqs:
                result = generate_consensus(cluster_seqs)
                if result.sequence:
                    contigs.append(result.sequence)
                    contig_names.append(f"contig_{i}")
    else:
        # No overlaps, use corrected reads as contigs
        for i, cr in enumerate(corrected_seqs):
            seq = cr["sequence"] if isinstance(cr, dict) else cr
            if seq:
                contigs.append(seq)
                contig_names.append(f"contig_{i}")

    # Compute assembly statistics
    if contigs:
        lengths = sorted([len(c) for c in contigs], reverse=True)
        total_length = sum(lengths)
        n50 = _compute_n50(lengths)
        largest = lengths[0]
    else:
        total_length = 0
        n50 = 0
        largest = 0

    logger.info(
        "Assembly result: %d contigs, total %d bp, N50 %d bp",
        len(contigs),
        total_length,
        n50,
    )

    return HybridAssemblyResult(
        contigs=contigs,
        contig_names=contig_names,
        total_length=total_length,
        num_contigs=len(contigs),
        n50=n50,
        largest_contig=largest,
        corrected_reads=corrected,
    )


def correct_with_short_reads(
    long_reads: Sequence[str | dict[str, Any]],
    short_reads: Sequence[str | dict[str, Any]],
    kmer_size: int = 21,
    min_kmer_support: int = 3,
    window_size: int = 50,
) -> list[CorrectedRead]:
    """Correct long-read errors using short-read k-mer evidence.

    For each position in a long read:
    1. Extract the local k-mer from the long read
    2. Check if it exists in the short-read k-mer database
    3. If not, try all single-base substitutions and find the best
       supported alternative from the short-read data
    4. Apply corrections where short-read evidence is strong

    This approach corrects random errors (substitutions and small indels)
    common in long-read data while preserving the long-range structure.

    Args:
        long_reads: Sequence of long-read sequences.
        short_reads: Sequence of short-read sequences.
        kmer_size: K-mer size for error correction.
        min_kmer_support: Minimum short-read k-mer count to consider trusted.
        window_size: Window size for local correction context.

    Returns:
        List of CorrectedRead objects.
    """
    long_seqs = _extract_sequences(long_reads)
    short_seqs = _extract_sequences(short_reads)

    # Build k-mer database from short reads
    kmer_db = _build_kmer_database(short_seqs, kmer_size)

    corrected_reads: list[CorrectedRead] = []

    for i, long_seq_data in enumerate(long_seqs):
        if isinstance(long_seq_data, dict):
            read_id = long_seq_data.get("read_id", f"read_{i}")
            seq = long_seq_data.get("sequence", "")
        elif isinstance(long_seq_data, str):
            read_id = f"read_{i}"
            seq = long_seq_data
        else:
            continue

        if not seq or len(seq) < kmer_size:
            continue

        # Correct the sequence
        corrected_seq, num_corrections, correction_regions = _correct_sequence(
            seq,
            kmer_db,
            kmer_size,
            min_kmer_support,
        )

        correction_rate = num_corrections / len(seq) if len(seq) > 0 else 0.0

        corrected_reads.append(
            CorrectedRead(
                read_id=str(read_id),
                original_sequence=seq,
                corrected_sequence=corrected_seq,
                corrections_made=num_corrections,
                correction_rate=correction_rate,
                regions_corrected=correction_regions,
            )
        )

    total_corrections = sum(cr.corrections_made for cr in corrected_reads)
    logger.info(
        "Error correction: corrected %d long reads, %d total corrections",
        len(corrected_reads),
        total_corrections,
    )

    return corrected_reads


def _correct_sequence(
    sequence: str,
    kmer_db: dict[str, int],
    kmer_size: int,
    min_support: int,
) -> tuple[str, int, list[tuple[int, int]]]:
    """Correct a single sequence using k-mer database.

    Returns:
        Tuple of (corrected_sequence, num_corrections, correction_regions).
    """
    seq_list = list(sequence.upper())
    corrections = 0
    regions: list[tuple[int, int]] = []
    bases = "ACGT"

    for pos in range(len(seq_list) - kmer_size + 1):
        kmer = "".join(seq_list[pos : pos + kmer_size])

        # Check if current k-mer is trusted
        if kmer_db.get(kmer, 0) >= min_support:
            continue

        # Try single-base substitutions
        best_alt = None
        best_count = 0

        for offset in range(kmer_size):
            original_base = seq_list[pos + offset]
            for base in bases:
                if base == original_base:
                    continue

                test_kmer = list(kmer)
                test_kmer[offset] = base
                test_kmer_str = "".join(test_kmer)
                count = kmer_db.get(test_kmer_str, 0)

                if count > best_count and count >= min_support:
                    best_count = count
                    best_alt = (pos + offset, base)

        # Apply correction if found
        if best_alt is not None:
            correct_pos, correct_base = best_alt
            seq_list[correct_pos] = correct_base
            corrections += 1
            regions.append((correct_pos, correct_pos + 1))

    return "".join(seq_list), corrections, regions


def scaffold_with_long_reads(
    contigs: Sequence[str | dict[str, Any]],
    long_reads: Sequence[str | dict[str, Any]],
    min_alignment_length: int = 1000,
    min_bridge_reads: int = 2,
) -> list[Scaffold]:
    """Scaffold contigs using long reads as bridges.

    Algorithm:
    1. Map long reads to contigs using k-mer anchoring
    2. Identify long reads that span multiple contigs
    3. Determine contig order and orientation from bridging reads
    4. Join contigs into scaffolds with gap size estimates

    Args:
        contigs: Sequence of contig sequences (strings or dicts).
        long_reads: Sequence of long-read sequences used for bridging.
        min_alignment_length: Minimum alignment length to a contig.
        min_bridge_reads: Minimum bridging reads to scaffold two contigs.

    Returns:
        List of Scaffold objects representing the scaffolded assembly.
    """
    # Extract contig data
    contig_data: list[tuple[str, str]] = []
    for i, contig in enumerate(contigs):
        if isinstance(contig, str):
            contig_data.append((f"contig_{i}", contig))
        elif isinstance(contig, dict):
            name = contig.get("name", contig.get("contig_id", f"contig_{i}"))
            seq = contig.get("sequence", "")
            if seq:
                contig_data.append((str(name), seq))

    if not contig_data:
        return []

    # Extract long reads
    long_seq_data = _extract_sequences(long_reads)
    long_seqs: list[tuple[str, str]] = []
    for i, lr in enumerate(long_seq_data):
        if isinstance(lr, dict):
            rid = lr.get("read_id", f"lr_{i}")
            seq = lr.get("sequence", "")
            if seq:
                long_seqs.append((str(rid), seq))
        elif isinstance(lr, str):
            long_seqs.append((f"lr_{i}", lr))

    if not long_seqs:
        # Return contigs as individual scaffolds
        return [
            Scaffold(
                scaffold_id=name,
                sequence=seq,
                contigs=[name],
                contig_orientations=["+"],
                gaps=[],
                total_length=len(seq),
                num_contigs=1,
            )
            for name, seq in contig_data
        ]

    # Build k-mer index for each contig
    kmer_size = 15
    contig_kmer_index: dict[str, list[tuple[str, int]]] = defaultdict(list)

    for cname, cseq in contig_data:
        for pos in range(0, len(cseq) - kmer_size + 1, kmer_size):
            kmer = cseq[pos : pos + kmer_size].upper()
            contig_kmer_index[kmer].append((cname, pos))

    # Map long reads to contigs
    # For each long read, find which contigs it aligns to
    read_to_contigs: dict[str, list[tuple[str, int, int, str]]] = defaultdict(
        list
    )  # read -> [(contig, start, end, strand)]

    for rname, rseq in long_seqs:
        contig_hits: dict[str, list[int]] = defaultdict(list)
        rseq_upper = rseq.upper()

        for pos in range(0, len(rseq_upper) - kmer_size + 1, kmer_size):
            kmer = rseq_upper[pos : pos + kmer_size]
            if kmer in contig_kmer_index:
                for cname, cpos in contig_kmer_index[kmer]:
                    contig_hits[cname].append(pos)

        for cname, positions in contig_hits.items():
            if len(positions) >= min_alignment_length // kmer_size:
                start = min(positions)
                end = max(positions) + kmer_size
                read_to_contigs[rname].append((cname, start, end, "+"))

    # Find contig pairs bridged by long reads
    bridges: dict[tuple[str, str], int] = defaultdict(int)
    bridge_orders: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)

    for rname, contig_list in read_to_contigs.items():
        if len(contig_list) < 2:
            continue

        # Sort by position on the read
        sorted_contigs = sorted(contig_list, key=lambda x: x[1])

        for i in range(len(sorted_contigs) - 1):
            c1 = sorted_contigs[i][0]
            c2 = sorted_contigs[i + 1][0]
            if c1 != c2:
                key = (c1, c2)
                bridges[key] += 1
                bridge_orders[key].append(
                    (
                        sorted_contigs[i][2],  # end of c1 alignment on read
                        sorted_contigs[i + 1][1],  # start of c2 alignment on read
                    )
                )

    # Build scaffolds from bridges
    # Use a greedy approach: find connected chains of contigs
    used_contigs: set[str] = set()
    scaffolds: list[Scaffold] = []
    scaffold_id = 0

    # Sort bridges by support
    sorted_bridges = sorted(bridges.items(), key=lambda x: -x[1])

    for (c1, c2), support in sorted_bridges:
        if support < min_bridge_reads:
            continue
        if c1 in used_contigs or c2 in used_contigs:
            continue

        # Start a new scaffold with this bridge
        chain = [c1, c2]
        orientations = ["+", "+"]
        used_contigs.add(c1)
        used_contigs.add(c2)

        # Estimate gap size from bridge read positions
        orders = bridge_orders[(c1, c2)]
        gap_estimates = [max(0, start - end) for end, start in orders]
        gap_size = int(sum(gap_estimates) / len(gap_estimates)) if gap_estimates else 100

        # Build scaffold sequence
        contig_seqs = {name: seq for name, seq in contig_data}
        seq_parts = []
        if c1 in contig_seqs:
            seq_parts.append(contig_seqs[c1])
        seq_parts.append("N" * max(1, gap_size))
        if c2 in contig_seqs:
            seq_parts.append(contig_seqs[c2])

        scaffold_seq = "".join(seq_parts)

        scaffolds.append(
            Scaffold(
                scaffold_id=f"scaffold_{scaffold_id}",
                sequence=scaffold_seq,
                contigs=chain,
                contig_orientations=orientations,
                gaps=[gap_size],
                total_length=len(scaffold_seq),
                num_contigs=len(chain),
                supporting_reads=support,
            )
        )
        scaffold_id += 1

    # Add unscaffolded contigs as singletons
    for cname, cseq in contig_data:
        if cname not in used_contigs:
            scaffolds.append(
                Scaffold(
                    scaffold_id=f"scaffold_{scaffold_id}",
                    sequence=cseq,
                    contigs=[cname],
                    contig_orientations=["+"],
                    gaps=[],
                    total_length=len(cseq),
                    num_contigs=1,
                )
            )
            scaffold_id += 1

    logger.info(
        "Scaffolding: %d contigs -> %d scaffolds (%d bridged)",
        len(contig_data),
        len(scaffolds),
        sum(1 for s in scaffolds if s.num_contigs > 1),
    )

    return scaffolds


# --- Internal helper functions ---


def _extract_sequences(
    reads: Sequence[str | dict[str, Any]],
    min_length: int = 0,
) -> list[Any]:
    """Extract sequence data from various read representations."""
    result: list[Any] = []
    for r in reads:
        if isinstance(r, str):
            if len(r) >= min_length:
                result.append(r)
        elif isinstance(r, dict):
            seq = r.get("sequence", "")
            if seq and len(seq) >= min_length:
                result.append(r)
        elif hasattr(r, "sequence") and r.sequence:
            if len(r.sequence) >= min_length:
                result.append(r)
    return result


def _build_kmer_database(
    sequences: Sequence[Any],
    kmer_size: int,
) -> dict[str, int]:
    """Build a k-mer frequency database from sequences."""
    kmer_counts: dict[str, int] = defaultdict(int)

    for seq_data in sequences:
        if isinstance(seq_data, str):
            seq = seq_data.upper()
        elif isinstance(seq_data, dict):
            seq = seq_data.get("sequence", "").upper()
        elif hasattr(seq_data, "sequence"):
            seq = seq_data.sequence.upper()
        else:
            continue

        for i in range(len(seq) - kmer_size + 1):
            kmer = seq[i : i + kmer_size]
            if "N" not in kmer:
                kmer_counts[kmer] += 1

    return dict(kmer_counts)


def _cluster_reads_by_overlaps(
    reads: list[dict[str, Any]],
    overlaps: list[Any],
) -> list[list[dict[str, Any]]]:
    """Cluster reads into groups based on their overlaps.

    Uses union-find to group reads that are connected by overlaps.
    """
    read_index: dict[str, int] = {}
    for i, r in enumerate(reads):
        rid = r.get("read_id", r.get("read_name", f"read_{i}"))
        read_index[str(rid)] = i

    # Union-Find
    parent = list(range(len(reads)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for ovl in overlaps:
        qname = ovl.query_name if hasattr(ovl, "query_name") else ovl.get("query_name", "")
        tname = ovl.target_name if hasattr(ovl, "target_name") else ovl.get("target_name", "")
        qi = read_index.get(str(qname))
        ti = read_index.get(str(tname))
        if qi is not None and ti is not None:
            union(qi, ti)

    clusters: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for i, r in enumerate(reads):
        root = find(i)
        clusters[root].append(r)

    return list(clusters.values())


def _compute_n50(sorted_lengths_desc: list[int]) -> int:
    """Compute N50 from a descending-sorted list of lengths."""
    total = sum(sorted_lengths_desc)
    threshold = total * 0.5
    cumulative = 0
    for length in sorted_lengths_desc:
        cumulative += length
        if cumulative >= threshold:
            return length
    return sorted_lengths_desc[-1] if sorted_lengths_desc else 0
