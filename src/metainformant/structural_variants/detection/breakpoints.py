"""Breakpoint detection and refinement for structural variants.

Provides algorithms for refining breakpoint positions to base-pair
resolution, detecting microhomology at breakpoint junctions, clustering
nearby breakpoints, and calculating breakpoint confidence scores.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]

logger = get_logger(__name__)


@dataclass
class Breakpoint:
    """A structural variant breakpoint.

    Attributes:
        chrom: Chromosome name.
        position: Genomic position (0-based).
        strand: Strand orientation ('+' or '-').
        confidence: Confidence score [0, 1].
        support: Number of supporting reads.
        microhomology: Microhomology sequence at the breakpoint, if any.
        inserted_sequence: Inserted sequence at the breakpoint, if any.
        read_names: Names of reads supporting this breakpoint.
    """

    chrom: str
    position: int
    strand: str = "+"
    confidence: float = 0.0
    support: int = 0
    microhomology: str = ""
    inserted_sequence: str = ""
    read_names: list[str] = field(default_factory=list)


@dataclass
class BreakpointPair:
    """A pair of breakpoints defining a structural variant.

    Attributes:
        bp1: First breakpoint.
        bp2: Second breakpoint.
        sv_type: Inferred structural variant type.
        confidence: Joint confidence score.
    """

    bp1: Breakpoint
    bp2: Breakpoint
    sv_type: str = "UNKNOWN"
    confidence: float = 0.0


def refine_breakpoints(
    variants: list[dict[str, Any]],
    reads: list[dict[str, Any]],
    window: int = 200,
    min_clip: int = 5,
) -> list[BreakpointPair]:
    """Refine breakpoint positions to base-pair resolution.

    For each variant, examines split reads near the reported breakpoint
    positions and determines the most likely exact breakpoint positions
    by consensus of soft-clip positions.

    Args:
        variants: List of variant dictionaries, each containing:
            - 'chrom' or 'chrom1': Chromosome of breakpoint 1
            - 'start' or 'breakpoint1': Approximate position of breakpoint 1
            - 'end' or 'breakpoint2': Approximate position of breakpoint 2
            - 'chrom2': Chromosome of breakpoint 2 (optional, defaults to chrom)
            - 'sv_type': Type of structural variant (optional)
        reads: List of read alignment dictionaries with:
            - 'chrom': Chromosome
            - 'pos': Alignment start position
            - 'cigar': CIGAR string
            - 'name': Read name
            - 'seq': Read sequence (optional, for microhomology detection)
        window: Window size around each breakpoint to search for refinement.
        min_clip: Minimum soft-clip length to consider.

    Returns:
        List of BreakpointPair objects with refined positions.
    """
    refined_pairs: list[BreakpointPair] = []

    for variant in variants:
        chrom1 = variant.get("chrom1", variant.get("chrom", ""))
        chrom2 = variant.get("chrom2", chrom1)
        bp1_approx = variant.get("breakpoint1", variant.get("start", 0))
        bp2_approx = variant.get("breakpoint2", variant.get("end", 0))
        sv_type = variant.get("sv_type", "UNKNOWN")

        # Collect soft-clip positions near breakpoint 1
        bp1_clips = _collect_clip_positions(reads, chrom1, bp1_approx, window, min_clip)
        # Collect soft-clip positions near breakpoint 2
        bp2_clips = _collect_clip_positions(reads, chrom2, bp2_approx, window, min_clip)

        # Refine breakpoint 1 position using consensus
        bp1_pos, bp1_support, bp1_names = _consensus_position(bp1_clips, bp1_approx)
        bp1_conf = calculate_breakpoint_confidence(
            {"support": bp1_support, "total_clips": len(bp1_clips), "position_std": _position_std(bp1_clips)}
        )

        # Refine breakpoint 2 position using consensus
        bp2_pos, bp2_support, bp2_names = _consensus_position(bp2_clips, bp2_approx)
        bp2_conf = calculate_breakpoint_confidence(
            {"support": bp2_support, "total_clips": len(bp2_clips), "position_std": _position_std(bp2_clips)}
        )

        # Detect microhomology at breakpoints
        bp1_mh = ""
        bp2_mh = ""
        for read in reads:
            seq = read.get("seq", "")
            if seq and read.get("chrom", "") == chrom1:
                mh = detect_microhomology(
                    {"chrom": chrom1, "position": bp1_pos},
                    seq,
                )
                if len(mh) > len(bp1_mh):
                    bp1_mh = mh
                break  # One sequence is enough for microhomology check

        bp1 = Breakpoint(
            chrom=chrom1,
            position=bp1_pos,
            strand="+",
            confidence=bp1_conf,
            support=bp1_support,
            microhomology=bp1_mh,
            read_names=bp1_names,
        )

        bp2 = Breakpoint(
            chrom=chrom2,
            position=bp2_pos,
            strand="+",
            confidence=bp2_conf,
            support=bp2_support,
            microhomology=bp2_mh,
            read_names=bp2_names,
        )

        joint_conf = (bp1_conf + bp2_conf) / 2.0

        refined_pairs.append(
            BreakpointPair(
                bp1=bp1,
                bp2=bp2,
                sv_type=sv_type if isinstance(sv_type, str) else str(sv_type),
                confidence=joint_conf,
            )
        )

    logger.info(f"Refined {len(refined_pairs)} breakpoint pairs from {len(variants)} variants")
    return refined_pairs


def detect_microhomology(
    breakpoint: dict[str, Any] | Breakpoint,
    flanking_seq: str,
    max_homology: int = 50,
) -> str:
    """Detect microhomology at a breakpoint junction.

    Examines the sequence around a breakpoint to find short regions of
    sequence identity (microhomology), which are indicative of the
    repair mechanism (e.g., NHEJ, MMEJ/alt-EJ).

    The algorithm slides a window from the breakpoint position in both
    directions, comparing bases for matches. The longest run of matching
    bases starting from the breakpoint is reported.

    Args:
        breakpoint: Breakpoint position info. Can be a dict with 'position'
            key or a Breakpoint object.
        flanking_seq: Sequence flanking the breakpoint. The breakpoint is
            assumed to be at the midpoint of this sequence.
        max_homology: Maximum microhomology length to search for.

    Returns:
        Microhomology sequence string (empty if none found).
    """
    if not flanking_seq:
        return ""

    if isinstance(breakpoint, Breakpoint):
        bp_pos = breakpoint.position
    else:
        bp_pos = breakpoint.get("position", 0)

    seq = flanking_seq.upper()
    mid = len(seq) // 2

    if mid <= 0 or mid >= len(seq):
        return ""

    # Compare bases extending from the breakpoint in both directions
    # Microhomology: identical bases on both sides of the junction
    homology: list[str] = []

    for i in range(min(max_homology, mid, len(seq) - mid)):
        left_idx = mid - 1 - i
        right_idx = mid + i

        if left_idx < 0 or right_idx >= len(seq):
            break

        if seq[left_idx] == seq[right_idx]:
            homology.append(seq[right_idx])
        else:
            break

    if not homology:
        return ""

    # Also check forward microhomology (same bases on each side of break)
    forward_mh: list[str] = []
    for i in range(min(max_homology, len(seq) - mid)):
        left_idx = mid - 1 - i
        right_idx = mid + i

        if left_idx < 0 or right_idx >= len(seq):
            break

        left_base = seq[left_idx]
        right_base = seq[right_idx]

        if left_base == right_base and left_base in "ACGT":
            forward_mh.append(left_base)
        else:
            break

    # Return the longer homology
    result = "".join(homology) if len(homology) >= len(forward_mh) else "".join(forward_mh)
    return result


def cluster_breakpoints(
    breakpoints: list[Breakpoint] | list[dict[str, Any]],
    max_distance: int = 100,
) -> list[list[Breakpoint]]:
    """Cluster nearby breakpoints using single-linkage clustering.

    Groups breakpoints that are on the same chromosome and within
    max_distance base pairs of each other. Returns clusters sorted
    by cluster size (largest first).

    Args:
        breakpoints: List of Breakpoint objects or dictionaries with
            'chrom' and 'position' keys.
        max_distance: Maximum distance between breakpoints to be
            assigned to the same cluster.

    Returns:
        List of breakpoint clusters, each cluster is a list of Breakpoint
        objects. Clusters are sorted by size (descending).
    """
    if not breakpoints:
        return []

    # Normalize to Breakpoint objects
    bp_list: list[Breakpoint] = []
    for bp in breakpoints:
        if isinstance(bp, Breakpoint):
            bp_list.append(bp)
        else:
            bp_list.append(
                Breakpoint(
                    chrom=bp.get("chrom", ""),
                    position=bp.get("position", 0),
                    strand=bp.get("strand", "+"),
                    confidence=bp.get("confidence", 0.0),
                    support=bp.get("support", 0),
                    read_names=bp.get("read_names", []),
                )
            )

    # Sort by chromosome and position
    bp_list.sort(key=lambda b: (b.chrom, b.position))

    # Single-linkage clustering
    clusters: list[list[Breakpoint]] = []
    current_cluster: list[Breakpoint] = [bp_list[0]]

    for bp in bp_list[1:]:
        # Check if this breakpoint is close enough to any in current cluster
        can_merge = False
        for existing in current_cluster:
            if bp.chrom == existing.chrom and abs(bp.position - existing.position) <= max_distance:
                can_merge = True
                break

        if can_merge:
            current_cluster.append(bp)
        else:
            clusters.append(current_cluster)
            current_cluster = [bp]

    clusters.append(current_cluster)

    # Sort clusters by size descending
    clusters.sort(key=len, reverse=True)

    logger.debug(f"Clustered {len(bp_list)} breakpoints into {len(clusters)} clusters")
    return clusters


def calculate_breakpoint_confidence(
    evidence: dict[str, Any] | Any,
) -> float:
    """Calculate confidence score for a breakpoint.

    Combines multiple evidence metrics into a single confidence score:
    - Number of supporting reads (more = higher confidence)
    - Positional consistency (lower std = higher confidence)
    - Ratio of supporting to total nearby split reads

    Args:
        evidence: Dictionary or object containing:
            - 'support': Number of reads supporting the exact position (int)
            - 'total_clips': Total split reads in the region (int, optional)
            - 'position_std': Standard deviation of clip positions (float, optional)
            - 'mapq_mean': Mean mapping quality of supporting reads (float, optional)

    Returns:
        Confidence score between 0 and 1.
    """
    if isinstance(evidence, dict):
        support = evidence.get("support", 0)
        total_clips = evidence.get("total_clips", support)
        position_std = evidence.get("position_std", 10.0)
        mapq_mean = evidence.get("mapq_mean", 30.0)
    else:
        support = getattr(evidence, "support", 0)
        total_clips = getattr(evidence, "total_clips", support)
        position_std = getattr(evidence, "position_std", 10.0)
        mapq_mean = getattr(evidence, "mapq_mean", 30.0)

    if support <= 0:
        return 0.0

    # Component 1: Support count (sigmoid scaling, saturates around 20)
    support_score = 1.0 - 1.0 / (1.0 + support / 5.0)

    # Component 2: Positional consistency (lower std is better)
    if position_std <= 0:
        position_score = 1.0
    else:
        position_score = 1.0 / (1.0 + position_std / 5.0)

    # Component 3: Fraction of local clips supporting this position
    if total_clips > 0:
        fraction_score = support / total_clips
    else:
        fraction_score = 0.5

    # Component 4: Mapping quality
    mapq_score = min(1.0, mapq_mean / 60.0)

    # Weighted combination
    confidence = (
        0.35 * support_score
        + 0.30 * position_score
        + 0.20 * fraction_score
        + 0.15 * mapq_score
    )

    return max(0.0, min(1.0, confidence))


def _collect_clip_positions(
    reads: list[dict[str, Any]],
    chrom: str,
    position: int,
    window: int,
    min_clip: int,
) -> list[tuple[int, str]]:
    """Collect soft-clip positions near a target position.

    Args:
        reads: List of read alignment dictionaries.
        chrom: Target chromosome.
        position: Target position.
        window: Window size around the position.
        min_clip: Minimum clip length.

    Returns:
        List of (clip_position, read_name) tuples.
    """
    clips: list[tuple[int, str]] = []

    for read in reads:
        if read.get("chrom", "") != chrom:
            continue

        read_pos = read.get("pos", 0)
        if abs(read_pos - position) > window + 500:
            continue

        cigar = read.get("cigar", "")
        if not cigar:
            continue

        name = read.get("name", "")
        parsed_clips = _parse_clips_from_cigar(cigar, read_pos, min_clip)

        for clip_pos, _clip_len in parsed_clips:
            if abs(clip_pos - position) <= window:
                clips.append((clip_pos, name))

    return clips


def _parse_clips_from_cigar(
    cigar: str,
    pos: int,
    min_clip: int,
) -> list[tuple[int, int]]:
    """Parse CIGAR to find soft-clip positions and lengths.

    Args:
        cigar: CIGAR string.
        pos: Alignment start position.
        min_clip: Minimum clip length.

    Returns:
        List of (clip_position, clip_length) tuples.
    """
    clips: list[tuple[int, int]] = []
    op_map = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}

    ops: list[tuple[int, int]] = []
    num_str = ""
    for char in cigar:
        if char.isdigit():
            num_str += char
        elif char in op_map:
            if num_str:
                ops.append((op_map[char], int(num_str)))
            num_str = ""

    ref_offset = 0
    for i, (op, length) in enumerate(ops):
        if op == 4 and length >= min_clip:  # Soft clip
            if i == 0:
                clips.append((pos, length))
            else:
                clips.append((pos + ref_offset, length))
        elif op in (0, 2, 3, 7, 8):  # Reference-consuming
            ref_offset += length

    return clips


def _consensus_position(
    clips: list[tuple[int, str]],
    fallback: int,
) -> tuple[int, int, list[str]]:
    """Find consensus breakpoint position from clip positions.

    Uses the mode (most frequent position) as the consensus. Falls back
    to the provided position if no clips are found.

    Args:
        clips: List of (position, read_name) tuples.
        fallback: Fallback position if no clips are available.

    Returns:
        Tuple of (consensus_position, support_count, read_names).
    """
    if not clips:
        return fallback, 0, []

    # Count position frequencies
    pos_counts: dict[int, list[str]] = {}
    for pos, name in clips:
        if pos not in pos_counts:
            pos_counts[pos] = []
        pos_counts[pos].append(name)

    # Find the position with most support
    best_pos = max(pos_counts, key=lambda p: len(pos_counts[p]))
    support = len(pos_counts[best_pos])
    names = pos_counts[best_pos]

    return best_pos, support, names


def _position_std(clips: list[tuple[int, str]]) -> float:
    """Calculate standard deviation of clip positions.

    Args:
        clips: List of (position, read_name) tuples.

    Returns:
        Standard deviation of positions. Returns 0 if fewer than 2 clips.
    """
    if len(clips) < 2:
        return 0.0

    positions = [c[0] for c in clips]

    if np is not None:
        return float(np.std(positions, ddof=1))

    mean = sum(positions) / len(positions)
    variance = sum((p - mean) ** 2 for p in positions) / (len(positions) - 1)
    return math.sqrt(variance)
