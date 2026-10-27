"""Breakpoint detection and refinement for structural variants.

Provides algorithms for refining breakpoint positions to base-pair
resolution, detecting microhomology at breakpoint junctions, clustering
nearby breakpoints, and calculating breakpoint confidence scores.
"""

from __future__ import annotations

import bisect
import math
from collections import Counter
from dataclasses import dataclass, field
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

np: Any  # None when numpy is unavailable (optional dependency)
try:
    import numpy as np
except ImportError:
    np = None

logger = get_logger(__name__)


@dataclass
class Breakpoint:
    """A structural variant breakpoint.

    Attributes:
        chrom: Chromosome name.
        position: Genomic position (0-based).
        strand: Strand orientation ('+' or '-'): the orientation of most of the
            reads supporting the breakpoint ('-' when more are reverse-complemented
            than not). '+' when no supporting read reports ``is_reverse``.
        confidence: Confidence score [0, 1].
        support: Number of supporting reads.
        microhomology: Microhomology sequence at the breakpoint, if any.
        inserted_sequence: Bases inserted between the two alignments of a split
            read at this breakpoint, as they read along the forward strand of the
            alignment at this breakpoint; empty when not derivable (see
            :func:`refine_breakpoints`).
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


# Longest microhomology reported per breakpoint (the default of :func:`detect_microhomology`).
_MAX_MICROHOMOLOGY = 50


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

    A read is located by the reference coordinate of its soft clip, not by where
    its alignment starts, so a read that spans the breakpoint from far away (a long
    alignment, or a long deletion or skip in its CIGAR) counts as long as the clip
    itself lies within ``window`` of the breakpoint.

    The reads that clip at the consensus position are the breakpoint's supporting
    reads. Their orientation, sequence and supplementary alignment give each
    breakpoint's ``strand``, ``microhomology`` and ``inserted_sequence``; both
    breakpoints of a variant are analysed the same way.

    Args:
        variants: List of variant dictionaries, each containing:
            - 'chrom' or 'chrom1': Chromosome of breakpoint 1
            - 'start' or 'breakpoint1': Approximate position of breakpoint 1
            - 'end' or 'breakpoint2': Approximate position of breakpoint 2
            - 'pos2': Position of breakpoint 2 on 'chrom2' (translocations; it
              takes the place of 'end', which then describes breakpoint 1 only).
              ``dataclasses.asdict`` of a StructuralVariant provides it.
            - 'chrom2': Chromosome of breakpoint 2 (optional; absent, None or ''
              means the same chromosome as breakpoint 1)
            - 'sv_type': Type of structural variant (optional)
        reads: List of read alignment dictionaries with:
            - 'chrom': Chromosome
            - 'pos': Alignment start position
            - 'cigar': CIGAR string
            - 'name': Read name
            - 'is_reverse': Whether the read is reverse-complemented (optional;
              gives the breakpoint strand)
            - 'seq': Read sequence as stored in the alignment (optional; used for
              microhomology and inserted-sequence detection)
            - 'sa_tag': SA tag of a read with exactly one supplementary
              alignment (optional; with 'seq' and 'is_reverse' it gives the
              inserted sequence)
        window: Window size around each breakpoint to search for refinement.
        min_clip: Minimum soft-clip length to consider.

    Returns:
        List of BreakpointPair objects with refined positions. A breakpoint with
        no soft clip in its window keeps its input position with ``support`` 0.
    """
    clip_index = _ClipIndex(reads, min_clip)
    refined_pairs: list[BreakpointPair] = []

    for variant in variants:
        chrom1 = variant.get("chrom1") or variant.get("chrom", "")
        chrom2 = variant.get("chrom2") or chrom1
        bp1_approx = _first_present(variant, ("breakpoint1", "start"))
        bp2_approx = _first_present(variant, ("breakpoint2", "pos2", "end"))
        sv_type = variant.get("sv_type", "UNKNOWN")

        bp1 = _refine_breakpoint(
            chrom1, bp1_approx, clip_index.near(chrom1, bp1_approx, window), chrom2, bp2_approx, window
        )
        bp2 = _refine_breakpoint(
            chrom2, bp2_approx, clip_index.near(chrom2, bp2_approx, window), chrom1, bp1_approx, window
        )

        joint_conf = (bp1.confidence + bp2.confidence) / 2.0

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
        support: int = evidence.get("support", 0)
        total_clips: int = evidence.get("total_clips", support)
        position_std: float = evidence.get("position_std", 10.0)
        mapq_mean: float = evidence.get("mapq_mean", 30.0)
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
    confidence = 0.35 * support_score + 0.30 * position_score + 0.20 * fraction_score + 0.15 * mapq_score

    return max(0.0, min(1.0, confidence))


@dataclass
class _ClipHit:
    """One soft-clipped junction of one read (internal).

    Attributes:
        position: Reference coordinate of the junction.
        name: Read name.
        clip_length: Length of the soft clip.
        query_offset: Offset of the junction in the read's ``seq`` (bases before it).
        order: ``(read index, clip index)``: the order in which the clips were read.
        read: The read dictionary the clip belongs to.
    """

    position: int
    name: str
    clip_length: int
    query_offset: int
    order: tuple[int, int]
    read: dict[str, Any]


class _ClipIndex:
    """Soft-clipped junctions of a read set, sorted by reference position per chromosome.

    Every read is parsed once and found by the coordinate of its clip, never by where
    its alignment starts: a long alignment (or a long deletion or skip in its CIGAR)
    puts the clip far from the read start, so the start position says nothing about
    which reads clip near a breakpoint.
    """

    def __init__(self, reads: list[dict[str, Any]], min_clip: int) -> None:
        by_chrom: dict[str, list[_ClipHit]] = {}
        for read_index, read in enumerate(reads):
            cigar = read.get("cigar", "")
            if not cigar:
                continue
            chrom = read.get("chrom", "")
            name = read.get("name", "")
            for clip_index, (clip_pos, clip_len, offset) in enumerate(
                _parse_clips_from_cigar(cigar, read.get("pos", 0), min_clip)
            ):
                by_chrom.setdefault(chrom, []).append(
                    _ClipHit(clip_pos, name, clip_len, offset, (read_index, clip_index), read)
                )

        self._hits: dict[str, list[_ClipHit]] = {}
        self._positions: dict[str, list[int]] = {}
        for chrom, hits in by_chrom.items():
            hits.sort(key=lambda hit: hit.position)
            self._hits[chrom] = hits
            self._positions[chrom] = [hit.position for hit in hits]

    def near(self, chrom: str, position: int, window: int) -> list[_ClipHit]:
        """Return the clips within ``window`` bp of ``position``, in the order they were read."""
        positions = self._positions.get(chrom)
        if not positions:
            return []
        low = bisect.bisect_left(positions, position - window)
        high = bisect.bisect_right(positions, position + window)
        return sorted(self._hits[chrom][low:high], key=lambda hit: hit.order)


def _first_present(variant: dict[str, Any], keys: tuple[str, ...]) -> int:
    """Return the first value among ``keys`` that is not None (0 when none is)."""
    for key in keys:
        value = variant.get(key)
        if value is not None:
            return int(value)
    return 0


def _refine_breakpoint(
    chrom: str,
    approx: int,
    hits: list[_ClipHit],
    other_chrom: str,
    other_approx: int,
    window: int,
) -> Breakpoint:
    """Build the refined breakpoint from the soft clips found near it.

    Args:
        chrom: Chromosome of this breakpoint.
        approx: Approximate position (kept when there is no clip).
        hits: Soft clips within the window around ``approx``.
        other_chrom: Chromosome of the variant's other breakpoint.
        other_approx: Approximate position of the other breakpoint.
        window: Search window (also the tolerance for locating the other breakpoint).

    Returns:
        The refined Breakpoint.
    """
    clips = [(hit.position, hit.name) for hit in hits]
    position, support, names = _consensus_position(clips, approx)
    confidence = calculate_breakpoint_confidence(
        {"support": support, "total_clips": len(clips), "position_std": _position_std(clips)}
    )
    supporters = [hit for hit in hits if hit.position == position] if support else []

    return Breakpoint(
        chrom=chrom,
        position=position,
        strand=_majority_strand(supporters),
        confidence=confidence,
        support=support,
        microhomology=_junction_microhomology(chrom, position, supporters),
        inserted_sequence=_split_read_insertion(supporters, other_chrom, other_approx, window),
        read_names=names,
    )


def _majority_strand(hits: list[_ClipHit]) -> str:
    """Orientation of most supporting reads: '-' when more are reverse-complemented than not.

    Reads that do not report ``is_reverse`` carry no orientation and are not counted;
    with no orientation at all (or a tie) the breakpoint keeps the default '+'.
    """
    reverse = forward = 0
    for hit in hits:
        flag = hit.read.get("is_reverse")
        if flag is None:
            continue
        if flag:
            reverse += 1
        else:
            forward += 1
    return "-" if reverse > forward else "+"


def _junction_microhomology(chrom: str, position: int, hits: list[_ClipHit]) -> str:
    """Longest microhomology among the supporting reads, each centred on its junction.

    :func:`detect_microhomology` assumes the breakpoint is at the midpoint of the
    sequence it is given, so each read contributes the window of its ``seq`` that is
    centred on its own clip junction (and not the whole read).
    """
    best = ""
    for hit in hits:
        seq = hit.read.get("seq")
        if not seq:
            continue
        half = min(_MAX_MICROHOMOLOGY, hit.query_offset, len(seq) - hit.query_offset)
        if half < 1:
            continue
        flank = seq[hit.query_offset - half : hit.query_offset + half]
        homology = detect_microhomology({"chrom": chrom, "position": position}, flank)
        if len(homology) > len(best):
            best = homology
    return best


def _split_read_insertion(hits: list[_ClipHit], other_chrom: str, other_approx: int, window: int) -> str:
    """Most common inserted sequence among the supporting reads, or '' when not derivable.

    Nothing but the read's own two alignments tells inserted bases from reference bases,
    so a sequence is reported only for a read that has its ``seq``, its ``is_reverse``
    flag and an ``sa_tag`` with exactly one supplementary alignment located at the other
    breakpoint. Ties go to the longer sequence, then to the first read.
    """
    found: Counter[str] = Counter()
    for hit in hits:
        inserted = _inserted_between_alignments(hit.read, other_chrom, other_approx, window)
        if inserted:
            found[inserted] += 1
    if not found:
        return ""
    return max(found.items(), key=lambda item: (item[1], len(item[0])))[0]


def _inserted_between_alignments(read: dict[str, Any], other_chrom: str, other_approx: int, window: int) -> str:
    """Bases of ``read`` between its primary and its single supplementary alignment.

    Both alignments cover an interval of the read (the CIGAR, leading clip to the end
    of the aligned part). Expressed in the orientation of the stored ``seq``, the bases
    between the two intervals are inserted; overlapping or abutting alignments
    insert nothing. An SA entry on the opposite strand describes the reverse-complemented
    read, so its interval is mirrored (its leading clip becomes its trailing clip).

    Args:
        read: Read dictionary with ``seq``, ``cigar``, ``is_reverse`` and ``sa_tag``.
        other_chrom: Chromosome the supplementary alignment must lie on.
        other_approx: Position the supplementary alignment must lie within ``window`` of.
        window: Positional tolerance in bp.

    Returns:
        The inserted sequence, or '' when it cannot be derived.
    """
    seq = read.get("seq", "")
    sa_tag = read.get("sa_tag", "")
    is_reverse = read.get("is_reverse")
    if not seq or not sa_tag or is_reverse is None:
        return ""

    primary = _query_extent(read.get("cigar", ""))
    entries = _parse_sa_alignments(sa_tag)
    if primary is None or len(entries) != 1 or sum(primary[:3]) != len(seq):
        return ""

    sa_chrom, sa_start, sa_strand, sa_cigar = entries[0]
    supplementary = _query_extent(sa_cigar)
    if supplementary is None or sum(supplementary[:3]) != len(seq):
        return ""
    sa_lead, sa_core, sa_trail, sa_ref_length = supplementary
    if sa_chrom != other_chrom or not (sa_start - window <= other_approx <= sa_start + sa_ref_length + window):
        return ""

    lead, core, _trail, _ref_length = primary
    primary_strand = "-" if is_reverse else "+"
    sa_begin = sa_lead if sa_strand == primary_strand else sa_trail
    sa_end = sa_begin + sa_core

    if sa_begin >= lead + core:
        return str(seq[lead + core : sa_begin])
    if lead >= sa_end:
        return str(seq[sa_end:lead])
    return ""


_CIGAR_LETTERS = "MIDNSHP=X"  # BAM operation codes 0-8


def _cigar_ops(cigar: str | Sequence[tuple[int, int]]) -> list[tuple[str, int]]:
    """Split a CIGAR into ``(operation letter, length)`` pairs.

    Args:
        cigar: A CIGAR string (``"50M30S"``) or a list of ``(BAM operation code, length)``
            tuples (``[(0, 50), (4, 30)]``).
    """
    if not isinstance(cigar, str):
        return [(_CIGAR_LETTERS[code], int(length)) for code, length in cigar]

    ops: list[tuple[str, int]] = []
    num_str = ""
    for char in cigar:
        if char.isdigit():
            num_str += char
        elif char in _CIGAR_LETTERS:
            if num_str:
                ops.append((char, int(num_str)))
            num_str = ""
    return ops


def _query_extent(cigar: str | Sequence[tuple[int, int]]) -> tuple[int, int, int, int] | None:
    """Return ``(leading clip, aligned query length, trailing clip, reference length)`` of a CIGAR.

    The leading and trailing clips count soft and hard clipped bases; the aligned query
    length counts M, I, = and X; the reference length counts M, D, N, = and X.
    Returns None for an empty or unparseable CIGAR.
    """
    ops = _cigar_ops(cigar)
    if not ops:
        return None
    lead = core = trail = ref_length = 0
    seen_aligned = False
    for op, length in ops:
        if op in "SH":
            if seen_aligned:
                trail += length
            else:
                lead += length
        elif op in "MIDN=X":
            if op in "MI=X":
                core += length
            if op in "MDN=X":
                ref_length += length
            seen_aligned = True
    return lead, core, trail, ref_length


def _parse_sa_alignments(sa_tag: str) -> list[tuple[str, int, str, str]]:
    """Parse an SA tag (``rname,pos,strand,CIGAR,mapQ,NM;...``) into ``(chrom, start, strand, CIGAR)``.

    ``start`` is the 0-based leftmost reference position. Entries without a numeric
    position or a CIGAR are skipped.
    """
    entries: list[tuple[str, int, str, str]] = []
    for entry in sa_tag.rstrip(";").split(";"):
        parts = entry.strip().split(",")
        if len(parts) < 4:
            continue
        try:
            start = int(parts[1]) - 1
        except ValueError:
            continue
        entries.append((parts[0], start, parts[2], parts[3]))
    return entries


def _parse_clips_from_cigar(
    cigar: str | Sequence[tuple[int, int]],
    pos: int,
    min_clip: int,
) -> list[tuple[int, int, int]]:
    """Parse CIGAR to find soft-clip junctions.

    Args:
        cigar: CIGAR string, or a list of ``(BAM operation code, length)`` tuples.
        pos: Alignment start position.
        min_clip: Minimum clip length.

    Returns:
        List of ``(clip_position, clip_length, query_offset)`` tuples: the reference
        coordinate of the junction between the clip and the aligned part, the clip
        length, and the offset of that junction in the read's stored sequence (hard
        clipped bases are not stored, so they are not counted).
    """
    clips: list[tuple[int, int, int]] = []
    ref_offset = 0
    query_offset = 0
    seen_aligned = False
    for op, length in _cigar_ops(cigar):
        if op == "S":
            if length >= min_clip:
                if seen_aligned:  # right clip: the junction lies where the aligned part ends
                    clips.append((pos + ref_offset, length, query_offset))
                else:  # left clip: the junction lies where the aligned part begins
                    clips.append((pos, length, query_offset + length))
            query_offset += length
        elif op in "MIDN=X":
            if op in "MDN=X":
                ref_offset += length
            if op in "MI=X":
                query_offset += length
            seen_aligned = True
        # H and P consume neither reference nor stored query bases

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
