"""Structural variant calling from aligned sequencing reads.

Detects structural variants (deletions, duplications, inversions,
translocations, insertions) from split-read and discordant read-pair
evidence in BAM/SAM alignment data.

Supports both pysam-backed BAM reading and plain dictionary-based
read representations for flexibility and testability.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]

try:
    from scipy import stats as scipy_stats
except ImportError:
    scipy_stats = None  # type: ignore[assignment]

logger = get_logger(__name__)

_ENV_PREFIX = "SV_"


class SVType(str, Enum):
    """Structural variant types."""

    DEL = "DEL"
    DUP = "DUP"
    INV = "INV"
    TRA = "TRA"
    INS = "INS"
    BND = "BND"
    UNKNOWN = "UNKNOWN"


@dataclass
class InsertSizeStats:
    """Insert size distribution statistics.

    Attributes:
        mean: Mean insert size.
        std: Standard deviation of insert size.
        median: Median insert size.
        mad: Median absolute deviation.
    """

    mean: float
    std: float
    median: float = 0.0
    mad: float = 0.0


@dataclass
class SVEvidence:
    """Evidence supporting a structural variant call.

    Attributes:
        split_reads: Number of supporting split reads.
        discordant_pairs: Number of supporting discordant read pairs.
        depth_support: Whether read depth supports the call.
        evidence_reads: List of read names providing evidence.
        breakpoint1: Estimated position of first breakpoint.
        breakpoint2: Estimated position of second breakpoint.
        chrom1: Chromosome of first breakpoint.
        chrom2: Chromosome of second breakpoint.
        strand1: Strand orientation at breakpoint 1.
        strand2: Strand orientation at breakpoint 2.
    """

    split_reads: int = 0
    discordant_pairs: int = 0
    depth_support: bool = False
    evidence_reads: list[str] = field(default_factory=list)
    breakpoint1: int = 0
    breakpoint2: int = 0
    chrom1: str = ""
    chrom2: str = ""
    strand1: str = "+"
    strand2: str = "+"


@dataclass
class StructuralVariant:
    """A called structural variant.

    Attributes:
        chrom: Primary chromosome.
        start: Start position (0-based).
        end: End position (0-based).
        sv_type: Type of structural variant.
        size: Size of the variant in base pairs.
        quality: Quality score.
        evidence: Supporting evidence.
        genotype: Genotype string (e.g., '0/1', '1/1').
        filter_status: Filter status ('PASS' or reason).
        chrom2: Secondary chromosome (for translocations).
        info: Additional info fields.
    """

    chrom: str
    start: int
    end: int
    sv_type: SVType
    size: int = 0
    quality: float = 0.0
    evidence: SVEvidence = field(default_factory=SVEvidence)
    genotype: str = "./."
    filter_status: str = "PASS"
    chrom2: str = ""
    info: dict[str, Any] = field(default_factory=dict)


def call_structural_variants(
    alignments: list[dict[str, Any]],
    min_support: int = 3,
    insert_size_stats: InsertSizeStats | None = None,
    min_mapq: int = 20,
    min_sv_size: int = 50,
) -> list[StructuralVariant]:
    """Call structural variants from alignment data.

    Combines split-read and discordant read-pair evidence to detect
    structural variants. Reads are first classified as split or discordant,
    then evidence is clustered by genomic position, and variant calls are
    made from clusters with sufficient support.

    Args:
        alignments: List of read alignment dictionaries, each containing:
            - 'name': Read name (str)
            - 'chrom': Reference chromosome (str)
            - 'pos': Alignment start position, 0-based (int)
            - 'mapq': Mapping quality (int)
            - 'cigar': CIGAR string (str) or list of (op, length) tuples
            - 'mate_chrom': Mate chromosome (str, optional)
            - 'mate_pos': Mate start position (int, optional)
            - 'insert_size': Observed insert size (int, optional)
            - 'is_reverse': Whether read is reverse-complemented (bool, optional)
            - 'mate_is_reverse': Whether mate is reverse-complemented (bool, optional)
            - 'is_supplementary': Whether this is a supplementary alignment (bool, optional)
            - 'sa_tag': SA tag for supplementary alignments (str, optional)
        min_support: Minimum number of supporting reads for a call.
        insert_size_stats: Pre-computed insert size distribution statistics.
            If None, will be estimated from the data.
        min_mapq: Minimum mapping quality for reads to be considered.
        min_sv_size: Minimum structural variant size in base pairs.

    Returns:
        List of StructuralVariant objects, sorted by chromosome and position.
    """
    env_min_support = os.environ.get(f"{_ENV_PREFIX}MIN_SUPPORT")
    if env_min_support is not None:
        min_support = int(env_min_support)

    # Filter low-quality reads
    filtered = [a for a in alignments if a.get("mapq", 0) >= min_mapq]
    logger.info(f"Processing {len(filtered)} reads (filtered from {len(alignments)} by MAPQ >= {min_mapq})")

    # Estimate insert size stats if not provided
    if insert_size_stats is None:
        insert_size_stats = _estimate_insert_size(filtered)

    # Detect split-read evidence
    split_evidence = detect_split_reads(filtered, min_clip=20)

    # Detect discordant pair evidence
    discordant_evidence = detect_discordant_pairs(filtered, insert_size_stats)

    # Combine and cluster evidence
    all_evidence = split_evidence + discordant_evidence
    if not all_evidence:
        logger.info("No structural variant evidence found")
        return []

    # Cluster evidence by proximity
    clusters = _cluster_evidence(all_evidence, max_distance=500)

    # Build variant calls from clusters
    variants: list[StructuralVariant] = []
    for cluster in clusters:
        total_support = sum(e.split_reads + e.discordant_pairs for e in cluster)
        if total_support < min_support:
            continue

        # Merge evidence in this cluster
        merged = _merge_evidence(cluster)

        # Classify SV type
        sv_type = classify_sv_type(merged)

        # Calculate size
        if merged.chrom1 == merged.chrom2 or not merged.chrom2:
            size = abs(merged.breakpoint2 - merged.breakpoint1)
        else:
            size = 0  # Inter-chromosomal

        if size < min_sv_size and sv_type != SVType.TRA:
            continue

        # Calculate quality
        quality = _calculate_sv_quality(merged, total_support)

        variant = StructuralVariant(
            chrom=merged.chrom1,
            start=min(merged.breakpoint1, merged.breakpoint2),
            end=max(merged.breakpoint1, merged.breakpoint2),
            sv_type=sv_type,
            size=size,
            quality=quality,
            evidence=merged,
            chrom2=merged.chrom2 if merged.chrom2 != merged.chrom1 else "",
        )
        variants.append(variant)

    # Sort by position
    variants.sort(key=lambda v: (v.chrom, v.start))

    # Genotype all variants
    for variant in variants:
        variant.genotype = genotype_sv(variant, filtered)

    logger.info(f"Called {len(variants)} structural variants from {len(all_evidence)} evidence items")
    return variants


def detect_split_reads(
    reads: list[dict[str, Any]],
    min_clip: int = 20,
) -> list[SVEvidence]:
    """Identify split-read evidence for structural variants.

    Examines CIGAR strings and supplementary alignment information to
    detect reads that span structural variant breakpoints. A split read
    has a soft-clipped portion of at least min_clip bases, indicating
    the read crosses a breakpoint.

    Args:
        reads: List of read alignment dictionaries (see call_structural_variants
            for format).
        min_clip: Minimum soft-clip length (bases) to consider as evidence
            for a breakpoint.

    Returns:
        List of SVEvidence objects, one per split-read breakpoint detected.
    """
    evidence_list: list[SVEvidence] = []

    for read in reads:
        cigar = read.get("cigar", "")
        if not cigar:
            continue

        chrom = read.get("chrom", "")
        pos = read.get("pos", 0)
        name = read.get("name", "")

        # Parse CIGAR for soft clips
        clips = _parse_cigar_clips(cigar, pos)

        for clip_pos, clip_len, clip_side in clips:
            if clip_len < min_clip:
                continue

            # Check for supplementary alignment info (SA tag)
            sa_tag = read.get("sa_tag", "")
            if sa_tag:
                sa_parts = _parse_sa_tag(sa_tag)
                for sa_chrom, sa_pos, sa_strand in sa_parts:
                    evidence = SVEvidence(
                        split_reads=1,
                        evidence_reads=[name],
                        breakpoint1=clip_pos,
                        breakpoint2=sa_pos,
                        chrom1=chrom,
                        chrom2=sa_chrom,
                        strand1="+" if not read.get("is_reverse", False) else "-",
                        strand2=sa_strand,
                    )
                    evidence_list.append(evidence)
            else:
                # No SA tag: record the clip position as a potential breakpoint
                evidence = SVEvidence(
                    split_reads=1,
                    evidence_reads=[name],
                    breakpoint1=clip_pos,
                    breakpoint2=clip_pos,  # Will be refined during clustering
                    chrom1=chrom,
                    chrom2=chrom,
                    strand1="+" if not read.get("is_reverse", False) else "-",
                    strand2="+" if not read.get("is_reverse", False) else "-",
                )
                evidence_list.append(evidence)

    logger.debug(f"Detected {len(evidence_list)} split-read evidence items from {len(reads)} reads")
    return evidence_list


def detect_discordant_pairs(
    pairs: list[dict[str, Any]],
    insert_size_stats: InsertSizeStats,
    n_std: float = 4.0,
) -> list[SVEvidence]:
    """Detect discordant read pair evidence for structural variants.

    A discordant read pair has an insert size significantly different from
    the expected distribution, or has reads mapping to different chromosomes,
    or has unexpected orientation patterns.

    Args:
        pairs: List of read alignment dictionaries (see call_structural_variants
            for format). Only reads with mate information are analyzed.
        insert_size_stats: Insert size distribution statistics.
        n_std: Number of standard deviations from the mean to classify
            a pair as discordant (default 4.0).

    Returns:
        List of SVEvidence objects for discordant pairs.
    """
    evidence_list: list[SVEvidence] = []

    min_isize = max(0, insert_size_stats.mean - n_std * insert_size_stats.std)
    max_isize = insert_size_stats.mean + n_std * insert_size_stats.std

    seen_pairs: set[str] = set()

    for read in pairs:
        name = read.get("name", "")
        if name in seen_pairs:
            continue

        mate_chrom = read.get("mate_chrom")
        if mate_chrom is None:
            continue

        chrom = read.get("chrom", "")
        pos = read.get("pos", 0)
        mate_pos = read.get("mate_pos", 0)
        insert_size = read.get("insert_size", 0)
        is_reverse = read.get("is_reverse", False)
        mate_is_reverse = read.get("mate_is_reverse", False)

        is_discordant = False
        disc_reason = ""

        # Inter-chromosomal
        if mate_chrom != chrom:
            is_discordant = True
            disc_reason = "inter_chrom"
        # Aberrant insert size
        elif abs(insert_size) > max_isize or (abs(insert_size) < min_isize and abs(insert_size) > 0):
            is_discordant = True
            disc_reason = "insert_size"
        # Same-strand orientation (expected: FR for Illumina)
        elif is_reverse == mate_is_reverse:
            is_discordant = True
            disc_reason = "orientation"
        # Read-pair orientation reversed (RF instead of FR)
        elif is_reverse and not mate_is_reverse and pos < mate_pos:
            is_discordant = True
            disc_reason = "rf_orientation"

        if is_discordant:
            seen_pairs.add(name)
            evidence = SVEvidence(
                discordant_pairs=1,
                evidence_reads=[name],
                breakpoint1=pos,
                breakpoint2=mate_pos,
                chrom1=chrom,
                chrom2=mate_chrom,
                strand1="-" if is_reverse else "+",
                strand2="-" if mate_is_reverse else "+",
            )
            evidence_list.append(evidence)

    logger.debug(f"Detected {len(evidence_list)} discordant pairs from {len(pairs)} reads")
    return evidence_list


def classify_sv_type(evidence: SVEvidence) -> SVType:
    """Classify the type of structural variant based on evidence patterns.

    Uses breakpoint positions, chromosome assignments, and strand orientations
    to determine the SV type:
        - DEL: Same chromosome, same strand, forward-reverse pair, distant breakpoints
        - DUP: Same chromosome, reverse-forward orientation
        - INV: Same chromosome, same strand orientation
        - TRA: Different chromosomes
        - INS: Same chromosome, close breakpoints with split-read evidence

    Args:
        evidence: SVEvidence object containing breakpoint and orientation info.

    Returns:
        SVType enum value.
    """
    chrom1 = evidence.chrom1
    chrom2 = evidence.chrom2
    strand1 = evidence.strand1
    strand2 = evidence.strand2
    bp1 = evidence.breakpoint1
    bp2 = evidence.breakpoint2

    # Inter-chromosomal: translocation
    if chrom2 and chrom1 != chrom2:
        return SVType.TRA

    distance = abs(bp2 - bp1)

    # Very close breakpoints with split-read support: insertion
    if distance < 50 and evidence.split_reads > 0:
        return SVType.INS

    # Same-strand orientation: inversion
    if strand1 == strand2:
        return SVType.INV

    # Determine orientation relative to position order
    if bp1 <= bp2:
        # Normal ordering
        if strand1 == "+" and strand2 == "-":
            # Forward-reverse at distant positions: deletion
            return SVType.DEL
        elif strand1 == "-" and strand2 == "+":
            # Reverse-forward: tandem duplication
            return SVType.DUP
    else:
        # Reversed ordering
        if strand1 == "-" and strand2 == "+":
            return SVType.DEL
        elif strand1 == "+" and strand2 == "-":
            return SVType.DUP

    # Default classification based on insert size
    if distance > 1000:
        return SVType.DEL

    return SVType.UNKNOWN


def genotype_sv(
    variant: StructuralVariant,
    reads: list[dict[str, Any]],
    min_mapq: int = 20,
) -> str:
    """Genotype a structural variant based on supporting and reference reads.

    Counts the number of reads supporting the variant (split reads and
    discordant pairs) versus reads supporting the reference allele
    (reads spanning the breakpoint region without split alignment).
    Uses a simple allele-fraction model:
        - AF < 0.15: homozygous reference (0/0)
        - 0.15 <= AF < 0.85: heterozygous (0/1)
        - AF >= 0.85: homozygous alternate (1/1)

    Args:
        variant: StructuralVariant to genotype.
        reads: List of read alignment dictionaries.
        min_mapq: Minimum mapping quality for counting reads.

    Returns:
        Genotype string: '0/0', '0/1', or '1/1'.
    """
    # Count reads supporting the SV
    sv_support = variant.evidence.split_reads + variant.evidence.discordant_pairs

    # Count reads spanning the breakpoint region that support reference
    ref_support = 0
    bp_start = variant.start - 100
    bp_end = variant.start + 100

    for read in reads:
        if read.get("mapq", 0) < min_mapq:
            continue
        if read.get("chrom", "") != variant.chrom:
            continue

        read_pos = read.get("pos", 0)
        read_end = read_pos + read.get("read_length", 150)

        # Read spans the breakpoint without being split
        if read_pos < bp_start and read_end > bp_end:
            name = read.get("name", "")
            if name not in variant.evidence.evidence_reads:
                ref_support += 1

    total = sv_support + ref_support
    if total == 0:
        return "./."

    allele_fraction = sv_support / total

    if allele_fraction < 0.15:
        return "0/0"
    elif allele_fraction >= 0.85:
        return "1/1"
    else:
        return "0/1"


def _estimate_insert_size(reads: list[dict[str, Any]]) -> InsertSizeStats:
    """Estimate insert size distribution from properly paired reads.

    Samples insert sizes from reads that appear properly paired (same
    chromosome, FR orientation) and computes robust statistics.

    Args:
        reads: List of read alignment dictionaries.

    Returns:
        InsertSizeStats with estimated distribution parameters.
    """
    insert_sizes: list[float] = []

    for read in reads:
        isize = read.get("insert_size", 0)
        if isize is None:
            continue
        isize = abs(isize)
        if 0 < isize < 10000:  # Filter extreme outliers
            chrom = read.get("chrom", "")
            mate_chrom = read.get("mate_chrom", chrom)
            if chrom == mate_chrom:
                insert_sizes.append(float(isize))

    if not insert_sizes:
        # Default to typical Illumina insert size
        logger.warning("Could not estimate insert size, using defaults (mean=400, std=100)")
        return InsertSizeStats(mean=400.0, std=100.0, median=400.0, mad=50.0)

    if np is not None:
        arr = np.array(insert_sizes)
        mean = float(np.mean(arr))
        std = float(np.std(arr, ddof=1)) if len(arr) > 1 else 100.0
        median = float(np.median(arr))
        mad = float(np.median(np.abs(arr - median)))
    else:
        insert_sizes.sort()
        n = len(insert_sizes)
        mean = sum(insert_sizes) / n
        median = insert_sizes[n // 2]
        variance = sum((x - mean) ** 2 for x in insert_sizes) / max(n - 1, 1)
        std = math.sqrt(variance)
        deviations = sorted(abs(x - median) for x in insert_sizes)
        mad = deviations[n // 2]

    return InsertSizeStats(mean=mean, std=std, median=median, mad=mad)


def _parse_cigar_clips(cigar: str | list[tuple[int, int]], pos: int) -> list[tuple[int, int, str]]:
    """Parse CIGAR string to find soft-clipped regions.

    Args:
        cigar: CIGAR string (e.g., '50M30S') or list of (operation, length) tuples.
        pos: Alignment start position.

    Returns:
        List of (clip_position, clip_length, side) tuples where side is
        'left' or 'right'.
    """
    clips: list[tuple[int, int, str]] = []

    if isinstance(cigar, list):
        ops = cigar
    else:
        ops = _parse_cigar_string(cigar)

    if not ops:
        return clips

    ref_offset = 0

    for i, (op, length) in enumerate(ops):
        if op == 4:  # S = soft clip
            if i == 0:
                # Left soft clip: breakpoint at alignment start
                clips.append((pos, length, "left"))
            else:
                # Right soft clip: breakpoint at current reference position
                clips.append((pos + ref_offset, length, "right"))
        elif op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume reference
            ref_offset += length
        # I (1), S (4), H (5), P (6) do not consume reference

    return clips


def _parse_cigar_string(cigar: str) -> list[tuple[int, int]]:
    """Parse a CIGAR string into (operation, length) tuples.

    CIGAR operations: M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8

    Args:
        cigar: CIGAR string (e.g., '50M30S').

    Returns:
        List of (operation_code, length) tuples.
    """
    ops: list[tuple[int, int]] = []
    op_map = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}

    num_str = ""
    for char in cigar:
        if char.isdigit():
            num_str += char
        elif char in op_map:
            if num_str:
                ops.append((op_map[char], int(num_str)))
            num_str = ""
        else:
            num_str = ""

    return ops


def _parse_sa_tag(sa_tag: str) -> list[tuple[str, int, str]]:
    """Parse SAM SA (supplementary alignment) tag.

    SA tag format: rname,pos,strand,CIGAR,mapQ,NM;

    Args:
        sa_tag: SA tag string.

    Returns:
        List of (chromosome, position, strand) tuples.
    """
    results: list[tuple[str, int, str]] = []

    entries = sa_tag.rstrip(";").split(";")
    for entry in entries:
        parts = entry.strip().split(",")
        if len(parts) >= 3:
            chrom = parts[0]
            try:
                pos = int(parts[1]) - 1  # Convert to 0-based
            except ValueError:
                continue
            strand = parts[2]
            results.append((chrom, pos, strand))

    return results


def _cluster_evidence(
    evidence_list: list[SVEvidence],
    max_distance: int = 500,
) -> list[list[SVEvidence]]:
    """Cluster SV evidence by genomic proximity.

    Groups evidence items whose breakpoints are within max_distance of each
    other on the same chromosome pair, using single-linkage clustering.

    Args:
        evidence_list: List of SVEvidence objects.
        max_distance: Maximum distance between breakpoints for clustering.

    Returns:
        List of evidence clusters (each cluster is a list of SVEvidence).
    """
    if not evidence_list:
        return []

    # Sort by chrom1, chrom2, breakpoint1
    sorted_ev = sorted(evidence_list, key=lambda e: (e.chrom1, e.chrom2, e.breakpoint1))

    clusters: list[list[SVEvidence]] = [[sorted_ev[0]]]

    for ev in sorted_ev[1:]:
        merged = False
        # Try to add to existing cluster (check last cluster first for efficiency)
        for cluster in reversed(clusters):
            rep = cluster[0]
            if (
                ev.chrom1 == rep.chrom1
                and ev.chrom2 == rep.chrom2
                and abs(ev.breakpoint1 - rep.breakpoint1) <= max_distance
                and abs(ev.breakpoint2 - rep.breakpoint2) <= max_distance
            ):
                cluster.append(ev)
                merged = True
                break

        if not merged:
            clusters.append([ev])

    return clusters


def _merge_evidence(cluster: list[SVEvidence]) -> SVEvidence:
    """Merge a cluster of evidence into a single SVEvidence.

    Computes consensus breakpoint positions (median), totals support counts,
    and combines read names.

    Args:
        cluster: List of SVEvidence items to merge.

    Returns:
        Merged SVEvidence with consensus breakpoints and combined support.
    """
    if len(cluster) == 1:
        return cluster[0]

    total_sr = sum(e.split_reads for e in cluster)
    total_dp = sum(e.discordant_pairs for e in cluster)
    all_reads: list[str] = []
    for e in cluster:
        all_reads.extend(e.evidence_reads)

    bp1_values = [e.breakpoint1 for e in cluster]
    bp2_values = [e.breakpoint2 for e in cluster]

    if np is not None:
        bp1 = int(np.median(bp1_values))
        bp2 = int(np.median(bp2_values))
    else:
        bp1_values.sort()
        bp2_values.sort()
        bp1 = bp1_values[len(bp1_values) // 2]
        bp2 = bp2_values[len(bp2_values) // 2]

    # Use most common strand orientation
    strand1_counts: dict[str, int] = {}
    strand2_counts: dict[str, int] = {}
    for e in cluster:
        strand1_counts[e.strand1] = strand1_counts.get(e.strand1, 0) + 1
        strand2_counts[e.strand2] = strand2_counts.get(e.strand2, 0) + 1

    strand1 = max(strand1_counts, key=strand1_counts.get)  # type: ignore[arg-type]
    strand2 = max(strand2_counts, key=strand2_counts.get)  # type: ignore[arg-type]

    return SVEvidence(
        split_reads=total_sr,
        discordant_pairs=total_dp,
        evidence_reads=list(set(all_reads)),
        breakpoint1=bp1,
        breakpoint2=bp2,
        chrom1=cluster[0].chrom1,
        chrom2=cluster[0].chrom2,
        strand1=strand1,
        strand2=strand2,
    )


def _calculate_sv_quality(evidence: SVEvidence, total_support: int) -> float:
    """Calculate quality score for a structural variant call.

    Quality is derived from: total supporting evidence, balance between
    split reads and discordant pairs, and uniqueness of supporting reads.

    Args:
        evidence: Merged SVEvidence.
        total_support: Total support count.

    Returns:
        Quality score (Phred-scaled, capped at 999).
    """
    if total_support <= 0:
        return 0.0

    # Base quality from support count
    base_qual = min(60.0, 10.0 * math.log10(total_support + 1) * 10)

    # Bonus for having both types of evidence
    if evidence.split_reads > 0 and evidence.discordant_pairs > 0:
        base_qual += 10.0

    # Bonus for many unique reads
    n_unique = len(set(evidence.evidence_reads))
    unique_bonus = min(20.0, n_unique * 2.0)

    quality = base_qual + unique_bonus
    return min(999.0, max(0.0, quality))
