"""Structural variant detection from long-read sequencing data.

Detects structural variants (SVs) including insertions, deletions, inversions,
duplications, and translocations from split/supplementary long-read alignments.
Supports phasing of SVs to haplotypes.

The approach uses split-read and discordant alignment signatures from
long-read mappers like minimap2, which produce supplementary alignments
(SA tag) for reads spanning structural variant breakpoints.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


@dataclass
class StructuralVariant:
    """A detected structural variant.

    Attributes:
        sv_type: Type of SV ('DEL', 'INS', 'INV', 'DUP', 'BND').
        chromosome: Reference chromosome for the SV.
        start: Start position (0-based).
        end: End position (0-based). For BND, this is the mate position.
        size: Size of the variant in base pairs.
        supporting_reads: Number of reads supporting the SV.
        read_names: Names of supporting reads.
        quality: Variant quality score.
        genotype: Genotype string (e.g., '0/1', '1/1').
        haplotype: Haplotype assignment (0, 1, or 2; 0=unphased).
        inserted_sequence: Inserted sequence (for INS type).
        mate_chromosome: Mate chromosome (for BND type).
        mate_position: Mate breakpoint position (for BND type).
        info: Additional variant information.
    """

    sv_type: str
    chromosome: str
    start: int
    end: int
    size: int = 0
    supporting_reads: int = 0
    read_names: list[str] = field(default_factory=list)
    quality: float = 0.0
    genotype: str = "./."
    haplotype: int = 0
    inserted_sequence: str = ""
    mate_chromosome: str = ""
    mate_position: int = 0
    info: dict[str, Any] = field(default_factory=dict)


def detect_sv_from_long_reads(
    alignments: Sequence[Any],
    min_size: int = 50,
    min_support: int = 2,
    min_mapping_quality: int = 20,
) -> list[StructuralVariant]:
    """Detect structural variants from long-read alignments.

    Uses three signature types:
    1. CIGAR-based: Large insertions/deletions in the CIGAR string
    2. Split-read: Supplementary alignments indicating breakpoints
    3. Discordant: Unexpected alignment orientations

    Args:
        alignments: Sequence of alignment objects (LongReadAlignment or dicts).
        min_size: Minimum SV size in base pairs.
        min_support: Minimum number of supporting reads.
        min_mapping_quality: Minimum mapping quality for supporting alignments.

    Returns:
        List of StructuralVariant objects, sorted by position.
    """
    raw_calls: list[StructuralVariant] = []

    for aln in alignments:
        mapq = _get_attr(aln, "mapping_quality", 0)
        if mapq < min_mapping_quality:
            continue

        is_unmapped = _get_attr(aln, "is_unmapped", False)
        if is_unmapped:
            continue

        read_name = _get_attr(aln, "read_name", "") or _get_attr(aln, "query_name", "")
        chrom = _get_attr(aln, "reference_name", "")
        ref_start = _get_attr(aln, "reference_start", 0)
        cigar = _get_attr(aln, "cigar_string", "")
        cigar_tuples = _get_attr(aln, "cigar_tuples", [])

        # 1. CIGAR-based SV detection
        cigar_svs = _detect_cigar_svs(chrom, ref_start, cigar, cigar_tuples, read_name, min_size)
        raw_calls.extend(cigar_svs)

        # 2. Split-read based SV detection from supplementary alignments
        sa_tag = ""
        tags = _get_attr(aln, "tags", {})
        if isinstance(tags, dict):
            sa_tag = str(tags.get("SA", ""))

        if sa_tag:
            split_svs = _detect_split_read_svs(
                chrom, ref_start, cigar, sa_tag, read_name,
                _get_attr(aln, "is_reverse", False),
                _get_attr(aln, "query_length", 0) or len(_get_attr(aln, "query_sequence", "")),
                min_size,
            )
            raw_calls.extend(split_svs)

    # Cluster and merge overlapping calls
    merged = _cluster_sv_calls(raw_calls, max_distance=500, min_support=min_support)

    # Sort by chromosome and position
    merged.sort(key=lambda sv: (sv.chromosome, sv.start))

    logger.info("Detected %d structural variants from %d alignments", len(merged), len(alignments))
    return merged


def _detect_cigar_svs(
    chrom: str,
    ref_start: int,
    cigar_string: str,
    cigar_tuples: list[tuple[int, int]],
    read_name: str,
    min_size: int,
) -> list[StructuralVariant]:
    """Detect SVs from large CIGAR operations (insertions and deletions)."""
    svs: list[StructuralVariant] = []

    if cigar_tuples:
        tuples = cigar_tuples
    elif cigar_string:
        tuples = _parse_cigar(cigar_string)
    else:
        return svs

    ref_pos = ref_start

    for op, length in tuples:
        if op == 2 and length >= min_size:  # D = deletion
            svs.append(StructuralVariant(
                sv_type="DEL",
                chromosome=chrom,
                start=ref_pos,
                end=ref_pos + length,
                size=length,
                supporting_reads=1,
                read_names=[read_name],
                quality=float(length),
                info={"source": "cigar"},
            ))
        elif op == 1 and length >= min_size:  # I = insertion
            svs.append(StructuralVariant(
                sv_type="INS",
                chromosome=chrom,
                start=ref_pos,
                end=ref_pos + 1,
                size=length,
                supporting_reads=1,
                read_names=[read_name],
                quality=float(length),
                info={"source": "cigar"},
            ))

        # Advance reference position
        if op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume reference
            ref_pos += length

    return svs


def _detect_split_read_svs(
    primary_chrom: str,
    primary_start: int,
    primary_cigar: str,
    sa_tag: str,
    read_name: str,
    is_reverse: bool,
    query_length: int,
    min_size: int,
) -> list[StructuralVariant]:
    """Detect SVs from supplementary alignments (split reads)."""
    svs: list[StructuralVariant] = []

    # Parse SA tag: chr,pos,strand,CIGAR,mapQ,NM;
    sa_entries = [e.strip() for e in sa_tag.split(";") if e.strip()]

    primary_ref_end = primary_start + _cigar_ref_length(primary_cigar)
    primary_strand = "-" if is_reverse else "+"

    for entry in sa_entries:
        parts = entry.split(",")
        if len(parts) < 6:
            continue

        sa_chrom = parts[0]
        try:
            sa_pos = int(parts[1]) - 1  # Convert to 0-based
        except ValueError:
            continue
        sa_strand = parts[2]
        sa_cigar = parts[3]
        try:
            sa_mapq = int(parts[4])
        except ValueError:
            sa_mapq = 0

        sa_ref_end = sa_pos + _cigar_ref_length(sa_cigar)

        # Classify SV type based on alignment signatures
        if sa_chrom == primary_chrom:
            # Same chromosome
            if sa_strand == primary_strand:
                # Same strand: deletion, insertion, or duplication
                gap = 0
                if sa_pos > primary_ref_end:
                    gap = sa_pos - primary_ref_end
                elif primary_start > sa_ref_end:
                    gap = primary_start - sa_ref_end

                if gap >= min_size:
                    # Large gap between alignments = deletion
                    del_start = min(primary_ref_end, sa_ref_end)
                    del_end = max(primary_start, sa_pos)
                    if del_end > del_start:
                        svs.append(StructuralVariant(
                            sv_type="DEL",
                            chromosome=primary_chrom,
                            start=del_start,
                            end=del_end,
                            size=del_end - del_start,
                            supporting_reads=1,
                            read_names=[read_name],
                            quality=float(sa_mapq),
                            info={"source": "split_read"},
                        ))

                # Check for tandem duplication (overlapping alignments)
                overlap = 0
                if sa_pos < primary_ref_end and sa_ref_end > primary_start:
                    overlap = min(primary_ref_end, sa_ref_end) - max(primary_start, sa_pos)

                if overlap >= min_size:
                    svs.append(StructuralVariant(
                        sv_type="DUP",
                        chromosome=primary_chrom,
                        start=max(primary_start, sa_pos),
                        end=min(primary_ref_end, sa_ref_end),
                        size=overlap,
                        supporting_reads=1,
                        read_names=[read_name],
                        quality=float(sa_mapq),
                        info={"source": "split_read"},
                    ))

            else:
                # Different strands on same chromosome = inversion
                inv_start = min(primary_start, sa_pos)
                inv_end = max(primary_ref_end, sa_ref_end)
                inv_size = inv_end - inv_start

                if inv_size >= min_size:
                    svs.append(StructuralVariant(
                        sv_type="INV",
                        chromosome=primary_chrom,
                        start=inv_start,
                        end=inv_end,
                        size=inv_size,
                        supporting_reads=1,
                        read_names=[read_name],
                        quality=float(sa_mapq),
                        info={"source": "split_read"},
                    ))
        else:
            # Different chromosomes = breakend / translocation
            svs.append(StructuralVariant(
                sv_type="BND",
                chromosome=primary_chrom,
                start=primary_ref_end,
                end=primary_ref_end + 1,
                size=0,
                supporting_reads=1,
                read_names=[read_name],
                quality=float(sa_mapq),
                mate_chromosome=sa_chrom,
                mate_position=sa_pos,
                info={"source": "split_read", "sa_strand": sa_strand},
            ))

    return svs


def detect_insertions(
    alignments: Sequence[Any],
    min_size: int = 50,
    min_support: int = 2,
) -> list[StructuralVariant]:
    """Detect insertions from long-read alignments.

    Specifically targets insertion events from CIGAR I operations and
    unaligned portions of split reads.

    Args:
        alignments: Sequence of alignment objects.
        min_size: Minimum insertion size.
        min_support: Minimum supporting read count.

    Returns:
        List of insertion StructuralVariant objects.
    """
    all_svs = detect_sv_from_long_reads(alignments, min_size=min_size, min_support=min_support)
    insertions = [sv for sv in all_svs if sv.sv_type == "INS"]
    logger.info("Detected %d insertions >= %d bp", len(insertions), min_size)
    return insertions


def detect_inversions(
    alignments: Sequence[Any],
    min_size: int = 50,
    min_support: int = 2,
) -> list[StructuralVariant]:
    """Detect inversions from supplementary alignments with strand switches.

    An inversion is identified when a read has supplementary alignments on
    the same chromosome but opposite strands, indicating a segment of the
    reference that is inverted in the sample.

    Args:
        alignments: Sequence of alignment objects.
        min_size: Minimum inversion size.
        min_support: Minimum supporting read count.

    Returns:
        List of inversion StructuralVariant objects.
    """
    all_svs = detect_sv_from_long_reads(alignments, min_size=min_size, min_support=min_support)
    inversions = [sv for sv in all_svs if sv.sv_type == "INV"]
    logger.info("Detected %d inversions >= %d bp", len(inversions), min_size)
    return inversions


def phase_structural_variants(
    variants: Sequence[StructuralVariant],
    haplotype_tags: dict[str, int],
) -> list[StructuralVariant]:
    """Phase structural variants to haplotypes using read-level haplotype assignments.

    Uses haplotype tag information (HP tag from haplotagged BAMs) to assign
    SVs to haplotype 1 or haplotype 2 based on the haplotype of their
    supporting reads.

    Args:
        variants: Sequence of StructuralVariant objects.
        haplotype_tags: Dictionary mapping read_name -> haplotype (1 or 2).

    Returns:
        List of StructuralVariant objects with haplotype field updated.
    """
    phased: list[StructuralVariant] = []

    for sv in variants:
        hp_counts: dict[int, int] = {1: 0, 2: 0}

        for read_name in sv.read_names:
            hp = haplotype_tags.get(read_name, 0)
            if hp in (1, 2):
                hp_counts[hp] += 1

        total_phased = hp_counts[1] + hp_counts[2]

        # Assign haplotype based on majority vote
        if total_phased == 0:
            assigned_hp = 0  # Unphased
        elif hp_counts[1] > hp_counts[2]:
            assigned_hp = 1
        elif hp_counts[2] > hp_counts[1]:
            assigned_hp = 2
        else:
            assigned_hp = 0  # Ambiguous

        # Determine genotype
        if assigned_hp == 0:
            genotype = "./."
        else:
            # If all supporting reads are on one haplotype, likely heterozygous
            other_hp = 2 if assigned_hp == 1 else 1
            if hp_counts[other_hp] == 0 and total_phased > 0:
                genotype = "0/1"
            elif hp_counts[1] > 0 and hp_counts[2] > 0:
                genotype = "1/1"
            else:
                genotype = "0/1"

        phased_sv = StructuralVariant(
            sv_type=sv.sv_type,
            chromosome=sv.chromosome,
            start=sv.start,
            end=sv.end,
            size=sv.size,
            supporting_reads=sv.supporting_reads,
            read_names=list(sv.read_names),
            quality=sv.quality,
            genotype=genotype,
            haplotype=assigned_hp,
            inserted_sequence=sv.inserted_sequence,
            mate_chromosome=sv.mate_chromosome,
            mate_position=sv.mate_position,
            info={
                **sv.info,
                "hp1_support": hp_counts[1],
                "hp2_support": hp_counts[2],
                "phase_confidence": max(hp_counts.values()) / total_phased if total_phased > 0 else 0.0,
            },
        )
        phased.append(phased_sv)

    phased_count = sum(1 for sv in phased if sv.haplotype > 0)
    logger.info("Phased %d/%d structural variants", phased_count, len(phased))

    return phased


# --- Internal helper functions ---


def _get_attr(obj: Any, attr: str, default: Any = None) -> Any:
    """Get attribute from an object or dictionary."""
    if isinstance(obj, dict):
        return obj.get(attr, default)
    return getattr(obj, attr, default)


def _parse_cigar(cigar_string: str) -> list[tuple[int, int]]:
    """Parse a CIGAR string into a list of (operation, length) tuples.

    Operations: M=0, I=1, D=2, N=3, S=4, H=5, P=6, =7, X=8
    """
    ops = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}
    tuples = []
    for match in re.finditer(r"(\d+)([MIDNSHP=X])", cigar_string):
        length = int(match.group(1))
        op = ops.get(match.group(2), 0)
        tuples.append((op, length))
    return tuples


def _cigar_ref_length(cigar_string: str) -> int:
    """Calculate reference-consuming length from a CIGAR string."""
    ref_consuming = {"M", "D", "N", "=", "X"}
    total = 0
    for match in re.finditer(r"(\d+)([MIDNSHP=X])", cigar_string):
        length = int(match.group(1))
        op = match.group(2)
        if op in ref_consuming:
            total += length
    return total


def _cluster_sv_calls(
    calls: list[StructuralVariant],
    max_distance: int = 500,
    min_support: int = 2,
) -> list[StructuralVariant]:
    """Cluster overlapping SV calls and merge them.

    Groups nearby SV calls of the same type and merges them into
    consensus calls with combined support.

    Args:
        calls: List of raw SV calls.
        max_distance: Maximum distance between breakpoints to cluster.
        min_support: Minimum read support after merging.

    Returns:
        List of merged SV calls.
    """
    if not calls:
        return []

    # Group by (sv_type, chromosome)
    groups: dict[tuple[str, str], list[StructuralVariant]] = defaultdict(list)
    for call in calls:
        groups[(call.sv_type, call.chromosome)].append(call)

    merged: list[StructuralVariant] = []

    for (sv_type, chrom), group_calls in groups.items():
        # Sort by start position
        group_calls.sort(key=lambda sv: sv.start)

        clusters: list[list[StructuralVariant]] = []
        current_cluster: list[StructuralVariant] = [group_calls[0]]

        for i in range(1, len(group_calls)):
            call = group_calls[i]
            cluster_end = max(sv.end for sv in current_cluster)
            cluster_start = min(sv.start for sv in current_cluster)

            # Check if this call overlaps or is near the current cluster
            if (abs(call.start - cluster_start) <= max_distance or
                    abs(call.start - cluster_end) <= max_distance):
                current_cluster.append(call)
            else:
                clusters.append(current_cluster)
                current_cluster = [call]

        clusters.append(current_cluster)

        # Merge each cluster into a single call
        for cluster in clusters:
            all_reads = []
            total_support = 0
            for sv in cluster:
                all_reads.extend(sv.read_names)
                total_support += sv.supporting_reads

            # De-duplicate read names
            unique_reads = list(set(all_reads))
            actual_support = len(unique_reads)

            if actual_support < min_support:
                continue

            # Use median breakpoints
            starts = sorted([sv.start for sv in cluster])
            ends = sorted([sv.end for sv in cluster])
            sizes = sorted([sv.size for sv in cluster])

            median_start = starts[len(starts) // 2]
            median_end = ends[len(ends) // 2]
            median_size = sizes[len(sizes) // 2]

            max_quality = max(sv.quality for sv in cluster)

            merged.append(StructuralVariant(
                sv_type=sv_type,
                chromosome=chrom,
                start=median_start,
                end=median_end,
                size=median_size,
                supporting_reads=actual_support,
                read_names=unique_reads,
                quality=max_quality,
                genotype="0/1" if actual_support >= min_support else "./.",
                info={
                    "cluster_size": len(cluster),
                    "sources": list(set(
                        sv.info.get("source", "unknown") for sv in cluster
                    )),
                },
            ))

    return merged
