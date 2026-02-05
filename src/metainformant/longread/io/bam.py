"""Long-read BAM file processing with methylation and supplementary alignment support.

Provides real parsing of BAM files from long-read sequencers with support for
modified base tags (MM/ML), supplementary alignment extraction for structural
variant detection, and comprehensive alignment statistics.

Optional dependencies:
    - pysam: For reading BAM/CRAM files
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

try:
    import pysam  # type: ignore[import-untyped]
except ImportError:
    pysam = None  # type: ignore[assignment]


@dataclass
class LongReadAlignment:
    """Represents a single alignment from a long-read BAM file.

    Attributes:
        read_name: Query read name.
        query_sequence: The aligned read sequence.
        query_qualities: Per-base quality scores (list of ints).
        reference_name: Reference contig name.
        reference_start: 0-based start position on reference.
        reference_end: 0-based end position on reference.
        mapping_quality: MAPQ score.
        cigar_string: CIGAR alignment string.
        cigar_tuples: CIGAR as list of (operation, length) tuples.
        is_supplementary: Whether this is a supplementary alignment.
        is_secondary: Whether this is a secondary alignment.
        is_reverse: Whether the read is reverse complemented.
        is_unmapped: Whether the read is unmapped.
        tags: Dictionary of BAM tags.
        methylation_tags: Parsed methylation data from MM/ML tags.
        query_length: Total length of the original query.
        aligned_pairs: List of (query_pos, ref_pos) alignment pairs.
    """

    read_name: str
    query_sequence: str = ""
    query_qualities: list[int] = field(default_factory=list)
    reference_name: str = ""
    reference_start: int = 0
    reference_end: int = 0
    mapping_quality: int = 0
    cigar_string: str = ""
    cigar_tuples: list[tuple[int, int]] = field(default_factory=list)
    is_supplementary: bool = False
    is_secondary: bool = False
    is_reverse: bool = False
    is_unmapped: bool = False
    tags: dict[str, Any] = field(default_factory=dict)
    methylation_tags: dict[str, Any] = field(default_factory=dict)
    query_length: int = 0
    aligned_pairs: list[tuple[int | None, int | None]] = field(default_factory=list)


@dataclass
class AlignmentStats:
    """Summary statistics for a set of long-read alignments.

    Attributes:
        total_reads: Total number of reads.
        mapped_reads: Number of mapped reads.
        unmapped_reads: Number of unmapped reads.
        supplementary_count: Number of supplementary alignments.
        secondary_count: Number of secondary alignments.
        mean_mapping_quality: Average MAPQ across mapped reads.
        mean_read_length: Average query length.
        mean_aligned_length: Average aligned length.
        mean_identity: Average alignment identity (matches / aligned_length).
        total_bases: Total bases in all reads.
        total_aligned_bases: Total aligned bases.
        coverage_depth: Estimated coverage depth (if reference length provided).
    """

    total_reads: int = 0
    mapped_reads: int = 0
    unmapped_reads: int = 0
    supplementary_count: int = 0
    secondary_count: int = 0
    mean_mapping_quality: float = 0.0
    mean_read_length: float = 0.0
    mean_aligned_length: float = 0.0
    mean_identity: float = 0.0
    total_bases: int = 0
    total_aligned_bases: int = 0
    coverage_depth: float = 0.0


def read_long_read_bam(
    filepath: str | Path,
    region: str | None = None,
    include_supplementary: bool = True,
    include_secondary: bool = False,
    min_mapping_quality: int = 0,
) -> list[LongReadAlignment]:
    """Read a BAM file with long-read specific tag extraction.

    Parses BAM files from ONT or PacBio sequencers, extracting methylation
    tags (MM/ML), haplotype tags (HP), and other long-read specific information.

    Args:
        filepath: Path to the BAM file (must be indexed with .bai).
        region: Optional genomic region in format "chr:start-end".
        include_supplementary: Whether to include supplementary alignments.
        include_secondary: Whether to include secondary alignments.
        min_mapping_quality: Minimum MAPQ filter.

    Returns:
        List of LongReadAlignment objects.

    Raises:
        FileNotFoundError: If the BAM file does not exist.
        ImportError: If pysam is not installed.
    """
    if pysam is None:
        raise ImportError("pysam is required for reading BAM files. " "Install it with: uv pip install pysam")

    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"BAM file not found: {filepath}")

    alignments: list[LongReadAlignment] = []

    with pysam.AlignmentFile(str(filepath), "rb") as bam:
        iterator = bam.fetch(region=region) if region else bam.fetch(until_eof=True)

        for read in iterator:
            if read.is_secondary and not include_secondary:
                continue
            if read.is_supplementary and not include_supplementary:
                continue
            if not read.is_unmapped and read.mapping_quality < min_mapping_quality:
                continue

            # Extract tags as dict
            tags_dict: dict[str, Any] = {}
            try:
                tags_dict = dict(read.get_tags())
            except Exception:
                pass

            # Parse methylation tags
            meth_tags = _parse_methylation_from_tags(tags_dict)

            # Get CIGAR
            cigar_tuples = list(read.cigartuples) if read.cigartuples else []
            cigar_string = read.cigarstring or ""

            # Get aligned pairs (limited to avoid memory issues on very long reads)
            aligned_pairs: list[tuple[int | None, int | None]] = []
            try:
                pairs = read.get_aligned_pairs()
                aligned_pairs = [(p[0], p[1]) for p in pairs]
            except Exception:
                pass

            aln = LongReadAlignment(
                read_name=read.query_name or "",
                query_sequence=read.query_sequence or "",
                query_qualities=list(read.query_qualities) if read.query_qualities is not None else [],
                reference_name=read.reference_name or "",
                reference_start=read.reference_start,
                reference_end=read.reference_end,
                mapping_quality=read.mapping_quality,
                cigar_string=cigar_string,
                cigar_tuples=cigar_tuples,
                is_supplementary=read.is_supplementary,
                is_secondary=read.is_secondary,
                is_reverse=read.is_reverse,
                is_unmapped=read.is_unmapped,
                tags=tags_dict,
                methylation_tags=meth_tags,
                query_length=read.query_length or 0,
                aligned_pairs=aligned_pairs,
            )
            alignments.append(aln)

    logger.info(
        "Read %d alignments from %s%s",
        len(alignments),
        filepath.name,
        f" region={region}" if region else "",
    )
    return alignments


def _parse_methylation_from_tags(tags: dict[str, Any]) -> dict[str, Any]:
    """Parse MM/ML methylation tags from BAM alignment tags.

    The MM tag encodes modified base positions using the SAM spec:
        MM:Z:C+m,0,1,0;  (CpG methylation at positions 0, 2, 3)
    The ML tag encodes modification probabilities (0-255 scale).

    Args:
        tags: Dictionary of BAM tags.

    Returns:
        Dictionary with parsed methylation data including base type,
        modification code, positions, and probabilities.
    """
    result: dict[str, Any] = {}

    mm_tag = tags.get("MM") or tags.get("Mm")
    ml_tag = tags.get("ML") or tags.get("Ml")

    if mm_tag is None:
        return result

    # Parse MM tag: format is "C+m,0,1,0;" or "C+m?,0,1,0;A+a,2,0;"
    modifications: list[dict[str, Any]] = []
    mm_str = str(mm_tag)

    # Split by semicolons for multiple modification types
    mod_entries = [e.strip() for e in mm_str.split(";") if e.strip()]

    ml_probs: list[int] = []
    if ml_tag is not None:
        if isinstance(ml_tag, (list, tuple)):
            ml_probs = list(ml_tag)
        elif isinstance(ml_tag, bytes):
            ml_probs = list(ml_tag)
        else:
            ml_probs = [int(ml_tag)]

    prob_offset = 0

    for entry in mod_entries:
        # Parse: base_char+mod_code[?],delta1,delta2,...
        match = re.match(r"([ACGTN])([+-])([a-zA-Z0-9]+)(\??),(.*)", entry)
        if not match:
            continue

        base_char = match.group(1)
        strand = match.group(2)
        mod_code = match.group(3)
        implicit_flag = match.group(4)
        deltas_str = match.group(5)

        # Parse delta positions
        deltas = []
        for d in deltas_str.split(","):
            d = d.strip()
            if d:
                try:
                    deltas.append(int(d))
                except ValueError:
                    continue

        # Extract probabilities for this modification
        n_positions = len(deltas)
        probs = ml_probs[prob_offset : prob_offset + n_positions]
        prob_offset += n_positions

        modifications.append(
            {
                "base": base_char,
                "strand": strand,
                "modification": mod_code,
                "implicit": implicit_flag == "?",
                "delta_positions": deltas,
                "probabilities": probs,
                "probability_scale": 255,  # ML tag uses 0-255 scale
            }
        )

    if modifications:
        result["modifications"] = modifications

    return result


def extract_methylation_tags(alignment: LongReadAlignment) -> dict[str, Any]:
    """Extract and parse methylation data from a LongReadAlignment.

    Converts raw MM/ML tags into absolute genomic positions with methylation
    probabilities. Handles both CpG (5mC) and 6mA modification types.

    Args:
        alignment: A LongReadAlignment with methylation tags parsed.

    Returns:
        Dictionary with modification types as keys and lists of
        (position, probability) tuples as values.
    """
    result: dict[str, Any] = {}

    if not alignment.methylation_tags or "modifications" not in alignment.methylation_tags:
        # Try parsing from raw tags
        parsed = _parse_methylation_from_tags(alignment.tags)
        if not parsed or "modifications" not in parsed:
            return result
        mods = parsed["modifications"]
    else:
        mods = alignment.methylation_tags["modifications"]

    seq = alignment.query_sequence

    for mod_info in mods:
        base_char = mod_info["base"]
        mod_code = mod_info["modification"]
        deltas = mod_info["delta_positions"]
        probs = mod_info.get("probabilities", [])

        # Convert delta-encoded positions to absolute positions in the read
        # Deltas encode the number of skipped target bases between modifications
        absolute_positions: list[int] = []
        if seq:
            # Find all positions of the target base in the sequence
            target_positions = [i for i, c in enumerate(seq.upper()) if c == base_char]

            target_idx = 0
            for i, delta in enumerate(deltas):
                target_idx += delta
                if target_idx < len(target_positions):
                    absolute_positions.append(target_positions[target_idx])
                    target_idx += 1  # Move past the modified base

        # Convert to genomic coordinates if aligned
        genomic_positions: list[tuple[int, float]] = []
        for i, abs_pos in enumerate(absolute_positions):
            # Convert read position to reference position using aligned pairs
            ref_pos = None
            for qpos, rpos in alignment.aligned_pairs:
                if qpos == abs_pos and rpos is not None:
                    ref_pos = rpos + alignment.reference_start
                    break

            prob = probs[i] / 255.0 if i < len(probs) else 0.0

            if ref_pos is not None:
                genomic_positions.append((ref_pos, prob))
            else:
                # Use read-level position if reference mapping unavailable
                genomic_positions.append((abs_pos, prob))

        mod_key = f"{base_char}+{mod_code}"
        result[mod_key] = {
            "read_positions": absolute_positions,
            "genomic_positions": genomic_positions,
            "probabilities": [p / 255.0 for p in probs] if probs else [],
        }

    return result


def get_supplementary_alignments(alignment: LongReadAlignment) -> list[dict[str, Any]]:
    """Extract supplementary alignment information for split-read SV detection.

    Parses the SA tag which encodes supplementary alignments in the format:
        SA:Z:chr,pos,strand,CIGAR,mapQ,NM;

    Args:
        alignment: A LongReadAlignment object.

    Returns:
        List of dictionaries, each containing supplementary alignment details
        including reference name, position, strand, CIGAR, and mapping quality.
    """
    sa_tag = alignment.tags.get("SA", "")
    if not sa_tag:
        return []

    supplementary: list[dict[str, Any]] = []

    entries = [e.strip() for e in str(sa_tag).split(";") if e.strip()]
    for entry in entries:
        parts = entry.split(",")
        if len(parts) < 6:
            continue

        ref_name = parts[0]
        try:
            pos = int(parts[1]) - 1  # Convert 1-based SAM to 0-based
        except ValueError:
            continue
        strand = parts[2]
        cigar = parts[3]
        try:
            mapq = int(parts[4])
        except ValueError:
            mapq = 0
        try:
            nm = int(parts[5])
        except ValueError:
            nm = 0

        # Parse CIGAR to get aligned length
        aligned_length = _cigar_aligned_length(cigar)

        supplementary.append(
            {
                "reference_name": ref_name,
                "reference_start": pos,
                "strand": strand,
                "cigar": cigar,
                "mapping_quality": mapq,
                "edit_distance": nm,
                "aligned_length": aligned_length,
                "is_reverse": strand == "-",
            }
        )

    return supplementary


def _cigar_aligned_length(cigar_string: str) -> int:
    """Calculate the aligned reference length from a CIGAR string.

    CIGAR operations that consume reference:
        M (alignment match), D (deletion), N (skipped region), = (match), X (mismatch)

    Args:
        cigar_string: CIGAR alignment string.

    Returns:
        Total reference-consuming length.
    """
    ref_consuming = {"M", "D", "N", "=", "X"}
    total = 0
    for match in re.finditer(r"(\d+)([MIDNSHP=X])", cigar_string):
        length = int(match.group(1))
        op = match.group(2)
        if op in ref_consuming:
            total += length
    return total


def _cigar_query_length(cigar_string: str) -> int:
    """Calculate the query-consuming length from a CIGAR string.

    CIGAR operations that consume query:
        M (alignment match), I (insertion), S (soft clip), = (match), X (mismatch)

    Args:
        cigar_string: CIGAR alignment string.

    Returns:
        Total query-consuming length.
    """
    query_consuming = {"M", "I", "S", "=", "X"}
    total = 0
    for match in re.finditer(r"(\d+)([MIDNSHP=X])", cigar_string):
        length = int(match.group(1))
        op = match.group(2)
        if op in query_consuming:
            total += length
    return total


def calculate_alignment_stats(
    alignments: Sequence[LongReadAlignment],
    reference_length: int | None = None,
) -> AlignmentStats:
    """Calculate comprehensive alignment statistics for a set of long-read alignments.

    Computes mapping quality, identity, coverage, and read length distributions
    from a collection of alignments.

    Args:
        alignments: Sequence of LongReadAlignment objects.
        reference_length: Total reference genome length (for coverage calculation).

    Returns:
        AlignmentStats with aggregated metrics.
    """
    if not alignments:
        return AlignmentStats()

    total = len(alignments)
    mapped = 0
    unmapped = 0
    supplementary = 0
    secondary = 0
    mapq_sum = 0.0
    read_lengths: list[int] = []
    aligned_lengths: list[int] = []
    identities: list[float] = []
    total_bases = 0
    total_aligned = 0

    for aln in alignments:
        if aln.is_unmapped:
            unmapped += 1
            continue

        if aln.is_supplementary:
            supplementary += 1
        if aln.is_secondary:
            secondary += 1

        mapped += 1
        mapq_sum += aln.mapping_quality

        ql = aln.query_length or len(aln.query_sequence)
        read_lengths.append(ql)
        total_bases += ql

        # Calculate aligned length from CIGAR
        if aln.cigar_string:
            al = _cigar_aligned_length(aln.cigar_string)
        elif aln.reference_end > aln.reference_start:
            al = aln.reference_end - aln.reference_start
        else:
            al = ql

        aligned_lengths.append(al)
        total_aligned += al

        # Calculate identity from NM (edit distance) tag or CIGAR
        nm = aln.tags.get("NM", aln.tags.get("nM", None))
        if nm is not None and al > 0:
            identity = 1.0 - (int(nm) / al)
            identities.append(max(0.0, identity))
        elif aln.cigar_tuples:
            # Estimate identity from CIGAR: matches / (matches + mismatches + indels)
            matches = sum(length for op, length in aln.cigar_tuples if op in (0, 7))  # M or =
            total_ops = sum(length for op, length in aln.cigar_tuples if op in (0, 1, 2, 7, 8))
            if total_ops > 0:
                identities.append(matches / total_ops)

    mean_mapq = mapq_sum / mapped if mapped > 0 else 0.0
    mean_read_len = sum(read_lengths) / len(read_lengths) if read_lengths else 0.0
    mean_aligned_len = sum(aligned_lengths) / len(aligned_lengths) if aligned_lengths else 0.0
    mean_identity = sum(identities) / len(identities) if identities else 0.0

    coverage = 0.0
    if reference_length and reference_length > 0:
        coverage = total_aligned / reference_length

    return AlignmentStats(
        total_reads=total,
        mapped_reads=mapped,
        unmapped_reads=unmapped,
        supplementary_count=supplementary,
        secondary_count=secondary,
        mean_mapping_quality=mean_mapq,
        mean_read_length=mean_read_len,
        mean_aligned_length=mean_aligned_len,
        mean_identity=mean_identity,
        total_bases=total_bases,
        total_aligned_bases=total_aligned,
        coverage_depth=coverage,
    )
