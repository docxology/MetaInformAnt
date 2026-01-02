"""FASTQ file processing and quality control for DNA sequences.

This module provides tools for reading, writing, and analyzing FASTQ files,
including quality score assessment, read filtering, and format conversion.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any, Dict, Iterator, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)


def read_fastq(path: str | Path) -> Dict[str, Tuple[str, str]]:
    """Read sequences and quality scores from a FASTQ file.

    Args:
        path: Path to FASTQ file (can be gzipped)

    Returns:
        Dictionary mapping read IDs to (sequence, quality_string) tuples

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid

    Example:
        >>> # Assuming test.fastq exists
        >>> reads = read_fastq("test.fastq")
        >>> isinstance(reads, dict)
        True
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"FASTQ file not found: {path}")

    reads = {}

    # Open file (handle gzip compression)
    opener = gzip.open if path.suffix == '.gz' else open
    mode = 'rt' if path.suffix == '.gz' else 'r'

    with opener(path, mode) as f:
        while True:
            # Read four lines: header, sequence, +, quality
            header_line = f.readline().strip()
            if not header_line:
                break  # End of file

            if not header_line.startswith('@'):
                raise ValueError(f"Invalid FASTQ format: expected '@' at line start, got '{header_line[:20]}...'")

            seq_line = f.readline().strip()
            plus_line = f.readline().strip()
            qual_line = f.readline().strip()

            if not seq_line or not qual_line:
                raise ValueError("Incomplete FASTQ record")

            if len(seq_line) != len(qual_line):
                raise ValueError(f"Sequence and quality string lengths don't match: {len(seq_line)} vs {len(qual_line)}")

            # Extract read ID (remove @ and take first word)
            read_id = header_line[1:].split()[0]

            reads[read_id] = (seq_line, qual_line)

    logger.info(f"Read {len(reads)} reads from {path}")
    return reads


def write_fastq(sequences: Dict[str, Tuple[str, str]], path: str | Path) -> None:
    """Write sequences and quality scores to a FASTQ file.

    Args:
        sequences: Dictionary mapping read IDs to (sequence, quality_string) tuples
        path: Output path for FASTQ file

    Example:
        >>> reads = {"read1": ("ATCG", "IIII")}
        >>> write_fastq(reads, "output.fastq")
    """
    path = Path(path)

    # Create parent directories if needed
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, 'w') as f:
        for read_id, (sequence, quality) in sequences.items():
            if len(sequence) != len(quality):
                raise ValueError(f"Sequence and quality lengths don't match for read {read_id}")

            f.write(f"@{read_id}\n")
            f.write(f"{sequence}\n")
            f.write("+\n")
            f.write(f"{quality}\n")

    logger.info(f"Wrote {len(sequences)} reads to {path}")


def assess_quality(fastq_path: str | Path) -> Dict[str, Any]:
    """Assess overall quality of reads in a FASTQ file.

    Args:
        fastq_path: Path to FASTQ file

    Returns:
        Dictionary with quality statistics

    Example:
        >>> # Assuming test.fastq exists
        >>> stats = assess_quality("test.fastq")
        >>> "mean_quality" in stats
        True
    """
    reads = read_fastq(fastq_path)

    if not reads:
        return {
            'total_reads': 0,
            'mean_quality': 0.0,
            'median_quality': 0.0,
            'quality_distribution': {},
            'gc_content': 0.0,
            'read_lengths': []
        }

    total_reads = len(reads)
    all_qualities = []
    gc_counts = 0
    total_bases = 0
    read_lengths = []

    for sequence, quality_string in reads.values():
        # Convert quality string to Phred scores
        qualities = [ord(c) - 33 for c in quality_string]  # Illumina 1.8+ encoding
        all_qualities.extend(qualities)

        # Calculate GC content
        seq_upper = sequence.upper()
        gc_counts += seq_upper.count('G') + seq_upper.count('C')
        total_bases += len(sequence)

        read_lengths.append(len(sequence))

    # Quality statistics
    mean_quality = sum(all_qualities) / len(all_qualities) if all_qualities else 0.0
    median_quality = sorted(all_qualities)[len(all_qualities) // 2] if all_qualities else 0.0

    # Quality distribution (bucketed)
    quality_distribution = {}
    for q in all_qualities:
        bucket = (q // 5) * 5  # Group by 5-point buckets
        quality_distribution[bucket] = quality_distribution.get(bucket, 0) + 1

    # GC content
    gc_content = (gc_counts / total_bases * 100) if total_bases > 0 else 0.0

    return {
        'total_reads': total_reads,
        'mean_quality': round(mean_quality, 2),
        'median_quality': median_quality,
        'quality_distribution': quality_distribution,
        'gc_content': round(gc_content, 2),
        'read_lengths': read_lengths,
        'min_length': min(read_lengths) if read_lengths else 0,
        'max_length': max(read_lengths) if read_lengths else 0,
        'mean_length': sum(read_lengths) / len(read_lengths) if read_lengths else 0.0
    }


def filter_reads(fastq_path: str | Path, min_quality: int = 20) -> Iterator[str]:
    """Filter reads based on minimum quality score.

    Args:
        fastq_path: Path to FASTQ file
        min_quality: Minimum average quality score (default: 20)

    Yields:
        FASTQ records (4 lines each) that pass quality filter

    Example:
        >>> # Assuming test.fastq exists
        >>> filtered = list(filter_reads("test.fastq", min_quality=25))
        >>> len(filtered) % 4 == 0  # FASTQ records come in groups of 4
        True
    """
    # Open file (handle gzip compression)
    path = Path(fastq_path)
    opener = gzip.open if path.suffix == '.gz' else open
    mode = 'rt' if path.suffix == '.gz' else 'r'

    with opener(path, mode) as f:
        while True:
            # Read complete FASTQ record
            lines = []
            for _ in range(4):
                line = f.readline()
                if not line:  # End of file
                    return
                lines.append(line.rstrip('\n'))

            if not lines[0].startswith('@'):
                continue  # Skip invalid records

            sequence = lines[1]
            quality_string = lines[3]

            if len(sequence) != len(quality_string):
                continue  # Skip malformed records

            # Calculate average quality
            qualities = [ord(c) - 33 for c in quality_string]  # Illumina 1.8+ encoding
            avg_quality = sum(qualities) / len(qualities)

            if avg_quality >= min_quality:
                # Yield complete FASTQ record
                yield '\n'.join(lines)


def convert_fastq_to_fasta(fastq_path: str | Path, fasta_path: str | Path) -> None:
    """Convert FASTQ file to FASTA format.

    Args:
        fastq_path: Input FASTQ file path
        fasta_path: Output FASTA file path

    Example:
        >>> # Assuming test.fastq exists
        >>> convert_fastq_to_fasta("test.fastq", "output.fasta")
    """
    reads = read_fastq(fastq_path)

    with open(fasta_path, 'w') as f:
        for read_id, (sequence, _) in reads.items():
            f.write(f">{read_id}\n")
            f.write(f"{sequence}\n")

    logger.info(f"Converted {len(reads)} reads from FASTQ to FASTA")


def trim_reads(fastq_path: str | Path, output_path: str | Path,
               min_length: int = 50, trim_5p: int = 0, trim_3p: int = 0) -> None:
    """Trim reads and filter by minimum length.

    Args:
        fastq_path: Input FASTQ file path
        output_path: Output FASTQ file path
        min_length: Minimum length after trimming (default: 50)
        trim_5p: Bases to trim from 5' end (default: 0)
        trim_3p: Bases to trim from 3' end (default: 0)

    Example:
        >>> # Assuming test.fastq exists
        >>> trim_reads("test.fastq", "trimmed.fastq", min_length=30)
    """
    reads = read_fastq(fastq_path)
    trimmed_reads = {}

    for read_id, (sequence, quality) in reads.items():
        # Trim sequence and quality
        start = trim_5p
        end = len(sequence) - trim_3p

        if end <= start or (end - start) < min_length:
            continue  # Skip reads that become too short

        trimmed_seq = sequence[start:end]
        trimmed_qual = quality[start:end]

        trimmed_reads[read_id] = (trimmed_seq, trimmed_qual)

    write_fastq(trimmed_reads, output_path)
    logger.info(f"Trimmed {len(trimmed_reads)} reads (from {len(reads)} original)")


def gc_content(seq: str) -> float:
    """Calculate GC content of a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        GC content as a fraction (0.0-1.0)

    Example:
        >>> gc_content("ATGC")
        0.5
        >>> gc_content("GGGG")
        1.0
        >>> gc_content("AAAA")
        0.0
    """
    if not seq:
        return 0.0

    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    # Exclude N's and other non-standard bases from total count
    valid_bases = seq.count('A') + seq.count('T') + seq.count('G') + seq.count('C')
    total_count = valid_bases

    return gc_count / total_count if total_count > 0 else 0.0


def average_phred_by_position(fastq_path: str | Path) -> Dict[int, float]:
    """Calculate average Phred quality score by position.

    Args:
        fastq_path: Path to FASTQ file

    Returns:
        Dictionary mapping position to average quality score

    Example:
        >>> # Assuming test.fastq exists
        >>> avg_qual = average_phred_by_position("test.fastq")
        >>> isinstance(avg_qual, dict)
        True
    """
    from collections import defaultdict

    fastq_path = Path(fastq_path)
    if not fastq_path.exists():
        raise FileNotFoundError(f"FASTQ file not found: {fastq_path}")

    position_qualities = defaultdict(list)

    # Open file (handle gzip compression)
    opener = gzip.open if fastq_path.suffix == '.gz' else open
    mode = 'rt' if fastq_path.suffix == '.gz' else 'r'

    with opener(fastq_path, mode) as f:
        line_num = 0
        for line in f:
            line_num += 1
            if line_num % 4 == 0:  # Quality line
                for pos, char in enumerate(line.strip()):
                    # Convert Phred character to quality score
                    quality_score = ord(char) - 33  # Illumina 1.8+ encoding
                    position_qualities[pos].append(quality_score)

    # Calculate averages
    return {pos: sum(qualities) / len(qualities)
            for pos, qualities in position_qualities.items()}


def iter_fastq(fastq_path: str | Path) -> Iterator[Tuple[str, str, str]]:
    """Iterate over FASTQ records yielding (header, sequence, quality).

    Args:
        fastq_path: Path to FASTQ file

    Yields:
        Tuples of (header, sequence, quality_string)

    Example:
        >>> # Assuming test.fastq exists
        >>> records = list(iter_fastq("test.fastq"))
        >>> len(records) >= 0
        True
    """
    fastq_path = Path(fastq_path)
    if not fastq_path.exists():
        raise FileNotFoundError(f"FASTQ file not found: {fastq_path}")

    # Open file (handle gzip compression)
    opener = gzip.open if fastq_path.suffix == '.gz' else open
    mode = 'rt' if fastq_path.suffix == '.gz' else 'r'

    with opener(fastq_path, mode) as f:
        while True:
            # Read four lines: header, sequence, +, quality
            header_line = f.readline().strip()
            if not header_line:  # End of file
                break

            seq_line = f.readline().strip()
            plus_line = f.readline().strip()  # Skip the + line
            qual_line = f.readline().strip()

            if not all([header_line, seq_line, plus_line, qual_line]):
                raise ValueError("Incomplete FASTQ record")

            # Extract read ID (everything after @ up to first space)
            read_id = header_line[1:].split()[0] if header_line.startswith('@') else header_line.split()[0]
            yield read_id, seq_line, qual_line


def calculate_per_base_quality(fastq_path: str | Path) -> Dict[int, Dict[str, float]]:
    """Calculate quality statistics for each base position.

    Args:
        fastq_path: Path to FASTQ file

    Returns:
        Dictionary mapping position to quality statistics

    Example:
        >>> # Assuming test.fastq exists
        >>> stats = calculate_per_base_quality("test.fastq")
        >>> isinstance(stats, dict)
        True
    """
    reads = read_fastq(fastq_path)

    if not reads:
        return {}

    # Assume all reads have same length (typical for Illumina)
    first_read_len = len(next(iter(reads.values()))[0])

    # Initialize statistics for each position
    position_stats = {}
    for pos in range(first_read_len):
        position_stats[pos] = {
            'mean_quality': 0.0,
            'median_quality': 0.0,
            'min_quality': float('inf'),
            'max_quality': 0.0,
            'qualities': []
        }

    # Collect quality scores for each position
    for sequence, quality_string in reads.values():
        qualities = [ord(c) - 33 for c in quality_string]

        for pos, q in enumerate(qualities):
            if pos < first_read_len:  # Safety check
                stats = position_stats[pos]
                stats['qualities'].append(q)
                stats['min_quality'] = min(stats['min_quality'], q)
                stats['max_quality'] = max(stats['max_quality'], q)

    # Calculate summary statistics
    for pos, stats in position_stats.items():
        qualities = stats['qualities']
        if qualities:
            stats['mean_quality'] = sum(qualities) / len(qualities)
            stats['median_quality'] = sorted(qualities)[len(qualities) // 2]
        else:
            stats['min_quality'] = 0.0

        # Remove raw qualities list to save memory
        del stats['qualities']

    return position_stats


def summarize_fastq(fastq_path: str | Path) -> Dict[str, int]:
    """Summarize contents of a FASTQ file.

    Args:
        fastq_path: Path to FASTQ file

    Returns:
        Dictionary with summary statistics
    """
    fastq_path = Path(fastq_path)
    if not fastq_path.exists():
        raise FileNotFoundError(f"FASTQ file not found: {fastq_path}")

    # Open file (handle gzip compression)
    opener = gzip.open if fastq_path.suffix == '.gz' else open
    mode = 'rt' if fastq_path.suffix == '.gz' else 'r'

    summary = {
        'total_reads': 0,
        'total_bases': 0,
        'min_length': float('inf'),
        'max_length': 0,
        'mean_length': 0.0
    }

    lengths = []

    with opener(fastq_path, mode) as f:
        line_num = 0
        for line in f:
            line_num += 1
            if line_num % 4 == 2:  # Sequence line
                seq = line.strip()
                seq_len = len(seq)
                lengths.append(seq_len)
                summary['total_reads'] += 1
                summary['total_bases'] += seq_len
                summary['min_length'] = min(summary['min_length'], seq_len)
                summary['max_length'] = max(summary['max_length'], seq_len)

    if lengths:
        summary['mean_length'] = sum(lengths) / len(lengths)
    else:
        summary['min_length'] = 0

    return summary



