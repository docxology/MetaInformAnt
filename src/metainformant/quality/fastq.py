"""FASTQ file quality analysis and processing.

This module provides comprehensive quality control for FASTQ sequencing files,
including per-base quality scores, sequence length distributions, adapter detection,
and various quality metrics.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Any, Iterator, Tuple, Optional
from collections import defaultdict, Counter
import gzip
import statistics

from metainformant.core import logging, errors, validation, io

logger = logging.get_logger(__name__)


class FastqRecord:
    """Represents a single FASTQ record."""

    def __init__(self, header: str, sequence: str, quality_header: str, quality: str):
        """Initialize a FASTQ record.

        Args:
            header: Read header line (starts with @)
            sequence: DNA sequence
            quality_header: Quality header line (starts with +)
            quality: Quality scores as ASCII characters
        """
        self.header = header
        self.sequence = sequence
        self.quality_header = quality_header
        self.quality = quality

        # Validate record structure
        if not header.startswith('@'):
            raise ValueError(f"Invalid FASTQ header: {header}")
        if not quality_header.startswith('+'):
            raise ValueError(f"Invalid FASTQ quality header: {quality_header}")
        if len(sequence) != len(quality):
            raise ValueError(f"Sequence and quality lengths don't match: {len(sequence)} != {len(quality)}")

    @property
    def name(self) -> str:
        """Get the read name from the header."""
        return self.header[1:].split()[0] if self.header.startswith('@') else ""

    @property
    def length(self) -> int:
        """Get the sequence length."""
        return len(self.sequence)

    def quality_scores(self) -> List[int]:
        """Convert ASCII quality scores to numeric values."""
        return [ord(c) - 33 for c in self.quality]  # Illumina 1.8+ encoding

    def mean_quality(self) -> float:
        """Calculate mean quality score."""
        return statistics.mean(self.quality_scores())

    def gc_content(self) -> float:
        """Calculate GC content percentage."""
        if not self.sequence:
            return 0.0
        gc_count = self.sequence.upper().count('G') + self.sequence.upper().count('C')
        return (gc_count / len(self.sequence)) * 100.0


def read_fastq_records(path: str | Path, max_records: int | None = None) -> Iterator[FastqRecord]:
    """Read FASTQ records from a file.

    Args:
        path: Path to FASTQ file (can be gzipped)
        max_records: Maximum number of records to read (None for all)

    Yields:
        FastqRecord objects

    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file format is invalid
    """
    path = validation.validate_path_exists(Path(path))

    logger.info(f"Reading FASTQ records from {path}")

    try:
        with io.open_text_auto(path, 'rt') as f:
            record_count = 0
            while max_records is None or record_count < max_records:
                # Read 4 lines for each record
                lines = []
                for i in range(4):
                    line = f.readline()
                    if not line:
                        if lines:  # Incomplete record
                            raise ValueError(f"Incomplete FASTQ record at end of file")
                        break  # End of file
                    lines.append(line.rstrip('\n\r'))

                if len(lines) == 0:
                    break  # End of file
                elif len(lines) != 4:
                    raise ValueError(f"Incomplete FASTQ record: got {len(lines)} lines, expected 4")

                record = FastqRecord(*lines)
                yield record
                record_count += 1

    except Exception as e:
        logger.error(f"Error reading FASTQ file {path}: {e}")
        raise errors.FileIOError(f"Failed to read FASTQ file: {e}") from e


def analyze_fastq_quality(fastq_path: str | Path, n_reads: int | None = None) -> Dict[str, Any]:
    """Perform comprehensive quality analysis on a FASTQ file.

    Args:
        fastq_path: Path to FASTQ file
        n_reads: Number of reads to analyze (None for all)

    Returns:
        Dictionary containing quality metrics
    """
    validation.validate_path_exists(Path(fastq_path))

    logger.info(f"Analyzing FASTQ quality for {fastq_path}")

    records = list(read_fastq_records(fastq_path, n_reads))
    if not records:
        raise ValueError("No reads found in FASTQ file")

    return {
        "basic_statistics": basic_statistics(records),
        "per_base_quality": per_base_quality(records),
        "per_sequence_quality": per_sequence_quality(records),
        "sequence_length_distribution": sequence_length_distribution(records),
        "gc_content_distribution": gc_content_distribution(records),
        "adapter_content": adapter_content(records),
        "overrepresented_sequences": overrepresented_sequences(records),
        "duplication_levels": duplication_levels(records),
        "n_content_per_position": n_content_per_position(records),
        "quality_score_distribution": quality_score_distribution(records),
    }


def basic_statistics(records: List[FastqRecord]) -> Dict[str, Any]:
    """Calculate basic statistics for FASTQ records.

    Args:
        records: List of FastqRecord objects

    Returns:
        Dictionary of basic statistics
    """
    if not records:
        return {}

    lengths = [r.length for r in records]
    qualities = [r.mean_quality() for r in records]
    gc_contents = [r.gc_content() for r in records]

    return {
        "total_reads": len(records),
        "total_bases": sum(lengths),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": statistics.mean(lengths),
        "median_length": statistics.median(lengths),
        "length_std": statistics.stdev(lengths) if len(lengths) > 1 else 0,
        "min_quality": min(qualities),
        "max_quality": max(qualities),
        "mean_quality": statistics.mean(qualities),
        "median_quality": statistics.median(qualities),
        "quality_std": statistics.stdev(qualities) if len(qualities) > 1 else 0,
        "min_gc": min(gc_contents),
        "max_gc": max(gc_contents),
        "mean_gc": statistics.mean(gc_contents),
        "median_gc": statistics.median(gc_contents),
        "gc_std": statistics.stdev(gc_contents) if len(gc_contents) > 1 else 0,
    }


def per_base_quality(records: List[FastqRecord]) -> Dict[str, Any]:
    """Calculate per-base quality statistics.

    Args:
        records: List of FastqRecord objects

    Returns:
        Dictionary with per-base quality statistics
    """
    if not records:
        return {}

    max_length = max(r.length for r in records)
    quality_matrix = []

    for record in records:
        scores = record.quality_scores()
        # Pad shorter reads with None
        padded = scores + [None] * (max_length - len(scores))
        quality_matrix.append(padded)

    # Calculate statistics for each position
    positions = []
    for pos in range(max_length):
        pos_qualities = [row[pos] for row in quality_matrix if row[pos] is not None]
        if pos_qualities:
            positions.append({
                "position": pos + 1,
                "mean": statistics.mean(pos_qualities),
                "median": statistics.median(pos_qualities),
                "lower_quartile": statistics.quantiles(pos_qualities, n=4)[0],
                "upper_quartile": statistics.quantiles(pos_qualities, n=4)[2],
                "min": min(pos_qualities),
                "max": max(pos_qualities),
                "count": len(pos_qualities),
            })

    return {"positions": positions}


def per_sequence_quality(records: List[FastqRecord]) -> Dict[str, Any]:
    """Calculate per-sequence quality score distribution.

    Args:
        records: List of FastqRecord objects

    Returns:
        Dictionary with quality score distribution
    """
    if not records:
        return {}

    mean_qualities = [r.mean_quality() for r in records]

    # Create histogram bins
    if mean_qualities:
        min_qual = min(mean_qualities)
        max_qual = max(mean_qualities)
        bin_width = 1.0
        bins = []

        for bin_start in range(int(min_qual), int(max_qual) + 2):
            count = sum(1 for q in mean_qualities if bin_start <= q < bin_start + bin_width)
            if count > 0:
                bins.append({
                    "bin_start": bin_start,
                    "bin_end": bin_start + bin_width,
                    "count": count,
                    "percentage": (count / len(mean_qualities)) * 100,
                })

    return {"bins": bins}


def sequence_length_distribution(records: List[FastqRecord]) -> Dict[str, Any]:
    """Calculate sequence length distribution.

    Args:
        records: List of FastqRecord objects

    Returns:
        Dictionary with length distribution
    """
    if not records:
        return {}

    lengths = [r.length for r in records]
    length_counts = Counter(lengths)

    distribution = []
    for length, count in sorted(length_counts.items()):
        distribution.append({
            "length": length,
            "count": count,
            "percentage": (count / len(records)) * 100,
        })

    return {"distribution": distribution}


def gc_content_distribution(records: List[FastqRecord]) -> Dict[str, Any]:
    """Calculate GC content distribution.

    Args:
        records: List of FastqRecord objects

    Returns:
        Dictionary with GC content distribution
    """
    if not records:
        return {}

    gc_contents = [r.gc_content() for r in records]

    # Create histogram bins (0-100% in 5% increments)
    bins = []
    for bin_start in range(0, 101, 5):
        bin_end = min(bin_start + 5, 100)
        count = sum(1 for gc in gc_contents if bin_start <= gc < bin_end)
        if count > 0:
            bins.append({
                "bin_start": bin_start,
                "bin_end": bin_end,
                "count": count,
                "percentage": (count / len(gc_contents)) * 100,
            })

    return {"bins": bins}


def adapter_content(records: List[FastqRecord], adapters: List[str] | None = None) -> Dict[str, Any]:
    """Detect adapter content in sequences.

    Args:
        records: List of FastqRecord objects
        adapters: List of adapter sequences to check (default: common Illumina adapters)

    Returns:
        Dictionary with adapter detection results
    """
    if not records:
        return {}

    # Default Illumina adapters
    if adapters is None:
        adapters = [
            "AGATCGGAAGAG",  # TruSeq Universal Adapter
            "GATCGGAAGAG",   # TruSeq Adapter, Read 1
            "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",  # TruSeq Adapter, Read 2
            "CTGTCTCTTAT",   # Nextera Transposase Sequence
        ]

    results = {}
    for adapter in adapters:
        adapter_name = f"adapter_{len(results) + 1}"
        adapter_length = len(adapter)

        positions = []
        for pos in range(0, 101, 10):  # Check every 10% of read length
            pos_bp = max(1, int((pos / 100) * max(r.length for r in records)))
            count = 0

            for record in records:
                if pos_bp + adapter_length <= record.length:
                    seq_slice = record.sequence[pos_bp:pos_bp + adapter_length]
                    if adapter in seq_slice:
                        count += 1

            positions.append({
                "position": pos,
                "count": count,
                "percentage": (count / len(records)) * 100 if records else 0,
            })

        results[adapter_name] = {
            "adapter_sequence": adapter,
            "positions": positions,
        }

    return {"adapters": results}


def overrepresented_sequences(records: List[FastqRecord], min_length: int = 20) -> Dict[str, Any]:
    """Find overrepresented sequences in the FASTQ file.

    Args:
        records: List of FastqRecord objects
        min_length: Minimum sequence length to consider

    Returns:
        Dictionary with overrepresented sequences
    """
    if not records:
        return {}

    # Extract sequences of minimum length
    sequences = []
    for record in records:
        if record.length >= min_length:
            sequences.append(record.sequence)

    if not sequences:
        return {"overrepresented": []}

    # Count sequence frequencies
    seq_counts = Counter(sequences)
    total_sequences = len(sequences)

    # Find sequences that appear more than expected by chance
    overrepresented = []
    for seq, count in seq_counts.most_common(10):  # Top 10
        percentage = (count / total_sequences) * 100
        if percentage > 1.0:  # More than 1% of reads
            overrepresented.append({
                "sequence": seq,
                "count": count,
                "percentage": percentage,
            })

    return {"overrepresented": overrepresented}


def duplication_levels(records: List[FastqRecord]) -> Dict[str, Any]:
    """Calculate sequence duplication levels.

    Args:
        records: List of FastqRecord objects

    Returns:
        Dictionary with duplication statistics
    """
    if not records:
        return {}

    sequences = [r.sequence for r in records]
    seq_counts = Counter(sequences)

    # Calculate duplication levels
    unique_sequences = len(seq_counts)
    total_sequences = len(sequences)

    # Group by duplication level
    duplication_bins = defaultdict(int)
    for count in seq_counts.values():
        if count <= 10:
            duplication_bins[str(count)] += 1
        else:
            duplication_bins["10+"] += 1

    return {
        "total_sequences": total_sequences,
        "unique_sequences": unique_sequences,
        "duplication_rate": (1 - unique_sequences / total_sequences) * 100,
        "duplication_bins": dict(duplication_bins),
    }


def n_content_per_position(records: List[FastqRecord]) -> Dict[str, Any]:
    """Calculate N content (unknown bases) per position.

    Args:
        records: List of FastqRecord objects

    Returns:
        Dictionary with N content per position
    """
    if not records:
        return {}

    max_length = max(r.length for r in records)
    n_counts = [0] * max_length

    for record in records:
        for i, base in enumerate(record.sequence):
            if base.upper() == 'N':
                if i < max_length:
                    n_counts[i] += 1

    positions = []
    for pos, count in enumerate(n_counts):
        positions.append({
            "position": pos + 1,
            "n_count": count,
            "n_percentage": (count / len(records)) * 100,
        })

    return {"positions": positions}


def quality_score_distribution(records: List[FastqRecord]) -> Dict[str, Any]:
    """Calculate overall quality score distribution.

    Args:
        records: List of FastqRecord objects

    Returns:
        Dictionary with quality score distribution
    """
    if not records:
        return {}

    # Collect all quality scores
    all_scores = []
    for record in records:
        all_scores.extend(record.quality_scores())

    if not all_scores:
        return {"distribution": []}

    # Create distribution
    score_counts = Counter(all_scores)
    distribution = []
    for score in sorted(score_counts.keys()):
        distribution.append({
            "quality_score": score,
            "count": score_counts[score],
            "percentage": (score_counts[score] / len(all_scores)) * 100,
        })

    return {"distribution": distribution}


def filter_reads(fastq_path: str | Path, output_path: str | Path,
                min_quality: float = 20.0, min_length: int | None = None,
                max_n_bases: int = 0) -> Dict[str, Any]:
    """Filter FASTQ reads based on quality criteria.

    Args:
        fastq_path: Input FASTQ file
        output_path: Output FASTQ file for filtered reads
        min_quality: Minimum mean quality score
        min_length: Minimum read length
        max_n_bases: Maximum number of N bases allowed

    Returns:
        Dictionary with filtering statistics
    """
    validation.validate_path_exists(Path(fastq_path))
    output_path = Path(output_path)

    logger.info(f"Filtering FASTQ reads from {fastq_path} to {output_path}")

    total_reads = 0
    passed_reads = 0

    with io.open_text_auto(output_path, 'w') as out_f:
        for record in read_fastq_records(fastq_path):
            total_reads += 1

            # Apply filters
            passes = True

            if record.mean_quality() < min_quality:
                passes = False

            if min_length is not None and record.length < min_length:
                passes = False

            if max_n_bases >= 0:
                n_count = record.sequence.upper().count('N')
                if n_count > max_n_bases:
                    passes = False

            if passes:
                # Write record
                out_f.write(f"{record.header}\n")
                out_f.write(f"{record.sequence}\n")
                out_f.write(f"{record.quality_header}\n")
                out_f.write(f"{record.quality}\n")
                passed_reads += 1

    stats = {
        "total_reads": total_reads,
        "passed_reads": passed_reads,
        "filtered_reads": total_reads - passed_reads,
        "pass_rate": (passed_reads / total_reads * 100) if total_reads > 0 else 0,
        "filters_applied": {
            "min_quality": min_quality,
            "min_length": min_length,
            "max_n_bases": max_n_bases,
        },
    }

    logger.info(f"Filtered {passed_reads}/{total_reads} reads ({stats['pass_rate']:.1f}%)")
    return stats

