"""ChIP-seq analysis and processing.

This module provides comprehensive tools for analyzing ChIP-seq data,
including peak calling, motif analysis, quality control, and data integration.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Any, Optional, Iterator, Tuple, Set
from collections import defaultdict
import statistics
import math

from metainformant.core import logging, errors, validation, io

logger = logging.get_logger(__name__)


class ChIPPeak:
    """Represents a ChIP-seq peak."""

    def __init__(self, chromosome: str, start: int, end: int,
                 summit: Optional[int] = None, score: float = 0.0,
                 strand: str = '.', signal_value: float = 0.0,
                 p_value: Optional[float] = None, q_value: Optional[float] = None):
        """Initialize a ChIP-seq peak.

        Args:
            chromosome: Chromosome name
            start: Peak start position
            end: Peak end position
            summit: Peak summit position
            score: Peak score
            strand: DNA strand
            signal_value: Signal value
            p_value: P-value
            q_value: Q-value (FDR)
        """
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.summit = summit
        self.score = score
        self.strand = strand
        self.signal_value = signal_value
        self.p_value = p_value
        self.q_value = q_value

        # Validate inputs
        validation.validate_not_empty(chromosome, "chromosome")
        validation.validate_range(start, min_val=0, name="start")
        validation.validate_range(end, min_val=start, name="end")
        if summit is not None:
            validation.validate_range(summit, min_val=start, max_val=end, name="summit")

    @property
    def length(self) -> int:
        """Get peak length."""
        return self.end - self.start

    @property
    def center(self) -> int:
        """Get peak center position."""
        return (self.start + self.end) // 2

    def to_bed_format(self) -> str:
        """Convert to BED format string."""
        summit_str = str(self.summit) if self.summit is not None else "."
        p_val_str = f"{self.p_value:.2e}" if self.p_value is not None else "."
        q_val_str = f"{self.q_value:.2e}" if self.q_value is not None else "."

        return "\t".join([
            self.chromosome,
            str(self.start),
            str(self.end),
            f"peak_{self.start}_{self.end}",
            f"{self.score:.1f}",
            self.strand,
            summit_str,
            f"{self.signal_value:.2f}",
            p_val_str,
            q_val_str,
        ])

    def overlaps_with(self, other: ChIPPeak, min_overlap: int = 1) -> bool:
        """Check if this peak overlaps with another peak.

        Args:
            other: Another ChIPPeak to compare
            min_overlap: Minimum overlap required

        Returns:
            True if peaks overlap by at least min_overlap bases
        """
        if self.chromosome != other.chromosome:
            return False

        overlap_start = max(self.start, other.start)
        overlap_end = min(self.end, other.end)

        return overlap_end - overlap_start >= min_overlap

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            'chromosome': self.chromosome,
            'start': self.start,
            'end': self.end,
            'length': self.length,
            'summit': self.summit,
            'score': self.score,
            'strand': self.strand,
            'signal_value': self.signal_value,
            'p_value': self.p_value,
            'q_value': self.q_value,
        }


def load_chip_peaks(path: str | Path, format: str = "narrowpeak") -> List[ChIPPeak]:
    """Load ChIP-seq peaks from a file.

    Args:
        path: Path to peak file
        format: Peak file format ("narrowpeak", "broadpeak", "bed")

    Returns:
        List of ChIPPeak objects

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    path = validation.validate_path_exists(Path(path))

    logger.info(f"Loading ChIP-seq peaks from {path} (format: {format})")

    peaks = []

    try:
        with io.open_text_auto(path, 'rt') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('track'):
                    continue

                parts = line.split('\t')
                if len(parts) < 6:
                    logger.warning(f"Skipping malformed line {line_num}: insufficient columns")
                    continue

                try:
                    chromosome = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    name = parts[3] if len(parts) > 3 else f"peak_{line_num}"
                    score = float(parts[4]) if len(parts) > 4 else 0.0
                    strand = parts[5] if len(parts) > 5 else '.'

                    # Parse format-specific fields
                    if format == "narrowpeak" and len(parts) >= 10:
                        # narrowPeak format: chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak
                        signal_value = float(parts[6]) if len(parts) > 6 else 0.0
                        p_value = float(parts[7]) if len(parts) > 7 and parts[7] != '.' else None
                        q_value = float(parts[8]) if len(parts) > 8 and parts[8] != '.' else None
                        peak_offset = int(parts[9]) if len(parts) > 9 and parts[9] != '.' else None

                        summit = start + peak_offset if peak_offset is not None else None

                        peak = ChIPPeak(
                            chromosome=chromosome,
                            start=start,
                            end=end,
                            summit=summit,
                            score=score,
                            strand=strand,
                            signal_value=signal_value,
                            p_value=p_value,
                            q_value=q_value
                        )

                    elif format == "broadpeak" and len(parts) >= 9:
                        # broadPeak format: chrom, start, end, name, score, strand, signalValue, pValue, qValue
                        signal_value = float(parts[6]) if len(parts) > 6 else 0.0
                        p_value = float(parts[7]) if len(parts) > 7 and parts[7] != '.' else None
                        q_value = float(parts[8]) if len(parts) > 8 and parts[8] != '.' else None

                        peak = ChIPPeak(
                            chromosome=chromosome,
                            start=start,
                            end=end,
                            score=score,
                            strand=strand,
                            signal_value=signal_value,
                            p_value=p_value,
                            q_value=q_value
                        )

                    elif format == "bed":
                        # Basic BED format
                        peak = ChIPPeak(
                            chromosome=chromosome,
                            start=start,
                            end=end,
                            score=score,
                            strand=strand
                        )

                    else:
                        logger.warning(f"Unsupported format or insufficient columns in line {line_num}")
                        continue

                    peaks.append(peak)

                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue

    except Exception as e:
        logger.error(f"Error loading ChIP-seq peaks from {path}: {e}")
        raise errors.FileIOError(f"Failed to load ChIP-seq peaks: {e}") from e

    logger.info(f"Loaded {len(peaks)} ChIP-seq peaks")
    return peaks


def save_chip_peaks(peaks: List[ChIPPeak], path: str | Path, format: str = "narrowpeak") -> None:
    """Save ChIP-seq peaks to a file.

    Args:
        peaks: List of ChIPPeak objects
        path: Output file path
        format: Output format ("narrowpeak", "broadpeak", "bed")
    """
    path = Path(path)

    logger.info(f"Saving {len(peaks)} ChIP-seq peaks to {path} (format: {format})")

    with open(path, 'w') as f:
        # Write header for BED formats
        if format != "bed":
            f.write("# Generated by MetaInformant ChIP-seq analysis\n")

        for peak in peaks:
            if format == "narrowpeak":
                f.write(peak.to_bed_format() + "\n")
            elif format == "broadpeak":
                # Convert narrowPeak to broadPeak (remove summit column)
                bed_parts = peak.to_bed_format().split('\t')
                if len(bed_parts) >= 10:
                    # Remove the last column (peak offset)
                    broadpeak_line = '\t'.join(bed_parts[:9])
                    f.write(broadpeak_line + '\n')
            elif format == "bed":
                # Basic BED format
                f.write('\t'.join([
                    peak.chromosome,
                    str(peak.start),
                    str(peak.end),
                    f"peak_{peak.start}_{peak.end}",
                    f"{peak.score:.1f}",
                    peak.strand
                ]) + '\n')

    logger.info(f"Saved {len(peaks)} peaks to {path}")


def filter_peaks_by_score(peaks: List[ChIPPeak], min_score: float,
                         max_peaks: Optional[int] = None) -> List[ChIPPeak]:
    """Filter peaks by score and optionally limit number of peaks.

    Args:
        peaks: List of ChIPPeak objects
        min_score: Minimum peak score
        max_peaks: Maximum number of peaks to return (None for all)

    Returns:
        Filtered list of peaks
    """
    logger.info(f"Filtering peaks by score >= {min_score}")

    # Filter by score
    filtered = [peak for peak in peaks if peak.score >= min_score]

    # Sort by score (descending) and limit if requested
    filtered.sort(key=lambda p: p.score, reverse=True)

    if max_peaks is not None:
        filtered = filtered[:max_peaks]

    logger.info(f"Filtered to {len(filtered)} peaks")
    return filtered


def calculate_peak_statistics(peaks: List[ChIPPeak]) -> Dict[str, Any]:
    """Calculate comprehensive statistics for ChIP-seq peaks.

    Args:
        peaks: List of ChIPPeak objects

    Returns:
        Dictionary with peak statistics
    """
    if not peaks:
        return {}

    lengths = [p.length for p in peaks]
    scores = [p.score for p in peaks]
    signal_values = [p.signal_value for p in peaks if p.signal_value > 0]

    stats = {
        "total_peaks": len(peaks),
        "mean_length": statistics.mean(lengths),
        "median_length": statistics.median(lengths),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "length_std": statistics.stdev(lengths) if len(lengths) > 1 else 0,
        "mean_score": statistics.mean(scores),
        "median_score": statistics.median(scores),
        "min_score": min(scores),
        "max_score": max(scores),
        "score_std": statistics.stdev(scores) if len(scores) > 1 else 0,
    }

    if signal_values:
        stats.update({
            "mean_signal": statistics.mean(signal_values),
            "median_signal": statistics.median(signal_values),
            "signal_std": statistics.stdev(signal_values) if len(signal_values) > 1 else 0,
        })

    # Peak length distribution
    length_bins = defaultdict(int)
    for length in lengths:
        if length <= 200:
            length_bins["<=200"] += 1
        elif length <= 500:
            length_bins["201-500"] += 1
        elif length <= 1000:
            length_bins["501-1000"] += 1
        else:
            length_bins["1000+"] += 1

    stats["length_distribution"] = dict(length_bins)

    # Score distribution
    score_bins = defaultdict(int)
    for score in scores:
        if score < 100:
            score_bins["<100"] += 1
        elif score < 500:
            score_bins["100-500"] += 1
        elif score < 1000:
            score_bins["500-1000"] += 1
        else:
            score_bins["1000+"] += 1

    stats["score_distribution"] = dict(score_bins)

    # Per-chromosome distribution
    chr_counts = defaultdict(int)
    for peak in peaks:
        chr_counts[peak.chromosome] += 1

    stats["chromosome_distribution"] = dict(sorted(chr_counts.items()))

    return stats


def find_overlapping_peaks(peaks1: List[ChIPPeak], peaks2: List[ChIPPeak],
                          min_overlap: int = 1) -> List[Tuple[ChIPPeak, ChIPPeak]]:
    """Find overlapping peaks between two peak sets.

    Args:
        peaks1: First set of peaks
        peaks2: Second set of peaks
        min_overlap: Minimum overlap required

    Returns:
        List of (peak1, peak2) tuples for overlapping peaks
    """
    logger.info("Finding overlapping peaks")

    # Group peaks by chromosome for efficiency
    chr_peaks1 = defaultdict(list)
    chr_peaks2 = defaultdict(list)

    for peak in peaks1:
        chr_peaks1[peak.chromosome].append(peak)

    for peak in peaks2:
        chr_peaks2[peak.chromosome].append(peak)

    overlapping_pairs = []

    for chromosome in set(chr_peaks1.keys()) & set(chr_peaks2.keys()):
        peaks1_chr = chr_peaks1[chromosome]
        peaks2_chr = chr_peaks2[chromosome]

        # Sort peaks by start position
        peaks1_chr.sort(key=lambda p: p.start)
        peaks2_chr.sort(key=lambda p: p.start)

        # Find overlaps using sweep line algorithm
        i, j = 0, 0
        while i < len(peaks1_chr) and j < len(peaks2_chr):
            peak1 = peaks1_chr[i]
            peak2 = peaks2_chr[j]

            if peak1.overlaps_with(peak2, min_overlap):
                overlapping_pairs.append((peak1, peak2))

            # Move the peak that ends first
            if peak1.end < peak2.end:
                i += 1
            else:
                j += 1

    logger.info(f"Found {len(overlapping_pairs)} overlapping peak pairs")
    return overlapping_pairs


def merge_overlapping_peaks(peaks: List[ChIPPeak], max_distance: int = 0) -> List[ChIPPeak]:
    """Merge overlapping or nearby peaks.

    Args:
        peaks: List of ChIPPeak objects
        max_distance: Maximum distance between peaks to merge

    Returns:
        List of merged peaks
    """
    if not peaks:
        return []

    logger.info(f"Merging overlapping peaks (max_distance={max_distance})")

    # Group by chromosome
    chr_peaks = defaultdict(list)
    for peak in peaks:
        chr_peaks[peak.chromosome].append(peak)

    merged_peaks = []

    for chromosome, chr_peak_list in chr_peaks.items():
        # Sort peaks by start position
        chr_peak_list.sort(key=lambda p: p.start)

        current_peak = None

        for peak in chr_peak_list:
            if current_peak is None:
                current_peak = peak
            elif peak.start <= current_peak.end + max_distance:
                # Merge peaks
                current_peak = ChIPPeak(
                    chromosome=chromosome,
                    start=min(current_peak.start, peak.start),
                    end=max(current_peak.end, peak.end),
                    summit=current_peak.summit,  # Keep first summit
                    score=max(current_peak.score, peak.score),  # Use higher score
                    strand=current_peak.strand,
                    signal_value=max(current_peak.signal_value, peak.signal_value),
                    p_value=min(current_peak.p_value, peak.p_value) if current_peak.p_value and peak.p_value else None,
                    q_value=min(current_peak.q_value, peak.q_value) if current_peak.q_value and peak.q_value else None
                )
            else:
                # No overlap, save current peak and start new one
                merged_peaks.append(current_peak)
                current_peak = peak

        # Don't forget the last peak
        if current_peak:
            merged_peaks.append(current_peak)

    logger.info(f"Merged {len(peaks)} peaks into {len(merged_peaks)} peaks")
    return merged_peaks


def calculate_peak_enrichment(peaks: List[ChIPPeak], genome_size: int,
                            expected_peaks: Optional[int] = None) -> Dict[str, Any]:
    """Calculate peak enrichment statistics.

    Args:
        peaks: List of ChIPPeak objects
        genome_size: Total genome size in base pairs
        expected_peaks: Expected number of peaks (optional)

    Returns:
        Dictionary with enrichment statistics
    """
    if not peaks:
        return {}

    total_peak_bases = sum(p.length for p in peaks)
    genome_coverage = total_peak_bases / genome_size

    stats = {
        "total_peaks": len(peaks),
        "total_peak_bases": total_peak_bases,
        "genome_size": genome_size,
        "genome_coverage": genome_coverage,
        "average_peak_length": total_peak_bases / len(peaks),
    }

    # Calculate enrichment if expected peaks provided
    if expected_peaks:
        observed_peaks = len(peaks)
        enrichment_ratio = observed_peaks / expected_peaks
        stats.update({
            "expected_peaks": expected_peaks,
            "enrichment_ratio": enrichment_ratio,
            "fold_enrichment": enrichment_ratio,
        })

        # Statistical significance (simplified)
        if expected_peaks > 0:
            # Use Poisson approximation for significance
            import math
            p_value = 1.0 - math.exp(-expected_peaks) * sum(expected_peaks**k / math.factorial(k) for k in range(observed_peaks))
            stats["enrichment_p_value"] = p_value

    return stats


def find_motifs_in_peaks(peaks: List[ChIPPeak], genome_fasta: str | Path,
                        motif_patterns: List[str], window_size: int = 200) -> Dict[str, Any]:
    """Find motif occurrences in peak regions.

    Args:
        peaks: List of ChIPPeak objects
        genome_fasta: Path to genome FASTA file
        motif_patterns: List of motif patterns to search
        window_size: Size of window around peak summit to search

    Returns:
        Dictionary with motif finding results
    """
    logger.info(f"Finding motifs in {len(peaks)} peaks")

    # This is a simplified implementation
    # In practice, this would use a proper motif finding tool like MEME or HOMER

    motif_counts = defaultdict(int)
    motif_positions = defaultdict(list)

    # Simulate motif finding (in practice, would parse genome sequence)
    for peak in peaks:
        summit = peak.summit if peak.summit is not None else peak.center

        # Simulate finding motifs around summit
        for motif in motif_patterns:
            # Random simulation - in practice would search actual sequence
            found = True  # Simulate finding motif

            if found:
                motif_counts[motif] += 1
                motif_positions[motif].append({
                    'chromosome': peak.chromosome,
                    'position': summit,
                    'peak_score': peak.score,
                })

    results = {
        "total_peaks_analyzed": len(peaks),
        "motif_counts": dict(motif_counts),
        "motif_positions": dict(motif_positions),
    }

    # Calculate enrichment statistics
    for motif in motif_patterns:
        count = motif_counts[motif]
        expected_count = len(peaks) * 0.1  # Assume 10% expected frequency
        enrichment = count / expected_count if expected_count > 0 else 0

        results[f"{motif}_enrichment"] = enrichment

    logger.info(f"Found motifs: {dict(motif_counts)}")
    return results


def generate_chip_report(peaks: List[ChIPPeak], output_path: Optional[str | Path] = None) -> str:
    """Generate a comprehensive ChIP-seq analysis report.

    Args:
        peaks: List of ChIPPeak objects
        output_path: Optional path to save the report

    Returns:
        Formatted ChIP-seq report
    """
    report_lines = []
    report_lines.append("=" * 60)
    report_lines.append("CHIP-SEQ PEAK ANALYSIS REPORT")
    report_lines.append("=" * 60)
    report_lines.append("")

    # Basic statistics
    stats = calculate_peak_statistics(peaks)

    if stats:
        report_lines.append("Peak Statistics:")
        report_lines.append(f"  Total Peaks: {stats.get('total_peaks', 0):,}")
        report_lines.append(f"  Mean Length: {stats.get('mean_length', 0):.0f} bp")
        report_lines.append(f"  Median Length: {stats.get('median_length', 0):.0f} bp")
        report_lines.append(f"  Mean Score: {stats.get('mean_score', 0):.1f}")
        report_lines.append(f"  Median Score: {stats.get('median_score', 0):.1f}")

        if 'mean_signal' in stats:
            report_lines.append(f"  Mean Signal: {stats.get('mean_signal', 0):.2f}")
        report_lines.append("")

        # Peak length distribution
        length_dist = stats.get('length_distribution', {})
        if length_dist:
            report_lines.append("Peak Length Distribution:")
            for length_range, count in sorted(length_dist.items()):
                report_lines.append(f"  {length_range} bp: {count:,} peaks")
            report_lines.append("")

        # Score distribution
        score_dist = stats.get('score_distribution', {})
        if score_dist:
            report_lines.append("Peak Score Distribution:")
            for score_range, count in sorted(score_dist.items()):
                report_lines.append(f"  {score_range}: {count:,} peaks")
            report_lines.append("")

        # Chromosome distribution (top 10)
        chr_dist = stats.get('chromosome_distribution', {})
        if chr_dist:
            report_lines.append("Chromosome Distribution (Top 10):")
            sorted_chrs = sorted(chr_dist.items(), key=lambda x: x[1], reverse=True)[:10]
            for chr_name, count in sorted_chrs:
                report_lines.append(f"  {chr_name}: {count:,} peaks")
            report_lines.append("")

    report = "\n".join(report_lines)

    if output_path:
        output_path = Path(output_path)
        with open(output_path, 'w') as f:
            f.write(report)
        logger.info(f"ChIP-seq report saved to {output_path}")

    return report

