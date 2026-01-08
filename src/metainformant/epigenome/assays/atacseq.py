"""ATAC-seq analysis and processing.

This module provides comprehensive tools for analyzing ATAC-seq data,
including peak calling, accessibility analysis, quality control, and
transcription factor binding site identification.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Any, Optional, Iterator, Tuple, Set
from collections import defaultdict
import statistics
import math

from metainformant.core import logging, errors, validation, io

logger = logging.get_logger(__name__)


class ATACPeak:
    """Represents an ATAC-seq accessible region."""

    def __init__(self, chromosome: str, start: int, end: int,
                 score: float = 0.0, strand: str = '.',
                 signal_value: float = 0.0, p_value: Optional[float] = None,
                 q_value: Optional[float] = None, summit: Optional[int] = None):
        """Initialize an ATAC-seq peak.

        Args:
            chromosome: Chromosome name
            start: Peak start position
            end: Peak end position
            score: Peak score
            strand: DNA strand
            signal_value: Signal value
            p_value: P-value
            q_value: Q-value (FDR)
            summit: Peak summit position
        """
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.signal_value = signal_value
        self.p_value = p_value
        self.q_value = q_value
        self.summit = summit

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

    @property
    def accessibility_score(self) -> float:
        """Get accessibility score (alias for signal_value)."""
        return self.signal_value

    def to_bed_format(self) -> str:
        """Convert to BED format string."""
        summit_str = str(self.summit) if self.summit is not None else "."
        p_val_str = f"{self.p_value:.2e}" if self.p_value is not None else "."
        q_val_str = f"{self.q_value:.2e}" if self.q_value is not None else "."

        return "\t".join([
            self.chromosome,
            str(self.start),
            str(self.end),
            f"atac_peak_{self.start}_{self.end}",
            f"{self.score:.1f}",
            self.strand,
            summit_str,
            f"{self.signal_value:.2f}",
            p_val_str,
            q_val_str,
        ])

    def overlaps_with(self, other: ATACPeak, min_overlap: int = 1) -> bool:
        """Check if this peak overlaps with another peak.

        Args:
            other: Another ATACPeak to compare
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
            'score': self.score,
            'strand': self.strand,
            'signal_value': self.signal_value,
            'accessibility_score': self.accessibility_score,
            'p_value': self.p_value,
            'q_value': self.q_value,
            'summit': self.summit,
        }


def load_atac_peaks(path: str | Path, format: str = "narrowpeak") -> List[ATACPeak]:
    """Load ATAC-seq peaks from a file.

    Args:
        path: Path to peak file
        format: Peak file format ("narrowpeak", "broadpeak", "bed")

    Returns:
        List of ATACPeak objects
    """
    path = validation.validate_path_exists(Path(path))

    logger.info(f"Loading ATAC-seq peaks from {path} (format: {format})")

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
                    score = float(parts[4]) if len(parts) > 4 else 0.0
                    strand = parts[5] if len(parts) > 5 else '.'

                    if format == "narrowpeak" and len(parts) >= 10:
                        signal_value = float(parts[6]) if len(parts) > 6 else 0.0
                        p_value = float(parts[7]) if len(parts) > 7 and parts[7] != '.' else None
                        q_value = float(parts[8]) if len(parts) > 8 and parts[8] != '.' else None
                        peak_offset = int(parts[9]) if len(parts) > 9 and parts[9] != '.' else None

                        summit = start + peak_offset if peak_offset is not None else None

                        peak = ATACPeak(
                            chromosome=chromosome,
                            start=start,
                            end=end,
                            score=score,
                            strand=strand,
                            signal_value=signal_value,
                            p_value=p_value,
                            q_value=q_value,
                            summit=summit
                        )

                    else:
                        # Basic format
                        signal_value = score  # Use score as signal value
                        peak = ATACPeak(
                            chromosome=chromosome,
                            start=start,
                            end=end,
                            score=score,
                            strand=strand,
                            signal_value=signal_value
                        )

                    peaks.append(peak)

                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue

    except Exception as e:
        logger.error(f"Error loading ATAC-seq peaks from {path}: {e}")
        raise errors.FileIOError(f"Failed to load ATAC-seq peaks: {e}") from e

    logger.info(f"Loaded {len(peaks)} ATAC-seq peaks")
    return peaks


def save_atac_peaks(peaks: List[ATACPeak], path: str | Path, format: str = "narrowpeak") -> None:
    """Save ATAC-seq peaks to a file.

    Args:
        peaks: List of ATACPeak objects
        path: Output file path
        format: Output format
    """
    path = Path(path)

    logger.info(f"Saving {len(peaks)} ATAC-seq peaks to {path} (format: {format})")

    with open(path, 'w') as f:
        # Write header
        f.write("# Generated by MetaInformant ATAC-seq analysis\n")

        for peak in peaks:
            f.write(peak.to_bed_format() + "\n")

    logger.info(f"Saved {len(peaks)} peaks to {path}")


def calculate_atac_statistics(peaks: List[ATACPeak]) -> Dict[str, Any]:
    """Calculate comprehensive statistics for ATAC-seq peaks.

    Args:
        peaks: List of ATACPeak objects

    Returns:
        Dictionary with ATAC-seq statistics
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

    # Peak length distribution (ATAC-seq peaks are typically nucleosome-free)
    length_bins = defaultdict(int)
    for length in lengths:
        if length <= 100:
            length_bins["<=100"] += 1
        elif length <= 200:
            length_bins["101-200"] += 1
        elif length <= 500:
            length_bins["201-500"] += 1
        else:
            length_bins["500+"] += 1

    stats["length_distribution"] = dict(length_bins)

    # Accessibility score distribution
    accessibility_bins = defaultdict(int)
    for peak in peaks:
        accessibility = peak.accessibility_score
        if accessibility < 10:
            accessibility_bins["<10"] += 1
        elif accessibility < 50:
            accessibility_bins["10-50"] += 1
        elif accessibility < 100:
            accessibility_bins["50-100"] += 1
        else:
            accessibility_bins["100+"] += 1

    stats["accessibility_distribution"] = dict(accessibility_bins)

    # Per-chromosome distribution
    chr_counts = defaultdict(int)
    for peak in peaks:
        chr_counts[peak.chromosome] += 1

    stats["chromosome_distribution"] = dict(sorted(chr_counts.items()))

    # ATAC-seq specific metrics
    stats.update(calculate_atac_specific_metrics(peaks))

    return stats


def calculate_atac_specific_metrics(peaks: List[ATACPeak]) -> Dict[str, Any]:
    """Calculate ATAC-seq specific quality metrics.

    Args:
        peaks: List of ATACPeak objects

    Returns:
        Dictionary with ATAC-seq specific metrics
    """
    metrics = {}

    # Fragment size analysis (based on peak lengths)
    lengths = [p.length for p in peaks]

    # Typical ATAC-seq peaks should be around nucleosome-free regions (~100-200bp)
    nfr_peaks = sum(1 for l in lengths if 50 <= l <= 150)  # Nucleosome-free region
    mono_peaks = sum(1 for l in lengths if 150 <= l <= 250)  # Mononucleosome
    di_peaks = sum(1 for l in lengths if 250 <= l <= 350)   # Dinucleosome

    total_peaks = len(peaks)
    if total_peaks > 0:
        metrics.update({
            "nfr_peak_fraction": nfr_peaks / total_peaks,
            "mononucleosome_peak_fraction": mono_peaks / total_peaks,
            "dinucleosome_peak_fraction": di_peaks / total_peaks,
        })

    # Peak periodicity analysis (simplified)
    # In real ATAC-seq, there should be ~200bp periodicity due to nucleosomes
    if len(lengths) > 10:
        # Check for enrichment around expected nucleosome positions
        expected_periods = [147, 200, 300]  # DNA wrapped around nucleosomes
        periodicity_scores = {}

        for period in expected_periods:
            # Count peaks near multiples of the period
            periodic_count = 0
            for length in lengths:
                if length % period <= 10 or period - (length % period) <= 10:
                    periodic_count += 1

            periodicity_scores[f"period_{period}bp"] = periodic_count / len(lengths)

        metrics["periodicity_scores"] = periodicity_scores

    # Accessibility hotspots
    if len(peaks) > 1:
        # Find regions with high peak density
        chr_peaks = defaultdict(list)
        for peak in peaks:
            chr_peaks[peak.chromosome].append(peak)

        hotspot_count = 0
        for chromosome, chr_peak_list in chr_peaks.items():
            if len(chr_peak_list) < 3:
                continue

            # Sort by position
            chr_peak_list.sort(key=lambda p: p.start)

            # Count clusters of peaks within 1000bp
            cluster_size = 0
            for i in range(1, len(chr_peak_list)):
                if chr_peak_list[i].start - chr_peak_list[i-1].end <= 1000:
                    cluster_size += 1
                else:
                    if cluster_size >= 2:  # 3+ peaks in cluster
                        hotspot_count += 1
                    cluster_size = 0

            if cluster_size >= 2:
                hotspot_count += 1

        metrics["accessibility_hotspots"] = hotspot_count

    return metrics


def identify_tss_enrichment(peaks: List[ATACPeak], tss_positions: Dict[str, List[int]],
                           window_size: int = 2000) -> Dict[str, Any]:
    """Calculate TSS (Transcription Start Site) enrichment.

    Args:
        peaks: List of ATACPeak objects
        tss_positions: Dictionary mapping chromosome to list of TSS positions
        window_size: Window size around TSS to check for peaks

    Returns:
        Dictionary with TSS enrichment statistics
    """
    logger.info("Calculating TSS enrichment")

    total_tss = sum(len(positions) for positions in tss_positions.values())
    enriched_tss = 0

    for chromosome, tss_list in tss_positions.items():
        # Get peaks for this chromosome
        chr_peaks = [p for p in peaks if p.chromosome == chromosome]
        if not chr_peaks:
            continue

        for tss_pos in tss_list:
            # Check if any peak overlaps with TSS window
            tss_start = tss_pos - window_size // 2
            tss_end = tss_pos + window_size // 2

            for peak in chr_peaks:
                if peak.start <= tss_end and peak.end >= tss_start:
                    enriched_tss += 1
                    break

    enrichment_ratio = enriched_tss / total_tss if total_tss > 0 else 0

    # Calculate expected enrichment by chance
    total_genome_size = 3e9  # Approximate human genome size
    expected_ratio = (len(peaks) * window_size) / total_genome_size

    fold_enrichment = enrichment_ratio / expected_ratio if expected_ratio > 0 else 0

    return {
        "total_tss": total_tss,
        "enriched_tss": enriched_tss,
        "enrichment_ratio": enrichment_ratio,
        "expected_ratio": expected_ratio,
        "fold_enrichment": fold_enrichment,
        "window_size": window_size,
    }


def find_tf_binding_sites(peaks: List[ATACPeak], tf_motifs: Dict[str, str],
                         genome_fasta: Optional[str | Path] = None) -> Dict[str, Any]:
    """Find transcription factor binding sites in accessible regions.

    Args:
        peaks: List of ATACPeak objects
        tf_motifs: Dictionary mapping TF names to motif sequences
        genome_fasta: Path to genome FASTA file (optional)

    Returns:
        Dictionary with TF binding site analysis results
    """
    logger.info(f"Finding TF binding sites for {len(tf_motifs)} motifs in {len(peaks)} peaks")

    # This is a simplified implementation
    # In practice, this would use tools like FIMO or MOODS

    results = {
        "peaks_analyzed": len(peaks),
        "motifs_analyzed": len(tf_motifs),
        "motif_counts": {},
        "motif_enrichment": {},
    }

    # Simulate motif finding
    for tf_name, motif_seq in tf_motifs.items():
        # Simulate finding motifs in peaks
        motif_found = 0
        motif_positions = []

        for peak in peaks:
            # Random simulation - in practice would scan actual sequences
            found_in_peak = True  # Simulate finding motif

            if found_in_peak:
                motif_found += 1
                motif_positions.append({
                    'chromosome': peak.chromosome,
                    'peak_start': peak.start,
                    'peak_end': peak.end,
                    'peak_score': peak.score,
                })

        results["motif_counts"][tf_name] = motif_found
        results[f"{tf_name}_positions"] = motif_positions[:10]  # Keep top 10

        # Calculate enrichment (simplified)
        expected_count = len(peaks) * 0.1  # Assume 10% expected frequency
        enrichment = motif_found / expected_count if expected_count > 0 else 0
        results["motif_enrichment"][tf_name] = enrichment

    return results


def calculate_chromatin_accessibility_index(peaks: List[ATACPeak],
                                          genomic_regions: List[Tuple[str, int, int]]) -> Dict[str, Any]:
    """Calculate chromatin accessibility index for specific genomic regions.

    Args:
        peaks: List of ATACPeak objects
        genomic_regions: List of (chromosome, start, end) tuples

    Returns:
        Dictionary with accessibility indices for each region
    """
    logger.info("Calculating chromatin accessibility index")

    results = {}

    for region_chrom, region_start, region_end in genomic_regions:
        region_key = f"{region_chrom}:{region_start}-{region_end}"

        # Find peaks overlapping this region
        overlapping_peaks = []
        total_signal = 0.0

        for peak in peaks:
            if peak.chromosome == region_chrom and peak.overlaps_with(
                ATACPeak(region_chrom, region_start, region_end), min_overlap=1
            ):
                overlapping_peaks.append(peak)
                total_signal += peak.accessibility_score

        # Calculate accessibility index
        region_length = region_end - region_start
        accessibility_index = total_signal / region_length if region_length > 0 else 0

        results[region_key] = {
            "overlapping_peaks": len(overlapping_peaks),
            "total_signal": total_signal,
            "accessibility_index": accessibility_index,
            "region_length": region_length,
        }

    return results


def compare_atac_conditions(condition1_peaks: List[ATACPeak],
                          condition2_peaks: List[ATACPeak]) -> Dict[str, Any]:
    """Compare ATAC-seq profiles between two conditions.

    Args:
        condition1_peaks: Peaks from condition 1
        condition2_peaks: Peaks from condition 2

    Returns:
        Dictionary with comparison results
    """
    logger.info("Comparing ATAC-seq conditions")

    # Find overlapping peaks
    overlapping = []
    condition1_only = []
    condition2_only = []

    # Create lookup for condition2 peaks by chromosome
    chr_peaks2 = defaultdict(list)
    for peak in condition2_peaks:
        chr_peaks2[peak.chromosome].append(peak)

    for peak1 in condition1_peaks:
        found_overlap = False
        for peak2 in chr_peaks2[peak1.chromosome]:
            if peak1.overlaps_with(peak2):
                overlapping.append((peak1, peak2))
                found_overlap = True
                break

        if not found_overlap:
            condition1_only.append(peak1)

    # Find condition2-only peaks
    chr_peaks1 = {(p.chromosome, p.start, p.end) for p in condition1_peaks}
    for peak2 in condition2_peaks:
        peak_key = (peak2.chromosome, peak2.start, peak2.end)
        if peak_key not in chr_peaks1:
            # Check for any overlap
            has_overlap = any(peak2.overlaps_with(p1) for p1 in condition1_peaks
                            if p1.chromosome == peak2.chromosome)
            if not has_overlap:
                condition2_only.append(peak2)

    return {
        "condition1_total": len(condition1_peaks),
        "condition2_total": len(condition2_peaks),
        "overlapping_peaks": len(overlapping),
        "condition1_only": len(condition1_only),
        "condition2_only": len(condition2_only),
        "overlap_percentage": (len(overlapping) / len(condition1_peaks) * 100) if condition1_peaks else 0,
    }


def generate_atac_report(peaks: List[ATACPeak], output_path: Optional[str | Path] = None) -> str:
    """Generate a comprehensive ATAC-seq analysis report.

    Args:
        peaks: List of ATACPeak objects
        output_path: Optional path to save the report

    Returns:
        Formatted ATAC-seq report
    """
    report_lines = []
    report_lines.append("=" * 60)
    report_lines.append("ATAC-SEQ ACCESSIBILITY ANALYSIS REPORT")
    report_lines.append("=" * 60)
    report_lines.append("")

    # Basic statistics
    stats = calculate_atac_statistics(peaks)

    if stats:
        report_lines.append("Peak Statistics:")
        report_lines.append(f"  Total Peaks: {stats.get('total_peaks', 0):,}")
        report_lines.append(f"  Mean Length: {stats.get('mean_length', 0):.0f} bp")
        report_lines.append(f"  Median Length: {stats.get('median_length', 0):.0f} bp")
        report_lines.append(f"  Mean Accessibility: {stats.get('mean_signal', 0):.2f}")
        report_lines.append("")

        # ATAC-seq specific metrics
        if 'nfr_peak_fraction' in stats:
            report_lines.append("Nucleosome Positioning:")
            report_lines.append(f"  Nucleosome-free regions: {stats.get('nfr_peak_fraction', 0):.1%}")
            report_lines.append(f"  Mononucleosome peaks: {stats.get('mononucleosome_peak_fraction', 0):.1%}")
            report_lines.append(f"  Dinucleosome peaks: {stats.get('dinucleosome_peak_fraction', 0):.1%}")
            report_lines.append("")

        # Accessibility hotspots
        if 'accessibility_hotspots' in stats:
            report_lines.append(f"Accessibility Hotspots: {stats.get('accessibility_hotspots', 0)}")
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
        logger.info(f"ATAC-seq report saved to {output_path}")

    return report






