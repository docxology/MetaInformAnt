"""DNA methylation analysis and processing.

This module provides comprehensive tools for analyzing DNA methylation data,
including loading methylation calls, calculating methylation levels,
differential methylation analysis, and methylation pattern recognition.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Any, Optional, Iterator, Tuple, Set
from collections import defaultdict
import statistics
import math

from metainformant.core import logging, errors, validation, io

logger = logging.get_logger(__name__)


class MethylationSite:
    """Represents a single CpG methylation site."""

    def __init__(self, chromosome: str, position: int, methylated_reads: int,
                 total_reads: int, strand: str = '+'):
        """Initialize a methylation site.

        Args:
            chromosome: Chromosome name
            position: Genomic position
            methylated_reads: Number of methylated reads
            total_reads: Total number of reads covering this site
            strand: DNA strand ('+' or '-')
        """
        self.chromosome = chromosome
        self.position = position
        self.methylated_reads = methylated_reads
        self.total_reads = total_reads
        self.strand = strand

        # Validate inputs
        validation.validate_not_empty(chromosome, "chromosome")
        validation.validate_range(position, min_val=0, name="position")
        validation.validate_range(methylated_reads, min_val=0, name="methylated_reads")
        validation.validate_range(total_reads, min_val=methylated_reads, name="total_reads")

    @property
    def methylation_level(self) -> float:
        """Calculate methylation level (beta value)."""
        if self.total_reads == 0:
            return 0.0
        return self.methylated_reads / self.total_reads

    @property
    def coverage(self) -> int:
        """Get read coverage."""
        return self.total_reads

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            'chromosome': self.chromosome,
            'position': self.position,
            'methylated_reads': self.methylated_reads,
            'total_reads': self.total_reads,
            'methylation_level': self.methylation_level,
            'coverage': self.coverage,
            'strand': self.strand,
        }


def load_methylation_bedgraph(path: str | Path, min_coverage: int = 1) -> Dict[str, List[MethylationSite]]:
    """Load methylation data from BEDgraph format.

    Args:
        path: Path to BEDgraph file
        min_coverage: Minimum read coverage to include

    Returns:
        Dictionary mapping chromosome names to lists of MethylationSite objects

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    path = validation.validate_path_exists(Path(path))

    logger.info(f"Loading methylation BEDgraph from {path}")

    methylation_data = defaultdict(list)

    try:
        with io.open_text_auto(path, 'rt') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('track'):
                    continue

                parts = line.split('\t')
                if len(parts) < 4:
                    logger.warning(f"Skipping malformed line {line_num}: {line}")
                    continue

                try:
                    chromosome = parts[0]
                    start_pos = int(parts[1])
                    end_pos = int(parts[2])
                    methylation_level = float(parts[3])

                    # BEDgraph format typically uses the start position
                    # Convert methylation level (0-1) to read counts (assuming coverage of 1 for simplicity)
                    # In practice, you might need additional coverage information
                    total_reads = max(1, int(1 / max(0.01, methylation_level)))  # Estimate coverage
                    methylated_reads = int(total_reads * methylation_level)

                    if total_reads >= min_coverage:
                        site = MethylationSite(
                            chromosome=chromosome,
                            position=start_pos,
                            methylated_reads=methylated_reads,
                            total_reads=total_reads
                        )
                        methylation_data[chromosome].append(site)

                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue

    except Exception as e:
        logger.error(f"Error loading BEDgraph file {path}: {e}")
        raise errors.FileIOError(f"Failed to load methylation BEDgraph: {e}") from e

    # Sort sites by position for each chromosome
    for chromosome in methylation_data:
        methylation_data[chromosome].sort(key=lambda x: x.position)

    total_sites = sum(len(sites) for sites in methylation_data.values())
    logger.info(f"Loaded {total_sites} methylation sites across {len(methylation_data)} chromosomes")

    return dict(methylation_data)


def load_methylation_cov(path: str | Path, min_coverage: int = 1) -> Dict[str, List[MethylationSite]]:
    """Load methylation data from Bismark coverage format.

    Args:
        path: Path to .cov file
        min_coverage: Minimum read coverage to include

    Returns:
        Dictionary mapping chromosome names to lists of MethylationSite objects
    """
    path = validation.validate_path_exists(Path(path))

    logger.info(f"Loading methylation coverage from {path}")

    methylation_data = defaultdict(list)

    try:
        with io.open_text_auto(path, 'rt') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 6:
                    logger.warning(f"Skipping malformed line {line_num}: {line}")
                    continue

                try:
                    chromosome = parts[0]
                    position = int(parts[1])
                    methylation_level = float(parts[3]) / 100.0  # Convert percentage to fraction
                    methylated_reads = int(parts[4])
                    total_reads = int(parts[5])

                    if total_reads >= min_coverage:
                        site = MethylationSite(
                            chromosome=chromosome,
                            position=position,
                            methylated_reads=methylated_reads,
                            total_reads=total_reads,
                            strand=parts[2] if len(parts) > 2 else '+'
                        )
                        methylation_data[chromosome].append(site)

                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue

    except Exception as e:
        logger.error(f"Error loading coverage file {path}: {e}")
        raise errors.FileIOError(f"Failed to load methylation coverage: {e}") from e

    # Sort sites by position
    for chromosome in methylation_data:
        methylation_data[chromosome].sort(key=lambda x: x.position)

    total_sites = sum(len(sites) for sites in methylation_data.values())
    logger.info(f"Loaded {total_sites} methylation sites across {len(methylation_data)} chromosomes")

    return dict(methylation_data)


def calculate_methylation_statistics(methylation_data: Dict[str, List[MethylationSite]]) -> Dict[str, Any]:
    """Calculate comprehensive methylation statistics.

    Args:
        methylation_data: Dictionary of methylation sites by chromosome

    Returns:
        Dictionary with methylation statistics
    """
    if not methylation_data:
        return {}

    all_sites = []
    for sites in methylation_data.values():
        all_sites.extend(sites)

    if not all_sites:
        return {}

    # Basic statistics
    methylation_levels = [site.methylation_level for site in all_sites]
    coverages = [site.coverage for site in all_sites]

    stats = {
        "total_sites": len(all_sites),
        "total_chromosomes": len(methylation_data),
        "mean_methylation": statistics.mean(methylation_levels),
        "median_methylation": statistics.median(methylation_levels),
        "methylation_std": statistics.stdev(methylation_levels) if len(methylation_levels) > 1 else 0,
        "mean_coverage": statistics.mean(coverages),
        "median_coverage": statistics.median(coverages),
        "min_coverage": min(coverages),
        "max_coverage": max(coverages),
    }

    # Coverage distribution
    coverage_bins = defaultdict(int)
    for coverage in coverages:
        if coverage <= 10:
            coverage_bins[str(coverage)] += 1
        elif coverage <= 50:
            coverage_bins["11-50"] += 1
        else:
            coverage_bins["50+"] += 1

    stats["coverage_distribution"] = dict(coverage_bins)

    # Methylation distribution
    methylation_bins = defaultdict(int)
    for level in methylation_levels:
        bin_key = f"{int(level * 10) / 10:.1f}-{int(level * 10) / 10 + 0.1:.1f}"
        methylation_bins[bin_key] += 1

    stats["methylation_distribution"] = dict(methylation_bins)

    # Per-chromosome statistics
    chromosome_stats = {}
    for chromosome, sites in methylation_data.items():
        if sites:
            levels = [s.methylation_level for s in sites]
            chromosome_stats[chromosome] = {
                "sites": len(sites),
                "mean_methylation": statistics.mean(levels),
                "median_methylation": statistics.median(levels),
            }

    stats["chromosome_stats"] = chromosome_stats

    return stats


def find_differentially_methylated_regions(site_data1: Dict[str, List[MethylationSite]],
                                          site_data2: Dict[str, List[MethylationSite]],
                                          delta_threshold: float = 0.2,
                                          p_value_threshold: float = 0.05,
                                          min_sites: int = 3) -> List[Dict[str, Any]]:
    """Find differentially methylated regions between two conditions.

    Args:
        site_data1: Methylation data for condition 1
        site_data2: Methylation data for condition 2
        delta_threshold: Minimum methylation difference
        p_value_threshold: Maximum p-value for significance
        min_sites: Minimum number of consecutive sites for DMR

    Returns:
        List of differentially methylated regions
    """
    logger.info("Finding differentially methylated regions")

    dmr_list = []

    # Get all chromosomes present in both datasets
    chromosomes = set(site_data1.keys()) & set(site_data2.keys())

    for chromosome in chromosomes:
        sites1 = site_data1[chromosome]
        sites2 = site_data2[chromosome]

        if not sites1 or not sites2:
            continue

        # Create position-to-site mapping
        pos_to_site1 = {site.position: site for site in sites1}
        pos_to_site2 = {site.position: site for site in sites2}

        # Find overlapping positions
        positions = sorted(set(pos_to_site1.keys()) & set(pos_to_site2.keys()))

        if len(positions) < min_sites:
            continue

        # Calculate differences and identify potential DMRs
        current_dmr = None

        for i, pos in enumerate(positions):
            site1 = pos_to_site1[pos]
            site2 = pos_to_site2[pos]

            delta = abs(site1.methylation_level - site2.methylation_level)

            # Simple proportion test for significance (could be improved with proper statistical test)
            if site1.coverage > 0 and site2.coverage > 0:
                # Use Fisher's exact test approximation for difference significance
                try:
                    # Calculate expected difference significance
                    p_diff = abs(site1.methylated_reads/site1.coverage - site2.methylated_reads/site2.coverage)
                    # Simplified p-value calculation (in practice, use proper statistical test)
                    se_diff = math.sqrt((site1.methylation_level * (1 - site1.methylation_level) / site1.coverage) +
                                      (site2.methylation_level * (1 - site2.methylation_level) / site2.coverage))
                    z_score = p_diff / se_diff if se_diff > 0 else 0
                    p_value = 2 * (1 - statistics.NormalDist().cdf(abs(z_score)))  # Two-tailed
                except:
                    p_value = 1.0
            else:
                p_value = 1.0

            if delta >= delta_threshold and p_value <= p_value_threshold:
                if current_dmr is None:
                    # Start new DMR
                    current_dmr = {
                        'chromosome': chromosome,
                        'start': pos,
                        'end': pos,
                        'sites': [pos],
                        'mean_delta': delta,
                        'mean_p_value': p_value,
                        'direction': 'hyper' if site1.methylation_level > site2.methylation_level else 'hypo'
                    }
                else:
                    # Extend current DMR
                    current_dmr['end'] = pos
                    current_dmr['sites'].append(pos)
                    current_dmr['mean_delta'] = (current_dmr['mean_delta'] * (len(current_dmr['sites']) - 1) + delta) / len(current_dmr['sites'])
                    current_dmr['mean_p_value'] = min(current_dmr['mean_p_value'], p_value)
            else:
                # End current DMR if it meets minimum size
                if current_dmr and len(current_dmr['sites']) >= min_sites:
                    dmr_list.append(current_dmr)
                current_dmr = None

        # Don't forget the last DMR
        if current_dmr and len(current_dmr['sites']) >= min_sites:
            dmr_list.append(current_dmr)

    logger.info(f"Found {len(dmr_list)} differentially methylated regions")
    return dmr_list


def identify_cpg_islands(methylation_data: Dict[str, List[MethylationSite]],
                        window_size: int = 200,
                        step_size: int = 50,
                        min_gc: float = 0.5,
                        min_cpg_ratio: float = 0.6) -> List[Dict[str, Any]]:
    """Identify CpG islands based on methylation patterns.

    Args:
        methylation_data: Dictionary of methylation sites by chromosome
        window_size: Size of sliding window in base pairs
        step_size: Step size for sliding window
        min_gc: Minimum GC content for CpG island
        min_cpg_ratio: Minimum CpG ratio (observed/expected)

    Returns:
        List of identified CpG islands
    """
    logger.info("Identifying CpG islands")

    cpg_islands = []

    for chromosome, sites in methylation_data.items():
        if len(sites) < 10:  # Need minimum sites for analysis
            continue

        # Sort sites by position
        sites.sort(key=lambda x: x.position)

        # Sliding window analysis
        max_pos = max(site.position for site in sites)
        min_pos = min(site.position for site in sites)

        for start_pos in range(min_pos, max_pos - window_size + 1, step_size):
            end_pos = start_pos + window_size

            # Get sites in window
            window_sites = [s for s in sites if start_pos <= s.position < end_pos]

            if len(window_sites) < 5:  # Need minimum sites in window
                continue

            # Calculate GC content and CpG ratio
            gc_count = 0
            cpg_count = 0
            total_bases = 0

            # Simplified calculation based on methylation levels
            # In practice, this would use actual sequence data
            avg_methylation = sum(s.methylation_level for s in window_sites) / len(window_sites)

            # Estimate GC content from methylation patterns
            # (highly methylated regions often have higher GC content)
            gc_content = min(1.0, 0.3 + avg_methylation * 0.4)  # Simplified estimation

            # Estimate CpG ratio
            cpg_ratio = avg_methylation * 1.2  # Simplified estimation

            if gc_content >= min_gc and cpg_ratio >= min_cpg_ratio:
                cpg_islands.append({
                    'chromosome': chromosome,
                    'start': start_pos,
                    'end': end_pos,
                    'sites': len(window_sites),
                    'avg_methylation': avg_methylation,
                    'gc_content': gc_content,
                    'cpg_ratio': cpg_ratio,
                })

    # Remove overlapping islands (keep the one with highest CpG ratio)
    cpg_islands.sort(key=lambda x: x['cpg_ratio'], reverse=True)

    filtered_islands = []
    for island in cpg_islands:
        # Check for overlap with already selected islands
        overlaps = False
        for selected in filtered_islands:
            if (island['chromosome'] == selected['chromosome'] and
                not (island['end'] <= selected['start'] or island['start'] >= selected['end'])):
                overlaps = True
                break

        if not overlaps:
            filtered_islands.append(island)

    logger.info(f"Identified {len(filtered_islands)} CpG islands")
    return filtered_islands


def calculate_methylation_entropy(methylation_data: Dict[str, List[MethylationSite]],
                                window_size: int = 1000) -> Dict[str, List[Tuple[int, float]]]:
    """Calculate methylation entropy across genomic regions.

    Args:
        methylation_data: Dictionary of methylation sites by chromosome
        window_size: Size of genomic windows for entropy calculation

    Returns:
        Dictionary mapping chromosomes to lists of (position, entropy) tuples
    """
    logger.info("Calculating methylation entropy")

    entropy_profiles = {}

    for chromosome, sites in methylation_data.items():
        if len(sites) < 10:
            continue

        sites.sort(key=lambda x: x.position)

        # Sliding window entropy calculation
        entropy_values = []
        max_pos = max(site.position for site in sites)

        for start_pos in range(sites[0].position, max_pos - window_size + 1, window_size // 2):
            end_pos = start_pos + window_size

            # Get sites in window
            window_sites = [s for s in sites if start_pos <= s.position < end_pos]

            if len(window_sites) < 5:
                continue

            # Calculate methylation distribution
            methylation_levels = [s.methylation_level for s in window_sites]

            # Create histogram bins
            bins = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
            hist = [0] * (len(bins) - 1)

            for level in methylation_levels:
                for i in range(len(bins) - 1):
                    if bins[i] <= level < bins[i + 1]:
                        hist[i] += 1
                        break

            # Calculate Shannon entropy
            total_sites = len(window_sites)
            entropy = 0.0

            for count in hist:
                if count > 0:
                    p = count / total_sites
                    entropy -= p * math.log2(p)

            entropy_values.append((start_pos + window_size // 2, entropy))

        entropy_profiles[chromosome] = entropy_values

    logger.info(f"Calculated entropy profiles for {len(entropy_profiles)} chromosomes")
    return entropy_profiles


def export_methylation_bedgraph(methylation_data: Dict[str, List[MethylationSite]],
                              output_path: str | Path) -> None:
    """Export methylation data to BEDgraph format.

    Args:
        output_path: Path for output BEDgraph file
        methylation_data: Dictionary of methylation sites by chromosome
    """
    output_path = Path(output_path)

    logger.info(f"Exporting methylation data to BEDgraph: {output_path}")

    with open(output_path, 'w') as f:
        # Write track header
        f.write('track type=bedGraph name="Methylation" description="DNA Methylation Levels" visibility=full color=200,100,0\n')

        for chromosome in sorted(methylation_data.keys()):
            sites = methylation_data[chromosome]
            sites.sort(key=lambda x: x.position)

            for site in sites:
                # BEDgraph format: chrom start end value
                f.write(f"{chromosome}\t{site.position}\t{site.position + 1}\t{site.methylation_level:.3f}\n")

    logger.info(f"Exported {sum(len(sites) for sites in methylation_data.values())} sites to {output_path}")


def load_cpg_table(path: str | Path) -> pd.DataFrame:
    """Load CpG methylation data from TSV file.

    Args:
        path: Path to TSV file with columns: chrom, pos, methylated, unmethylated

    Returns:
        DataFrame with methylation data
    """
    import pandas as pd

    df = pd.read_csv(path, sep='\t')
    required_cols = ['chrom', 'pos', 'methylated', 'unmethylated']

    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV file must contain columns: {required_cols}")

    return df


def compute_beta_values(df: pd.DataFrame) -> pd.DataFrame:
    """Compute beta values from methylated/unmethylated counts.

    Args:
        df: DataFrame with 'methylated' and 'unmethylated' columns

    Returns:
        DataFrame with added 'beta' column
    """
    if 'methylated' not in df.columns or 'unmethylated' not in df.columns:
        raise ValueError("DataFrame must contain 'methylated' and 'unmethylated' columns")

    df_copy = df.copy()
    total_reads = df_copy['methylated'] + df_copy['unmethylated']
    df_copy['beta'] = df_copy['methylated'] / total_reads
    df_copy['beta'] = df_copy['beta'].fillna(0.0)  # Handle division by zero

    return df_copy


def summarize_beta_by_chromosome(df: pd.DataFrame) -> pd.DataFrame:
    """Summarize beta values by chromosome.

    Args:
        df: DataFrame with 'chrom' and 'beta' columns

    Returns:
        DataFrame with summary statistics by chromosome
    """
    if 'chrom' not in df.columns or 'beta' not in df.columns:
        raise ValueError("DataFrame must contain 'chrom' and 'beta' columns")

    return df.groupby('chrom')['beta'].agg(['mean', 'std', 'min', 'max', 'count'])


def generate_methylation_report(methylation_data: Dict[str, List[MethylationSite]],
                               output_path: Optional[str | Path] = None) -> str:
    """Generate a comprehensive methylation analysis report.

    Args:
        methylation_data: Dictionary of methylation sites by chromosome
        output_path: Optional path to save the report

    Returns:
        Formatted methylation report
    """
    report_lines = []
    report_lines.append("=" * 60)
    report_lines.append("DNA METHYLATION ANALYSIS REPORT")
    report_lines.append("=" * 60)
    report_lines.append("")

    # Basic statistics
    stats = calculate_methylation_statistics(methylation_data)

    if stats:
        report_lines.append("Global Statistics:")
        report_lines.append(f"  Total Sites: {stats.get('total_sites', 0):,}")
        report_lines.append(f"  Chromosomes: {stats.get('total_chromosomes', 0)}")
        report_lines.append(f"  Mean Methylation: {stats.get('mean_methylation', 0):.3f}")
        report_lines.append(f"  Median Methylation: {stats.get('median_methylation', 0):.3f}")
        report_lines.append(f"  Mean Coverage: {stats.get('mean_coverage', 0):.1f}")
        report_lines.append("")

        # Coverage distribution
        coverage_dist = stats.get('coverage_distribution', {})
        if coverage_dist:
            report_lines.append("Coverage Distribution:")
            for cov_range, count in sorted(coverage_dist.items()):
                report_lines.append(f"  {cov_range}x: {count:,} sites")
            report_lines.append("")

        # Top chromosomes
        chr_stats = stats.get('chromosome_stats', {})
        if chr_stats:
            report_lines.append("Per-Chromosome Statistics (Top 5):")
            sorted_chrs = sorted(chr_stats.items(), key=lambda x: x[1]['sites'], reverse=True)[:5]
            for chr_name, chr_data in sorted_chrs:
                report_lines.append(f"  {chr_name}: {chr_data['sites']:,} sites, {chr_data['mean_methylation']:.3f} mean methylation")
            report_lines.append("")

    report = "\n".join(report_lines)

    if output_path:
        output_path = Path(output_path)
        with open(output_path, 'w') as f:
            f.write(report)
        logger.info(f"Methylation report saved to {output_path}")

    return report

