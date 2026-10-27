"""DNA methylation analysis and processing.

This module provides comprehensive tools for analyzing DNA methylation data,
including loading methylation calls, calculating methylation levels,
differential methylation analysis, and methylation pattern recognition.
"""

from __future__ import annotations

import math
import statistics
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Optional, Tuple

from metainformant.core import io
from metainformant.core.data import validation
from metainformant.core.utils import errors, logging

if TYPE_CHECKING:
    import pandas as pd

logger = logging.get_logger(__name__)


class MethylationSite:
    """Represents a single CpG methylation site."""

    def __init__(self, chromosome: str, position: int, methylated_reads: int, total_reads: int, strand: str = "+"):
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
            "chromosome": self.chromosome,
            "position": self.position,
            "methylated_reads": self.methylated_reads,
            "total_reads": self.total_reads,
            "methylation_level": self.methylation_level,
            "coverage": self.coverage,
            "strand": self.strand,
        }


_BEDGRAPH_LEVEL_SCALES = ("auto", "fraction", "percent")


def load_methylation_bedgraph(
    path: str | Path, min_coverage: int = 1, level_scale: str = "auto"
) -> Dict[str, List[MethylationSite]]:
    """Load methylation data from BEDgraph format.

    Column 4 is the methylation level. BEDgraph carries levels, not read counts, and tools write them either as a
    fraction (0-1) or as a percentage (0-100); ``level_scale`` says which. The default ``"auto"`` reads the file as
    percentages when any level in it is above 1 (a warning is logged) and as fractions otherwise. A line whose
    level is outside the range of the scale (negative, above 1 for fractions, above 100 for percentages, or not a
    number) is skipped with a warning that names the line; it no longer aborts the whole load.

    Args:
        path: Path to BEDgraph file
        min_coverage: Minimum read coverage to include
        level_scale: How to read column 4: ``"fraction"`` (0-1), ``"percent"`` (0-100) or ``"auto"``. Force
            ``"fraction"`` for a fraction file that may hold a stray value above 1.

    Returns:
        Dictionary mapping chromosome names to lists of MethylationSite objects

    Raises:
        ValidationError: If the file doesn't exist
        ValueError: If level_scale is not one of the values above
    """
    if level_scale not in _BEDGRAPH_LEVEL_SCALES:
        raise ValueError(f"level_scale must be one of {', '.join(_BEDGRAPH_LEVEL_SCALES)}, got {level_scale!r}")

    path = validation.validate_path_exists(Path(path))

    logger.info(f"Loading methylation BEDgraph from {path}")

    methylation_data: Dict[str, List[MethylationSite]] = defaultdict(list)

    try:
        records: List[Tuple[int, str, int, float]] = []
        with io.open_text_auto(path, "rt") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue

                parts = line.split("\t")
                if len(parts) < 4:
                    logger.warning(f"Skipping malformed line {line_num}: {line}")
                    continue

                try:
                    chromosome = parts[0]
                    start_pos = int(parts[1])
                    int(parts[2])
                    methylation_level = float(parts[3])
                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue

                records.append((line_num, chromosome, start_pos, methylation_level))

        # The scale is a property of the file, so it is settled once all the levels have been read
        percent = level_scale == "percent"
        if level_scale == "auto":
            highest = max((level for _, _, _, level in records if not math.isnan(level)), default=0.0)
            percent = highest > 1.0
            if percent:
                logger.warning(
                    f"Methylation levels in {path} reach {highest:g}, above 1: reading column 4 as percentages "
                    "(0-100); pass level_scale='fraction' to read it as a fraction"
                )
        upper_bound = 100.0 if percent else 1.0

        for line_num, chromosome, start_pos, methylation_level in records:
            if not 0.0 <= methylation_level <= upper_bound:
                logger.warning(
                    f"Skipping line {line_num}: methylation level {methylation_level:g} is outside "
                    f"0-{upper_bound:g} ({'percent' if percent else 'fraction'} scale)"
                )
                continue

            fraction = methylation_level / 100.0 if percent else methylation_level

            # BEDgraph format typically uses the start position
            # BEDgraph carries a level, not read counts. Estimate
            # counts at a fixed 100-read precision so the level survives
            # the count round trip (int truncation at total_reads=1
            # previously collapsed any level > 0.5 to 0.0).
            total_reads = max(min_coverage, 100)
            methylated_reads = int(round(fraction * total_reads))

            methylation_data[chromosome].append(
                MethylationSite(
                    chromosome=chromosome,
                    position=start_pos,
                    methylated_reads=methylated_reads,
                    total_reads=total_reads,
                )
            )

    except Exception as e:
        logger.error(f"Error loading BEDgraph file {path}: {e}")
        raise errors.IOError(f"Failed to load methylation BEDgraph: {e}") from e

    # Sort sites by position for each chromosome
    for chromosome in methylation_data:
        methylation_data[chromosome].sort(key=lambda x: x.position)

    total_sites = sum(len(sites) for sites in methylation_data.values())
    logger.info(f"Loaded {total_sites} methylation sites across {len(methylation_data)} chromosomes")

    return dict(methylation_data)


def load_methylation_cov(path: str | Path, min_coverage: int = 1) -> Dict[str, List[MethylationSite]]:
    """Load methylation data from a Bismark coverage (``.cov``) file.

    Reads the tab-separated, header-less layout of a Bismark coverage file (``*.bismark.cov``, as written by
    ``bismark2bedGraph``), gzip-compressed or not::

        chromosome  start  end  methylation_percentage  count_methylated  count_unmethylated

    ``start`` becomes the site position (as written in the file), and the coverage of a site is
    ``count_methylated + count_unmethylated``: the line ``chr1  100  100  75.0  15  5`` is 15 methylated reads out
    of 20, a methylation level of 0.75. The methylation level is computed from the counts; the percentage column is
    only checked to be a number. A coverage file records no strand, so every site keeps the default strand. Blank
    lines and lines starting with ``#`` are skipped, and so are malformed ones (logged as warnings), for example
    lines in an older layout that put a strand symbol in the third column.

    Args:
        path: Path to .cov file
        min_coverage: Minimum read coverage (methylated + unmethylated reads) to include

    Returns:
        Dictionary mapping chromosome names to lists of MethylationSite objects
    """
    path = validation.validate_path_exists(Path(path))

    logger.info(f"Loading methylation coverage from {path}")

    methylation_data = defaultdict(list)

    try:
        with io.open_text_auto(path, "rt") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split("\t")
                if len(parts) < 6:
                    logger.warning(f"Skipping malformed line {line_num}: {line}")
                    continue

                try:
                    chromosome = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    float(parts[3])  # Validate percentage column parses as a number
                    methylated_reads = int(parts[4])
                    unmethylated_reads = int(parts[5])

                    if start < 0 or end < start:
                        raise ValueError(f"invalid coordinates start={start}, end={end}")
                    if methylated_reads < 0 or unmethylated_reads < 0:
                        raise ValueError("read counts must not be negative")
                    total_reads = methylated_reads + unmethylated_reads

                    if total_reads >= min_coverage:
                        site = MethylationSite(
                            chromosome=chromosome,
                            position=start,
                            methylated_reads=methylated_reads,
                            total_reads=total_reads,
                        )
                        methylation_data[chromosome].append(site)

                except (ValueError, IndexError) as e:
                    logger.warning(
                        f"Error parsing line {line_num} (expected Bismark .cov columns chromosome, start, end, "
                        f"percentage, methylated, unmethylated): {e}"
                    )
                    continue

    except Exception as e:
        logger.error(f"Error loading coverage file {path}: {e}")
        raise errors.IOError(f"Failed to load methylation coverage: {e}") from e

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

    stats: Dict[str, Any] = {
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
    coverage_bins: Dict[str, int] = defaultdict(int)
    for coverage in coverages:
        if coverage <= 10:
            coverage_bins[str(coverage)] += 1
        elif coverage <= 50:
            coverage_bins["11-50"] += 1
        else:
            coverage_bins["50+"] += 1

    stats["coverage_distribution"] = dict(coverage_bins)

    # Methylation distribution
    methylation_bins: Dict[str, int] = defaultdict(int)
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


def find_differentially_methylated_regions(
    site_data1: Dict[str, List[MethylationSite]],
    site_data2: Dict[str, List[MethylationSite]],
    delta_threshold: float = 0.2,
    p_value_threshold: float = 0.05,
    min_sites: int = 3,
    max_gap: Optional[int] = None,
) -> List[Dict[str, Any]]:
    """Find differentially methylated regions between two conditions.

    A position is significant when its methylation difference is at least ``delta_threshold`` and the two-sided
    p-value of the test that both conditions have the same methylation proportion (from their methylated and
    total read counts) is at most ``p_value_threshold``. The test is a pooled two-proportion z-test with continuity
    correction when every expected count of the methylated/unmethylated by condition table is at least 10, and
    Fisher's exact test otherwise (low coverage, or levels at 0% or 100%), so sites that are fully methylated in one
    condition and unmethylated in the other are tested rather than discarded. Runs of consecutive significant
    positions with at least ``min_sites`` positions become regions.

    Args:
        site_data1: Methylation data for condition 1
        site_data2: Methylation data for condition 2
        delta_threshold: Minimum methylation difference
        p_value_threshold: Maximum p-value for significance
        min_sites: Minimum number of consecutive sites for DMR
        max_gap: Maximum genomic gap (bp) between consecutive DMR sites.
            ``None`` (default) keeps the historical behaviour of chaining
            any consecutive significant sites regardless of distance.

    Returns:
        List of differentially methylated regions, ordered by chromosome name and position. Each is a dict with
        ``chromosome``, ``start``, ``end``, ``sites`` (positions), ``mean_delta``, ``mean_p_value`` (the mean of the
        per-site p-values), ``min_p_value`` (the smallest of them) and ``direction``.
    """
    logger.info("Finding differentially methylated regions")

    dmr_list = []

    # Get all chromosomes present in both datasets, in a fixed order (set order varies from run to run)
    chromosomes = sorted(set(site_data1.keys()) & set(site_data2.keys()))

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
        current_dmr: Optional[Dict[str, Any]] = None

        for i, pos in enumerate(positions):
            site1 = pos_to_site1[pos]
            site2 = pos_to_site2[pos]

            delta = abs(site1.methylation_level - site2.methylation_level)

            # The test is only needed where the difference is large enough to matter
            p_value = 1.0
            significant = delta >= delta_threshold
            if significant:
                p_value = _two_proportion_p_value(
                    site1.methylated_reads, site1.coverage, site2.methylated_reads, site2.coverage
                )
                significant = p_value <= p_value_threshold

            if significant:
                gap_ok = max_gap is None or current_dmr is None or pos - current_dmr["end"] <= max_gap
                if current_dmr is not None and not gap_ok:
                    # Flush the current DMR before starting a new one
                    if len(current_dmr["sites"]) >= min_sites:
                        dmr_list.append(_finish_dmr(current_dmr))
                    current_dmr = None
                if current_dmr is None:
                    # Start new DMR
                    current_dmr = {
                        "chromosome": chromosome,
                        "start": pos,
                        "end": pos,
                        "sites": [pos],
                        "mean_delta": delta,
                        "mean_p_value": p_value,
                        "min_p_value": p_value,
                        "direction": "hyper" if site1.methylation_level > site2.methylation_level else "hypo",
                        "_p_values": [p_value],
                    }
                else:
                    # Extend current DMR
                    current_dmr["end"] = pos
                    current_dmr["sites"].append(pos)
                    current_dmr["mean_delta"] = (
                        current_dmr["mean_delta"] * (len(current_dmr["sites"]) - 1) + delta
                    ) / len(current_dmr["sites"])
                    current_dmr["_p_values"].append(p_value)
            else:
                # End current DMR if it meets minimum size
                if current_dmr and len(current_dmr["sites"]) >= min_sites:
                    dmr_list.append(_finish_dmr(current_dmr))
                current_dmr = None

        # Don't forget the last DMR
        if current_dmr and len(current_dmr["sites"]) >= min_sites:
            dmr_list.append(_finish_dmr(current_dmr))

    logger.info(f"Found {len(dmr_list)} differentially methylated regions")
    return dmr_list


# Smallest expected cell count of the 2x2 table above which the normal approximation is used instead of the exact test
_NORMAL_APPROXIMATION_MIN_EXPECTED = 10.0


def _finish_dmr(dmr: Dict[str, Any]) -> Dict[str, Any]:
    """Replace a DMR's running list of per-site p-values by their mean and minimum."""
    p_values = dmr.pop("_p_values")
    dmr["mean_p_value"] = sum(p_values) / len(p_values)
    dmr["min_p_value"] = min(p_values)
    return dmr


def _log_choose(n: int, k: int) -> float:
    """Natural logarithm of the binomial coefficient ``n`` choose ``k``."""
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _fisher_exact_two_sided(methylated1: int, total1: int, methylated2: int, total2: int) -> float:
    """Two-sided Fisher exact p-value that two samples share one methylation proportion.

    The 2x2 table is methylated and unmethylated read counts by sample. The p-value is the total probability, under
    the hypergeometric distribution with the margins fixed, of the tables that are no more likely than the observed
    one (with a relative tolerance of 1e-7 for ties, as R's ``fisher.test``).

    Args:
        methylated1: Methylated reads in sample 1.
        total1: Total reads in sample 1.
        methylated2: Methylated reads in sample 2.
        total2: Total reads in sample 2.

    Returns:
        The p-value, in [0, 1].
    """
    methylated = methylated1 + methylated2
    log_total_ways = _log_choose(total1 + total2, methylated)

    def log_probability(methylated_in_1: int) -> float:
        return _log_choose(total1, methylated_in_1) + _log_choose(total2, methylated - methylated_in_1) - log_total_ways

    cutoff = log_probability(methylated1) + 1e-7
    p_value = 0.0
    for methylated_in_1 in range(max(0, methylated - total2), min(total1, methylated) + 1):
        log_p = log_probability(methylated_in_1)
        if log_p <= cutoff:
            p_value += math.exp(log_p)
    return min(1.0, p_value)


def _two_proportion_p_value(methylated1: int, total1: int, methylated2: int, total2: int) -> float:
    """Two-sided p-value that two samples share one methylation proportion, from read counts.

    Uses a continuity-corrected z-test with the pooled proportion when every expected cell count of the
    methylated/unmethylated by sample table is at least 10, and Fisher's exact test otherwise. The pooled standard
    error is the one under the null
    hypothesis, so it is zero only when every read in both samples agrees (no methylated reads at all, or no
    unmethylated ones), where the samples cannot differ and the p-value is 1. The standard error built from each
    sample's own proportion is zero as soon as both proportions are 0 or 1, however different they are, which is
    why it must not be used for the test.

    Args:
        methylated1: Methylated reads in sample 1.
        total1: Total reads in sample 1.
        methylated2: Methylated reads in sample 2.
        total2: Total reads in sample 2.

    Returns:
        The p-value, in [0, 1]; 1.0 when either sample has no reads.
    """
    if total1 <= 0 or total2 <= 0:
        return 1.0

    methylated = methylated1 + methylated2
    total = total1 + total2
    unmethylated = total - methylated
    if methylated == 0 or unmethylated == 0:
        return 1.0

    # Smallest cell expectation of the 2x2 table: smallest row total x smallest column total / table total
    if min(total1, total2) * min(methylated, unmethylated) / total >= _NORMAL_APPROXIMATION_MIN_EXPECTED:
        pooled = methylated / total
        standard_error = math.sqrt(pooled * (1.0 - pooled) * (1.0 / total1 + 1.0 / total2))
        # Yates' continuity correction: read counts are discrete, and without it the z-test is liberal next to
        # the exact test it hands over to at the smallest expected counts
        difference = abs(methylated1 / total1 - methylated2 / total2) - 0.5 * (1.0 / total1 + 1.0 / total2)
        z_score = max(0.0, difference) / standard_error
        return math.erfc(z_score / math.sqrt(2.0))

    return _fisher_exact_two_sided(methylated1, total1, methylated2, total2)


def _cpg_statistics(sequence: str) -> Tuple[float, float, int]:
    """GC fraction, observed/expected CpG ratio and CpG count of an upper-case DNA sequence.

    Observed/expected CpG is ``CpG count * length / (C count * G count)`` (Gardiner-Garden & Frommer, 1987),
    taken as 0.0 when the sequence has no C or no G.
    """
    length = len(sequence)
    c_count = sequence.count("C")
    g_count = sequence.count("G")
    cpg_count = sequence.count("CG")
    gc_content = (c_count + g_count) / length if length else 0.0
    cpg_ratio = cpg_count * length / (c_count * g_count) if c_count and g_count else 0.0
    return gc_content, cpg_ratio, cpg_count


def find_cpg_islands(
    sequences: Mapping[str, str],
    window_size: int = 200,
    step_size: int = 50,
    min_gc: float = 0.5,
    min_cpg_ratio: float = 0.6,
    min_length: int = 200,
) -> List[Dict[str, Any]]:
    """Find CpG islands in DNA sequence with the Gardiner-Garden & Frommer criteria.

    A CpG island is a stretch of at least ``min_length`` bp (200) whose GC fraction is at least ``min_gc`` (0.5)
    and whose observed/expected CpG ratio, ``CpG count * length / (C count * G count)``, is at least
    ``min_cpg_ratio`` (0.6). Each chromosome is scanned with windows of ``window_size`` bp every ``step_size`` bp
    (the last window is aligned to the end of the sequence); windows that meet the GC and CpG criteria are merged
    when they overlap or touch, and a merged stretch is reported when it is at least ``min_length`` long and still
    meets both criteria as a whole. Island boundaries are therefore resolved to about ``window_size`` bp (a window
    that is partly island still qualifies), finer with a smaller ``window_size`` and ``step_size``. Sequence is
    upper-cased, and bases other than C and G (including N) count as neither.

    Args:
        sequences: Mapping of chromosome name to its DNA sequence.
        window_size: Scanning window size in bp.
        step_size: Distance between window starts in bp.
        min_gc: Minimum GC fraction.
        min_cpg_ratio: Minimum observed/expected CpG ratio.
        min_length: Minimum island length in bp.

    Returns:
        Islands in the order of ``sequences`` and then by start. Each is a dict with ``chromosome``, ``start`` and
        ``end`` (0-based, end exclusive), ``length``, ``gc_content``, ``cpg_ratio`` and ``cpg_count``.

    Raises:
        ValueError: If ``window_size`` or ``step_size`` is below 1.
    """
    if window_size < 1 or step_size < 1:
        raise ValueError(f"window_size and step_size must be at least 1, got {window_size} and {step_size}")

    islands: List[Dict[str, Any]] = []

    for chromosome, raw_sequence in sequences.items():
        sequence = raw_sequence.upper()
        if len(sequence) < window_size:
            continue

        starts = list(range(0, len(sequence) - window_size + 1, step_size))
        if starts[-1] != len(sequence) - window_size:
            starts.append(len(sequence) - window_size)

        # Merge qualifying windows that overlap or touch into candidate stretches
        stretches: List[List[int]] = []
        for start in starts:
            gc_content, cpg_ratio, _ = _cpg_statistics(sequence[start : start + window_size])
            if gc_content < min_gc or cpg_ratio < min_cpg_ratio:
                continue
            if stretches and start <= stretches[-1][1]:
                stretches[-1][1] = start + window_size
            else:
                stretches.append([start, start + window_size])

        for start, end in stretches:
            if end - start < min_length:
                continue
            gc_content, cpg_ratio, cpg_count = _cpg_statistics(sequence[start:end])
            if gc_content >= min_gc and cpg_ratio >= min_cpg_ratio:
                islands.append(
                    {
                        "chromosome": chromosome,
                        "start": start,
                        "end": end,
                        "length": end - start,
                        "gc_content": gc_content,
                        "cpg_ratio": cpg_ratio,
                        "cpg_count": cpg_count,
                    }
                )

    return islands


def identify_cpg_islands(
    methylation_data: Dict[str, List[MethylationSite]],
    window_size: int = 200,
    step_size: int = 50,
    min_gc: float = 0.5,
    min_cpg_ratio: float = 0.6,
    *,
    sequences: Optional[Mapping[str, str]] = None,
    min_length: int = 200,
) -> List[Dict[str, Any]]:
    """Identify CpG islands in DNA sequence and summarize the methylation sites inside them.

    Islands are called from the sequence with :func:`find_cpg_islands` (Gardiner-Garden & Frommer: at least
    ``min_length`` bp, GC fraction at least ``min_gc``, observed/expected CpG at least ``min_cpg_ratio``). Methylation
    levels cannot stand in for sequence: GC content and CpG depletion are properties of the DNA, so the sequence is
    required. The methylation sites are only summarized per island.

    Args:
        methylation_data: Dictionary of methylation sites by chromosome (may be empty or lack chromosomes: the
            islands are then reported with no sites)
        window_size: Size of sliding window in base pairs
        step_size: Step size for sliding window
        min_gc: Minimum GC content for CpG island
        min_cpg_ratio: Minimum CpG ratio (observed/expected)
        sequences: Mapping of chromosome name to DNA sequence (required, keyword-only). Site positions are
            compared with the 0-based sequence coordinates, so shift 1-based positions (such as Bismark's) first
            if the difference matters.
        min_length: Minimum island length in base pairs

    Returns:
        List of identified CpG islands: the dicts of :func:`find_cpg_islands` plus ``sites`` (number of methylation
        sites with ``start <= position < end``) and ``avg_methylation`` (their mean level, ``None`` without sites)

    Raises:
        ValueError: If ``sequences`` is not given, or ``window_size`` or ``step_size`` is below 1.
    """
    if sequences is None:
        raise ValueError(
            "identify_cpg_islands needs the DNA sequence: pass sequences={chromosome: sequence}. CpG islands are "
            "defined by the GC content and CpG depletion of the sequence, which methylation levels cannot supply."
        )

    logger.info("Identifying CpG islands")

    islands = find_cpg_islands(
        sequences,
        window_size=window_size,
        step_size=step_size,
        min_gc=min_gc,
        min_cpg_ratio=min_cpg_ratio,
        min_length=min_length,
    )

    sites_by_chromosome: Dict[str, Tuple[List[int], List[float]]] = {}
    for island in islands:
        chromosome = island["chromosome"]
        if chromosome not in sites_by_chromosome:
            ordered = sorted(methylation_data.get(chromosome, []), key=lambda site: site.position)
            sites_by_chromosome[chromosome] = (
                [site.position for site in ordered],
                [site.methylation_level for site in ordered],
            )
        positions, levels = sites_by_chromosome[chromosome]
        inside = levels[bisect_left(positions, island["start"]) : bisect_left(positions, island["end"])]
        island["sites"] = len(inside)
        island["avg_methylation"] = sum(inside) / len(inside) if inside else None

    logger.info(f"Identified {len(islands)} CpG islands")
    return islands


def calculate_methylation_entropy(
    methylation_data: Dict[str, List[MethylationSite]], window_size: int = 1000
) -> Dict[str, List[Tuple[int, float]]]:
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

        sites = sorted(sites, key=lambda x: x.position)

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
                    # The last bin is closed on the right, so a fully methylated site (1.0) is counted
                    in_last_bin = i == len(bins) - 2 and level == bins[i + 1]
                    if bins[i] <= level < bins[i + 1] or in_last_bin:
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


def export_methylation_bedgraph(methylation_data: Dict[str, List[MethylationSite]], output_path: str | Path) -> None:
    """Export methylation data to BEDgraph format.

    Args:
        output_path: Path for output BEDgraph file
        methylation_data: Dictionary of methylation sites by chromosome
    """
    output_path = Path(output_path)

    logger.info(f"Exporting methylation data to BEDgraph: {output_path}")

    with open(output_path, "w") as f:
        # Write track header
        f.write(
            'track type=bedGraph name="Methylation" '
            'description="DNA Methylation Levels" visibility=full '
            "color=200,100,0\n"
        )

        for chromosome in sorted(methylation_data.keys()):
            sites = methylation_data[chromosome]
            sites = sorted(sites, key=lambda x: x.position)

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

    df = pd.read_csv(path, sep="\t")
    required_cols = ["chrom", "pos", "methylated", "unmethylated"]

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
    if "methylated" not in df.columns or "unmethylated" not in df.columns:
        raise ValueError("DataFrame must contain 'methylated' and 'unmethylated' columns")

    df_copy = df.copy()
    total_reads = df_copy["methylated"] + df_copy["unmethylated"]
    df_copy["beta"] = df_copy["methylated"] / total_reads
    df_copy["beta"] = df_copy["beta"].fillna(0.0)  # Handle division by zero

    return df_copy


def summarize_beta_by_chromosome(df: pd.DataFrame) -> pd.DataFrame:
    """Summarize beta values by chromosome.

    Args:
        df: DataFrame with 'chrom' and 'beta' columns

    Returns:
        DataFrame with summary statistics by chromosome
    """
    if "chrom" not in df.columns or "beta" not in df.columns:
        raise ValueError("DataFrame must contain 'chrom' and 'beta' columns")

    return df.groupby("chrom")["beta"].agg(["mean", "std", "min", "max", "count"])


def generate_methylation_report(
    methylation_data: Dict[str, List[MethylationSite]], output_path: Optional[str | Path] = None
) -> str:
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
        coverage_dist = stats.get("coverage_distribution", {})
        if coverage_dist:
            report_lines.append("Coverage Distribution:")
            for cov_range, count in sorted(coverage_dist.items()):
                report_lines.append(f"  {cov_range}x: {count:,} sites")
            report_lines.append("")

        # Top chromosomes
        chr_stats = stats.get("chromosome_stats", {})
        if chr_stats:
            report_lines.append("Per-Chromosome Statistics (Top 5):")
            sorted_chrs = sorted(chr_stats.items(), key=lambda x: x[1]["sites"], reverse=True)[:5]
            for chr_name, chr_data in sorted_chrs:
                report_lines.append(
                    f"  {chr_name}: {chr_data['sites']:,} sites, {chr_data['mean_methylation']:.3f} mean methylation"
                )
            report_lines.append("")

    report = "\n".join(report_lines)

    if output_path:
        output_path = Path(output_path)
        with open(output_path, "w") as f:
            f.write(report)
        logger.info(f"Methylation report saved to {output_path}")

    return report
