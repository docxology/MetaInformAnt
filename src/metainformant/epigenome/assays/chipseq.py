"""ChIP-seq analysis and processing.

This module provides comprehensive tools for analyzing ChIP-seq data,
including peak calling, motif analysis, quality control, and data integration.
"""

from __future__ import annotations

import math
import re
import statistics
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

from metainformant.core import io
from metainformant.core.data import validation
from metainformant.core.sequence import reverse_complement
from metainformant.core.utils import errors, logging
from metainformant.epigenome.assays._peak_format import (
    MIN_COLUMNS,
    format_peak_line,
    is_skipped_line,
    normalise_peak_format,
    parse_peak_fields,
)

logger = logging.get_logger(__name__)


class ChIPPeak:
    """Represents a ChIP-seq peak.

    ``p_value`` and ``q_value`` are probabilities (between 0 and 1, or ``None``). Peak files store them as
    ``-log10`` values, as the ENCODE specification and MACS2 do; :func:`load_chip_peaks` and
    :func:`save_chip_peaks` convert at the file boundary, so a q-value of 0.01 is the column value 2.
    """

    def __init__(
        self,
        chromosome: str,
        start: int,
        end: int,
        summit: Optional[int] = None,
        score: float = 0.0,
        strand: str = ".",
        signal_value: float = 0.0,
        p_value: Optional[float] = None,
        q_value: Optional[float] = None,
    ):
        """Initialize a ChIP-seq peak.

        Args:
            chromosome: Chromosome name
            start: Peak start position
            end: Peak end position
            summit: Peak summit position
            score: Peak score
            strand: DNA strand
            signal_value: Signal value
            p_value: P-value, as a probability between 0 and 1
            q_value: Q-value (FDR), as a probability between 0 and 1
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

    def to_bed_format(self, format: str = "narrowpeak") -> str:
        """Convert to one line of a standard peak file.

        The columns follow the ENCODE layouts that :func:`load_chip_peaks` reads: ``narrowpeak`` (default,
        ten columns ending in the summit as an offset from ``start``), ``broadpeak`` (nine columns) or
        ``bed`` (six columns). The p-value and q-value are probabilities and are written as ``-log10`` values,
        as the specification requires (``1e-05`` is written ``5``); a probability of exactly 0 is written as the
        ``-log10`` of the smallest positive float (about 323.3). A missing p-value, q-value or summit is written
        as ``-1``.

        Args:
            format: Peak layout to write.

        Returns:
            The tab-separated line, without a trailing newline.

        Raises:
            ValueError: If the layout is not supported, or the p-value or q-value is not between 0 and 1
                (only the ``narrowpeak`` and ``broadpeak`` layouts write them).
        """
        return format_peak_line(
            name=f"peak_{self.start}_{self.end}",
            chromosome=self.chromosome,
            start=self.start,
            end=self.end,
            score=self.score,
            strand=self.strand,
            signal_value=self.signal_value,
            p_value=self.p_value,
            q_value=self.q_value,
            summit=self.summit,
            format=format,
        )

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
            "chromosome": self.chromosome,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "summit": self.summit,
            "score": self.score,
            "strand": self.strand,
            "signal_value": self.signal_value,
            "p_value": self.p_value,
            "q_value": self.q_value,
        }


def load_chip_peaks(path: str | Path, format: str = "narrowpeak") -> List[ChIPPeak]:
    """Load ChIP-seq peaks from a file.

    The columns are read in the ENCODE layouts that :func:`save_chip_peaks` writes: ``narrowpeak`` has ten
    columns (``chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak``, with ``peak`` the
    summit as an offset from ``start``), ``broadpeak`` the first nine and ``bed`` the first six (only the first
    three are required). ``-1`` or ``.`` in the p-value, q-value or summit column means "not available" and
    gives ``None``.

    The pValue and qValue columns hold ``-log10`` values, as the ENCODE specification and MACS2 write them, and
    the peaks carry probabilities: a column value ``v`` becomes ``p_value`` or ``q_value`` ``10 ** -v``, so a
    column of 2 is a q-value of 0.01 and a MACS2 value of 9.1 is about 7.9e-10. A column above about 323.6
    underflows to a probability of 0.0. Any negative column value other than ``-1`` (or ``nan``) is not a
    ``-log10`` probability, and the line is skipped with a warning.

    Blank, ``#``, ``track`` and ``browser`` lines are ignored. A line with too few columns for the layout, or
    with a value that does not parse, is skipped with a warning.

    Args:
        path: Path to peak file (optionally gzip-compressed)
        format: Peak file layout ("narrowpeak", "broadpeak", "bed"), in any letter case

    Returns:
        List of ChIPPeak objects

    Raises:
        ValidationError: If the file doesn't exist
        ValueError: If the format is not supported
    """
    fmt = normalise_peak_format(format)
    path = validation.validate_path_exists(Path(path))

    logger.info(f"Loading ChIP-seq peaks from {path} (format: {fmt})")

    peaks = []
    min_columns = MIN_COLUMNS[fmt]

    try:
        with io.open_text_auto(path, "rt") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if is_skipped_line(line):
                    continue

                parts = line.split("\t")
                if len(parts) < min_columns:
                    logger.warning(
                        f"Skipping line {line_num}: {fmt} needs at least {min_columns} columns, found {len(parts)}"
                    )
                    continue

                try:
                    fields = parse_peak_fields(parts, fmt)
                    peak = ChIPPeak(
                        chromosome=fields.chromosome,
                        start=fields.start,
                        end=fields.end,
                        summit=fields.summit,
                        score=fields.score,
                        strand=fields.strand,
                        signal_value=fields.signal_value or 0.0,
                        p_value=fields.p_value,
                        q_value=fields.q_value,
                    )
                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue

                peaks.append(peak)

    except Exception as e:
        logger.error(f"Error loading ChIP-seq peaks from {path}: {e}")
        raise errors.IOError(f"Failed to load ChIP-seq peaks: {e}") from e

    logger.info(f"Loaded {len(peaks)} ChIP-seq peaks")
    return peaks


def save_chip_peaks(peaks: List[ChIPPeak], path: str | Path, format: str = "narrowpeak") -> None:
    """Save ChIP-seq peaks to a file that :func:`load_chip_peaks` reads back.

    The columns follow the ENCODE layouts (see :meth:`ChIPPeak.to_bed_format`): ten columns for ``narrowpeak``,
    nine for ``broadpeak`` and six for ``bed``. ``narrowpeak`` and ``broadpeak`` files start with a ``#``
    comment line, which the loaders skip. The name column (``peak_<start>_<end>``) is regenerated on every
    write; peaks carry no name. A ``broadpeak`` or ``bed`` file does not hold the summit, and a ``bed`` file
    holds neither the signal value nor the p- and q-values.

    The peaks' p-values and q-values are probabilities and are written as ``-log10`` values, as the ENCODE
    specification requires (a q-value of 0.01 is written ``2``); a probability of exactly 0 is written as the
    ``-log10`` of the smallest positive float (about 323.3), so it reads back as that float rather than 0.0.
    Probabilities lose a few units in the last place in the round trip, and values below about 1e-308 lose
    precision as the float underflows.

    Every peak is formatted before the file is opened, so a failure leaves no partial file and an existing file
    is not overwritten.

    Args:
        peaks: List of ChIPPeak objects
        path: Output file path
        format: Output layout ("narrowpeak", "broadpeak", "bed"), in any letter case

    Raises:
        ValueError: If the format is not supported, or a peak's p-value or q-value is not between 0 and 1
            (nothing is written)
    """
    fmt = normalise_peak_format(format)
    path = Path(path)

    # Format every peak first: an invalid p- or q-value raises ValueError here, before the file is touched
    lines = [peak.to_bed_format(fmt) for peak in peaks]

    logger.info(f"Saving {len(peaks)} ChIP-seq peaks to {path} (format: {fmt})")

    with open(path, "w") as f:
        if fmt != "bed":
            f.write("# Generated by MetaInformant ChIP-seq analysis\n")

        for line in lines:
            f.write(line + "\n")

    logger.info(f"Saved {len(peaks)} peaks to {path}")


def filter_peaks_by_score(peaks: List[ChIPPeak], min_score: float, max_peaks: Optional[int] = None) -> List[ChIPPeak]:
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

    stats: Dict[str, Any] = {
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
        stats.update(
            {
                "mean_signal": statistics.mean(signal_values),
                "median_signal": statistics.median(signal_values),
                "signal_std": statistics.stdev(signal_values) if len(signal_values) > 1 else 0,
            }
        )

    # Peak length distribution
    length_bins: Dict[str, int] = defaultdict(int)
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
    score_bins: Dict[str, int] = defaultdict(int)
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
    chr_counts: Dict[str, int] = defaultdict(int)
    for peak in peaks:
        chr_counts[peak.chromosome] += 1

    stats["chromosome_distribution"] = dict(sorted(chr_counts.items()))

    return stats


def find_overlapping_peaks(
    peaks1: List[ChIPPeak], peaks2: List[ChIPPeak], min_overlap: int = 1
) -> List[Tuple[ChIPPeak, ChIPPeak]]:
    """Find every pair of overlapping peaks between two peak sets.

    A pair ``(peak1, peak2)`` is returned when both peaks are on the same chromosome and
    ``peak1.overlaps_with(peak2, min_overlap)`` holds. Peaks within one set may overlap each other (or be
    nested), and a peak that overlaps several peaks of the other set appears in several pairs. The result is
    ordered by chromosome (in order of first appearance in ``peaks1``), then by ``peak1`` start, then by
    ``peak2`` start, so it does not depend on hash seeds.

    Args:
        peaks1: First set of peaks
        peaks2: Second set of peaks
        min_overlap: Minimum overlap required, in bases

    Returns:
        List of (peak1, peak2) tuples for overlapping peaks
    """
    logger.info("Finding overlapping peaks")

    # Group peaks by chromosome for efficiency
    chr_peaks1: Dict[str, List[ChIPPeak]] = defaultdict(list)
    chr_peaks2: Dict[str, List[ChIPPeak]] = defaultdict(list)

    for peak in peaks1:
        chr_peaks1[peak.chromosome].append(peak)

    for peak in peaks2:
        chr_peaks2[peak.chromosome].append(peak)

    overlapping_pairs: List[Tuple[ChIPPeak, ChIPPeak]] = []

    for chromosome, peaks1_chr in chr_peaks1.items():
        peaks2_chr = chr_peaks2.get(chromosome)
        if not peaks2_chr:
            continue

        # Sort peaks2 by start and keep the running maximum of the ends. A peak2 can overlap a given peak1 only if
        # it starts at or before peak1.end - min_overlap and some peak2 up to it ends at or after
        # peak1.start + min_overlap, so two bisections bound the candidates; each is then checked exactly.
        peaks2_chr.sort(key=lambda p: (p.start, p.end))
        starts = [p.start for p in peaks2_chr]
        running_end: List[int] = []
        highest_end = peaks2_chr[0].end
        for peak in peaks2_chr:
            highest_end = max(highest_end, peak.end)
            running_end.append(highest_end)

        peaks1_chr.sort(key=lambda p: (p.start, p.end))
        for peak1 in peaks1_chr:
            first = bisect_left(running_end, peak1.start + min_overlap)
            last = bisect_right(starts, peak1.end - min_overlap)
            for peak2 in peaks2_chr[first:last]:
                if peak1.overlaps_with(peak2, min_overlap):
                    overlapping_pairs.append((peak1, peak2))

    logger.info(f"Found {len(overlapping_pairs)} overlapping peak pairs")
    return overlapping_pairs


def _lowest_probability(first: Optional[float], second: Optional[float]) -> Optional[float]:
    """The lower of two probabilities, or ``None`` when either is missing.

    A probability of 0.0 is a value, not a gap: peak files hold ``-log10`` values, so the strongest peaks (a column
    above about 323) are read as exactly 0.0, and a truth test on them would discard the best significance.
    """
    return min(first, second) if first is not None and second is not None else None


def merge_overlapping_peaks(peaks: List[ChIPPeak], max_distance: int = 0) -> List[ChIPPeak]:
    """Merge overlapping or nearby peaks.

    A merged peak takes the lower p-value and q-value of the two peaks it joins (``None`` when either has none).

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
                    p_value=_lowest_probability(current_peak.p_value, peak.p_value),
                    q_value=_lowest_probability(current_peak.q_value, peak.q_value),
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


def _poisson_upper_tail(observed: int, expected: float) -> Tuple[float, float]:
    """Return ``P(X >= observed)`` for ``X ~ Poisson(expected)`` as ``(p, log10(p))``.

    The sum is carried out on the probability mass function evaluated in log space, so it cannot overflow for
    large counts and does not lose a tiny tail to the cancellation in ``1 - P(X < observed)``. Above the mean the
    tail is summed directly from the mass at ``observed`` (the terms fall off geometrically); at or below the mean
    the tail is at least about one half and is taken as ``1 - P(X < observed)``, where nothing cancels. ``p``
    loses precision below about 1e-308 and underflows to 0.0 below about 5e-324, while ``log10(p)`` stays finite.

    Args:
        observed: Observed count.
        expected: Poisson mean (the expected count).

    Returns:
        The upper-tail probability and its base-10 logarithm.
    """
    if observed <= 0 or not math.isfinite(expected):
        return 1.0, 0.0
    if expected <= 0:
        return 0.0, -math.inf

    log_expected = math.log(expected)

    def log_pmf(count: int) -> float:
        return count * log_expected - expected - math.lgamma(count + 1)

    if observed > expected:
        # pmf(k + 1) = pmf(k) * expected / (k + 1), a ratio below one from the mean upward
        term = 1.0
        series = 1.0
        count = observed
        while term >= series * 1e-17:
            count += 1
            term *= expected / count
            series += term
        log_p = log_pmf(observed) + math.log(series)
        return math.exp(log_p), log_p / math.log(10.0)

    # observed <= expected: add the lower tail P(X <= observed - 1) from its largest term downwards,
    # pmf(k - 1) = pmf(k) * k / expected, a ratio below one under the mean
    term = math.exp(log_pmf(observed - 1))
    lower = term
    count = observed - 1
    while count > 0 and term > lower * 1e-17:
        term *= count / expected
        lower += term
        count -= 1
    p_value = min(1.0, max(0.0, 1.0 - lower))
    return p_value, (math.log10(p_value) if p_value > 0.0 else -math.inf)


def calculate_peak_enrichment(
    peaks: List[ChIPPeak], genome_size: int, expected_peaks: Optional[float] = None
) -> Dict[str, Any]:
    """Calculate peak enrichment statistics.

    When ``expected_peaks`` is positive the result also holds the observed/expected ratio and the Poisson
    probability of seeing at least the observed number of peaks (``enrichment_p_value``), together with its
    base-10 logarithm (``enrichment_log10_p_value``), which stays finite where the probability itself underflows
    to 0.0. The probability is computed in log space, so any number of peaks is handled.

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
        stats.update(
            {
                "expected_peaks": expected_peaks,
                "enrichment_ratio": enrichment_ratio,
                "fold_enrichment": enrichment_ratio,
            }
        )

        # Poisson probability of at least the observed number of peaks, computed in log space
        if expected_peaks > 0:
            p_value, log10_p_value = _poisson_upper_tail(observed_peaks, float(expected_peaks))
            stats["enrichment_p_value"] = p_value
            stats["enrichment_log10_p_value"] = log10_p_value

    return stats


# IUPAC nucleotide code table. Kept local because the epigenome domain must
# not import from the dna domain (scripts/quality/check_module_boundaries.py).
_IUPAC_BASES: Dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}

_IUPAC_COMPLEMENT: Dict[str, str] = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
    "N": "N",
}


def _iupac_pattern_to_regex(pattern: str) -> "re.Pattern[str]":
    """Compile an IUPAC nucleotide pattern into a regular expression.

    Args:
        pattern: Motif pattern using IUPAC nucleotide codes.

    Returns:
        Compiled case-insensitive regular expression.

    Raises:
        ValueError: If the pattern is empty or contains a non-IUPAC character.
    """
    if not pattern:
        # An empty regular expression matches at every position, which would report a site everywhere
        raise ValueError("Motif pattern must not be empty")
    parts: List[str] = []
    for char in pattern.upper():
        bases = _IUPAC_BASES.get(char)
        if bases is None:
            raise ValueError(f"Invalid IUPAC motif character {char!r} in pattern {pattern!r}")
        parts.append(f"[{bases}]" if len(bases) > 1 else bases)
    return re.compile("".join(parts))


def _reverse_complement_iupac(pattern: str) -> str:
    """Return the reverse complement of an IUPAC motif pattern."""
    return "".join(_IUPAC_COMPLEMENT[char] for char in reversed(pattern.upper()))


def _iter_matches(regex: "re.Pattern[str]", text: str) -> Iterator[Tuple[int, str]]:
    """Yield ``(offset, matched text)`` for every match of ``regex`` in ``text``, overlapping ones included.

    ``re.finditer`` resumes after the end of a match, so it skips an occurrence that starts inside the previous
    one (``ATAT`` in ``ATATAT`` is found at 0 and 2, not only at 0). Resuming one base after the start of each
    match reports them all. ``regex`` must not be able to match the empty string.
    """
    position = 0
    while True:
        match = regex.search(text, position)
        if match is None:
            return
        yield match.start(), match.group(0)
        position = match.start() + 1


def _read_fasta_sequences(path: str | Path) -> Dict[str, str]:
    """Read a FASTA file into a mapping of sequence ID to sequence.

    The sequence ID is the first whitespace-delimited word of the header,
    the same convention as the canonical FASTA reader in the dna domain.

    Args:
        path: Path to the (optionally gzipped) FASTA file.

    Returns:
        Dict mapping sequence IDs to upper-cased sequence strings.
    """
    sequences: Dict[str, str] = {}
    current_id: Optional[str] = None
    chunks: List[str] = []

    with io.open_text_auto(Path(path)) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(chunks)
                current_id = line[1:].split()[0]
                chunks = []
            elif current_id is not None:
                chunks.append(line.upper())

    if current_id is not None:
        sequences[current_id] = "".join(chunks)

    return sequences


def find_motifs_in_peaks(
    peaks: List[ChIPPeak],
    genome_fasta: str | Path,
    motif_patterns: List[str],
    window_size: int = 200,
) -> Dict[str, Any]:
    """Find motif occurrences in peak regions.

    Loads the genome FASTA once, then for every peak scans the
    ``window_size`` window centered on the peak summit (clipped to the peak
    boundaries) for each IUPAC motif pattern on both strands. Every occurrence
    is reported, including ones that overlap an earlier occurrence of the same
    pattern. Peaks whose chromosome is absent from the FASTA and windows
    without any match contribute no occurrences -- results are never
    fabricated.

    Args:
        peaks: List of ChIPPeak objects
        genome_fasta: Path to genome FASTA file
        motif_patterns: IUPAC motif patterns to search (e.g. ``"TTGACAW"``)
        window_size: Size of the window around the peak summit to search

    Returns:
        Dictionary with motif scanning results:

        - ``total_peaks_analyzed``: number of input peaks
        - ``motif_counts``: per pattern, the number of peaks with >= 1 match
        - ``match_counts``: per pattern, the total number of matches
        - ``motif_positions``: per pattern, a list of match dicts containing
          ``chromosome``, ``position`` (0-based genomic start of the match),
          ``strand`` ('+' or '-'), ``sequence`` (the matched sequence),
          ``peak_start``, ``peak_end``, and ``peak_score``
    """
    logger.info(f"Scanning {len(peaks)} peaks for {len(motif_patterns)} motif patterns")

    regexes = {pattern: _iupac_pattern_to_regex(pattern) for pattern in motif_patterns}
    rc_regexes = {pattern: _iupac_pattern_to_regex(_reverse_complement_iupac(pattern)) for pattern in motif_patterns}

    genome: Dict[str, str] = _read_fasta_sequences(genome_fasta) if peaks else {}

    motif_counts: Dict[str, int] = defaultdict(int)
    match_counts: Dict[str, int] = defaultdict(int)
    motif_positions: Dict[str, List[Dict[str, Any]]] = defaultdict(list)

    half_window = max(window_size, 0) // 2

    for peak in peaks:
        chromosome = genome.get(peak.chromosome)
        if chromosome is None:
            continue  # Chromosome missing from the FASTA: no matches, honestly

        summit = peak.summit if peak.summit is not None else peak.center
        window_start = max(peak.start, summit - half_window)
        window_end = min(peak.end, summit + half_window)
        if window_end - window_start < 1:
            continue

        window = chromosome[window_start:window_end].upper()

        for pattern in motif_patterns:
            matches: List[Dict[str, Any]] = []

            for offset, matched in _iter_matches(regexes[pattern], window):
                matches.append(
                    {
                        "chromosome": peak.chromosome,
                        "position": window_start + offset,
                        "strand": "+",
                        "sequence": matched,
                        "peak_start": peak.start,
                        "peak_end": peak.end,
                        "peak_score": peak.score,
                    }
                )

            # Minus-strand occurrences: an interval of the forward window
            # whose reverse complement matches the motif. Searching the
            # forward window with the reverse-complemented pattern keeps
            # positions in forward genomic coordinates.
            for offset, matched in _iter_matches(rc_regexes[pattern], window):
                matches.append(
                    {
                        "chromosome": peak.chromosome,
                        "position": window_start + offset,
                        "strand": "-",
                        "sequence": reverse_complement(matched),
                        "peak_start": peak.start,
                        "peak_end": peak.end,
                        "peak_score": peak.score,
                    }
                )

            if matches:
                motif_counts[pattern] += 1
                match_counts[pattern] += len(matches)
                motif_positions[pattern].extend(matches)

    results: Dict[str, Any] = {
        "total_peaks_analyzed": len(peaks),
        "motif_counts": dict(motif_counts),
        "match_counts": dict(match_counts),
        "motif_positions": dict(motif_positions),
    }

    logger.info(f"Motif scan matches per pattern: {dict(match_counts)}")
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

        if "mean_signal" in stats:
            report_lines.append(f"  Mean Signal: {stats.get('mean_signal', 0):.2f}")
        report_lines.append("")

        # Peak length distribution
        length_dist = stats.get("length_distribution", {})
        if length_dist:
            report_lines.append("Peak Length Distribution:")
            for length_range, count in sorted(length_dist.items()):
                report_lines.append(f"  {length_range} bp: {count:,} peaks")
            report_lines.append("")

        # Score distribution
        score_dist = stats.get("score_distribution", {})
        if score_dist:
            report_lines.append("Peak Score Distribution:")
            for score_range, count in sorted(score_dist.items()):
                report_lines.append(f"  {score_range}: {count:,} peaks")
            report_lines.append("")

        # Chromosome distribution (top 10)
        chr_dist = stats.get("chromosome_distribution", {})
        if chr_dist:
            report_lines.append("Chromosome Distribution (Top 10):")
            sorted_chrs = sorted(chr_dist.items(), key=lambda x: x[1], reverse=True)[:10]
            for chr_name, count in sorted_chrs:
                report_lines.append(f"  {chr_name}: {count:,} peaks")
            report_lines.append("")

    report = "\n".join(report_lines)

    if output_path:
        output_path = Path(output_path)
        with open(output_path, "w") as f:
            f.write(report)
        logger.info(f"ChIP-seq report saved to {output_path}")

    return report
