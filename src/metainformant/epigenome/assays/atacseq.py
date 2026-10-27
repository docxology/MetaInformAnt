"""ATAC-seq analysis and processing.

This module provides comprehensive tools for analyzing ATAC-seq data,
including peak calling, accessibility analysis, quality control, and
transcription factor binding site identification.
"""

from __future__ import annotations

import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

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
from metainformant.epigenome.assays.chipseq import (
    _IUPAC_BASES,
    _iter_matches,
    _iupac_pattern_to_regex,
    _read_fasta_sequences,
    _reverse_complement_iupac,
)

logger = logging.get_logger(__name__)


class ATACPeak:
    """Represents an ATAC-seq accessible region.

    ``p_value`` and ``q_value`` are probabilities (between 0 and 1, or ``None``). Peak files store them as
    ``-log10`` values, as the ENCODE specification and MACS2 do; :func:`load_atac_peaks` and
    :func:`save_atac_peaks` convert at the file boundary, so a q-value of 0.01 is the column value 2.
    """

    def __init__(
        self,
        chromosome: str,
        start: int,
        end: int,
        score: float = 0.0,
        strand: str = ".",
        signal_value: float = 0.0,
        p_value: Optional[float] = None,
        q_value: Optional[float] = None,
        summit: Optional[int] = None,
    ):
        """Initialize an ATAC-seq peak.

        Args:
            chromosome: Chromosome name
            start: Peak start position
            end: Peak end position
            score: Peak score
            strand: DNA strand
            signal_value: Signal value
            p_value: P-value, as a probability between 0 and 1
            q_value: Q-value (FDR), as a probability between 0 and 1
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

    def to_bed_format(self, format: str = "narrowpeak") -> str:
        """Convert to one line of a standard peak file.

        The columns follow the ENCODE layouts that :func:`load_atac_peaks` reads: ``narrowpeak`` (default,
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
            name=f"atac_peak_{self.start}_{self.end}",
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
            "chromosome": self.chromosome,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "score": self.score,
            "strand": self.strand,
            "signal_value": self.signal_value,
            "accessibility_score": self.accessibility_score,
            "p_value": self.p_value,
            "q_value": self.q_value,
            "summit": self.summit,
        }


def load_atac_peaks(path: str | Path, format: str = "narrowpeak") -> List[ATACPeak]:
    """Load ATAC-seq peaks from a file.

    The columns are read in the ENCODE layouts that :func:`save_atac_peaks` writes: ``narrowpeak`` has ten
    columns (``chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak``, with ``peak`` the
    summit as an offset from ``start``), ``broadpeak`` the first nine and ``bed`` the first six (only the first
    three are required). ``-1`` or ``.`` in the p-value, q-value or summit column means "not available" and
    gives ``None``.

    The pValue and qValue columns hold ``-log10`` values, as the ENCODE specification and MACS2 write them, and
    the peaks carry probabilities: a column value ``v`` becomes ``p_value`` or ``q_value`` ``10 ** -v``, so a
    column of 2 is a q-value of 0.01 and a MACS2 value of 9.1 is about 7.9e-10. A column above about 323.6
    underflows to a probability of 0.0. Any negative column value other than ``-1`` (or ``nan``) is not a
    ``-log10`` probability, and the line is skipped with a warning.

    A line with at least six columns but fewer than the requested layout needs is read as BED, and a BED line
    carries no signal column, so its score stands in as the signal value. Blank, ``#``, ``track`` and
    ``browser`` lines are ignored; a line with too few columns or a value that does not parse is skipped with a
    warning.

    Args:
        path: Path to peak file (optionally gzip-compressed)
        format: Peak file layout ("narrowpeak", "broadpeak", "bed"), in any letter case

    Returns:
        List of ATACPeak objects

    Raises:
        ValidationError: If the file doesn't exist
        ValueError: If the format is not supported
    """
    fmt = normalise_peak_format(format)
    path = validation.validate_path_exists(Path(path))

    logger.info(f"Loading ATAC-seq peaks from {path} (format: {fmt})")

    peaks = []
    # A short line is still read as BED, as long as it holds the six BED columns
    minimum_columns = MIN_COLUMNS["bed"] if fmt == "bed" else 6

    try:
        with io.open_text_auto(path, "rt") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if is_skipped_line(line):
                    continue

                parts = line.split("\t")
                if len(parts) < minimum_columns:
                    logger.warning(f"Skipping malformed line {line_num}: insufficient columns")
                    continue

                line_format = fmt if len(parts) >= MIN_COLUMNS[fmt] else "bed"

                try:
                    fields = parse_peak_fields(parts, line_format)
                    peak = ATACPeak(
                        chromosome=fields.chromosome,
                        start=fields.start,
                        end=fields.end,
                        score=fields.score,
                        strand=fields.strand,
                        # Basic (BED) lines carry no signal column: the score stands in as the signal value
                        signal_value=fields.signal_value if fields.signal_value is not None else fields.score,
                        p_value=fields.p_value,
                        q_value=fields.q_value,
                        summit=fields.summit,
                    )
                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue

                peaks.append(peak)

    except Exception as e:
        logger.error(f"Error loading ATAC-seq peaks from {path}: {e}")
        raise errors.IOError(f"Failed to load ATAC-seq peaks: {e}") from e

    logger.info(f"Loaded {len(peaks)} ATAC-seq peaks")
    return peaks


def save_atac_peaks(peaks: List[ATACPeak], path: str | Path, format: str = "narrowpeak") -> None:
    """Save ATAC-seq peaks to a file that :func:`load_atac_peaks` reads back.

    The columns follow the ENCODE layouts (see :meth:`ATACPeak.to_bed_format`): ten columns for ``narrowpeak``,
    nine for ``broadpeak`` and six for ``bed``. ``narrowpeak`` and ``broadpeak`` files start with a ``#``
    comment line, which the loaders skip. The name column (``atac_peak_<start>_<end>``) is regenerated on every
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
        peaks: List of ATACPeak objects
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

    logger.info(f"Saving {len(peaks)} ATAC-seq peaks to {path} (format: {fmt})")

    with open(path, "w") as f:
        if fmt != "bed":
            f.write("# Generated by MetaInformant ATAC-seq analysis\n")

        for line in lines:
            f.write(line + "\n")

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

    # Peak length distribution (ATAC-seq peaks are typically nucleosome-free)
    length_bins: Dict[str, int] = defaultdict(int)
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
    accessibility_bins: Dict[str, int] = defaultdict(int)
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
    chr_counts: Dict[str, int] = defaultdict(int)
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
    metrics: Dict[str, Any] = {}

    # Fragment size analysis (based on peak lengths)
    lengths = [p.length for p in peaks]

    # Typical ATAC-seq peaks should be around nucleosome-free regions (~100-200bp)
    nfr_peaks = sum(1 for length in lengths if 50 <= length <= 150)  # Nucleosome-free region
    mono_peaks = sum(1 for length in lengths if 150 <= length <= 250)  # Mononucleosome
    di_peaks = sum(1 for length in lengths if 250 <= length <= 350)  # Dinucleosome

    total_peaks = len(peaks)
    if total_peaks > 0:
        metrics.update(
            {
                "nfr_peak_fraction": nfr_peaks / total_peaks,
                "mononucleosome_peak_fraction": mono_peaks / total_peaks,
                "dinucleosome_peak_fraction": di_peaks / total_peaks,
            }
        )

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
                if chr_peak_list[i].start - chr_peak_list[i - 1].end <= 1000:
                    cluster_size += 1
                else:
                    if cluster_size >= 2:  # 3+ peaks in cluster
                        hotspot_count += 1
                    cluster_size = 0

            if cluster_size >= 2:
                hotspot_count += 1

        metrics["accessibility_hotspots"] = hotspot_count

    return metrics


def identify_tss_enrichment(
    peaks: List[ATACPeak], tss_positions: Dict[str, List[int]], window_size: int = 2000
) -> Dict[str, Any]:
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


def _match_probability(pattern: str, frequencies: Dict[str, float]) -> float:
    """Probability that a window of independent bases matches an IUPAC pattern.

    Args:
        pattern: Motif pattern using IUPAC nucleotide codes.
        frequencies: Probability of each of ``A``, ``C``, ``G`` and ``T``.

    Returns:
        Product over the pattern positions of the summed frequency of the bases that position allows.
    """
    probability = 1.0
    for char in pattern.upper():
        probability *= sum(frequencies[base] for base in _IUPAC_BASES[char])
    return probability


def find_tf_binding_sites(
    peaks: List[ATACPeak],
    tf_motifs: Dict[str, str],
    genome_fasta: Optional[str | Path] = None,
    peak_sequences: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """Find transcription factor binding sites in accessible regions.

    Scans the DNA of every peak for each motif on both strands and reports every occurrence, including ones that
    overlap each other. Each motif is a consensus string in IUPAC nucleotide codes (``"CCGCGNGGNGGCAG"``,
    ``"TGASTCA"``); a position matches when the base is one the code allows. The DNA must be supplied, either as a
    FASTA file covering the peaks' chromosomes (``genome_fasta``; a peak whose chromosome is missing from it is
    not scanned) or directly as one sequence per peak (``peak_sequences``, in the order of ``peaks``, each exactly
    as long as its peak and starting at ``peak.start``). Nothing is scanned or reported without sequence.

    Args:
        peaks: List of ATACPeak objects
        tf_motifs: Dictionary mapping TF names to IUPAC consensus motifs
        genome_fasta: Path to a FASTA file holding the sequence of the peaks' chromosomes
        peak_sequences: One DNA sequence per peak, in the order of ``peaks``

    Returns:
        Dictionary with the results:

        - ``peaks_analyzed``: number of input peaks; ``peaks_scanned``: how many of them had sequence
        - ``motifs_analyzed``: number of motifs
        - ``motif_counts``: per TF, the number of peaks with at least one site
        - ``site_counts``: per TF, the total number of sites (a palindromic motif counts once per strand)
        - ``expected_sites``: per TF, the number of sites expected by chance in the scanned sequence, from its
          base composition (independent bases); ``None`` when it cannot be computed
        - ``motif_enrichment``: per TF, ``site_counts / expected_sites``, or ``None`` when ``expected_sites`` is
          missing or zero
        - ``<TF>_positions``: per TF, a list of site dicts with ``chromosome``, ``position`` (0-based forward
          coordinate of the site start), ``end``, ``strand`` (``"+"`` or ``"-"``), ``sequence`` (the site read
          in the motif's orientation), ``peak_start``, ``peak_end`` and ``peak_score``

    Raises:
        ValueError: If neither ``genome_fasta`` nor ``peak_sequences`` is given, if both are, if a motif is empty
            or not IUPAC, or if ``peak_sequences`` does not match ``peaks`` in number or length.
    """
    if genome_fasta is not None and peak_sequences is not None:
        raise ValueError("Pass either genome_fasta or peak_sequences, not both")

    sequences: List[Optional[str]] = []
    if peak_sequences is not None:
        if len(peak_sequences) != len(peaks):
            raise ValueError(
                f"peak_sequences holds {len(peak_sequences)} sequences for {len(peaks)} peaks; "
                "give one per peak, in order"
            )
        for peak, given_sequence in zip(peaks, peak_sequences):
            if len(given_sequence) != peak.length:
                raise ValueError(
                    f"The sequence given for peak {peak.chromosome}:{peak.start}-{peak.end} is "
                    f"{len(given_sequence)} bp, but the peak is {peak.length} bp long"
                )
            sequences.append(given_sequence.upper())
    elif genome_fasta is not None:
        genome = _read_fasta_sequences(genome_fasta) if peaks else {}
        for peak in peaks:
            chromosome_sequence = genome.get(peak.chromosome)
            sequences.append(None if chromosome_sequence is None else chromosome_sequence[peak.start : peak.end])
    else:
        raise ValueError(
            "find_tf_binding_sites needs the DNA sequence of the peaks: pass genome_fasta (a FASTA file covering the "
            "peaks' chromosomes) or peak_sequences (one sequence per peak). There is nothing to scan without it."
        )

    logger.info(f"Finding TF binding sites for {len(tf_motifs)} motifs in {len(peaks)} peaks")

    forward = {tf_name: _iupac_pattern_to_regex(motif) for tf_name, motif in tf_motifs.items()}
    reverse = {
        tf_name: _iupac_pattern_to_regex(_reverse_complement_iupac(motif)) for tf_name, motif in tf_motifs.items()
    }

    composition = {base: 0 for base in "ACGT"}
    peaks_scanned = 0
    windows: Dict[str, int] = {tf_name: 0 for tf_name in tf_motifs}
    peaks_with_site: Dict[str, int] = {tf_name: 0 for tf_name in tf_motifs}
    sites: Dict[str, List[Dict[str, Any]]] = {tf_name: [] for tf_name in tf_motifs}

    for peak, sequence in zip(peaks, sequences):
        if sequence is None:
            continue  # Chromosome missing from the FASTA: nothing to scan, nothing reported
        peaks_scanned += 1
        for base in composition:
            composition[base] += sequence.count(base)

        for tf_name, motif in tf_motifs.items():
            windows[tf_name] += max(0, len(sequence) - len(motif) + 1)

            peak_sites: List[Dict[str, Any]] = []
            for strand, regex in (("+", forward[tf_name]), ("-", reverse[tf_name])):
                for offset, matched in _iter_matches(regex, sequence):
                    peak_sites.append(
                        {
                            "chromosome": peak.chromosome,
                            "position": peak.start + offset,
                            "end": peak.start + offset + len(motif),
                            "strand": strand,
                            # A minus-strand site is read in the motif's orientation, as in find_motifs_in_peaks
                            "sequence": matched if strand == "+" else reverse_complement(matched),
                            "peak_start": peak.start,
                            "peak_end": peak.end,
                            "peak_score": peak.score,
                        }
                    )
            if peak_sites:
                peaks_with_site[tf_name] += 1
                peak_sites.sort(key=lambda site: (site["position"], site["strand"]))
                sites[tf_name].extend(peak_sites)

    # Expected sites by chance: every motif-length window of every scanned sequence matches with the probability
    # given by the sequence's own base composition (independent bases), once per strand.
    total_bases = sum(composition.values())
    frequencies = {base: count / total_bases for base, count in composition.items()} if total_bases else {}

    expected_sites: Dict[str, Optional[float]] = {}
    motif_enrichment: Dict[str, Optional[float]] = {}
    for tf_name, motif in tf_motifs.items():
        expected: Optional[float] = None
        if frequencies and windows[tf_name] > 0:
            per_window = _match_probability(motif, frequencies) + _match_probability(
                _reverse_complement_iupac(motif), frequencies
            )
            expected = windows[tf_name] * per_window
        expected_sites[tf_name] = expected
        motif_enrichment[tf_name] = len(sites[tf_name]) / expected if expected else None

    results: Dict[str, Any] = {
        "peaks_analyzed": len(peaks),
        "peaks_scanned": peaks_scanned,
        "motifs_analyzed": len(tf_motifs),
        "motif_counts": peaks_with_site,
        "site_counts": {tf_name: len(tf_sites) for tf_name, tf_sites in sites.items()},
        "expected_sites": expected_sites,
        "motif_enrichment": motif_enrichment,
    }
    for tf_name, tf_sites in sites.items():
        results[f"{tf_name}_positions"] = tf_sites

    return results


def calculate_chromatin_accessibility_index(
    peaks: List[ATACPeak], genomic_regions: List[Tuple[str, int, int]]
) -> Dict[str, Any]:
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


def compare_atac_conditions(condition1_peaks: List[ATACPeak], condition2_peaks: List[ATACPeak]) -> Dict[str, Any]:
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
            has_overlap = any(peak2.overlaps_with(p1) for p1 in condition1_peaks if p1.chromosome == peak2.chromosome)
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
        if "nfr_peak_fraction" in stats:
            report_lines.append("Nucleosome Positioning:")
            report_lines.append(f"  Nucleosome-free regions: {stats.get('nfr_peak_fraction', 0):.1%}")
            report_lines.append(f"  Mononucleosome peaks: {stats.get('mononucleosome_peak_fraction', 0):.1%}")
            report_lines.append(f"  Dinucleosome peaks: {stats.get('dinucleosome_peak_fraction', 0):.1%}")
            report_lines.append("")

        # Accessibility hotspots
        if "accessibility_hotspots" in stats:
            report_lines.append(f"Accessibility Hotspots: {stats.get('accessibility_hotspots', 0)}")
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
        logger.info(f"ATAC-seq report saved to {output_path}")

    return report
