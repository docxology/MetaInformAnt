# ChIP-seq Analysis

Chromatin immunoprecipitation sequencing analysis for histone modifications and transcription factor binding, implemented in `metainformant.epigenome.assays.chipseq`. Includes peak loading and saving, score filtering, statistics, overlap analysis, merging, genome-level enrichment, and motif scanning.

## Key Concepts

**ChIP-seq peaks** represent genomic regions enriched for a specific protein-DNA interaction (histone mark or transcription factor binding). Each peak has a location, a score, a signal value, and optional p-value and q-value (FDR) estimates.

**Peak overlaps** identify shared regulatory elements between experiments or conditions by finding genomic intervals that intersect.

**Peak enrichment** summarizes how much of the genome the peaks cover and, when an expected peak count is supplied, the observed/expected ratio and the Poisson probability of seeing at least that many peaks. It does not test enrichment in annotated genomic features such as promoters.

## Data Model

### `ChIPPeak`

Class representing a single ChIP-seq peak (a plain class, not a dataclass). The constructor requires a non-empty `chromosome`, `start >= 0`, `end >= start`, and, when given, `start <= summit <= end`.

| Field | Type | Description |
|-------|------|-------------|
| `chromosome` | `str` | Chromosome name |
| `start` | `int` | Start coordinate |
| `end` | `int` | End coordinate |
| `summit` | `int \| None` | Position of maximum signal (default `None`) |
| `score` | `float` | Peak score (default `0.0`) |
| `strand` | `str` | Strand (default `"."`) |
| `signal_value` | `float` | Signal value (default `0.0`) |
| `p_value` | `float \| None` | P-value (default `None`) |
| `q_value` | `float \| None` | FDR-adjusted significance (default `None`) |

Read-only properties: `length` (`end - start`) and `center`. Methods: `to_bed_format()`, `overlaps_with(other, min_overlap=1)`, and `to_dict()`. There is no peak name or fold-enrichment field.

## Function Reference

### `load_chip_peaks(path, format="narrowpeak") -> List[ChIPPeak]`

Load peaks from a tab-delimited peak file. `format` is `"narrowpeak"`, `"broadpeak"` or `"bed"`, in lower case: any other value (including `"narrowPeak"`) matches no parser and returns an empty list with warnings. narrowPeak lines need all ten columns (the tenth is the summit offset from `start`) and broadPeak lines nine; `"bed"` reads the first six columns. Comment (`#`) and `track` lines are ignored; lines with too few columns or unparseable values are skipped with a warning. A line that parses but breaks the `ChIPPeak` constraints (for example `end < start`, or a summit outside `start`--`end`) is not skipped: it aborts the whole load with `metainformant.core.utils.errors.IOError`.

### `save_chip_peaks(peaks, path, format="narrowpeak") -> None`

Write peaks to a BED-like text file: ten columns for `"narrowpeak"`, nine for `"broadpeak"`, six for `"bed"`; any other `format` value writes only a `#` header line and no peaks. After the strand column, the narrowpeak and broadpeak lines hold the summit coordinate, signal value and p-value (narrowpeak then adds the q-value, which broadpeak leaves out). `load_chip_peaks` reads signal value, p-value, q-value and then the summit offset, so a saved narrowpeak/broadpeak file cannot be reloaded faithfully.

### `filter_peaks_by_score(peaks, min_score, max_peaks=None) -> List[ChIPPeak]`

Keep peaks with `score >= min_score`, sorted by score from highest to lowest and optionally truncated to the `max_peaks` best. Filtering on q-value, signal or chromosome is done with a list comprehension over the peak attributes.

### `calculate_peak_statistics(peaks) -> Dict`

Summary statistics, or `{}` for an empty list: `total_peaks`; `mean_length`, `median_length`, `min_length`, `max_length`, `length_std`; `mean_score`, `median_score`, `min_score`, `max_score`, `score_std`; `mean_signal`, `median_signal`, `signal_std` (computed over the peaks with a positive `signal_value`, and present only when there is one); binned `length_distribution` and `score_distribution`; and `chromosome_distribution`.

### `find_overlapping_peaks(peaks1, peaks2, min_overlap=1) -> List[Tuple[ChIPPeak, ChIPPeak]]`

Find peak pairs on the same chromosome that share at least `min_overlap` bp. Returns `(peak1, peak2)` tuples; the overlap length is not included, so compute it from the coordinates. The sweep assumes peaks within each input list do not overlap one another and can miss pairs otherwise, so merge each list first with `merge_overlapping_peaks` when in doubt.

### `merge_overlapping_peaks(peaks, max_distance=0) -> List[ChIPPeak]`

Merge peaks that overlap, or lie within `max_distance` bp of each other, into single peaks per chromosome. A merged peak keeps the summit and strand of the first (lowest-`start`) peak, the highest score and signal value, and the lowest p-value and q-value when both peaks have a non-zero one (otherwise `None`).

### `calculate_peak_enrichment(peaks, genome_size, expected_peaks=None) -> Dict`

Genome-level enrichment summary, or `{}` for an empty list. Always returns `total_peaks`, `total_peak_bases`, `genome_size`, `genome_coverage` (summed peak length divided by `genome_size`) and `average_peak_length`. When `expected_peaks` is a positive number it also returns `expected_peaks`, `enrichment_ratio` (observed over expected peaks, repeated as `fold_enrichment`) and `enrichment_p_value` (Poisson probability of at least the observed number of peaks).

### `find_motifs_in_peaks(peaks, genome_fasta, motif_patterns, window_size=200) -> Dict`

Scan peaks for IUPAC motif patterns (for example `"TTGACAW"`) in a genome FASTA file. For each peak, the `window_size` bp window centered on the summit (the peak center when there is no summit, clipped to the peak) is searched on both strands. Peaks on chromosomes absent from the FASTA contribute no matches. Returns `total_peaks_analyzed`; `motif_counts` (peaks with at least one match, per pattern); `match_counts` (total matches per pattern); and `motif_positions` (per pattern, a list of dicts with `chromosome`, `position`, `strand`, `sequence`, `peak_start`, `peak_end` and `peak_score`). A pattern with no match at all is absent from all three dicts. There is no de novo motif discovery.

## Usage Examples

```python-snippet
from metainformant.epigenome.assays.chipseq import (
    calculate_peak_enrichment,
    calculate_peak_statistics,
    filter_peaks_by_score,
    find_overlapping_peaks,
    load_chip_peaks,
    merge_overlapping_peaks,
)

# Load peaks, keep the highest-scoring ones, then filter on q-value
peaks = load_chip_peaks("h3k27ac_peaks.narrowPeak")
top_peaks = filter_peaks_by_score(peaks, min_score=100, max_peaks=1000)
high_quality = [p for p in top_peaks if p.q_value is not None and p.q_value <= 0.01]

# Summary statistics
stats = calculate_peak_statistics(high_quality)
print(f"Total peaks: {stats['total_peaks']}, median length: {stats['median_length']}")

# Find overlaps between two marks (merge each set first so no peaks overlap within a set)
h3k4me3 = load_chip_peaks("h3k4me3_peaks.narrowPeak")
overlaps = find_overlapping_peaks(
    merge_overlapping_peaks(high_quality),
    merge_overlapping_peaks(h3k4me3),
)
print(f"Overlapping peak pairs: {len(overlaps)}")

# Genome-level enrichment against an expected number of peaks
enrichment = calculate_peak_enrichment(
    high_quality, genome_size=3_000_000_000, expected_peaks=500
)
print(f"Fold enrichment: {enrichment['fold_enrichment']:.2f}")
```

## Configuration

No dedicated configuration file or environment prefix is used by the ChIP-seq module; call the functions above directly from Python.

## Optional Dependencies

None. The module uses only the standard library and `metainformant.core`.

## Related Modules

- `metainformant.epigenome.assays.methylation` -- DNA methylation
- `metainformant.epigenome.assays.atacseq` -- chromatin accessibility
- `metainformant.epigenome.peak_calling` -- de novo peak calling from signal data
- `metainformant.epigenome.chromatin_state` -- chromatin state learning
