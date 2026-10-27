# ChIP-seq Analysis

Chromatin immunoprecipitation sequencing analysis for histone modifications and transcription factor binding, implemented in `metainformant.epigenome.assays.chipseq`. Includes peak loading and saving, score filtering, statistics, overlap analysis, merging, genome-level enrichment, and motif scanning.

## Key Concepts

**ChIP-seq peaks** represent genomic regions enriched for a specific protein-DNA interaction (histone mark or transcription factor binding). Each peak has a location, a score, a signal value, and optional p-value and q-value (FDR) estimates. The p-value and q-value are probabilities between 0 and 1; peak files store them as `-log10` values (see [Peak File Layouts](#peak-file-layouts)).

**Peak overlaps** identify shared regulatory elements between experiments or conditions by finding genomic intervals that intersect.

**Peak enrichment** summarizes how much of the genome the peaks cover and, when an expected peak count is supplied, the observed/expected ratio and the Poisson probability of seeing at least that many peaks. It does not test enrichment in annotated genomic features such as promoters.

## Data Model

### `ChIPPeak`

Class representing a single ChIP-seq peak (a plain class, not a dataclass). The constructor requires a non-empty `chromosome`, `start >= 0`, `end >= start`, and, when given, `start <= summit <= end`. The p-value and q-value are not range-checked when a peak is built; saving a narrowPeak or broadPeak file rejects a value outside 0--1.

| Field | Type | Description |
|-------|------|-------------|
| `chromosome` | `str` | Chromosome name |
| `start` | `int` | Start coordinate |
| `end` | `int` | End coordinate |
| `summit` | `int \| None` | Absolute position of maximum signal (default `None`) |
| `score` | `float` | Peak score (default `0.0`) |
| `strand` | `str` | Strand (default `"."`) |
| `signal_value` | `float` | Signal value (default `0.0`) |
| `p_value` | `float \| None` | P-value as a probability (default `None`) |
| `q_value` | `float \| None` | FDR-adjusted significance as a probability (default `None`) |

Read-only properties: `length` (`end - start`) and `center`. Methods: `to_bed_format(format="narrowpeak")` (one line of a peak file in the given layout, named `peak_<start>_<end>`), `overlaps_with(other, min_overlap=1)`, and `to_dict()`. There is no peak name or fold-enrichment field.

## Peak File Layouts

`load_chip_peaks` and `save_chip_peaks`, and their ATAC-seq counterparts `load_atac_peaks` and `save_atac_peaks`, read and write the three ENCODE/UCSC tab-separated layouts. The columns are defined once, in `assays/_peak_format.py`, so a file written with one layout is read back with the same layout.

| Layout | Columns |
|--------|---------|
| `narrowpeak` | 10: `chrom`, `chromStart`, `chromEnd`, `name`, `score`, `strand`, `signalValue`, `pValue`, `qValue`, `peak` |
| `broadpeak` | The first 9 narrowPeak columns (no `peak`) |
| `bed` | The first 6 columns: `chrom`, `chromStart`, `chromEnd`, `name`, `score`, `strand`; reading needs only the first 3 |

- **Summit.** The narrowPeak `peak` column is the summit as an offset from `chromStart`; the `summit` attribute is the absolute coordinate (`start + offset`). broadPeak and BED files hold no summit, and a BED file holds neither the signal value nor the p- and q-values.
- **Missing values.** A missing p-value, q-value or summit is written as `-1`, and `-1` or `.` in those columns is read back as `None`.
- **p- and q-values.** Peaks carry probabilities, while the files store `-log10` values, as the ENCODE specification and MACS2 do. Reading turns a column value `v` into `10 ** -v` (a column of 2 is 0.01, and a MACS2 value of 9.1 is about 7.9e-10); writing turns a probability `p` into `-log10(p)`. A probability of exactly 0 is written as about 323.306, the `-log10` of the smallest positive float, and reads back as that float (5e-324) rather than 0.0; a column above about 323.6 underflows to 0.0. A round trip loses a few units in the last place, and values below about 1e-308 lose precision as the float underflows.
- **Invalid values.** A negative column other than `-1` (or `nan`) is not a `-log10` probability, so the loader skips the line with a warning. A p-value or q-value outside 0--1 raises `ValueError` when saving (the narrowPeak and broadPeak layouts, which write them).
- **Format names.** `format` is case-insensitive (`"narrowPeak"` works), and an unsupported name raises `ValueError`.
- **Saving.** Every line is formatted before the file is opened, so a failure leaves no partial file and does not overwrite an existing one. narrowPeak and broadPeak files start with a `#` comment line, which loaders skip. The name column is regenerated on every write, because peaks carry no name. Loaders read gzip files by their `.gz` suffix; saving always writes plain text.

## Function Reference

### `load_chip_peaks(path, format="narrowpeak") -> List[ChIPPeak]`

Load peaks from a tab-delimited peak file in one of the layouts above. `format` is `"narrowpeak"`, `"broadpeak"` or `"bed"` in any letter case; any other value raises `ValueError`, and a missing file raises `metainformant.core.utils.errors.ValidationError`. narrowPeak lines need all ten columns, broadPeak lines nine and BED lines three. Blank, `#`, `track` and `browser` lines are ignored; a line with too few columns, or with a value that does not parse (including a negative p or q column other than `-1`), is skipped with a warning. A line that parses but breaks the `ChIPPeak` constraints (for example `end < start`, or a summit outside `start`--`end`) is not skipped: it aborts the whole load with `metainformant.core.utils.errors.IOError`.

### `save_chip_peaks(peaks, path, format="narrowpeak") -> None`

Write peaks to a text file that `load_chip_peaks` reads back, in the layout given by `format` (case-insensitive; an unsupported value raises `ValueError`). The summit is written as an offset from `start`, so a narrowPeak file restores it. See [Peak File Layouts](#peak-file-layouts) for the columns, the `-log10` conversion, the `-1` placeholders and the all-or-nothing write.

### `filter_peaks_by_score(peaks, min_score, max_peaks=None) -> List[ChIPPeak]`

Keep peaks with `score >= min_score`, sorted by score from highest to lowest and optionally truncated to the `max_peaks` best. Filtering on q-value, signal or chromosome is done with a list comprehension over the peak attributes.

### `calculate_peak_statistics(peaks) -> Dict`

Summary statistics, or `{}` for an empty list: `total_peaks`; `mean_length`, `median_length`, `min_length`, `max_length`, `length_std`; `mean_score`, `median_score`, `min_score`, `max_score`, `score_std`; `mean_signal`, `median_signal`, `signal_std` (computed over the peaks with a positive `signal_value`, and present only when there is one); binned `length_distribution` and `score_distribution`; and `chromosome_distribution`.

### `find_overlapping_peaks(peaks1, peaks2, min_overlap=1) -> List[Tuple[ChIPPeak, ChIPPeak]]`

Find every pair of peaks, one from each list, that are on the same chromosome and share at least `min_overlap` bp. Returns `(peak1, peak2)` tuples; the overlap length is not included, so compute it from the coordinates. Peaks within one list may overlap each other or be nested, and a peak that overlaps several peaks of the other list appears in several pairs. The result is ordered by chromosome (in order of first appearance in `peaks1`), then by `peak1` start, then by `peak2` start, so it does not depend on hash seeds.

### `merge_overlapping_peaks(peaks, max_distance=0) -> List[ChIPPeak]`

Merge peaks that overlap, touch, or lie within `max_distance` bp of each other into single peaks, per chromosome. The result lists chromosomes in order of first appearance, each sorted by start. A merged peak keeps the summit and strand of the first (lowest-`start`) peak, the highest score and signal value, and the lowest p-value and q-value among the peaks it joins (`None` as soon as one of them has none). A p-value or q-value of 0.0 is a value, not a missing one.

### `calculate_peak_enrichment(peaks, genome_size, expected_peaks=None) -> Dict`

Genome-level enrichment summary, or `{}` for an empty list. Always returns `total_peaks`, `total_peak_bases`, `genome_size`, `genome_coverage` (summed peak length divided by `genome_size`) and `average_peak_length`. When `expected_peaks` is a positive number it also returns `expected_peaks`, `enrichment_ratio` (observed over expected peaks, repeated as `fold_enrichment`), `enrichment_p_value` (Poisson probability of at least the observed number of peaks) and `enrichment_log10_p_value` (its base-10 logarithm). The tail is summed in log space, so any number of peaks is handled without overflow. The p-value loses precision below about 1e-308 and is 0.0 below about 5e-324, while `enrichment_log10_p_value` stays finite.

### `find_motifs_in_peaks(peaks, genome_fasta, motif_patterns, window_size=200) -> Dict`

Scan peaks for IUPAC motif patterns (for example `"TTGACAW"`, in any letter case) in a genome FASTA file. For each peak, the `window_size` bp window centered on the summit (the peak center when there is no summit, clipped to the peak) is searched on both strands. Every occurrence is reported, including one that overlaps an earlier occurrence of the same pattern (`ATAT` is found at offsets 0 and 2 of `ATATAT`). Peaks on chromosomes absent from the FASTA contribute no matches. An empty pattern, or one with a character outside the IUPAC nucleotide codes, raises `ValueError` before the FASTA is read.

Returns `total_peaks_analyzed`; `motif_counts` (peaks with at least one match, per pattern); `match_counts` (total matches per pattern); and `motif_positions` (per pattern, a list of dicts with `chromosome`, `position` (0-based start on the forward strand), `strand`, `sequence` (the match in the motif's orientation), `peak_start`, `peak_end` and `peak_score`). A pattern with no match at all is absent from all three dicts. There is no de novo motif discovery.

## Usage Examples

```python-snippet
from metainformant.epigenome.assays.chipseq import (
    calculate_peak_enrichment,
    calculate_peak_statistics,
    filter_peaks_by_score,
    find_overlapping_peaks,
    load_chip_peaks,
    merge_overlapping_peaks,
    save_chip_peaks,
)

# Load peaks (q-values arrive as probabilities), keep the highest-scoring ones, then filter on q-value
peaks = load_chip_peaks("h3k27ac_peaks.narrowPeak")
top_peaks = filter_peaks_by_score(peaks, min_score=100, max_peaks=1000)
high_quality = [p for p in top_peaks if p.q_value is not None and p.q_value <= 0.01]

# Summary statistics
stats = calculate_peak_statistics(high_quality)
print(f"Total peaks: {stats['total_peaks']}, median length: {stats['median_length']}")

# Find overlaps between two marks (peaks may overlap within a set)
h3k4me3 = load_chip_peaks("h3k4me3_peaks.narrowPeak")
overlaps = find_overlapping_peaks(high_quality, h3k4me3)
print(f"Overlapping peak pairs: {len(overlaps)}")

# Merge peaks within 200 bp of each other and save them as narrowPeak
merged = merge_overlapping_peaks(high_quality, max_distance=200)
save_chip_peaks(merged, "h3k27ac_merged.narrowPeak")

# Genome-level enrichment against an expected number of peaks
enrichment = calculate_peak_enrichment(
    high_quality, genome_size=3_000_000_000, expected_peaks=500
)
print(f"Fold enrichment: {enrichment['fold_enrichment']:.2f}")
print(f"log10 p-value: {enrichment['enrichment_log10_p_value']:.1f}")
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
