# Peak Calling

De novo peak detection from epigenomic signal data, implemented in `metainformant.epigenome.peak_calling.peak_detection`: Poisson enrichment models, broad-peak calling, peak merging, quality filtering, FRiP calculation, and differential peak analysis. The caller works on plain Python lists with one signal value per genomic bin and returns peaks as dicts whose `start` and `end` are bin indices.

## Key Concepts

**Narrow peaks** are sharp, well-defined signal enrichments typically produced by transcription factor ChIP-seq or ATAC-seq. Detected using local Poisson enrichment against a background estimate.

**Broad peaks** span larger genomic regions characteristic of histone modifications like H3K36me3 or H3K27me3. Detected using a two-pass approach: seed regions passing a stringent threshold are extended through flanking regions passing a relaxed threshold.

**FRiP (Fraction of Reads in Peaks)** is a key quality metric measuring the proportion of total reads falling within called peaks. ENCODE guidelines suggest FRiP >= 1% for a usable experiment and >= 5% for a good one.

**Differential peaks** identify regions where signal intensity changes between conditions using the log2 fold change of mean signal and a statistical test.

## Function Reference

### `call_peaks_simple(signal, control=None, threshold=5.0, min_length=100, merge_distance=50, chrom="chr_unknown") -> List[Dict]`

Call narrow peaks using local Poisson enrichment testing. Each bin is tested against the larger of the global mean and the local mean (the maximum over windows of +/-500, +/-2500 and +/-5000 bins) of the control track, or of `signal` itself when no control is given.

**Parameters:**
- `signal`: List of signal values, one per genomic bin (must not be empty).
- `control`: Optional control/input signal of the same length (`ValueError` if the lengths differ). If `None`, the background is estimated from `signal`.
- `threshold`: Minimum -log10(p-value) for a bin to count as enriched (default 5.0, i.e. p <= 1e-5).
- `min_length`: Minimum peak width in bins (default 100, so short signals return no peaks).
- `merge_distance`: Merge enriched runs separated by at most this many bins (default 50).
- `chrom`: Chromosome the signal belongs to; it is written to the `chrom` key of every peak (`ValueError` if empty). The default `"chr_unknown"` is only a placeholder for a signal with no chromosome. Peak coordinates stay bin indices into `signal`; they equal genomic positions when each value covers one base and the track starts at position 0 of `chrom`. Set `chrom` when peaks from several chromosomes are merged, filtered against a blacklist or compared with `differential_peaks`, which all match peaks by their `chrom` key.

**Returns** a list of dicts with `chrom` (the `chrom` argument), `start`, `end` (bin indices, `end` exclusive), `summit`, `score` (-log10 of the best p-value in the peak, 300 when that p-value is 0), `p_value`, `fold_enrichment`, `signal_value` (maximum signal in the peak) and `q_value` (Benjamini-Hochberg across the returned peaks).

### `call_peaks_broad(signal, control=None, p_threshold=0.1, broad_cutoff=0.1, min_length=200, chrom="chr_unknown") -> List[Dict]`

Call broad peaks using a two-pass seed-and-extend strategy. Seeds are runs of bins with Poisson p <= `p_threshold / 10`; each seed is extended left and right while the p-value stays <= `broad_cutoff`. Regions shorter than `min_length` bins are dropped and overlapping domains are merged. `chrom` works as in `call_peaks_simple`.

Returns the same keys as `call_peaks_simple`, plus `peak_type` (`"broad"`).

### `merge_peaks(peaks, distance=0) -> List[Dict]`

Merge peaks whose gap is at most `distance` bins, one chromosome at a time. Peaks are grouped by their `chrom` key (`"chr_unknown"` for a peak without one) and sorted by `start` within each chromosome; peaks on different chromosomes are never merged. The result lists the chromosomes in order of first appearance. A merged peak keeps the summit and `signal_value` of the strongest constituent, the lowest `p_value` and `q_value` of its members (the q-value is not re-estimated for the merged set), and the highest `score` and `fold_enrichment`; every other key, including `chrom`, comes from the first peak.

### `filter_peaks(peaks, min_fold=2.0, max_qvalue=0.05, blacklist_regions=None) -> List[Dict]`

Keep peaks with `fold_enrichment >= min_fold` and `q_value <= max_qvalue` (a peak without `q_value` is judged on its `p_value`), then remove peaks that overlap a blacklist region. Peaks without fold enrichment or without either significance field are removed with a counted warning. `blacklist_regions` is a list of dicts with `chrom`, `start` and `end`; they are matched against each peak's `chrom` key.

### `compute_frip(reads_in_peaks, total_reads) -> Dict`

Calculate the Fraction of Reads in Peaks from two counts. The function does not count reads itself: `reads_in_peaks` is the number of reads overlapping the called peaks, obtained elsewhere. Raises `ValueError` if `total_reads` is not positive or `reads_in_peaks` is negative or larger than `total_reads`.

**Returns** dict with:
- `frip`: Fraction of reads in peaks (0.0--1.0).
- `reads_in_peaks`: Reads overlapping peaks.
- `total_reads`: Total library size.
- `reads_outside_peaks`: `total_reads - reads_in_peaks`.
- `percent_in_peaks`: `frip` as a percentage.
- `quality_assessment`: `"good"` for FRiP >= 0.05, `"acceptable"` for FRiP >= 0.01, otherwise `"low"`.

### `differential_peaks(peaks_a, peaks_b, signals_a, signals_b) -> List[Dict]`

Compare measured signal over the union of each peak pair. Each A peak selects the B peak with greatest overlap, and each B peak selects the A peak with greatest overlap; ties use start, end and then list order. A pair selected from both sides is emitted once. A peak overlapping anything is `shared`; `unique_A` and `unique_B` mean no partner exists.

Tracks may be lists for a single chromosome or chromosome-to-track dictionaries for several chromosomes. Unequal track lengths, invalid coordinates, empty/out-of-track regions, negative regional mean signals and ambiguous multi-chromosome plain tracks raise `ValueError`.

Every record, including a unique peak, uses both measured tracks. `log2_fold_change` is `log2((mean_a + 1e-10) / (mean_b + 1e-10))`. Its two-sided Welch p-value (normal approximation without SciPy) ranks regions; autocorrelated bases are not independent biological replicates. A one-base region has `NaN` p; equal constant signals give p=1 and unequal constant signals give p=0. `direction` is `up_in_A` or `up_in_B` only when absolute log2 fold change is at least 1 and p < 0.05, otherwise `unchanged`.

Returns `chrom`, `start`, `end`, `log2_fold_change`, `p_value`, `direction`, `status`, `mean_signal_a` and `mean_signal_b`. There is no fold-change threshold argument.

## Usage Examples

```python
from metainformant.epigenome.peak_calling.peak_detection import (
    call_peaks_broad,
    call_peaks_simple,
    compute_frip,
    differential_peaks,
    filter_peaks as filter_called_peaks,
)

# Narrow peak calling: one signal value per bin, with a matching control track
signal = [1.0] * 400 + [30.0] * 60 + [1.0] * 400
control = [1.0] * len(signal)
peaks = call_peaks_simple(signal, control, threshold=5.0, min_length=20, merge_distance=5, chrom="chr1")

# Broad peak calling for histone marks
broad_signal = [1.0] * 600 + [12.0] * 200 + [1.0] * 600
broad_peaks = call_peaks_broad(
    broad_signal, p_threshold=0.05, broad_cutoff=0.1, min_length=100, chrom="chr1"
)

# Quality filtering
filtered = filter_called_peaks(peaks, min_fold=2.0, max_qvalue=0.05)

# FRiP quality control (reads_in_peaks is counted elsewhere, e.g. from the alignment file)
frip = compute_frip(reads_in_peaks=150_000, total_reads=10_000_000)
print(f"FRiP: {frip['frip']:.3f}, quality: {frip['quality_assessment']}")

# Differential peaks between two conditions (the signal tracks must have equal length)
treated_signal = signal
untreated_signal = [1.0] * 400 + [10.0] * 60 + [1.0] * 400
untreated_peaks = call_peaks_simple(
    untreated_signal, control, threshold=5.0, min_length=20, merge_distance=5, chrom="chr1"
)
diff = differential_peaks(peaks, untreated_peaks, treated_signal, untreated_signal)
for region in diff:
    print(region["status"], region["direction"], region["log2_fold_change"])
```

## Configuration

No dedicated configuration file or environment prefix is used by the peak calling module; call the functions above directly from Python.

## Related Modules

- `metainformant.epigenome.assays.chipseq` -- ChIP-seq peak data
- `metainformant.epigenome.assays.atacseq` -- ATAC-seq peak data
- `metainformant.epigenome.chromatin_state` -- chromatin state learning from histone mark signal
- `metainformant.epigenome.workflow` -- end-to-end analysis
