# Peak Calling

De novo peak detection from epigenomic signal data using Poisson enrichment models, broad-peak calling, peak merging, quality filtering, FRiP calculation, and differential peak analysis.

## Key Concepts

**Narrow peaks** are sharp, well-defined signal enrichments typically produced by transcription factor ChIP-seq or ATAC-seq. Detected using local Poisson enrichment against a background estimate.

**Broad peaks** span larger genomic regions characteristic of histone modifications like H3K36me3 or H3K27me3. Detected using a two-pass approach: seed regions passing a stringent threshold are extended through flanking regions passing a relaxed threshold.

**FRiP (Fraction of Reads in Peaks)** is a key quality metric measuring the proportion of total reads falling within called peaks. ENCODE guidelines recommend FRiP >= 1%.

**Differential peaks** identify regions where signal intensity changes significantly between conditions using fold-change and statistical testing.

## Function Reference

### `call_peaks_simple(signal, background=None, p_threshold=1e-5, min_length=100, merge_distance=50) -> List[Dict]`

Call narrow peaks using local Poisson enrichment testing.

**Parameters:**
- `signal`: List of signal values per genomic bin.
- `background`: Optional background/input signal for comparison. If None, uses a local background estimate.
- `p_threshold`: Poisson p-value threshold (default 1e-5).
- `min_length`: Minimum peak width in bins.
- `merge_distance`: Merge peaks within this distance (bins).

**Returns** list of dicts with `start`, `end`, `summit`, `score`, `p_value`, `fold_enrichment`.

### `call_peaks_broad(signal, background=None, seed_threshold=1e-5, extend_threshold=1e-3, min_length=500, max_gap=200) -> List[Dict]`

Call broad peaks using a two-pass seed-and-extend strategy. Seeds are identified at `seed_threshold`, then extended through flanking regions passing `extend_threshold`. Gaps up to `max_gap` bins are bridged.

### `merge_peaks(peaks, max_gap=0) -> List[Dict]`

Merge overlapping or nearby peaks. Scores and fold enrichment are combined (max or weighted average).

### `filter_peaks(peaks, min_fold=None, max_q_value=None, blacklist_regions=None) -> List[Dict]`

Filter peaks by fold enrichment, q-value, and/or blacklist regions. Blacklist regions are genomic intervals known to produce artifacts.

### `compute_frip(peaks, total_reads) -> Dict`

Calculate the Fraction of Reads in Peaks.

**Returns** dict with:
- `frip`: Fraction of reads in peaks (0.0--1.0).
- `reads_in_peaks`: Total reads overlapping peaks.
- `total_reads`: Total library size.
- `passes_encode`: Whether FRiP >= 0.01.

### `differential_peaks(peaks_a, signal_a, peaks_b, signal_b, min_fold_change=2.0) -> List[Dict]`

Identify peaks with significantly different signal intensity between two conditions. Uses log2 fold change and significance testing.

Returns dicts with `chrom`, `start`, `end`, `log2fc`, `p_value`, `direction`.

## Usage Examples

```python
from metainformant.epigenome import (
    call_peaks_simple,
    call_peaks_broad,
    filter_peaks as filter_called_peaks,
    compute_frip,
    differential_peaks,
)

# Narrow peak calling
signal = [0.5, 0.8, 1.2, 5.0, 8.0, 12.0, 7.0, 3.0, 1.0, 0.5]
background = [0.5] * 10
peaks = call_peaks_simple(signal, background, p_threshold=1e-5)

# Broad peak calling for histone marks
broad_signal = [0.5, 1.0, 2.0, 2.5, 3.0, 2.8, 2.5, 2.0, 1.5, 0.5]
broad_peaks = call_peaks_broad(broad_signal, seed_threshold=1e-5,
                                extend_threshold=1e-3)

# Quality filtering
filtered = filter_called_peaks(peaks, min_fold=2.0, max_q_value=0.05)

# FRiP quality control
frip = compute_frip(filtered, total_reads=10_000_000)
print(f"FRiP: {frip['frip']:.3f}, ENCODE pass: {frip['passes_encode']}")

# Differential peaks between conditions
diff = differential_peaks(treatment_peaks, treatment_signal,
                          control_peaks, control_signal,
                          min_fold_change=2.0)
```

## Configuration

Environment variable prefix: `EPI_`

## Related Modules

- `metainformant.epigenome.chipseq` -- ChIP-seq peak data
- `metainformant.epigenome.atacseq` -- ATAC-seq peak data
- `metainformant.epigenome.chromatin_state` -- state annotation from peaks
- `metainformant.epigenome.workflow` -- end-to-end analysis
