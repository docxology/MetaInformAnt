# Peak Calling

Signal-based peak detection for ChIP-seq and ATAC-seq data, including narrow peak calling with Poisson statistics, broad domain detection, summit refinement, and differential peak analysis.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `peak_detection` module |
| `peak_detection.py` | Peak calling algorithms, merging, filtering, FRiP metrics, differential analysis |

## Key Functions

| Function | Description |
|----------|-------------|
| `call_peaks_simple()` | Narrow peak calling using Poisson significance testing |
| `call_peaks_broad()` | Broad domain detection for diffuse histone marks |
| `peak_summit_refinement()` | Refine peak summits to single-base resolution |
| `merge_peaks()` | Merge overlapping or nearby peaks |
| `filter_peaks()` | Filter peaks by significance, fold-change, or size |
| `compute_frip()` | Calculate Fraction of Reads in Peaks quality metric |
| `compute_local_lambda()` | Estimate local background signal level |
| `differential_peaks()` | Identify differential peaks between conditions |

## Usage

```python
from metainformant.epigenome.peak_calling import peak_detection

peaks = peak_detection.call_peaks_simple(signal, control, p_threshold=1e-5)
broad = peak_detection.call_peaks_broad(signal, control)
merged = peak_detection.merge_peaks(peaks, max_gap=500)
frip = peak_detection.compute_frip(reads_in_peaks=50000, total_reads=1000000)
diff = peak_detection.differential_peaks(peaks_a, peaks_b)
```
