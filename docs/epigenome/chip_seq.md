# ChIP-seq Analysis

Chromatin immunoprecipitation sequencing analysis for histone modifications and transcription factor binding. Includes peak loading, filtering, statistics, overlap analysis, merging, enrichment, and motif finding.

## Key Concepts

**ChIP-seq peaks** represent genomic regions enriched for a specific protein-DNA interaction (histone mark or transcription factor binding). Each peak has a location, signal intensity (score), fold enrichment, and statistical significance (q-value).

**Peak overlaps** identify shared regulatory elements between experiments or conditions by finding genomic intervals that intersect.

**Enrichment analysis** assesses whether peaks are overrepresented in specific genomic features (promoters, enhancers, gene bodies) relative to a background model.

## Data Model

### `ChIPPeak`

Dataclass representing a single ChIP-seq peak.

| Field | Type | Description |
|-------|------|-------------|
| `chrom` | `str` | Chromosome |
| `start` | `int` | Start coordinate |
| `end` | `int` | End coordinate |
| `name` | `str` | Peak identifier |
| `score` | `float` | Signal intensity |
| `fold_enrichment` | `float` | Fold enrichment over background |
| `q_value` | `float` | FDR-adjusted significance |
| `summit` | `int` | Position of maximum signal |

## Function Reference

### `load_peaks(file_path, format="narrowPeak") -> List[ChIPPeak]`

Load peaks from ENCODE narrowPeak or broadPeak format files.

### `save_peaks(peaks, file_path, format="narrowPeak") -> None`

Write peaks to a BED-like file.

### `filter_peaks(peaks, min_score=None, max_q_value=None, min_fold=None, chromosomes=None) -> List[ChIPPeak]`

Filter peaks by quality thresholds: minimum score, maximum q-value, minimum fold enrichment, or specific chromosomes.

### `peak_statistics(peaks) -> Dict`

Summary statistics: total peaks, mean/median width, total coverage, mean fold enrichment, mean q-value, chromosome distribution.

### `find_overlapping_peaks(peaks_a, peaks_b, min_overlap=1) -> List[Tuple]`

Find peaks that overlap between two peak sets. Returns tuples of (peak_a, peak_b, overlap_bp).

### `merge_peaks(peaks, max_gap=0) -> List[ChIPPeak]`

Merge nearby or overlapping peaks within `max_gap` bp into single peaks, combining scores and enrichment values.

### `calculate_enrichment(peaks, features, genome_size) -> Dict`

Calculate fold enrichment of peaks in genomic features (e.g., promoters, exons) compared to random expectation. Uses Fisher's exact test when scipy is available.

### `find_motifs(peaks, sequences, motif_length=8, top_n=10) -> List[Dict]`

Discover enriched sequence motifs within peak regions using k-mer frequency analysis. Returns motifs with their enrichment scores and p-values.

## Usage Examples

```python
from metainformant.epigenome import (
    ChIPPeak,
    load_peaks,
    filter_peaks,
    peak_statistics,
    find_overlapping_peaks,
    calculate_enrichment,
)

# Load and filter peaks
peaks = load_peaks("h3k27ac_peaks.narrowPeak")
high_quality = filter_peaks(peaks, max_q_value=0.01, min_fold=2.0)

# Summary statistics
stats = peak_statistics(high_quality)
print(f"Total peaks: {stats['total_peaks']}, median width: {stats['median_width']}")

# Find overlaps between two marks
h3k4me3 = load_peaks("h3k4me3_peaks.narrowPeak")
overlaps = find_overlapping_peaks(high_quality, h3k4me3)
print(f"Overlapping peaks: {len(overlaps)}")

# Enrichment in promoters
features = {"promoter": [("chr1", 1000, 3000), ("chr1", 50000, 52000)]}
enrichment = calculate_enrichment(high_quality, features, genome_size=3_000_000_000)
```

## Configuration

Environment variable prefix: `EPI_`

## Optional Dependencies

- `scipy` -- Fisher's exact test for enrichment analysis

## Related Modules

- `metainformant.epigenome.methylation` -- DNA methylation
- `metainformant.epigenome.atacseq` -- chromatin accessibility
- `metainformant.epigenome.peak_calling` -- de novo peak calling from signal data
- `metainformant.epigenome.chromatin_state` -- chromatin state learning
