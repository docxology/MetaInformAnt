# Epigenome Analysis

Genomic track loading, manipulation, and statistics for BED, BedGraph, and BigWig formats used in epigenomic data analysis.

## Contents

| File | Purpose |
|------|---------|
| `tracks.py` | GenomicTrack data model, BED/BedGraph/BigWig I/O, track merging and statistics |

## Key Functions

| Function | Description |
|----------|-------------|
| `load_bed_track()` | Parse BED file into a GenomicTrack object |
| `load_bedgraph_track()` | Parse BedGraph signal file into GenomicTrack |
| `load_bigwig_track()` | Read BigWig continuous signal track |
| `save_bed_track()` | Write GenomicTrack features to BED format |
| `save_bedgraph_track()` | Export GenomicTrack as BedGraph |
| `merge_tracks()` | Merge multiple tracks with `union`, `intersection`, `mean` or `max` (an unknown operation raises `ValueError`); the input tracks are not modified |
| `calculate_track_statistics()` | Compute coverage, mean signal, feature counts per track |
| `extract_track_region()` | Subset track to a specific chromosome:start-end region |
| `compare_tracks()` | Statistical comparison of two genomic tracks |
| `generate_track_report()` | Generate summary report for a track |

## Merging tracks

`merge_tracks(tracks, operation="union")` returns a new `GenomicTrack` and leaves the inputs unchanged; chromosomes are processed in sorted order.

| Operation | Result |
|-----------|--------|
| `union` | Features of any track that overlap or touch are merged into one feature spanning them, with the highest value |
| `mean`, `max` | Features that overlap or touch, directly or through a chain of such features, form a group that becomes one feature spanning it, with the mean or the highest value of its features |
| `intersection` | Only the regions covered by every track are kept, so a chromosome missing from any track has none. Within a track the value over a region is the highest value among its features covering it, and an intersection region takes the lowest of those values across the tracks. Touching regions with the same value are joined, and features keep only `start`, `end` and `value` |

Any other operation raises `ValueError`.

## Usage

```python
from metainformant.epigenome.analysis.tracks import (
    calculate_track_statistics,
    extract_track_region,
    load_bed_track,
    merge_tracks,
)

track = load_bed_track("peaks.bed", name="H3K27ac")
stats = calculate_track_statistics(track)
region = extract_track_region(track, "chr1", 1000000, 2000000)

# Regions covered by both tracks, valued at the lower of the two
other = load_bed_track("other_peaks.bed", name="H3K4me3")
shared = merge_tracks([track, other], operation="intersection")
```
