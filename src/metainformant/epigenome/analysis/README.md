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
| `merge_tracks()` | Merge multiple tracks with union/intersection operations |
| `calculate_track_statistics()` | Compute coverage, mean signal, feature counts per track |
| `extract_track_region()` | Subset track to a specific chromosome:start-end region |
| `compare_tracks()` | Statistical comparison of two genomic tracks |
| `generate_track_report()` | Generate summary report for a track |

## Usage

```python
from metainformant.epigenome.analysis.tracks import load_bed_track, calculate_track_statistics

track = load_bed_track("peaks.bed", name="H3K27ac")
stats = calculate_track_statistics(track)
region = extract_track_region(track, "chr1", 1000000, 2000000)
```
