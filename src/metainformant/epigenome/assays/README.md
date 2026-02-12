# Epigenome Assays

Data loading, peak analysis, and quality metrics for the three core epigenomic assay types: ATAC-seq, ChIP-seq, and DNA methylation.

## Contents

| File | Purpose |
|------|---------|
| `atacseq.py` | ATAC-seq peak loading (narrowPeak), TSS enrichment, chromatin accessibility |
| `chipseq.py` | ChIP-seq peak loading, filtering, enrichment, motif discovery in peaks |
| `methylation.py` | Methylation BedGraph/cov parsing, DMR detection, CpG island identification |

## Key Functions

| Function | Description |
|----------|-------------|
| `load_atac_peaks()` | Parse ATAC-seq narrowPeak file into ATACPeak list |
| `calculate_atac_statistics()` | Peak width distribution, signal stats, FRiP |
| `identify_tss_enrichment()` | Score ATAC signal enrichment around TSS sites |
| `find_tf_binding_sites()` | Identify transcription factor motifs within peaks |
| `load_chip_peaks()` | Parse ChIP-seq narrowPeak file into ChIPPeak list |
| `filter_peaks_by_score()` | Threshold-based peak filtering with optional top-N |
| `calculate_peak_enrichment()` | Fold enrichment of peaks in genomic regions |
| `find_motifs_in_peaks()` | Motif scanning within ChIP-seq peak sequences |
| `load_methylation_bedgraph()` | Parse methylation BedGraph with coverage filtering |
| `find_differentially_methylated_regions()` | Detect DMRs between conditions |
| `identify_cpg_islands()` | Locate CpG islands by GC content and CpG ratio |
| `calculate_methylation_entropy()` | Shannon entropy of methylation levels per region |

## Usage

```python
from metainformant.epigenome.assays.chipseq import load_chip_peaks, calculate_peak_statistics
from metainformant.epigenome.assays.methylation import load_methylation_bedgraph

peaks = load_chip_peaks("H3K4me3.narrowPeak")
stats = calculate_peak_statistics(peaks)
meth = load_methylation_bedgraph("bismark.bedGraph", min_coverage=5)
```
