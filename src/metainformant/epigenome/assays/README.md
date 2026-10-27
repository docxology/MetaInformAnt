# Epigenome Assays

Data loading, peak analysis, and quality metrics for the three core epigenomic assay types: ATAC-seq, ChIP-seq, and DNA methylation.

## Contents

| File | Purpose |
|------|---------|
| `_peak_format.py` | Internal narrowPeak / broadPeak / BED reader and writer shared by `atacseq.py` and `chipseq.py`: the column layouts, the `-log10` conversion of p- and q-values, and the layout names |
| `atacseq.py` | ATAC-seq peak loading and saving, TSS enrichment, chromatin accessibility, TF motif scan of peak sequence |
| `chipseq.py` | ChIP-seq peak loading and saving, filtering, overlap and merging, enrichment, IUPAC motif scan of peak windows |
| `methylation.py` | Methylation BedGraph/cov parsing, DMR detection, sequence-based CpG island calling |

## Key Functions

| Function | Description |
|----------|-------------|
| `load_atac_peaks()` | Parse an ATAC-seq narrowPeak, broadPeak or BED file into an ATACPeak list (p- and q-value columns hold `-log10` values and are read as probabilities) |
| `save_atac_peaks()` | Write ATACPeak objects in the same layouts; every line is formatted before the file is opened |
| `calculate_atac_statistics()` | Peak width distribution and signal statistics |
| `identify_tss_enrichment()` | Score ATAC signal enrichment around TSS sites |
| `find_tf_binding_sites()` | IUPAC consensus motif scan of peak sequence on both strands (sequence from `genome_fasta` or `peak_sequences`, which is required); site counts, expected sites from base composition and enrichment per TF. No PWM scoring |
| `load_chip_peaks()` | Parse a ChIP-seq narrowPeak, broadPeak or BED file into a ChIPPeak list |
| `save_chip_peaks()` | Write ChIPPeak objects in the same layouts |
| `filter_peaks_by_score()` | Threshold-based peak filtering with optional top-N |
| `find_overlapping_peaks()` | Every overlapping pair of peaks between two sets, found exactly per chromosome (nested peaks included) and returned in a deterministic order |
| `merge_overlapping_peaks()` | Merge overlapping or nearby peaks per chromosome, keeping the lowest p- and q-values |
| `calculate_peak_enrichment()` | Genome coverage of the peaks and, given an expected peak count, fold enrichment and the Poisson upper-tail p-value (with its log10, which stays finite when the p-value underflows) |
| `find_motifs_in_peaks()` | IUPAC motif scan of peak windows against a genome FASTA (both strands, overlapping occurrences included; real sequence matching) |
| `load_methylation_bedgraph()` | Parse methylation BedGraph with levels written as fractions or percentages (`level_scale`) |
| `find_differentially_methylated_regions()` | Detect DMRs between conditions (per-site Yates-corrected z-test or Fisher's exact test on the read counts) |
| `find_cpg_islands()` | CpG island calling from DNA sequence with the Gardiner-Garden & Frommer criteria |
| `identify_cpg_islands()` | `find_cpg_islands` plus the number and mean level of the methylation sites inside each island (sequence required) |
| `calculate_methylation_entropy()` | Shannon entropy of methylation levels per region |

## Usage

```python
from metainformant.epigenome.assays.chipseq import load_chip_peaks, calculate_peak_statistics
from metainformant.epigenome.assays.methylation import load_methylation_bedgraph

peaks = load_chip_peaks("H3K4me3.narrowPeak")
stats = calculate_peak_statistics(peaks)
meth = load_methylation_bedgraph("bismark.bedGraph", min_coverage=5)
```
