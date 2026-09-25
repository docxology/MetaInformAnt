# ATAC-seq Analysis

Assay for Transposase-Accessible Chromatin sequencing analysis. Includes peak management, fragment size distribution, nucleosome positioning, TSS enrichment, transcription factor binding site detection, and multi-condition comparisons.

## Key Concepts

**ATAC-seq** uses hyperactive Tn5 transposase to probe open chromatin. Fragment sizes reflect nucleosome organization: sub-nucleosomal fragments (<147 bp) indicate nucleosome-free regions (NFRs), while mono-nucleosomal fragments (~147--294 bp) indicate positioned nucleosomes.

**TSS enrichment** measures signal accumulation around transcription start sites, serving as a key quality control metric. High TSS enrichment indicates successful capture of regulatory elements.

**Nucleosome-free regions (NFRs)** are the primary signal in ATAC-seq and correspond to active regulatory elements including promoters, enhancers, and transcription factor binding sites.

## Data Model

### `ATACPeak`

Class representing an ATAC-seq accessible region.

| Field | Type | Description |
|-------|------|-------------|
| `chromosome` | `str` | Chromosome name |
| `start` | `int` | Peak start position |
| `end` | `int` | Peak end position |
| `score` | `float` | Peak score |
| `strand` | `str` | DNA strand |
| `signal_value` | `float` | Signal value |
| `p_value` | `float \| None` | P-value |
| `q_value` | `float \| None` | Q-value (FDR) |
| `summit` | `int \| None` | Peak summit position |

## Function Reference

### `load_atac_peaks(path, format="narrowpeak") -> List[ATACPeak]`

Load ATAC-seq peaks from narrowPeak, broadPeak, or BED format.

### `save_atac_peaks(peaks, path, format="narrowpeak") -> None`

Write peaks to disk in the specified format.

### `calculate_atac_statistics(peaks) -> Dict`

Compute summary statistics for a peak set: total count, length distribution (mean/median/min/max/std), score distribution, and per-chromosome signal breakdown.

### `calculate_atac_specific_metrics(peaks) -> Dict`

ATAC-seq-specific QC computed from peak lengths: `nfr_peak_fraction` (50--150 bp), `mononucleosome_peak_fraction` (150--250 bp), `dinucleosome_peak_fraction` (250--350 bp), and periodicity scores around expected nucleosome periods (147/200/300 bp). Takes a peak list, not raw fragment sizes.

### `identify_tss_enrichment(peaks, tss_positions, window_size=2000) -> Dict`

Measure peak density around annotated transcription start sites. `tss_positions` maps chromosome names to lists of TSS positions. Returns `total_tss`, `enriched_tss`, `enrichment_ratio`, `expected_ratio`, `fold_enrichment`, and `window_size`.

### `find_tf_binding_sites(peaks, tf_motifs, genome_fasta=None) -> Dict`

Scan accessible regions for transcription factor binding motifs. `tf_motifs` maps TF names to motif sequences; pass `genome_fasta` to analyze real sequence context.

### `calculate_chromatin_accessibility_index(peaks, genomic_regions) -> Dict`

Accessibility index for specific genomic regions given as `(chromosome, start, end)` tuples.

### `compare_atac_conditions(condition1_peaks, condition2_peaks) -> Dict`

Compare accessibility between two conditions. Returns `condition1_total`, `condition2_total`, `overlapping_peaks`, condition-only counts, and `overlap_percentage`.

## Usage Examples

```python-snippet
from metainformant.epigenome.assays.atacseq import (
    ATACPeak,
    load_atac_peaks,
    calculate_atac_statistics,
    calculate_atac_specific_metrics,
    identify_tss_enrichment,
    find_tf_binding_sites,
    compare_atac_conditions,
)

# Load ATAC-seq peaks
peaks = load_atac_peaks("atac_peaks.narrowpeak")
stats = calculate_atac_statistics(peaks)
print(f"Accessible regions: {stats['total_peaks']}")

# Nucleosome-positioning QC from peak lengths
metrics = calculate_atac_specific_metrics(peaks)
print(f"NFR peak fraction: {metrics['nfr_peak_fraction']:.2%}")
print(f"Mono-nucleosomal: {metrics['mononucleosome_peak_fraction']:.2%}")

# TSS enrichment QC (chromosome -> TSS positions)
tss_positions = {"chr1": [1000, 50000], "chr2": [25000]}
tss = identify_tss_enrichment(peaks, tss_positions, window_size=2000)
print(f"Fold enrichment: {tss['fold_enrichment']:.2f}")

# TF binding site scan
tf_motifs = {"CTCF": "CCGCGNGGNGGCAG", "SOX2": "CATTGTT"}
binding = find_tf_binding_sites(peaks, tf_motifs)

# Condition comparison
diff = compare_atac_conditions(treatment_peaks, control_peaks)
print(f"Overlapping peaks: {diff['overlapping_peaks']}")
print(f"Overlap percentage: {diff['overlap_percentage']:.1f}%")
```

## Configuration

No dedicated configuration file or environment prefix is used by the ATAC-seq module; call the functions above directly from Python.

## Related Modules

- `metainformant.epigenome.assays.chipseq` -- histone modification peaks
- `metainformant.epigenome.peak_calling` -- de novo peak calling
- `metainformant.epigenome.chromatin_state` -- chromatin state annotation
- `metainformant.epigenome.workflow` -- integrated analysis pipelines
