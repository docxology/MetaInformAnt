# ATAC-seq Analysis

Assay for Transposase-Accessible Chromatin sequencing analysis. Includes peak management, fragment size distribution, nucleosome positioning, TSS enrichment, transcription factor binding site detection, and multi-condition comparisons.

## Key Concepts

**ATAC-seq** uses hyperactive Tn5 transposase to probe open chromatin. Fragment sizes reflect nucleosome organization: sub-nucleosomal fragments (<147 bp) indicate nucleosome-free regions (NFRs), while mono-nucleosomal fragments (~147--294 bp) indicate positioned nucleosomes.

**TSS enrichment** measures signal accumulation around transcription start sites, serving as a key quality control metric. High TSS enrichment indicates successful capture of regulatory elements.

**Nucleosome-free regions (NFRs)** are the primary signal in ATAC-seq and correspond to active regulatory elements including promoters, enhancers, and transcription factor binding sites.

## Data Model

### `ATACPeak`

Dataclass representing an ATAC-seq accessible region.

| Field | Type | Description |
|-------|------|-------------|
| `chrom` | `str` | Chromosome |
| `start` | `int` | Start coordinate |
| `end` | `int` | End coordinate |
| `name` | `str` | Peak identifier |
| `score` | `float` | Accessibility score |
| `fold_enrichment` | `float` | Fold enrichment over background |
| `q_value` | `float` | FDR-adjusted significance |
| `summit` | `int` | Position of maximum accessibility |

## Function Reference

### `load_peaks(file_path, format="narrowPeak") -> List[ATACPeak]`

Load ATAC-seq peaks from narrowPeak or BED format.

### `save_peaks(peaks, file_path, format="narrowPeak") -> None`

Write peaks to disk in the specified format.

### `peak_statistics(peaks) -> Dict`

Compute summary statistics for a peak set including total count, coverage, width distribution, and per-chromosome breakdown.

### `calculate_nucleosome_fractions(fragment_sizes) -> Dict`

Classify fragments into nucleosome-positioning categories:
- `nfr_fraction`: nucleosome-free (< 147 bp)
- `mono_nucleosome`: mono-nucleosomal (147--294 bp)
- `di_nucleosome`: di-nucleosomal (294--441 bp)
- `tri_nucleosome`: tri-nucleosomal (> 441 bp)

### `tss_enrichment(peaks, tss_positions, window=2000) -> Dict`

Calculate TSS enrichment score by measuring peak density around annotated transcription start sites within the specified window. Returns enrichment score, total TSS tested, and overlapping count.

### `find_tf_binding_sites(peaks, motif_database, sequences) -> List[Dict]`

Scan peak sequences for transcription factor binding motifs. Each result includes the TF name, motif match score, position, and strand.

### `compare_conditions(peaks_a, peaks_b, labels=None) -> Dict`

Compare accessibility between two conditions. Returns shared peaks, condition-specific peaks, differential statistics, and a Jaccard similarity index.

## Usage Examples

```python
from metainformant.epigenome import (
    ATACPeak,
    load_peaks as load_atac_peaks,
    peak_statistics as atac_stats,
    calculate_nucleosome_fractions,
    tss_enrichment,
    compare_conditions,
)

# Load ATAC-seq peaks
peaks = load_atac_peaks("atac_peaks.narrowPeak")
stats = atac_stats(peaks)
print(f"Accessible regions: {stats['total_peaks']}")

# Fragment size analysis
fragments = [120, 80, 200, 350, 160, 95, 180, 250, 400]
fractions = calculate_nucleosome_fractions(fragments)
print(f"NFR fraction: {fractions['nfr_fraction']:.2%}")
print(f"Mono-nucleosomal: {fractions['mono_nucleosome']:.2%}")

# TSS enrichment QC
tss_list = [("chr1", 1000, "+"), ("chr1", 50000, "-")]
tss_score = tss_enrichment(peaks, tss_list, window=2000)
print(f"TSS enrichment: {tss_score['enrichment_score']:.2f}")

# Condition comparison
diff = compare_conditions(treatment_peaks, control_peaks,
                          labels=["treatment", "control"])
print(f"Jaccard similarity: {diff['jaccard']:.3f}")
```

## Configuration

Environment variable prefix: `EPI_`

## Related Modules

- `metainformant.epigenome.chipseq` -- histone modification peaks
- `metainformant.epigenome.peak_calling` -- de novo peak calling
- `metainformant.epigenome.chromatin_state` -- chromatin state annotation
- `metainformant.epigenome.workflow` -- integrated analysis pipelines
