### Epigenome: Overview

Epigenetic modification analysis including DNA methylation patterns, histone modification (ChIP-seq) peaks, chromatin accessibility and chromatin states.

## Key Components

### DNA Methylation Analysis
Implemented in `metainformant.epigenome.assays.methylation`.
- **Methylation data loading**: Read per-site methylation calls from BEDgraph and coverage files or count tables (`load_methylation_bedgraph`, `load_methylation_cov`, `load_cpg_table`). Calling methylation from reads or alignments is not implemented in this module
- **Differential methylation**: Compare methylation patterns between conditions (`find_differentially_methylated_regions`)
- **CpG island screening**: Window scan that estimates CpG-island-like regions from methylation levels (`identify_cpg_islands`); it does not read sequence

### Chromatin Accessibility
- **ATAC-seq analysis**: Load, summarize and compare accessible-region peaks (`metainformant.epigenome.assays.atacseq`)
- **Peak calling**: Statistical detection of enriched regions from signal tracks with a Poisson model (`metainformant.epigenome.peak_calling.peak_detection`)
- **Motif scanning**: Search peak windows for IUPAC motif patterns in a genome FASTA (`assays.chipseq.find_motifs_in_peaks`)
- **Differential accessibility**: Compare peak sets between conditions (`atacseq.compare_atac_conditions`, `peak_detection.differential_peaks`)

### Histone Modifications
- **ChIP-seq peak analysis**: Load, filter, summarize, overlap and merge peaks (`metainformant.epigenome.assays.chipseq`)
- **Chromatin state learning**: Discover and label chromatin states from histone mark signal (`metainformant.epigenome.chromatin_state`)

## Integration with Genomics

Epigenetic data integrates with genomic annotation in two ways within this module:
- **Multi-assay integration**: `workflow.workflow.integrate_epigenome_results` relates methylation sites to ChIP-seq and ATAC-seq peaks
- **Functional annotation**: `chromatin_state.state_learning.interpret_states` labels chromatin states from histone mark patterns, and `compute_state_enrichment` tests the states against annotated bins

Not implemented here: connecting methylation or accessibility to gene expression (the `metainformant.multiomics` module accepts an epigenome table through `analysis.integration.integrate_omics_data(..., epigenome_data=...)`), linking to 3D genome organization, and cross-species comparison.

## Analysis Workflows

### Methylation Analysis Pipeline
```python-snippet
from metainformant.epigenome.assays import methylation

# Load methylation data (Bismark coverage files: chromosome, start, end,
# percentage, methylated reads, unmethylated reads)
treatment_data = methylation.load_methylation_cov("treatment.cov")
control_data = methylation.load_methylation_cov("control.cov")

# Summarize methylation levels and coverage
stats = methylation.calculate_methylation_statistics(treatment_data)
print(stats["total_sites"], stats["mean_methylation"])

# Find differentially methylated regions
differential = methylation.find_differentially_methylated_regions(
    treatment_data,
    control_data,
    delta_threshold=0.2,
    p_value_threshold=0.05,
)
```

### Chromatin Accessibility Pipeline
```python-snippet
from metainformant.epigenome.assays import atacseq

# Load ATAC-seq peaks (narrowPeak format)
control_peaks = atacseq.load_atac_peaks("control_peaks.narrowpeak")
treated_peaks = atacseq.load_atac_peaks("treated_peaks.narrowpeak")

# Summary statistics, including nucleosome fractions derived from peak lengths
stats = atacseq.calculate_atac_statistics(treated_peaks)
print(stats["total_peaks"], stats["nfr_peak_fraction"])

# Accessibility index for regions of interest, given as (chromosome, start, end)
regions = [("chr1", 1000, 5000), ("chr2", 200, 4000)]
accessibility = atacseq.calculate_chromatin_accessibility_index(treated_peaks, regions)

# Compare accessibility between conditions
comparison = atacseq.compare_atac_conditions(control_peaks, treated_peaks)
print(comparison["overlapping_peaks"], comparison["overlap_percentage"])
```

## Planned Extensions

Future epigenetic analysis modules will include:
- **DNA hydroxymethylation**: 5hmC analysis and mapping
- **RNA modifications**: epitranscriptomics analysis
- **3D genome structure**: Hi-C and chromatin conformation analysis
