# SPEC: Epigenome Scripts

Orchestration of methylation, ChIP-seq, and ATAC-seq analysis pipelines.

## Pipelines

### 1. Methylation Calling
Processing bedGraph files to identify differentially methylated regions (DMRs).
- `analyze_methylation.py`: Standard workflow for DMR detection.

### 2. Chromatin Accessibility
Processing ATAC-seq and ChIP-seq data.
- `peak_calling.py`: Wrapper around standard peak callers with MetaInformAnt output standards.

## Standards

- **Output Management**: All genomic tracks must be saved to `output/epigenome/tracks/`.
- **Integration**: Use `metainformant.epigenome.assays` for domain-specific logic.
