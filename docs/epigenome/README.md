# EPIGENOME

## Overview
Epigenome analysis module for METAINFORMANT.

## Contents
- **[index.md](index.md)**: module overview
- **[methylation.md](methylation.md)**: methylation levels, DMRs and CpG islands
- **[chip_seq.md](chip_seq.md)**: ChIP-seq peaks, overlap and enrichment
- **[atac_seq.md](atac_seq.md)**: ATAC-seq peaks and TF binding-site scans
- **[peak_calling.md](peak_calling.md)**: peak calling, merging and filtering
- **[chromatin.md](chromatin.md)**: chromatin states
- **[workflow.md](workflow.md)**: end-to-end workflows and integration

## Structure

```mermaid
graph TD
    subgraph "Epigenome Module"
        AS[assays/] --> |methylation.py| ME[Methylation Levels + DMR]
        AS --> |chipseq.py| CH[ChIP-seq Peak Handling]
        AS --> |atacseq.py| AT[ATAC-seq Accessibility]

        PK[peak_calling/] --> |peak_detection.py| PD[Poisson Peak Calling]

        CS[chromatin_state/] --> |state_learning.py| SL[ChromHMM-style GMM]

        AN[analysis/] --> |tracks.py| TR[Genomic Track Analysis]

        WF[workflow/] --> |workflow.py| WO[EpigenomeConfig Pipeline]

        VZ[visualization/] --> |visualization.py| VP[Epigenomic Plots]
    end
```

## Usage
Import module:
```python-snippet
from metainformant.epigenome import ...
```
