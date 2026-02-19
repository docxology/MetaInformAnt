# EPIGENOME

## Overview
Epigenome analysis module for METAINFORMANT.

## 📦 Contents
- **[analysis/](analysis/)**
- **[assays/](assays/)**
- **[visualization/](visualization/)**
- **[workflow/](workflow/)**
- `[__init__.py](__init__.py)`

## 📊 Structure

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
```python
from metainformant.epigenome import ...
```
