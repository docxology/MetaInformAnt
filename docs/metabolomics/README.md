# Metabolomics

## Overview

Metabolomics analysis module for METAINFORMANT. Covers mass spectrometry data processing, metabolite identification, pathway mapping, and metabolite-gene integration.

## Contents

- **io/** - Mass spectrometry file reading (mzML, mzXML, CSV), format conversion
- **analysis/** - Metabolite identification, peak quantification, normalization
- **pathways/** - KEGG/Reactome pathway mapping, metabolite set enrichment
- **visualization/** - Volcano plots, PCA ordination, concentration heatmaps

## Architecture

```mermaid
graph TD
    subgraph "Metabolomics Module"
        IO[io/] --> |formats.py| FMT[mzML / mzXML / CSV I/O]

        AN[analysis/] --> |identification.py| ID[Metabolite Identification]
        AN --> |quantification.py| QT[Peak Quantification]
        AN --> |normalization.py| NM[Normalization & Scaling]

        PW[pathways/] --> |mapping.py| PM[Pathway Mapping]
        PW --> |enrichment.py| PE[Metabolite Set Enrichment]

        VZ[visualization/] --> |plots.py| PL[Volcano / PCA / Heatmap]
    end

    IO --> AN
    AN --> PW
    AN --> VZ
```

## Usage

```python
from metainformant.metabolomics.io import formats
from metainformant.metabolomics.analysis import identification, quantification, normalization
from metainformant.metabolomics.pathways import mapping, enrichment
from metainformant.metabolomics.visualization import plots
```
