# Spatial Transcriptomics

## Overview

Spatial transcriptomics analysis module for METAINFORMANT. Covers platform I/O (Visium, MERFISH, Xenium), spatial statistics, cell-cell communication, deconvolution, and scRNA-seq integration.

## Contents

- **io/** - Platform loaders for Visium, MERFISH, and Xenium data formats
- **analysis/** - Spatial autocorrelation, clustering, neighborhood analysis, deconvolution
- **communication/** - Ligand-receptor interaction scoring and communication networks
- **deconvolution/** - Advanced cell type deconvolution with reference profiles and niche ID
- **integration/** - scRNA-seq to spatial mapping, label transfer, gene imputation
- **niche/** - Tissue niche identification via spatial smoothing and K-Means clustering
- **visualization/** - Spatial scatter plots, tissue overlays, expression maps

## Architecture

```mermaid
graph TD
    subgraph "Spatial Module"
        IO[io/] --> |visium.py| VI[10x Visium Loader]
        IO --> |merfish.py| MF[MERFISH Loader]
        IO --> |xenium.py| XN[10x Xenium Loader]

        AN[analysis/] --> |autocorrelation.py| AC[Moran's I / Geary's C]
        AN --> |clustering.py| CL[Leiden / Spatial Domains]
        AN --> |neighborhood.py| NB[Neighborhood Analysis]
        AN --> |deconvolution.py| DC[Core Deconvolution]

        CM[communication/] --> |cell_communication.py| CC[Ligand-Receptor Scoring]

        DV[deconvolution/] --> |spatial_deconvolution.py| SD[Advanced Deconvolution]

        IG[integration/] --> |scrna_mapping.py| SM[scRNA-seq Mapping]

        VZ[visualization/] --> |plots.py| PL[Spatial Plots]
    end

    IO --> AN
    IO --> CM
    SM --> DV
```

## Usage

```python
from metainformant.spatial.io import visium, merfish, xenium
from metainformant.spatial.analysis import autocorrelation, clustering, neighborhood
from metainformant.spatial.communication import cell_communication
from metainformant.spatial.deconvolution import spatial_deconvolution
from metainformant.spatial.integration import scrna_mapping
```
