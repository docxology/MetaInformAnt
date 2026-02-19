# SINGLECELL

## Overview
Single-cell analysis module for METAINFORMANT.

## 📦 Contents
- **[analysis/](analysis/)**
- **[data/](data/)**
- **[visualization/](visualization/)**
- `[__init__.py](__init__.py)`

## 📊 Structure

```mermaid
graph TD
    subgraph "Single-Cell Module"
        DA[data/] --> |preprocessing.py| SC[SingleCellData + QC]
        DA --> |integration.py| BI[Batch Correction: BBKNN, MNN]

        AN[analysis/] --> |clustering.py| CL[Leiden, Louvain, K-means]
        AN --> |dimensionality.py| DR[PCA, UMAP, t-SNE]
        AN --> |trajectory.py| TJ[Diffusion Pseudotime]

        CT[celltyping/] --> |annotation.py| CA[Cell Type Annotation]
        DE[differential/] --> |expression.py| DX[Differential Expression]
        VE[velocity/] --> |rna_velocity.py| RV[RNA Velocity]

        VZ[visualization/] --> |visualization.py| VP[UMAP, Trajectory Plots]
    end
```

## Usage
Import module:
```python
from metainformant.singlecell import ...
```
