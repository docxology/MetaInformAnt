# Specification: spatial

## 🎯 Scope

Spatial transcriptomics analysis: platform I/O (Visium, MERFISH, Xenium), spatial statistics, cell-cell communication, deconvolution, and scRNA-seq integration.

## 🧱 Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures

- **Sub-packages**: io, analysis, communication, deconvolution, integration, visualization
- **Key Concepts**: `SpatialDataset`, `TissuePosition`, Moran's I, ligand-receptor scoring

## 🔌 API Definition

### Exports

- `__init__.py`
