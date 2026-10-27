# Agent Directives: docs/singlecell

## Role

Documentation for the single-cell RNA-seq analysis module.

## Module Scope

Preprocessing (filtering, normalization), dimensionality reduction, clustering (Leiden, Louvain, k-means, hierarchical), differential expression, trajectory inference, RNA velocity, and cell-type annotation.

## Key Source Files

| Path | Description |
|------|-------------|
| `src/metainformant/singlecell/analysis/` | Clustering, dimensionality reduction, neighbor graph and trajectory inference |
| `src/metainformant/singlecell/data/` | `SingleCellData`, loading, QC, normalization and batch integration |
| `src/metainformant/singlecell/differential/` | Differential expression |
| `src/metainformant/singlecell/velocity/` | RNA velocity |
| `src/metainformant/singlecell/celltyping/` | Cell-type annotation |
| `src/metainformant/singlecell/visualization/` | Plotting functions |
| `src/metainformant/singlecell/io.py` | Dataset download, Salmon alevin, and the not-implemented R and pySCENIC helpers |

## Documentation Rules

- Every call quoted in a `python` or `python-snippet` block must bind to the real signature; run the documentation gate (`scripts/verify_documentation_code.py`) after editing.
- State behaviour that depends on optional packages (Leiden needs `networkx`, `igraph` and `leidenalg`; Louvain needs `python-louvain`; UMAP needs `umap-learn`).
- Do not describe the `run_*` helpers of `io.py` for R and pySCENIC tools as working wrappers: they raise `NotImplementedError`.
