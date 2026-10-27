# Spatial Transcriptomics Index

High‑level map.

## Getting started

- **[GETTING_STARTED.md](GETTING_STARTED.md)** — first‑time demo
- **[README.md](README.md)** — concise entry point and source links
- **[io.md](io.md)** — loaders per platform, their required files and errors, `coordinate_order`

## Core guides

- **[ARCHITECTURE.md](ARCHITECTURE.md)** — module design + data‑flow
- **[CONFIGURATION.md](CONFIGURATION.md)** — keyword arguments and defaults (there are no `spatial.*` config keys), seeds
- **[SPEC.md](SPEC.md)** — interface and boundary policy
- **[PERFORMANCE.md](PERFORMANCE.md)** — benchmarks & tuning
- **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** — error → fix table

## Analysis

- **[autocorrelation.md](autocorrelation.md)** — Moran, LISA, Geary, Getis‑Ord, variogram, spatially variable genes
- **[clustering.md](clustering.md)** — spatially‑constrained Leiden / Louvain / k‑means, graph options, spatial domains and the gap statistic that chooses their number
- **[neighborhood.md](neighborhood.md)** — Ripley's K, interaction index, neighborhood enrichment
- **[deconvolution.md](deconvolution.md)** — NNLS / NMF reference‑based deconvolution, marker‑based validation
- **[communication.md](communication.md)** — ligand‑receptor scoring (genes matched by exact name)
- **[niche.md](niche.md)** — niche discovery + annotation

## Visualisation

- **[visualization.md](visualization.md)** — matplotlib plots, `coordinate_order`

## Integration

- **[INTEGRATION.md](INTEGRATION.md)** — scRNA-seq reference mapping, label transfer and gene imputation

## API

- **[CAPABILITIES.md](CAPABILITIES.md)** — function and class catalogue by submodule
