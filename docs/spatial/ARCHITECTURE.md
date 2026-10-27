# Architecture: Spatial

The spatial module is a set of independent sub-packages that exchange plain arrays and small
dataclasses. It needs no AnnData object, has no configuration layer (see
[CONFIGURATION.md](CONFIGURATION.md)) and keeps no global state: every setting is a keyword
argument and every random draw comes from a local generator.

## Sub-packages

1. **I/O** (`spatial.io`) — one module per platform. `visium.load_visium`, `xenium.load_xenium` and
   `merfish.load_merfish` return `SpatialDataset`, `XeniumDataset` and `MERFISHDataset`
   (expression matrix, `n × 2` coordinates, gene names, identifiers). See [io.md](io.md).
2. **Analysis** (`spatial.analysis`) — algorithms on arrays:
   - `autocorrelation` — Moran's I, Geary's C, local Moran's I (LISA), Getis-Ord G*, variogram,
     spatial weights
   - `clustering` — spatial graphs, Leiden / Louvain clustering with a built-in greedy fallback,
     `spatial_cluster`, `spatial_domains`
   - `neighborhood` — Ripley's K, interaction matrix, neighborhood enrichment, niche detection,
     spatial ligand-receptor scores
   - `deconvolution` — NNLS / NMF deconvolution against reference profiles
   - `spatially_variable_genes` — per-gene Moran's I ranking
3. **Communication** (`spatial.communication`) — cell-type level ligand-receptor scoring, networks and
   patterns.
4. **Deconvolution** (`spatial.deconvolution`) — dictionary-based reference profiles, cell-type
   mapping, validation and niche identification.
5. **Niche** (`spatial.niche`) — `identify_niches` on cell-type proportions.
6. **Integration** (`spatial.integration`) — `scrna_mapping`: label transfer from scRNA-seq,
   correlation mapping and gene imputation. See [INTEGRATION.md](INTEGRATION.md).
7. **Visualization** (`spatial.visualization`) — matplotlib plots that save to `output_path` and
   return the `Figure`. See [visualization.md](visualization.md).

There is no batch correction, section registration, interactive (plotly or bokeh) plotting or YAML
configuration in the module.

## Data flow

```
Loader  →  dataset (expression, coordinates)
              ↓
   analysis / communication / deconvolution / niche / integration
   (arrays in, result dataclasses or dictionaries out)
              ↓
   visualization (coordinates + values, or a dataset; figure files out)
```

- Analysis functions are pure with respect to their inputs: they return a result object (for example
  `MoransIResult`, `SpatialClusterResult`, `DeconvolutionResult`) and modify nothing in place.
- Spatial graphs (`build_spatial_graph`, `spatial_weights_matrix`) are SciPy sparse matrices that you
  pass on explicitly (for example to `compute_interaction_matrix`); nothing is cached on the dataset.
- The coordinate column order is a property of the dataset class (`coordinate_order`: `"row_col"` for
  Visium, `"xy"` for Xenium and MERFISH); plotting functions read it from a dataset, or take it as an
  argument for bare arrays.

## Dependencies between sub-packages

The sub-packages hardly import each other. There are three links: `analysis.spatially_variable_genes`
uses `analysis.autocorrelation`; `integration.scrna_mapping.map_scrna_to_spatial` imports
`analysis.deconvolution.create_reference_profiles` inside the function when `method='correlation'`; and
`communication` and `deconvolution` both import the private helper module `spatial/_gene_names.py`.
That module holds the gene-name matching they share: it splits a matrix, DataFrame or dataset into
matrix and names, and looks up each ligand, receptor or marker by exact name (see
[communication.md](communication.md#cell-type-level-scoring)). It imports without NumPy (the two helpers
that need it load it when called) and is not part of the public API.

No analysis code imports `visualization`, and `visualization` reads datasets only through attributes
(`expression`, `coordinates`, `gene_names`, `coordinate_order`). The numerical packages are optional at
import time: a missing NumPy, SciPy, scikit-learn or matplotlib leaves the name unbound and the function
that needs it raises an `ImportError` with an install hint. `h5py`, Pillow, `leidenalg` and
`python-louvain` are optional extras.
