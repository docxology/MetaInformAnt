# Specification: spatial

## Scope

Spatial transcriptomics analysis: platform I/O (Visium, MERFISH, Xenium), spatial statistics, cell-cell communication, deconvolution, and scRNA-seq integration.

## Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## Data Structures

- **Sub-packages**: io, analysis, communication, deconvolution, integration, niche, visualization
- **Key Concepts**: `SpatialDataset`, `XeniumDataset`, `MERFISHDataset`, `TissuePosition`, Moran's I, ligand-receptor scoring

## API Definition

### Exports

- `__init__.py`

### Conventions

- **Coordinates.** Each dataset class declares the column order of its `coordinates` in a `coordinate_order` class constant: `"row_col"` for `SpatialDataset` (Visium, `[pixel_row, pixel_col]`) and `"xy"` for `XeniumDataset` and `MERFISHDataset` (`[x, y]`). The plot functions take a `coordinate_order` argument (`"row_col"` or `"xy"`, otherwise `ValueError`); the ones that take a dataset default to the dataset's own value, the ones that take bare arrays default to `"row_col"`.
- **Loaders.** `load_visium` filters the matrix, barcodes, tissue positions and coordinates to the spots found in both the matrix and the positions file, in matrix order, and raises `ValueError` when no barcode matches. `load_xenium` and `load_merfish` raise `FileNotFoundError` (missing matrix, `cells.csv(.gz)` or `cell_metadata.csv`/`.tsv`) or `ValueError` (missing id or coordinate columns, or a cell of the matrix without a row); no cell receives placeholder `[0, 0]` coordinates. MERFISH tables take their delimiter from the file extension (tab for `.tsv` and `.tab`, comma otherwise); `load_merfish` finds `cell_by_gene.csv`, `cell_by_gene.tsv` or `counts.csv`, `cell_metadata.csv` or `.tsv` and, with `load_transcripts=True`, only `detected_transcripts.csv` or `transcripts.csv`.
- **Gene names.** `compute_ligand_receptor_interactions` and `validate_deconvolution(spatial_markers=...)` match ligand, receptor and marker genes to expression columns by exact, case-sensitive name. The names come from `gene_names=` or from a dataset or DataFrame (explicit names win); a plain matrix without names raises `ValueError`, as does a needed gene that labels more than one column. Ligand-receptor pairs with an absent gene are skipped and counted in `summary` (`n_pairs_in_database`, `n_pairs_scored`, `n_pairs_skipped`), and the built-in database is written in human HGNC symbols. The shared helpers are in the private module `spatial/_gene_names.py`.
- **Ripley's K.** `ripley_k` applies Ripley's isotropic edge correction at every number of points, with the bounding box of `points` as the window (`edge_correction="isotropic"`, the default; `"none"` gives every pair weight 1; anything else raises `ValueError`). `area` sets only the intensity.
- **Domains.** `spatial_domains(n_domains=None)` chooses between 2 and 15 domains with the gap statistic (Tibshirani, Walther and Hastie, 2001) and records `metadata["n_domains_selection"]` as `"gap_statistic"` or `"user"`.
- **Deconvolution validation.** `validate_deconvolution` correlates, across spots, the estimated proportion of each cell type with the mean expression of its markers (`correlation_method="pearson"` by default, `"spearman"` on request) and needs `expression` and gene names whenever markers are given. A type it cannot score is listed in `marker_correlations_skipped` with the reason, never reported as a number.
- **Spatial graphs.** `build_spatial_graph(method="knn")` weighs an edge `1 / (1 + distance)` and symmetrises with the entrywise maximum; `"delaunay"` and `"radius"` edges weigh 1, and `method="radius"` requires `radius`. `spatial_cluster` accepts `radius` for `graph_method="radius"` and raises `ValueError` without it.
- **Clustering fallbacks.** `leiden_clustering` and `louvain_clustering` use `leidenalg` / `python-louvain` and raise `ImportError` when those are not importable, unless `allow_greedy_fallback=True` selects the built-in greedy modularity routine (with a `UserWarning`). `spatial_cluster(..., allow_greedy_fallback=...)` records what ran: `SpatialClusterResult.method` is `"greedy_modularity"` after a fallback, `metadata["requested_method"]` the method asked for. `leiden_clustering(n_iterations=-1)` runs until convergence, `0` raises `ValueError`.
- **Undefined statistics.** Moran's I, Geary's C, local Moran's I, Getis-Ord G* and the variogram return NaN (with `status` on the global results) where they are undefined (constant values, a weights matrix without edges, too few observations), never `I = 0`, `C = 1` or `p = 1`; the spatially variable genes' Benjamini-Hochberg correction runs over the genes that have a p-value. Permutation p-values (`neighborhood_enrichment`, `compute_ligand_receptor_interactions`) are `(b + 1) / (N + 1)`, never 0. Unknown quality or volume values in the Xenium and MERFISH loaders are NaN, and `read_transcripts` refuses a quality filter on a file without a quality column.
- **Randomness.** Every random draw uses a local generator; numpy's global random state is never read or advanced. Most functions take `seed` (default 42; `random_state` for `identify_niches`); `compute_ligand_receptor_interactions` and `niche_identification` default to `seed=None`, which draws from operating-system entropy and is not reproducible.
