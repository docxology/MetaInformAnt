# Troubleshooting: Spatial

## Common errors & fixes

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `FileNotFoundError: Visium directory not found: ...` | wrong `path` | pass the Space Ranger output directory, the one that contains `filtered_feature_bc_matrix*` and `spatial/` |
| `FileNotFoundError: No expression matrix found in ...` (Visium) | neither `filtered_feature_bc_matrix.h5` nor `filtered_feature_bc_matrix/` in `path` | check the directory layout in [io.md](io.md) |
| `FileNotFoundError: No tissue_positions file found in .../spatial` | no `tissue_positions.csv` or `tissue_positions_list.csv` | restore the `spatial/` folder of the Space Ranger output |
| `ValueError: No barcodes match between the expression matrix (...) and the tissue positions (...)` (Visium) | none of the matrix barcodes has a tissue position: the matrix and `spatial/` come from different runs, or `in_tissue_only=True` dropped every spot of the positions file | use the matrix and the `spatial/` folder of one Space Ranger run; the message gives the counts and the first barcodes of both sides, and names `in_tissue_only` when it is set |
| `FileNotFoundError: No cells.csv(.gz) found in ...` (Xenium) | the cell centroids are required, so the loader does not guess them | add `cells.csv.gz` (or `cells.csv`) from the same Xenium run |
| `ValueError: Required columns (cell_id, x_centroid, y_centroid) not found in ...` (Xenium) | `cells.csv(.gz)` has other column names | rename the columns, or build the dataset yourself |
| `ValueError: N of M cells in the feature matrix have no row in cells.csv...` (Xenium) | `cells.csv(.gz)` and the feature matrix come from different runs or cell sets | use files of the same run; the error lists the first missing cell ids |
| `FileNotFoundError: No cell-by-gene expression file found in ...` (MERFISH) | none of `cell_by_gene.csv`, `cell_by_gene.tsv` or `counts.csv` is in `path` | add the matrix under one of those names (cells as rows, genes as columns; the delimiter follows the extension: tab for `.tsv`, comma otherwise) |
| `FileNotFoundError: No cell_metadata.csv/.tsv found in ...` (MERFISH) | the cell coordinates are required | add `cell_metadata.csv` (or `cell_metadata.tsv`) with a cell id column and `center_x`/`center_y` |
| `load_merfish(..., load_transcripts=True)` gives a dataset whose `transcript_spots` is `None` and logs the warning `No detected_transcripts file found; skipping transcript spots` | the loader only looks for `detected_transcripts.csv` and `transcripts.csv` in `path` (not for `.tsv` variants or other names) | rename the file, or read it yourself with `load_transcript_spots(file)`, which accepts any path and takes its delimiter from the extension, and assign the result to `dataset.transcript_spots` |
| `ValueError: No coordinate columns (x, y) found in ...` (MERFISH) | `cell_metadata.csv` has other column names | see the accepted names in [io.md](io.md) |
| `ValueError: N of M cells in the expression matrix have no row in cell_metadata.csv ...` (MERFISH) | some matrix cells have no metadata row | make the two files describe the same cells |
| The plotted tissue looks transposed or mirrored | an `n × 2` array in `[x, y]` order (Xenium, MERFISH) was plotted with the default `coordinate_order='row_col'` | pass `coordinate_order='xy'` (or `dataset.coordinate_order`); see [visualization.md](visualization.md#coordinate-order) |
| `ValueError: Unknown coordinate_order '...': use 'row_col' or 'xy'` | a typo in `coordinate_order` | use `'row_col'` or `'xy'` |
| `ValueError: graph_method='radius' needs the radius argument` | `spatial_cluster(graph_method='radius')` without `radius` | pass `radius=` in the units of the coordinates |
| `ValueError: n_iterations must be positive, or negative to iterate until convergence` | `leiden_clustering(n_iterations=0)` | use `-1` (until convergence) or a positive count |
| `spatial_domains` returns more domains than expected | with `n_domains=None` the gap statistic can prefer a few extra domains when each domain has only tens of spots | pass `n_domains`; `result.metadata['n_domains_selection']` is `'gap_statistic'` when the number was chosen and `'user'` when you gave it (see [clustering.md](clustering.md#choosing-the-number-of-domains)) |
| `ValueError: Unknown edge_correction '...': use 'isotropic' or 'none'` | a typo in `ripley_k(edge_correction=...)` | use `'isotropic'` (default) or `'none'` |
| `ValueError: radii must be a non-empty array of non-negative values` | `ripley_k` got an empty `radii` or a negative radius | pass at least one radius, all ≥ 0 |
| `ValueError: the points have no extent in x or in y, so there is no study window to correct for ...` | `ripley_k` with the default isotropic correction on points that lie on one horizontal or vertical line (or all coincide) | pass `edge_correction='none'` for such a pattern (see [neighborhood.md](neighborhood.md)) |
| `ValueError: No shared genes between spatial and reference data` | no overlap between the gene names passed to `deconvolve_spots` | harmonise the identifiers (Ensembl vs symbols) |
| `ValueError: gene_names is required: genes are matched to expression columns by name ...` | `compute_ligand_receptor_interactions` was given a plain matrix without `gene_names`, or `validate_deconvolution` was given `spatial_markers` and a plain matrix without names | pass `gene_names=` with one name per column, or pass a dataset or DataFrame that carries them (see [communication.md](communication.md#cell-type-level-scoring)) |
| `ValueError: None of the N ligand-receptor pairs has both of its genes among the M gene names (...)` | genes are matched by exact, case-sensitive name and the built-in database is written in human HGNC symbols, so mouse symbols, Ensembl ids or lower-case names match nothing | pass `lr_database=` written in the gene-name convention of the data, or rename the genes |
| `ValueError: N ligand/receptor name(s) label more than one expression column ...` (or `N marker name(s) ...` from `validate_deconvolution`) | a gene needed for scoring appears under several columns, so no column is the right one | make `gene_names` unique first, for example by collapsing the duplicate columns |
| `ValueError: gene_names has N entries but the expression matrix has M columns` (or `cell_types has N entries but the expression matrix has M cells`) | the names or labels do not line up with the matrix | pass one name per column and one cell-type label per row |
| `ValueError: spatial_markers needs the spatial expression matrix to correlate against ...` | `validate_deconvolution(spatial_markers=...)` without `expression` | pass `expression=` (spots × genes, in the order of the proportions) and `gene_names=`, or a dataset carrying both |
| `ValueError: expression has N rows but the deconvolution has M spots` | the `expression` of `validate_deconvolution` is not the matrix the proportions were estimated from | pass the same spots × genes matrix, in the same spot order |
| `ValueError: Unknown correlation_method '...': use 'pearson' or 'spearman'` | a typo in `validate_deconvolution(correlation_method=...)`; it is checked even when no markers are given | use `'pearson'` (default) or `'spearman'` (needs SciPy) |
| A cell type is missing from `marker_correlations` | its correlation cannot be computed; the reason is in `marker_correlations_skipped`: the type is not in the result, none of its markers is among `gene_names`, or its proportion or marker expression is constant across spots | compare the marker names with `gene_names` (exact, case-sensitive) and check the proportions; see [deconvolution.md](deconvolution.md) |
| `ImportError: h5py is required ...` | an `.h5` count matrix without `h5py` | `uv pip install h5py` |
| `ImportError: ... is required: uv pip install ...` | NumPy, SciPy, scikit-learn, matplotlib or Pillow is missing | install the named package |
| `MemoryError` or a very slow `build_spatial_graph(method='radius')` | a radius much larger than the spot spacing creates a huge number of neighbour pairs | use `method='knn'`, a smaller radius, or fewer spots |
| Leiden or Louvain results without `leidenalg` or `python-louvain` installed | the built-in greedy modularity routine ran; the result does not record it | install the packages (see [clustering.md](clustering.md#which-implementation-runs)) |
| P-values or niches differ between runs | `compute_ligand_receptor_interactions` and `niche_identification` default to `seed=None` | pass an integer `seed` (see [CONFIGURATION.md](CONFIGURATION.md#seeds)) |

## Step‑by‑step sanity check

```python
from metainformant.spatial.io.visium import load_visium

dataset = load_visium('sample/')
print(dataset.n_spots, dataset.n_genes)          # shape: spots x genes
print(dataset.coordinates[:5])                   # first 5 [pixel_row, pixel_col] pairs
print(dataset.platform, dataset.coordinate_order)  # 'visium' 'row_col'
```

`path` is the Space Ranger output directory (the one that holds `spatial/` and the matrix):

```
sample/
 ├─ filtered_feature_bc_matrix/       (or filtered_feature_bc_matrix.h5)
 │   ├─ barcodes.tsv.gz
 │   ├─ features.tsv.gz
 │   └─ matrix.mtx.gz
 └─ spatial/
     ├─ tissue_positions.csv          (or tissue_positions_list.csv)
     └─ scalefactors_json.json        (optional)
```

`load_visium` keeps only the spots that appear in both the matrix and the tissue-position records (with
`in_tissue_only=True`, the default, only spots marked in tissue), so `expression`, `barcodes`,
`tissue_positions` and `coordinates` always describe the same spots in the order of the matrix, and `n_spots`
counts those spots. When no barcode of the matrix has a position it raises `ValueError` instead of returning an
unfiltered matrix next to empty coordinates; check that the matrix and `spatial/` come from the same run.

## Deconvolution looks wrong

`deconvolve_spots` in `spatial.analysis.deconvolution` is a per-spot non-negative least-squares
solve against fixed reference profiles; it has no training loop, so there is no loss to plateau.
Typical causes of poor results:

1. **Unclean reference** — filter low-quality cells and normalise the single-cell reference before
   building profiles with `create_reference_profiles`.
2. **Mismatched gene names** — the reference and spatial genes must use the same identifiers
   (for example both symbols or both Ensembl ids); otherwise the shared-gene intersection is empty
   (`ValueError`) or tiny.
3. **Batch effects** — the spatial and reference data should be on a comparable scale (the same
   normalisation) before deconvolution.

Per-spot fit quality is in `deconv.residuals`; see [Deconvolution](deconvolution.md).

## Deeper diagnostics

The module logs through the standard `logging` module (one logger per module, for example
`metainformant.spatial.analysis.clustering`). Raise a logger to `DEBUG` or `INFO` to see what ran:

```python
import logging
logging.getLogger('metainformant.spatial.analysis.clustering').setLevel('INFO')
```

The INFO line `leidenalg unavailable (...); using greedy modularity fallback` shows that the fallback
clustering ran.

## Getting help

Include: platform (OS, Python, metainformant version), loader used, dataset size
(spots × genes), the full traceback, and a small sample of the input files if possible.
