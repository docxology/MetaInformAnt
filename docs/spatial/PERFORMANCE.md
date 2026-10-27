# Performance: Spatial

## Benchmarks (Apple silicon arm64, 14 logical cores, Python 3.12.15, NumPy 2.5.4, SciPy 1.18.1, scikit-learn 1.9.1, metainformant 0.4.0)

Measured on synthetic data: one section of 5 000 spots × 2 000 genes (four spatial
domains, SciPy CSR, ~4.5 % non-zero). Each figure is a median of 3–5 runs and the ranges
span three passes. The machine was shared with other jobs and timings moved by up to 4×
between passes, so the ranges below are orders of magnitude; re-measure on your own data.
`leidenalg` is not installed in the benchmark environment, so the Leiden row measures the
built-in greedy modularity fallback.

| Operation | Time | Notes |
|-----------|------|-------|
| `load_visium()` | 0.09–0.14 s | MEX matrix, `load_image=False` |
| `build_spatial_graph(method='knn', n_neighbors=6)` | 0.02–0.03 s | scikit-learn `NearestNeighbors`, then a Python loop over the edges |
| `build_spatial_graph(method='radius', radius=150)` | ~0.02 s | KD-tree pair search |
| `build_spatial_graph(method='delaunay')` | 0.06–0.09 s | |
| `spatial_cluster(method='leiden', resolution=0.8)` | 9–11 s | greedy-modularity fallback; includes densifying the matrix, z-scoring, PCA (15 components) and both graphs |
| `spatial_cluster(method='kmeans', n_clusters=8)` | 2–3 s | PCA + KMeans on the concatenated features |
| `spatial_weights_matrix(method='knn', k=6)` | 0.01–0.06 s | |
| `morans_i()` (single gene) | 1–2 ms | analytic statistic, no permutations |
| `local_morans_i()` (single gene) | ~5 ms | analytic |
| `detect_spatially_variable_genes()` (200 genes, dense input) | 0.25–1 s | Python loop, one `morans_i` per gene |
| `deconvolve_spots(method='nnls')` (20 types, 2 000 genes) | 3.5–15 s | one NNLS solve per spot |
| `deconvolve_spots(method='nmf')` (`alpha=0`) | 5–15 s | same per-spot solver |
| `plot_spatial_scatter()` | 0.15–0.3 s | 8 categorical groups; includes saving the PNG |

## Scaling to larger sections

Same synthetic data with 2 000 genes and 20 cell types, from one pass (medians of 3 runs,
except a single run for the 100 k NNLS row):

| Spots | `build_spatial_graph` (knn) | weights + `morans_i` (one gene) | `deconvolve_spots` NNLS |
|-------|-----------------------------|---------------------------------|-------------------------|
| 5 k   | 0.02 s | 0.01 s | 3.5 s |
| 20 k  | 0.07 s | 0.07 s | 15 s |
| 100 k | 0.6 s  | 0.4 s  | 100 s |

The graphs are SciPy sparse matrices (at most 2 × `n_spots` × `n_neighbors` stored
entries for `'knn'`), so graph memory grows linearly with the number of spots; for
`'radius'` it grows with the number of neighbour pairs inside the radius. Dense memory is
the larger concern: `spatial_cluster` and `deconvolve_spots` convert the expression matrix
to a dense float64 array (`n_spots × n_genes × 8` bytes, e.g. 100 000 × 2 000 = 1.6 GB),
and `spatial_cluster` makes another dense copy when it z-scores the genes.

## Tuning knobs

- **Gene subset** — pass highly variable genes (a few thousand) rather than the full
  matrix to `spatial_cluster` and `deconvolve_spots`; both densify their input, so memory
  and time scale with the number of genes.
- **Graph density** — `n_neighbors` (default 6) sets the size of the KNN graphs in
  `build_spatial_graph` and `spatial_cluster`; smaller values give sparser graphs.
- **Permutation and simulation counts** — `neighborhood_enrichment(n_permutations=1000)`
  and `ripley_k(n_simulations=99)` are the only Monte Carlo counts you can set; lower them
  for exploration. (`compute_ligand_receptor_interactions` also runs a permutation test,
  fixed at 100 permutations for every ligand–receptor pair and source/target cell-type
  pair.) Moran's I, Geary's C, local Moran's I and Getis–Ord G* are analytic and have no
  such setting.
- **Deconvolution** — `method='nnls'` is an exact per-spot solver (no epochs or early
  stopping); cost grows linearly with the number of spots, roughly linearly with the
  number of genes and faster than linearly with the number of cell types. Restrict the
  genes to markers or to the genes shared with the reference.
- **Device** — the module is CPU-only; there is no device option.

## Bottleneck checklist

1. **Loader I/O** — `load_visium` and `load_xenium` read the whole count matrix into
   memory (MEX through `scipy.io.mmread`, or HDF5 through h5py) and keep it as CSR;
   `load_merfish` parses `cell_by_gene.csv` with Python's `csv` module into a dense
   `float32` array (`n_cells × n_genes × 4` bytes, but the rows are first held as Python
   lists of floats, so the peak is roughly ten times that), so it is the most
   memory-hungry loader. Keep the data on fast local storage.
2. **Neighbour graph** — KNN construction loops over spots in Python, and the radius
   method collects every neighbour pair inside the radius in a Python set, so a radius much
   larger than the spot spacing costs a lot of memory and time. Prefer `'knn'` or a small
   radius.
3. **Densification** — `spatial_cluster` and `deconvolve_spots` convert the matrix to dense
   (see above); subset the genes first.
4. **Per-spot and per-gene loops** — `deconvolve_spots` solves one NNLS problem per spot and
   `detect_spatially_variable_genes` calls `morans_i` once per gene. Both are single-threaded
   Python loops; nothing in the module uses `joblib` or has an `n_jobs` option, so split the
   spots or genes across processes yourself to use several cores.

## Profiling recipe

```bash
python -m cProfile -o spatial.prof myscript.py
snakeviz spatial.prof
```

Look for accumulated time in:
- `spatial/io/visium.py`, `xenium.py`, `merfish.py` — file parsing
- `spatial/analysis/clustering.py` — graph construction, PCA and community detection
- `spatial/analysis/autocorrelation.py` — spatial weights and Moran's I
- `spatial/analysis/deconvolution.py` — the per-spot NNLS loop

## Known limits

- A radius graph with a large radius on many spots holds every neighbour pair in a Python
  set and may exceed RAM; use `method='knn'` or a smaller radius.
- `detect_spatially_variable_genes` accepts SciPy sparse input and densifies one gene at a
  time, so memory stays low, but it computes one Moran's I per gene: pass a gene subset
  for large panels.
- The module is CPU-only: there is no GPU backend and no PyTorch or FAISS dependency.
