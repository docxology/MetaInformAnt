# Performance: Spatial

## Benchmarks (Intel i7‑12700K, 32 GB, Python 3.12, metainformant 0.8.0)

All tests on a single Visium section (~5 200 spots × 36 500 genes, sparse matrix
~4 % fill).

| Operation | Median | 95th pct | Peak RAM | Notes |
|-----------|--------|----------|----------|-------|
| `load_visium()` | 2.1 s | 2.4 s | 180 MB | includes Space Ranger JSON parse |
| `neighbors(method='knn')` | 180 ms | 210 ms | 50 MB | sklearn KDTree, expression PCA‑50 |
| `neighbors(method='radius', r=100)` | 420 ms | 480 ms | 55 MB | radius search |
| Spatial Leiden (res=0.8) | 320 ms | 380 ms | 80 MB | igraph partition |
| Moran's I (single gene, 999 perm) | 350 ms | 400 ms | 10 MB | permutation loop |
| Local Moran (single gene, 99 perm) | 200 ms | 240 ms | 10 MB | |
| Stereoscope (5 k spots × 20 types, 400 epoch) | 6.2 min | — | 450 MB (CPU) | 90 s on CUDA A100 |
| Tangram (same) | 1.4 min | — | 380 MB | CUDA 45 s |
| `spatial_scatter()` (matplotlib) | 120 ms | 140 ms | 20 MB | figure render only |

## Scaling to larger sections

| Spots | Graph (knn) | Moran (999 perm) | Stereoscope (CPU) |
|-------|-------------|------------------|-------------------|
| 5 k   | 0.18 s      | 0.35 s           | 6 min |
| 20 k  | 0.75 s      | 1.4 s            | 25 min |
| 100 k | 4.2 s       | 7.8 s            | 2 h + |

Memory scales linearly with N_spots for dense graphs (~N² for radius method).

## Tuning knobs

- **Graph sparsification** — reduce `n_neighbors` to 10 for very large N; speeds
  clustering ×1.5 with minor accuracy loss.
- **Permutation count** — for discovery use `n_permutations=99`; for publication
  use `999`. Analytic Moran's I (no permutation) is available but less accurate
  for small N.
- **Deconvolution epochs** — early stopping monitors ELBO; typical convergence
  plateau at 300–400 epochs. Reduce to 200 for prototyping.
- **Device** — CUDA accelerates Stereoscope×8; Tangram×4 (memory bound).

## Bottleneck checklist

1. **Loader I/O** — raw matrix HDF5 read dominates initial load; place data on NVMe.
2. **Neighbour graph** — radius method O(N²); prefer KNN for N > 20 k.
3. **Deconvolution** — GPU memory can OOM with >100 cell types; lower `batch_size`.
4. **Permutation tests** — auto‑parallelise across genes via `joblib` if you have
   many cores: `morans_i(..., n_jobs=8)`.

## Profiling recipe

```bash
python -m cProfile -o spatial.prof myscript.py
snakeviz spatial.prof
```

Look for accumulated time in:
- `spatial/io/loader.py` — file parsing
- `spatial/analysis/neighbors.py` — scikit‑learn tree build
- `spatial/analysis/autocorrelation.py` — permutation loop

## Known limits

- Radius graph with N > 50 k spots may exceed RAM; use `method='knn'`.
- CUDA support requires `pytorch` + `faiss‑gpu`; CPU fallback works everywhere.
