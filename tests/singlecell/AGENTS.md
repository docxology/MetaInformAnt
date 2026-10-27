# AGENTS.md — `MetaInformAnt/tests/singlecell/`

tests for the singlecell module.
Files (verified 2026-10-27): __init__.py, test_singlecell_basic.py, test_singlecell_celltyping.py,
test_singlecell_clustering_diffusion_contracts.py, test_singlecell_clustering_trajectory.py, test_singlecell_data_doublet.py,
test_singlecell_differential.py, test_singlecell_diffusion_pseudotime.py, test_singlecell_dimensionality.py,
test_singlecell_docstring_exceptions.py, test_singlecell_hvg_seurat_v3.py, test_singlecell_input_contracts.py,
test_singlecell_leiden_fallback.py, test_singlecell_preprocessing.py, test_singlecell_preprocessing_consistency.py,
test_singlecell_randomness.py, test_singlecell_velocity.py, test_singlecell_visualization.py (17 test files plus `__init__.py`).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/singlecell -p no:cacheprovider` (verified 2026-10-27: 434 passed, 7 skipped).
  Every skip needs the Leiden/Louvain stack (`python-igraph`, `leidenalg`, `python-louvain`), which the venv does not install; the
  Leiden-available classes of `test_singlecell_leiden_fallback.py` have therefore not been executed in this venv.
- 2026-10-27 files: `test_singlecell_diffusion_pseudotime.py` (DPT against a curved trajectory, scale invariance, the paper's
  brute-force definition and the ARPACK path), `test_singlecell_input_contracts.py` (`inplace` contract of `select_hvgs`,
  `compute_neighbors`, `compute_diffusion_map`), `test_singlecell_leiden_fallback.py` (no silent k-means),
  `test_singlecell_preprocessing_consistency.py` (per-gene QC refresh, `uns` isolation, size-factor direction) and
  `test_singlecell_randomness.py` (no use of numpy's global random state).
- `test_singlecell_hvg_seurat_v3.py` pins the `seurat_v3` flavor of `select_hvgs`: the local-regression smoother against
  statsmodels' lowess (degree 1) and a brute-force loop (degrees 1 and 2), the statistic against an explicit per-gene loop of
  the documented formula, and the selection against genes planted with extra-Poisson variance. scanpy is not the oracle
  because its `seurat_v3` needs scikit-misc, which the venv does not install.
- `test_singlecell_docstring_exceptions.py` scans every `Raises:` section of the package with `ast`: a function whose checks
  raise `ValidationError` (not a `ValueError` or `TypeError` subclass) may not document `ValueError` or `TypeError` unless a
  library it calls raises it (two listed cases). It also pins what `load_count_matrix` raises for a missing file, an unknown
  format and a Matrix Market load without gene and barcode files.
- `test_singlecell_basic.py::TestIoUtilities` pins the `io.py` contract: the nine helpers for R and pySCENIC tools raise
  `NotImplementedError` with a message that names the function and creates no file or directory, and the offline paths of
  `fetch_atlas_datasets` (unknown names skipped, ArrayExpress accessions get only their directory).
- The 434 passed and 7 skipped add up to the 441 tests that `pytest --collect-only` finds in the directory.
Repo-wide policy: see the repository-root `AGENTS.md`.