# AGENTS.md — `MetaInformAnt/tests/spatial/`

tests for the spatial module.
Files (verified 2026-10-27): __init__.py, test_spatial_analysis.py, test_spatial_autocorrelation.py,
test_spatial_communication.py, test_spatial_communication_gene_names.py, test_spatial_deconvolution_advanced.py,
test_spatial_deconvolution_validation.py, test_spatial_domains_deconv_niche.py, test_spatial_domains_gap_statistic.py,
test_spatial_graph_clustering_contracts.py, test_spatial_integration.py, test_spatial_io.py,
test_spatial_loader_consistency.py, test_spatial_loaders_require_coordinates.py, test_spatial_neighborhood.py,
test_spatial_optional_dependencies.py, test_spatial_plot_coordinate_order.py, test_spatial_randomness.py,
test_spatial_ripley_edge_correction.py, test_spatial_svg.py, test_spatial_svg_sparse.py,
test_spatial_visualization.py (21 test files plus `__init__.py`).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/spatial -p no:cacheprovider` (verified 2026-09-05: 178 passed pre-review, 192 passed post-review with 14 added tests; 237 passed on 2026-10-26 with 45 further tests for sparse weight formats, correlation-method validation and optional-dependency fallbacks; 262 passed on 2026-10-26 with the 25 tests of `test_spatial_svg_sparse.py`).
- 396 passed and 1 skipped on 2026-10-27 (`env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/spatial -p no:cacheprovider`). The skip is
  the one test that needs `python-igraph` and `leidenalg`, which the venv does not install.
- The four 2026-10-27 files pin the graph, loader, plotting and randomness contracts: kNN edge weights depend on the
  endpoint distance only (not on mutual neighbors), `spatial_cluster(graph_method="radius")` takes a `radius`, the Xenium and
  MERFISH loaders refuse placeholder `[0, 0]` coordinates, every plot that draws coordinates honours `coordinate_order`
  (Visium `row_col`, Xenium and MERFISH `xy`), and no analysis path reads or advances numpy's global random state.
- The five later 2026-10-27 files pin the numerical contracts of the second review round, each against an independent
  computation: `test_spatial_communication_gene_names.py` (ligand-receptor genes matched to columns by exact name, a planted
  pair found, missing names or no present pair raise), `test_spatial_deconvolution_validation.py` (marker correlations equal
  scipy's Pearson and Spearman values, and change with the markers), `test_spatial_loader_consistency.py` (`load_visium`
  raises when no barcode matches; MERFISH reads `.csv` and `.tsv`), `test_spatial_ripley_edge_correction.py` (isotropic edge
  correction at every n, against a brute-force integration of the circle) and `test_spatial_domains_gap_statistic.py`
  (`spatial_domains` recovers 2 to 6 planted domains with the gap statistic).
  The 396 passed and 1 skipped add up to the 397 tests that `pytest --collect-only` finds in the directory.
- `test_spatial_optional_dependencies.py` runs a child interpreter with numpy, scipy, matplotlib and PIL blocked, to pin the
  runtime fallback (`None`-bound names, `ImportError` guidance) behind the `if not TYPE_CHECKING:` guards in `src/metainformant/spatial/`.
- `test_spatial_svg_sparse.py` checks that `detect_spatially_variable_genes` gives the dense result for SciPy sparse input
  (every sparse format and both orientations, plus the CSR matrix from the real `load_visium` reader) and uses `tracemalloc`
  to check that a whole-slide matrix is not densified.
Repo-wide policy: see the repository-root `AGENTS.md`.
