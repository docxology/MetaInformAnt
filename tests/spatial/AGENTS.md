# AGENTS.md — `MetaInformAnt/tests/spatial/`

tests for the spatial module.
Files (verified 2026-10-26): __init__.py, test_spatial_analysis.py, test_spatial_autocorrelation.py,
test_spatial_communication.py, test_spatial_deconvolution_advanced.py, test_spatial_domains_deconv_niche.py,
test_spatial_integration.py, test_spatial_io.py, test_spatial_neighborhood.py,
test_spatial_optional_dependencies.py, test_spatial_svg.py, test_spatial_svg_sparse.py,
test_spatial_visualization.py (12 test files plus `__init__.py`).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/spatial -p no:cacheprovider` (verified 2026-09-05: 178 passed pre-review, 192 passed post-review with 14 added tests; 237 passed on 2026-10-26 with 45 further tests for sparse weight formats, correlation-method validation and optional-dependency fallbacks; 262 passed on 2026-10-26 with the 25 tests of `test_spatial_svg_sparse.py`).
- `test_spatial_optional_dependencies.py` runs a child interpreter with numpy, scipy, matplotlib and PIL blocked, to pin the
  runtime fallback (`None`-bound names, `ImportError` guidance) behind the `if not TYPE_CHECKING:` guards in `src/metainformant/spatial/`.
- `test_spatial_svg_sparse.py` checks that `detect_spatially_variable_genes` gives the dense result for SciPy sparse input
  (every sparse format and both orientations, plus the CSR matrix from the real `load_visium` reader) and uses `tracemalloc`
  to check that a whole-slide matrix is not densified.
Repo-wide policy: see the repository-root `AGENTS.md`.
