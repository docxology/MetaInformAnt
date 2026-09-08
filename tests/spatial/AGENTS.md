# AGENTS.md — `MetaInformAnt/tests/spatial/`

tests for the spatial module.
Files (verified 2026-09-05): __init__.py, test_spatial_analysis.py, test_spatial_autocorrelation.py,
test_spatial_communication.py, test_spatial_deconvolution_advanced.py, test_spatial_domains_deconv_niche.py,
test_spatial_integration.py, test_spatial_io.py, test_spatial_neighborhood.py, test_spatial_svg.py,
test_spatial_visualization.py (10 test files plus `__init__.py`).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/spatial -p no:cacheprovider` (verified 2026-09-05: 178 passed pre-review, 192 passed post-review with 14 added tests).
Repo-wide policy: see the repository-root `AGENTS.md`.