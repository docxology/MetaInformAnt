# AGENTS.md — `MetaInformAnt/tests/spatial/`

tests for the spatial module.
Files (verified 2026-08-29): __init__.py, test_spatial_analysis.py, test_spatial_autocorrelation.py, test_spatial_communication.py, test_spatial_deconvolution_advanced.py, test_spatial_integration.py, test_spatial_io.py, test_spatial_neighborhood.py … (9 files).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/spatial` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.