# AGENTS.md — `MetaInformAnt/tests/singlecell/`

tests for the singlecell module.
Files (verified 2026-08-29): __init__.py, test_singlecell_basic.py, test_singlecell_celltyping.py, test_singlecell_differential.py, test_singlecell_dimensionality.py, test_singlecell_preprocessing.py, test_singlecell_velocity.py, test_singlecell_visualization.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/singlecell` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.