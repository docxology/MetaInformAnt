# AGENTS.md — `MetaInformAnt/tests/visualization/`

tests for the visualization module.
Files (verified 2026-08-29): __init__.py, test_cross_species.py, test_visualization.py, test_visualization_animations.py, test_visualization_basic.py, test_visualization_comprehensive.py, test_visualization_dimred.py, test_visualization_expression.py … (17 files).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/visualization` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.