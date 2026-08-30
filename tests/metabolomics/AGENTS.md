# AGENTS.md — `MetaInformAnt/tests/metabolomics/`

tests for the metabolomics module.
Files (verified 2026-08-29): __init__.py, test_metabolomics_identification.py, test_metabolomics_io.py, test_metabolomics_pathways.py, test_metabolomics_visualization.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/metabolomics` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.