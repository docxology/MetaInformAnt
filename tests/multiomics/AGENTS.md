# AGENTS.md — `MetaInformAnt/tests/multiomics/`

tests for the multiomics module.
Files (verified 2026-08-29): __init__.py, test_multiomics_comprehensive.py, test_multiomics_integration.py, test_multiomics_methods.py, test_multiomics_pathways.py, test_multiomics_sample_mapping.py, test_multiomics_survival.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/multiomics` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.