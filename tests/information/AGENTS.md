# AGENTS.md — `MetaInformAnt/tests/information/`

tests for the information-theory module.
Files (verified 2026-08-29): __init__.py, test_information_comprehensive.py, test_information_geometry_decomposition.py, test_information_integration.py, test_information_new_modules.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/information` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.