# AGENTS.md — `MetaInformAnt/tests/ecology/`

tests for the ecology module.
Files (verified 2026-08-29): __init__.py, test_ecology_basic.py, test_ecology_comprehensive.py, test_ecology_functional.py, test_ecology_macroecology.py, test_ecology_ordination.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/ecology` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.