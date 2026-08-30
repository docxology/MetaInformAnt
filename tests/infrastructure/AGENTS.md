# AGENTS.md — `MetaInformAnt/tests/infrastructure/`

tests for repo infrastructure behavior.
Files (verified 2026-08-29): __init__.py, test_build.py, test_cli.py, test_dependency_verifier.py, test_repo_structure.py, test_test_infrastructure.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/infrastructure` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.