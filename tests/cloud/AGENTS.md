# AGENTS.md — `MetaInformAnt/tests/cloud/`

tests for the cloud/GCP worker scripts (mirror of `src/metainformant/cloud/`).
Files (verified 2026-08-29): __init__.py, test_cloud.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/cloud` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.