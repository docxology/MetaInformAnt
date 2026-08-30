# AGENTS.md — `MetaInformAnt/tests/fixtures/`

shared pytest fixtures.
Files (verified 2026-08-29): __init__.py.


## Children (documented on disk)
- gwas/

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/fixtures` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.