# AGENTS.md — `MetaInformAnt/tests/scripts/`

tests exercising `scripts/` workflow tooling.
Files (verified 2026-08-29): test_verify_utils.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/scripts` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.