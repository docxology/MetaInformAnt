# AGENTS.md — `MetaInformAnt/tests/cstmm_fixture/`

shared CSTMM test fixtures.
Files (verified 2026-08-29): (none — fixture/support dir).


## Children (documented on disk)
- csca/, cstmm/, curate/, merge/, metadata/, orthogroups/

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/cstmm_fixture` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.