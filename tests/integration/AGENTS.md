# AGENTS.md — `MetaInformAnt/tests/integration/`

cross-module integration tests.
Files (verified 2026-08-29): __init__.py, test_eqtl_integration.py, test_eqtl_scripts.py, test_integration_comprehensive.py, test_rna_gwas_handoff.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/integration` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.