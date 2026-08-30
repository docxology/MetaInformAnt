# AGENTS.md — `MetaInformAnt/tests/other/`

miscellaneous tests not tied to a single module.
Files (verified 2026-08-29): __init__.py, test_domain_modules.py, test_examples.py, test_import_verification.py, test_orchestrators.py, test_ortholog_generation.py, test_tissue_config.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/other` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.