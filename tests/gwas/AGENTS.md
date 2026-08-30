# AGENTS.md — `MetaInformAnt/tests/gwas/`

tests for the GWAS module.
Files (verified 2026-08-29): __init__.py, test_gwas_annotation.py, test_gwas_association.py, test_gwas_benchmarking.py, test_gwas_calling.py, test_gwas_config.py, test_gwas_config_pbarbatus.py, test_gwas_contact_policy.py … (49 files).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/gwas` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.