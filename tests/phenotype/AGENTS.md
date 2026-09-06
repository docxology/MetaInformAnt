# AGENTS.md — `MetaInformAnt/tests/phenotype/`

tests for the phenotype module.
Files (verified 2026-09-05): __init__.py, test_antwiki_records.py, test_mappings.py, test_multivariate.py, test_phenotype_basic.py, test_phenotype_behavior.py, test_phenotype_comprehensive.py, test_phenotype_integration.py, test_phenotype_life_course.py, test_phenotype_modules.py, test_phenotype_morphological.py, test_phenotype_scraper.py, test_phenotype_workflow.py, test_phewas_regression.py, test_plots.py, test_statistical.py (16 files).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/phenotype -p no:cacheprovider` (verified 2026-09-05).
Repo-wide policy: see the repository-root `AGENTS.md`.