# AGENTS.md — `MetaInformAnt/tests/pharmacogenomics/`

tests for the pharmacogenomics module.
Files (verified 2026-09-05): __init__.py, test_pharmacogenomics_alleles.py, test_pharmacogenomics_annotations.py, test_pharmacogenomics_clinical.py, test_pharmacogenomics_metabolism.py, test_pharmacogenomics_visualization.py, test_drug_interactions_depth.py, test_reporting_depth.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/pharmacogenomics`
Repo-wide policy: see the repository-root `AGENTS.md`.