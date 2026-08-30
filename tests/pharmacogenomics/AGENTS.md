# AGENTS.md — `MetaInformAnt/tests/pharmacogenomics/`

tests for the pharmacogenomics module.
Files (verified 2026-08-29): __init__.py, test_pharmacogenomics_alleles.py, test_pharmacogenomics_annotations.py, test_pharmacogenomics_clinical.py, test_pharmacogenomics_metabolism.py, test_pharmacogenomics_visualization.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/pharmacogenomics` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.