# AGENTS.md — `MetaInformAnt/tests/structural_variants/`

tests for the structural-variants module.
Files (verified 2026-09-05): __init__.py, test_structural_variants.py, test_structural_variants_detection_advanced.py, test_structural_variants_population.py, test_quality_filter_depth.py, test_functional_impact_depth.py, test_merge_depth.py, test_detection_depth.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/structural_variants -p no:cacheprovider` (verified 2026-09-05: 226 passed).
Repo-wide policy: see the repository-root `AGENTS.md`.