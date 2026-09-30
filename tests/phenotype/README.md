# phenotype tests

pytest suite for the `phenotype` domain of METAINFORMANT. Tests import from `src/metainformant/phenotype` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_multivariate.py`
- `test_phenotype_basic.py`
- `test_phenotype_behavior.py`
- `test_phenotype_comprehensive.py`
- `test_phenotype_integration.py`
- `test_phenotype_life_course.py`
- `test_phenotype_modules.py`
- `test_phenotype_morphological.py`
- `test_phenotype_scraper.py`
- `test_phenotype_workflow.py`
- `test_plots.py`

Run from the repo root: `pytest tests/phenotype/ -v`. Tests follow the real-implementation policy (no mocks).
