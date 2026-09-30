# other tests

Cross-cutting tests that don't belong to one domain: domain-module imports, examples validation, import verification, orchestrators, ortholog generation, tissue config.

## Files

- `__init__.py`
- `test_domain_modules.py`
- `test_examples.py`
- `test_import_verification.py`
- `test_orchestrators.py`
- `test_ortholog_generation.py`
- `test_tissue_config.py`

Run from the repo root: `pytest tests/other/ -v`. Tests follow the real-implementation policy (no mocks).
