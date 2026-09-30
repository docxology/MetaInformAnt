# multiomics tests

pytest suite for the `multiomics` domain of METAINFORMANT. Tests import from `src/metainformant/multiomics` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_multiomics_comprehensive.py`
- `test_multiomics_integration.py`
- `test_multiomics_methods.py`
- `test_multiomics_pathways.py`
- `test_multiomics_sample_mapping.py`
- `test_multiomics_survival.py`

Run from the repo root: `pytest tests/multiomics/ -v`. Tests follow the real-implementation policy (no mocks).
