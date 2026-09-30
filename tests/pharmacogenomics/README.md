# pharmacogenomics tests

pytest suite for the `pharmacogenomics` domain of METAINFORMANT. Tests import from `src/metainformant/pharmacogenomics` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_pharmacogenomics_alleles.py`
- `test_pharmacogenomics_annotations.py`
- `test_pharmacogenomics_clinical.py`
- `test_pharmacogenomics_metabolism.py`
- `test_pharmacogenomics_visualization.py`

Run from the repo root: `pytest tests/pharmacogenomics/ -v`. Tests follow the real-implementation policy (no mocks).
