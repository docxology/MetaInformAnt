# metabolomics tests

pytest suite for the `metabolomics` domain of METAINFORMANT. Tests import from `src/metainformant/metabolomics` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_metabolomics_identification.py`
- `test_metabolomics_io.py`
- `test_metabolomics_pathways.py`
- `test_metabolomics_visualization.py`

Run from the repo root: `pytest tests/metabolomics/ -v`. Tests follow the real-implementation policy (no mocks).
