# spatial tests

pytest suite for the `spatial` domain of METAINFORMANT. Tests import from `src/metainformant/spatial` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_spatial_analysis.py`
- `test_spatial_autocorrelation.py`
- `test_spatial_communication.py`
- `test_spatial_deconvolution_advanced.py`
- `test_spatial_integration.py`
- `test_spatial_io.py`
- `test_spatial_neighborhood.py`
- `test_spatial_visualization.py`

Run from the repo root: `pytest tests/spatial/ -v`. Tests follow the real-implementation policy (no mocks).
