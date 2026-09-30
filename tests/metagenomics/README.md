# metagenomics tests

pytest suite for the `metagenomics` domain of METAINFORMANT. Tests import from `src/metainformant/metagenomics` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_metagenomics_amplicon.py`
- `test_metagenomics_comparative.py`
- `test_metagenomics_diversity.py`
- `test_metagenomics_functional.py`
- `test_metagenomics_shotgun.py`
- `test_metagenomics_visualization.py`

Run from the repo root: `pytest tests/metagenomics/ -v`. Tests follow the real-implementation policy (no mocks).
