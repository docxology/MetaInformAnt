# gwas tests

pytest suite for the `gwas` domain of METAINFORMANT. Tests import from `src/metainformant/gwas` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_gwas_annotation.py`
- `test_gwas_association.py`
- `test_gwas_benchmarking.py`
- `test_gwas_calling.py`
- `test_gwas_config.py`
- `test_gwas_config_pbarbatus.py`
- `test_gwas_contact_policy.py`
- `test_gwas_correction.py`
- `test_gwas_download.py`
- `test_gwas_download_contracts.py`
- `test_gwas_end_to_end.py`

Run from the repo root: `pytest tests/gwas/ -v`. Tests follow the real-implementation policy (no mocks).
