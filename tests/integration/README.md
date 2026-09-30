# integration tests

Cross-domain integration tests: eQTL integration/scripts, RNA→GWAS handoff, comprehensive integration.

## Files

- `__init__.py`
- `test_eqtl_integration.py`
- `test_eqtl_scripts.py`
- `test_integration_comprehensive.py`
- `test_rna_gwas_handoff.py`

Run from the repo root: `pytest tests/integration/ -v`. Tests follow the real-implementation policy (no mocks).
