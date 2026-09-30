# AGENTS.md — `MetaInformAnt/tests/gwas`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `gwas` domain of METAINFORMANT. Tests import from `src/metainformant/gwas` and follow the repo's real-implementation policy.

## Layout

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
- `test_gwas_genome.py`
- `test_gwas_heritability.py`
- … (+35 more test modules)

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/gwas/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
