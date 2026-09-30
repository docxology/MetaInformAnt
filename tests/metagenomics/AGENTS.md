# AGENTS.md — `MetaInformAnt/tests/metagenomics`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `metagenomics` domain of METAINFORMANT. Tests import from `src/metainformant/metagenomics` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_metagenomics_amplicon.py`
- `test_metagenomics_comparative.py`
- `test_metagenomics_diversity.py`
- `test_metagenomics_functional.py`
- `test_metagenomics_shotgun.py`
- `test_metagenomics_visualization.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/metagenomics/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
