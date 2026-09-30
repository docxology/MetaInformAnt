# AGENTS.md — `MetaInformAnt/tests/visualization`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `visualization` domain of METAINFORMANT. Tests import from `src/metainformant/visualization` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_cross_species.py`
- `test_visualization.py`
- `test_visualization_animations.py`
- `test_visualization_basic.py`
- `test_visualization_comprehensive.py`
- `test_visualization_dimred.py`
- `test_visualization_expression.py`
- `test_visualization_genomics.py`
- `test_visualization_information.py`
- `test_visualization_multidim.py`
- `test_visualization_networks.py`
- `test_visualization_phylo.py`
- `test_visualization_quality.py`
- … (+3 more test modules)

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/visualization/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
