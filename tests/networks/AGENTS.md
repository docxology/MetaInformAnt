# AGENTS.md — `MetaInformAnt/tests/networks`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `networks` domain of METAINFORMANT. Tests import from `src/metainformant/networks` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_networks_community.py`
- `test_networks_comprehensive.py`
- `test_networks_graph.py`
- `test_networks_pathway.py`
- `test_networks_ppi.py`
- `test_networks_regulatory.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/networks/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
