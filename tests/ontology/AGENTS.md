# AGENTS.md — `MetaInformAnt/tests/ontology`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `ontology` domain of METAINFORMANT. Tests import from `src/metainformant/ontology` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_ontology_api.py`
- `test_ontology_comprehensive.py`
- `test_ontology_enrichment.py`
- `test_ontology_go_basic.py`
- `test_ontology_obo_parser.py`
- `test_ontology_query.py`
- `test_ontology_serialization.py`
- `test_ontology_serialize.py`
- `test_ontology_types.py`
- `test_ontology_visualization.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/ontology/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
