# AGENTS.md — `MetaInformAnt/tests/data/ontology/GO_v3`

Verified against disk 2026-08-30 (doc-realization fleet pass). Gene-Ontology pipeline fixture: numbered step scripts (1-4) that build a gene-to-GO annotation summary from UniProt IDs, plus a `.keep` marker.

## Layout

- `.keep`
- `1_uniprot_ID_extract.py`
- `2_genetogo.py`
- `3_genetogotoanno.py`
- `4_genetogo_summary.py`
- `PAI.md`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/data/ontology/GO_v3/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
- Scripts are numbered pipeline steps (1_uniprot_ID_extract → 4_genetogo_summary); run in order.
