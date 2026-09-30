# AGENTS.md — `MetaInformAnt/tests/cstmm_fixture/cstmm/speciesB`

Verified against disk 2026-08-30 (doc-realization fleet pass). Fixture data for the CSTMM test pipeline stage `speciesB` (per-species TSVs and/or PDFs consumed by CSTMM tests).

## Layout

- `speciesB_cstmm_counts.tsv`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/cstmm_fixture/cstmm/speciesB/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
- Data is tiny, deterministic fixture data — safe to inspect, never edit by hand.
