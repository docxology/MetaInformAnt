# merge tests

Fixture data for the CSTMM test pipeline stage `merge` (per-species TSVs and/or PDFs consumed by CSTMM tests).

## Files

- `metadata.tsv`

Subdirectories: `speciesA`, `speciesB`.

Run from the repo root: `pytest tests/cstmm_fixture/merge/ -v`. Tests follow the real-implementation policy (no mocks).
