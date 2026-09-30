# tables tests

Fixture data for the CSTMM test pipeline stage `tables` (per-species TSVs and/or PDFs consumed by CSTMM tests).

## Files

- `speciesB.metadata.tsv`
- `speciesB.no.tc.tsv`
- `speciesB.uncorrected.tc.tsv`

Run from the repo root: `pytest tests/cstmm_fixture/curate/speciesB/tables/ -v`. Tests follow the real-implementation policy (no mocks).
