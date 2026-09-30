# tables tests

Fixture data for the CSTMM test pipeline stage `tables` (per-species TSVs and/or PDFs consumed by CSTMM tests).

## Files

- `speciesA.metadata.tsv`
- `speciesA.no.tc.tsv`
- `speciesA.uncorrected.tc.tsv`

Run from the repo root: `pytest tests/cstmm_fixture/curate/speciesA/tables/ -v`. Tests follow the real-implementation policy (no mocks).
