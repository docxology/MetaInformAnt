# cstmm tests

Fixture data for the CSTMM test pipeline stage `cstmm` (per-species TSVs and/or PDFs consumed by CSTMM tests).

## Files

- `cstmm_exclusion.pdf`
- `cstmm_mean_expression_boxplot.pdf`
- `cstmm_orthogroup_genecount.tsv`
- `metadata.tsv`

Subdirectories: `speciesA`, `speciesB`.

Run from the repo root: `pytest tests/cstmm_fixture/cstmm/ -v`. Tests follow the real-implementation policy (no mocks).
