# csca tests

Fixture data for the CSTMM test pipeline stage `csca` (per-species TSVs and/or PDFs consumed by CSTMM tests).

## Files

- `csca_SVA_dendrogram.pdf`
- `csca_SVA_heatmap.pdf`
- `csca_averaged_summary.pdf`
- `csca_boxplot.pdf`
- `csca_color_averaged.tsv`
- `csca_color_unaveraged.tsv`
- `csca_exclusion.pdf`
- `csca_group_cor_scatter.pdf`
- `csca_ortholog_averaged.corrected.tsv`
- `csca_ortholog_averaged.imputed.corrected.tsv`
- `csca_ortholog_averaged.imputed.uncorrected.tsv`
- `csca_ortholog_averaged.uncorrected.tsv`

Run from the repo root: `pytest tests/cstmm_fixture/csca/ -v`. Tests follow the real-implementation policy (no mocks).
