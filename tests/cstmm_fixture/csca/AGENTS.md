# AGENTS.md — `MetaInformAnt/tests/cstmm_fixture/csca`

Verified against disk 2026-08-30 (doc-realization fleet pass). Fixture data for the CSTMM test pipeline stage `csca` (per-species TSVs and/or PDFs consumed by CSTMM tests).

## Layout

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
- `csca_ortholog_unaveraged.corrected.tsv`
- `csca_ortholog_unaveraged.imputed.corrected.tsv`
- … (+7 more test modules)

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/cstmm_fixture/csca/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
- Data is tiny, deterministic fixture data — safe to inspect, never edit by hand.
