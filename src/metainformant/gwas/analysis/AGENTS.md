# AGENTS.md — `MetaInformAnt/src/metainformant/gwas/analysis`

Source module under the METAINFORMANT bioinformatics toolkit (`gwas` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `alignment.py` — Read alignment and BAM processing for GWAS analysis.
- `annotation.py` — SNP-to-gene annotation for GWAS results.
- `association.py` — GWAS association testing utilities.
- `benchmarking.py` — GWAS compute-time benchmarking and runtime extrapolation.
- `calling.py` — Variant calling integration for GWAS analysis.
- `correction.py` — GWAS multiple testing correction utilities.
- `enrichment.py` — Gene-set enrichment analysis for GWAS hits.
- `eqtl.py` — eQTL Analysis Wrapper.
- `gene_annotation_api.py` — Gene region annotation using live NCBI E-utilities for Apis mellifera.
- `heritability.py` — SNP heritability estimation and partitioning.
- `hwe.py` — Hardy-Weinberg Equilibrium (HWE) analysis for GWAS QC.
- `ld_decay.py` — LD decay analysis for population genomics and GWAS QC.
- `ld_pruning.py` — Linkage disequilibrium (LD) pruning for GWAS.
- `ldsr.py` — LD Score Regression (LDSR) — SNP-heritability and intercept estimation.
- `mixed_model.py` — Mixed linear model (MLM) for GWAS using the EMMA algorithm.
- `power.py` — GWAS statistical power estimation, convergence, and saturation analysis.
- `prs.py` — Polygenic Risk Score (PRS) construction and validation.
- `quality.py` — GWAS quality control and data filtering utilities.
- `strain_analysis.py` — Strain-specific variant analysis for population-structured GWAS.
- `structure.py` — GWAS population structure analysis utilities.
- `summary_stats.py` — GWAS summary statistics output utilities.
- `utils.py` — Shared utility functions for GWAS analysis.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-gwas-analysis` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
