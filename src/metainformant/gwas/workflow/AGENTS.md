# AGENTS.md — `src/metainformant/gwas/workflow`

Sub-package `workflow` of the `gwas` domain module
(verified 2026-08-29: 7 files, 0 subdirs).
Part of the METAINFORMANT package; keep public API in sync with docs and
README. `__pycache__/` is generated.
Repo-wide policy: see the repository-root `AGENTS.md`.

## Rules

- Pass `layout="sample_major"` to `_extract_genotype_matrix()` for anything that came from `parse_vcf_full()`; never let the shape decide when the producer is known.
- Keep failed tests as rows with NaN statistics, take plots, corrections, the top hit and fine-mapping from the rows with a p-value only, and report the count as `n_tests_failed`.
- Attach `maf` to results with `_attach_maf()` and pass the genotypes to the writers; never fill the MAF column or panel from another quantity (such as an effect size).
