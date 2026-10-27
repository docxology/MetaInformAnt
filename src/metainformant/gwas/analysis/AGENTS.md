# AGENTS.md — `src/metainformant/gwas/analysis`

Sub-package `analysis` of the `gwas` domain module
(verified 2026-08-29: 25 files, 0 subdirs).
Part of the METAINFORMANT package; keep public API in sync with docs and
README. `__pycache__/` is generated.
Repo-wide policy: see the repository-root `AGENTS.md`.

## Rules

- A test that cannot be fitted returns NaN for its statistics, with `status`, `error` and `converged=False`; it never returns p = 1 or an effect of 0. Do not raise for it, and do not fill the NaN.
- Anything that counts or corrects tests takes only the p-values that are present: use `correction.as_p_value()` to decide, report the rest as `n_excluded` (corrections) or `n_tests_failed` (summaries), and keep outputs aligned with the input.
- Writers print missing numeric values as `NA`.
- `apply_qc_filters()` raises `ValueError` when the genotype rows differ in length or the number of variant records differs from the number of genotype columns.
- Results that reach the writers carry `maf` (`summary_stats.minor_allele_frequency()` computes it); an existing finite `maf` or `MAF` wins.
- Code that samples takes a `seed` and draws from `random.Random(seed)`, never from the global state (`calculate_ld_decay`, `compute_ld_decay`, `ld_decay_plot`, `subsample_convergence`, `ldsr_plot`).
- `run_eqtl_analysis()` scans all pairs unless both `variant_positions` and `gene_positions` are given; with them it keeps `|position - TSS| <= cis_window` on the same chromosome and adds `distance`. A `cis_window` without positions raises.
