# AGENTS.md — `src/metainformant/gwas/finemapping`

Sub-package `finemapping` of the `gwas` domain module
(verified 2026-08-29: 5 files, 0 subdirs).
Part of the METAINFORMANT package; keep public API in sync with docs and
README. `__pycache__/` is generated.
Repo-wide policy: see the repository-root `AGENTS.md`.

## Public API

| Module | Functions |
|--------|-----------|
| `credible_sets.py` | `compute_credible_set()`, `compute_log_bayes_factors()`, `compute_bayes_factors()`, `log_sum_exp()`, `susie_regression()`, `colocalization()`, `colocalization_from_log_bayes_factors()`, `conditional_analysis()`, `annotate_credible_set()` |
| `colocalization.py` | `multi_trait_coloc()`, `eqtl_coloc()`, `compute_clpp()`, `regional_coloc()` |
| `eqtl.py` | `cis_eqtl_scan()`, `trans_eqtl_scan()`, `conditional_eqtl()`, `eqtl_effect_sizes()`, `eqtl_summary_stats()`, `load_transcriptome_variants()` |

The package `__init__.py` exports the three modules, not their functions. `metainformant.gwas` re-exports `eqtl_coloc`, `multi_trait_coloc`, `compute_clpp`, `regional_coloc`, `cis_eqtl_scan`, `trans_eqtl_scan`, `conditional_eqtl`, `eqtl_effect_sizes` and `eqtl_summary_stats`; its `compute_credible_set` is the association-result wrapper in `visualization/interactive/finemapping.py`, not `credible_sets.compute_credible_set`.

## Rules

- Every Bayes factor is a Wakefield (2009) ABF from `compute_log_bayes_factors()`. `V` and `W` must share a scale: pass `standard_errors=` (`V = SE**2`) or `sample_size=` with `maf=` and `case_fraction=` as available. Without either, the z-score scale applies (`variance_model == "z_scale"`, `V = 1`, prior SD 5.2). Never pair a beta-scale `prior_variance` with `V = 1`.
- `standard_errors` and `sample_size` are mutually exclusive. The default `W` is coloc's: `0.15**2` quantitative, `0.2**2` when `case_fraction` marks a case-control trait.
- Keep Bayes factors, PIPs and colocalization posteriors in log space and use `log_sum_exp()`; exponentiate only after normalising. The linear-scale `bayes_factors` returned by `compute_credible_set()` and `compute_bayes_factors()` are capped at `exp(500)`.
- A function that computes Bayes factors from Z-scores accepts the same scale arguments and reports the `variance_model` and `prior_variance` it used (`variance_model` and `prior_variance` from `compute_credible_set()`, `bayes_factor_model` from `colocalization()`).
- `colocalization_from_log_bayes_factors()` takes per-variant log Bayes factors from any source; `colocalization()` computes them from Z-scores first. The `i != j` sum behind `H3` is built from leave-one-out log sums, not from `sum(BF1) * sum(BF2) - sum(BF1 * BF2)`, which cancels catastrophically when one variant dominates both traits.
- Dict-returning functions report invalid input as `{"status": "error", "message": ...}`. `compute_log_bayes_factors()` and `compute_bayes_factors()` raise `ValueError`.
- `multi_trait_coloc()` and `regional_coloc()` take per-trait scale inputs as mappings keyed by trait name, and a trait in no mapping is scored on the z-score scale. `multi_trait_coloc()` returns an error for a key that is not a trait of `z_scores`. `regional_coloc()` slices per-variant inputs with each window and returns an error when one does not hold one value per variant of the region.

## Tests

Real-implementation tests with independent oracles live in `tests/gwas/`: `test_gwas_finemapping_abf.py` (Bayes factors, credible sets, colocalization), `test_gwas_finemapping_statistics.py` (PIPs, SuSiE, conditional analysis), `test_gwas_interactive_credible_set.py` and `test_gwas_visualization_finemapping.py` (plotting wrappers).
