# Specification: analysis

## 🎯 Scope
GWAS analysis modules.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 23 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `alignment.py`
- `annotation.py`
- `association.py`
- `benchmarking.py`
- `calling.py`
- `correction.py`
- `enrichment.py`
- `eqtl.py`
- `gene_annotation_api.py`
- ...

## Behavioral Contracts

The contracts below are stated in full in [README.md](README.md) and in the gwas [SPEC.md](../SPEC.md).

- **Failed tests** are returned, not raised: NaN for `beta`, `se`, `t_stat` (`z_stat` for the logistic test), `p_value` and `r_squared` (`odds_ratio` for the logistic test), plus `status`, `error` and `converged: False`. A linear test with a constant phenotype fails with an error message; an exact fit with a nonzero genotype effect reports t of plus or minus infinity and p = 0; zero residual variance with a zero effect fails. Mixed-model rows carry `status`.
- **Corrections** (`fdr_correction`, `bonferroni_correction`, `qvalue_estimation`, `adjust_p_values`) leave NaN or `None` p-values out of the family, keep their outputs aligned with the input, and report `n_excluded` in the dictionaries of `fdr_correction` and `bonferroni_correction`; `as_p_value` decides what is missing.
- **Summaries**: `create_results_summary` counts only the tests that ran in `n_variants_tested` and reports `n_tests_failed`; `compute_comprehensive_summary` also reports `n_tests_failed`; `per_chromosome_summary`, `sign_test`, `subsample_convergence` and `jackknife_se` skip failed tests; the writers print `NA`.
- **QC**: `apply_qc_filters` raises `ValueError` when the genotype rows differ in length or the number of variant records differs from the number of genotype columns.
- **MAF**: `summary_stats.minor_allele_frequency` is public; results carry `maf`, and the workflows write it to `summary_statistics.tsv`.
- **Randomness**: `calculate_ld_decay`, `compute_ld_decay`, `ld_decay_plot`, `subsample_convergence` and `ldsr_plot` take `seed` (default 42) and use their own generator.
- **eQTL**: `run_eqtl_analysis` takes `variant_positions` and `gene_positions`, restricts the scan to variants on the transcript's chromosome with `|position - TSS| <= cis_window`, and adds a `distance` column; a `cis_window` without positions raises `ValueError`.
