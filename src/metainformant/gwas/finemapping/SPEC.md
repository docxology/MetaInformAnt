# Specification: finemapping

## 🎯 Scope
GWAS statistical fine-mapping submodule: Wakefield approximate Bayes factors, credible sets, SuSiE regression, coloc-style colocalization, conditional analysis, and eQTL scanning.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 4 Python modules
- **Key Concepts**: Wakefield approximate Bayes factors accumulated in log space, posterior inclusion probabilities (PIPs), credible sets, and the five coloc hypotheses H0 to H4. Most results are plain dicts with a `status` key.

## 🔌 API Definition
### Exports
- `__init__.py`
- `colocalization.py`
- `credible_sets.py`
- `eqtl.py`

### Bayes-factor model

`log ABF = (log(V / (V + W)) + z^2 * W / (V + W)) / 2` per variant (Wakefield 2009), where `z` is the Z-score, `V` the sampling variance of the effect estimate and `W` the prior variance of the true effect. Exactly one source fixes `V`:

| Source | `V` | `variance_model` |
|--------|-----|------------------|
| `standard_errors` (scalar or per variant) | `SE**2` | `"standard_errors"` |
| `sample_size` (scalar or per variant), optional `maf`, optional `case_fraction` `s` | `1 / (N * 2f(1 - f) * s(1 - s))` for a phenotype of unit variance, where the `2f(1 - f)` factor is 1 without `maf` and the `s(1 - s)` factor is 1 without `case_fraction` | `"sample_size"` |
| neither | 1 | `"z_scale"` |

`prior_variance` is `W`. It defaults to `DEFAULT_PRIOR_VARIANCE_QUANTITATIVE` (`0.15**2`), to `DEFAULT_PRIOR_VARIANCE_CASE_CONTROL` (`0.2**2`) when `case_fraction` is given together with `standard_errors` or `sample_size`, and to `DEFAULT_PRIOR_VARIANCE_Z_SCALE` (`5.2**2`) on the z-score scale. The module constants `DEFAULT_PRIOR_SD_QUANTITATIVE`, `DEFAULT_PRIOR_SD_CASE_CONTROL`, `DEFAULT_PRIOR_SD_Z_SCALE`, `VARIANCE_MODEL_STANDARD_ERRORS`, `VARIANCE_MODEL_SAMPLE_SIZE` and `VARIANCE_MODEL_Z_SCALE` are public.

`standard_errors` and `sample_size` are mutually exclusive. `case_fraction` must lie in (0, 1), `maf` (used with `sample_size`) must lie in (0, 1), `prior_variance`, `standard_errors` and `sample_size` must be positive, and every input must be finite.

### credible_sets

| Function | Signature | Returns |
|----------|-----------|---------|
| `compute_credible_set` | `(z_scores, ld_matrix=None, coverage=0.95, prior="uniform", *, prior_variance=None, standard_errors=None, sample_size=None, maf=None, case_fraction=None)` | `status`, `snps_in_set`, `pips`, `coverage_achieved`, `n_snps`, `bayes_factors`, `log_bayes_factors`, `variance_model`, `prior_variance` |
| `compute_log_bayes_factors` | `(z_scores, prior_variance=None, *, standard_errors=None, sample_size=None, maf=None, case_fraction=None)` | `list[float]` of natural-log ABFs; raises `ValueError` |
| `compute_bayes_factors` | same arguments | `list[float]` of ABFs capped at `exp(500)`; raises `ValueError` |
| `log_sum_exp` | `(values)` | `float`; `-inf` for empty input or all `-inf` |
| `susie_regression` | `(X, y, L=10, max_iter=100, tol=1e-3)` | `status`, `alpha`, `mu`, `sigma2`, `pip`, `credible_sets`, `converged`, `elbo` |
| `colocalization` | `(z_scores_1, z_scores_2, ld_matrix=None, prior_p1=1e-4, prior_p2=1e-4, prior_p12=1e-5, *, prior_variance_1=None, prior_variance_2=None, standard_errors_1=None, standard_errors_2=None, sample_size_1=None, sample_size_2=None, maf=None, case_fraction_1=None, case_fraction_2=None)` | `status`, `PP_H0` to `PP_H4`, `summary`, `snp_pp_h4`, `bayes_factor_model` |
| `colocalization_from_log_bayes_factors` | `(log_bayes_factors_1, log_bayes_factors_2, prior_p1=1e-4, prior_p2=1e-4, prior_p12=1e-5)` | `status`, `PP_H0` to `PP_H4`, `summary`, `snp_pp_h4` |
| `conditional_analysis` | `(z_scores, ld_matrix, max_signals=10, p_threshold=5e-8)` | `list[dict]` with `index`, `z_score`, `z_conditional`, `p_value`, `signal_number` |
| `annotate_credible_set` | `(credible_set, annotations=None)` | `status`, `credible_set`, `annotated_snps`, `enrichment`, `n_annotated` |

`compute_credible_set` accepts `prior` of `"uniform"` or `"distance"` and `coverage` in (0, 1]. Its `ld_matrix` holds LD correlations r (signed entries in [-1, 1], ones on the diagonal, square with one row per Z-score), not r-squared: the adjustment multiplies each variant's Bayes factor by `1 - r**2`, where r is its largest absolute correlation with a variant that has a larger Bayes factor (the factor never drops below 0.01). `snp_pp_h4` holds each variant's probability of being the shared causal variant given H4 and sums to 1. `bayes_factor_model` maps `"trait_1"` and `"trait_2"` to `{"variance_model", "prior_variance"}`.

### colocalization

| Function | Signature | Returns |
|----------|-----------|---------|
| `multi_trait_coloc` | `(z_scores, ld_matrix=None, *, prior_variances=None, standard_errors=None, sample_sizes=None, maf=None, case_fractions=None)` | `status`, `pairwise`, `multi_trait_pp_shared`, `n_traits`, `n_variants`, `summary` |
| `eqtl_coloc` | `(gwas_z, eqtl_z, gene_id, ld_matrix=None, *, gwas_prior_variance=None, eqtl_prior_variance=None, gwas_standard_errors=None, eqtl_standard_errors=None, gwas_sample_size=None, eqtl_sample_size=None, maf=None, gwas_case_fraction=None)` | `status`, `gene_id`, `coloc_result`, `pp_h4`, `is_colocalized`, `lead_variant_gwas`, `lead_variant_eqtl`, `shared_variant`, `interpretation` |
| `compute_clpp` | `(pip_1, pip_2)` | `status`, `clpp`, `per_variant_clpp`, `max_clpp_index`, `max_clpp_value`, `is_colocalized`, `n_variants` |
| `regional_coloc` | `(z_scores_dict, region, ld_matrix=None, *, prior_variances=None, standard_errors=None, sample_sizes=None, maf=None, case_fractions=None)` | `status`, `best_window`, `all_windows`, `n_windows`, `n_failed_windows`, `region` |

The per-trait arguments of `multi_trait_coloc` and `regional_coloc` are mappings keyed by trait name; `maf` is shared by all traits. `eqtl_coloc` treats the eQTL trait as quantitative.

`regional_coloc` picks `best_window` among the windows that succeeded and counts the others in `n_failed_windows`; a failed window stays in `all_windows` with its own `status` and `message`. When every window fails, the result is `{"status": "error", "message": ..., "all_windows": ..., "n_windows": ..., "n_failed_windows": ...}` with no `best_window`.

### eqtl

| Function | Signature | Returns |
|----------|-----------|---------|
| `cis_eqtl_scan` | `(expression_matrix, genotype_matrix, gene_positions, variant_positions, cis_window=1_000_000, maf_threshold=0.05)` | DataFrame: `gene_id`, `variant_id`, `beta`, `se`, `pvalue`, `distance`, `maf`, `status`, `error`; raises `ValueError` for a repeated id in `gene_positions`, `variant_positions` or the genotype index |
| `trans_eqtl_scan` | `(expression_matrix, genotype_matrix, gene_positions, variant_positions, cis_window=1_000_000, maf_threshold=0.05, pvalue_threshold=1e-5)` | DataFrame of the variants with `pvalue` below the threshold: `gene_id`, `variant_id`, `beta`, `se`, `pvalue`, `is_trans`, `maf` |
| `conditional_eqtl` | `(expression_vector, genotype_matrix, lead_variant_id, test_variant_ids=None)` | DataFrame: `variant_id`, `conditioned_on`, `beta`, `se`, `pvalue`, `status`, `error` |
| `eqtl_effect_sizes` | `(expression_matrix, genotype_matrix, associations)` | DataFrame: `gene_id`, `variant_id`, `beta`, `se`, `pvalue`, `r_squared`, `status`, `error`, `mean_expr_gt0` to `mean_expr_gt2` where a genotype class is present |
| `eqtl_summary_stats` | `(cis_results, fdr_threshold=0.05)` | `n_tests`, `n_tests_failed`, `n_egenes`, `n_eqtls` and, when there are tests, `fdr_threshold`, `mean_effect_size` |
| `write_eqtl_results` | `(results, output_path)` | `Path` of the tab-separated file; every missing value is written as `NA` |
| `load_transcriptome_variants` | `(vcf_path)` | `(genotype_matrix, variant_positions)` from a VCF or VCF.gz |

A fit that cannot be made (fewer than 3 complete observations, a constant genotype, or a zero standard error) has NaN `beta`, `se` and `pvalue` (and NaN `r_squared` in `eqtl_effect_sizes`), `status == "error"` and the reason in `error`, in `cis_eqtl_scan`, `conditional_eqtl` and `eqtl_effect_sizes`. NaN is missing, never p = 1. `eqtl_summary_stats` counts in `n_tests` the rows that have a p-value, runs Benjamini-Hochberg over them and reports the rest as `n_tests_failed`. `trans_eqtl_scan` never reports a failed fit, because NaN is not below the threshold.

### Error contract

Functions that return a dict report invalid input as `{"status": "error", "message": ...}`. `compute_log_bayes_factors` and `compute_bayes_factors` raise `ValueError`.
