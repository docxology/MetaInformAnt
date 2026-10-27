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

`compute_credible_set` accepts `prior` of `"uniform"` or `"distance"` and `coverage` in (0, 1]. `snp_pp_h4` holds each variant's probability of being the shared causal variant given H4 and sums to 1. `bayes_factor_model` maps `"trait_1"` and `"trait_2"` to `{"variance_model", "prior_variance"}`.

### colocalization

| Function | Signature | Returns |
|----------|-----------|---------|
| `multi_trait_coloc` | `(z_scores, ld_matrix=None, *, prior_variances=None, standard_errors=None, sample_sizes=None, maf=None, case_fractions=None)` | `status`, `pairwise`, `multi_trait_pp_shared`, `n_traits`, `n_variants`, `summary` |
| `eqtl_coloc` | `(gwas_z, eqtl_z, gene_id, ld_matrix=None, *, gwas_prior_variance=None, eqtl_prior_variance=None, gwas_standard_errors=None, eqtl_standard_errors=None, gwas_sample_size=None, eqtl_sample_size=None, maf=None, gwas_case_fraction=None)` | `status`, `gene_id`, `coloc_result`, `pp_h4`, `is_colocalized`, `lead_variant_gwas`, `lead_variant_eqtl`, `shared_variant`, `interpretation` |
| `compute_clpp` | `(pip_1, pip_2)` | `status`, `clpp`, `per_variant_clpp`, `max_clpp_index`, `max_clpp_value`, `is_colocalized`, `n_variants` |
| `regional_coloc` | `(z_scores_dict, region, ld_matrix=None, *, prior_variances=None, standard_errors=None, sample_sizes=None, maf=None, case_fractions=None)` | `status`, `best_window`, `all_windows`, `n_windows`, `region` |

The per-trait arguments of `multi_trait_coloc` and `regional_coloc` are mappings keyed by trait name; `maf` is shared by all traits. `eqtl_coloc` treats the eQTL trait as quantitative.

### eqtl

`cis_eqtl_scan`, `trans_eqtl_scan`, `conditional_eqtl`, `eqtl_effect_sizes`, `eqtl_summary_stats` and `load_transcriptome_variants`.

### Error contract

Functions that return a dict report invalid input as `{"status": "error", "message": ...}`. `compute_log_bayes_factors` and `compute_bayes_factors` raise `ValueError`.
