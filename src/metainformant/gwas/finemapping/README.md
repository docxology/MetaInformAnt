# GWAS Fine-Mapping

Statistical fine-mapping methods for identifying causal variants from GWAS signals, including credible set analysis, SuSiE regression, colocalization, and eQTL integration.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `credible_sets`, `colocalization`, `eqtl` |
| `credible_sets.py` | Wakefield ABF credible sets, log Bayes factors, SuSiE regression, pairwise colocalization, conditional analysis |
| `colocalization.py` | Multi-trait colocalization, GWAS-eQTL coloc, CLPP, regional coloc |
| `eqtl.py` | cis/trans eQTL scanning, conditional eQTL, effect size estimation, summary statistics, result writer, transcriptome VCF loading |

## Key Functions

| Function | Description |
|----------|-------------|
| `compute_credible_set()` | Wakefield ABF credible set from Z-scores with PIP calculation |
| `compute_log_bayes_factors()` | Natural-log Wakefield approximate Bayes factors, one per variant |
| `compute_bayes_factors()` | The same Bayes factors on the linear scale, capped at `exp(500)` |
| `log_sum_exp()` | `log(sum(exp(v)))` without overflow or underflow |
| `susie_regression()` | Sum of Single Effects regression for multi-causal fine-mapping |
| `colocalization()` | Pairwise trait colocalization (coloc-style) from Z-scores |
| `colocalization_from_log_bayes_factors()` | Coloc posteriors from per-variant log Bayes factors (`coloc.bf_bf`) |
| `conditional_analysis()` | Stepwise conditional association analysis |
| `annotate_credible_set()` | Functional annotation enrichment of credible sets |
| `multi_trait_coloc()` | Multi-trait colocalization extending pairwise to N traits |
| `eqtl_coloc()` | GWAS-eQTL colocalization with gene-level summaries |
| `compute_clpp()` | Colocalization posterior probability from individual PIPs |
| `regional_coloc()` | Colocalization across a genomic region with automatic window selection; reports `n_failed_windows` |
| `cis_eqtl_scan()` | Scan for cis-eQTLs within a window around each gene; a failed fit has NaN statistics, `status` and `error` |
| `trans_eqtl_scan()` | Genome-wide trans-eQTL association testing |
| `conditional_eqtl()` | Test the remaining variants after conditioning on a lead variant |
| `eqtl_effect_sizes()` | Effect size, standard error, `r_squared` and mean expression per genotype class for given gene-variant pairs |
| `eqtl_summary_stats()` | Test counts, eGenes and eQTLs from a `cis_eqtl_scan()` table (Benjamini-Hochberg over the tests with a p-value) |
| `write_eqtl_results()` | Write an eQTL table as tab-separated text with `NA` for missing values |
| `load_transcriptome_variants()` | Read a (multi-sample) VCF into the genotype matrix and variant positions the scans take |

## Approximate Bayes Factors

Credible sets and colocalization are built on Wakefield (2009) approximate Bayes factors. The effect estimate is modelled as `beta_hat ~ N(beta, V)` and the true effect as `beta ~ N(0, W)`, which gives

```
ABF = sqrt(V / (V + W)) * exp(z^2 * W / (2 * (V + W))),   z = beta_hat / sqrt(V)
```

`V` and `W` must be on the same scale, so Z-scores alone do not fix an ABF. The scale of `V` comes from exactly one of three sources:

| Input | `V` | `variance_model` |
|-------|-----|------------------|
| `standard_errors=` (scalar or one per variant) | `SE^2` | `"standard_errors"` |
| `sample_size=` (scalar or one per variant), with optional `maf=` and `case_fraction=` | Modelled as in coloc's `Var.data` and `Var.data.cc`: `1 / (2 N f (1 - f))` for a quantitative trait with unit phenotype variance and `1 / (2 N f (1 - f) s (1 - s))` for a case-control trait, with `f = maf` and `s = case_fraction`. Without `maf` the genotype is standardised and the `2 f (1 - f)` factor is 1. | `"sample_size"` |
| neither | 1, the z-score scale (a coarse fallback; prefer one of the two sources above) | `"z_scale"` |

`standard_errors` and `sample_size` are mutually exclusive.

The prior variance `W` is `prior_variance=` and defaults to coloc's values: `0.15**2` for a quantitative trait (in phenotype-SD units; pass `prior_variance=(0.15 * sd_y) ** 2` for raw units), `0.2**2` on the log-odds scale when `case_fraction` marks a case-control trait, and `5.2**2` on the z-score scale, where `W` is the prior variance of the noncentrality.

Bayes factors, PIPs and colocalization posteriors are accumulated in log space, because a correctly scaled ABF reaches `exp(z^2 / 2)` and overflows double precision for `|z|` above about 38. The linear-scale `bayes_factors` returned by `compute_credible_set()` and `compute_bayes_factors()` are capped at `exp(500)` so they stay finite; PIPs use the uncapped logs.

### Return values

| Function | Result |
|----------|--------|
| `compute_credible_set()` | Dict with `status`, `snps_in_set`, `pips`, `coverage_achieved`, `n_snps`, `bayes_factors`, `log_bayes_factors`, `variance_model`, `prior_variance` |
| `colocalization()` | Dict with `status`, `PP_H0` to `PP_H4`, `summary`, `snp_pp_h4` (per-variant probability of being the shared causal variant given H4, summing to 1), `bayes_factor_model` (per trait: `variance_model`, `prior_variance`) |
| `colocalization_from_log_bayes_factors()` | The same dict without `bayes_factor_model` |
| `multi_trait_coloc()` | Dict with `status`, `pairwise`, `multi_trait_pp_shared`, `n_traits`, `n_variants`, `summary` |
| `eqtl_coloc()` | Dict with `status`, `gene_id`, `coloc_result`, `pp_h4`, `is_colocalized`, `lead_variant_gwas`, `lead_variant_eqtl`, `shared_variant`, `interpretation` |
| `regional_coloc()` | Dict with `status`, `best_window` (the best of the windows that succeeded), `all_windows`, `n_windows`, `n_failed_windows` and `region`. When every window fails, `status` is `"error"` with a `message` and there is no `best_window` |

Functions that return a dict report invalid input as `{"status": "error", "message": ...}`. `compute_log_bayes_factors()` and `compute_bayes_factors()` raise `ValueError` for a non-positive `prior_variance`, for both `standard_errors` and `sample_size`, and for non-finite, out-of-range or wrongly sized inputs.

`compute_credible_set()` takes an optional `ld_matrix` of LD correlations r (signed entries in [-1, 1], ones on the diagonal), not r-squared: the adjustment squares it. Each variant's Bayes factor is multiplied by `1 - r**2`, where r is its largest absolute correlation with a variant that has a larger Bayes factor (the factor never drops below 0.01).

`multi_trait_coloc()` and `regional_coloc()` take the scale inputs per trait as mappings keyed by trait name (`standard_errors=`, `sample_sizes=`, `prior_variances=`, `case_fractions=`) plus a shared `maf=`; a trait that appears in no mapping is scored on the z-score scale. `eqtl_coloc()` takes `gwas_*` and `eqtl_*` keyword arguments for the same inputs. `metainformant.gwas.compute_credible_set` is a different function: the association-result wrapper in `metainformant.gwas.visualization.interactive.finemapping`.

## Usage

```python
from metainformant.gwas.finemapping import credible_sets, colocalization, eqtl
from metainformant.gwas.finemapping.credible_sets import (
    colocalization_from_log_bayes_factors,
    compute_log_bayes_factors,
)

z_scores = [0.4, -1.1, 6.8, 5.9, 0.2]
standard_errors = [0.05] * 5

# Credible set with Bayes factors scaled by the standard errors
cs = credible_sets.compute_credible_set(z_scores, coverage=0.95, standard_errors=standard_errors)
print(cs["snps_in_set"], cs["variance_model"])  # [2] standard_errors

# Without standard errors, model V from the sample size, MAF and case fraction
cs = credible_sets.compute_credible_set(z_scores, sample_size=20_000, maf=0.2, case_fraction=0.4)

# Per-variant log Bayes factors, combined into coloc posteriors
log_bf_1 = compute_log_bayes_factors(z_scores, standard_errors=standard_errors)
log_bf_2 = compute_log_bayes_factors([0.1, -0.7, 5.2, 4.8, 0.3], standard_errors=0.04)
coloc = colocalization_from_log_bayes_factors(log_bf_1, log_bf_2)

# Multi-trait colocalization with the standard errors given per trait
multi = colocalization.multi_trait_coloc(
    {"gwas": z_scores, "eqtl": [0.1, -0.7, 5.2, 4.8, 0.3]},
    standard_errors={"gwas": 0.05, "eqtl": 0.04},
)

# Multi-causal fine-mapping from a genotype matrix (samples x variants) and phenotypes
susie = credible_sets.susie_regression(genotype_matrix, phenotypes, L=5)

results = eqtl.cis_eqtl_scan(expression, genotypes, gene_pos, var_pos)
```

## eQTL Scans

`cis_eqtl_scan(expression_matrix, genotype_matrix, gene_positions, variant_positions, cis_window=1_000_000, maf_threshold=0.05)` takes expression (genes x samples), genotype dosages (variants x samples) and the position tables `[gene_id, chrom, tss_position]` and `[variant_id, chrom, position]`. It returns one row per tested gene-variant pair with `gene_id`, `variant_id`, `beta`, `se`, `pvalue`, `distance` (variant position minus TSS), `maf`, `status` and `error`.

- **Unique identifiers**: `cis_eqtl_scan()` raises `ValueError` when a gene id repeats in `gene_positions`, or a variant id repeats in `variant_positions` or in the index of `genotype_matrix`. A repeated id would test one pair twice or run a gene against another variant's genotypes. `metainformant.eqtl.workflow.synthetic.create_synthetic_genotypes()` builds its variant ids as `var_<full gene id>_<i>`, so they are unique whenever the gene ids are.
- **Failed fits**: a fit that cannot be made (fewer than 3 complete observations, a constant genotype, or a zero standard error because the expression is constant or fits the genotype exactly) has NaN for `beta`, `se` and `pvalue`, `status == "error"` and the reason in `error`. `conditional_eqtl()` (columns `variant_id`, `conditioned_on`, `beta`, `se`, `pvalue`, `status`, `error`) and `eqtl_effect_sizes()` follow the same rule, and `eqtl_effect_sizes()` also sets `r_squared` to NaN. NaN is missing, not an effect of 0 or a p-value of 1. A failed fit is never below the threshold of `trans_eqtl_scan()`, so that scan does not report it.
- **Summaries**: `eqtl_summary_stats()` counts in `n_tests` only the tests that have a p-value and reports the others as `n_tests_failed`; the Benjamini-Hochberg family is the tests that ran. It also returns `n_egenes`, `n_eqtls` and, when there are tests, `fdr_threshold` and `mean_effect_size`.
- **Writing**: `write_eqtl_results(results, output_path)` writes the table of any of the scans tab-separated, creates the parent directories and writes every missing value, such as the NaN of a failed fit, as `NA`.
