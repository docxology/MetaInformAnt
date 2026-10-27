# GWAS Analysis

Statistical methods for GWAS: association tests, multiple testing correction, VCF quality control, population structure, LD pruning, heritability estimation, and variant annotation.

## Contents

| File | Purpose |
|------|---------|
| `association.py` | Linear (Student's t) and logistic (Wald) association tests on complete cases |
| `mixed_model.py` | EMMA mixed-model association with kinship correction |
| `correction.py` | Bonferroni, FDR, genomic control, q-value corrections |
| `quality.py` | VCF parsing, QUAL/MAF/missingness/HWE filtering, haplodiploidy QC, filtered-VCF output |
| `structure.py` | PCA, kinship matrices (VanRaden, IBS, Yang), population clustering |
| `ld_pruning.py` | LD-based SNP pruning by r-squared threshold |
| `ld_decay.py` | Pairwise r-squared against distance in 20 distance bins, the distance at which the mean falls to 0.5, LD decay plot |
| `ldsr.py` | LD scores, LD score regression and its plot |
| `power.py` | Power calculations, subsampling convergence, block-jackknife standard errors |
| `heritability.py` | REML heritability estimation with chromosome partitioning |
| `eqtl.py` | Linear eQTL scan of genotype columns against expression columns, optionally restricted to cis pairs |
| `calling.py` | Variant calling via bcftools, GATK, FreeBayes; bcftools filtering; VCF merging/indexing |
| `annotation.py` | Variant-to-gene annotation and functional location classification |
| `summary_stats.py` | Summary statistics TSV, significant hit extraction, JSON results summary with lambda GC |
| `utils.py` | Shared analysis utilities and helper functions |

## Key Functions

| Function | Description |
|----------|-------------|
| `association_test_linear()` | Linear regression GWAS for quantitative traits |
| `association_test_logistic()` | Logistic regression GWAS for case-control traits |
| `association_test_mixed()` | EMMA mixed-model with kinship correction |
| `bonferroni_correction()` | Bonferroni p-value adjustment |
| `fdr_correction()` | Benjamini-Hochberg FDR correction |
| `qvalue_estimation()` | Storey q-values with an estimated or given pi0 |
| `adjust_p_values()` | One entry point for the Bonferroni, FDR, genomic control and q-value adjustments |
| `as_p_value()` | A p-value as a finite float, or `None` when it is missing (NaN, infinite, `None`, boolean or non-numeric) |
| `parse_vcf_full()` | Complete VCF parsing into genotype/variant dicts |
| `apply_qc_filters()` | Apply QUAL, MAF, missingness and HWE filters to VCF data |
| `filter_variants_bcftools()` | `bcftools filter` that keeps sites meeting `min_qual`, `max_missing` and `min_depth` |
| `compute_pca()` | Principal component analysis of genotype matrix |
| `compute_kinship_matrix()` | Kinship matrix by VanRaden, IBS, or Yang method |
| `ld_prune()` | LD pruning with sliding window and r-squared threshold |
| `estimate_heritability()` | REML narrow-sense heritability from kinship + phenotypes |
| `run_eqtl_analysis()` | Linear eQTL scan, all-vs-all or restricted to cis pairs |
| `minor_allele_frequency()` | Minor allele frequency of one variant's 0/1/2 dosages |
| `write_summary_statistics()` | Summary statistics TSV with Benjamini-Hochberg `Q_FDR` |
| `write_significant_hits()` | The same table restricted to variants below a p-value threshold |
| `create_results_summary()` | JSON summary: tests run and failed, lambda GC, significance counts, top hits |

## Association Tests

`association_test_linear(genotypes, phenotypes, covariates=None)` and `association_test_logistic(genotypes, phenotypes, covariates=None, max_iter=100)` test one variant across samples. Both analyse complete cases: a sample is dropped when its genotype is negative or non-finite (`-1` marks a missing call) or when a covariate value is non-finite. The linear test also drops samples whose phenotype is missing (`NaN`, `None` or non-numeric), and the logistic test raises `ValueError` for a phenotype that is not 0 or 1. Results report `n_samples` (the samples analysed), `n_missing_dropped` and `converged`.

- **Linear**: p-values come from Student's t distribution with `n - p` degrees of freedom, where `n` is the number of complete cases and `p` the number of fitted parameters including the intercept (a normal approximation when SciPy is unavailable).
- **Logistic**: `converged` is the fit's convergence flag; a fit that did not converge (for example under perfect separation) keeps `status: "success"` but its estimate and p-value are unreliable.

### Failed tests

A test that cannot be fitted is returned, not raised. It carries NaN for its statistics, `status`, `converged=False` and the reason in `error`:

| Test | NaN statistics | `status` of a failed test |
|------|----------------|---------------------------|
| Linear | `beta`, `se`, `t_stat`, `p_value`, `r_squared` | `"error"` |
| Logistic | `beta`, `se`, `z_stat`, `p_value`, `odds_ratio` | `"failed"` (fewer than 2 cases or 2 controls among the complete cases) or `"error"` (the statsmodels fit failed) |
| Mixed model | `beta`, `se`, `t_stat`, `p_value` | `"error"` |

NaN means missing: it is neither an effect of 0 nor a p-value of 1. The corrections, summaries and writers below leave a failed test out, and the writers print it as `NA`. Every row of `run_mixed_model_gwas()` carries `status`.

A linear test fails when fewer than 3 complete cases remain, the genotype is constant, the design matrix is singular or has no residual degrees of freedom, or the phenotype is constant (the error message says so). Zero residual variance is otherwise read for what it means: an exact fit with a nonzero genotype effect reports `se` 0, `t_stat` of plus or minus infinity and `p_value` 0, and zero residual variance with a zero effect is a failed test.

## Multiple-Testing Correction

`fdr_correction()`, `bonferroni_correction()`, `qvalue_estimation()` and `adjust_p_values()` leave a missing p-value (NaN or `None`, the p-value of a failed test) out of the family: the number of tests `m` counts only the p-values that are present. Outputs stay aligned with the input list and hold NaN at a missing entry, and `significant` is False there. The dictionaries of `fdr_correction()` and `bonferroni_correction()` report `n_tests` (p-values present) and `n_excluded` (p-values missing). `genomic_control()` estimates lambda from the p-values in (0, 1] and passes the others through unchanged, and `adjust_p_values(method="genomic_control")` does the same.

## Summaries

- `create_results_summary()` counts in `n_variants_tested` only the tests that ran (those with a p-value) and reports the others as `n_tests_failed`. The inflation factor, significance counts and top hits come from the tests that ran.
- `compute_comprehensive_summary()` reports `n_variants` for the tests that ran and, likewise, `n_tests_failed`.
- `per_chromosome_summary()`, `sign_test()`, `subsample_convergence()` and `jackknife_se()` skip failed tests.

## Quality Control

`apply_qc_filters(vcf_input, qc_config)` accepts a VCF path or the dict from `parse_vcf_full()`. It applies the QUAL, minor allele frequency, missing-data and Hardy-Weinberg filters in that order. `qc_config` keys:

| Key | Default | Meaning |
|-----|---------|---------|
| `min_maf` | 0.01 | Minimum minor allele frequency |
| `max_missing` | 0.1 | Maximum proportion of missing genotypes |
| `min_hwe_p` | 1e-6 | Minimum Hardy-Weinberg p-value |
| `min_qual` | none | Minimum QUAL, inclusive; absent or `None` disables the QUAL filter |

When `min_qual` is set, a variant whose QUAL is missing (`.` in the VCF, `None` once parsed) cannot meet the threshold and is removed, whatever the threshold value. The result's `qual_filter` entry holds `min_qual`, `n_removed` and `n_missing_qual` (the part of `n_removed` with a missing QUAL); leave `min_qual` unset for a VCF without QUAL values. A non-numeric `min_qual` raises `ValueError`.

The genotypes are sample-major (one row per sample). `apply_qc_filters()` raises `ValueError` when those rows differ in length or when the number of variant records differs from the number of genotype columns, whatever `min_qual` is.

`write_filtered_vcf()` writes a missing QUAL, ID, FILTER or INFO as `.`. `filter_variants_bcftools(vcf_path, output_vcf, filters)` joins `QUAL>=min_qual`, `F_MISSING<=max_missing` and `DP>=min_depth` (for the keys present in `filters`) with `&&` and passes them to `bcftools filter -i`, so sites that meet every threshold are kept.

## Summary Statistics

`write_summary_statistics(results, variant_info, output_path, *, genotypes=None)` writes a tab-separated table with the columns `CHR POS SNP REF ALT BETA SE P Q_FDR N MAF`; `write_significant_hits()` writes the same columns for results with `p_value` below `threshold` (default 5e-8), sorted by p-value.

- A missing or non-finite numeric value is written as `NA`, and a missing text field as `.`. The BETA, SE, P and Q_FDR of a failed test are `NA`.
- `Q_FDR` is the result's own `q_value`, `q_fdr` or `Q_FDR` when present. Otherwise it is the Benjamini-Hochberg q-value over the p-values of `results`, so pass every test rather than a subset (`write_significant_hits()` also derives it over the full `results` before selecting hits).
- `MAF` is the result's own `maf` or `MAF` when present. Otherwise it is computed from the matching row of `genotypes` (variants x samples; dosages 0, 1 and 2, any other code treated as missing), and it is `NA` when no genotype is called or `genotypes` is not given.
- `N` is the result's `n_samples` or `N`.

`minor_allele_frequency(genotype_row)` in `summary_stats` is the public helper behind the computed column: it returns the minor allele frequency of one variant's 0/1/2 dosages, or `None` when no genotype is called. `execute_gwas_workflow()` and `run_gwas()` give every association result a `maf` (the frequency of the genotype row the test ran on; a finite `maf` or `MAF` already present is kept, and a variant without a called genotype gets none), and all three workflows pass the genotypes to the writers, so `summary_statistics.tsv` carries a MAF column instead of `NA`.

## Randomness

Functions that draw a random sample take a `seed` (default 42) and use their own generator: the same inputs and seed give the same result, and the global `random` state is neither read nor advanced.

| Function | What the seed draws |
|----------|---------------------|
| `calculate_ld_decay(..., seed=42)` (`summary_stats`) | The variant pairs |
| `compute_ld_decay(..., seed=42)` (`ld_decay`) | The sampled variant pairs, when there are more than `n_pairs_max` |
| `ld_decay_plot(..., seed=42)` | The 5000 pairs scattered, when there are more |
| `subsample_convergence(..., seed=42)` | The variant subsamples |
| `ldsr_plot(..., seed=42)` | The 3000 variants scattered, when there are more |

## eQTL Scan

`run_eqtl_analysis(genotype_matrix, expression_matrix, covariates=None, cis_window=None, max_transcripts=None, *, variant_positions=None, gene_positions=None)` tests every genotype column (variant) against every expression column (transcript) with `association_test_linear()` and returns one row per test.

- Without position tables the scan is all-vs-all.
- With `variant_positions` (columns `variant_id`, `chrom`, `position`) and `gene_positions` (columns `gene_id`, `chrom`, `tss_position`), given together, only cis pairs are tested: the variant is on the transcript's chromosome and `|position - TSS| <= cis_window` (default `DEFAULT_CIS_WINDOW`, 1 000 000 bp, on either side). The result gains a `distance` column, the variant position minus the TSS in bp. A variant or transcript with no row in its table is not tested, and a warning counts them.
- `ValueError` is raised for a `cis_window` without both tables (a window without positions would otherwise be ignored), for only one table, a negative window, a position table that lacks a required column or repeats an identifier, and `max_transcripts` below 1.

## Usage

```python
from metainformant.gwas.analysis.association import association_test_linear
from metainformant.gwas.analysis.quality import apply_qc_filters, parse_vcf_full
from metainformant.gwas.analysis.summary_stats import write_summary_statistics

vcf_data = parse_vcf_full("variants.vcf")
qc = apply_qc_filters(vcf_data, {"min_maf": 0.05, "max_missing": 0.1, "min_qual": 30})
filtered = qc["filtered_data"]
print(qc["num_variants_before"], qc["num_variants_after"], qc.get("qual_filter"))

genotypes = [0, 1, 2, 1, 0, 2, -1, 1]  # -1 marks a missing call
phenotypes = [0.1, 0.9, 2.1, 1.2, -0.2, 1.8, 0.5, 1.1]
variant_info = [{"chrom": "1", "pos": 1042, "id": "rs1", "ref": "A", "alt": ["G"]}]

result = association_test_linear(genotypes, phenotypes)
write_summary_statistics([result], variant_info, "output/gwas/summary_statistics.tsv", genotypes=[genotypes])
```
