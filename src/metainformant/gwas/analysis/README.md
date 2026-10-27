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
| `heritability.py` | REML heritability estimation with chromosome partitioning |
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
| `parse_vcf_full()` | Complete VCF parsing into genotype/variant dicts |
| `apply_qc_filters()` | Apply QUAL, MAF, missingness and HWE filters to VCF data |
| `filter_variants_bcftools()` | `bcftools filter` that keeps sites meeting `min_qual`, `max_missing` and `min_depth` |
| `compute_pca()` | Principal component analysis of genotype matrix |
| `compute_kinship_matrix()` | Kinship matrix by VanRaden, IBS, or Yang method |
| `ld_prune()` | LD pruning with sliding window and r-squared threshold |
| `estimate_heritability()` | REML narrow-sense heritability from kinship + phenotypes |
| `write_summary_statistics()` | Summary statistics TSV with Benjamini-Hochberg `Q_FDR` |
| `write_significant_hits()` | The same table restricted to variants below a p-value threshold |

## Association Tests

`association_test_linear(genotypes, phenotypes, covariates=None)` and `association_test_logistic(genotypes, phenotypes, covariates=None, max_iter=100)` test one variant across samples. Both analyse complete cases: a sample is dropped when its genotype is negative or non-finite (`-1` marks a missing call) or when a covariate value is non-finite. The linear test also drops samples whose phenotype is missing (`NaN`, `None` or non-numeric), and the logistic test raises `ValueError` for a phenotype that is not 0 or 1. Results report `n_samples` (the samples analysed), `n_missing_dropped` and `converged`.

- **Linear**: p-values come from Student's t distribution with `n - p` degrees of freedom, where `n` is the number of complete cases and `p` the number of fitted parameters including the intercept (a normal approximation when SciPy is unavailable). A degenerate design (fewer than 3 complete cases, a constant genotype or a singular design matrix) returns `status: "error"` with `p_value` 1.0 and an `error` message.
- **Logistic**: `status` is `"success"`, `"failed"` (fewer than 2 cases or 2 controls among the complete cases) or `"error"` (the statsmodels fit failed). `converged` is the fit's convergence flag; a fit that did not converge (for example under perfect separation) keeps `status: "success"` but its estimate and p-value are unreliable.

## Quality Control

`apply_qc_filters(vcf_input, qc_config)` accepts a VCF path or the dict from `parse_vcf_full()`. It applies the QUAL, minor allele frequency, missing-data and Hardy-Weinberg filters in that order. `qc_config` keys:

| Key | Default | Meaning |
|-----|---------|---------|
| `min_maf` | 0.01 | Minimum minor allele frequency |
| `max_missing` | 0.1 | Maximum proportion of missing genotypes |
| `min_hwe_p` | 1e-6 | Minimum Hardy-Weinberg p-value |
| `min_qual` | none | Minimum QUAL, inclusive; absent or `None` disables the QUAL filter |

When `min_qual` is set, a variant whose QUAL is missing (`.` in the VCF, `None` once parsed) cannot meet the threshold and is removed, whatever the threshold value. The result's `qual_filter` entry holds `min_qual`, `n_removed` and `n_missing_qual` (the part of `n_removed` with a missing QUAL); leave `min_qual` unset for a VCF without QUAL values. A non-numeric `min_qual` raises `ValueError`.

`write_filtered_vcf()` writes a missing QUAL, ID, FILTER or INFO as `.`. `filter_variants_bcftools(vcf_path, output_vcf, filters)` joins `QUAL>=min_qual`, `F_MISSING<=max_missing` and `DP>=min_depth` (for the keys present in `filters`) with `&&` and passes them to `bcftools filter -i`, so sites that meet every threshold are kept.

## Summary Statistics

`write_summary_statistics(results, variant_info, output_path, *, genotypes=None)` writes a tab-separated table with the columns `CHR POS SNP REF ALT BETA SE P Q_FDR N MAF`; `write_significant_hits()` writes the same columns for results with `p_value` below `threshold` (default 5e-8), sorted by p-value.

- A missing or non-finite numeric value is written as `NA`, and a missing text field as `.`.
- `Q_FDR` is the result's own `q_value`, `q_fdr` or `Q_FDR` when present. Otherwise it is the Benjamini-Hochberg q-value over the p-values of `results`, so pass every test rather than a subset (`write_significant_hits()` also derives it over the full `results` before selecting hits).
- `MAF` is the result's own `maf` or `MAF` when present. Otherwise it is computed from the matching row of `genotypes` (variants x samples; dosages 0, 1 and 2, any other code treated as missing), and it is `NA` when no genotype is called or `genotypes` is not given.
- `N` is the result's `n_samples` or `N`.

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
