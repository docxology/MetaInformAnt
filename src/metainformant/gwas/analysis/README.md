# GWAS Analysis

Statistical methods for GWAS: association tests, multiple testing correction, VCF quality control, population structure, LD pruning, heritability estimation, and variant annotation.

## Contents

| File | Purpose |
|------|---------|
| `association.py` | Linear and logistic regression association tests |
| `mixed_model.py` | EMMA mixed-model association with kinship correction |
| `correction.py` | Bonferroni, FDR, genomic control, q-value corrections |
| `quality.py` | VCF parsing, MAF/missingness/HWE filtering, haplodiploidy QC |
| `structure.py` | PCA, kinship matrices (VanRaden, IBS, Yang), population clustering |
| `ld_pruning.py` | LD-based SNP pruning by r-squared threshold |
| `heritability.py` | REML heritability estimation with chromosome partitioning |
| `calling.py` | Variant calling via bcftools, GATK, FreeBayes; VCF merging/indexing |
| `annotation.py` | Variant-to-gene annotation and functional location classification |
| `summary_stats.py` | Summary statistics output, significant hit extraction, lambda GC |

## Key Functions

| Function | Description |
|----------|-------------|
| `association_test_linear()` | Linear regression GWAS for quantitative traits |
| `association_test_logistic()` | Logistic regression GWAS for case-control traits |
| `association_test_mixed()` | EMMA mixed-model with kinship correction |
| `bonferroni_correction()` | Bonferroni p-value adjustment |
| `fdr_correction()` | Benjamini-Hochberg FDR correction |
| `parse_vcf_full()` | Complete VCF parsing into genotype/variant dicts |
| `apply_qc_filters()` | Apply MAF, missingness, HWE filters to VCF data |
| `compute_pca()` | Principal component analysis of genotype matrix |
| `compute_kinship_matrix()` | Kinship matrix by VanRaden, IBS, or Yang method |
| `ld_prune()` | LD pruning with sliding window and r-squared threshold |
| `estimate_heritability()` | REML narrow-sense heritability from kinship + phenotypes |

## Usage

```python
from metainformant.gwas.analysis.quality import parse_vcf_full, apply_qc_filters

vcf_data = parse_vcf_full("variants.vcf")
filtered = apply_qc_filters(vcf_data, maf=0.05, missing=0.1)
```
