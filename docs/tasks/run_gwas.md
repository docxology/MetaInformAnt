# GWAS Quick Reference

Complete association testing workflow from QC to fine-mapping.

## When to Use

Use `run_gwas` for genome-wide association studies on SNP genotype data—not for eQTL colocalization (see `run_eqtl`) or rare variant aggregation (use `run_skat`). Supports linear, logistic, and Cox models.

## Table of Contents

- [End-to-End Pipeline](#end-to-end-pipeline)
- [End-to-End Pipeline](#end-to-end-pipeline)
- [Models](#models)
- [Output Files](#output-files)
- [Flags](#flags)
- [Interpreting Results](#interpreting-results)
- [Expected Output](#expected-output)
- [Common Pitfalls](#common-pitfalls)

---

```python-snippet
from metainformant.gwas import run_gwas

# One-call workflow: QC, PCA/kinship, association, multiple-testing correction, plots
results = run_gwas(
    vcf_path="data/genotypes.vcf.gz",
    phenotype_path="data/traits.tsv",
    config={"model": "linear", "trait": "height"},  # model: linear | logistic | mixed
    output_dir="output/gwas",
)

# View top hits
top_hits = sorted(results["results"]["association_results"], key=lambda r: r["p_value"])[:10]
for hit in top_hits:
    print(hit["variant_id"], hit["chrom"], hit["pos"], hit["beta"], hit["p_value"])
```

## End-to-End Pipeline

```bash
# 1. QC & filter VCF
python3 scripts/gwas/qc/filter_variants.py \
  --input raw.vcf.gz \
  --output filtered.vcf.gz \
  --maf 0.05 \
  --missingness 0.1

# 2. Population structure (PCA)
python3 scripts/gwas/structure/pca.py \
  --vcf filtered.vcf.gz \
  --output pca.tsv \
  --components 10

# 3. Association testing
python3 scripts/gwas/association/run_all.py \
  --vcf filtered.vcf.gz \
  --pheno data/traits.tsv \
  --covariates pca.tsv \
  --method linear  # or logistic for binary traits

# 4. Multiple testing correction
python3 scripts/gwas/qc/correct_pvalues.py \
  --input gwas_results.tsv \
  --method fdr  # or bonferroni, holm
```

## Models

| Trait Type | Method | Implementation |
|------------|--------|---------------|
| Quantitative (height, weight) | Linear regression | `--model linear` (SAIGE, PLINK) |
| Binary (disease yes/no) | Logistic regression | `--model logistic` |
| Survival (time-to-event) | Cox proportional hazards | `--model cox` |

## Output Files

```
output/gwas/
 results.tsv # Raw association statistics
 manhattan.png # Manhattan plot
 qq.png # QQ plot (inflation check)
 top_hits.tsv # Significant SNPs (p < 5e-8)
 lambdas.tsv # Genomic inflation factors
```

## Flags

```bash
--maf MIN               # Minor allele frequency filter (default: 0.01)
--geno MISSING          # Missingness filter per SNP (default: 0.1)
--mind MISSING          # Missingness filter per sample (default: 0.1)
--hwe P                 # Hardy-Weinberg equilibrium p-value threshold
--threads N             # Parallel threads (default: auto)
--covariates FILE       # Covariate matrix (PCA, batch, sex, age, etc.)
--method NAME           # linear | logistic | cox
--output-prefix PREFIX  # Output file prefix
```

## Interpreting Results

**Genomic inflation (λ):**
- λ ≈ 1.0 → Well-controlled population structure
- λ > 1.05 → Residual confounding, add more PCs

**Manhattan plot:**
- Each dot = SNP; y-axis = -log10(p-value)
- Genome-wide significance threshold: 5e-8
- Chromosome bands alternate colors

**QQ plot:**
- Points follow diagonal = well-calibrated p-values
- Inflation at top = true associations
- Whole-plot lift = population stratification

## Advanced Examples

### Credible set from summary statistics (Wakefield approximate Bayes factors)
```python-snippet
import numpy as np

from metainformant.gwas.finemapping.credible_sets import compute_credible_set

# Z-scores and standard errors for the variants of one locus, in the same order. In practice they come from
# the association results (z = beta / se); here a synthetic locus of 100 variants has one causal variant (index 40)
rng = np.random.default_rng(0)
z = rng.normal(size=100)
z[40] = 8.0
se = np.full(100, 0.05)

# standard_errors fixes the scale of the Bayes factors (V = se**2); the prior variance defaults to 0.15**2
credible = compute_credible_set(z.tolist(), coverage=0.95, standard_errors=se.tolist())
print(
    f"{credible['n_snps']} variant(s) in the 95% credible set: {credible['snps_in_set']} "
    f"({credible['variance_model']}, prior variance {credible['prior_variance']:.4f})"
)
```
Expected output:
```
1 variant(s) in the 95% credible set: [40] (standard_errors, prior variance 0.0225)
```
Without `standard_errors` (or a `sample_size`, with `maf` and, for a case-control trait, `case_fraction`), the Bayes factors fall back to the Z-score scale and `variance_model` reads `z_scale`.

### Fine-mapping with SuSiE (Sum of Single Effects)
```python-snippet
from metainformant.gwas.finemapping.credible_sets import susie_regression

# X: genotype matrix for the locus (n_samples x n_variants); y: phenotype vector
# Fit up to 10 causal effects
fit = susie_regression(X, y, L=10)

print(f"Fitted {len(fit['credible_sets'])} single-effect components (converged: {fit['converged']})")
for i, cs in enumerate(fit["credible_sets"], start=1):
    print(f"CS{i}: {cs['size']} variants ({cs['coverage_achieved']:.0%} coverage)")
```
Expected output:
```
Fitted 10 single-effect components (converged: True)
CS1: 1 variants (100% coverage)
CS2: 9 variants (95% coverage)
...
```

### Colocalization with eQTL (COLOC)
```python-snippet
import numpy as np

from metainformant.gwas.finemapping.colocalization import eqtl_coloc

# Z-scores and standard errors for the same variants (same order) at one locus. In practice gwas_z and gwas_se
# come from the GWAS summary statistics and eqtl_z and eqtl_se from the GTEx eQTL for the gene; here a synthetic
# locus of 100 variants shares one causal variant (index 40)
rng = np.random.default_rng(0)
gwas_z, eqtl_z = rng.normal(size=100), rng.normal(size=100)
gwas_z[40], eqtl_z[40] = 8.0, 7.0
gwas_se, eqtl_se = np.full(100, 0.05), np.full(100, 0.10)

result = eqtl_coloc(
    gwas_z=gwas_z.tolist(),
    eqtl_z=eqtl_z.tolist(),
    gene_id="AMEL_004123",
    gwas_standard_errors=gwas_se.tolist(),  # the standard errors fix the scale of each Bayes factor
    eqtl_standard_errors=eqtl_se.tolist(),
)
print(
    f"{result['gene_id']}: PP4 (shared causal) = {result['pp_h4']:.2f}, colocalized: {result['is_colocalized']}, "
    f"shared variant: {result['shared_variant']}"
)
```
Expected output:
```
AMEL_004123: PP4 (shared causal) = 1.00, colocalized: True, shared variant: 40
```
`result["coloc_result"]` holds `PP_H0` to `PP_H4`, `snp_pp_h4` (the per-variant probability of being the shared causal variant given H4) and `bayes_factor_model`, which names the variance model and prior variance used for each trait.

## Expected Output

### Association test console log
```
[2026-04-26 09:00:01] Loading VCF: filtered.vcf.gz (12.3M variants)
[2026-04-26 09:02:45] Loading phenotypes: traits.tsv (n=8,432)
[2026-04-26 09:03:12] Computing principal components...
[2026-04-26 09:05:44] PC1 explains 4.2% variance, PC2: 2.1%
[2026-04-26 09:06:00] Starting linear association (SAIGE)
[2026-04-26 09:20:11] Chromosome 1 complete: 1,234,567 tests
[2026-04-26 10:45:33] Chromosome 16 complete: 987,654 tests
[2026-04-26 11:30:01] All chromosomes done: 12,456,789 tests total
[2026-04-26 11:32:22] Multiple testing correction (FDR < 0.05)
[2026-04-26 11:33:44] Significant hits: 142 (p < 5e-8), 1,203 (FDR < 0.05)
[2026-04-26 11:35:01] Writing results: gwas_results.tsv
[2026-04-26 11:36:12] Generating Manhattan plot...
[2026-04-26 11:38:45] Done: manhattan.png, qq.png
```

## Common Pitfalls

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| `SAIGE convergence failed` | Null model didn't converge (rare phenotype or extreme imbalance) | Increase `--step-size`, check phenotype distribution; if binary trait has < 1% cases, use `--firth` or Fisher's exact test |
| `No significant hits` despite strong phenotype | Underpowered (n < 500), high genomic inflation (λ > 1.1) | Increase sample size; add more PCs; remove outlier samples; verify phenotype measurement quality |
| `MemoryError` on VCF load | 100k+ samples × millions of SNPs → >100 GB RAM | Use BGEN format (Plink2), downsample SNPs (--thin 10), or use REGENIE (two-step, lower memory) |
| `Sample mismatch` between VCF and phenotype file | Different sample IDs, ordering, or subset | Run `plink --bgen ... --sample ... --make-bed` to align; check `--keep` filters; ensure consistent sample ordering |
| Manhattan plot shows spike at chr1 (or chrX) | Residual population stratification not fully corrected | Add more PCs (top 20), use linear mixed model (SAIGE/BOLT-LMM), or remove ancestry outliers from PCA |
| `LD matrix singular` | Too many SNPs relative to samples for LD calculation | Prune SNPs first (`plink --indep-pairwise 50 5 0.2`), or use shrinkage estimator (`--ld-shrink 0.5`) |

---

**Related:** [Full GWAS docs](../gwas/index.md) | Fine-mapping guide | [eQTL integration](../eqtl/) | [Multi-omics](../multiomics/)
