# GWAS Quick Reference

Complete association testing workflow from QC to fine-mapping.

## When to Use

Use `run_gwas` for genome-wide association studies on SNP genotype data—not for eQTL colocalization (see `run_eqtl`) or rare variant aggregation (use `run_skat`). Supports linear, logistic, and Cox models.

## Table of Contents

- [Minimal Working Example](#minimal-working-example)
- [End-to-End Pipeline](#end-to-end-pipeline)
- [Models](#models)
- [Output Files](#output-files)
- [Flags](#flags)
- [Interpreting Results](#interpreting-results)
- [Expected Output](#expected-output)
- [Common Pitfalls](#common-pitfalls)

---

```python-snippet
from metainformant.gwas import run_association

# One-line association
results = run_association(
    vcf="data/genotypes.vcf.gz",
    phenotypes="data/traits.tsv",
    covariates="data/covariates.tsv"
)

# View top hits
results.top_hits(n=10).to_csv("gwas_hits.tsv", sep='	')
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

### REGENIE (efficient mixed-model for large cohorts)
```python-snippet
from metainformant.gwas import regenie

# Step 1: Build genetic relationship matrix (GRM)
grm = regenie.build_grm(
    vcf="cohort.vcf.gz",
    sample_sheet="samples.tsv",
    chunk_size=10000
)

# Step 2: Association testing with covariates
results = regenie.associate(
    grm=grm,
    phenotype="phenotypes.tsv",
    covariates="covariates.tsv",
    model="linear"
)
results.to_csv("regenie_association.tsv", sep='\t')
```
Expected output:
```
[regenie] Step 1: Building GRM from 18374 variants... done (12m34s)
[regenie] Step 2: Testing 5 traits on 10 PCs + age + sex...
[regenie]   Trait 1/5: 12,456,789 tests, 23 hits (p < 5e-8)
[regenie]   Trait 2/5: 12,456,789 tests, 0 hits
```

### Fine-mapping with SuSiE (Summarized Individual Effect)
```python-snippet
from metainformant.gwas.finemapping import susie

# Load summary statistics and LD matrix
sumstats = susie.load_sumstats("gwas_results.tsv")
ld_matrix = susie.compute_ld("filtered.vcf.gz", sample_size=10000)

# Run SuSiE with 10 causal configurations
credible_sets = susie.fine_map(
    sumstats=sumstats,
    ld=ld_matrix,
    n_causal=10
)
credible_sets.to_csv("susie_credible_sets.tsv", sep='\t')
print(f"Found {len(credible_sets)} credible sets, 95% CS coverage")
```
Expected output:
```
Found 18 credible sets, 95% CS coverage
Lead SNP chr6:28538123 (HLA region) → CS size = 12 variants
```

### Colocalization with eQTL (COLOC)
```python-snippet
from metainformant.gwas.coloc import coloc

# GWAS summary statistics per gene
gwas_per_gene = coloc.summarize_by_gene(
    gwas="gwas_results.tsv",
    gtf="annotation.gtf",
    window_kb=500
)

# Colocalize with eQTL from GTEx
 coloc_results = coloc.run(
    gwas_signal=gwas_per_gene,
    eqtl_file="gtex_liver_eqtl.tsv",
    sample_overlap=0.3
 )
 coloc_results.to_csv("coloc_hits.tsv", sep='\t')
```
Expected output:
```
Testing 1,234 genes near GWAS hits...
PP4 (shared causal) > 0.8 for 23 genes
Top hit: AMEL_004123 (PP4=0.96, log10BF=12.4)
```

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

### Summary statistics (top 10 hits)
```
> head gwas_hits.tsv
SNP            CHR   POS        REF   ALT   BETA   SE    P      P_BOLT_LMM
rs123456789    1     1456673   A     G     0.234  0.032 1.2e-08 8.7e-09
rs987654321    6     28538123  C     T    -0.187  0.028 3.4e-08 2.1e-08
...
Genomic inflation λ: 1.02
LD Score intercept: 1.01 (s.e. 0.03)
```

### Fine-mapping credible set output
```
> head susie_credible_sets.tsv
credible_set  variant         rs_id      chromosome  position  posterior  in_CS
CS1           chr1:1456673    rs123456   1           1456673   0.87       True
CS1           chr1:1456721    rs123457   1           1456721   0.08       True
CS2           chr6:28538123   rs987654   6           28538123  0.94       True
...
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
