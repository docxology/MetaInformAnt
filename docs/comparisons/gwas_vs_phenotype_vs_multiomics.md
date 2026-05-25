# GWAS vs Phenotype vs Multi-Omics: Study Design & Analysis Comparison

## Overview

This guide compares three complementary approaches for connecting genetic variation to biological phenotypes: Genome-Wide Association Studies (GWAS), phenotype analysis, and multi-omics integration. Understanding their strengths, limitations, and appropriate use cases is critical for study design and interpretation.

## Module Summary

| Module | Core Purpose | Data Type | Sample Size Typical | Primary Output |
|--------|--------------|-----------|--------------------|----------------|
| **[gwas](../gwas/)** | Identify genetic variants associated with traits | Genotypes (VCF), phenotypes | 1,000–1,000,000+ (population-scale) | Association statistics (p-values, effect sizes) |
| **[phenotype](../phenotype/)** | Characterize and analyze phenotypic traits | Morphological, behavioral, chemical measurements | 10–10,000 (individual-level) | Trait profiles, measurements, indices |
| **[multiomics](../multiomics/)** | Integrate multiple omic layers (DNA+RNA+Protein) | Multi-modal (genomics, transcriptomics, proteomics) | 10–1,000 (coordinated multi-omic) | Integrated components, cross-omic associations |
| **[eqtl](../eqtl/)** *(cross-cutting)* | Bridge GWAS variants with gene expression | GWAS summary stats + RNA expression | 100–10,000 (expression cohorts) | eQTLs, colocalization, mediated effects |

> **Note**: eQTL analysis is a cross-cutting integration pipeline (logic in `gwas/finemapping/` and `multiomics/`) — it leverages both GWAS and RNA data but is not a standalone module.

## 1. Study Design Comparison

### Research Question Alignment

| Question Type | Best Module(s) | Rationale |
|---------------|----------------|-----------|
| **"Which SNPs are associated with height?"** | `gwas` | Direct variant-trait association via regression |
| **"What is the heritability of diabetes?"** | `gwas` + `phenotype` | SNP-heritability from GWAS; phenotype definition |
| **"How do gene expression changes mediate genetic effects?"** | `eqtl` + `multiomics` | Links DNA→RNA→phenotype via mediation |
| **"What morphological traits define this species?"** | `phenotype` | Direct measurement without genotypes |
| **"Which multi-omic signature predicts disease?"** | `multiomics` | Integrates DNA+RNA+Protein for prediction |
| **"Do genetic variants affect gene expression?"** | `eqtl` | cis/trans eQTL scanning |
| **"What pathways are enriched among associated genes?"** | `gwas` → `ontology` | Post-GWAS enrichment analysis |

### Study Scale Requirements

| Scale Dimension | gwas | phenotype | multiomics | eQtl |
|-----------------|------|-----------|------------|------|
| **Samples needed** | 1,000 – 1,000,000+ (large cohorts) | 10 – 10,000 (phenotyping possible) | 10 – 1,000 (expensive multi-omic) | 100 – 10,000 (expression cohorts) |
| **Variants tested** | 10⁵ – 10⁷ (genome-wide) | N/A | 10⁴ – 10⁶ (subset) | Matched to expression SNPs |
| **Phenotype measurements** | 1 – 20 (typical) | 1 – 1,000+ (comprehensive) | 1 – 20 (harmonized) | Same as GWAS phenotype |
| **Cost** | Low-medium ($50–200/sample genotyping) | Medium ($100–500/individual phenotyping) | Very high ($1,000–10,000/sample multi-omic) | High ($1,000–5,000/sample RNA-seq + genotypes) |
| **Time to complete** | Months–years | Months–years | 1–2 years | 6–12 months |

> **Rule of thumb**: GWAS scales cheaply with samples via array genotyping or imputation; multi-omics costs rise steeply with depth.

### Data Type & Format Comparison

| Data Aspect | gwas | phenotype | multiomics |
|-------------|------|-----------|------------|
| **Genotype source** | VCF, BGEN, PLINK .bed/.bim/.fam | Optional (only if integrated) | VCF required (DNA layer) |
| **Phenotype format** | Single column (continuous) or binary | Full measurement tables (multi-trait) | Single or multi-trait |
| **Covariates** | Age, sex, PC1–PC10 (population structure) | Missing value imputation, batch | Same as GWAS + omics-specific |
| **Metadata** | Cohort, ancestry, relatedness | Specimen ID, measurement protocol | Sample ID, tissue, batch, platform |
| **Input files** | `samples.vcf.gz`, `phenotypes.tsv`, `covariates.tsv` | `measurements.tsv`, `metadata.json` | `genomics.vcf`, `transcriptomics.tsv`, `proteomics.tsv` |

## 2. Statistical Power & Sensitivity

| Aspect | gwas | phenotype | multiomics |
|--------|------|-----------|------------|
| **Power driver** | Sample size * minor allele frequency | Measurement precision + n | Both sample size AND multi-omic depth |
| **Effect detectable** | OR 1.05–1.2 (common variants, large n) | Cohen's d 0.2–0.8 | Integrated effect size depends on modality |
| **Multiple testing** | Genome-wide Bonferroni (~5×10⁻⁸) | Bonferroni per trait or FDR | Per-feature across multi-view |
| **Confounders** | Population stratification, cryptic relatedness | Measurement error, observer bias | Batch effects across platforms |
| **Replication needed** | Independent cohort (high) | Biological replicates (medium) | Cross-modal validation (high) |

### Power Calculator Guidelines

- **GWAS**: Need 80% power to detect OR=1.1 at MAF=0.2 → ~50,000 samples
- **Phenotype**: Effect size d=0.5, α=0.05 → ~64 samples/group
- **Multi-omics**: For joint PCA to capture cross-modal signal → at least 50–100 samples with all omics measured

## 3. Integration Strategies

### How do these modules work together?

```
┌─────────────────────────────────────────────────────────────────────┐
│                       Analysis Flow                                 │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│   ┌─────────────┐        ┌──────────────┐        ┌─────────────┐  │
│   │ Phenotyping │ ────── │   GWAS       │ ────── │  eQTL       │  │
│   │ (module:    │        │ (module:     │        │ (pipeline:  │  │
│   │ phenotype)  │        │ gwas)        │        │ gwas+multi) │  │
│   └─────────────┘        └──────────────┘        └──────┬──────┘  │
│                                                          │        │
│                                   ┌──────────────────────▼────────┐│
│                                   │  Multi-Omic Integration         ││
│                                   │  (module: multiomics)           ││
│                                   │  • Joint PCA/CCA/NMF             ││
│                                   │  • Multi-view clustering          ││
│                                   │  • Cross-omics biomarker disc.    ││
│                                   └───────────────────────────────────┘│
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

### Integration Entry Points

| From → To | Connection Method | Example |
|-----------|------------------|---------|
| phenotype → gwas | Phenotype file input to `gwas.workflow` | `phenotypes.tsv` → trait column |
| gwas → multiomics | Variant file → `multiomics.from_dna_variants()` | VCF → DNA layer |
| multiomics → gwas | Summary stats → fine-mapping | LD matrix from multi-omics |
| gwas → eqtl | Summary statistics → colocalization | coloc, eCAVIAR |
| phenotype → multiomics | Trait as multi-omic outcome | Multi-view regression |

## 4. Comparison Tables

### Study Design Dimension

| Criterion | gwas | phenotype | multiomics |
|-----------|------|-----------|------------|
| **Hypothesis** | Variant → trait association | Trait characterization | Multi-causal layer discovery |
| **Directionality** | DNA→Phenotype (observational) | Measurement only | Bi-directional across omics |
| **Confounder control** | PCA, kinship, mixed models | Replication + standardization | Batch correction across platforms |
| **Causal inference** | Limited (observational) | Descriptive only | Mediation analysis possible |
| **Sample origin** | Population cohorts or biobanks | Field/lab measurements | Coordinated multi-omic sampling |
| **Ethical constraints** | Consent for genetics + data sharing | Standard human/animal ethics | Multi-tiered consent (genomics + transcriptomics) |

### Data Requirements Table

| Requirement | gwas | phenotype | multiomics |
|-------------|------|-----------|------------|
| **Genotype data** | **Mandatory** (VCF/BGEN) | Optional | Mandatory (DNA layer) |
| **Phenotype data** | **Mandatory** (single trait column) | **Mandatory** (full table) | Optional (can use as outcome) |
| **RNA expression** | Optional (for eQTL) | Optional | Optional (transcriptome layer) |
| **Protein data** | Never | Rarely | Optional (proteome layer) |
| **Metadata (age/sex)** | Strongly recommended | Recommended | Recommended |
| **Batch variables** | Optional (but critical if present) | Critical for measurement | Absolutely critical |
| **Sample overlap** | Same samples needed | Independent OK | Must be same individuals |

### Statistical Method Comparison

| Method Category | gwas | phenotype | multiomics |
|-----------------|------|-----------|------------|
| **Association** | Linear/logistic regression, MLM | Correlations, MANOVA | Multi-view PLS, MOFA |
| **Dim. reduction** | PCA (for confounders only) | PCA (for summary) | Joint PCA, CCA, NMF (core) |
| **Feature selection** | Clumping, thresholding | Correlation filtering | Multi-view feature selection |
| **Multiple test** | Bonferroni, FDR, genomic control | FDR per trait | Multi-omic FDR (weighted) |
| **Effect size** | β, OR, SE, p-value | Cohen's d, eta² | Multi-view loadings |
| **Replication** | Independent cohort split | Technical/biological replicates | Cross-validation + external validation |

### Visualization Type Comparison

| Plot Type | gwas | phenotype | multiomics |
|-----------|------|-----------|------------|
| **Manhattan** | ✓ Primary association visualization | ✗ | ✗ |
| **Q-Q** | ✓ Inflation diagnostic | ✗ | ✗ |
| **Regional** | ✓ Locus zoom | ✗ | ✗ |
| **Scatter matrix** | ✗ | ✓ Correlation matrix | ✓ Cross-omics loading scatter |
| **Heatmap** | ✗ | ✓ Trait correlations | ✓ Multi-omic feature matrix |
| **PCA** | ✓ Population structure | ✓ Sample summary | ✓ Joint PCA (multi-view) |
| **Volcano** | ✓ Effect vs significance | ✓ DE among phenotypes | ✓ Multi-omic signature |
| **Network** | ✗ | ✗ | ✓ Multi-layer network (omics interactions) |

## 5. When to Choose Each Approach

### Choose **GWAS** when:

- ✓ You have **large sample size** (>1,000) with genotype data
- ✓ Your primary goal is identifying **genetic loci** affecting a trait
- ✓ Phenotype is **well-defined** (binary or quantitative)
- ✓ You can control for population structure (need ancestry PCs or kinship)
- ✓ Cost is constrained (genotyping arrays cheaper than multi-omic)
- ✓ You want to leverage **public summary statistics** (UK Biobank, etc.)

**Avoid GWAS when**:
- ✗ Your phenotype is not easily quantifiable or has high measurement error
- ✗ Sample size is small (<1,000) — underpowered for genome-wide tests
- ✗ You need mechanistic insights beyond association
- ✗ Genetic architecture is rare variants only (need sequencing, larger n)

---

### Choose **Phenotype Analysis** when:

- ✓ You are **characterizing traits** without genetic data (or genetics comes later)
- ✓ Detailed **morphological, behavioral, or chemical** profiling is the goal
- ✓ High-precision measurements possible (calipers, spectrometers, tracking)
- ✓ You need to define **derived indices** (allometric ratios, shape indices)
- ✓ Building **phenotype-genotype maps** for downstream integration

**Avoid phenotype-only** when:
- ✗ You already have genotypes and want to test associations directly (use GWAS)
- ✗ Phenotype measurement is too noisy or subjective (unless well-controlled)
- ✗ Traits are purely genetic (monogenic) — consider DNA diagnostics instead

---

### Choose **Multi-Omics Integration** when:

- ✓ You have **coordinated multi-layer data** (DNA + RNA + Protein) on same samples
- ✓ Goal is **systems-level understanding** (pathways, networks, multi-omic signatures)
- ✓ Sufficient sample size (≥50–100) across all omics layers
- ✓ You can afford deep profiling costs (or have existing biobank with multi-omics)
- ✓ Want to identify **cross-omics interactions** (e.g., variant → expression → protein)

**Avoid multi-omics** when:
- ✗ Only one or two omic layers available (use single-module: dna or rna)
- ✗ Sample size is tiny (<20 per omic) — overfitting risk
- ✗ Missing data is extensive across layers
- ✗ Cost constraints preclude deep profiling

---

### Choose **eQTL Integration** when:

- ✓ You have **both GWAS and RNA-seq** data from overlapping samples
- ✓ Goal is to **link genetic associations to gene expression**
- ✓ Want to **mechanistically interpret** GWAS hits (which gene is causal?)
- ✓ Have sufficient expression sample size (≥100) for cis-eQTL power
- ✓ Interested in **colocalization** (same variant drives both)

**Note**: eQTL leverages `gwas` and `multiomics` modules; scripts live in `scripts/eqtl/`, core logic in `gwas.finemapping.eqtl`.

## 6. Sample Size & Statistical Power Estimation

### GWAS Sample Size Calculator

For a binary trait with minor allele frequency (MAF) and desired power (80%) at α=5×10⁻⁸:

```
Required samples ≈ (Z_α/2 + Z_β)² × (1 + (1-MAF)²/MAF²) / (OR × log(OR))²
```

Practical rules:

| MAF | Detectable OR (80% power) | Required n |
|-----|---------------------------|------------|
| 0.40 | 1.10 | ~100,000 |
| 0.20 | 1.15 | ~60,000 |
| 0.10 | 1.25 | ~30,000 |
| 0.05 | 1.40 | ~15,000 |
| 0.01 | 2.00 | ~5,000 |

### Phenotype Effect Size Detection

For continuous trait with standard deviation σ:

| Desired detectable d | Required n/group | Total n |
|---------------------|------------------|---------|
| 0.3 (small) | 280 | 560 |
| 0.5 (medium) | 64 | 128 |
| 0.8 (large) | 21 | 42 |

### Multi-Omics Integration Sample Estimation

Because multi-omics reduces feature space via joint factorization:

| Model | Required samples | Notes |
|-------|-----------------|-------|
| Joint PCA | ≥ 5 × max(num_features) | Rule of thumb: 50–100 minimum |
| CCA | ≥ max(n_dna, n_rna) | Balanced sample sizes critical |
| NMF | ≥ 2 × (comp1 + comp2) | Less stringent but <100 risky |
| Multi-view clustering | ≥ 30 per view | Each omic layer needs samples |

## 7. Input Format Differences

### gwas: Standard Pipeline Input

```yaml
# config/gwas/gwas_template.yaml
work_dir: output/gwas/work
threads: 8

genome:
  accession: GCF_000001405.40
  dest_dir: output/gwas/genome

variants:
  vcf_files:
    - data/variants/cohort.vcf.gz

qc:
  min_maf: 0.01
  max_missing: 0.05
  hwe_pval: 1e-6

samples:
  phenotype_file: data/phenotypes.tsv  # columns: sample_id, phenotype, covariates
  covariates_file: data/covariates.tsv

association:
  model: linear    # or logistic
  trait: height
  covariates: [age, sex, PC1, PC2]

correction:
  method: bonferroni  # or fdr, genomic_control
```

**Phenotype file format** (`phenotypes.tsv`):
```
sample_id  height  bmi  disease_status
S001       172.3   22.1 0
S002       165.1   24.5 1
...
```

### phenotype: Measurement Table Format

```tsv
# measurements.tsv
specimen_id  measurement_type  value  unit  method  replicate  notes
sp001        head_length       1.85   mm    caliper  1          healthy
sp001        head_width        1.62   mm    caliper  1          healthy
sp001        thorax_length     2.10   mm    caliper  1          healthy
sp002        head_length       1.92   mm    caliper  1          healthy
...
```

**Metadata** (`metadata.json`):
```json
{
  "species": "Apis mellifera",
  "colony": "colony_A",
  "date": "2025-03-15",
  "observer": "researcher_1",
  "environment": "lab"
}
```

### multiomics: Container Format

```python-snippet
from metainformant.multiomics import integration

# Create container
multi = integration.MultiOmicsData()

# Add omics layers
dna_df = pd.read_csv("genotypes.tsv")      # samples × SNPs
rna_df = pd.read_csv("expression.tsv")    # samples × genes
prot_df = pd.read_csv("proteomics.tsv")   # samples × proteins

multi.add_omics("genome", dna_df, feature_type="snps")
multi.add_omics("transcriptome", rna_df, feature_type="genes")
multi.add_omics("proteome", prot_df, feature_type="proteins")

# Alignment happens automatically
assert multi.samples_aligned == True

# Run joint analysis
joint_pca = integration.joint_pca(multi, n_components=10)
cca_result = integration.cca_integration(multi, view1="genome", view2="transcriptome")
```

## 8. Output Interpretation

### gwas: Association Statistics

Output file: `output/gwas/results/association.tsv`

```
# Columns:
chrom   pos     snp_id   ref   alt   n_samples   beta   se     p_value   p_bonferroni
1       10583   rs12345  A     G     5000        0.12   0.02   1.2e-08   2.4e-02
...
```

**Interpretation**:
- `beta` > 0 → minor allele associated with increased trait
- `p_value < 5e-8` → genome-wide significant
- `p_bonferroni < 0.05` → after Bonferroni correction
- Manhattan/Q-Q plots visualize distribution

### phenotype: Trait Summaries

Output: `output/phenotype/measurements_summary.json`

```json
{
  "n_specimens": 250,
  "measurements": {
    "head_length": {"mean": 1.85, "sd": 0.12, "unit": "mm"},
    "head_width": {"mean": 1.62, "sd": 0.08, "unit": "mm"}
  },
  "indices": {
    "cephalic_index": {"mean": 87.6, "sd": 4.2}
  },
  "outliers": ["sp037", "sp089"]
}
```

**Interpretation**:
- Descriptive statistics (mean, SD, median, IQR)
- Outlier detection (specimens deviating >3 SD)
- Correlation matrix between traits

### multiomics: Integrated Results

Output: `output/multiomics/integration_results.h5`

**(HDF5 structure)**:
```
/
├─ joint_pca/
│  ├─ scores/ (samples × components)
│  ├─ loadings/ (features × components)
│  └─ variance_explained/ (per component)
├─ cca/
│  ├─ canonical_correlations/ (vector)
│  ├─ dna_loadings/ (SNPs × components)
│  └─ rna_loadings/ (genes × components)
├─ clusters/ (sample cluster assignments per omic)
└─ biomarkers/ (ranked cross-omic feature pairs)
```

**Interpretation**:
- `joint_pca.scores`: sample embeddings in shared space
- `cca.canonical_correlations`: strength of cross-omic relationship
- `biomarkers`: top integrated feature pairs (e.g., SNP+gene)

## 9. Performance Comparison

### Computational Complexity

| Operation | gwas | phenotype | multiomics |
|-----------|------|-----------|------------|
| **Association testing** | O(n_samples × n_snps) | O(n × p²) for correlations | O(n × (p1 + p2)) for joint PCA |
| **Memory footprint** | O(n_snps) for BED | O(n × p) measurements | O(n × (p1+p2+p3)) — largest |
| **I/O pattern** | Random access (plink) | Sequential read (CSV) | Block reads (HDF5) |
| **Bottleneck** | RAM (kinship matrix) | CPU (nested loops) | Memory (feature matrices) |

**Rough benchmarks** (10,000 samples):

| Task | gwas (linear model) | phenotype (correlation) | multiomics (joint PCA) |
|------|---------------------|------------------------|------------------------|
| **100k SNPs** | 2–5 min (8 threads) | N/A | N/A |
| **1000 measurements** | N/A | 1–2 min | N/A |
| **10k features × 3 omics** | N/A | N/A | 10–30 min (memory 32GB+) |

## 10. Example Commands / Code Snippets

### Equivalent: Running Each Module

#### GWAS: Association testing

```bash
# CLI
uv run python -m metainformant gwas run --config config/gwas/gwas_template.yaml

# Python API
from metainformant.gwas.workflow import run_gwas_workflow
results = run_gwas_workflow(config="config/gwas/gwas_template.yaml")
```

#### Phenotype: Trait analysis

```bash
uv run python -m metainformant phenotype analyze --measurements data/pheno.tsv
```

```python
from metainformant.phenotype import MorphometricProfile
profile = MorphometricProfile.from_tsv("measurements.tsv")
ci = profile.calculate_index("cephalic_index", "head_width", "head_length")
```

#### Multi-Omics: Joint integration

```bash
uv run python -m metainformant multiomics integrate \
  --dna data/genotypes.tsv \
  --rna data/expression.tsv \
  --method joint_pca
```

```python-snippet
from metainformant.multiomics import integration
multi = integration.MultiOmicsData()
multi.add_omics("dna", dna_df)
multi.add_omics("rna", rna_df)
result = integration.joint_pca(multi, n_components=20)
```

### Cross-Module Pipeline Example: GWAS → eQTL → Multi-omics

```python-snippet
# Step 1: Run GWAS
from metainformant.gwas.workflow import run_gwas_workflow
gwas_results = run_gwas_workflow(config="config/gwas/gwas.yaml")

# Step 2: Extract lead SNPs for eQTL
from metainformant.gwas.finemapping import credible_sets
credible = credible_sets.extract_credible_set(gwas_results, prob=0.95)

# Step 3: Run eQTL on those SNPs
from metainformant.gwas.finemapping.eqtl import run_eqtl
eqtl_results = run_eqtl(
    genotypes="cohort.vcf.gz",
    expression="expression.tsv",
    snps_to_test=credible.snp_ids
)

# Step 4: Integrate with proteomics (multiomics)
from metainformant.multiomics import integration
multi = integration.MultiOmicsData()
multi.add_omics("gwas", gwas_results[credible.snp_ids])
multi.add_omics("eqtl", eqtl_results)
multi.add_omics("proteome", protein_df)

# Joint factor analysis to find multi-omic axis
joint = integration.joint_pca(multi)
print(f"Top component explains {joint.variance[0]:.1%} integrated variance")
```

## 11. Common Pitfalls & Gotchas

| Pitfall | GWAS | Phenotype | Multiomics |
|---------|------|-----------|------------|
| **Small n, large p** | Yes — need n >> p for stable | Less critical | Severe — must regularize |
| **Batch effects** | Genotyping batch (fix with PC) | Measurement batch (critical) | Omics platform batch (critical) |
| **Population structure** | Major confounder (use PCA/kinship) | Minor (unless across pops) | Same as GWAS component |
| **Missing data** | Imputation possible | Often problematic | Hard — require complete cases or impute per-layer |
| **Multiple testing** | Genome-wide correction essential | FDR per-trait OK | Bonferroni or FDR across all features |
| **Reproducibility** | Summary statistics shareable | Raw measurements needed | Data volume often prohibitive to share |

## Decision Matrix: Which Tool First?

| Your Question | Start Here | Next Steps |
|---------------|------------|------------|
| "What genetic variants affect my trait?" | `gwas` | → If interested in mechanism: add eQTL/multiomics |
| "How do I characterize these phenotypes?" | `phenotype` | → If have genotypes: integrate with GWAS |
| "Can I find multi-omic biomarkers?" | `multiomics` | → Needs all omics layers collected |
| "Which gene mediates a GWAS hit?" | `eqtl` | ← Uses `gwas` results + `rna` expression |
| "Do these traits group into syndromes?" | `phenotype` | → Correlation + clustering |
| "What pathways are affected?" | `gwas` → `ontology` | Post-GWAS enrichment |

## 12. Special Considerations

### Ethical & Privacy

- **GWAS**: Human genetics requires IRB consent; summary statistics can often be shared, but individual genotypes cannot
- **Phenotype**: May include sensitive traits (behavioral, medical); de-identification essential
- **Multi-omics**: Highest sensitivity (genome + transcriptome + proteome + maybe metabolites); consortium-level data use agreements often required

### Cost Comparison (per sample, approximate 2026)

| Module | Genotyping | RNA-seq | Proteomics | Metabolomics | **Total** |
|--------|------------|---------|------------|--------------|-----------|
| gwas only | $50–150 | — | — | — | **$50–150** |
| phenotype only | — | — | — | — | **$100–500** (measurement labor/instrument) |
| rna only (bulk) | — | $50–200 | — | — | **$50–200** |
| singlecell | — | $200–1,000 | — | — | **$200–1,000** |
| multiomics (DNA+RNA) | $50–150 | $50–200 | — | — | **$100–350** |
| multiomics (DNA+RNA+Protein) | $50–150 | $50–200 | $200–500 | — | **$300–850** |
| multiomics (full) | $50–150 | $50–200 | $200–500 | $100–300 | **$400–1,150** |

### Reproducibility Checklist

**For GWAS**: ✓ Software version pinning (plink2, BOLT-LMM), ✓ Genome build (GRCh38/hg38), ✓ QC thresholds, ✓ Imputation panel (if used), ✓ Covariate list

**For Phenotype**: ✓ Measurement protocol, ✓ Calibration curves, ✓ Observer(s) documented, ✓ Raw data archived, ✓ Codebook

**For Multi-omics**: ✓ All of the above per omic layer, ✓ Batch correction parameters, ✓ Integration method hyperparameters, ✓ Sample matching log, ✓ Cross-modal validation split

## Related Documentation

- **[gwas/index.md](../gwas/index.md)** — Complete GWAS module guide
- **[phenotype/index.md](../phenotype/index.md)** — Phenotype analysis documentation
- **[multiomics/index.md](../multiomics/index.md)** — Multi-omic integration guide
- **[eqtl/README.md](../eqtl/README.md)** — eQTL integration pipeline (bridges GWAS+RNA)
- **[COMPARISON_GUIDES.md](../COMPARISON_GUIDES.md)** — Master comparison index
