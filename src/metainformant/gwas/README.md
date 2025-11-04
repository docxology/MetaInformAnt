# GWAS Module

Genome-Wide Association Studies (GWAS) module for identifying genetic variants associated with phenotypic traits.

## Overview

This module provides comprehensive GWAS functionality integrated with METAINFORMANT's bioinformatics ecosystem.

## Components

- **quality.py** - Quality control filters (MAF, missingness, HWE, quality scores)
- **structure.py** - Population structure analysis (PCA, kinship matrices)
- **association.py** - Association testing (linear/logistic regression)
- **correction.py** - Multiple testing correction (Bonferroni, FDR, genomic control)
- **visualization.py** - Result visualization (Manhattan plots, Q-Q plots)
- **calling.py** - Variant calling integration (bcftools, GATK) ✅ Integrated into workflow
- **download.py** - Data acquisition (reference genomes ✅, variant databases ⚠️ placeholder)
- **config.py** - Configuration management and validation
- **workflow.py** - Complete workflow orchestration

## Quick Start

```python
from metainformant.gwas import load_gwas_config, execute_gwas_workflow

# Load and execute workflow
config = load_gwas_config("config/gwas_config.yaml")
results = execute_gwas_workflow(config)
```

## CLI Usage

```bash
python -m metainformant gwas run --config config/gwas_config.yaml
```

## Documentation

Complete documentation available in [`docs/gwas/`](../../docs/gwas/):
- [Overview](../../docs/gwas/index.md)
- [Workflow Guide](../../docs/gwas/workflow.md)
- [Configuration Reference](../../docs/gwas/config.md)

## Integration with Other Modules

### With DNA Variants Module
```python
from metainformant.dna import variants
from metainformant.gwas import parse_vcf_full, apply_qc_filters

# VCF parsing and analysis workflow
# Parse VCF file containing variant data
vcf_data = parse_vcf_full("variants.vcf")

# Use DNA module for additional variant operations
vcf_info = variants.parse_vcf("variants.vcf")

# Apply quality control filters
qc_variants = apply_qc_filters(
    vcf_data,
    maf_threshold=0.05,
    missing_threshold=0.1,
    hwe_threshold=1e-6
)

# Extract variant information for association testing
genotype_matrix = extract_genotypes(qc_variants)
variant_positions = extract_positions(qc_variants)
```

### With DNA Population Module
```python
from metainformant.dna import population
from metainformant.gwas import compute_pca, compute_kinship

# Population genetics integration
# Load sequences for population analysis
sequences = population.read_sequences("populations.fasta")

# Calculate population diversity metrics
nuc_diversity = population.nucleotide_diversity(sequences)
tajima_d = population.tajimas_D(sequences)

# Use population structure for GWAS
# PCA for population stratification
pca_result = compute_pca(genotype_matrix, n_components=10)
pca_covariates = pca_result["components"]

# Kinship matrix for relatedness correction
kinship_matrix = compute_kinship(genotype_matrix)

# Include population genetics metrics in association testing
association_results = association_test_linear(
    genotype_matrix,
    phenotypes,
    covariates=pca_covariates,  # Population structure covariates
    kinship=kinship_matrix  # Relatedness matrix
)
```

### With Math Popgen Module
```python
from metainformant.math import popgen_stats
from metainformant.gwas import association_test_linear

# Theoretical population genetics calculations
# Calculate expected allele frequencies
allele_freqs = calculate_allele_frequencies(genotype_matrix)

# Use theoretical models for interpretation
# Expected heterozygosity
expected_het = popgen_stats.expected_heterozygosity(allele_freqs)

# Hardy-Weinberg equilibrium
hwe_pvalues = popgen_stats.hardy_weinberg_test(genotype_matrix)

# Fst between populations
fst_values = popgen_stats.fst_calculation(pop1_genotypes, pop2_genotypes)

# Use theoretical expectations for filtering
# Filter variants deviating from HWE
hwe_filtered = genotype_matrix[hwe_pvalues > 1e-6]

# Association testing with theoretical context
association_results = association_test_linear(
    hwe_filtered,
    phenotypes,
    covariates=None
)
```

### With ML Regression Module
```python
from metainformant.ml import BiologicalRegressor, cross_validate_biological
from metainformant.gwas import association_test_linear, association_test_logistic

# Association testing using ML regression models
# Linear regression for continuous traits
linear_results = association_test_linear(
    genotype_matrix,
    phenotypes,
    covariates=covariates
)

# Logistic regression for binary traits
logistic_results = association_test_logistic(
    genotype_matrix,
    binary_phenotypes,
    covariates=covariates
)

# Use ML module for advanced regression
regressor = BiologicalRegressor(algorithm="lasso", alpha=0.01, random_state=42)
# Fit regressor for each variant
for variant_idx in range(genotype_matrix.shape[1]):
    variant_genotypes = genotype_matrix[:, variant_idx]
    regressor.fit(variant_genotypes.reshape(-1, 1), phenotypes)
    predictions = regressor.predict(variant_genotypes.reshape(-1, 1))

# Cross-validation for model assessment
cv_results = cross_validate_biological(
    regressor,
    genotype_matrix[:, 0:1],  # Single variant
    phenotypes,
    cv=5
)
```

### With Visualization Module
```python
from metainformant.gwas import association_test_linear, parse_vcf_full
from metainformant.visualization import manhattan_plot, qq_plot, regional_plot

# Comprehensive GWAS visualization
# Run association testing
association_results = association_test_linear(
    genotype_matrix,
    phenotypes,
    covariates=covariates
)

# Manhattan plot for genome-wide results
ax = manhattan_plot(
    association_results["chromosome"],
    association_results["position"],
    association_results["p_value"],
    title="GWAS Manhattan Plot",
    xlabel="Chromosome",
    ylabel="-log10(p-value)"
)

# Q-Q plot for p-value distribution
ax = qq_plot(
    association_results["p_value"],
    title="GWAS Q-Q Plot",
    xlabel="Expected -log10(p)",
    ylabel="Observed -log10(p)"
)

# Regional plot for specific locus
significant_variant = association_results.loc[
    association_results["p_value"].idxmin()
]
ax = regional_plot(
    chromosome=significant_variant["chromosome"],
    start=significant_variant["position"] - 500000,
    end=significant_variant["position"] + 500000,
    association_results=association_results,
    title="Regional Association Plot"
)

# Effect size visualization
ax = scatter_plot(
    association_results["beta"],
    -np.log10(association_results["p_value"]),
    xlabel="Effect Size (β)",
    ylabel="-log10(p-value)",
    title="Effect Size vs Significance"
)
```

### With Quality Module
```python
from metainformant.gwas import apply_qc_filters, parse_vcf_full
from metainformant.quality import analyze_fastq_quality
from metainformant.dna.fastq import iter_fastq

# Comprehensive QC workflow for GWAS
# Quality assessment before variant calling
genomic_reads = list(iter_fastq("genomic_reads.fastq"))
quality_stats = analyze_fastq_quality(genomic_reads)

# Parse VCF after variant calling
vcf_data = parse_vcf_full("variants.vcf")

# Apply GWAS-specific QC filters
qc_variants = apply_qc_filters(
    vcf_data,
    maf_threshold=0.05,  # Minor allele frequency
    missing_threshold=0.1,  # Missing data threshold
    hwe_threshold=1e-6  # Hardy-Weinberg equilibrium
)

# Generate QC report
qc_report = {
    "initial_variants": len(vcf_data),
    "after_maf_filter": len(qc_variants[qc_variants["maf"] >= 0.05]),
    "after_missing_filter": len(qc_variants[qc_variants["missing_rate"] <= 0.1]),
    "after_hwe_filter": len(qc_variants[qc_variants["hwe_pvalue"] > 1e-6])
}

# Use quality-filtered variants for association testing
genotype_matrix = extract_genotypes(qc_variants)
association_results = association_test_linear(
    genotype_matrix,
    phenotypes,
    covariates=covariates
)
```

## Scientific Methods

Based on established GWAS methods:
- Price et al. (2006) - PCA for population stratification
- VanRaden (2008) - Kinship matrix computation
- Benjamini & Hochberg (1995) - FDR correction
- Yang et al. (2011) - GCTA methods

## Related Modules

- **DNA**: Variant parsing and population genetics
- **Math**: Theoretical population genetics
- **ML**: Regression and classification models
- **Visualization**: Manhattan and Q-Q plots

