# GWAS Visualization Gallery

This document showcases all available visualization types in the METAINFORMANT GWAS module, organized by category.

## Quick Start

```python
from metainformant.gwas import generate_all_plots

# Generate all visualizations at once
results = generate_all_plots(
    association_results="output/gwas/results/association_results.tsv",
    output_dir="output/gwas/plots",
    pca_file="output/gwas/structure/pca_components.tsv",
    kinship_file="output/gwas/structure/kinship_matrix.tsv",
    significance_threshold=5e-8
)

print(f"Generated {results['num_plots_generated']} plots")
```

## Visualization Categories

### 1. Genome-Wide Visualizations

#### Manhattan Plot
**Module**: `visualization_genome.py`  
**Function**: `manhattan_plot()`

Standard chromosome-wide association plot optimized for millions of SNPs.

**Features**:
- Intelligent point thinning for performance
- Alternating chromosome colors
- Genome-wide and suggestive significance lines
- Handles up to 50,000 points per chromosome

**Use Cases**:
- Primary GWAS results visualization
- Identifying genome-wide significant loci
- Publication figures

**Example**:
```python
from metainformant.gwas.visualization_genome import manhattan_plot

result = manhattan_plot(
    results="results.tsv",
    output_path="plots/manhattan.png",
    significance_threshold=5e-8,
    max_points_per_chrom=50000
)
```

#### Circular Manhattan Plot
**Function**: `circular_manhattan_plot()`

Circular genome-wide view with radial p-values.

**Use Cases**:
- Compact genome-wide overview
- Identifying genome-wide patterns
- Alternative to linear Manhattan plot

#### Chromosome Ideogram
**Function**: `chromosome_ideogram()`

Chromosome map with significant loci marked.

**Use Cases**:
- Visualizing distribution of hits
- Quick overview of significant regions
- Identifying clustered associations

---

### 2. Statistical Diagnostics

#### Q-Q Plot
**Module**: `visualization_statistical.py`  
**Function**: `qq_plot()`

Quantile-quantile plot with 95% CI and genomic inflation factor.

**Features**:
- 95% confidence intervals
- Lambda GC calculation and display
- Automatic inflation status assessment

**Interpretation**:
- Points on diagonal: well-calibrated
- Deviation at right: true associations
- Global inflation: population stratification or relatedness

**Example**:
```python
from metainformant.gwas.visualization_statistical import qq_plot

result = qq_plot(
    results="results.tsv",
    output_path="plots/qq_plot.png",
    show_ci=True,
    show_lambda_gc=True
)
```

#### Stratified Q-Q Plot
**Function**: `qq_plot_stratified()`

Q-Q plot separated by MAF bins.

**Use Cases**:
- Checking if inflation varies by frequency
- Identifying frequency-specific issues
- Quality control by allele frequency

#### Lambda GC by Chromosome
**Function**: `lambda_gc_plot()`

Genomic inflation factor per chromosome.

**Use Cases**:
- Identifying chromosome-specific issues
- Quality control across genome
- Detecting localized problems

#### Volcano Plot
**Function**: `volcano_plot()`

Effect size vs significance scatter plot.

**Use Cases**:
- Identifying large-effect variants
- Combining significance with biological relevance
- Publication figures

#### Power Plot
**Function**: `power_plot()`

Statistical power curves at different sample sizes.

**Use Cases**:
- Study design and planning
- Post-hoc power analysis
- Sample size estimation

---

### 3. Regional/Locus-Specific

#### Regional Association Plot
**Module**: `visualization_regional.py`  
**Function**: `regional_plot()`

Detailed association view for genomic window.

**Use Cases**:
- Fine-mapping significant loci
- Identifying lead SNPs
- Examining local LD structure

**Example**:
```python
from metainformant.gwas.visualization_regional import regional_plot

result = regional_plot(
    results="results.tsv",
    output_path="plots/region_chr1.png",
    chrom="chr1",
    start=1000000,
    end=2000000,
    lead_snp_pos=1500000
)
```

#### Regional LD Plot
**Function**: `regional_ld_plot()`

Linkage disequilibrium around lead SNP.

**Note**: Requires external LD calculation (PLINK recommended).

#### Gene Annotation Plot
**Function**: `gene_annotation_plot()`

Association with gene structures overlaid.

**Use Cases**:
- Identifying variants in genes
- Functional annotation
- Publication figures

#### Recombination Rate Plot
**Function**: `recombination_rate_plot()`

Association with recombination hotspots.

**Note**: Requires recombination map data.

---

### 4. Population Structure

#### PCA Plot
**Module**: `visualization_population.py`  
**Function**: `pca_plot()`

Principal component scatter plot (2D or 3D).

**Features**:
- Optional phenotype coloring
- 2D and 3D modes
- Sample outlier identification

**Use Cases**:
- Visualizing population stratification
- Quality control (outlier detection)
- Adjusting for ancestry

**Example**:
```python
from metainformant.gwas.visualization_population import pca_plot

result = pca_plot(
    pca_file="structure/pca_components.tsv",
    output_path="plots/pca.png",
    pc1=1,
    pc2=2,
    phenotype_file="phenotypes.tsv"
)
```

#### PCA Scree Plot
**Function**: `pca_scree_plot()`

Variance explained by each PC.

**Use Cases**:
- Determining number of PCs to adjust for
- Understanding population structure complexity

#### Kinship Heatmap
**Function**: `kinship_heatmap()`

Pairwise genetic relationships.

**Use Cases**:
- Identifying related samples
- Mixed model adjustment
- Quality control

#### Admixture Plot
**Function**: `admixture_plot()`

Ancestry proportions (stacked bar plot).

**Note**: Requires ADMIXTURE software output.

#### Population Tree
**Function**: `population_tree()`

Hierarchical clustering dendrogram.

**Use Cases**:
- Visualizing population relationships
- Identifying genetic clusters

---

### 5. Variant Properties

#### MAF Distribution
**Module**: `visualization_variants.py`  
**Function**: `maf_distribution()`

Allele frequency spectrum histogram.

**Features**:
- Rare variant highlighting (MAF < 0.05)
- Summary statistics

**Use Cases**:
- Quality control
- Understanding variant frequency distribution
- Assessing QC filter effects

**Example**:
```python
from metainformant.gwas.visualization_variants import maf_distribution

result = maf_distribution(
    results="results.tsv",
    output_path="plots/maf_dist.png",
    bins=50
)
```

#### Variant Density Plot
**Function**: `variant_density_plot()`

SNP density across genome.

**Use Cases**:
- Identifying high/low density regions
- Quality control
- Visualizing variant distribution

#### HWE Deviation Plot
**Function**: `hwe_deviation_plot()`

Hardy-Weinberg equilibrium p-value distribution.

**Use Cases**:
- Quality control
- Identifying genotyping errors
- Detecting population substructure

#### Missingness Plot
**Function**: `missingness_plot()`

Missing data patterns.

**Note**: Requires VCF parsing implementation.

#### Transition/Transversion Ratio
**Function**: `transition_transversion_plot()`

Ts/Tv ratio visualization and interpretation.

**Interpretation**:
- Expected Ts/Tv: 2.0-2.1 (whole genome)
- Low Ts/Tv: possible quality issues
- Automated quality assessment

---

### 6. Effect Sizes and Functional Annotation

#### Forest Plot
**Module**: `visualization_effects.py`  
**Function**: `effect_size_forest_plot()`

Effect estimates with 95% confidence intervals.

**Features**:
- Top N significant variants
- Confidence intervals
- Directional coloring

**Use Cases**:
- Reporting significant associations
- Comparing effect sizes
- Publication figures

**Example**:
```python
from metainformant.gwas.visualization_effects import effect_size_forest_plot

result = effect_size_forest_plot(
    results="results.tsv",
    output_path="plots/forest.png",
    top_n=20,
    significance_threshold=5e-8
)
```

#### Effect Direction Plot
**Function**: `effect_direction_plot()`

Distribution of positive/negative effects.

**Use Cases**:
- Checking for directional bias
- Quality control
- Understanding effect patterns

#### Functional Enrichment Plot
**Function**: `functional_enrichment_plot()`

Variant distribution by functional category.

**Note**: Requires variant annotation (ANNOVAR, VEP, SnpEff).

#### Allelic Series Plot
**Function**: `allelic_series_plot()`

Multiple alleles at locus with dose-response.

**Note**: For complex loci with multiple functional alleles.

---

### 7. Multi-Trait and Cross-Cohort Comparisons

#### Miami Plot
**Module**: `visualization_comparison.py`  
**Function**: `miami_plot()`

Back-to-back Manhattan plots for two traits.

**Use Cases**:
- Identifying pleiotropic loci
- Comparing two phenotypes
- Publication figures for trait pairs

**Example**:
```python
from metainformant.gwas.visualization_comparison import miami_plot

result = miami_plot(
    results1="trait1_results.tsv",
    results2="trait2_results.tsv",
    output_path="plots/miami.png",
    trait1_name="Height",
    trait2_name="BMI"
)
```

#### Multi-Trait Manhattan
**Function**: `multi_trait_manhattan()`

Stacked Manhattan plots for multiple phenotypes.

**Use Cases**:
- Comparing many traits simultaneously
- Identifying shared associations
- Comprehensive trait analysis

#### Cross-Cohort Forest Plot
**Function**: `cross_cohort_forest()`

Meta-analysis forest plot.

**Note**: For comparing variant effects across studies.

#### Concordance Plot
**Function**: `concordance_plot()`

Discovery vs replication scatter plot.

**Features**:
- Correlation calculation
- Best-fit line
- Concordance statistics

**Use Cases**:
- Validating replication
- Meta-analysis quality control
- Assessing consistency

---

## Comprehensive Visualization Suite

Generate all plots at once:

```python
from metainformant.gwas import generate_all_plots

results = generate_all_plots(
    association_results="output/gwas/results/association_results.tsv",
    output_dir="output/gwas/plots",
    pca_file="output/gwas/structure/pca_components.tsv",
    kinship_file="output/gwas/structure/kinship_matrix.tsv",
    vcf_file="data/variants.vcf.gz",
    significance_threshold=5e-8
)

# Check results
print(f"Successful: {results['num_plots_generated']}")
print(f"Skipped: {results['num_plots_skipped']}")
print(f"Errors: {results['num_errors']}")

# List generated plots
for plot_name in results['successful_plots']:
    print(f"  ✓ {plot_name}")
```

## Workflow Integration

Enable comprehensive plots in workflow:

```yaml
# config/gwas/gwas_amellifera.yaml
output:
  results_dir: output/gwas/amellifera/results
  plots_dir: output/gwas/amellifera/plots
  comprehensive_plots: true  # Enable all visualizations
  significance_threshold: 5e-8
```

Or programmatically:

```python
from metainformant.gwas import execute_gwas_workflow, load_gwas_config

config = load_gwas_config("config/gwas/gwas_amellifera.yaml")
config.output["comprehensive_plots"] = True

results = execute_gwas_workflow(config)
```

## Performance Notes

### Large-Scale Data

For genome-scale data (millions of SNPs):

1. **Manhattan Plot**: Uses intelligent thinning (keeps all significant, samples non-significant)
2. **Q-Q Plot**: Handles millions of points efficiently
3. **Circular Manhattan**: May be slow with >5M SNPs
4. **Kinship Heatmap**: Auto-downsamples if >200 samples

### Optimization Tips

- Use `max_points_per_chrom` parameter for Manhattan plots
- Generate regional plots only for significant loci
- Use summary statistics for large cohorts
- Consider rasterization for publication figures

## Customization

All visualization functions accept standard matplotlib parameters:

```python
# Custom colors
manhattan_plot(
    results="results.tsv",
    output_path="plots/manhattan.png",
    chrom_colors=["#FF6B6B", "#4ECDC4", "#45B7D1"]
)

# Custom size
qq_plot(
    results="results.tsv",
    output_path="plots/qq.png",
    title="My Custom QQ Plot"
)
```

## Exporting for Publication

All plots support high-resolution export:

- Default DPI: 300 (publication quality)
- Formats: PNG (default), SVG, PDF
- Vector graphics for scalability

```python
# Save as vector graphic
import matplotlib.pyplot as plt

# After generating plot, save as PDF
plt.savefig("figure.pdf", dpi=300, bbox_inches="tight")
```

## Troubleshooting

### Common Issues

**Issue**: "No variants found"
- Check input file format (TSV with headers)
- Ensure CHROM, POS, p_value columns exist

**Issue**: "Plot generation skipped"
- Some plots require additional data (PCA, kinship, VCF)
- Check log messages for specific requirements

**Issue**: "Memory error with large datasets"
- Reduce `max_points_per_chrom` for Manhattan plots
- Use regional plots instead of genome-wide
- Increase system memory or use HPC

### Getting Help

For issues or questions:
- Check log files for detailed error messages
- See individual function docstrings for parameters
- Refer to examples in `tests/gwas/` directory

## Future Enhancements

Planned visualization additions:
- Interactive plots (Plotly/Bokeh)
- Gene track annotations
- LD heatmaps
- Conditional analysis plots
- PheWAS plots
- GWAS catalog integration

---

## Summary

**Total Visualization Functions**: 30+

**Categories**:
- Genome-wide: 4 functions
- Statistical: 5 functions
- Regional: 4 functions
- Population: 5 functions
- Variants: 5 functions
- Effects: 4 functions
- Comparison: 4 functions

**Quick Reference**: Use `generate_all_plots()` to create comprehensive visualization suite automatically.



