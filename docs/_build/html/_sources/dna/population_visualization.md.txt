# DNA: Population Genetics Visualization

Comprehensive visualization functions for population genetics analysis.

**Functions**: `plot_diversity_comparison`, `plot_tajimas_d_comparison`, `plot_fst_comparison`, `plot_pca_results`, `plot_kinship_matrix`, `plot_site_frequency_spectrum`, `plot_neutrality_test_summary`, `plot_demographic_comparison`, `plot_summary_statistics_grid`, `plot_linkage_disequilibrium_decay`, `plot_allele_frequency_spectrum`, `plot_pairwise_distance_distribution`, `plot_heterozygosity_distribution`, `plot_statistic_distribution`, `plot_pi_vs_theta`, `plot_statistic_correlation_matrix`, `plot_regression_analysis`, `plot_statistic_relationships`, `plot_fst_matrix`, `plot_admixture_plot`, `plot_population_structure_tree`, `plot_isolation_by_distance`, `plot_three_population_f3`, `plot_neutrality_test_suite`, `plot_hardy_weinberg_test`, `plot_bootstrap_distribution`, `plot_permutation_test`, `plot_outlier_detection`.

## Overview

The `dna.population_viz` module provides publication-quality visualizations for population genetics analysis, including diversity metrics, neutrality tests, population structure, and demographic models.

## Functions

### `plot_diversity_comparison(diversity_values, ...)`

Plot comparison of nucleotide diversity (π) across scenarios.

**Parameters**:
- `diversity_values`: Dictionary mapping scenario names to π values
- `output_path`: Optional path to save figure
- `title`: Plot title
- `figsize`: Figure size

**Returns**: matplotlib Figure object

**Example**:
```python
from metainformant.dna.population_viz import plot_diversity_comparison

diversity = {
    "Neutral": 0.01,
    "High Diversity": 0.05,
    "Low Diversity": 0.001,
}

plot_diversity_comparison(
    diversity,
    output_path="output/diversity_comparison.png",
    title="Nucleotide Diversity Comparison"
)
```

### `plot_tajimas_d_comparison(tajimas_d_values, ...)`

Plot comparison of Tajima's D across scenarios with color coding.

**Parameters**:
- `tajimas_d_values`: Dictionary mapping scenario names to Tajima's D values
- `output_path`: Optional path to save figure
- `title`: Plot title
- `figsize`: Figure size

**Returns**: matplotlib Figure object

**Features**:
- Red bars for negative D (population expansion)
- Blue bars for positive D (balancing selection)
- Reference lines at D = -2, 0, 2

**Example**:
```python
from metainformant.dna.population_viz import plot_tajimas_d_comparison

tajimas_d = {
    "Neutral": 0.1,
    "Bottleneck": -1.5,
    "Expansion": -1.2,
}

plot_tajimas_d_comparison(
    tajimas_d,
    output_path="output/tajimas_d.png"
)
```

### `plot_fst_comparison(fst_values, ...)`

Plot comparison of Fst values with interpretation thresholds.

**Parameters**:
- `fst_values`: Dictionary mapping comparison names to Fst values
- `output_path`: Optional path to save figure
- `title`: Plot title
- `figsize`: Figure size

**Returns**: matplotlib Figure object

**Features**:
- Reference lines for Fst interpretation (0.05 = low, 0.15 = moderate, 0.25 = high)

**Example**:
```python
from metainformant.dna.population_viz import plot_fst_comparison

fst_values = {
    "Pop1 vs Pop2": 0.2,
    "Pop1 vs Pop3": 0.05,
    "Pop2 vs Pop3": 0.3,
}

plot_fst_comparison(fst_values, output_path="output/fst.png")
```

### `plot_pca_results(pca_result, ...)`

Plot PCA results with multiple views.

**Parameters**:
- `pca_result`: Dictionary with PCA results (from `compute_pca`)
- `output_path`: Optional path to save figure
- `n_components`: Number of components to plot
- `title`: Plot title
- `figsize`: Figure size

**Returns**: matplotlib Figure object

**Plots include**:
1. PC1 vs PC2 scatter plot
2. Explained variance per component
3. Cumulative explained variance

**Example**:
```python
from metainformant.dna.population_viz import plot_pca_results
from metainformant.gwas.structure import compute_pca

pca_result = compute_pca(genotype_matrix, n_components=10)
plot_pca_results(
    pca_result,
    output_path="output/pca.png",
    n_components=10
)
```

### `plot_kinship_matrix(kinship_result, ...)`

Plot kinship matrix as heatmap.

**Parameters**:
- `kinship_result`: Dictionary with kinship results (from `compute_kinship_matrix`)
- `output_path`: Optional path to save figure
- `max_samples`: Maximum number of samples to plot (for large matrices)
- `title`: Plot title
- `figsize`: Figure size

**Returns**: matplotlib Figure object

**Example**:
```python
from metainformant.dna.population_viz import plot_kinship_matrix
from metainformant.gwas.structure import compute_kinship_matrix

kinship_result = compute_kinship_matrix(genotype_matrix, method="vanraden")
plot_kinship_matrix(
    kinship_result,
    output_path="output/kinship.png",
    max_samples=100
)
```

### `plot_site_frequency_spectrum(sfs, ...)`

Plot site frequency spectrum.

**Parameters**:
- `sfs`: Site frequency spectrum (list of counts per frequency bin)
- `output_path`: Optional path to save figure
- `title`: Plot title
- `figsize`: Figure size

**Returns**: matplotlib Figure object

**Example**:
```python
from metainformant.dna.population_viz import plot_site_frequency_spectrum
from metainformant.simulation.popgen import generate_site_frequency_spectrum

sfs = generate_site_frequency_spectrum(
    sample_size=30,
    n_sites=100,
    theta=0.01
)
plot_site_frequency_spectrum(
    sfs,
    output_path="output/sfs.png"
)
```

### `plot_neutrality_test_summary(neutrality_results, ...)`

Plot comprehensive summary of neutrality tests across scenarios.

**Parameters**:
- `neutrality_results`: Dictionary mapping scenario names to neutrality test results
- `output_path`: Optional path to save figure
- `title`: Plot title
- `figsize`: Figure size

**Returns**: matplotlib Figure object

**Plots include**:
1. Tajima's D by scenario
2. π/θ ratio by scenario
3. Nucleotide diversity by scenario
4. Segregating sites by scenario

**Example**:
```python
from metainformant.dna.population_viz import plot_neutrality_test_summary
from metainformant.dna.population_analysis import neutrality_test_suite

neutrality_data = {
    "Neutral": neutrality_test_suite(neutral_seqs),
    "Bottleneck": neutrality_test_suite(bottleneck_seqs),
    "Expansion": neutrality_test_suite(expansion_seqs),
}

plot_neutrality_test_summary(
    neutrality_data,
    output_path="output/neutrality_summary.png"
)
```

### `plot_demographic_comparison(demographic_results, ...)`

Plot comparison of demographic models.

**Parameters**:
- `demographic_results`: Dictionary with demographic model results
- `output_path`: Optional path to save figure
- `title`: Plot title
- `figsize`: Figure size

**Returns**: matplotlib Figure object

**Plots include**:
1. Estimated effective population size (Ne)
2. Observed diversity (π)

**Example**:
```python
from metainformant.dna.population_viz import plot_demographic_comparison

demographic_results = {
    "bottleneck": {
        "estimated_ne": 15.0,
        "observed_diversity": 0.025,
    },
    "expansion": {
        "estimated_ne": 250.0,
        "observed_diversity": 0.005,
    },
}

plot_demographic_comparison(
    demographic_results,
    output_path="output/demographic_comparison.png"
)
```

### `plot_summary_statistics_grid(summary_stats, ...)`

Plot grid of summary statistics across scenarios.

**Parameters**:
- `summary_stats`: Dictionary mapping scenario names to summary statistics
- `output_path`: Optional path to save figure
- `title`: Plot title
- `figsize`: Figure size

**Returns**: matplotlib Figure object

**Plots include**:
1. Nucleotide diversity (π)
2. Segregating sites (S)
3. Watterson's theta (θ_W)
4. Tajima's D
5. Sample size
6. Sequence length

**Example**:
```python
from metainformant.dna.population_viz import plot_summary_statistics_grid
from metainformant.dna.population_analysis import calculate_summary_statistics

summary_stats = {
    "Neutral": calculate_summary_statistics(sequences=neutral_seqs),
    "High Diversity": calculate_summary_statistics(sequences=high_div_seqs),
    "Low Diversity": calculate_summary_statistics(sequences=low_div_seqs),
}

plot_summary_statistics_grid(
    summary_stats,
    output_path="output/summary_grid.png"
)
```

### `plot_linkage_disequilibrium_decay(ld_values, distances, ...)`

Plot linkage disequilibrium decay with distance.

**Parameters**:
- `ld_values`: r² values at different distances
- `distances`: Optional distances (if None, uses indices)
- `output_path`: Optional path to save figure
- `title`: Plot title
- `figsize`: Figure size

**Returns**: matplotlib Figure object

**Example**:
```python
from metainformant.dna.population_viz import plot_linkage_disequilibrium_decay

r_squared_values = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05]
distances = [0, 1, 2, 3, 4, 5]  # In units (e.g., kb)

plot_linkage_disequilibrium_decay(
    r_squared_values,
    distances=distances,
    output_path="output/ld_decay.png"
)
```

## Integration Examples

### Complete Workflow with Visualizations

```python
from metainformant.dna.population_analysis import (
    calculate_summary_statistics,
    compare_populations,
    neutrality_test_suite,
)
from metainformant.dna.population_viz import (
    plot_diversity_comparison,
    plot_neutrality_test_summary,
    plot_summary_statistics_grid,
)

# Analyze data
neutral_stats = calculate_summary_statistics(sequences=neutral_seqs)
neutrality = neutrality_test_suite(neutral_seqs)

# Create visualizations
plot_diversity_comparison(
    {"Neutral": neutral_stats["nucleotide_diversity"]},
    output_path="output/diversity.png"
)

plot_neutrality_test_summary(
    {"Neutral": neutrality},
    output_path="output/neutrality.png"
)
```

### Batch Visualization Generation

```python
from metainformant.dna.population_viz import (
    plot_diversity_comparison,
    plot_tajimas_d_comparison,
    plot_fst_comparison,
)

# Generate all comparison plots
scenarios = ["neutral", "high_diversity", "low_diversity", "bottleneck", "expansion"]
diversity_values = {s: get_diversity(s) for s in scenarios}
tajimas_d_values = {s: get_tajimas_d(s) for s in scenarios}

plot_diversity_comparison(
    diversity_values,
    output_path="output/diversity_comparison.png"
)

plot_tajimas_d_comparison(
    tajimas_d_values,
    output_path="output/tajimas_d_comparison.png"
)
```

## Best Practices

### 1. High-Resolution Output

All functions support saving at 300 DPI for publication:

```python
plot_diversity_comparison(
    diversity_values,
    output_path="output/diversity.png"  # Saved at 300 DPI
)
```

### 2. Custom Styling

Functions return matplotlib Figure objects for further customization:

```python
fig = plot_diversity_comparison(diversity_values)
ax = fig.axes[0]
ax.set_ylabel("Diversity (π)", fontsize=14)
fig.savefig("output/custom.png", dpi=300)
```

### 3. Combining Multiple Plots

Use matplotlib subplots to combine visualizations:

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Plot 1: Diversity
ax1 = axes[0]
# ... customization

# Plot 2: Tajima's D
ax2 = axes[1]
# ... customization

plt.tight_layout()
fig.savefig("output/combined.png", dpi=300)
```

## See Also

- [`dna.population`](population.md) - Core population genetics functions
- [`dna.population_analysis`](population_workflows.md) - Analysis workflows
- [`visualization.plots`](../visualization/plots.md) - General plotting utilities


## Additional Visualization Functions

### Distribution Plots

#### `plot_allele_frequency_spectrum(sfs, ...)`

Plot allele frequency spectrum (SFS) showing the distribution of allele frequencies.

**Parameters**:
- `sfs`: Site frequency spectrum (counts of sites at each frequency)
- `unfolded`: Whether SFS is unfolded (requires ancestral state)
- `output_path`: Optional path to save figure

#### `plot_pairwise_distance_distribution(pairwise_distances, ...)`

Plot distribution of pairwise nucleotide differences.

**Parameters**:
- `pairwise_distances`: Pairwise distance values
- `bins`: Number of histogram bins (default: 30)

#### `plot_heterozygosity_distribution(heterozygosity_values, ...)`

Plot distribution of observed heterozygosity across sites.

**Parameters**:
- `heterozygosity_values`: Per-site heterozygosity values
- `bins`: Number of histogram bins (default: 30)

#### `plot_statistic_distribution(statistic_values, ...)`

Plot distribution of statistics across scenarios with multiple plot types.

**Parameters**:
- `statistic_values`: Dictionary mapping scenario names to statistic values
- `plot_type`: Type of plot ("histogram", "violin", "box")

### Correlation and Relationship Plots

#### `plot_pi_vs_theta(pi_values, theta_values, ...)`

Plot scatter of π vs θ_W with expected 1:1 line and regression analysis.

**Parameters**:
- `pi_values`: Nucleotide diversity (π) values
- `theta_values`: Watterson's theta (θ_W) values

#### `plot_statistic_correlation_matrix(statistics, ...)`

Plot correlation heatmap between statistics with hierarchical clustering.

**Parameters**:
- `statistics`: Dictionary mapping statistic names to values

#### `plot_regression_analysis(x_values, y_values, ...)`

Plot regression analysis with confidence bands and residual analysis.

**Parameters**:
- `x_values`: Independent variable values
- `y_values`: Dependent variable values
- `x_label`: Label for x-axis
- `y_label`: Label for y-axis

#### `plot_statistic_relationships(statistics, ...)`

Plot multi-panel relationship plots between statistics.

**Parameters**:
- `statistics`: Dictionary mapping statistic names to values

### Multi-Population Visualizations

#### `plot_fst_matrix(fst_matrix, population_names, ...)`

Plot heatmap of pairwise Fst values with hierarchical clustering dendrogram.

**Parameters**:
- `fst_matrix`: 2D array or nested dict of Fst values
- `population_names`: Names of populations

#### `plot_admixture_plot(ancestry_proportions, ...)`

Plot ancestry proportions as stacked bar plot or heatmap.

**Parameters**:
- `ancestry_proportions`: List of ancestry proportion lists per individual
- `population_names`: Names of populations/components

#### `plot_population_structure_tree(fst_matrix, ...)`

Plot neighbor-joining or UPGMA tree based on Fst or genetic distances.

**Parameters**:
- `fst_matrix`: Fst matrix or distance matrix
- `population_names`: Names of populations

#### `plot_isolation_by_distance(fst_values, geographic_distances, ...)`

Plot Fst vs geographic distance with Mantel test visualization.

**Parameters**:
- `fst_values`: Pairwise Fst values
- `geographic_distances`: Corresponding geographic distances

#### `plot_three_population_f3(f3_values, ...)`

Plot F3-statistic for admixture testing.

**Parameters**:
- `f3_values`: Dictionary of F3 values (key format: "Pop1,Pop2,Pop3")

### Statistical Test Visualizations

#### `plot_neutrality_test_suite(test_results, ...)`

Plot comprehensive panel of all neutrality tests with significance indicators.

**Parameters**:
- `test_results`: Dictionary mapping test names to result dictionaries

#### `plot_hardy_weinberg_test(hwe_results, ...)`

Plot Hardy-Weinberg test results with chi-square and p-value distributions.

**Parameters**:
- `hwe_results`: HWE test results (single dict or list of dicts)

#### `plot_bootstrap_distribution(bootstrap_values, observed_value, ...)`

Plot bootstrap distribution with confidence interval visualization.

**Parameters**:
- `bootstrap_values`: Bootstrap replicate values
- `observed_value`: Observed statistic value
- `confidence_level`: Confidence level for CI

#### `plot_permutation_test(permuted_values, observed_value, ...)`

Plot permuted vs observed distribution with p-value calculation.

**Parameters**:
- `permuted_values`: Permuted statistic values
- `observed_value`: Observed statistic value
- `p_value`: Optional p-value to display

#### `plot_outlier_detection(statistic_values, outlier_indices, ...)`

Plot Manhattan plot of statistics with outliers highlighted.

**Parameters**:
- `statistic_values`: Statistic values per locus
- `outlier_indices`: Indices of outlier loci
- `positions`: Genomic positions (if None, uses indices)
- `fdr_threshold`: FDR threshold line

## Complete Function List

The module now includes **28 visualization functions** covering:

1. **Basic Comparisons**: Diversity, Tajima's D, Fst, summary statistics
2. **Population Structure**: PCA, kinship matrices, population trees
3. **Distribution Plots**: Allele frequencies, pairwise distances, heterozygosity
4. **Correlation Analysis**: π vs θ, correlation matrices, regression plots
5. **Multi-Population**: Fst matrices, admixture plots, isolation by distance
6. **Statistical Tests**: Neutrality test suites, HWE tests, bootstrap/permutation plots
7. **Specialized**: Site frequency spectra, linkage disequilibrium decay, demographic models

All functions follow a consistent API with optional `output_path` for saving figures and customizable titles and figure sizes.
