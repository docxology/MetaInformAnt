# DNA: Population Genetics Visualization

Visualization functions for population genetics analysis.

**Functions**: `plot_fst_matrix`, `plot_tajima_d_distribution`, `plot_selection_statistics`, `plot_population_diversity`, `plot_ld_decay`, `plot_population_structure`, `plot_demographic_history`, `create_population_summary_plot`, `plot_mutation_spectrum`, `plot_allele_frequency_spectrum`, `plot_bootstrap_distribution`, `plot_demographic_comparison`, `plot_diversity_comparison`, `plot_fst_comparison`, `plot_hardy_weinberg_test`, `plot_heterozygosity_distribution`, `plot_kinship_matrix`, `plot_linkage_disequilibrium_decay`, `plot_neutrality_test_suite`, `plot_neutrality_test_summary`, `plot_outlier_detection`, `plot_pca_results`, `plot_permutation_test`, `plot_pi_vs_theta`, `plot_site_frequency_spectrum`, `plot_statistic_correlation_matrix`, `plot_statistic_distribution`, `plot_summary_statistics_grid`, `plot_tajimas_d_comparison`.

## Overview

The `metainformant.dna.population.visualization` module provides matplotlib figures (seaborn heatmaps for three of them) for population genetics analysis, including diversity metrics, neutrality tests, population structure, and demographic summaries. It re-exports the 29 functions defined in `visualization_core` (11 core plots) and `visualization_stats` (18 statistical comparison plots), so one import path reaches all of them:

```python
from metainformant.dna.population.visualization import plot_diversity_comparison, plot_pca_results
```

### Common behavior

- **Saving**: pass `output_file` to save the figure at 300 DPI with a tight bounding box. Thirteen functions in `visualization_stats` also accept `output_path` as an alias; `plot_population_diversity` accepts it but ignores it, so `output_file` is the argument that works everywhere. The output directory is not created for you.
- **Returns**: the `matplotlib.figure.Figure`, already closed in pyplot so it is not displayed (save it, or adjust it through `fig.axes`), or `None` when matplotlib (or seaborn / SciPy where noted) is not installed. Most functions also return `None` for empty input.
- **Fixed titles and sizes**: none of the functions accept `title` or `figsize`; titles, labels and sizes are set inside each plot.
- **Optional packages**: `plot_fst_matrix`, `plot_kinship_matrix` and `plot_statistic_correlation_matrix` need seaborn; `plot_bootstrap_distribution` needs SciPy, and `plot_ld_decay` adds its exponential fit only when SciPy is installed.

## Functions

### `plot_diversity_comparison(diversity_data, output_file=None, output_path=None)`

Bar chart comparing nucleotide diversity (π) across scenarios or populations, with the value written above each bar.

**Parameters**:
- `diversity_data`: Dictionary mapping scenario names to π values
- `output_file`: Optional path to save the figure (alias `output_path`)

**Returns**: matplotlib Figure object

**Example**:
```python-snippet
from pathlib import Path

from metainformant.dna.population.visualization import plot_diversity_comparison

Path("output").mkdir(exist_ok=True)  # the plot functions do not create directories

diversity = {
    "Neutral": 0.01,
    "High Diversity": 0.05,
    "Low Diversity": 0.001,
}

fig = plot_diversity_comparison(
    diversity,
    output_file="output/diversity_comparison.png"
)
```

### `plot_tajimas_d_comparison(tajima_d_values, output_file=None, output_path=None)`

Bar chart comparing Tajima's D across scenarios with color coding.

**Parameters**:
- `tajima_d_values`: Dictionary mapping scenario names to Tajima's D values
- `output_file`: Optional path to save the figure (alias `output_path`)

**Returns**: matplotlib Figure object

**Features**:
- Red bars for D below -2, green bars for D above 2, blue bars in between
- Reference lines at D = 0 (solid), D = 2 (labelled "Balancing selection") and D = -2 (labelled "Positive selection")

**Example**:
```python-snippet
from pathlib import Path

from metainformant.dna.population.visualization import plot_tajimas_d_comparison

Path("output").mkdir(exist_ok=True)

tajimas_d = {
    "Neutral": 0.1,
    "Bottleneck": -1.5,
    "Expansion": -2.4,
    "Balanced": 2.3,
}

fig = plot_tajimas_d_comparison(
    tajimas_d,
    output_file="output/tajimas_d.png"
)
```

### `plot_fst_comparison(fst_data, output_file=None, output_path=None)`

Bar chart comparing Fst values across loci or population pairs.

**Parameters**:
- `fst_data`: Dictionary mapping locus or comparison names to Fst values
- `output_file`: Optional path to save the figure (alias `output_path`)

**Returns**: matplotlib Figure object

**Features**:
- A dashed reference line at Fst = 0.05, labelled "Significance threshold (0.05)" in the legend; it is a fixed line, not a computed test

**Example**:
```python-snippet
from pathlib import Path

from metainformant.dna.population.visualization import plot_fst_comparison

Path("output").mkdir(exist_ok=True)

fst_values = {
    "Pop1 vs Pop2": 0.2,
    "Pop1 vs Pop3": 0.05,
    "Pop2 vs Pop3": 0.3,
}

fig = plot_fst_comparison(fst_values, output_file="output/fst.png")
```

### `plot_pca_results(pca_result, sample_names=None, output_file=None, output_path=None, n_components=None)`

Plot the first two principal components.

**Parameters**:
- `pca_result`: Result dictionary from `compute_pca` (coordinates under `pcs` or `coordinates`, explained variance as fractions of the total under `explained_variance_ratio` or `explained_variance`; a `labels` or `cluster_labels` entry colors the points), or a `(coordinates, explained_variance, labels)` tuple. A bare coordinate array is not accepted: it returns `None` with a logged warning
- `sample_names`: Optional names drawn next to the points
- `output_file`: Optional path to save the figure (alias `output_path`)
- `n_components`: Accepted for compatibility but unused

**Returns**: matplotlib Figure object, or `None` when no coordinates are found or they are not two-dimensional with at least two columns. A result dictionary whose `status` is `"failed"` raises `ValidationError`.

**Plot**: a single PC1 versus PC2 scatter plot. When the explained variance is available, the axis labels show the percentage of variance of each component.

**Example**:
```python-snippet
import random
from pathlib import Path

from metainformant.dna.population.visualization import plot_pca_results
from metainformant.gwas.analysis.structure import compute_pca
from metainformant.simulation.models.popgen import generate_genotype_matrix

Path("output").mkdir(exist_ok=True)

# Individuals x sites matrix of 0/1/2 copies of the alternate allele
genotype_matrix = generate_genotype_matrix(n_individuals=60, n_sites=200, rng=random.Random(1))

pca_result = compute_pca(genotype_matrix, n_components=4)
fig = plot_pca_results(
    pca_result,
    output_file="output/pca.png"
)
```

### `plot_kinship_matrix(kinship_matrix, sample_names=None, output_file=None, output_path=None, max_samples=None)`

Plot a kinship matrix as a heatmap (requires seaborn).

**Parameters**:
- `kinship_matrix`: Result dictionary from `compute_kinship_matrix` (a `status` of `"failed"` raises `ValidationError`) or the kinship matrix itself as an array
- `sample_names`: Optional axis labels (default `Sample_0`, `Sample_1`, ...)
- `output_file`: Optional path to save the figure (alias `output_path`)
- `max_samples`: Plot only the first `max_samples` samples; applied only when a result dictionary is passed

**Returns**: matplotlib Figure object

**Example**:
```python-snippet
import random
from pathlib import Path

from metainformant.dna.population.visualization import plot_kinship_matrix
from metainformant.gwas.analysis.structure import compute_kinship_matrix
from metainformant.simulation.models.popgen import generate_genotype_matrix

Path("output").mkdir(exist_ok=True)

genotype_matrix = generate_genotype_matrix(n_individuals=60, n_sites=200, rng=random.Random(1))

kinship_result = compute_kinship_matrix(genotype_matrix, method="vanraden")
fig = plot_kinship_matrix(
    kinship_result,
    output_file="output/kinship.png",
    max_samples=30
)
```

### `plot_site_frequency_spectrum(sfs_data, output_file=None, output_path=None)`

Bar chart of a site frequency spectrum.

**Parameters**:
- `sfs_data`: List of site counts, one per frequency class starting with class 1, as returned by `generate_site_frequency_spectrum`
- `output_file`: Optional path to save the figure (alias `output_path`)

**Returns**: matplotlib Figure object

**Example**:
```python-snippet
import random
from pathlib import Path

from metainformant.dna.population.visualization import plot_site_frequency_spectrum
from metainformant.simulation.models.popgen import generate_site_frequency_spectrum

Path("output").mkdir(exist_ok=True)

sfs = generate_site_frequency_spectrum(
    sample_size=30,
    n_sites=100,
    rng=random.Random(1)
)
fig = plot_site_frequency_spectrum(
    sfs,
    output_file="output/sfs.png"
)
```

### `plot_neutrality_test_summary(summary, output_file=None, output_path=None)`

A 2x2 grid with one single-bar panel per neutrality test statistic.

**Parameters**:
- `summary`: Dictionary mapping test names to a number, or to a dictionary with `statistic` (or `value`) and an optional `p_value`, which is added to the panel title. Only the first four entries are drawn.
- `output_file`: Optional path to save the figure (alias `output_path`)

**Returns**: matplotlib Figure object

`neutrality_test_suite` returns a flat dictionary (`tajima_d`, `tajima_d_p_value`, `fu_li_d_star`, ...), not this shape, so pick out the statistics to plot.

**Example**:
```python-snippet
import random
from pathlib import Path

from metainformant.dna.population.analysis import neutrality_test_suite
from metainformant.dna.population.visualization import plot_neutrality_test_summary
from metainformant.simulation.models.popgen import generate_population_sequences

Path("output").mkdir(exist_ok=True)

sequences = generate_population_sequences(
    n_sequences=20, sequence_length=500, nucleotide_diversity=0.01, rng=random.Random(1)
)
suite = neutrality_test_suite(sequences)

fig = plot_neutrality_test_summary(
    {
        "Tajima's D": {"statistic": suite["tajima_d"], "p_value": suite["tajima_d_p_value"]},
        "Fu and Li's D*": {"statistic": suite["fu_li_d_star"], "p_value": suite["fu_li_d_star_p_value"]},
        "Fu and Li's F*": {"statistic": suite["fu_li_f_star"], "p_value": suite["fu_li_f_star_p_value"]},
        "Fay and Wu's H": {"statistic": suite["fay_wu_h"]},  # no p-value without an outgroup
    },
    output_file="output/neutrality_summary.png"
)
```

### `plot_demographic_comparison(pop1_data, pop2_data=None, output_file=None, output_path=None)`

Two bar-chart panels, "Effective Population Size" and "Observed Diversity", one bar per label.

**Parameters**:
- `pop1_data`: Dictionary mapping labels to dictionaries with `estimated_ne` and `observed_diversity`
- `pop2_data`: Optional second dictionary of the same shape. When given, the bars of both dictionaries are drawn side by side and labelled "Population 1" and "Population 2"; the tick labels come from `pop1_data`.
- `output_file`: Optional path to save the figure (alias `output_path`)

**Returns**: matplotlib Figure object

**Example**:
```python-snippet
from pathlib import Path

from metainformant.dna.population.visualization import plot_demographic_comparison

Path("output").mkdir(exist_ok=True)

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

fig = plot_demographic_comparison(
    demographic_results,
    output_file="output/demographic_comparison.png"
)
```

### `plot_summary_statistics_grid(stats_dict, output_file=None, output_path=None)`

A 3x2 grid with one bar-chart panel per summary statistic, one bar per scenario.

**Parameters**:
- `stats_dict`: Dictionary mapping scenario names to `{metric: value}` dictionaries; the values must be numbers
- `output_file`: Optional path to save the figure (alias `output_path`)

**Returns**: matplotlib Figure object

Metrics are sorted alphabetically and only the first six get a panel. A full `calculate_summary_statistics` result has nine metrics, so select the ones to plot, as below.

**Example**:
```python-snippet
import random
from pathlib import Path

from metainformant.dna.population.analysis import calculate_summary_statistics
from metainformant.dna.population.visualization import plot_summary_statistics_grid
from metainformant.simulation.models.popgen import generate_population_sequences

Path("output").mkdir(exist_ok=True)

METRICS = ("nucleotide_diversity", "segregating_sites", "watterson_theta", "tajima_d", "sample_size", "sequence_length")

def summarize(diversity, seed):
    sequences = generate_population_sequences(
        n_sequences=20, sequence_length=500, nucleotide_diversity=diversity, rng=random.Random(seed)
    )
    stats = calculate_summary_statistics(sequences=sequences)
    return {metric: stats[metric] for metric in METRICS}

summary_stats = {
    "Neutral": summarize(0.01, 1),
    "High Diversity": summarize(0.05, 2),
    "Low Diversity": summarize(0.002, 3),
}

fig = plot_summary_statistics_grid(
    summary_stats,
    output_file="output/summary_grid.png"
)
```

### `plot_linkage_disequilibrium_decay(ld_data, distances=None, output_file=None, output_path=None)`

Line plot of linkage disequilibrium (r²) against distance.

**Parameters**:
- `ld_data`: r² values, or `(distance, r²)` tuples such as those returned by `calculate_ld_decay`
- `distances`: Distances for plain r² values (default: 0, 1, 2, ...)
- `output_file`: Optional path to save the figure (alias `output_path`)

**Returns**: matplotlib Figure object, or `None` for empty input

**Example**:
```python-snippet
import random
from pathlib import Path

from metainformant.dna.population.analysis import calculate_ld_decay
from metainformant.dna.population.visualization import plot_linkage_disequilibrium_decay
from metainformant.simulation.models.popgen import generate_population_sequences

Path("output").mkdir(exist_ok=True)

r_squared_values = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05]
distances = [0, 1, 2, 3, 4, 5]  # In units (e.g., kb)

fig = plot_linkage_disequilibrium_decay(
    r_squared_values,
    distances=distances,
    output_file="output/ld_decay.png"
)

# Or pass (distance, r²) tuples straight from calculate_ld_decay
sequences = generate_population_sequences(
    n_sequences=20, sequence_length=500, nucleotide_diversity=0.01, rng=random.Random(1)
)
fig = plot_linkage_disequilibrium_decay(
    calculate_ld_decay(sequences, max_distance=50),
    output_file="output/ld_decay_sequences.png"
)
```

## Integration Examples

### Complete Workflow with Visualizations

```python-snippet
import random
from pathlib import Path

from metainformant.dna.population.analysis import calculate_summary_statistics, neutrality_test_suite
from metainformant.dna.population.visualization import (
    plot_diversity_comparison,
    plot_neutrality_test_summary,
)
from metainformant.simulation.models.popgen import generate_population_sequences

Path("output").mkdir(exist_ok=True)

neutral_seqs = generate_population_sequences(
    n_sequences=20, sequence_length=500, nucleotide_diversity=0.01, rng=random.Random(1)
)

# Analyze data
neutral_stats = calculate_summary_statistics(sequences=neutral_seqs)
neutrality = neutrality_test_suite(neutral_seqs)

# Create visualizations
plot_diversity_comparison(
    {"Neutral": neutral_stats["nucleotide_diversity"]},
    output_file="output/diversity.png"
)

plot_neutrality_test_summary(
    {
        "Tajima's D": {"statistic": neutrality["tajima_d"], "p_value": neutrality["tajima_d_p_value"]},
        "Fu and Li's D*": {"statistic": neutrality["fu_li_d_star"], "p_value": neutrality["fu_li_d_star_p_value"]},
    },
    output_file="output/neutrality.png"
)
```

### Batch Visualization Generation

```python-snippet
import random
from pathlib import Path

from metainformant.dna.population import hudson_fst, nucleotide_diversity, tajimas_d
from metainformant.dna.population.visualization import (
    plot_diversity_comparison,
    plot_tajimas_d_comparison,
    plot_fst_comparison,
)
from metainformant.simulation.models.popgen import (
    generate_population_sequences,
    generate_two_populations,
    simulate_bottleneck_population,
    simulate_population_expansion,
)

Path("output").mkdir(exist_ok=True)
rng = random.Random(1)

# Simulate the scenarios to compare
scenarios = {
    "neutral": generate_population_sequences(n_sequences=20, sequence_length=500, nucleotide_diversity=0.01, rng=rng),
    "high_diversity": generate_population_sequences(n_sequences=20, sequence_length=500, nucleotide_diversity=0.05, rng=rng),
    "low_diversity": generate_population_sequences(n_sequences=20, sequence_length=500, nucleotide_diversity=0.002, rng=rng),
    "bottleneck": simulate_bottleneck_population(n_sequences=20, sequence_length=500, bottleneck_size=5, rng=rng),
    "expansion": simulate_population_expansion(n_sequences=20, sequence_length=500, expansion_factor=10.0, rng=rng),
}
diversity_values = {name: nucleotide_diversity(seqs) for name, seqs in scenarios.items()}
tajimas_d_values = {name: tajimas_d(seqs) for name, seqs in scenarios.items()}

fst_values = {}
for target in (0.05, 0.2, 0.5):
    pop1, pop2 = generate_two_populations(n_pop1=10, n_pop2=10, sequence_length=500, fst=target, n_sites=250, rng=rng)
    fst_values[f"target {target}"] = hudson_fst(pop1, pop2)

# Generate all comparison plots
plot_diversity_comparison(
    diversity_values,
    output_file="output/diversity_comparison.png"
)

plot_tajimas_d_comparison(
    tajimas_d_values,
    output_file="output/tajimas_d_comparison.png"
)

plot_fst_comparison(
    fst_values,
    output_file="output/fst_comparison.png"
)
```

## Best Practices

### 1. High-Resolution Output

All functions support saving at 300 DPI for publication:

```python
from pathlib import Path

from metainformant.dna.population.visualization import plot_diversity_comparison

Path("output").mkdir(exist_ok=True)
diversity_values = {"Neutral": 0.01, "High Diversity": 0.05}

plot_diversity_comparison(
    diversity_values,
    output_file="output/diversity.png"  # Saved at 300 DPI
)
```

### 2. Custom Styling

Functions return matplotlib Figure objects for further customization:

```python
from pathlib import Path

from metainformant.dna.population.visualization import plot_diversity_comparison

Path("output").mkdir(exist_ok=True)
diversity_values = {"Neutral": 0.01, "High Diversity": 0.05}

fig = plot_diversity_comparison(diversity_values)
ax = fig.axes[0]
ax.set_ylabel("Diversity (π)", fontsize=14)
fig.savefig("output/custom.png", dpi=300)
```

### 3. Combining Multiple Plots

The plot functions draw on figures of their own and take no `ax` argument. To combine them, save each plot and place the images side by side with matplotlib:

```python
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

from metainformant.dna.population.visualization import plot_diversity_comparison, plot_tajimas_d_comparison

Path("output").mkdir(exist_ok=True)

# Plot 1: Diversity
plot_diversity_comparison({"Neutral": 0.01, "High Diversity": 0.05}, output_file="output/diversity.png")

# Plot 2: Tajima's D
plot_tajimas_d_comparison({"Neutral": 0.1, "Bottleneck": -1.5}, output_file="output/tajimas_d.png")

fig, axes = plt.subplots(1, 2, figsize=(15, 6))
for ax, path in zip(axes, ["output/diversity.png", "output/tajimas_d.png"]):
    ax.imshow(mpimg.imread(path))
    ax.axis("off")

plt.tight_layout()
fig.savefig("output/combined.png", dpi=300)
```

## Additional Visualization Functions

### Comparison Plots

#### `plot_population_diversity(diversity_data, output_file=None, output_path=None)`

Bar chart of π per population with the value (four decimals) above each bar. It draws the same kind of chart as `plot_diversity_comparison`. `output_path` is accepted but ignored; use `output_file`.

**Parameters**:
- `diversity_data`: Dictionary mapping population names to π values

#### `create_population_summary_plot(population_data, output_file=None)`

One bar-chart panel per metric, stacked vertically, one bar per population.

**Parameters**:
- `population_data`: Dictionary mapping population names to `{metric: value}` dictionaries; a metric missing for a population is plotted as 0

```python-snippet
from pathlib import Path

from metainformant.dna.population.visualization import create_population_summary_plot, plot_population_diversity

Path("output").mkdir(exist_ok=True)

fig = plot_population_diversity(
    {"pop1": 0.01, "pop2": 0.015, "pop3": 0.008},
    output_file="output/population_diversity.png"
)
fig = create_population_summary_plot(
    {
        "pop1": {"diversity": 0.01, "fst": 0.05},
        "pop2": {"diversity": 0.015, "fst": 0.03},
    },
    output_file="output/population_summary.png"
)
```

### Distribution Plots

#### `plot_allele_frequency_spectrum(allele_frequencies, output_file=None)`

Histogram of allele frequencies in ten equal bins from 0 to 1.

**Parameters**:
- `allele_frequencies`: List of allele frequencies between 0 and 1 (for example from `dna.population.allele_frequencies`)

#### `plot_heterozygosity_distribution(het_data, output_file=None)`

Histogram (20 bins) of heterozygosity values.

**Parameters**:
- `het_data`: List of per-site heterozygosity values, or a dictionary mapping population names to such lists (overlaid histograms with a legend)

#### `plot_statistic_distribution(stat_values, stat_name="Statistic", output_file=None, plot_type="histogram")`

Histogram (20 bins) of a statistic. For a single list, dashed lines mark the mean and the median.

**Parameters**:
- `stat_values`: List of values, or a dictionary mapping names to lists of values (overlaid histograms)
- `stat_name`: Axis label and title text
- `plot_type`: Accepted but unused; the plot is always a histogram

#### `plot_tajima_d_distribution(tajima_d_values, output_file=None)`

Histogram (30 bins) of Tajima's D values with a dashed line at D = 0 labelled "Neutral expectation".

**Parameters**:
- `tajima_d_values`: List of Tajima's D values

#### `plot_mutation_spectrum(mutation_counts, output_file=None)`

Bar chart of mutation counts by type, with the count above each bar.

**Parameters**:
- `mutation_counts`: Dictionary mapping mutation types (such as `"A>G"`) to counts

```python-snippet
import random

from metainformant.dna.population import allele_frequencies, nucleotide_diversity, tajimas_d
from metainformant.dna.population.visualization import (
    plot_allele_frequency_spectrum,
    plot_heterozygosity_distribution,
    plot_mutation_spectrum,
    plot_statistic_distribution,
    plot_tajima_d_distribution,
)
from metainformant.dna.variation.mutations import calculate_substitution_matrix
from metainformant.simulation.models.popgen import generate_genotype_matrix, simulate_bottleneck_population
from metainformant.simulation.models.sequences import generate_random_dna, mutate_sequence

genotypes = generate_genotype_matrix(n_individuals=100, n_sites=300, rng=random.Random(1))
by_site = list(zip(*genotypes))  # one row per site, as allele_frequencies expects

fig = plot_allele_frequency_spectrum(allele_frequencies(by_site))

# Observed heterozygosity per site (genotype 1 is heterozygous)
het = [sum(1 for genotype in site if genotype == 1) / len(site) for site in by_site]
fig = plot_heterozygosity_distribution(het)
fig = plot_heterozygosity_distribution({"Pop A": het[:150], "Pop B": het[150:]})

# Distribution of a statistic over simulation replicates
replicates = [
    simulate_bottleneck_population(n_sequences=20, sequence_length=500, bottleneck_size=5, rng=random.Random(seed))
    for seed in range(40)
]
fig = plot_statistic_distribution([nucleotide_diversity(seqs) for seqs in replicates], stat_name="Nucleotide diversity (π)")
fig = plot_tajima_d_distribution([tajimas_d(seqs) for seqs in replicates])

# Mutation spectrum between an ancestor and a mutated copy
ancestor = generate_random_dna(1000, rng=random.Random(2))
derived = mutate_sequence(ancestor, n_mut=100, rng=random.Random(3))
matrix = calculate_substitution_matrix(ancestor, derived)
fig = plot_mutation_spectrum({f"{a}>{b}": count for (a, b), count in sorted(matrix.items())})
```

### Correlation and Relationship Plots

#### `plot_pi_vs_theta(pi_values, theta_values, locus_names=None, output_file=None)`

Scatter plot of π (y-axis) against Watterson's θ (x-axis) with a dashed 1:1 line. `locus_names` is accepted but unused. Put both on the same per-site scale: `wattersons_theta` returns S / a_n for the whole alignment, so divide it by the sequence length.

**Parameters**:
- `pi_values`: Nucleotide diversity (π) values
- `theta_values`: Watterson's theta (θ_W) values, same length as `pi_values`

#### `plot_statistic_correlation_matrix(stats_data, output_file=None)`

Annotated heatmap (seaborn) of the correlations between statistics.

**Parameters**:
- `stats_data`: Dictionary mapping statistic names to equal-length lists of values

```python-snippet
import random

from metainformant.dna.population import nucleotide_diversity, segregating_sites, tajimas_d, wattersons_theta
from metainformant.dna.population.visualization import plot_pi_vs_theta, plot_statistic_correlation_matrix
from metainformant.simulation.models.popgen import generate_population_sequences

sequences = generate_population_sequences(
    n_sequences=20, sequence_length=1000, nucleotide_diversity=0.01, rng=random.Random(1)
)

# Per-window statistics over 100 bp windows
windows = [[seq[start:start + 100] for seq in sequences] for start in range(0, 1000, 100)]
pi = [nucleotide_diversity(window) for window in windows]
theta = [wattersons_theta(window) / 100 for window in windows]  # per site, like pi

# theta_W sits well above pi here: the simulated genealogy is star-like (negative Tajima's D)
fig = plot_pi_vs_theta(pi, theta)

fig = plot_statistic_correlation_matrix(
    {
        "pi": pi,
        "theta_w": theta,
        "segregating_sites": [segregating_sites(window) for window in windows],
        "tajimas_d": [tajimas_d(window) for window in windows],
    }
)
```

### Population Structure and Multi-Population Plots

#### `plot_population_structure(pca_coords, cluster_labels=None, output_file=None)`

Scatter plot of PCA coordinates: two-dimensional for two components, three-dimensional (first three components) for three or more. Integer `cluster_labels` color the points and add a "Population Cluster" color bar.

**Parameters**:
- `pca_coords`: NumPy array of shape (n_samples, n_components); fewer than two columns returns `None`
- `cluster_labels`: Optional list of integer cluster assignments

#### `plot_fst_matrix(data, pop_names=None, output_file=None)`

Annotated heatmap (seaborn) of pairwise Fst values.

**Parameters**:
- `data`: Dictionary mapping population names to lists of aligned sequences (pairwise values are computed with `calculate_fst`, the heterozygosity-based estimator rather than `hudson_fst`), or a square Fst matrix as an array or nested lists
- `pop_names`: Names for the rows and columns of a matrix (default `Pop1`, `Pop2`, ...)

```python-snippet
import random

import numpy as np

from metainformant.dna.population.visualization import plot_fst_matrix, plot_population_structure
from metainformant.gwas.analysis.structure import compute_pca
from metainformant.simulation.models.popgen import generate_genotype_matrix, generate_two_populations

# Two simulated groups with different allele frequencies, stacked into one matrix
group_a = generate_genotype_matrix(n_individuals=30, n_sites=200, rng=random.Random(1))
group_b = generate_genotype_matrix(n_individuals=30, n_sites=200, rng=random.Random(2))
pca_result = compute_pca(group_a + group_b, n_components=3)
fig = plot_population_structure(np.array(pca_result["pcs"]), cluster_labels=[0] * 30 + [1] * 30)

# Pairwise Fst from sequences, or from a precomputed matrix
pop1, pop2 = generate_two_populations(
    n_pop1=10, n_pop2=10, sequence_length=500, fst=0.3, n_sites=250, rng=random.Random(3)
)
fig = plot_fst_matrix({"Pop1": pop1, "Pop2": pop2, "Mixed": pop1[:5] + pop2[:5]})
fig = plot_fst_matrix([[0, 0.1, 0.3], [0.1, 0, 0.2], [0.3, 0.2, 0]], pop_names=["A", "B", "C"])
```

### Statistical Test Visualizations

#### `plot_neutrality_test_suite(results, output_file=None, output_path=None)`

Bar chart with one bar per test statistic, red when negative and blue otherwise, and a line at zero.

**Parameters**:
- `results`: Dictionary mapping test names to a number, or to a dictionary with `statistic` (or `value`)

#### `plot_hardy_weinberg_test(results, output_file=None)`

Bars of -log10(p-value) per locus with a dashed line at p = 0.05.

**Parameters**:
- `results`: List of dictionaries, one per locus, with a `p_value` (default 1.0) and an optional `locus` name (default `Locus_<index>`)

#### `plot_bootstrap_distribution(bootstrap_values, observed_value=None, output_file=None)`

Density histogram (30 bins) of bootstrap replicates with a kernel density estimate, and vertical lines for the observed value (when given), the 2.5th and 97.5th percentiles (a 95% interval), the mean and the median. Requires SciPy.

**Parameters**:
- `bootstrap_values`: List of bootstrap replicate values
- `observed_value`: Observed statistic value to mark on the plot

#### `plot_permutation_test(results, observed_value=None, p_value=None, output_file=None, output_path=None)`

Histogram of permuted statistics with a dashed line at the observed value. The p-value is only shown in the legend when you pass it; it is not computed.

**Parameters**:
- `results`: List of permuted values, or a dictionary with `permuted_values` (or `null_distribution`), `observed_value` and `p_value`
- `observed_value`: Observed statistic value (overrides the dictionary entry)
- `p_value`: p-value to display (overrides the dictionary entry)

#### `plot_outlier_detection(results, outliers=None, output_file=None, output_path=None)`

Scatter plot of statistic values by observation index, with the indices listed in `outliers` highlighted in red.

**Parameters**:
- `results`: List of numbers, or of dictionaries with `value` (or `statistic`)
- `outliers`: Indices of the observations to highlight

#### `plot_selection_statistics(statistics, output_file=None)`

One stacked panel per statistic, drawn as a line with markers along the sequence, with a dashed zero line.

**Parameters**:
- `statistics`: Dictionary mapping statistic names to lists of values (for example windowed Tajima's D)

```python-snippet
import random

import numpy as np

from metainformant.dna.population import nucleotide_diversity, tajimas_d
from metainformant.dna.population.analysis import neutrality_test_suite
from metainformant.dna.population.visualization import (
    plot_bootstrap_distribution,
    plot_hardy_weinberg_test,
    plot_neutrality_test_suite,
    plot_outlier_detection,
    plot_permutation_test,
    plot_selection_statistics,
)
from metainformant.simulation.models.popgen import generate_population_sequences

rng = random.Random(1)
sequences = generate_population_sequences(n_sequences=20, sequence_length=500, nucleotide_diversity=0.01, rng=rng)

suite = neutrality_test_suite(sequences)
fig = plot_neutrality_test_suite(
    {
        "Tajima's D": suite["tajima_d"],
        "Fu and Li's D*": suite["fu_li_d_star"],
        "Fu and Li's F*": suite["fu_li_f_star"],
        "Fay and Wu's H": suite["fay_wu_h"],
    }
)

# One dictionary per locus
fig = plot_hardy_weinberg_test(
    [
        {"locus": "L1", "p_value": 0.40},
        {"locus": "L2", "p_value": 0.01},
        {"locus": "L3", "p_value": 0.0004},
    ]
)

# Bootstrap of nucleotide diversity: resample the sequences with replacement
bootstrap = [nucleotide_diversity(rng.choices(sequences, k=len(sequences))) for _ in range(200)]
fig = plot_bootstrap_distribution(bootstrap, observed_value=nucleotide_diversity(sequences))

# Permutation test on a difference of means; the p-value is yours to compute
np_rng = np.random.default_rng(1)
group_a = np_rng.normal(0.30, 0.05, 40)
group_b = np_rng.normal(0.34, 0.05, 40)
observed = group_b.mean() - group_a.mean()
pooled = np.concatenate([group_a, group_b])
permuted = []
for _ in range(1000):
    np_rng.shuffle(pooled)
    permuted.append(pooled[40:].mean() - pooled[:40].mean())
p_value = (sum(abs(value) >= abs(observed) for value in permuted) + 1) / (len(permuted) + 1)
fig = plot_permutation_test(permuted, observed_value=observed, p_value=p_value)

# Outliers among per-locus values
values = np_rng.normal(0, 1, 200)
values[[20, 120]] = [4.5, -4.2]
fig = plot_outlier_detection(values.tolist(), outliers=[i for i, v in enumerate(values) if abs(v) > 4])

# Statistics along the sequence in 100 bp windows
windows = [[seq[start:start + 100] for seq in sequences] for start in range(0, 500, 100)]
fig = plot_selection_statistics(
    {
        "tajima_d": [tajimas_d(window) for window in windows],
        "nucleotide_diversity": [nucleotide_diversity(window) for window in windows],
    }
)
```

### Specialized Plots

#### `plot_ld_decay(ld_data, output_file=None)`

Line plot of LD (r²) against distance in bp. With more than three points and SciPy installed, a dashed exponential-decay fit `a * exp(-b * distance)` is added and `b` is shown in the legend.

**Parameters**:
- `ld_data`: List of `(distance, r²)` tuples, as returned by `calculate_ld_decay`; an empty list returns `None`

#### `plot_demographic_history(ne_estimates, generations, output_file=None)`

Effective population size against generation on a logarithmic y-axis, with each estimate written next to its point.

**Parameters**:
- `ne_estimates`: Effective population size estimates
- `generations`: Generation numbers, one per estimate

```python-snippet
import random

from metainformant.dna.population.analysis import calculate_ld_decay
from metainformant.dna.population.visualization import plot_demographic_history, plot_ld_decay
from metainformant.simulation.models.popgen import generate_population_sequences

sequences = generate_population_sequences(
    n_sequences=20, sequence_length=500, nucleotide_diversity=0.01, rng=random.Random(1)
)
fig = plot_ld_decay(calculate_ld_decay(sequences, max_distance=100))

fig = plot_demographic_history([1000, 1200, 800, 1500], [0, 1000, 2000, 3000])
```

## Complete Function List

The module provides **29 visualization functions** (11 in `visualization_core`, 18 in `visualization_stats`), all importable from `metainformant.dna.population.visualization`:

1. **Basic Comparisons**: Diversity (`plot_diversity_comparison`, `plot_population_diversity`), Tajima's D, Fst, summary statistics (`plot_summary_statistics_grid`, `create_population_summary_plot`)
2. **Population Structure**: PCA (`plot_pca_results`, `plot_population_structure`), kinship matrices
3. **Distribution Plots**: Allele frequencies, heterozygosity, statistic and Tajima's D distributions, mutation spectrum
4. **Correlation Analysis**: π vs θ, correlation matrices
5. **Multi-Population**: Fst matrices
6. **Statistical Tests**: Neutrality test suites and summaries, HWE tests, bootstrap and permutation plots, outlier detection, selection statistics along a sequence
7. **Specialized**: Site frequency spectra, linkage disequilibrium decay (`plot_linkage_disequilibrium_decay`, `plot_ld_decay`), demographic comparison and history

All functions follow a consistent API: the first argument is the data, `output_file` saves the figure at 300 DPI, and the return value is the matplotlib figure. Titles and figure sizes are fixed.

## See Also

- [`dna.population`](population.md) - Core population genetics functions
- [Population workflows](population_workflows.md) - Analysis workflows
- [`visualization.plots`](../visualization/plots.md) - General plotting utilities
