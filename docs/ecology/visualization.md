# Ecology Visualization

Visualization functions for ecological data including species abundance distributions, community composition, ordination plots, diversity comparisons, distance heatmaps, and ecological networks.

## Key Concepts

All plot functions return a matplotlib `Axes` object that can be further customized. Figures can be saved to disk via the `output_path` parameter. Interactive dashboards require Plotly.

## Function Reference

### `plot_species_abundance_distribution(abundance_data, species_names=None, ax=None, output_path=None, figsize=(10,6), **kwargs) -> Axes`

Rank-abundance curve (Whittaker plot). Pass `log_scale=True` in kwargs for log y-axis. Top 5 species are labeled when names are provided.

### `plot_diversity_accumulation_curve(diversity_data, ax=None, output_path=None, figsize=(10,6)) -> Axes`

Species accumulation / rarefaction curve with optional 95% confidence intervals. Input is a list of dicts with `sample_size`, `species_count`, and optional `confidence_lower`/`confidence_upper`.

### `plot_community_composition(community_matrix, species_names=None, sample_names=None, ax=None, output_path=None, figsize=(12,8)) -> Axes`

Stacked bar chart of relative abundances across samples.

### `plot_beta_diversity_ordination(ordination_coords, sample_groups=None, ax=None, output_path=None, figsize=(8,6)) -> Axes`

Scatter plot of PCoA/NMDS coordinates. Points colored by `sample_groups` if provided, otherwise by sample index with a colorbar.

### `plot_diversity_indices_comparison(diversity_indices, index_names=None, ax=None, output_path=None, figsize=(10,6)) -> Axes`

Grouped bar chart comparing multiple diversity indices across samples. Input is a dict mapping index names to arrays.

### `plot_rank_abundance_curve_comparison(abundance_datasets, ax=None, output_path=None, figsize=(10,6)) -> Axes`

Overlay multiple rank-abundance curves for cross-community comparison.

### `plot_biodiversity_rarefaction(rarefaction_data, ax=None, output_path=None, figsize=(10,6)) -> Axes`

Rarefaction curves for multiple communities. Input is a dict mapping community names to lists of `{sample_size, species_count}` dicts.

### `plot_ecological_distance_heatmap(distance_matrix, sample_names=None, ax=None, output_path=None, figsize=(10,8)) -> Axes`

Heatmap of pairwise ecological distances. Labels shown for 20 or fewer samples.

### `plot_ecological_network(interaction_matrix, species_names=None, ax=None, output_path=None, figsize=(10,8), **kwargs) -> Axes`

Spring-layout network of species interactions. Requires `networkx`. Pass `threshold` in kwargs to control minimum edge weight.

### `create_interactive_ecology_dashboard(ecology_data, output_path=None, **kwargs) -> plotly.Figure`

Interactive Plotly dashboard. Requires `plotly`. Pass `diversity_indices` key in the data dict. Saves as HTML when `output_path` is given.

## Usage Examples

```python
import numpy as np
from metainformant.ecology.visualization.visualization import (
    plot_species_abundance_distribution,
    plot_beta_diversity_ordination,
    plot_ecological_distance_heatmap,
)

# Rank-abundance curve
abundances = np.array([45, 23, 12, 8, 5, 3, 2, 1, 1])
ax = plot_species_abundance_distribution(abundances, log_scale=True)

# Ordination plot with groups
coords = np.array([[0.1, 0.2], [-0.3, 0.1], [0.2, -0.4], [-0.1, 0.3]])
groups = np.array(["forest", "forest", "grassland", "grassland"])
ax = plot_beta_diversity_ordination(coords, groups, output_path="output/ordination.png")

# Distance heatmap
dm = np.array([[0, 0.3, 0.7], [0.3, 0, 0.5], [0.7, 0.5, 0]])
ax = plot_ecological_distance_heatmap(dm, sample_names=["A", "B", "C"])
```

## Optional Dependencies

- `seaborn` -- enhanced styling (optional)
- `plotly` -- interactive dashboard (required for `create_interactive_ecology_dashboard`)
- `networkx` -- ecological network plots (required for `plot_ecological_network`)

## Configuration

Environment variable prefix: `ECO_`

## Related Modules

- `metainformant.ecology.community` -- data sources for diversity plots
- `metainformant.ecology.ordination` -- coordinates for ordination plots
- `metainformant.visualization` -- general-purpose plot library
