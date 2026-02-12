# Ecology Visualization

Plotting functions for ecological community analysis including species abundance, diversity, ordination, and network visualizations.

## Contents

| File | Purpose |
|------|---------|
| `visualization.py` | Ecological plots: abundance distributions, ordination, diversity curves |

## Key Functions

| Function | Description |
|----------|-------------|
| `plot_species_abundance_distribution()` | Species abundance histogram or rank-abundance curve |
| `plot_diversity_accumulation_curve()` | Rarefaction or species accumulation visualization |
| `plot_community_composition()` | Stacked bar chart of community composition |
| `plot_beta_diversity_ordination()` | PCoA/NMDS ordination scatter plot |
| `plot_diversity_indices_comparison()` | Multi-panel comparison of diversity metrics |
| `plot_ecological_network()` | Bipartite or interaction network diagram |
| `plot_rank_abundance_curve_comparison()` | Overlaid rank-abundance curves for multiple sites |
| `plot_biodiversity_rarefaction()` | Rarefaction curves with confidence intervals |
| `plot_ecological_distance_heatmap()` | Distance or similarity matrix heatmap |
| `create_interactive_ecology_dashboard()` | Interactive Plotly dashboard for exploration |

## Usage

```python
from metainformant.ecology.visualization.visualization import (
    plot_species_abundance_distribution,
    plot_beta_diversity_ordination,
)

plot_species_abundance_distribution(abundances, output_path="output/sad.png")
plot_beta_diversity_ordination(coords, groups=labels, output_path="output/ord.png")
```
