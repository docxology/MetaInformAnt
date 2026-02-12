# Phenotype Visualization

Plotting functions for phenotypic trait distributions, life course trajectories, morphological measurements, and behavioral patterns.

## Contents

| File | Purpose |
|------|---------|
| `visualization.py` | Trait distribution, correlation, trajectory, heritability, and network plots |

## Key Functions

| Function | Description |
|----------|-------------|
| `plot_trait_distribution()` | Histogram or density plot for trait values |
| `plot_trait_correlation_matrix()` | Correlation heatmap across multiple traits |
| `plot_life_course_trajectory()` | Longitudinal trait trajectories over age |
| `plot_morphological_measurements()` | Body measurement comparisons across groups |
| `plot_behavioral_patterns()` | Temporal behavioral pattern visualization |
| `plot_phenotype_pca()` | PCA of multi-trait phenotype space |
| `plot_trait_heritability()` | Heritability estimates with confidence intervals |
| `plot_life_history_comparison()` | Life history trait comparison across species |
| `plot_phenotype_network()` | Network of trait correlations |
| `create_interactive_phenotype_browser()` | Interactive Plotly phenotype explorer |

## Usage

```python
from metainformant.phenotype.visualization.visualization import (
    plot_trait_distribution,
    plot_life_course_trajectory,
)

plot_trait_distribution(values, trait_name="body_size", output_path="output/trait.png")
plot_life_course_trajectory(sequences, output_path="output/trajectory.png")
```
