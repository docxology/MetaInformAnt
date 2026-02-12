# Phenotype Visualization

Comprehensive visualization functions for phenotypic data including trait
distributions, life course trajectories, morphological measurements,
behavioural patterns, PCA, heritability, networks, and interactive browsers.

## Key Concepts

All plot functions follow a consistent API pattern: they accept data, optional
matplotlib Axes for composition, an optional output path for saving, and a
configurable figure size. Every function returns a matplotlib Axes object for
further customisation.

Optional dependencies: **seaborn** enhances histograms and heatmaps;
**plotly** enables the interactive phenotype browser; **networkx** is required
for phenotype network plots.

## Function Reference

### plot_trait_distribution

```python
def plot_trait_distribution(
    trait_values: np.ndarray,
    trait_name: str = "Trait",
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 6),
    **kwargs,
) -> Axes
```

Histogram with kernel density estimate overlay. Uses seaborn `histplot` when
available, falls back to matplotlib histogram with `scipy.stats.gaussian_kde`.

### plot_trait_correlation_matrix

```python
def plot_trait_correlation_matrix(
    trait_data: pd.DataFrame,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Axes
```

Lower-triangle heatmap of pairwise Pearson correlations. Annotated with
correlation values when seaborn is available.

### plot_life_course_trajectory

```python
def plot_life_course_trajectory(
    life_events: List[Dict[str, Any]],
    individual_id: str | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes
```

Timeline plot of life events. Each event dictionary must contain `age` and
`event_type` keys. Events are colour-coded by type with annotated labels.

### plot_morphological_measurements

```python
def plot_morphological_measurements(
    measurements: Dict[str, np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes
```

Box plots comparing distributions of multiple named measurements. Uses
seaborn `boxplot` or colour-coded matplotlib patches.

### plot_behavioral_patterns

```python
def plot_behavioral_patterns(
    behavioral_data: pd.DataFrame,
    time_column: str = "time",
    behavior_column: str = "behavior",
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes
```

Overlaid histograms showing temporal frequency distributions of distinct
behaviour categories.

### plot_phenotype_pca

```python
def plot_phenotype_pca(
    phenotype_data: np.ndarray,
    trait_names: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 6),
    **kwargs,
) -> Axes
```

PCA scatter plot (PC1 vs PC2) with explained variance annotations. Requires
scikit-learn for PCA computation.

### plot_trait_heritability

```python
def plot_trait_heritability(
    heritability_estimates: Dict[str, float],
    confidence_intervals: Dict[str, Tuple[float, float]] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes
```

Bar chart of heritability (h-squared) estimates per trait with optional error
bars for confidence intervals. Y-axis is clamped to [0, 1].

### plot_life_history_comparison

```python
def plot_life_history_comparison(
    species_data: Dict[str, List[Dict[str, Any]]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes
```

Cumulative event curves across species for comparative life history analysis.

### plot_phenotype_network

```python
def plot_phenotype_network(
    phenotype_correlations: np.ndarray,
    phenotype_names: List[str],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Axes
```

Spring-layout network graph where edges connect phenotypes with correlations
exceeding a threshold (default 0.5, configurable via `threshold` kwarg).
Requires networkx.

### create_interactive_phenotype_browser

```python
def create_interactive_phenotype_browser(
    phenotype_data: pd.DataFrame,
    *,
    output_path: str | Path | None = None,
    **kwargs,
) -> Any
```

Parallel coordinates plot via Plotly for interactive exploration. Saves to
HTML when output_path is given. Requires plotly.

## Usage Example

```python
import numpy as np
import pandas as pd
from metainformant.phenotype.visualization.visualization import (
    plot_trait_distribution,
    plot_trait_correlation_matrix,
    plot_phenotype_pca,
)

# Single trait distribution
values = np.random.normal(10.0, 2.0, 200)
ax = plot_trait_distribution(values, trait_name="Body Length (mm)")

# Correlation matrix
df = pd.DataFrame({
    "length": np.random.normal(10, 2, 100),
    "width": np.random.normal(5, 1, 100),
    "mass": np.random.normal(15, 3, 100),
})
ax = plot_trait_correlation_matrix(df, output_path="output/phenotype/corr.png")

# PCA
data = np.random.randn(50, 5)
ax = plot_phenotype_pca(data, output_path="output/phenotype/pca.png")
```

## Related Modules

- `metainformant.phenotype.analysis.life_course` -- life course data sources
- `metainformant.phenotype.morphological` -- morphometric data classes
- `metainformant.visualization` -- generic 80+ plot types
