# VISUALIZATION

## Overview
Visualization and plotting utilities module for METAINFORMANT.

## 📦 Contents
- **[analysis/](analysis/)** — Statistical, quality, and dimensionality reduction plots
- **[genomics/](genomics/)** — Genomics-specific visualizations (expression, trees, networks)
- **[plots/](plots/)** — Core plotting functions (basic, specialized, animations)
- **[config/](config/)** — Color palettes and theme configuration
- **[dashboards/](dashboards/)** — Composite and interactive dashboard layouts
- **[interactive_dashboards/](interactive_dashboards/)** — Interactive web-based dashboards

## 📊 Structure

```mermaid
graph TD
    subgraph "Visualization Module"
        P[plots/] --> |basic.py| BA[Line, Scatter, Bar, Heatmap]
        P --> |general.py| GN[Volcano, Manhattan, PCA, QQ]
        P --> |specialized.py| SP[Domain-Specific Plots]
        P --> |multidim.py| MD[Multi-Dimensional Plots]
        P --> |animations.py| AN[Animated Plots]

        G[genomics/] --> |expression.py| EX[Expression Heatmaps]
        G --> |genomics.py| GM[Manhattan, Volcano, Ideogram]
        G --> |networks.py| NW[Network Graphs]
        G --> |trees.py| TR[Phylogenetic Trees]

        A[analysis/] --> |statistical.py| ST[Histogram, Boxplot, ROC]
        A --> |dimred.py| DR[PCA, UMAP, t-SNE Plots]
        A --> |quality.py| QA[Quality Assessment Plots]

        D[dashboards/] --> |composite.py| CM[Multi-Panel Layouts]
        D --> |interactive.py| IN[Interactive Plotly Plots]

        CF[config/] --> |themes.py| TH[Theme Management]
        CF --> |palettes.py| PL[Color Palettes]
    end

    CF --> P
    CF --> G
    P --> D
    G --> D
```

## Usage
Import module:
```python
from metainformant.visualization import ...
```
