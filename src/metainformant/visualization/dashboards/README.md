# Visualization Dashboards

Composite multi-panel figures and interactive Plotly-based plots. Provides a generic panel builder for publication-ready figures, genomic overview dashboards, QC summary panels, and interactive scatter/heatmap/volcano/Manhattan plots with static matplotlib fallback.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `composite` and `interactive` modules |
| `composite.py` | Multi-panel figure builder, genomic overview, QC summary dashboards |
| `interactive.py` | Interactive Plotly plots with matplotlib fallback |

## Key Functions

| Function | Description |
|----------|-------------|
| `composite.multi_panel()` | Arrange multiple plots into a single gridded figure |
| `composite.genomic_overview()` | Create genomic overview dashboard with multiple panels |
| `composite.qc_summary()` | Generate QC summary dashboard for sample metrics |
| `interactive.interactive_scatter()` | Interactive scatter plot via Plotly |
| `interactive.interactive_heatmap()` | Interactive heatmap with hover tooltips |
| `interactive.interactive_volcano()` | Interactive volcano plot for DE results |
| `interactive.interactive_manhattan()` | Interactive Manhattan plot for GWAS results |

## Usage

```python
from metainformant.visualization.dashboards import composite, interactive

fig = composite.multi_panel(panels, ncols=2, title="Overview")
scatter = interactive.interactive_scatter(df, x="PC1", y="PC2", color="cluster")
```
