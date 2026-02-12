# Visualization Interactive Dashboards

JSON-serializable interactive dashboard data structures for web-based visualization. Generates plot specifications for scatter plots, heatmaps, genome browser tracks, volcano plots, and composed dashboards, with standalone HTML export using embedded SVG and vanilla JavaScript.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `dashboards` module |
| `dashboards.py` | JSON-serializable plot builders, HTML export, dashboard composition |

## Key Functions

| Function | Description |
|----------|-------------|
| `create_interactive_scatter()` | Generate scatter plot data structure with hover info |
| `create_interactive_heatmap()` | Generate heatmap data with hierarchical clustering order |
| `create_genome_browser_track()` | Generate genome browser track data for genomic regions |
| `create_interactive_volcano()` | Generate volcano plot data for differential expression |
| `export_to_html()` | Export dashboard to standalone HTML file |
| `create_dashboard()` | Compose multiple plot specifications into a dashboard |

## Usage

```python
from metainformant.visualization.interactive_dashboards import dashboards

scatter = dashboards.create_interactive_scatter(x_vals, y_vals, labels=names)
dashboard = dashboards.create_dashboard([scatter, heatmap], title="Analysis")
dashboards.export_to_html(dashboard, "output/dashboard.html")
```
