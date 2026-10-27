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
| `interactive.interactive_scatter()` | Interactive scatter plot via Plotly; without Plotly, a numeric `color` column becomes a continuous colour scale and any other column one coloured series per category, with a legend |
| `interactive.interactive_heatmap()` | Interactive heatmap with hover tooltips |
| `interactive.interactive_volcano()` | Interactive volcano plot for DE results |
| `interactive.interactive_manhattan()` | Interactive Manhattan plot for GWAS results |

## Usage

```python
from metainformant.visualization.dashboards import composite, interactive

fig = composite.multi_panel(panels, ncols=2, title="Overview")
scatter = interactive.interactive_scatter(df, x="PC1", y="PC2", color="cluster")
```

## Output Paths and the Plotly Fallback

The four `interactive_*` functions return a Plotly `Figure` when Plotly is installed and a matplotlib `Axes` otherwise. What `output_path` does depends on that:

| Plotly | `output_path` | Result |
|--------|---------------|--------|
| installed | any path | Written as an HTML file (a suffix other than `.html` is replaced by `.html`); missing parent directories are created |
| missing | ends in `.html` or `.htm` (any letter case) | `ImportError` naming Plotly, raised before any figure, file or directory is created |
| missing | any other path | Saved as a static image through matplotlib (300 dpi, tight bounding box, format taken from the suffix); missing parent directories are created |
| either | omitted | Nothing is saved |

```python
import numpy as np
import pandas as pd

from metainformant.visualization.dashboards import interactive

rng = np.random.default_rng(0)
df = pd.DataFrame({"PC1": rng.normal(size=30), "PC2": rng.normal(size=30), "cluster": rng.choice(["a", "b"], 30)})

try:
    # With Plotly: writes output/dashboards/scatter.html and returns a Plotly Figure
    fig = interactive.interactive_scatter(df, x="PC1", y="PC2", color="cluster", output_path="output/dashboards/scatter.html")
except ImportError:
    # Without Plotly the HTML path is refused; an image path gives the static fallback (a matplotlib Axes)
    ax = interactive.interactive_scatter(df, x="PC1", y="PC2", color="cluster", output_path="output/dashboards/scatter.png")
```

Plotly is part of the `visualization` extra (`uv sync --extra visualization`).
