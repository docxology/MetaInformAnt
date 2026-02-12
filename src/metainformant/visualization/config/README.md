# Visualization Config

Predefined color palettes and themes for consistent, publication-ready biological data visualization. Includes colorblind-safe chromosome palettes, expression gradients, significance coloring, and scoped matplotlib theme application.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `palettes` and `themes` modules |
| `palettes.py` | Color palettes for chromosomes, expression, significance, and categories |
| `themes.py` | Matplotlib theme definitions (publication, presentation) with context manager |

## Key Functions

| Function | Description |
|----------|-------------|
| `palettes.chromosome_palette()` | Get colors for chromosome visualization |
| `palettes.categorical()` | Generate colorblind-safe categorical palette (Wong palette) |
| `palettes.expression_gradient()` | Create diverging colormap for expression data |
| `palettes.significance_palette()` | Map significance levels to colors |
| `palettes.significance_color()` | Get color for a specific p-value |
| `palettes.heatmap_cmap()` | Get named colormap for heatmaps |
| `palettes.alternating_pair()` | Generate alternating color pairs |
| `themes.list_themes()` | List all available theme names |
| `themes.get_theme()` | Retrieve a theme's parameter dictionary |
| `themes.apply_theme()` | Apply a theme globally to matplotlib |
| `themes.theme()` | Context manager for scoped theme application |
| `themes.register_theme()` | Register a custom theme |
| `themes.reset_theme()` | Reset matplotlib to default styling |

## Usage

```python
from metainformant.visualization.config import palettes, themes

colors = palettes.chromosome_palette()
with themes.theme("publication"):
    fig, ax = plt.subplots()
```
