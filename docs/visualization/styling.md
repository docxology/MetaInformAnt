# Styling and Customization

Utilities for creating publication-quality figures with consistent styling, color palettes, and formatting.
Themes live in `metainformant.visualization.config.themes` and color palettes in
`metainformant.visualization.config.palettes`.

## Style Presets

### `apply_theme(name)`

Apply a named theme globally (it persists until `reset_theme()` is called or another theme is applied).

**Available themes:** `'publication'`, `'presentation'`, `'poster'`, `'dark'`, `'minimal'`, `'colorblind'`.
`list_themes()` returns the registered names and `get_theme(name)` returns a theme's matplotlib
rc-parameters.

**Example:**
```python-snippet
from metainformant.visualization.config.themes import apply_theme

apply_theme("publication")
```

To override individual settings, register a modified copy of a theme and apply it:

```python-snippet
from metainformant.visualization.config.themes import apply_theme, get_theme, register_theme

register_theme("publication-12pt", {**get_theme("publication"), "font.size": 12})
apply_theme("publication-12pt")
```

## Color Palettes

### `categorical(n, *, palette='wong')`

Get `n` colorblind-safe categorical colors as a list of hex strings.

**Available palettes:**
- `'wong'`: Wong (2011) colorblind-safe palette, 8 colors (default)
- `'tol'`: Tol bright palette, 7 colors
- `'ibm'`: IBM colorblind-safe palette, 5 colors

For `n` larger than the chosen palette, colors fall back to matplotlib's `tab20`.

**Example:**
```python-snippet
from metainformant.visualization.config.palettes import categorical

colors = categorical(5, palette='wong')
```

The same module provides `chromosome_palette`, `expression_gradient`, `heatmap_cmap`,
`significance_palette` / `significance_color` and `alternating_pair` for common bioinformatics uses.

## Figure Sizes

Figure size is part of each theme (the `figure.figsize` rc-parameter):

- `'publication'`: (7, 5)
- `'presentation'`: (10, 7)
- `'poster'`: (14, 10)
- `'dark'`: (8, 6)
- `'minimal'`: (7, 5)
- `'colorblind'`: (7, 5)

## Font Management

Fonts are set through a theme's rc-parameters (`font.family`, `font.size`, `axes.titlesize`, ...).

### `register_theme(name, params)`

Register a custom theme from a dictionary of matplotlib rc-parameters (see the example above).

### `theme(name)`

Context manager that applies a theme only inside a `with` block and restores the previous settings afterwards.

### `reset_theme()`

Reset matplotlib's rc-parameters to their defaults.
