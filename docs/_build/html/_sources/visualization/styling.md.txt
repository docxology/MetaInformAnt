# Styling and Customization

Utilities for creating publication-quality figures with consistent styling, color palettes, and formatting.

## Style Presets

### `apply_publication_style(**overrides)`

Apply publication-quality style settings.

**Example:**
```python
from metainformant.visualization.style import apply_publication_style

apply_publication_style(font.size=12)
```

## Color Palettes

### `get_color_palette(name='primary')`

Get a color palette by name.

**Available palettes:**
- `'primary'`: Standard scientific colors
- `'accessible'`: Accessible color schemes
- `'colorblind'`: Colorblind-friendly palette
- `'viridis_like'`: Viridis-like sequential colors

**Example:**
```python
from metainformant.visualization.style import get_color_palette

colors = get_color_palette('colorblind')
```

## Figure Sizes

### `get_figure_size(size='medium')`

Get a figure size by name.

**Available sizes:**
- `'small'`: (4, 3)
- `'medium'`: (6, 4.5)
- `'large'`: (8, 6)
- `'wide'`: (10, 4)
- `'tall'`: (6, 8)
- `'square'`: (6, 6)
- `'presentation'`: (10, 7.5)
- `'poster'`: (12, 9)

## Font Management

### `set_font_family(family='sans-serif')`

Set the default font family.

### `set_font_size(size=10)`

Set the default font size.

### `reset_style()`

Reset matplotlib style to default.

