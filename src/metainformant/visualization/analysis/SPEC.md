# SPEC: Statistical Analysis Visualization

## Standards

- **Color Palettes**: Uses color-blind friendly palettes by default (e.g., Viridis, Colorcet).
- **Labeling**: Every plot must have axis labels and a title if applicable.
- **Return Type**: Consistent return of `matplotlib.axes.Axes`.

## Implementation Details

- **Scaling**: Automatic scaling is handled by matplotlib, but can be overridden via `kwargs`.
- **Performance**: Large datasets (e.g., single-cell counts) are optimized using sub-sampling or hexbinning where appropriate.
