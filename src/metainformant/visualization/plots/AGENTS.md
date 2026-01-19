# AI Agents: Plotting and Visualizations

Function index and technical details for the `visualization.plots` module.

## Function Index

### `metainformant.visualization.plots.basic`
- **`lineplot(x: np.ndarray, y: np.ndarray) -> Axes`**: Renders a line plot.
- **`scatter_plot(x: np.ndarray, y: np.ndarray) -> Axes`**: Renders a scatter plot.
- **`heatmap(data: np.ndarray) -> Axes`**: Renders a heatmap.

### `metainformant.visualization.plots.specialized`
- **`plot_venn_diagram(sets: dict[str, set]) -> Axes`**: Generates Venn diagrams.
- **`plot_sankey_diagram(data: pd.DataFrame) -> Axes`**: Generates Sankey diagrams.

### `metainformant.visualization.plots.animations`
- **`animate_time_series(data: np.ndarray) -> tuple[Figure, FuncAnimation]`**: Animates time-series data.
- **`animate_evolution(generations: list[np.ndarray]) -> tuple[Figure, FuncAnimation]`**: Animates evolutionary changes.
