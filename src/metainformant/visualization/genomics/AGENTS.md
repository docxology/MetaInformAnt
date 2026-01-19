# AI Agents: Genomics Visualization

Function index and technical details for the `visualization.genomics` module.

## Function Index

### `metainformant.visualization.genomics.genomics`
- **`manhattan_plot(results: pd.DataFrame, significance_threshold: float = 5e-8) -> Axes`**: Generates a Manhattan plot.
- **`volcano_plot(results: pd.DataFrame, p_col: str, lfc_col: str) -> Axes`**: Generates a volcano plot.
- **`regional_plot(results: pd.DataFrame, chrom: str, start: int, end: int) -> Axes`**: Focused view of a specific genomic region.

### `metainformant.visualization.genomics.trees`
- **`plot_phylo_tree(tree: Any, title: str = "Phylogenetic Tree") -> Axes`**: Renders phylogenetic trees.
- **`circular_tree_plot(tree: Any) -> Axes`**: Circular representation of trees.

### `metainformant.visualization.genomics.expression`
- **`plot_expression_heatmap(counts_matrix: pd.DataFrame) -> Axes`**: Visualizes gene expression patterns.
