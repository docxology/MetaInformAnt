# AI Agents in Single-Cell Genomics Development

This document outlines AI assistance in developing METAINFORMANT's single-cell RNA sequencing analysis capabilities.

## AI Contributions

### Single-Cell Architecture
**Code Assistant Agent** designed:
- Comprehensive single-cell analysis framework
- AnnData-compatible data structure
- Preprocessing pipeline organization
- Integration with scanpy ecosystem

### Analysis Components
**Code Assistant Agent** contributed to:
- Quality control and preprocessing utilities
- Dimensionality reduction algorithms
- Clustering and trajectory analysis
- Visualization and integration methods

### Quality Assurance
**Documentation Agent** assisted with:
- Single-cell analysis documentation
- API reference generation for preprocessing functions
- Usage examples and best practices
- Integration guides for scRNA-seq workflows

## Development Approach

- **Modular Design**: AI helped design flexible single-cell modules
- **scanpy Compatibility**: Established AnnData-compatible interfaces
- **Performance Optimization**: Efficient algorithms for large single-cell datasets
- **Extensibility**: Framework for adding new single-cell analysis methods

## Quality Assurance

- Human oversight ensures biological accuracy and relevance
- AI assistance accelerates development while maintaining standards
- Comprehensive testing validates single-cell functionality

This single-cell infrastructure provides a solid foundation for scRNA-seq analysis workflows.

## Complete Function Signatures

### Preprocessing (`preprocessing.py`)
- `load_count_matrix(filepath: str | Path, format: str = "h5ad", **kwargs) -> SingleCellData`
- `calculate_qc_metrics(data: SingleCellData) -> SingleCellData`
- `filter_cells(data: SingleCellData, min_counts: int | None = None, max_counts: int | None = None, min_genes: int | None = None, max_genes: int | None = None, max_mito_percent: float | None = None) -> SingleCellData`
- `filter_genes(data: SingleCellData, min_cells: int | None = None, max_cells: int | None = None) -> SingleCellData`
- `normalize_counts(data: SingleCellData, target_sum: float | None = None, normalize_method: str = "total") -> SingleCellData`
- `log_transform(data: SingleCellData, base: float = np.e) -> SingleCellData`
- `scale_data(data: SingleCellData, zero_center: bool = True, max_value: float | None = None) -> SingleCellData`

### Clustering (`clustering.py`)
- `leiden_clustering(data: SingleCellData, resolution: float = 1.0, n_neighbors: int = 15, random_state: int | None = None) -> SingleCellData`
- `louvain_clustering(data: SingleCellData, resolution: float = 1.0, n_neighbors: int = 15, random_state: int | None = None) -> SingleCellData`
- `kmeans_clustering(data: SingleCellData, n_clusters: int = 10, random_state: int | None = None) -> SingleCellData`
- `hierarchical_clustering(data: SingleCellData, n_clusters: int = 10, linkage: str = "ward") -> SingleCellData`
- `find_marker_genes(data: SingleCellData, groupby: str, method: str = "t-test", n_genes: int = 100) -> pd.DataFrame`
- `compute_cluster_composition(data: SingleCellData, groupby: str, cluster_col: str = "cluster") -> pd.DataFrame`
- `compute_cluster_silhouette(data: SingleCellData, cluster_col: str = "cluster") -> Dict[str, float]`

### Dimensionality Reduction (`dimensionality.py`)
- `pca_reduction(data: SingleCellData, n_components: int = 50, random_state: int | None = None) -> SingleCellData`
- `tsne_reduction(data: SingleCellData, n_components: int = 2, perplexity: float = 30.0, random_state: int | None = None) -> SingleCellData`
- `umap_reduction(data: SingleCellData, n_components: int = 2, n_neighbors: int = 15, random_state: int | None = None) -> SingleCellData`
- `diffusion_map_reduction(data: SingleCellData, n_components: int = 10) -> SingleCellData`

### Trajectory Inference (`trajectory.py`)
- `compute_diffusion_pseudotime(data: SingleCellData, root_cell: int | None = None) -> SingleCellData`
- `dpt_trajectory(data: SingleCellData, root_cell: int | None = None) -> SingleCellData`
- `paga_trajectory(data: SingleCellData, groups: str) -> SingleCellData`
- `slingshot_trajectory(data: SingleCellData, start_cluster: str, end_clusters: List[str]) -> SingleCellData`

### Integration (`integration.py`)
- `bbknn_integration(data: SingleCellData, batch_key: str, n_neighbors: int = 15) -> SingleCellData`
- `harmony_integration(data: SingleCellData, batch_key: str, n_components: int = 50) -> SingleCellData`
- `scanorama_integration(data_list: List[SingleCellData], batch_key: str) -> SingleCellData`
- `mnn_integration(data_list: List[SingleCellData], batch_key: str) -> SingleCellData`

### Visualization (`visualization.py`)
- `plot_umap(data: SingleCellData, color: str | None = None, **kwargs) -> matplotlib.figure.Figure`
- `plot_tsne(data: SingleCellData, color: str | None = None, **kwargs) -> matplotlib.figure.Figure`
- `plot_pca(data: SingleCellData, color: str | None = None, **kwargs) -> matplotlib.figure.Figure`
- `plot_trajectory(data: SingleCellData, trajectory_key: str, **kwargs) -> matplotlib.figure.Figure`
- `plot_marker_expression(data: SingleCellData, marker_genes: List[str], **kwargs) -> matplotlib.figure.Figure`

## Complete Function Signatures

### Preprocessing (`preprocessing.py`)
- `load_count_matrix(filepath: str | Path, format: str = "h5ad", **kwargs) -> SingleCellData`
- `calculate_qc_metrics(data: SingleCellData) -> SingleCellData`
- `filter_cells(data: SingleCellData, min_counts: int | None = None, max_counts: int | None = None, min_genes: int | None = None, max_genes: int | None = None, max_mito_percent: float | None = None) -> SingleCellData`
- `filter_genes(data: SingleCellData, min_cells: int | None = None, max_cells: int | None = None) -> SingleCellData`
- `normalize_counts(data: SingleCellData, target_sum: float | None = None, normalize_method: str = "total") -> SingleCellData`
- `log_transform(data: SingleCellData, base: float = np.e) -> SingleCellData`
- `scale_data(data: SingleCellData, zero_center: bool = True, max_value: float | None = None) -> SingleCellData`
- `identify_highly_variable_genes(data: SingleCellData, n_top_genes: int = 2000, flavor: str = "seurat") -> SingleCellData`
- `remove_batch_effects(data: SingleCellData, batch_key: str, method: str = "regress_out") -> SingleCellData`

### Clustering (`clustering.py`)
- `leiden_clustering(data: SingleCellData, resolution: float = 1.0, n_neighbors: int = 15, random_state: int | None = None, use_weights: bool = True) -> SingleCellData`
- `louvain_clustering(data: SingleCellData, resolution: float = 1.0, n_neighbors: int = 15, random_state: int | None = None, use_weights: bool = True) -> SingleCellData`
- `kmeans_clustering(data: SingleCellData, n_clusters: int = 10, random_state: int | None = None, n_init: int = 10) -> SingleCellData`
- `hierarchical_clustering(data: SingleCellData, n_clusters: int = 10, linkage_method: str = "ward", metric: str = "euclidean") -> SingleCellData`
- `find_marker_genes(data: SingleCellData, groupby: str, method: str = "t-test", n_genes: int = 100) -> pd.DataFrame`
- `compute_cluster_composition(data: SingleCellData, groupby: str, cluster_col: str = "cluster") -> pd.DataFrame`
- `compute_cluster_silhouette(data: SingleCellData, cluster_col: str = "cluster") -> Dict[str, float]`
- `evaluate_clustering_performance(data: SingleCellData, cluster_col: str = "cluster", ground_truth_col: Optional[str] = None) -> Dict[str, Any]`

### Dimensionality Reduction (`dimensionality.py`)
- `pca_reduction(data: SingleCellData, n_components: int = 50, random_state: int | None = None, scale_data: bool = True) -> SingleCellData`
- `tsne_reduction(data: SingleCellData, n_components: int = 2, perplexity: float = 30.0, random_state: int | None = None, learning_rate: float = 200.0, max_iter: int = 1000) -> SingleCellData`
- `umap_reduction(data: SingleCellData, n_components: int = 2, n_neighbors: int = 15, min_dist: float = 0.1, random_state: int | None = None, metric: str = "euclidean") -> SingleCellData`
- `diffusion_map_reduction(data: SingleCellData, n_components: int = 10, n_neighbors: int = 15, alpha: float = 1.0) -> SingleCellData`
- `mds_reduction(data: SingleCellData, n_components: int = 2, metric: bool = True, random_state: int | None = None) -> SingleCellData`
- `ica_reduction(data: SingleCellData, n_components: int = 10, random_state: int | None = None, max_iter: int = 1000) -> SingleCellData`
- `factor_analysis_reduction(data: SingleCellData, n_components: int = 10, random_state: int | None = None) -> SingleCellData`
- `compute_dimensionality_metrics(data: SingleCellData, embedding_cols: List[str]) -> Dict[str, Any]`

### Trajectory Inference (`trajectory.py`)
- `compute_diffusion_pseudotime(data: SingleCellData, root_cell: int | None = None, n_components: int = 10) -> SingleCellData`
- `dpt_trajectory(data: SingleCellData, root_cell: int | None = None) -> SingleCellData`
- `paga_trajectory(data: SingleCellData, groups: str) -> SingleCellData`
- `slingshot_trajectory(data: SingleCellData, start_cluster: str, end_clusters: List[str]) -> SingleCellData`
- `compute_pseudotime_from_dimensionality_reduction(data: SingleCellData, dim_red_cols: List[str], root_cell: int | None = None) -> SingleCellData`
- `find_trajectory_branches(data: SingleCellData, pseudotime_col: str, min_branch_size: int = 10) -> Dict[str, Any]`
- `compute_trajectory_entropy(data: SingleCellData, pseudotime_col: str, window_size: int = 100) -> Dict[str, Any]`

### Integration (`integration.py`)
- `bbknn_integration(data: SingleCellData, batch_key: str, n_neighbors: int = 15) -> SingleCellData`
- `harmony_integration(data: SingleCellData, batch_key: str, n_components: int = 50) -> SingleCellData`
- `scanorama_integration(data_list: List[SingleCellData], batch_key: str) -> SingleCellData`
- `mnn_integration(data_list: List[SingleCellData], batch_key: str) -> SingleCellData`
- `combat_integration(data: SingleCellData, batch_key: str, covariates: Optional[List[str]] = None) -> SingleCellData`
- `integrate_multiple_batches(data_list: List[SingleCellData], integration_method: str = "scanorama", **kwargs) -> SingleCellData`
- `evaluate_integration_quality(integrated_data: SingleCellData, batch_key: str, cluster_key: Optional[str] = None) -> Dict[str, Any]`

### Visualization (`visualization.py`)
- `plot_umap(data: SingleCellData, color: Optional[str] = None, **kwargs) -> Any`
- `plot_tsne(data: SingleCellData, color: Optional[str] = None, **kwargs) -> Any`
- `plot_pca(data: SingleCellData, color: Optional[str] = None, n_components: int = 2, **kwargs) -> Any`
- `plot_trajectory(data: SingleCellData, trajectory_key: str, color_by_pseudotime: bool = True, **kwargs) -> Any`
- `plot_marker_expression(data: SingleCellData, marker_genes: List[str], method: str = 'dotplot', **kwargs) -> Any`
- `plot_qc_metrics(data: SingleCellData, **kwargs) -> Any`
- `plot_cluster_comparison(data: SingleCellData, cluster_cols: List[str], embedding_cols: Optional[List[str]] = None, **kwargs) -> Any`
