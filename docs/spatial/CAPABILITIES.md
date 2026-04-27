# Spatial Module Capabilities

## Quick Reference

### Functions by Submodule

| Function | Submodule | Purpose |
|----------|-----------|---------|
| `gearys_c()` | `analysis.autocorrelation` | Compute Geary's C spatial autocorrelation statistic. |
| `getis_ord_g()` | `analysis.autocorrelation` | Compute Getis-Ord G* statistic for hot/cold spot detection. |
| `local_morans_i()` | `analysis.autocorrelation` | Compute Local Moran's I (LISA) statistics. |
| `morans_i()` | `analysis.autocorrelation` | Compute Moran's I spatial autocorrelation statistic. |
| `spatial_variogram()` | `analysis.autocorrelation` | Compute the empirical spatial variogram (semivariogram). |
| `spatial_weights_matrix()` | `analysis.autocorrelation` | Build a spatial weights matrix from coordinates. |
| `_compute_modularity()` | `analysis.clustering` | Compute Newman-Girvan modularity Q for a partition. |
| `_estimate_n_clusters()` | `analysis.clustering` | Estimate optimal number of clusters using the elbow method with spatial features. |
| `_greedy_modularity_clustering()` | `analysis.clustering` | Greedy modularity maximization as fallback clustering. |
| `build_spatial_graph()` | `analysis.clustering` | Build a spatial neighborhood graph from coordinates. |
| `leiden_clustering()` | `analysis.clustering` | Leiden community detection on a graph. |
| `louvain_clustering()` | `analysis.clustering` | Louvain community detection on a graph. |
| `spatial_cluster()` | `analysis.clustering` | Spatially-aware clustering combining expression and spatial proximity. |
| `spatial_domains()` | `analysis.clustering` | Identify spatial domains using a BayesSpace-inspired iterative approach. |
| `_nmf_deconvolution()` | `analysis.deconvolution` | NMF-based deconvolution. |
| `create_reference_profiles()` | `analysis.deconvolution` | Build cell type reference expression profiles from scRNA-seq data. |
| `deconvolve_spots()` | `analysis.deconvolution` | Deconvolve spatial spots to estimate cell type composition. |
| `enrichment_score()` | `analysis.deconvolution` | Compute cell type enrichment score per spot. |
| `estimate_cell_fractions()` | `analysis.deconvolution` | Normalize deconvolution weights to cell type fractions (proportions summing to 1). |
| `nnls_deconvolution()` | `analysis.deconvolution` | Non-negative least squares deconvolution for a single spot or batch. |
| `compute_interaction_matrix()` | `analysis.neighborhood` | Compute pairwise cell type interaction scores from a spatial graph. |
| `ligand_receptor_spatial()` | `analysis.neighborhood` | Spatial ligand-receptor interaction analysis. |
| `neighborhood_enrichment()` | `analysis.neighborhood` | Compute cell type neighborhood enrichment (co-localization analysis). |
| `niche_detection()` | `analysis.neighborhood` | Identify cellular niches based on local cell type composition. |
| `ripley_k()` | `analysis.neighborhood` | Compute Ripley's K function for spatial point pattern analysis. |
| `_gene_name_to_index()` | `communication.cell_communication` | Convert gene name to index. |
| `_pairwise_euclidean()` | `communication.cell_communication` | Compute pairwise Euclidean distance matrix without scipy. |
| `_simple_nmf()` | `communication.cell_communication` | Simple multiplicative update NMF implementation. |
| `build_communication_network()` | `communication.cell_communication` | Build cell-cell communication network from interaction results. |
| `communication_pattern_analysis()` | `communication.cell_communication` | Identify communication patterns using NMF on interaction matrix. |
| `compute_ligand_receptor_interactions()` | `communication.cell_communication` | Compute ligand-receptor interaction scores between cell types. |
| `default_lr_database()` | `communication.cell_communication` | Return built-in ligand-receptor pair database. |
| `spatial_interaction_score()` | `communication.cell_communication` | Score cell-cell communication considering spatial proximity. |
| `_kmeans()` | `deconvolution.spatial_deconvolution` | Simple k-means clustering implementation. |
| `_normalize_rows()` | `deconvolution.spatial_deconvolution` | Normalize each row of a matrix to sum to 1. |
| `_pairwise_euclidean()` | `deconvolution.spatial_deconvolution` | Compute pairwise Euclidean distance matrix without scipy. |
| `_regression_deconvolution()` | `deconvolution.spatial_deconvolution` | Ordinary least squares deconvolution with non-negativity clipping. |
| `build_reference_profiles()` | `deconvolution.spatial_deconvolution` | Build cell-type reference profiles from single-cell data. |
| `deconvolve_spots()` | `deconvolution.spatial_deconvolution` | Deconvolve spatial spots into cell type proportions. |
| `niche_identification()` | `deconvolution.spatial_deconvolution` | Identify tissue niches from cell type composition and spatial proximity. |
| `spatial_cell_type_mapping()` | `deconvolution.spatial_deconvolution` | Map cell type proportions to spatial coordinates for visualization. |
| `validate_deconvolution()` | `deconvolution.spatial_deconvolution` | Validate spatial deconvolution using marker gene expression agreement. |
| `_intersect_genes()` | `integration.scrna_mapping` | Find shared genes between spatial and reference datasets. |
| `_to_dense()` | `integration.scrna_mapping` | Convert sparse matrix to dense numpy array. |
| `anchor_based_transfer()` | `integration.scrna_mapping` | Anchor-based label transfer from scRNA-seq to spatial data. |
| `correlation_mapping()` | `integration.scrna_mapping` | Correlation-based mapping of spatial spots to cell type profiles. |
| `impute_spatial_genes()` | `integration.scrna_mapping` | Impute unmeasured genes in spatial data using scRNA-seq reference. |
| `map_scrna_to_spatial()` | `integration.scrna_mapping` | Map scRNA-seq cell type annotations to spatial spots. |
| `_find_column()` | `io.merfish` | Find the first matching column name from a list of candidates. |
| `aggregate_to_cells()` | `io.merfish` | Aggregate individual transcript spots to cell-level expression counts. |
| `load_merfish()` | `io.merfish` | Load a complete MERFISH dataset from a directory. |
| `load_transcript_spots()` | `io.merfish` | Load individual transcript spot coordinates from a MERFISH detected transcripts file. |
| `parse_cell_metadata()` | `io.merfish` | Parse MERFISH cell metadata CSV file. |
| `_read_h5_matrix()` | `io.visium` | Read a 10x HDF5 filtered feature-barcode matrix. |
| `_read_mex_matrix()` | `io.visium` | Read a Market Exchange (MEX) format matrix directory. |
| `_read_scale_factors()` | `io.visium` | Read Visium scale factors JSON. |
| `create_spatial_dataset()` | `io.visium` | Create a unified SpatialDataset from components. |
| `filter_tissue_spots()` | `io.visium` | Filter tissue position records to keep only spots overlapping tissue. |
| `load_visium()` | `io.visium` | Load a 10x Visium Spatial Gene Expression dataset. |
| `read_spatial_image()` | `io.visium` | Load an H&E tissue image as a numpy array. |
| `read_tissue_positions()` | `io.visium` | Parse a Visium tissue_positions.csv file. |
| `_find_col()` | `io.xenium` | Find first matching column name from candidates. |
| `_load_boundaries_parquet()` | `io.xenium` | Load cell boundaries from a parquet file using pandas. |
| `_open_maybe_gzipped()` | `io.xenium` | Open a file that may or may not be gzipped. |
| `_read_xenium_h5()` | `io.xenium` | Read Xenium HDF5 feature matrix. |
| `_read_xenium_mex()` | `io.xenium` | Read Xenium MEX format feature matrix. |
| `load_cell_boundaries()` | `io.xenium` | Load cell segmentation polygon boundaries. |
| `load_xenium()` | `io.xenium` | Load a complete 10x Xenium dataset from a directory. |
| `read_cell_features()` | `io.xenium` | Read a Xenium cell-level feature matrix. |
| `read_transcripts()` | `io.xenium` | Read per-transcript coordinates from a Xenium transcripts file. |
| `_kmeans()` | `niche.identification` | Simple K-Means clustering. |
| `_spatial_smooth()` | `niche.identification` | Smooth values over spatial neighbors. |
| `identify_niches()` | `niche.identification` | Identify tissue niches from local cell type composition. |
| `_ensure_plotting_deps()` | `visualization.plots` | Check that matplotlib is available. |
| `_save_figure()` | `visualization.plots` | Save a matplotlib figure to disk, creating directories as needed. |
| `plot_cell_type_map()` | `visualization.plots` | Plot spatial distribution of cell types. |
| `plot_deconvolution_pie()` | `visualization.plots` | Plot pie charts per spatial spot showing cell type fractions. |
| `plot_gene_expression_map()` | `visualization.plots` | Plot spatial expression map for a single gene. |
| `plot_interaction_heatmap()` | `visualization.plots` | Plot a cell type interaction heatmap. |
| `plot_neighborhood_graph()` | `visualization.plots` | Plot the spatial neighborhood graph on tissue coordinates. |
| `plot_spatial_autocorrelation()` | `visualization.plots` | Plot LISA cluster map showing spatial autocorrelation patterns. |
| `plot_spatial_scatter()` | `visualization.plots` | Create a spatial scatter plot colored by continuous or categorical values. |
| `plot_tissue_overlay()` | `visualization.plots` | Overlay expression values on a tissue H&E image. |

---

## analysis.autocorrelation

### Functions

#### `spatial_weights_matrix()`

**Signature**: `spatial_weights_matrix(coordinates, method, k, bandwidth)`

Build a spatial weights matrix from coordinates.

Args:
    coordinates: Spatial coordinates (n x 2).
    method: Weights construction method:
        - "knn": K-nearest neighbors with binary weights.
        - "distance": Inverse distance weighting within bandwidth.
        - "binary": Binary weights within bandwidth.
    k: Number of neighbors for KNN method.
    bandwidth: Distance threshold for distance/binary methods.
        If None, uses median nearest-neighbor distance * 1.5.
    row_standardize: If True, normalize rows to sum to 1.

Returns:
    Sparse weights matrix (n x n) in CSR format.

#### `morans_i()`

**Signature**: `morans_i(values, weights)`

Compute Moran's I spatial autocorrelation statistic.

Moran's I measures the overall spatial autocorrelation of a variable.
It is the spatial analog of Pearson's correlation coefficient.

Formula:
    I = (n / S0) * (z^T W z) / (z^T z)
where z = x - mean(x), S0 = sum of all weights.

Args:
    values: Observed values (length n).
    weights: Spatial weights matrix (n x n), sparse or dense.

Returns:
    MoransIResult with I statistic, expected value, variance, z-score, p-value.

#### `gearys_c()`

**Signature**: `gearys_c(values, weights)`

Compute Geary's C spatial autocorrelation statistic.

Geary's C is based on squared differences between neighboring observations,
making it more sensitive to local spatial patterns than Moran's I.

Formula:
    C = ((n-1) / (2 * S0)) * (sum_ij w_ij (x_i - x_j)^2) / (sum_i (x_i - x_bar)^2)

Args:
    values: Observed values (length n).
    weights: Spatial weights matrix (n x n).

Returns:
    GearyCResult with C statistic, z-score, and p-value.

#### `local_morans_i()`

**Signature**: `local_morans_i(values, weights)`

Compute Local Moran's I (LISA) statistics.

Local Indicators of Spatial Association (LISA) decompose Moran's I into
per-observation contributions, identifying local clusters and spatial outliers.

Formula for observation i:
    I_i = z_i * sum_j(w_ij * z_j) / (sum z^2 / n)

Classification:
    - HH: High value surrounded by high values (hot spot).
    - LL: Low value surrounded by low values (cold spot).
    - HL: High value surrounded by low values (outlier).
    - LH: Low value surrounded by high values (outlier).
    - NS: Not significant at the given significance level.

Args:
    values: Observed values (length n).
    weights: Spatial weights matrix (n x n). Should be row-standardized.
    significance: Significance level for cluster classification.

Returns:
    LocalMoransResult with local I values, z-scores, p-values, and cluster labels.

#### `getis_ord_g()`

**Signature**: `getis_ord_g(values, weights)`

Compute Getis-Ord G* statistic for hot/cold spot detection.

G* measures the concentration of high or low values in the neighborhood
of each observation. Unlike LISA, G* includes the observation itself.

Formula:
    G*_i = (sum_j w_ij * x_j - x_bar * sum_j w_ij) /
           (S * sqrt((n * sum_j w_ij^2 - (sum_j w_ij)^2) / (n-1)))

where S is the standard deviation of all values.

Args:
    values: Observed values (length n).
    weights: Spatial weights matrix (n x n). Should NOT be row-standardized.
    significance: Significance level for hot/cold spot identification.

Returns:
    GetisOrdResult with G* statistics and hot/cold spot identification.

#### `spatial_variogram()`

**Signature**: `spatial_variogram(values, coordinates, n_bins)`

Compute the empirical spatial variogram (semivariogram).

The semivariogram gamma(h) measures the average squared difference
between pairs of observations separated by distance h:

    gamma(h) = (1 / 2|N(h)|) * sum_{(i,j) in N(h)} (x_i - x_j)^2

where N(h) is the set of pairs at distance h (within a bin).

Also estimates variogram model parameters:
- Nugget: discontinuity at the origin (measurement error + micro-scale variation).
- Sill: the plateau level (total variance).
- Range: distance at which the sill is reached.

Args:
    values: Observed values (length n).
    coordinates: Spatial coordinates (n x 2).
    n_bins: Number of distance bins.
    max_distance: Maximum distance to consider. If None, uses half the
        maximum pairwise distance.

Returns:
    VariogramResult with semivariance values and estimated parameters.

### Classes

#### `MoransIResult`

Result of Moran's I spatial autocorrelation test.

Attributes:
    I: Moran's I statistic. Range roughly [-1, 1].
        Positive = positive spatial autocorrelation (similar values cluster).
        Near 0 = random spatial pattern.
        Negative = negative spatial autocorrelation (checkerboard pattern).
    expected_I: Expected I under null hypothesis of no autocorrelation.
    variance_I: Variance of I under normality assumption.
    z_score: Standardized Z-score: (I - E[I]) / sqrt(Var[I]).
    p_value: Two-sided p-value from normal approximation.
    n: Number of observations.

#### `GearyCResult`

Result of Geary's C spatial autocorrelation test.

Attributes:
    C: Geary's C statistic. Range [0, 2].
        C < 1 = positive spatial autocorrelation.
        C = 1 = no spatial autocorrelation.
        C > 1 = negative spatial autocorrelation.
    expected_C: Expected C under null hypothesis (= 1).
    variance_C: Variance of C under normality assumption.
    z_score: Standardized Z-score.
    p_value: Two-sided p-value.
    n: Number of observations.

#### `LocalMoransResult`

Result of Local Moran's I (LISA) analysis.

Attributes:
    local_I: Local Moran's I values per observation (length n).
    expected_I: Expected local I under null.
    z_scores: Z-scores per observation.
    p_values: P-values per observation.
    cluster_labels: LISA cluster classification per observation:
        "HH" = High-High (hot spot), "LL" = Low-Low (cold spot),
        "HL" = High-Low (spatial outlier), "LH" = Low-High (spatial outlier),
        "NS" = Not significant.
    significance_level: Alpha level used for classification.

#### `GetisOrdResult`

Result of Getis-Ord G* statistic analysis.

Attributes:
    g_star: G* statistic per observation (length n). Z-score scale.
    p_values: P-values per observation.
    hot_spots: Boolean mask for significant hot spots.
    cold_spots: Boolean mask for significant cold spots.
    significance_level: Alpha level used.

#### `VariogramResult`

Result of spatial variogram/semivariogram analysis.

Attributes:
    bin_centers: Distance bin centers (length n_bins).
    semivariance: Estimated semivariance at each bin.
    n_pairs: Number of point pairs in each bin.
    nugget: Estimated nugget (semivariance at distance 0).
    sill: Estimated sill (plateau semivariance).
    range_param: Estimated range (distance at which sill is reached).

## analysis.clustering

### Functions

#### `build_spatial_graph()`

**Signature**: `build_spatial_graph(coordinates, method, n_neighbors, radius)`

Build a spatial neighborhood graph from coordinates.

Args:
    coordinates: Array of shape (n, 2) with spatial coordinates.
    method: Graph construction method:
        - "knn": K-nearest neighbors graph.
        - "delaunay": Delaunay triangulation graph.
        - "radius": Fixed-radius neighborhood graph.
    n_neighbors: Number of neighbors for KNN (default 6 for hexagonal Visium grid).
    radius: Radius for radius-based graph (required if method="radius").

Returns:
    Sparse adjacency matrix (scipy CSR) of shape (n, n).

Raises:
    ImportError: If scipy or sklearn is not installed.
    ValueError: If method is "radius" but no radius is provided.

#### `leiden_clustering()`

**Signature**: `leiden_clustering(adjacency_matrix, resolution, n_iterations, seed)`

Leiden community detection on a graph.

Implements the Leiden algorithm using a modularity-optimization approach.
Falls back to a greedy modularity maximization if the leidenalg package
is not available.

Args:
    adjacency_matrix: Sparse adjacency matrix (n x n).
    resolution: Resolution parameter (higher = more clusters).
    n_iterations: Number of iterations (-1 for until convergence).
    seed: Random seed for reproducibility.

Returns:
    Tuple of (labels, modularity) where labels is an integer array
    and modularity is the partition quality score.

Raises:
    ImportError: If neither leidenalg nor required fallback deps are available.

#### `louvain_clustering()`

**Signature**: `louvain_clustering(adjacency_matrix, resolution, seed)`

Louvain community detection on a graph.

Implements the Louvain method for modularity-based community detection.
Uses community_louvain if available, otherwise falls back to a greedy approach.

Args:
    adjacency_matrix: Sparse adjacency matrix (n x n).
    resolution: Resolution parameter (higher = more clusters).
    seed: Random seed.

Returns:
    Tuple of (labels, modularity).

#### `_greedy_modularity_clustering()`

**Signature**: `_greedy_modularity_clustering(adj, resolution, seed, max_iter)`

Greedy modularity maximization as fallback clustering.

Implements a simplified agglomerative modularity optimization:
1. Start with each node in its own community.
2. Repeatedly move nodes to the community that maximizes modularity gain.
3. Stop when no improvement is found.

This is a simplified version of the Louvain first phase.

Args:
    adj: Sparse adjacency matrix.
    resolution: Resolution parameter.
    seed: Random seed.
    max_iter: Maximum iterations.

Returns:
    Tuple of (labels, modularity).

#### `_compute_modularity()`

**Signature**: `_compute_modularity(adj, labels, m)`

Compute Newman-Girvan modularity Q for a partition.

Q = (1/2m) * sum_ij [ A_ij - k_i*k_j/(2m) ] * delta(c_i, c_j)

Args:
    adj: Sparse adjacency matrix.
    labels: Cluster label array.
    m: Total edge weight (sum of all edges).

Returns:
    Modularity score in [-1, 1].

#### `spatial_cluster()`

**Signature**: `spatial_cluster(expression, coordinates, n_clusters, method)`

Spatially-aware clustering combining expression and spatial proximity.

Constructs a joint graph that balances transcriptomic similarity with spatial
proximity, then applies community detection or KMeans.

For graph-based methods (leiden, louvain):
    1. Build expression similarity graph (KNN in PCA space).
    2. Build spatial proximity graph.
    3. Combine: A_combined = (1-w)*A_expr + w*A_spatial.
    4. Run community detection on combined graph.

For KMeans:
    1. PCA on expression.
    2. Concatenate scaled PCA coordinates with scaled spatial coordinates.
    3. Run KMeans on concatenated features.

Args:
    expression: Expression matrix (n_spots x n_genes), dense or sparse.
    coordinates: Spatial coordinates (n_spots x 2).
    n_clusters: Number of clusters (required for KMeans, optional for graph methods).
    method: Clustering algorithm ("leiden", "louvain", "kmeans").
    n_neighbors: Number of spatial neighbors.
    graph_method: How to build the spatial graph.
    resolution: Resolution parameter for graph methods.
    spatial_weight: Weight for spatial graph in [0, 1]. 0 = expression only, 1 = spatial only.
    n_pcs: Number of PCA components for expression.
    seed: Random seed.

Returns:
    SpatialClusterResult with labels and metadata.

#### `spatial_domains()`

**Signature**: `spatial_domains(expression, coordinates, n_domains)`

Identify spatial domains using a BayesSpace-inspired iterative approach.

Algorithm:
1. PCA on expression data.
2. Initialize domains via KMeans on PCA + spatial features.
3. Iteratively refine assignments by incorporating spatial neighbors:
   a. For each spot, compute mean expression profile of its spatial neighbors.
   b. Blend spot's own profile with its neighborhood mean.
   c. Re-cluster the blended profiles.
4. Repeat until convergence or max_iterations.

This is a simplified version of the BayesSpace spatial smoothing approach,
using iterative neighborhood averaging instead of full MCMC.

Args:
    expression: Expression matrix (n x genes), dense or sparse.
    coordinates: Spatial coordinates (n x 2).
    n_domains: Number of spatial domains. If None, auto-selects via gap statistic heuristic.
    n_pcs: Number of PCA components.
    n_neighbors: Spatial neighbors for smoothing.
    max_iterations: Maximum refinement iterations.
    seed: Random seed.

Returns:
    SpatialClusterResult with domain labels.

#### `_estimate_n_clusters()`

**Signature**: `_estimate_n_clusters(expr_pca, coords, max_k, seed)`

Estimate optimal number of clusters using the elbow method with spatial features.

Uses the second derivative of the within-cluster sum of squares (WCSS)
curve to find the elbow point.

Args:
    expr_pca: PCA-transformed expression (n x pcs).
    coords: Spatial coordinates (n x 2).
    max_k: Maximum number of clusters to try.
    seed: Random seed.

Returns:
    Estimated optimal number of clusters (minimum 2).

### Classes

#### `SpatialClusterResult`

Result of spatial clustering.

Attributes:
    labels: Cluster label per spot/cell (integer array of length n).
    n_clusters: Number of clusters found.
    method: Clustering method used.
    modularity: Modularity score (for graph-based methods).
    spatial_graph: Adjacency matrix used for clustering.
    metadata: Additional result metadata.

## analysis.deconvolution

### Functions

#### `create_reference_profiles()`

**Signature**: `create_reference_profiles(scrna_data, cell_type_labels)`

Build cell type reference expression profiles from scRNA-seq data.

For each cell type, computes the average (or median) expression profile
across all cells of that type.

Args:
    scrna_data: scRNA-seq expression matrix (n_cells x n_genes), dense or sparse.
    cell_type_labels: Array or list of cell type labels (length n_cells).
    gene_names: List of gene names (length n_genes). If None, uses indices.
    method: Aggregation method ("mean" or "median").

Returns:
    Tuple of (reference_profiles, cell_type_names, gene_names) where
    reference_profiles is (n_types x n_genes).

#### `nnls_deconvolution()`

**Signature**: `nnls_deconvolution(bulk_expression, reference_signatures)`

Non-negative least squares deconvolution for a single spot or batch.

Solves: min ||bulk - reference^T * x||^2  subject to x >= 0

For each spot, finds the non-negative weights that best reconstruct the
observed expression from the reference cell type signatures.

Args:
    bulk_expression: Expression vector(s). If 1D (n_genes,), single spot.
        If 2D (n_spots x n_genes), batch of spots.
    reference_signatures: Reference profiles matrix (n_types x n_genes).

Returns:
    Tuple of (weights, residuals). weights has shape matching input:
    (n_types,) for single spot or (n_spots, n_types) for batch.
    residuals is the per-spot/vector fitting residual norm.

Raises:
    ImportError: If scipy is not installed.

#### `estimate_cell_fractions()`

**Signature**: `estimate_cell_fractions(deconvolution_result)`

Normalize deconvolution weights to cell type fractions (proportions summing to 1).

Args:
    deconvolution_result: Either a DeconvolutionResult object or a raw weights
        matrix (n_spots x n_types).

Returns:
    Normalized fractions array (n_spots x n_types) where each row sums to 1.

#### `enrichment_score()`

**Signature**: `enrichment_score(observed, expected)`

Compute cell type enrichment score per spot.

Enrichment is the log2 fold change of observed fractions vs expected
(background) fractions, with a pseudocount for numerical stability.

enrichment = log2((observed + epsilon) / (expected + epsilon))

Args:
    observed: Observed cell type fractions (n_spots x n_types) or (n_types,).
    expected: Expected (background) cell type fractions (n_types,).

Returns:
    Enrichment scores with same shape as observed.

#### `deconvolve_spots()`

**Signature**: `deconvolve_spots(spatial_expression, reference_profiles, method)`

Deconvolve spatial spots to estimate cell type composition.

Supports multiple deconvolution methods:
- "nnls": Non-negative least squares (fast, robust).
- "nmf": Non-negative matrix factorization based approach.

When gene_names are provided for both spatial and reference data, the function
automatically intersects to shared genes.

Args:
    spatial_expression: Spatial expression matrix (n_spots x n_genes), dense or sparse.
    reference_profiles: Reference cell type profiles (n_types x n_genes).
    method: Deconvolution method ("nnls" or "nmf").
    cell_type_names: Names of cell types (length n_types).
    gene_names: Deprecated, use spatial_gene_names and reference_gene_names.
    spatial_gene_names: Gene names for spatial data.
    reference_gene_names: Gene names for reference profiles.
    alpha: Regularization parameter (for NMF).

Returns:
    DeconvolutionResult with weights, fractions, and residuals.

#### `_nmf_deconvolution()`

**Signature**: `_nmf_deconvolution(spatial, reference, alpha)`

NMF-based deconvolution.

Uses the reference profiles as a fixed basis (W) and solves for the
coefficient matrix (H) that reconstructs the spatial expression.

spatial ~= H @ reference (where H is n_spots x n_types)

Args:
    spatial: Expression matrix (n_spots x n_genes).
    reference: Reference profiles (n_types x n_genes).
    alpha: Regularization strength.

Returns:
    Tuple of (weights, residuals).

### Classes

#### `DeconvolutionResult`

Result of cell type deconvolution.

Attributes:
    weights: Raw deconvolution weights matrix (n_spots x n_types).
    fractions: Normalized cell type fractions (n_spots x n_types), rows sum to 1.
    cell_type_names: List of cell type names.
    residuals: Per-spot fitting residuals.
    method: Deconvolution method used.
    metadata: Additional result metadata.

**Methods**:

| Method | Purpose |
|--------|---------|
| `n_spots()` | Number of spatial spots. |
| `n_types()` | Number of cell types. |

## analysis.neighborhood

### Functions

#### `neighborhood_enrichment()`

**Signature**: `neighborhood_enrichment(cell_types, coordinates, radius)`

Compute cell type neighborhood enrichment (co-localization analysis).

For each pair of cell types (A, B), counts how often type-A cells
are neighbors of type-B cells compared to a random permutation baseline.
Returns Z-score enrichment: positive means co-localized, negative means avoided.

Algorithm:
1. Build spatial neighbor graph (KNN or radius-based).
2. Count observed pairwise type interactions.
3. Permute cell type labels N times to build null distribution.
4. Compute Z-score: (observed - mean_null) / std_null.

Args:
    cell_types: Array of cell type labels (length n).
    coordinates: Spatial coordinates (n x 2).
    radius: If specified, use radius-based neighbors instead of KNN.
    n_neighbors: Number of neighbors for KNN (ignored if radius is set).
    n_permutations: Number of permutations for null distribution.
    seed: Random seed.

Returns:
    NeighborhoodEnrichmentResult with enrichment matrix and statistics.

#### `compute_interaction_matrix()`

**Signature**: `compute_interaction_matrix(cell_types, spatial_graph)`

Compute pairwise cell type interaction scores from a spatial graph.

For each pair of cell types, computes the interaction score as the number
of edges between them, optionally normalized by the product of their
frequencies (to account for abundance).

Args:
    cell_types: Array of cell type labels (length n).
    spatial_graph: Sparse adjacency matrix (n x n).
    normalize_by_type_frequency: If True, normalize by expected frequency.

Returns:
    InteractionResult with interaction matrix.

#### `ligand_receptor_spatial()`

**Signature**: `ligand_receptor_spatial(expression, lr_pairs, coordinates)`

Spatial ligand-receptor interaction analysis.

For each ligand-receptor pair, computes a spatial interaction score
that measures the co-expression of the ligand in one spot and the
receptor in neighboring spots.

Score for pair (L, R) = mean over all spots i of:
    expression(L, i) * mean(expression(R, neighbors(i)))

Args:
    expression: Expression matrix (n_spots x n_genes), dense or sparse.
    lr_pairs: List of (ligand_gene, receptor_gene) tuples.
    coordinates: Spatial coordinates (n_spots x 2).
    gene_names: Gene name list (length n_genes) to map names to columns.
    radius: Radius for neighbor definition (if None, uses KNN).
    n_neighbors: Number of neighbors for KNN.

Returns:
    Dictionary with keys:
        - "scores": Dict mapping (ligand, receptor) -> float interaction score.
        - "per_spot_scores": Dict mapping (ligand, receptor) -> array of per-spot scores.
        - "n_pairs_tested": Number of pairs with both genes present.

#### `niche_detection()`

**Signature**: `niche_detection(cell_types, coordinates, n_niches)`

Identify cellular niches based on local cell type composition.

Algorithm:
1. Build spatial neighbor graph.
2. For each cell, compute the cell type composition of its neighborhood.
3. Cluster these neighborhood composition vectors to define niches.

Args:
    cell_types: Array of cell type labels (length n).
    coordinates: Spatial coordinates (n x 2).
    n_niches: Number of niches to identify.
    n_neighbors: Number of spatial neighbors to consider.
    seed: Random seed.

Returns:
    NicheResult with niche labels and compositions.

#### `ripley_k()`

**Signature**: `ripley_k(points, radii, area)`

Compute Ripley's K function for spatial point pattern analysis.

Ripley's K(r) counts the expected number of points within distance r
of a typical point, normalized by intensity. Under Complete Spatial
Randomness (CSR), K(r) = pi * r^2.

Also computes Besag's L-function: L(r) = sqrt(K(r)/pi) - r,
which is centered at 0 under CSR. Positive L(r) indicates clustering,
negative indicates regularity/inhibition.

Uses Ripley's isotropic edge correction.

Args:
    points: Point coordinates (n x 2).
    radii: Array of radii at which to evaluate K.
    area: Total study area (e.g., bounding box area).
    n_simulations: Number of CSR simulations for confidence envelope.
    seed: Random seed.

Returns:
    RipleyKResult with K values, L values, and CSR envelope.

### Classes

#### `NeighborhoodEnrichmentResult`

Result of neighborhood enrichment analysis.

Attributes:
    enrichment_matrix: Enrichment Z-scores (n_types x n_types).
        Positive = co-localized more than expected, negative = avoided.
    count_matrix: Observed interaction counts (n_types x n_types).
    expected_matrix: Expected interaction counts under random spatial arrangement.
    p_values: P-values from permutation testing (n_types x n_types).
    cell_type_names: List of cell type names.

#### `InteractionResult`

Result of pairwise cell type interaction analysis.

Attributes:
    interaction_matrix: Interaction scores (n_types x n_types).
    cell_type_names: List of cell type names.
    method: Scoring method used.

#### `NicheResult`

Result of cellular niche detection.

Attributes:
    niche_labels: Niche assignment per cell/spot (length n).
    niche_compositions: Cell type composition per niche (n_niches x n_types).
    n_niches: Number of niches found.
    cell_type_names: List of cell type names.

#### `RipleyKResult`

Result of Ripley's K function analysis.

Attributes:
    radii: Array of evaluation radii.
    k_values: K(r) values at each radius.
    l_values: L(r) = sqrt(K(r)/pi) - r (Besag's L-function, centered).
    csr_envelope_lower: Lower bound of CSR envelope (from simulations).
    csr_envelope_upper: Upper bound of CSR envelope.
    n_points: Number of points analyzed.
    area: Study area.

## communication.cell_communication

### Functions

#### `compute_ligand_receptor_interactions()`

**Signature**: `compute_ligand_receptor_interactions(expression, cell_types, lr_database)`

Compute ligand-receptor interaction scores between cell types.

For each ligand-receptor pair, computes an interaction score between
every source (ligand-expressing) and target (receptor-expressing) cell
type pair. The score is the product of mean ligand expression in the
source type and mean receptor expression in the target type, normalized
by background expression.

Statistical significance is assessed by permutation testing: cell type
labels are shuffled to build a null distribution and a p-value is
computed for each interaction.

Args:
    expression: Expression matrix (n_cells x n_genes) as a numpy array
        or list of lists.
    cell_types: Cell type label for each cell (length n_cells).
    lr_database: Ligand-receptor pair database. Dictionary with key
        ``"pairs"`` mapping to a list of dicts, each with ``"ligand"``
        (gene name) and ``"receptor"`` (gene name). If None, uses the
        built-in database from ``default_lr_database()``.

Returns:
    Dictionary with keys:
        - ``interactions``: List of interaction dicts, each containing
          ``ligand``, ``receptor``, ``source_type``, ``target_type``,
          ``score`` (float), ``p_value`` (float).
        - ``n_significant``: Number of interactions with p < 0.05.
        - ``summary``: Dictionary with total pairs tested, unique
          cell types, unique ligands, unique receptors.

Raises:
    ImportError: If numpy is not available.

#### `spatial_interaction_score()`

**Signature**: `spatial_interaction_score(expression, coordinates, lr_pairs, max_distance)`

Score cell-cell communication considering spatial proximity.

Only counts ligand-receptor interactions between cells that are
within ``max_distance`` of each other. Applies an exponential distance
decay to weight interactions by proximity.

Args:
    expression: Expression matrix (n_cells x n_genes).
    coordinates: (x, y) coordinates for each cell.
    lr_pairs: List of ligand-receptor pair dicts, each with
        ``"ligand_idx"`` (int, gene index) and ``"receptor_idx"``
        (int, gene index).
    max_distance: Maximum Euclidean distance between cells for an
        interaction to be considered.

Returns:
    Dictionary with keys:
        - ``spatial_scores``: List of dicts with ``ligand_idx``,
          ``receptor_idx``, ``score``, ``n_interacting_pairs``.
        - ``distance_decay``: The decay constant used (1 / max_distance).
        - ``significant_pairs``: Number of pairs with score > 0.

Raises:
    ImportError: If numpy is not available.

#### `build_communication_network()`

**Signature**: `build_communication_network(interactions, min_score)`

Build cell-cell communication network from interaction results.

Constructs a directed graph where nodes are cell types and edges
represent significant ligand-receptor interactions. Edge weights are
the sum of interaction scores between each pair of cell types.

Args:
    interactions: List of interaction dicts as returned by
        ``compute_ligand_receptor_interactions``. Each must have
        ``source_type``, ``target_type``, ``score``.
    min_score: Minimum interaction score to include an edge.

Returns:
    Dictionary with keys:
        - ``adjacency_matrix``: 2D list (n_types x n_types) of edge
          weights.
        - ``cell_types``: Sorted list of unique cell type names.
        - ``edge_list``: List of dicts with ``source``, ``target``,
          ``weight``, ``n_interactions``.
        - ``hub_types``: List of cell types with highest total outgoing
          interaction weight (top 3).
        - ``pathway_summary``: Dictionary summarizing interactions by
          ligand-receptor pair.

Raises:
    ValueError: If interactions list is empty.

#### `default_lr_database()`

**Signature**: `default_lr_database()`

Return built-in ligand-receptor pair database.

Provides a curated subset of approximately 200 ligand-receptor pairs
covering major signaling pathways including chemokines, growth factors,
Wnt, Notch, Hedgehog, TGF-beta, interleukins, and adhesion molecules.
Gene names follow HGNC nomenclature.

Returns:
    Dictionary with key ``"pairs"`` containing a list of dicts, each
    with ``"ligand"`` and ``"receptor"`` gene name strings.

#### `communication_pattern_analysis()`

**Signature**: `communication_pattern_analysis(interactions, n_patterns)`

Identify communication patterns using NMF on interaction matrix.

Decomposes the cell-type interaction matrix into a small number of
communication patterns using non-negative matrix factorization. Each
pattern represents a group of correlated ligand-receptor interactions.

Args:
    interactions: Interaction results from
        ``compute_ligand_receptor_interactions``. Must contain
        ``interactions`` list with ``source_type``, ``target_type``,
        ``ligand``, ``receptor``, ``score``.
    n_patterns: Number of communication patterns to discover.

Returns:
    Dictionary with keys:
        - ``patterns``: 2D list (n_patterns x n_lr_pairs) of pattern
          weights.
        - ``pattern_loadings``: 2D list (n_type_pairs x n_patterns) of
          how strongly each cell-type pair participates in each pattern.
        - ``dominant_pathways_per_pattern``: Dictionary mapping pattern
          index to list of top ligand-receptor pairs.

Raises:
    ImportError: If numpy is not available.

#### `_gene_name_to_index()`

**Signature**: `_gene_name_to_index(gene_name, n_genes)`

Convert gene name to index.

If the gene name is numeric, uses it directly as an index. Otherwise,
hashes the name to a deterministic index within the gene space.

Args:
    gene_name: Gene name string or numeric string.
    n_genes: Total number of genes.

Returns:
    Gene index, or None if invalid.

#### `_pairwise_euclidean()`

**Signature**: `_pairwise_euclidean(coords)`

Compute pairwise Euclidean distance matrix without scipy.

Args:
    coords: Coordinate array (n_points x n_dims).

Returns:
    Distance matrix (n_points x n_points).

#### `_simple_nmf()`

**Signature**: `_simple_nmf(matrix, n_components, max_iter, tol)`

Simple multiplicative update NMF implementation.

Factorizes ``matrix ~ W @ H`` where W, H >= 0.

Args:
    matrix: Non-negative input matrix (m x n).
    n_components: Number of components (k).
    max_iter: Maximum iterations.
    tol: Convergence tolerance.

Returns:
    Tuple of (W, H) factor matrices.

## deconvolution.spatial_deconvolution

### Functions

#### `deconvolve_spots()`

**Signature**: `deconvolve_spots(spatial_counts, reference_profiles, method)`

Deconvolve spatial spots into cell type proportions.

Estimates the fraction of each cell type present in every spatial spot
by solving a constrained optimization problem. For NNLS, solves
``min ||counts - R @ w||^2`` subject to ``w >= 0`` for each spot, where
R is the reference profile matrix.

For the regression method, uses ordinary least squares with negative
values clipped to zero and results renormalized.

Args:
    spatial_counts: Expression count matrix (n_spots x n_genes). Can be a
        numpy array, list of lists, or any array-like.
    reference_profiles: Dictionary mapping cell type names to their
        reference expression profiles (each a list of length n_genes).
    method: Deconvolution algorithm. One of ``"nnls"`` (non-negative least
        squares, recommended) or ``"regression"`` (clipped OLS).

Returns:
    Dictionary with keys:
        - ``proportions_matrix``: Array of shape (n_spots, n_types) with
          normalized cell type proportions (rows sum to 1).
        - ``cell_types``: List of cell type names in column order.
        - ``confidence_scores``: Per-spot fitting confidence (1 - normalized
          residual), shape (n_spots,).
        - ``spots_summary``: Dictionary with ``n_spots``, ``n_types``,
          ``mean_confidence``, ``dominant_types`` count per cell type.

Raises:
    ImportError: If numpy is not available.
    ValueError: If method is unrecognized or dimensions mismatch.

#### `build_reference_profiles()`

**Signature**: `build_reference_profiles(sc_expression, cell_types, n_markers)`

Build cell-type reference profiles from single-cell data.

Computes per-cell-type mean expression profiles and selects top marker
genes per type based on fold-change over other types. The returned
profiles are restricted to the union of selected marker genes for
dimensionality reduction.

Args:
    sc_expression: Single-cell expression matrix (n_cells x n_genes).
        Can be dense array or list of lists.
    cell_types: Cell type labels for each cell (length n_cells).
    n_markers: Number of top marker genes to select per cell type.
        The final profile uses the union of all selected markers.

Returns:
    Dictionary with keys:
        - ``profiles``: Dictionary mapping cell type name to expression
          profile (list of floats, length = number of selected genes).
        - ``gene_indices``: Indices of selected marker genes in the
          original gene space.
        - ``n_markers_per_type``: Number of markers selected per type.
        - ``unique_types``: Sorted list of unique cell types.

Raises:
    ImportError: If numpy is not available.
    ValueError: If cell_types length does not match n_cells.

#### `spatial_cell_type_mapping()`

**Signature**: `spatial_cell_type_mapping(proportions, coordinates, cell_types)`

Map cell type proportions to spatial coordinates for visualization.

Assigns a dominant cell type to each spot and computes the mixing entropy
(Shannon entropy of the proportion vector) to quantify how mixed each
spot is between cell types.

Args:
    proportions: Cell type proportions matrix (n_spots x n_types).
        Rows should sum to approximately 1.
    coordinates: List of (x, y) coordinate tuples for each spot.
    cell_types: List of cell type names corresponding to columns.

Returns:
    Dictionary with keys:
        - ``spatial_map``: List of dicts, one per spot, each containing
          ``x``, ``y``, ``dominant_type``, ``proportions`` dict, ``entropy``.
        - ``dominant_type_per_spot``: List of dominant cell type names.
        - ``mixing_entropy``: Array of per-spot Shannon entropy values.

Raises:
    ImportError: If numpy is not available.
    ValueError: If dimensions are inconsistent.

#### `validate_deconvolution()`

**Signature**: `validate_deconvolution(estimated, spatial_markers)`

Validate spatial deconvolution using marker gene expression agreement.

Checks consistency of estimated cell type proportions against known marker
gene expression. For each cell type with known markers, computes the
correlation between estimated proportion and the mean expression of its
markers across spots.

When no marker information is provided, performs internal consistency
checks: proportion normalization, confidence distribution, and dominant
type diversity.

Args:
    estimated: Deconvolution result dictionary as returned by
        ``deconvolve_spots``. Must contain ``proportions_matrix``,
        ``cell_types``, and ``confidence_scores``.
    spatial_markers: Optional dictionary mapping cell type names to lists
        of marker gene names. Used for external validation against
        observed marker expression.

Returns:
    Dictionary with keys:
        - ``normalization_check``: Whether all rows sum to ~1.
        - ``mean_confidence``: Mean confidence across spots.
        - ``confidence_distribution``: Dictionary with ``min``, ``max``,
          ``median``, ``std`` of confidence scores.
        - ``type_diversity``: Number of distinct dominant types observed.
        - ``marker_correlations``: Dictionary mapping cell type to Pearson
          correlation with marker expression (only if spatial_markers
          provided and numpy available).
        - ``overall_validity``: Boolean indicating if basic checks pass.

Raises:
    ImportError: If numpy is not available.

#### `niche_identification()`

**Signature**: `niche_identification(proportions, coordinates, n_niches)`

Identify tissue niches from cell type composition and spatial proximity.

Groups spatial spots into niches (neighborhoods with similar cell type
composition) using k-means clustering on the proportions matrix, weighted
by spatial proximity. Spatial coherence is measured as the fraction of
each spot's spatial neighbors that share the same niche label.

Args:
    proportions: Cell type proportions matrix (n_spots x n_types).
    coordinates: List of (x, y) coordinate tuples for each spot.
    n_niches: Number of niches to identify.

Returns:
    Dictionary with keys:
        - ``niche_labels``: Array of niche assignments (length n_spots),
          integers from 0 to n_niches-1.
        - ``niche_compositions``: Dictionary mapping niche label to mean
          cell type proportion vector.
        - ``spatial_coherence``: Float between 0 and 1 indicating how
          spatially contiguous the niches are.

Raises:
    ImportError: If numpy is not available.
    ValueError: If n_niches exceeds number of spots.

#### `_normalize_rows()`

**Signature**: `_normalize_rows(matrix)`

Normalize each row of a matrix to sum to 1.

Rows with zero sum are left as zeros to avoid division by zero.

Args:
    matrix: 2D numpy array (n_rows x n_cols).

Returns:
    Row-normalized array with same shape.

#### `_regression_deconvolution()`

**Signature**: `_regression_deconvolution(counts, ref_matrix)`

Ordinary least squares deconvolution with non-negativity clipping.

Solves the normal equations ``(R^T R) w = R^T b`` for each spot,
clips negative weights to zero, and returns the weights and residuals.

Args:
    counts: Expression matrix (n_spots x n_genes).
    ref_matrix: Reference matrix (n_genes x n_types).

Returns:
    Tuple of (weights, residuals) arrays.

#### `_kmeans()`

**Signature**: `_kmeans(data, k, max_iter)`

Simple k-means clustering implementation.

Uses k-means++ initialization and Lloyd's algorithm.

Args:
    data: Feature matrix (n_samples x n_features).
    k: Number of clusters.
    max_iter: Maximum number of iterations.

Returns:
    Array of cluster labels (n_samples,).

#### `_pairwise_euclidean()`

**Signature**: `_pairwise_euclidean(coords)`

Compute pairwise Euclidean distance matrix without scipy.

Args:
    coords: Coordinate array (n_points x n_dims).

Returns:
    Distance matrix (n_points x n_points).

## integration.scrna_mapping

### Functions

#### `_to_dense()`

**Signature**: `_to_dense(matrix)`

Convert sparse matrix to dense numpy array.

#### `_intersect_genes()`

**Signature**: `_intersect_genes(spatial_genes, ref_genes)`

Find shared genes between spatial and reference datasets.

Returns:
    Tuple of (spatial_indices, ref_indices, shared_gene_names).

#### `correlation_mapping()`

**Signature**: `correlation_mapping(spatial_expression, scrna_profiles)`

Correlation-based mapping of spatial spots to cell type profiles.

For each spatial spot, computes the correlation between its expression
vector and each reference cell type profile, then assigns the most
correlated type.

Args:
    spatial_expression: Spatial expression matrix (n_spots x n_genes).
    scrna_profiles: Reference cell type profiles (n_types x n_genes).
        Each row is the mean expression for one cell type.
    cell_type_names: Names of cell types. If None, uses indices.
    correlation_method: "pearson" or "spearman".

Returns:
    MappingResult with predicted labels and probabilities.

#### `anchor_based_transfer()`

**Signature**: `anchor_based_transfer(spatial_data, scrna_data, anchors)`

Anchor-based label transfer from scRNA-seq to spatial data.

Inspired by Seurat's anchor-based integration approach:
1. Project both datasets into shared PCA space using shared genes.
2. Find mutual nearest neighbors (MNNs) as anchors between datasets.
3. Transfer labels through anchor-weighted voting.

If pre-computed anchors are provided (as index pairs), uses those directly.
Otherwise, computes MNN anchors.

Args:
    spatial_data: Spatial expression matrix (n_spots x n_spatial_genes).
    scrna_data: scRNA-seq expression matrix (n_cells x n_scrna_genes).
    anchors: Either:
        - Array of (spatial_idx, scrna_idx) pairs, shape (n_anchors, 2).
        - None to auto-compute anchors via MNN.
    scrna_labels: Cell type labels for scRNA-seq cells (length n_cells).
        Required for label transfer.
    cell_type_names: List of possible cell types.
    spatial_genes: Gene names for spatial data.
    scrna_genes: Gene names for scRNA-seq data.
    n_pcs: Number of PCA components for embedding.
    k_anchor: Number of nearest neighbors for MNN anchor finding.
    seed: Random seed.

Returns:
    MappingResult with transferred labels.

#### `map_scrna_to_spatial()`

**Signature**: `map_scrna_to_spatial(scrna_data, spatial_data, method)`

Map scRNA-seq cell type annotations to spatial spots.

Unified interface for label transfer from scRNA-seq reference to
spatial transcriptomics data.

Args:
    scrna_data: scRNA-seq expression matrix (n_cells x n_genes).
    spatial_data: Spatial expression matrix (n_spots x n_genes).
    method: Mapping method ("correlation" or "anchor").
    scrna_labels: Cell type labels for scRNA-seq cells.
    scrna_genes: Gene names for scRNA-seq data.
    spatial_genes: Gene names for spatial data.
    cell_type_names: List of cell type names.
    **kwargs: Additional arguments passed to the specific method.

Returns:
    MappingResult with predicted spatial labels.

#### `impute_spatial_genes()`

**Signature**: `impute_spatial_genes(spatial_data, scrna_data, genes)`

Impute unmeasured genes in spatial data using scRNA-seq reference.

For genes not present in the spatial panel but measured in scRNA-seq:
1. Project both datasets into shared PCA space (using overlapping genes).
2. For each spatial spot, find K nearest scRNA-seq cells in PCA space.
3. Impute target gene expression as weighted average of neighbor values.

Args:
    spatial_data: Spatial expression matrix (n_spots x n_spatial_genes).
    scrna_data: scRNA-seq expression matrix (n_cells x n_scrna_genes).
    genes: List of gene names to impute (must be present in scrna_genes).
    spatial_genes: Gene names for spatial data.
    scrna_genes: Gene names for scRNA-seq data.
    n_neighbors: Number of scRNA-seq neighbors for imputation.
    n_pcs: Number of PCA components.
    seed: Random seed.

Returns:
    ImputationResult with imputed expression for requested genes.

### Classes

#### `MappingResult`

Result of scRNA-seq to spatial mapping.

Attributes:
    predicted_labels: Predicted cell type label per spatial spot (length n_spots).
    prediction_scores: Confidence score per spot (length n_spots).
    label_probabilities: Probability matrix (n_spots x n_types) for each cell type.
    cell_type_names: List of cell type names.
    method: Mapping method used.
    metadata: Additional result metadata.

#### `ImputationResult`

Result of spatial gene imputation.

Attributes:
    imputed_expression: Imputed expression matrix (n_spots x n_imputed_genes).
    gene_names: Names of imputed genes.
    confidence: Per-gene confidence scores.
    method: Imputation method used.

## io.merfish

### Functions

#### `parse_cell_metadata()`

**Signature**: `parse_cell_metadata(metadata_file)`

Parse MERFISH cell metadata CSV file.

Reads cell centroid positions, volumes, FOVs, and any additional columns.
The file should have columns including: cell_id (or CellID), center_x, center_y
(or x, y), and optionally center_z, volume, fov.

Args:
    metadata_file: Path to cell_metadata.csv.

Returns:
    List of CellMetadata records.

Raises:
    FileNotFoundError: If the metadata file does not exist.
    ValueError: If required columns are missing.

#### `load_transcript_spots()`

**Signature**: `load_transcript_spots(spots_file)`

Load individual transcript spot coordinates from a MERFISH detected transcripts file.

The file should have columns: gene, x (or global_x), y (or global_y),
and optionally z, cell_id, quality.

Args:
    spots_file: Path to detected_transcripts.csv.

Returns:
    List of TranscriptSpot records.

Raises:
    FileNotFoundError: If the spots file does not exist.

#### `aggregate_to_cells()`

**Signature**: `aggregate_to_cells(transcript_spots, cell_boundaries)`

Aggregate individual transcript spots to cell-level expression counts.

For each cell, counts the number of transcripts per gene that fall within
the cell's assignment. If transcripts have pre-assigned cell_id, uses that;
otherwise uses nearest-centroid assignment.

Args:
    transcript_spots: List of TranscriptSpot records.
    cell_boundaries: List of CellMetadata records (used for cell IDs and positions).
    unassigned_label: Label for unassigned transcripts (excluded from output).

Returns:
    Tuple of (expression_matrix, cell_ids, gene_names) where expression_matrix
    is a numpy array of shape (n_cells, n_genes) with integer counts.

Raises:
    ImportError: If numpy is not installed.

#### `load_merfish()`

**Signature**: `load_merfish(path)`

Load a complete MERFISH dataset from a directory.

Expects:
    path/
        cell_by_gene.csv: Cell-level expression matrix (cells as rows, genes as columns).
        cell_metadata.csv: Cell positions and metadata.
        detected_transcripts.csv (optional): Per-transcript coordinates.

Args:
    path: Path to the MERFISH output directory.
    load_transcripts: If True, also load individual transcript spots.

Returns:
    MERFISHDataset with expression, coordinates, and metadata.

Raises:
    FileNotFoundError: If required files are missing.
    ImportError: If numpy is not installed.

#### `_find_column()`

**Signature**: `_find_column(headers, candidates)`

Find the first matching column name from a list of candidates.

### Classes

#### `CellMetadata`

Metadata for a single MERFISH cell.

Attributes:
    cell_id: Unique cell identifier.
    x: Cell centroid x-coordinate.
    y: Cell centroid y-coordinate.
    z: Cell centroid z-coordinate (layer/z-plane).
    volume: Cell volume (in cubic microns).
    fov: Field of view index.
    extra: Additional metadata key-value pairs.

#### `TranscriptSpot`

A single detected transcript location.

Attributes:
    gene: Gene name.
    x: Transcript x-coordinate.
    y: Transcript y-coordinate.
    z: Transcript z-coordinate.
    cell_id: Cell assignment (if segmented), or empty.
    quality: Detection quality/confidence score.

#### `MERFISHDataset`

Complete MERFISH spatial dataset.

Attributes:
    expression: Cell-by-gene expression matrix (numpy array, cells x genes).
    coordinates: Cell centroid coordinates array of shape (n_cells, 2).
    cell_ids: List of cell identifier strings.
    gene_names: List of gene names.
    cell_metadata: List of CellMetadata records.
    transcript_spots: Optional list of individual TranscriptSpot records.
    metadata: Additional dataset metadata.

**Methods**:

| Method | Purpose |
|--------|---------|
| `n_cells()` | Number of cells. |
| `n_genes()` | Number of genes. |

## io.visium

### Functions

#### `read_tissue_positions()`

**Signature**: `read_tissue_positions(positions_file)`

Parse a Visium tissue_positions.csv file.

Supports both Space Ranger v1 (tissue_positions_list.csv, no header)
and v2 (tissue_positions.csv, with header) formats.

Args:
    positions_file: Path to tissue_positions.csv or tissue_positions_list.csv.

Returns:
    List of TissuePosition records.

Raises:
    FileNotFoundError: If the positions file does not exist.
    ValueError: If the file format is unrecognized.

#### `read_spatial_image()`

**Signature**: `read_spatial_image(image_path)`

Load an H&E tissue image as a numpy array.

Args:
    image_path: Path to the tissue image (PNG/JPEG).

Returns:
    Numpy array of shape (height, width, channels) with values in [0, 255] uint8.

Raises:
    FileNotFoundError: If image file does not exist.
    ImportError: If Pillow (PIL) is not installed.

#### `filter_tissue_spots()`

**Signature**: `filter_tissue_spots(positions, in_tissue_only)`

Filter tissue position records to keep only spots overlapping tissue.

Args:
    positions: List of TissuePosition records.
    in_tissue_only: If True, keep only spots with in_tissue=True.
        If False, return all positions unchanged.

Returns:
    Filtered list of TissuePosition records.

#### `_read_mex_matrix()`

**Signature**: `_read_mex_matrix(matrix_dir)`

Read a Market Exchange (MEX) format matrix directory.

Expected files: matrix.mtx.gz or matrix.mtx, features.tsv.gz or genes.tsv.gz,
barcodes.tsv.gz or barcodes.tsv.

Returns:
    Tuple of (sparse_matrix, barcodes, gene_names, gene_ids).

#### `_read_h5_matrix()`

**Signature**: `_read_h5_matrix(h5_path)`

Read a 10x HDF5 filtered feature-barcode matrix.

Returns:
    Tuple of (sparse_matrix, barcodes, gene_names, gene_ids).

#### `_read_scale_factors()`

**Signature**: `_read_scale_factors(spatial_dir)`

Read Visium scale factors JSON.

#### `load_visium()`

**Signature**: `load_visium(path)`

Load a 10x Visium Spatial Gene Expression dataset.

Expects the standard Space Ranger output directory structure:
    path/
        filtered_feature_bc_matrix/ (or filtered_feature_bc_matrix.h5)
        spatial/
            tissue_positions_list.csv (or tissue_positions.csv)
            tissue_hires_image.png
            tissue_lowres_image.png
            scalefactors_json.json

Args:
    path: Path to the Space Ranger output directory.
    in_tissue_only: If True, keep only spots overlapping tissue.
    load_image: If True, load the tissue image.
    image_resolution: Which resolution image to load ("hires" or "lowres").

Returns:
    SpatialDataset with expression matrix, coordinates, and optionally image.

Raises:
    FileNotFoundError: If required files are missing.

#### `create_spatial_dataset()`

**Signature**: `create_spatial_dataset(matrix, positions, image, gene_names, gene_ids, scale_factors, platform)`

Create a unified SpatialDataset from components.

Args:
    matrix: Expression matrix (spots x genes) as numpy array or scipy sparse.
    positions: List of TissuePosition records (one per row in matrix).
    image: Optional tissue image array.
    gene_names: List of gene names. Defaults to G0, G1, ...
    gene_ids: List of gene IDs.
    scale_factors: Optional scale factor dictionary.
    platform: Originating platform name.

Returns:
    SpatialDataset instance.

### Classes

#### `TissuePosition`

A single Visium spot position on the tissue.

Attributes:
    barcode: Spot barcode string.
    in_tissue: Whether the spot overlaps tissue (1) or not (0).
    array_row: Row index in the Visium array grid.
    array_col: Column index in the Visium array grid.
    pixel_row: Row pixel coordinate in the full-resolution image.
    pixel_col: Column pixel coordinate in the full-resolution image.

#### `SpatialDataset`

Unified spatial transcriptomics dataset.

Attributes:
    expression: Gene expression matrix (spots x genes) as numpy array or scipy sparse.
    coordinates: Spot coordinates array of shape (n_spots, 2) [row, col in pixels].
    barcodes: List of spot barcodes.
    gene_names: List of gene names/symbols.
    gene_ids: List of gene IDs (e.g., Ensembl).
    tissue_positions: Full TissuePosition records per spot.
    image: Optional tissue image as numpy array (H, W, C).
    scale_factors: Dictionary of Visium scale factors.
    platform: Originating platform (visium, merfish, xenium).
    metadata: Additional metadata dictionary.

**Methods**:

| Method | Purpose |
|--------|---------|
| `n_spots()` | Number of spots/cells. |
| `n_genes()` | Number of genes. |

## io.xenium

### Functions

#### `_open_maybe_gzipped()`

**Signature**: `_open_maybe_gzipped(filepath, mode)`

Open a file that may or may not be gzipped.

#### `read_cell_features()`

**Signature**: `read_cell_features(feature_matrix_path)`

Read a Xenium cell-level feature matrix.

Supports both MEX directory format and HDF5 format.

Args:
    feature_matrix_path: Path to cell_feature_matrix/ directory or .h5 file.

Returns:
    Tuple of (expression_matrix, cell_ids, gene_names, gene_ids).
    Matrix is (n_cells, n_genes) in CSR sparse format.

Raises:
    FileNotFoundError: If the path does not exist.
    ImportError: If required dependencies are missing.

#### `_read_xenium_mex()`

**Signature**: `_read_xenium_mex(matrix_dir)`

Read Xenium MEX format feature matrix.

#### `_read_xenium_h5()`

**Signature**: `_read_xenium_h5(h5_path)`

Read Xenium HDF5 feature matrix.

#### `read_transcripts()`

**Signature**: `read_transcripts(transcripts_path)`

Read per-transcript coordinates from a Xenium transcripts file.

Args:
    transcripts_path: Path to transcripts.csv or transcripts.csv.gz.
    min_quality: Minimum Phred quality score to keep (default 20 = Q20).

Returns:
    List of XeniumTranscript records passing the quality filter.

Raises:
    FileNotFoundError: If the transcripts file does not exist.

#### `load_cell_boundaries()`

**Signature**: `load_cell_boundaries(boundaries_path)`

Load cell segmentation polygon boundaries.

Supports CSV/CSV.GZ format where each row has cell_id, vertex_x, vertex_y,
with multiple rows per cell forming a polygon.

Args:
    boundaries_path: Path to cell_boundaries.csv(.gz) or cell_boundaries.parquet.

Returns:
    List of CellBoundary records with polygon vertices.

Raises:
    FileNotFoundError: If the boundaries file does not exist.

#### `_load_boundaries_parquet()`

**Signature**: `_load_boundaries_parquet(parquet_path)`

Load cell boundaries from a parquet file using pandas.

#### `load_xenium()`

**Signature**: `load_xenium(path)`

Load a complete 10x Xenium dataset from a directory.

Expects the Xenium Ranger output structure:
    path/
        cell_feature_matrix/ (or cell_feature_matrix.h5)
        cells.csv.gz: Cell centroid positions
        transcripts.csv.gz (optional): Per-transcript coordinates
        cell_boundaries.csv.gz (optional): Cell polygon boundaries

Args:
    path: Path to the Xenium output directory.
    load_transcripts: If True, load individual transcript coordinates.
    load_boundaries: If True, load cell segmentation boundaries.
    min_transcript_quality: Minimum quality for transcript filtering.

Returns:
    XeniumDataset with expression matrix, coordinates, and optional extras.

Raises:
    FileNotFoundError: If required files are missing.

#### `_find_col()`

**Signature**: `_find_col(headers, candidates)`

Find first matching column name from candidates.

### Classes

#### `XeniumTranscript`

A single Xenium detected transcript.

Attributes:
    transcript_id: Unique transcript identifier.
    gene: Gene name.
    x: X-coordinate in microns.
    y: Y-coordinate in microns.
    z: Z-coordinate in microns.
    cell_id: Assigned cell ID (0 or empty if unassigned).
    quality_value: Phred-scaled quality score.
    overlaps_nucleus: Whether the transcript overlaps a segmented nucleus.

#### `CellBoundary`

Polygon boundary for a single cell.

Attributes:
    cell_id: Unique cell identifier.
    vertices: List of (x, y) coordinate tuples forming the polygon boundary.

#### `XeniumDataset`

Complete 10x Xenium spatial dataset.

Attributes:
    expression: Cell-by-gene expression matrix (numpy or scipy sparse).
    coordinates: Cell centroid coordinates (n_cells, 2).
    cell_ids: List of cell identifier strings.
    gene_names: List of gene names.
    gene_ids: List of gene IDs.
    transcripts: Optional list of individual transcript records.
    cell_boundaries: Optional list of cell boundary polygons.
    metadata: Additional dataset metadata.

**Methods**:

| Method | Purpose |
|--------|---------|
| `n_cells()` | Number of cells. |
| `n_genes()` | Number of genes. |

## niche.identification

### Functions

#### `identify_niches()`

**Signature**: `identify_niches(cell_type_proportions, coordinates, cell_type_names, n_niches, n_neighbors, spatial_weight, random_state)`

Identify tissue niches from local cell type composition.

1. Build a neighborhood graph from spatial coordinates.
2. Smooth cell type proportions over the spatial neighborhood.
3. Cluster smoothed proportions using K-Means to define niches.
4. Characterize each niche by its mean composition and diversity.

Args:
    cell_type_proportions: 2D array (spots × cell_types), rows sum to 1.
    coordinates: 2D array (spots × 2) of spatial coordinates.
    cell_type_names: Optional names for cell types.
    n_niches: Number of niches to identify.
    n_neighbors: Neighbors for spatial smoothing.
    spatial_weight: Weight for spatial smoothing (0 = no smoothing, 1 = full).
    random_state: Random seed.

Returns:
    NicheResult with niche assignments and characterizations.

#### `_spatial_smooth()`

**Signature**: `_spatial_smooth(values, coords, k, weight)`

Smooth values over spatial neighbors.

Args:
    values: 2D array (spots × features).
    coords: 2D array (spots × 2).
    k: Number of nearest neighbors.
    weight: Smoothing weight (0-1).

Returns:
    Smoothed values.

#### `_kmeans()`

**Signature**: `_kmeans(data, k, rng, max_iter)`

Simple K-Means clustering.

Args:
    data: 2D array (samples × features).
    k: Number of clusters.
    rng: Random state.
    max_iter: Maximum iterations.

Returns:
    1D array of cluster labels.

### Classes

#### `NicheResult`

Result of spatial niche identification.

Attributes:
    niche_labels: Per-spot niche assignment.
    n_niches: Number of identified niches.
    compositions: Per-niche mean cell type composition (niche × cell_type).
    cell_type_names: Names of cell types.
    niche_sizes: Number of spots per niche.
    niche_diversity: Shannon diversity of cell types within each niche.

## visualization.plots

### Functions

#### `_ensure_plotting_deps()`

**Signature**: `_ensure_plotting_deps()`

Check that matplotlib is available.

#### `_save_figure()`

**Signature**: `_save_figure(fig, output_path, dpi)`

Save a matplotlib figure to disk, creating directories as needed.

#### `plot_spatial_scatter()`

**Signature**: `plot_spatial_scatter(coordinates, values, output_path)`

Create a spatial scatter plot colored by continuous or categorical values.

Args:
    coordinates: Spatial coordinates (n x 2).
    values: Values for coloring (length n). Can be numeric or categorical.
    output_path: Path to save the plot image.
    cmap: Matplotlib colormap name.
    title: Plot title.
    point_size: Scatter point size.
    alpha: Point transparency.
    figsize: Figure dimensions (width, height) in inches.
    colorbar_label: Label for the colorbar.
    vmin: Minimum value for color scale.
    vmax: Maximum value for color scale.

Returns:
    Matplotlib Figure object.

#### `plot_tissue_overlay()`

**Signature**: `plot_tissue_overlay(coordinates, values, tissue_image, output_path)`

Overlay expression values on a tissue H&E image.

Args:
    coordinates: Spatial coordinates (n x 2) in pixel space.
    values: Expression values or cluster labels (length n).
    tissue_image: Tissue image array (H, W, C).
    output_path: Path to save the plot.
    cmap: Colormap for expression overlay.
    title: Plot title.
    point_size: Scatter point size.
    alpha: Overlay transparency.
    figsize: Figure dimensions.
    scale_factor: Scale factor to convert coordinates to image pixel space.

Returns:
    Matplotlib Figure object.

#### `plot_gene_expression_map()`

**Signature**: `plot_gene_expression_map(spatial_data, gene, output_path)`

Plot spatial expression map for a single gene.

Args:
    spatial_data: A SpatialDataset (or compatible object with .expression,
        .coordinates, .gene_names attributes).
    gene: Gene name to plot.
    output_path: Path to save the plot.
    cmap: Colormap for expression.
    figsize: Figure dimensions.

Returns:
    Matplotlib Figure object.

Raises:
    ValueError: If gene is not found in the dataset.

#### `plot_cell_type_map()`

**Signature**: `plot_cell_type_map(spatial_data, cell_types, output_path)`

Plot spatial distribution of cell types.

Args:
    spatial_data: A SpatialDataset or compatible object with .coordinates attribute.
    cell_types: Cell type labels (length n_spots), string or integer.
    output_path: Path to save the plot.
    figsize: Figure dimensions.
    point_size: Scatter point size.

Returns:
    Matplotlib Figure object.

#### `plot_neighborhood_graph()`

**Signature**: `plot_neighborhood_graph(coordinates, spatial_graph, output_path)`

Plot the spatial neighborhood graph on tissue coordinates.

Draws edges between connected spots and colors nodes by cluster or expression.

Args:
    coordinates: Spatial coordinates (n x 2).
    spatial_graph: Sparse adjacency matrix (n x n).
    output_path: Path to save the plot.
    node_colors: Values to color nodes by (length n). If None, uniform color.
    title: Plot title.
    figsize: Figure dimensions.
    node_size: Node marker size.
    edge_alpha: Edge line transparency.
    edge_width: Edge line width.

Returns:
    Matplotlib Figure object.

#### `plot_interaction_heatmap()`

**Signature**: `plot_interaction_heatmap(interaction_matrix, output_path)`

Plot a cell type interaction heatmap.

Args:
    interaction_matrix: Interaction scores (n_types x n_types).
    output_path: Path to save the plot.
    cell_type_names: Labels for rows/columns.
    title: Plot title.
    cmap: Colormap.
    figsize: Figure dimensions.
    annot: If True, annotate cells with values.

Returns:
    Matplotlib Figure object.

#### `plot_deconvolution_pie()`

**Signature**: `plot_deconvolution_pie(coordinates, fractions, output_path)`

Plot pie charts per spatial spot showing cell type fractions.

Each spot is represented as a small pie chart showing the estimated
cell type proportions from deconvolution.

Args:
    coordinates: Spatial coordinates (n_spots x 2).
    fractions: Cell type fractions (n_spots x n_types), rows sum to 1.
    output_path: Path to save the plot.
    cell_type_names: Names of cell types.
    title: Plot title.
    figsize: Figure dimensions.
    pie_radius: Radius of each pie chart. If None, auto-calculated.
    min_fraction: Minimum fraction to display (smaller merged into "other").

Returns:
    Matplotlib Figure object.

#### `plot_spatial_autocorrelation()`

**Signature**: `plot_spatial_autocorrelation(coordinates, local_scores, output_path)`

Plot LISA cluster map showing spatial autocorrelation patterns.

Colors spots by their LISA cluster classification:
HH (red), LL (blue), HL (pink), LH (lightblue), NS (gray).

Args:
    coordinates: Spatial coordinates (n x 2).
    local_scores: Local Moran's I values (length n) or similar local statistic.
    output_path: Path to save the plot.
    cluster_labels: LISA cluster labels ("HH", "LL", "HL", "LH", "NS").
        If None, uses local_scores as continuous values.
    title: Plot title.
    figsize: Figure dimensions.
    point_size: Scatter point size.

Returns:
    Matplotlib Figure object.

## Data Structures

| Class | Purpose | Key Attributes |
|-------|---------|----------------|
| Dataset | Container for input data | `.data`, `.metadata` |
| Result | Analysis output | `.values`, `.stats` |
| Config | Settings object | `.params`, `.paths` |
