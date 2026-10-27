"""Spatial neighborhood analysis for spatial transcriptomics.

Implements cell type co-localization analysis, interaction scoring,
spatial ligand-receptor analysis, niche detection, and Ripley's K function
for spatial point pattern statistics.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal, cast

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependencies. Type checkers model them as installed; at runtime a missing
# dependency leaves the name bound to None and the functions that need it raise ImportError.
try:
    import numpy as np
except ImportError:
    if not TYPE_CHECKING:
        np = None

try:
    from scipy import sparse as sp_sparse
    from scipy.spatial import KDTree
except ImportError:
    if not TYPE_CHECKING:
        sp_sparse = None
        KDTree = None

try:
    from sklearn.cluster import KMeans
except ImportError:
    KMeans = None


@dataclass
class NeighborhoodEnrichmentResult:
    """Result of neighborhood enrichment analysis.

    Attributes:
        enrichment_matrix: Enrichment Z-scores (n_types x n_types).
            Positive = co-localized more than expected, negative = avoided.
        count_matrix: Observed interaction counts (n_types x n_types).
        expected_matrix: Expected interaction counts under random spatial arrangement.
        p_values: P-values from permutation testing (n_types x n_types).
        cell_type_names: List of cell type names.
    """

    enrichment_matrix: Any  # np.ndarray (n_types, n_types)
    count_matrix: Any  # np.ndarray (n_types, n_types)
    expected_matrix: Any  # np.ndarray (n_types, n_types)
    p_values: Any  # np.ndarray (n_types, n_types)
    cell_type_names: list[str]


@dataclass
class InteractionResult:
    """Result of pairwise cell type interaction analysis.

    Attributes:
        interaction_matrix: Interaction scores (n_types x n_types).
        cell_type_names: List of cell type names.
        method: Scoring method used.
    """

    interaction_matrix: Any  # np.ndarray (n_types, n_types)
    cell_type_names: list[str]
    method: str


@dataclass
class NicheResult:
    """Result of cellular niche detection.

    Attributes:
        niche_labels: Niche assignment per cell/spot (length n).
        niche_compositions: Cell type composition per niche (n_niches x n_types).
        n_niches: Number of niches found.
        cell_type_names: List of cell type names.
    """

    niche_labels: Any  # np.ndarray (n,)
    niche_compositions: Any  # np.ndarray (n_niches, n_types)
    n_niches: int
    cell_type_names: list[str]


@dataclass
class RipleyKResult:
    """Result of Ripley's K function analysis.

    Attributes:
        radii: Array of evaluation radii.
        k_values: K(r) values at each radius.
        l_values: L(r) = sqrt(K(r)/pi) - r (Besag's L-function, centered).
        csr_envelope_lower: Lower bound of CSR envelope (from simulations).
        csr_envelope_upper: Upper bound of CSR envelope.
        n_points: Number of points analyzed.
        area: Study area.
    """

    radii: Any  # np.ndarray
    k_values: Any  # np.ndarray
    l_values: Any  # np.ndarray
    csr_envelope_lower: Any | None = None  # np.ndarray or None
    csr_envelope_upper: Any | None = None  # np.ndarray or None
    n_points: int = 0
    area: float = 0.0


def neighborhood_enrichment(
    cell_types: Any,
    coordinates: Any,
    radius: float | None = None,
    *,
    n_neighbors: int = 6,
    n_permutations: int = 1000,
    seed: int = 42,
) -> NeighborhoodEnrichmentResult:
    """Compute cell type neighborhood enrichment (co-localization analysis).

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
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if KDTree is None:
        raise ImportError("scipy.spatial is required: uv pip install scipy")

    rng = np.random.RandomState(seed)
    types = np.asarray(cell_types)
    coords = np.asarray(coordinates, dtype=np.float64)
    n = len(types)

    unique_types = sorted(set(types.tolist()))
    type_to_idx = {t: i for i, t in enumerate(unique_types)}
    n_types = len(unique_types)
    type_indices = np.array([type_to_idx[t] for t in types], dtype=np.int32)

    # Build neighbor lists
    tree = KDTree(coords)
    if radius is not None:
        neighbor_pairs = tree.query_pairs(r=radius)
        neighbors_dict: dict[int, list[int]] = {i: [] for i in range(n)}
        for i, j in neighbor_pairs:
            neighbors_dict[i].append(j)
            neighbors_dict[j].append(i)
    else:
        k = min(n_neighbors + 1, n)
        _, indices = tree.query(coords, k=k)
        neighbors_dict = {}
        for i in range(n):
            neighbors_dict[i] = [int(j) for j in indices[i, 1:] if j != i]

    # Count observed interactions
    observed = np.zeros((n_types, n_types), dtype=np.float64)
    for i in range(n):
        ti = type_indices[i]
        for j in neighbors_dict[i]:
            tj = type_indices[j]
            observed[ti, tj] += 1.0

    # Permutation test for null distribution
    perm_counts = np.zeros((n_permutations, n_types, n_types), dtype=np.float64)

    for p in range(n_permutations):
        perm_labels = rng.permutation(type_indices)
        for i in range(n):
            ti = perm_labels[i]
            for j in neighbors_dict[i]:
                tj = perm_labels[j]
                perm_counts[p, ti, tj] += 1.0

    # Compute statistics
    expected = perm_counts.mean(axis=0)
    std_null = perm_counts.std(axis=0)
    std_null[std_null == 0] = 1.0  # avoid division by zero

    enrichment = (observed - expected) / std_null

    # P-values: fraction of permutations with count >= observed
    p_values = np.zeros((n_types, n_types), dtype=np.float64)
    for ti in range(n_types):
        for tj in range(n_types):
            p_values[ti, tj] = np.mean(perm_counts[:, ti, tj] >= observed[ti, tj])

    type_names = [str(t) for t in unique_types]
    logger.info(
        f"Neighborhood enrichment: {n_types} types, "
        f"{n_permutations} permutations, max enrichment={enrichment.max():.2f}"
    )

    return NeighborhoodEnrichmentResult(
        enrichment_matrix=enrichment,
        count_matrix=observed,
        expected_matrix=expected,
        p_values=p_values,
        cell_type_names=type_names,
    )


def compute_interaction_matrix(
    cell_types: Any,
    spatial_graph: Any,
    *,
    normalize_by_type_frequency: bool = True,
) -> InteractionResult:
    """Compute pairwise cell type interaction scores from a spatial graph.

    For each pair of cell types, computes the interaction score as the number
    of edges between them, optionally normalized by the product of their
    frequencies (to account for abundance).

    Args:
        cell_types: Array of cell type labels (length n).
        spatial_graph: Sparse adjacency matrix (n x n).
        normalize_by_type_frequency: If True, normalize by expected frequency.

    Returns:
        InteractionResult with interaction matrix.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if sp_sparse is None:
        raise ImportError("scipy is required: uv pip install scipy")

    types = np.asarray(cell_types)
    adj = sp_sparse.csr_matrix(spatial_graph)
    n = len(types)

    unique_types = sorted(set(types.tolist()))
    type_to_idx = {t: i for i, t in enumerate(unique_types)}
    n_types = len(unique_types)
    type_indices = np.array([type_to_idx[t] for t in types], dtype=np.int32)

    # Count interactions
    interaction = np.zeros((n_types, n_types), dtype=np.float64)
    sources, targets = adj.nonzero()
    weights = np.array(adj[sources, targets]).flatten()

    for s, t, w in zip(sources, targets, weights):
        ti = type_indices[s]
        tj = type_indices[t]
        interaction[ti, tj] += w

    # Symmetrize (for undirected graphs, should already be symmetric)
    interaction = (interaction + interaction.T) / 2.0

    method = "raw_count"
    if normalize_by_type_frequency:
        # Compute type frequencies
        type_counts = np.array([np.sum(type_indices == i) for i in range(n_types)], dtype=np.float64)
        freq = type_counts / n
        # Expected interaction under random mixing
        total_edges = adj.sum() / 2.0
        expected = np.outer(freq, freq) * total_edges * 2.0
        expected[expected == 0] = 1.0
        interaction = interaction / expected
        method = "frequency_normalized"

    type_names = [str(t) for t in unique_types]
    logger.info(f"Interaction matrix: {n_types} types, method={method}")

    return InteractionResult(
        interaction_matrix=interaction,
        cell_type_names=type_names,
        method=method,
    )


def ligand_receptor_spatial(
    expression: Any,
    lr_pairs: list[tuple[str, str]],
    coordinates: Any,
    *,
    gene_names: list[str] | None = None,
    radius: float | None = None,
    n_neighbors: int = 6,
) -> dict[str, Any]:
    """Spatial ligand-receptor interaction analysis.

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
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if KDTree is None:
        raise ImportError("scipy.spatial is required: uv pip install scipy")

    if sp_sparse is not None and sp_sparse.issparse(expression):
        # scipy-stubs narrow issparse() to a method-less SparseABC; every sparse format provides toarray()
        expr = cast("sp_sparse.csr_matrix", expression).toarray()
    else:
        expr = np.asarray(expression, dtype=np.float64)

    coords = np.asarray(coordinates, dtype=np.float64)
    n = expr.shape[0]

    if gene_names is None:
        gene_names = [str(i) for i in range(expr.shape[1])]

    gene_to_col = {g: i for i, g in enumerate(gene_names)}

    # Build neighbor lists
    tree = KDTree(coords)
    if radius is not None:
        pairs_set = tree.query_pairs(r=radius)
        neighbors_dict: dict[int, list[int]] = {i: [] for i in range(n)}
        for i, j in pairs_set:
            neighbors_dict[i].append(j)
            neighbors_dict[j].append(i)
    else:
        k = min(n_neighbors + 1, n)
        _, indices = tree.query(coords, k=k)
        neighbors_dict = {}
        for i in range(n):
            neighbors_dict[i] = [int(j) for j in indices[i, 1:] if j != i]

    scores: dict[tuple[str, str], float] = {}
    per_spot_scores: dict[tuple[str, str], Any] = {}
    n_tested = 0

    for ligand, receptor in lr_pairs:
        if ligand not in gene_to_col or receptor not in gene_to_col:
            logger.debug(f"Skipping L-R pair ({ligand}, {receptor}): gene(s) not in expression data")
            continue

        l_col = gene_to_col[ligand]
        r_col = gene_to_col[receptor]

        ligand_expr = expr[:, l_col]
        receptor_expr = expr[:, r_col]

        spot_scores = np.zeros(n, dtype=np.float64)
        for i in range(n):
            nbrs = neighbors_dict[i]
            if len(nbrs) == 0:
                continue
            neighbor_receptor_mean = receptor_expr[nbrs].mean()
            spot_scores[i] = ligand_expr[i] * neighbor_receptor_mean

        pair_key = (ligand, receptor)
        scores[pair_key] = float(spot_scores.mean())
        per_spot_scores[pair_key] = spot_scores
        n_tested += 1

    logger.info(f"Ligand-receptor analysis: {n_tested}/{len(lr_pairs)} pairs tested")
    return {
        "scores": scores,
        "per_spot_scores": per_spot_scores,
        "n_pairs_tested": n_tested,
    }


def niche_detection(
    cell_types: Any,
    coordinates: Any,
    n_niches: int = 5,
    *,
    n_neighbors: int = 10,
    seed: int = 42,
) -> NicheResult:
    """Identify cellular niches based on local cell type composition.

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
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if KDTree is None:
        raise ImportError("scipy.spatial is required: uv pip install scipy")
    if KMeans is None:
        raise ImportError("scikit-learn is required for niche detection: uv pip install scikit-learn")

    types = np.asarray(cell_types)
    coords = np.asarray(coordinates, dtype=np.float64)
    n = len(types)

    unique_types = sorted(set(types.tolist()))
    type_to_idx = {t: i for i, t in enumerate(unique_types)}
    n_types = len(unique_types)

    # Build KNN
    k = min(n_neighbors + 1, n)
    tree = KDTree(coords)
    _, indices = tree.query(coords, k=k)

    # Compute neighborhood composition for each cell
    compositions = np.zeros((n, n_types), dtype=np.float64)
    for i in range(n):
        for j_idx in range(k):
            j = indices[i, j_idx]
            tj = type_to_idx[types[j]]
            compositions[i, tj] += 1.0
        # Normalize to fractions
        total = compositions[i, :].sum()
        if total > 0:
            compositions[i, :] /= total

    # Cluster compositions to find niches
    n_niches_actual = min(n_niches, n)
    km = KMeans(n_clusters=n_niches_actual, random_state=seed, n_init=10)
    niche_labels = km.fit_predict(compositions)

    # Compute per-niche composition
    niche_compositions = np.zeros((n_niches_actual, n_types), dtype=np.float64)
    for ni in range(n_niches_actual):
        mask = niche_labels == ni
        if mask.any():
            niche_compositions[ni, :] = compositions[mask, :].mean(axis=0)

    type_names = [str(t) for t in unique_types]
    logger.info(f"Niche detection: {n_niches_actual} niches from {n} cells, {n_types} types")

    return NicheResult(
        niche_labels=niche_labels.astype(np.int32),
        niche_compositions=niche_compositions,
        n_niches=n_niches_actual,
        cell_type_names=type_names,
    )


def _isotropic_edge_weights(x: Any, y: Any, radius: Any, window: tuple[float, float, float, float]) -> Any:
    """Ripley's isotropic edge-correction weights for circles about points inside a rectangular window.

    The weight of a circle of radius ``radius`` centred at ``(x, y)`` is 1 / (share of its circumference that
    lies inside the window). The share is exact for a rectangle: the circle leaves through a side at
    distance ``d`` when ``radius > d``, over an arc of half-width ``arccos(d / radius)`` centred on that
    side's outward normal, and arcs of opposite sides never meet, so the escaped length is the sum of the
    four arcs minus the overlaps of the four pairs of adjacent sides (corner regions).

    Args:
        x, y: Coordinates of the circle centres (arrays of equal length, inside the window).
        radius: Circle radii (array of the same length).
        window: ``(x_min, x_max, y_min, y_max)``.

    Returns:
        Array of weights, 1 where the circle lies entirely inside the window.
    """
    x_min, x_max, y_min, y_max = window
    x, y, radius = np.asarray(x, float), np.asarray(y, float), np.asarray(radius, float)

    # Distances to the right, top, left and bottom sides: outward normals at 0, 90, 180 and 270 degrees
    side_distances = (x_max - x, y_max - y, x - x_min, y - y_min)
    half_widths = []
    with np.errstate(divide="ignore", invalid="ignore"):
        for distance in side_distances:
            ratio = np.clip(distance / radius, -1.0, 1.0)
            half_widths.append(np.where(radius > distance, np.arccos(ratio), 0.0))

    escaped = 2.0 * sum(half_widths)
    for k in range(4):
        escaped = escaped - np.maximum(0.0, half_widths[k] + half_widths[(k + 1) % 4] - np.pi / 2.0)
    inside = 1.0 - escaped / (2.0 * np.pi)
    with np.errstate(divide="ignore"):
        return np.where(inside > 0.0, 1.0 / inside, np.inf)


def _pairs_within(point_set: Any, max_radius: float, use_tree: bool) -> tuple[Any, Any]:
    """Index arrays ``(i, j)`` of every unordered pair of points with ``i < j`` at distance <= ``max_radius``."""
    if use_tree:
        pairs = KDTree(point_set).query_pairs(r=max_radius, output_type="ndarray")
        return pairs[:, 0], pairs[:, 1]

    first: list[Any] = []
    second: list[Any] = []
    for a in range(len(point_set) - 1):
        distances = np.linalg.norm(point_set[a + 1 :] - point_set[a], axis=1)
        hits = np.flatnonzero(distances <= max_radius) + a + 1
        first.append(np.full(len(hits), a))
        second.append(hits)
    if not first:
        return np.empty(0, dtype=int), np.empty(0, dtype=int)
    return np.concatenate(first), np.concatenate(second)


def _ripley_pair_sums(
    point_set: Any,
    radii: Any,
    window: tuple[float, float, float, float],
    edge_correction: str,
    use_tree: bool,
) -> Any:
    """For each radius r, the sum over ordered pairs (i, j), i != j, with d_ij <= r of the pair's weight."""
    i, j = _pairs_within(point_set, float(radii.max()), use_tree)
    if len(i) == 0:
        return np.zeros(len(radii), dtype=np.float64)

    distance = np.linalg.norm(point_set[i] - point_set[j], axis=1)
    if edge_correction == "isotropic":
        # the ordered pair (i, j) is weighted by the circle about i, and (j, i) by the circle about j
        weights = _isotropic_edge_weights(point_set[i, 0], point_set[i, 1], distance, window)
        weights = weights + _isotropic_edge_weights(point_set[j, 0], point_set[j, 1], distance, window)
    else:
        weights = np.full(len(distance), 2.0)

    order = np.argsort(distance)
    cumulative = np.cumsum(weights[order])
    n_within = np.searchsorted(distance[order], radii, side="right")
    return np.where(n_within > 0, cumulative[np.maximum(n_within - 1, 0)], 0.0)


def ripley_k(
    points: Any,
    radii: Any,
    area: float,
    *,
    n_simulations: int = 99,
    seed: int = 42,
    edge_correction: Literal["isotropic", "none"] = "isotropic",
) -> RipleyKResult:
    """Compute Ripley's K function for spatial point pattern analysis.

    Ripley's K(r) counts the expected number of points within distance r
    of a typical point, normalized by intensity. Under Complete Spatial
    Randomness (CSR), K(r) = pi * r^2.

    Also computes Besag's L-function: L(r) = sqrt(K(r)/pi) - r,
    which is centered at 0 under CSR. Positive L(r) indicates clustering,
    negative indicates regularity/inhibition.

    By default Ripley's isotropic edge correction is applied, for every number of points: each ordered pair
    (i, j) at distance d_ij <= r counts with weight 1 / (share of the circle of radius d_ij about point i
    that lies inside the study window), the exact share for a rectangle, so that
    ``K(r) = area / (n (n - 1)) * sum of the weights``. The window is the bounding box of ``points``.
    ``edge_correction="none"`` counts every pair with weight 1 instead, which is biased low near the
    window boundary. Keep the radii well below half the shorter side of the window (a quarter is a
    common rule): the weights grow quickly beyond that. Memory grows with the number of point pairs
    within the largest radius.

    Args:
        points: Point coordinates (n x 2).
        radii: Array of radii at which to evaluate K.
        area: Total study area (e.g., bounding box area). It sets the intensity only; the edge-correction
            window is always the bounding box of ``points``.
        n_simulations: Number of CSR simulations (uniform in the bounding box, with the same edge
            correction) for the confidence envelope.
        seed: Random seed.
        edge_correction: ``"isotropic"`` (default) or ``"none"``.

    Returns:
        RipleyKResult with K values, L values, and CSR envelope.

    Raises:
        ValueError: If ``edge_correction`` is unknown, ``radii`` is empty or negative, or the isotropic
            correction is requested for points with no extent in x or in y (no window to correct for).
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if edge_correction not in ("isotropic", "none"):
        raise ValueError(f"Unknown edge_correction {edge_correction!r}: use 'isotropic' or 'none'")

    rng = np.random.RandomState(seed)
    pts = np.asarray(points, dtype=np.float64)
    r_arr = np.asarray(radii, dtype=np.float64)
    n = pts.shape[0]

    if r_arr.size == 0 or float(r_arr.min()) < 0.0:
        raise ValueError("radii must be a non-empty array of non-negative values")

    if n < 2:
        return RipleyKResult(
            radii=r_arr,
            k_values=np.zeros_like(r_arr),
            l_values=np.zeros_like(r_arr),
            n_points=n,
            area=area,
        )

    # Bounding box: the window for edge correction and for the CSR simulations
    x_min, y_min = pts.min(axis=0)
    x_max, y_max = pts.max(axis=0)
    if edge_correction == "isotropic" and not (x_max > x_min and y_max > y_min):
        raise ValueError(
            "the points have no extent in x or in y, so there is no study window to correct for; "
            "use edge_correction='none' for such a pattern"
        )
    window = (float(x_min), float(x_max), float(y_min), float(y_max))
    use_tree = KDTree is not None

    def _compute_k(point_set: Any, n_pts: int) -> Any:
        """K(r) = area / (n (n - 1)) * sum of the pair weights."""
        return (area / (n_pts * (n_pts - 1))) * _ripley_pair_sums(point_set, r_arr, window, edge_correction, use_tree)

    # Observed K
    k_observed = _compute_k(pts, n)

    # L-function: L(r) = sqrt(K(r)/pi) - r
    l_observed = np.sqrt(np.maximum(k_observed, 0) / np.pi) - r_arr

    # CSR envelope via simulation
    k_simulations = np.zeros((n_simulations, len(r_arr)), dtype=np.float64)
    for sim in range(n_simulations):
        # Generate CSR points in bounding box
        sim_x = rng.uniform(x_min, x_max, n)
        sim_y = rng.uniform(y_min, y_max, n)
        sim_pts = np.column_stack([sim_x, sim_y])
        k_simulations[sim, :] = _compute_k(sim_pts, n)

    l_simulations = np.sqrt(np.maximum(k_simulations, 0) / np.pi) - r_arr

    # 95% envelope
    envelope_lower = np.percentile(l_simulations, 2.5, axis=0)
    envelope_upper = np.percentile(l_simulations, 97.5, axis=0)

    logger.info(f"Ripley's K: {n} points, {len(r_arr)} radii, " f"max L={l_observed.max():.3f}")

    return RipleyKResult(
        radii=r_arr,
        k_values=k_observed,
        l_values=l_observed,
        csr_envelope_lower=envelope_lower,
        csr_envelope_upper=envelope_upper,
        n_points=n,
        area=area,
    )
