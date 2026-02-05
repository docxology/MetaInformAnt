"""Spatial autocorrelation statistics for spatial transcriptomics.

Implements classical spatial statistics including Moran's I, Geary's C,
Local Moran's I (LISA), Getis-Ord G*, and spatial variograms for measuring
spatial patterns in gene expression and other continuous variables.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependencies
try:
    import numpy as np
    from numpy.typing import NDArray
except ImportError:
    np = None  # type: ignore[assignment]
    NDArray = None  # type: ignore[assignment,misc]

try:
    from scipy import sparse as sp_sparse
    from scipy.spatial import KDTree
    from scipy.stats import norm
except ImportError:
    sp_sparse = None  # type: ignore[assignment]
    KDTree = None  # type: ignore[assignment,misc]
    norm = None  # type: ignore[assignment]


@dataclass
class MoransIResult:
    """Result of Moran's I spatial autocorrelation test.

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
    """

    I: float
    expected_I: float
    variance_I: float
    z_score: float
    p_value: float
    n: int


@dataclass
class GearyCResult:
    """Result of Geary's C spatial autocorrelation test.

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
    """

    C: float
    expected_C: float
    variance_C: float
    z_score: float
    p_value: float
    n: int


@dataclass
class LocalMoransResult:
    """Result of Local Moran's I (LISA) analysis.

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
    """

    local_I: Any  # np.ndarray (n,)
    expected_I: float
    z_scores: Any  # np.ndarray (n,)
    p_values: Any  # np.ndarray (n,)
    cluster_labels: list[str]
    significance_level: float


@dataclass
class GetisOrdResult:
    """Result of Getis-Ord G* statistic analysis.

    Attributes:
        g_star: G* statistic per observation (length n). Z-score scale.
        p_values: P-values per observation.
        hot_spots: Boolean mask for significant hot spots.
        cold_spots: Boolean mask for significant cold spots.
        significance_level: Alpha level used.
    """

    g_star: Any  # np.ndarray (n,)
    p_values: Any  # np.ndarray (n,)
    hot_spots: Any  # np.ndarray (n,) bool
    cold_spots: Any  # np.ndarray (n,) bool
    significance_level: float


@dataclass
class VariogramResult:
    """Result of spatial variogram/semivariogram analysis.

    Attributes:
        bin_centers: Distance bin centers (length n_bins).
        semivariance: Estimated semivariance at each bin.
        n_pairs: Number of point pairs in each bin.
        nugget: Estimated nugget (semivariance at distance 0).
        sill: Estimated sill (plateau semivariance).
        range_param: Estimated range (distance at which sill is reached).
    """

    bin_centers: Any  # np.ndarray (n_bins,)
    semivariance: Any  # np.ndarray (n_bins,)
    n_pairs: Any  # np.ndarray (n_bins,) int
    nugget: float = 0.0
    sill: float = 0.0
    range_param: float = 0.0


def spatial_weights_matrix(
    coordinates: Any,
    method: Literal["knn", "distance", "binary"] = "knn",
    k: int = 6,
    bandwidth: float | None = None,
    *,
    row_standardize: bool = True,
) -> Any:
    """Build a spatial weights matrix from coordinates.

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
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if sp_sparse is None:
        raise ImportError("scipy is required: uv pip install scipy")
    if KDTree is None:
        raise ImportError("scipy.spatial is required: uv pip install scipy")

    coords = np.asarray(coordinates, dtype=np.float64)
    n = coords.shape[0]
    tree = KDTree(coords)

    if method == "knn":
        k_query = min(k + 1, n)
        distances, indices = tree.query(coords, k=k_query)

        rows: list[int] = []
        cols: list[int] = []
        vals: list[float] = []

        for i in range(n):
            for j_idx in range(1, k_query):
                j = int(indices[i, j_idx])
                rows.append(i)
                cols.append(j)
                vals.append(1.0)

        W = sp_sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
        # Symmetrize
        W = W.maximum(W.T)

    elif method == "distance":
        if bandwidth is None:
            # Use median of nearest-neighbor distance * 1.5
            dists, _ = tree.query(coords, k=2)
            bandwidth = float(np.median(dists[:, 1])) * 1.5

        pairs = tree.query_pairs(r=bandwidth)
        rows = []
        cols = []
        vals = []

        for i, j in pairs:
            d = np.linalg.norm(coords[i] - coords[j])
            if d > 0:
                w = 1.0 / d
                rows.extend([i, j])
                cols.extend([j, i])
                vals.extend([w, w])

        W = sp_sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))

    elif method == "binary":
        if bandwidth is None:
            dists, _ = tree.query(coords, k=2)
            bandwidth = float(np.median(dists[:, 1])) * 1.5

        pairs = tree.query_pairs(r=bandwidth)
        rows = []
        cols = []
        vals = []

        for i, j in pairs:
            rows.extend([i, j])
            cols.extend([j, i])
            vals.extend([1.0, 1.0])

        W = sp_sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))

    else:
        raise ValueError(f"Unknown method: {method}. Use 'knn', 'distance', or 'binary'.")

    if row_standardize:
        row_sums = np.array(W.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1.0
        diag_inv = sp_sparse.diags(1.0 / row_sums)
        W = diag_inv @ W

    logger.info(f"Built spatial weights ({method}): {n} observations, {W.nnz} non-zero entries")
    return W


def morans_i(
    values: Any,
    weights: Any,
) -> MoransIResult:
    """Compute Moran's I spatial autocorrelation statistic.

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
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if norm is None:
        raise ImportError("scipy.stats is required: uv pip install scipy")

    x = np.asarray(values, dtype=np.float64).flatten()
    n = len(x)

    if sp_sparse is not None and sp_sparse.issparse(weights):
        W = weights
    else:
        W = sp_sparse.csr_matrix(np.asarray(weights, dtype=np.float64))

    # Center values
    x_mean = x.mean()
    z = x - x_mean

    # Sum of weights
    S0 = float(W.sum())
    if S0 == 0:
        return MoransIResult(I=0.0, expected_I=0.0, variance_I=0.0, z_score=0.0, p_value=1.0, n=n)

    # Moran's I
    numerator = float(z @ (W @ z))
    denominator = float(z @ z)

    if denominator == 0:
        return MoransIResult(I=0.0, expected_I=0.0, variance_I=0.0, z_score=0.0, p_value=1.0, n=n)

    I = (n / S0) * (numerator / denominator)

    # Expected value under null
    E_I = -1.0 / (n - 1)

    # Variance under normality assumption (randomization)
    S1 = float(((W + W.T).multiply(W + W.T)).sum()) / 2.0
    S2_vec = np.array(W.sum(axis=0)).flatten() + np.array(W.sum(axis=1)).flatten()
    S2 = float((S2_vec**2).sum())

    # Kurtosis term
    b2 = (n * np.sum(z**4)) / (np.sum(z**2) ** 2)

    A = n * ((n**2 - 3 * n + 3) * S1 - n * S2 + 3 * S0**2)
    B = b2 * ((n**2 - n) * S1 - 2 * n * S2 + 6 * S0**2)
    C = (n - 1) * (n - 2) * (n - 3) * S0**2

    if C == 0:
        Var_I = 0.0
    else:
        Var_I = (A - B) / C - E_I**2

    # Z-score and p-value
    if Var_I > 0:
        z_score = (I - E_I) / np.sqrt(Var_I)
        p_value = float(2 * norm.sf(abs(z_score)))
    else:
        z_score = 0.0
        p_value = 1.0

    logger.info(f"Moran's I={I:.4f}, z={z_score:.2f}, p={p_value:.4e}")
    return MoransIResult(
        I=I,
        expected_I=E_I,
        variance_I=Var_I,
        z_score=z_score,
        p_value=p_value,
        n=n,
    )


def gearys_c(
    values: Any,
    weights: Any,
) -> GearyCResult:
    """Compute Geary's C spatial autocorrelation statistic.

    Geary's C is based on squared differences between neighboring observations,
    making it more sensitive to local spatial patterns than Moran's I.

    Formula:
        C = ((n-1) / (2 * S0)) * (sum_ij w_ij (x_i - x_j)^2) / (sum_i (x_i - x_bar)^2)

    Args:
        values: Observed values (length n).
        weights: Spatial weights matrix (n x n).

    Returns:
        GearyCResult with C statistic, z-score, and p-value.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if norm is None:
        raise ImportError("scipy.stats is required: uv pip install scipy")

    x = np.asarray(values, dtype=np.float64).flatten()
    n = len(x)

    if sp_sparse is not None and sp_sparse.issparse(weights):
        W = weights
    else:
        W = sp_sparse.csr_matrix(np.asarray(weights, dtype=np.float64))

    x_mean = x.mean()
    z = x - x_mean
    S0 = float(W.sum())

    if S0 == 0:
        return GearyCResult(C=1.0, expected_C=1.0, variance_C=0.0, z_score=0.0, p_value=1.0, n=n)

    # Numerator: sum of w_ij * (x_i - x_j)^2
    sources, targets = W.nonzero()
    w_vals = np.array(W[sources, targets]).flatten()
    sq_diffs = (x[sources] - x[targets]) ** 2
    numerator = float(np.sum(w_vals * sq_diffs))

    # Denominator: sum of (x_i - x_bar)^2
    denominator = float(np.sum(z**2))

    if denominator == 0:
        return GearyCResult(C=1.0, expected_C=1.0, variance_C=0.0, z_score=0.0, p_value=1.0, n=n)

    C = ((n - 1) / (2.0 * S0)) * (numerator / denominator)

    # Expected C under null
    E_C = 1.0

    # Variance under normality assumption
    S1 = float(((W + W.T).multiply(W + W.T)).sum()) / 2.0
    S2_vec = np.array(W.sum(axis=0)).flatten() + np.array(W.sum(axis=1)).flatten()
    S2 = float((S2_vec**2).sum())

    n_f = float(n)
    if n > 3:
        Var_C = ((2 * S1 + S2) * (n_f - 1) - 4 * S0**2) / (2 * (n_f + 1) * S0**2)
    else:
        Var_C = 0.0

    if Var_C > 0:
        z_score = (C - E_C) / np.sqrt(Var_C)
        p_value = float(2 * norm.sf(abs(z_score)))
    else:
        z_score = 0.0
        p_value = 1.0

    logger.info(f"Geary's C={C:.4f}, z={z_score:.2f}, p={p_value:.4e}")
    return GearyCResult(
        C=C,
        expected_C=E_C,
        variance_C=Var_C,
        z_score=z_score,
        p_value=p_value,
        n=n,
    )


def local_morans_i(
    values: Any,
    weights: Any,
    *,
    significance: float = 0.05,
) -> LocalMoransResult:
    """Compute Local Moran's I (LISA) statistics.

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
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if norm is None:
        raise ImportError("scipy.stats is required: uv pip install scipy")

    x = np.asarray(values, dtype=np.float64).flatten()
    n = len(x)

    if sp_sparse is not None and sp_sparse.issparse(weights):
        W = weights
    else:
        W = sp_sparse.csr_matrix(np.asarray(weights, dtype=np.float64))

    x_mean = x.mean()
    z = x - x_mean

    # Variance of z
    m2 = float(np.sum(z**2) / n)
    if m2 == 0:
        return LocalMoransResult(
            local_I=np.zeros(n),
            expected_I=-1.0 / (n - 1),
            z_scores=np.zeros(n),
            p_values=np.ones(n),
            cluster_labels=["NS"] * n,
            significance_level=significance,
        )

    # Local Moran's I for each observation
    # I_i = (z_i / m2) * sum_j(w_ij * z_j)
    Wz = np.array(W @ z).flatten()
    local_I = (z / m2) * Wz

    # Expected value
    E_Ii = -1.0 / (n - 1)

    # Variance of local I (under randomization)
    # Simplified: use conditional permutation variance approximation
    w_i_sq = np.array(W.multiply(W).sum(axis=1)).flatten()
    w_i_sum = np.array(W.sum(axis=1)).flatten()

    b2 = (np.sum(z**4) / n) / (m2**2)

    Var_Ii = np.zeros(n, dtype=np.float64)
    for i in range(n):
        wii_sq = w_i_sq[i]
        wi_sum = w_i_sum[i]
        Var_Ii[i] = wii_sq * (n - b2) / (n - 1) + wi_sum**2 * (2 * b2 - n) / ((n - 1) * (n - 2)) - E_Ii**2

    Var_Ii = np.maximum(Var_Ii, 1e-20)

    # Z-scores and p-values
    z_scores = (local_I - E_Ii) / np.sqrt(Var_Ii)
    p_values = 2 * norm.sf(np.abs(z_scores))

    # Spatial lag (mean of neighbors)
    spatial_lag = Wz  # already W @ z

    # Classify clusters
    cluster_labels: list[str] = []
    for i in range(n):
        if p_values[i] > significance:
            cluster_labels.append("NS")
        elif z[i] > 0 and spatial_lag[i] > 0:
            cluster_labels.append("HH")
        elif z[i] < 0 and spatial_lag[i] < 0:
            cluster_labels.append("LL")
        elif z[i] > 0 and spatial_lag[i] < 0:
            cluster_labels.append("HL")
        elif z[i] < 0 and spatial_lag[i] > 0:
            cluster_labels.append("LH")
        else:
            cluster_labels.append("NS")

    n_sig = sum(1 for l in cluster_labels if l != "NS")
    logger.info(f"Local Moran's I: {n_sig}/{n} significant at alpha={significance}")

    return LocalMoransResult(
        local_I=local_I,
        expected_I=E_Ii,
        z_scores=z_scores,
        p_values=p_values,
        cluster_labels=cluster_labels,
        significance_level=significance,
    )


def getis_ord_g(
    values: Any,
    weights: Any,
    *,
    significance: float = 0.05,
) -> GetisOrdResult:
    """Compute Getis-Ord G* statistic for hot/cold spot detection.

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
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if norm is None:
        raise ImportError("scipy.stats is required: uv pip install scipy")

    x = np.asarray(values, dtype=np.float64).flatten()
    n = len(x)

    if sp_sparse is not None and sp_sparse.issparse(weights):
        W = weights
    else:
        W = sp_sparse.csr_matrix(np.asarray(weights, dtype=np.float64))

    x_bar = x.mean()
    S = np.sqrt(np.sum((x - x_bar) ** 2) / n)

    if S == 0:
        return GetisOrdResult(
            g_star=np.zeros(n),
            p_values=np.ones(n),
            hot_spots=np.zeros(n, dtype=bool),
            cold_spots=np.zeros(n, dtype=bool),
            significance_level=significance,
        )

    # Include self-weight (G* vs G)
    # Add identity to weights for G*
    W_star = W + sp_sparse.eye(n)

    g_star = np.zeros(n, dtype=np.float64)

    W_star_sum = np.array(W_star.sum(axis=1)).flatten()
    W_star_sq_sum = np.array(W_star.multiply(W_star).sum(axis=1)).flatten()
    Wx = np.array(W_star @ x).flatten()

    for i in range(n):
        wi_sum = W_star_sum[i]
        wi_sq_sum = W_star_sq_sum[i]

        numerator = Wx[i] - x_bar * wi_sum
        denominator_inner = (n * wi_sq_sum - wi_sum**2) / (n - 1)

        if denominator_inner > 0:
            denominator = S * np.sqrt(denominator_inner)
            g_star[i] = numerator / denominator
        else:
            g_star[i] = 0.0

    # P-values (two-sided normal)
    p_values = 2 * norm.sf(np.abs(g_star))

    # Identify hot/cold spots
    hot_spots = (g_star > 0) & (p_values < significance)
    cold_spots = (g_star < 0) & (p_values < significance)

    n_hot = int(hot_spots.sum())
    n_cold = int(cold_spots.sum())
    logger.info(f"Getis-Ord G*: {n_hot} hot spots, {n_cold} cold spots at alpha={significance}")

    return GetisOrdResult(
        g_star=g_star,
        p_values=p_values,
        hot_spots=hot_spots,
        cold_spots=cold_spots,
        significance_level=significance,
    )


def spatial_variogram(
    values: Any,
    coordinates: Any,
    n_bins: int = 20,
    *,
    max_distance: float | None = None,
) -> VariogramResult:
    """Compute the empirical spatial variogram (semivariogram).

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
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")

    x = np.asarray(values, dtype=np.float64).flatten()
    coords = np.asarray(coordinates, dtype=np.float64)
    n = len(x)

    if n < 2:
        return VariogramResult(
            bin_centers=np.array([]),
            semivariance=np.array([]),
            n_pairs=np.array([], dtype=np.int64),
        )

    # Compute pairwise distances and squared differences
    # For large datasets, use KDTree for efficiency
    if KDTree is not None and n > 500:
        tree = KDTree(coords)

        if max_distance is None:
            # Estimate max distance as half the diagonal
            extent = coords.max(axis=0) - coords.min(axis=0)
            max_distance = float(np.linalg.norm(extent)) / 2.0

        bin_edges = np.linspace(0, max_distance, n_bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
        semivariance = np.zeros(n_bins, dtype=np.float64)
        pair_counts = np.zeros(n_bins, dtype=np.int64)

        # Query pairs within max_distance
        pairs = tree.query_pairs(r=max_distance)
        for i, j in pairs:
            d = np.linalg.norm(coords[i] - coords[j])
            sq_diff = (x[i] - x[j]) ** 2
            bin_idx = int(np.searchsorted(bin_edges[1:], d))
            if 0 <= bin_idx < n_bins:
                semivariance[bin_idx] += sq_diff
                pair_counts[bin_idx] += 1

    else:
        # Brute force for small datasets
        all_dists: list[float] = []
        all_sq_diffs: list[float] = []

        for i in range(n):
            for j in range(i + 1, n):
                d = float(np.linalg.norm(coords[i] - coords[j]))
                sq_diff = (x[i] - x[j]) ** 2
                all_dists.append(d)
                all_sq_diffs.append(sq_diff)

        dists_arr = np.array(all_dists)
        sq_diffs_arr = np.array(all_sq_diffs)

        if max_distance is None:
            max_distance = float(dists_arr.max()) / 2.0

        bin_edges = np.linspace(0, max_distance, n_bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
        semivariance = np.zeros(n_bins, dtype=np.float64)
        pair_counts = np.zeros(n_bins, dtype=np.int64)

        for k in range(len(dists_arr)):
            d = dists_arr[k]
            if d > max_distance:
                continue
            bin_idx = int(np.searchsorted(bin_edges[1:], d))
            if 0 <= bin_idx < n_bins:
                semivariance[bin_idx] += sq_diffs_arr[k]
                pair_counts[bin_idx] += 1

    # Average semivariance per bin: gamma(h) = (1/2N(h)) * sum (x_i - x_j)^2
    valid = pair_counts > 0
    semivariance[valid] = semivariance[valid] / (2.0 * pair_counts[valid])

    # Estimate variogram parameters
    nugget = float(semivariance[valid][0]) if valid.any() else 0.0
    sill = float(np.max(semivariance[valid])) if valid.any() else 0.0

    # Range: first bin where semivariance reaches ~95% of sill
    range_param = 0.0
    if sill > 0 and valid.any():
        threshold = 0.95 * sill
        for bi in range(n_bins):
            if valid[bi] and semivariance[bi] >= threshold:
                range_param = float(bin_centers[bi])
                break

    logger.info(f"Variogram: nugget={nugget:.4f}, sill={sill:.4f}, " f"range={range_param:.2f}, {n_bins} bins")

    return VariogramResult(
        bin_centers=bin_centers,
        semivariance=semivariance,
        n_pairs=pair_counts,
        nugget=nugget,
        sill=sill,
        range_param=range_param,
    )
