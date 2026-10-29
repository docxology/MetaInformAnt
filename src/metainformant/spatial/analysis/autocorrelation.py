"""Spatial autocorrelation statistics for spatial transcriptomics.

Implements classical spatial statistics including Moran's I, Geary's C,
Local Moran's I (LISA), Getis-Ord G*, and spatial variograms for measuring
spatial patterns in gene expression and other continuous variables.

A statistic that is undefined for the input (constant values, a weights matrix without edges, too few
observations for a variance) is NaN, never a stand-in such as I = 0, C = 1 or p = 1: those would read as
a measured "no autocorrelation". The result then says why (``status`` on the global results; NaN entries on
the per-observation ones) and its p-value is NaN.
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
    from scipy.stats import norm
except ImportError:
    if not TYPE_CHECKING:
        sp_sparse = None
        KDTree = None
        norm = None


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
        status: ``"ok"``, or why the statistic (``I``, NaN) or its variance, z-score and p-value (NaN) are undefined.
    """

    I: float  # noqa: E741 - public API uses Moran's I notation.
    expected_I: float
    variance_I: float
    z_score: float
    p_value: float
    n: int
    status: str = "ok"


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
        status: ``"ok"``, or why the statistic (``C``, NaN) or its variance, z-score and p-value (NaN) are undefined.
    """

    C: float
    expected_C: float
    variance_C: float
    z_score: float
    p_value: float
    n: int
    status: str = "ok"


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
            "NS" = Not significant (also for observations whose z-score and p-value are NaN).
        significance_level: Alpha level used for classification.
        status: ``"ok"``, or why the statistics of every observation are undefined (NaN).
    """

    local_I: Any  # np.ndarray (n,)
    expected_I: float
    z_scores: Any  # np.ndarray (n,)
    p_values: Any  # np.ndarray (n,)
    cluster_labels: list[str]
    significance_level: float
    status: str = "ok"


@dataclass
class GetisOrdResult:
    """Result of Getis-Ord G* statistic analysis.

    Attributes:
        g_star: G* statistic per observation (length n). Z-score scale.
        p_values: P-values per observation.
        hot_spots: Boolean mask for significant hot spots.
        cold_spots: Boolean mask for significant cold spots.
        significance_level: Alpha level used.
        status: ``"ok"``, or why G* and the p-values of every observation are undefined (NaN).
    """

    g_star: Any  # np.ndarray (n,)
    p_values: Any  # np.ndarray (n,)
    hot_spots: Any  # np.ndarray (n,) bool
    cold_spots: Any  # np.ndarray (n,) bool
    significance_level: float
    status: str = "ok"


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
    nugget: float = float("nan")
    sill: float = float("nan")
    range_param: float = float("nan")


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
        W = (diag_inv @ W).tocsr()

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
        MoransIResult with I statistic, expected value, variance, z-score, p-value. ``I`` is NaN (and the
        variance, z-score and p-value too) with fewer than 2 observations, a weights matrix without edges,
        or constant values; the variance, z-score and p-value are NaN with fewer than 4 observations or a
        variance that is not positive. ``status`` says which.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if norm is None:
        raise ImportError("scipy.stats is required: uv pip install scipy")

    x = np.asarray(values, dtype=np.float64).flatten()
    n = len(x)

    if sp_sparse is not None and sp_sparse.issparse(weights):
        # scipy-stubs narrow issparse() to a method-less SparseABC; every sparse format provides these methods
        W = cast("sp_sparse.csr_matrix", weights)
    else:
        W = sp_sparse.csr_matrix(np.asarray(weights, dtype=np.float64))

    # What is returned when the statistic cannot be computed: NaN, and why
    expected_i = -1.0 / (n - 1) if n > 1 else float("nan")

    def undefined(reason: str) -> MoransIResult:
        nan = float("nan")
        return MoransIResult(I=nan, expected_I=expected_i, variance_I=nan, z_score=nan, p_value=nan, n=n, status=reason)

    if n < 2:
        return undefined("fewer than 2 observations")

    # Center values
    x_mean = x.mean()
    z = x - x_mean

    # Sum of weights
    S0 = float(W.sum())
    if S0 == 0:
        return undefined("the weights matrix has no edges (all weights are zero)")

    # Moran's I
    numerator = float(z @ (W @ z))
    denominator = float(z @ z)

    if denominator == 0:
        return undefined("the values are constant, so Moran's I (0 / 0) is undefined")

    moran_i = (n / S0) * (numerator / denominator)

    # Expected value under null
    E_I = expected_i

    # Variance under normality assumption (randomization)
    S1 = float(((W + W.T).multiply(W + W.T)).sum()) / 2.0
    S2_vec = np.array(W.sum(axis=0)).flatten() + np.array(W.sum(axis=1)).flatten()
    S2 = float((S2_vec**2).sum())

    # Kurtosis term
    b2 = (n * np.sum(z**4)) / (np.sum(z**2) ** 2)

    A = n * ((n**2 - 3 * n + 3) * S1 - n * S2 + 3 * S0**2)
    B = b2 * ((n**2 - n) * S1 - 2 * n * S2 + 6 * S0**2)
    C = (n - 1) * (n - 2) * (n - 3) * S0**2

    status = "ok"
    if C == 0:
        Var_I = float("nan")
        status = "the variance of Moran's I is undefined for fewer than 4 observations"
    else:
        Var_I = (A - B) / C - E_I**2

    # Z-score and p-value (NaN when the variance is missing or not positive)
    if Var_I > 0:
        z_score = (moran_i - E_I) / np.sqrt(Var_I)
        p_value = float(2 * norm.sf(abs(z_score)))
    else:
        z_score = float("nan")
        p_value = float("nan")
        if status == "ok":
            status = "the variance of Moran's I is not positive, so its z-score and p-value are undefined"

    logger.info(f"Moran's I={moran_i:.4f}, z={z_score:.2f}, p={p_value:.4e}")
    return MoransIResult(
        I=moran_i,
        expected_I=E_I,
        variance_I=Var_I,
        z_score=z_score,
        p_value=p_value,
        n=n,
        status=status,
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
        weights: Spatial weights matrix (n x n), dense or sparse in any scipy format.

    Returns:
        GearyCResult with C statistic, z-score, and p-value. ``C`` is NaN (and the variance, z-score and
        p-value too) with fewer than 2 observations, a weights matrix without edges, or constant values; the
        variance, z-score and p-value are NaN with fewer than 4 observations or a variance that is not
        positive. ``status`` says which.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if norm is None:
        raise ImportError("scipy.stats is required: uv pip install scipy")

    x = np.asarray(values, dtype=np.float64).flatten()
    n = len(x)

    if sp_sparse is not None and sp_sparse.issparse(weights):
        W = cast("sp_sparse.csr_matrix", weights)
    else:
        W = sp_sparse.csr_matrix(np.asarray(weights, dtype=np.float64))

    # What is returned when the statistic cannot be computed: NaN, and why (the null expectation of C is 1)
    def undefined(reason: str) -> GearyCResult:
        nan = float("nan")
        return GearyCResult(C=nan, expected_C=1.0, variance_C=nan, z_score=nan, p_value=nan, n=n, status=reason)

    if n < 2:
        return undefined("fewer than 2 observations")

    x_mean = x.mean()
    z = x - x_mean
    S0 = float(W.sum())

    if S0 == 0:
        return undefined("the weights matrix has no edges (all weights are zero)")

    # Numerator: sum of w_ij * (x_i - x_j)^2. Read the stored entries as (row, col, value)
    # triplets of a CSR copy instead of indexing W[sources, targets], which raises for COO, DIA
    # and BSR weights and does not return a plain array for LIL and DOK weights.
    W_coo = W.tocsr().tocoo()
    sources, targets, w_vals = W_coo.row, W_coo.col, W_coo.data
    sq_diffs = (x[sources] - x[targets]) ** 2
    numerator = float(np.sum(w_vals * sq_diffs))

    # Denominator: sum of (x_i - x_bar)^2
    denominator = float(np.sum(z**2))

    if denominator == 0:
        return undefined("the values are constant, so Geary's C (0 / 0) is undefined")

    C = ((n - 1) / (2.0 * S0)) * (numerator / denominator)

    # Expected C under null
    E_C = 1.0

    # Variance under normality assumption
    S1 = float(((W + W.T).multiply(W + W.T)).sum()) / 2.0
    S2_vec = np.array(W.sum(axis=0)).flatten() + np.array(W.sum(axis=1)).flatten()
    S2 = float((S2_vec**2).sum())

    n_f = float(n)
    status = "ok"
    if n > 3:
        Var_C = ((2 * S1 + S2) * (n_f - 1) - 4 * S0**2) / (2 * (n_f + 1) * S0**2)
    else:
        Var_C = float("nan")
        status = "the variance of Geary's C is undefined for fewer than 4 observations"

    # Z-score and p-value (NaN when the variance is missing or not positive)
    if Var_C > 0:
        z_score = (C - E_C) / np.sqrt(Var_C)
        p_value = float(2 * norm.sf(abs(z_score)))
    else:
        z_score = float("nan")
        p_value = float("nan")
        if status == "ok":
            status = "the variance of Geary's C is not positive, so its z-score and p-value are undefined"

    logger.info(f"Geary's C={C:.4f}, z={z_score:.2f}, p={p_value:.4e}")
    return GearyCResult(
        C=C,
        expected_C=E_C,
        variance_C=Var_C,
        z_score=z_score,
        p_value=p_value,
        n=n,
        status=status,
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
        LocalMoransResult with local I values, z-scores, p-values, and cluster labels. The z-scores and
        p-values come from a normal approximation with a simplified variance (see the code comment). They are
        NaN for an observation without a positive variance (one with no neighbor), and everything is NaN
        (``status`` says why) with fewer than 3 observations or constant values; such observations are "NS".
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if norm is None:
        raise ImportError("scipy.stats is required: uv pip install scipy")

    x = np.asarray(values, dtype=np.float64).flatten()
    n = len(x)

    if sp_sparse is not None and sp_sparse.issparse(weights):
        W = cast("sp_sparse.csr_matrix", weights)
    else:
        W = sp_sparse.csr_matrix(np.asarray(weights, dtype=np.float64))

    # What is returned when the statistics cannot be computed: NaN for every observation, and why
    expected_i = -1.0 / (n - 1) if n > 1 else float("nan")

    def undefined(reason: str) -> LocalMoransResult:
        return LocalMoransResult(
            local_I=np.full(n, np.nan),
            expected_I=expected_i,
            z_scores=np.full(n, np.nan),
            p_values=np.full(n, np.nan),
            cluster_labels=["NS"] * n,
            significance_level=significance,
            status=reason,
        )

    if n < 3:
        return undefined("fewer than 3 observations")

    x_mean = x.mean()
    z = x - x_mean

    # Variance of z
    m2 = float(np.sum(z**2) / n)
    if m2 == 0:
        return undefined("the values are constant, so the local statistics (0 / 0) are undefined")

    # Local Moran's I for each observation
    # I_i = (z_i / m2) * sum_j(w_ij * z_j)
    Wz = np.array(W @ z).flatten()
    local_I = (z / m2) * Wz

    # Expected value
    E_Ii = expected_i

    # Variance of local I: a normal approximation with a simplified variance. It keeps Anselin's (1995)
    # randomization form, but uses (sum_j w_ij)^2 where his cross-product term has 2 * sum_{k != h} w_ik w_ih.
    w_i_sq = np.array(W.multiply(W).sum(axis=1)).flatten()
    w_i_sum = np.array(W.sum(axis=1)).flatten()

    b2 = (np.sum(z**4) / n) / (m2**2)

    Var_Ii = (w_i_sq * (n - b2) / (n - 1) + w_i_sum**2 * (2 * b2 - n) / ((n - 1) * (n - 2)) - E_Ii**2).astype(
        np.float64
    )

    # Z-scores and p-values exist only where the variance is positive (an observation without any neighbor has
    # none); elsewhere they are NaN and the observation is "NS"
    z_scores = np.full(n, np.nan)
    p_values = np.full(n, np.nan)
    positive = Var_Ii > 0
    z_scores[positive] = (local_I[positive] - E_Ii) / np.sqrt(Var_Ii[positive])
    p_values[positive] = 2 * norm.sf(np.abs(z_scores[positive]))

    # Spatial lag (mean of neighbors)
    spatial_lag = Wz  # already W @ z

    # Classify clusters
    cluster_labels: list[str] = []
    for i in range(n):
        if not p_values[i] <= significance:  # also true for a NaN p-value
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

    n_sig = sum(1 for label in cluster_labels if label != "NS")
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
        GetisOrdResult with G* statistics and hot/cold spot identification. G* and the p-value are NaN for an
        observation that weighs all observations equally (its weights have no spread to standardize by), and
        for every observation (``status`` says why) with fewer than 2 observations or constant values; such
        observations are neither hot nor cold spots.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if norm is None:
        raise ImportError("scipy.stats is required: uv pip install scipy")

    x = np.asarray(values, dtype=np.float64).flatten()
    n = len(x)

    if sp_sparse is not None and sp_sparse.issparse(weights):
        W = cast("sp_sparse.csr_matrix", weights)
    else:
        W = sp_sparse.csr_matrix(np.asarray(weights, dtype=np.float64))

    # What is returned when G* cannot be computed: NaN for every observation, and why
    def undefined(reason: str) -> GetisOrdResult:
        return GetisOrdResult(
            g_star=np.full(n, np.nan),
            p_values=np.full(n, np.nan),
            hot_spots=np.zeros(n, dtype=bool),
            cold_spots=np.zeros(n, dtype=bool),
            significance_level=significance,
            status=reason,
        )

    if n < 2:
        return undefined("fewer than 2 observations")

    x_bar = x.mean()
    S = np.sqrt(np.sum((x - x_bar) ** 2) / n)

    if S == 0:
        return undefined("the values are constant, so G* (0 / 0) is undefined")

    # Include self-weight (G* vs G)
    # Add identity to weights for G*
    W_star = W + sp_sparse.eye(n)

    W_star_sum = np.array(W_star.sum(axis=1)).flatten()
    W_star_sq_sum = np.array(W_star.multiply(W_star).sum(axis=1)).flatten()
    Wx = np.array(W_star @ x).flatten()

    # G* is a z-score only where the denominator is positive; elsewhere (an observation that weighs all n
    # observations equally, so its weights have no spread) it is undefined: NaN, not 0
    g_star = np.full(n, np.nan)
    numerator = Wx - x_bar * W_star_sum
    denominator_inner = (n * W_star_sq_sum - W_star_sum**2) / (n - 1)
    defined = denominator_inner > 0
    g_star[defined] = numerator[defined] / (S * np.sqrt(denominator_inner[defined]))

    # P-values (two-sided normal); NaN where G* is undefined
    p_values = 2 * norm.sf(np.abs(g_star))

    # Identify hot/cold spots (a NaN comparison is False, so undefined observations are neither)
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
        VariogramResult with semivariance values and estimated parameters. A bin without any pair has a
        NaN semivariance (``n_pairs`` is 0 there); ``nugget``, ``sill`` and ``range_param`` are NaN when they
        cannot be estimated (no pair at all, or a sill of zero for the range), and for fewer than 2 points.
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
            nugget=float("nan"),
            sill=float("nan"),
            range_param=float("nan"),
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

    # Average semivariance per bin: gamma(h) = (1/2N(h)) * sum (x_i - x_j)^2. A bin without pairs has no
    # semivariance: NaN, not a zero that would read as identical values
    valid = pair_counts > 0
    semivariance[valid] = semivariance[valid] / (2.0 * pair_counts[valid])
    semivariance[~valid] = np.nan

    # Estimate variogram parameters (NaN when no bin has a pair, or when the sill is not reached/positive)
    nugget = float(semivariance[valid][0]) if valid.any() else float("nan")
    sill = float(np.max(semivariance[valid])) if valid.any() else float("nan")

    # Range: first bin where semivariance reaches ~95% of sill
    range_param = float("nan")
    if sill > 0:
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
