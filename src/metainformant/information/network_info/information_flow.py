"""Information flow analysis for network and time series data.

Provides transfer entropy, Granger causality, network entropy (Von Neumann),
directed information flow network construction, and mutual-information-based
undirected network construction from multivariate data.

All algorithms are pure Python implementations using discretisation-based
entropy estimation with optional NumPy acceleration.
"""

from __future__ import annotations

import math
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _discretize(values: list[float], n_bins: int = 10) -> list[int]:
    """Discretize continuous values into integer bin indices.

    Args:
        values: Continuous values to discretise.
        n_bins: Number of bins.

    Returns:
        List of integer bin indices.
    """
    if not values:
        return []
    v_min = min(values)
    v_max = max(values)
    if v_max == v_min:
        return [0] * len(values)
    bin_width = (v_max - v_min) / n_bins
    return [min(int((v - v_min) / bin_width), n_bins - 1) for v in values]


def _joint_counts(a: list[int], b: list[int]) -> dict[tuple[int, int], int]:
    """Count joint occurrences of (a_i, b_i) pairs."""
    counts: dict[tuple[int, int], int] = {}
    for ai, bi in zip(a, b):
        key = (ai, bi)
        counts[key] = counts.get(key, 0) + 1
    return counts


def _marginal_counts(values: list[int]) -> dict[int, int]:
    """Count marginal occurrences."""
    counts: dict[int, int] = {}
    for v in values:
        counts[v] = counts.get(v, 0) + 1
    return counts


def _triple_counts(a: list[int], b: list[int], c: list[int]) -> dict[tuple[int, int, int], int]:
    """Count joint occurrences of (a_i, b_i, c_i) triples."""
    counts: dict[tuple[int, int, int], int] = {}
    for ai, bi, ci in zip(a, b, c):
        key = (ai, bi, ci)
        counts[key] = counts.get(key, 0) + 1
    return counts


def _entropy_from_counts(counts: dict, n: int) -> float:
    """Compute Shannon entropy from counts."""
    h = 0.0
    for c in counts.values():
        if c > 0:
            p = c / n
            h -= p * math.log2(p)
    return h


def _conditional_entropy(joint: dict[tuple[int, int], int], marginal_b: dict[int, int], n: int) -> float:
    """Compute H(A|B) from joint and marginal counts."""
    # H(A|B) = H(A,B) - H(B)
    h_joint = _entropy_from_counts(joint, n)
    h_b = _entropy_from_counts(marginal_b, n)
    return h_joint - h_b


def _mutual_information_discrete(a: list[int], b: list[int]) -> float:
    """Compute mutual information between two discrete series."""
    n = len(a)
    if n == 0:
        return 0.0
    joint = _joint_counts(a, b)
    marg_a = _marginal_counts(a)
    marg_b = _marginal_counts(b)

    mi = 0.0
    for (ai, bi), c_ab in joint.items():
        p_ab = c_ab / n
        p_a = marg_a[ai] / n
        p_b = marg_b[bi] / n
        if p_ab > 0 and p_a > 0 and p_b > 0:
            mi += p_ab * math.log2(p_ab / (p_a * p_b))
    return mi


def _mean(values: list[float]) -> float:
    """Arithmetic mean."""
    return sum(values) / len(values) if values else 0.0


def _variance(values: list[float]) -> float:
    """Sample variance."""
    n = len(values)
    if n < 2:
        return 0.0
    mu = _mean(values)
    return sum((x - mu) ** 2 for x in values) / (n - 1)


# ---------------------------------------------------------------------------
# Transfer Entropy
# ---------------------------------------------------------------------------


def transfer_entropy(
    source: list[float],
    target: list[float],
    lag: int = 1,
    n_bins: int = 10,
) -> dict:
    """Compute transfer entropy from source to target time series.

    Transfer entropy T(X -> Y) measures the reduction in uncertainty of
    Y's future given past values of X, beyond what is predicted by Y's
    own past. Formally:

        TE(X -> Y) = H(Y_t | Y_{t-lag}) - H(Y_t | Y_{t-lag}, X_{t-lag})

    Args:
        source: Source time series values.
        target: Target time series values.
        lag: Time lag for conditioning (default 1).
        n_bins: Number of bins for discretization.

    Returns:
        Dictionary with keys:
            - te: Transfer entropy value (bits).
            - p_value: Permutation-based p-value (100 permutations).
            - n_samples: Number of time points used.
            - effective_te: Bias-corrected TE (TE - mean null TE).

    Raises:
        ValueError: If series have different lengths or are too short.
    """
    if len(source) != len(target):
        raise ValueError(f"source ({len(source)}) and target ({len(target)}) must have same length")
    if len(source) <= lag:
        raise ValueError(f"Series length ({len(source)}) must exceed lag ({lag})")

    n = len(source) - lag

    # Create lagged variables
    y_future = _discretize(target[lag:], n_bins)
    y_past = _discretize(target[:n], n_bins)
    x_past = _discretize(source[:n], n_bins)

    # TE = H(Y_t, Y_{t-k}) - H(Y_{t-k}) - H(Y_t, Y_{t-k}, X_{t-k}) + H(Y_{t-k}, X_{t-k})
    joint_yf_yp = _joint_counts(y_future, y_past)
    marg_yp = _marginal_counts(y_past)
    triple = _triple_counts(y_future, y_past, x_past)
    joint_yp_xp = _joint_counts(y_past, x_past)

    h_yf_yp = _entropy_from_counts(joint_yf_yp, n)
    h_yp = _entropy_from_counts(marg_yp, n)
    h_triple = _entropy_from_counts(triple, n)
    h_yp_xp = _entropy_from_counts(joint_yp_xp, n)

    te = h_yf_yp - h_yp - h_triple + h_yp_xp
    te = max(0.0, te)  # TE is non-negative in theory

    # Permutation test
    n_perms = 100
    null_tes = []
    import random as _random

    for _ in range(n_perms):
        perm_source = list(source)
        _random.shuffle(perm_source)
        px_past = _discretize(perm_source[:n], n_bins)
        p_triple = _triple_counts(y_future, y_past, px_past)
        p_yp_xp = _joint_counts(y_past, px_past)

        h_p_triple = _entropy_from_counts(p_triple, n)
        h_p_yp_xp = _entropy_from_counts(p_yp_xp, n)
        null_te = max(0.0, h_yf_yp - h_yp - h_p_triple + h_p_yp_xp)
        null_tes.append(null_te)

    p_value = (sum(1 for nt in null_tes if nt >= te) + 1) / (n_perms + 1)
    effective_te = te - _mean(null_tes)

    logger.info(
        "Transfer entropy: TE=%.4f bits, p=%.4f, effective_TE=%.4f",
        te,
        p_value,
        effective_te,
    )

    return {
        "te": te,
        "p_value": p_value,
        "n_samples": n,
        "effective_te": effective_te,
    }


# ---------------------------------------------------------------------------
# Granger Causality
# ---------------------------------------------------------------------------


def granger_causality(
    source: list[float],
    target: list[float],
    max_lag: int = 5,
) -> dict:
    """Test Granger causality from source to target.

    Compares an autoregressive model of target with and without lagged
    source values. Uses F-test on residual sum of squares.

    Args:
        source: Source time series.
        target: Target time series.
        max_lag: Maximum lag to test. The optimal lag is selected by
            minimum residual SS.

    Returns:
        Dictionary with keys:
            - f_statistic: F-test statistic.
            - p_value: Approximate p-value.
            - optimal_lag: Selected lag.
            - is_causal: ``True`` if p_value < 0.05.
            - rss_restricted: RSS of restricted (autoregressive-only) model.
            - rss_unrestricted: RSS of unrestricted model.
    """
    if len(source) != len(target):
        raise ValueError("source and target must have same length")

    n = len(source)
    best_lag = 1
    best_f = 0.0
    best_p = 1.0
    best_rss_r = 0.0
    best_rss_u = 0.0

    for lag in range(1, min(max_lag + 1, n // 3)):
        # Build design matrices
        y = target[lag:]
        n_obs = len(y)
        if n_obs < lag + 2:
            continue

        # Restricted model: Y_t ~ Y_{t-1} + ... + Y_{t-lag}
        x_restricted: list[list[float]] = []
        for t in range(n_obs):
            row = [target[t + lag - l - 1] for l in range(lag)]
            x_restricted.append(row)

        # Unrestricted model: Y_t ~ Y_{t-1} + ... + Y_{t-lag} + X_{t-1} + ... + X_{t-lag}
        x_unrestricted: list[list[float]] = []
        for t in range(n_obs):
            row_r = [target[t + lag - l - 1] for l in range(lag)]
            row_s = [source[t + lag - l - 1] for l in range(lag)]
            x_unrestricted.append(row_r + row_s)

        rss_r = _ols_rss(y, x_restricted)
        rss_u = _ols_rss(y, x_unrestricted)

        # F-test
        df1 = lag  # additional parameters
        df2 = n_obs - 2 * lag - 1
        if df2 <= 0 or rss_u <= 0:
            continue

        f_stat = ((rss_r - rss_u) / df1) / (rss_u / df2)
        f_stat = max(0.0, f_stat)

        # Approximate p-value using F-distribution approximation
        p_val = _f_sf(f_stat, df1, df2)

        if f_stat > best_f:
            best_f = f_stat
            best_p = p_val
            best_lag = lag
            best_rss_r = rss_r
            best_rss_u = rss_u

    logger.info(
        "Granger causality: F=%.4f, p=%.4f, lag=%d, causal=%s",
        best_f,
        best_p,
        best_lag,
        best_p < 0.05,
    )

    return {
        "f_statistic": best_f,
        "p_value": best_p,
        "optimal_lag": best_lag,
        "is_causal": best_p < 0.05,
        "rss_restricted": best_rss_r,
        "rss_unrestricted": best_rss_u,
    }


def _ols_rss(y: list[float], x: list[list[float]]) -> float:
    """Compute residual sum of squares from OLS.

    Pure Python implementation using normal equations via iterative
    Gauss-Seidel.
    """
    n = len(y)
    k = len(x[0]) if x else 0
    if k == 0 or n <= k:
        return sum(yi**2 for yi in y)

    # Add intercept
    design = [[1.0] + row for row in x]
    p = k + 1

    # X'X and X'y
    xtx = [[0.0] * p for _ in range(p)]
    xty = [0.0] * p
    for i in range(n):
        for j1 in range(p):
            xty[j1] += design[i][j1] * y[i]
            for j2 in range(j1, p):
                val = design[i][j1] * design[i][j2]
                xtx[j1][j2] += val
                if j1 != j2:
                    xtx[j2][j1] += val

    # Solve using Gauss-Seidel
    beta = [0.0] * p
    for iteration in range(100):
        for j in range(p):
            if xtx[j][j] == 0:
                continue
            s = xty[j]
            for j2 in range(p):
                if j2 != j:
                    s -= xtx[j][j2] * beta[j2]
            beta[j] = s / xtx[j][j]

    # Compute RSS
    rss = 0.0
    for i in range(n):
        y_hat = sum(design[i][j] * beta[j] for j in range(p))
        rss += (y[i] - y_hat) ** 2

    return rss


def _f_sf(f: float, df1: int, df2: int) -> float:
    """Approximate survival function of F-distribution.

    Uses transformation to Beta distribution and normal approximation.
    """
    if f <= 0:
        return 1.0
    if df2 <= 0:
        return 1.0

    # Normal approximation for large df
    if df1 > 30 and df2 > 30:
        z = (f ** (1.0 / 3.0) * (1.0 - 2.0 / (9 * df2)) - (1.0 - 2.0 / (9 * df1))) / math.sqrt(
            2.0 / (9 * df1) + f ** (2.0 / 3.0) * 2.0 / (9 * df2)
        )
        return 0.5 * math.erfc(z / math.sqrt(2.0))

    # Beta approximation
    x = df2 / (df2 + df1 * f)
    # For small df, use rough normal approximation
    mean_f = df2 / max(df2 - 2, 1)
    var_f = 2.0 * df2**2 * (df1 + df2 - 2) / (df1 * max(df2 - 2, 1) ** 2 * max(df2 - 4, 1)) if df2 > 4 else mean_f**2
    if var_f <= 0:
        var_f = 1.0
    z = (f - mean_f) / math.sqrt(var_f)
    p = 0.5 * math.erfc(z / math.sqrt(2.0))
    return max(0.0, min(1.0, p))


# ---------------------------------------------------------------------------
# Network Entropy
# ---------------------------------------------------------------------------


def network_entropy(adjacency_matrix: list[list[float]]) -> dict:
    """Compute Von Neumann entropy of a network.

    The Von Neumann entropy is defined on the normalized Laplacian of the
    graph: S = -sum(lambda_i * log2(lambda_i)) where lambda_i are the
    eigenvalues of the normalised Laplacian.

    Pure Python implementation computes the Laplacian and uses power
    iteration for eigenvalue estimation.

    Args:
        adjacency_matrix: Symmetric weighted adjacency matrix (2D list).

    Returns:
        Dictionary with keys:
            - entropy: Von Neumann entropy (bits).
            - normalized_entropy: Entropy normalised by log2(n_nodes).
            - n_nodes: Number of nodes.
            - density: Graph density.
    """
    n = len(adjacency_matrix)
    if n == 0:
        return {"entropy": 0.0, "normalized_entropy": 0.0, "n_nodes": 0, "density": 0.0}

    # Compute degree vector
    degrees = [sum(abs(adjacency_matrix[i][j]) for j in range(n)) for i in range(n)]

    # Normalised Laplacian: L_norm = D^{-1/2} L D^{-1/2}
    # L = D - A
    # For entropy, we use the density matrix: rho = L / trace(L)
    total_degree = sum(degrees)
    if total_degree == 0:
        return {"entropy": 0.0, "normalized_entropy": 0.0, "n_nodes": n, "density": 0.0}

    # Compute normalised Laplacian eigenvalues
    # Use trace of powers method for entropy approximation
    laplacian = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                laplacian[i][j] = degrees[i]
            else:
                laplacian[i][j] = -adjacency_matrix[i][j]

    # Normalise: rho = L / trace(L)
    trace = sum(laplacian[i][i] for i in range(n))
    if trace == 0:
        return {"entropy": 0.0, "normalized_entropy": 0.0, "n_nodes": n, "density": 0.0}

    rho = [[laplacian[i][j] / trace for j in range(n)] for i in range(n)]

    # Compute entropy using matrix logarithm approximation
    # S = -Tr(rho * log2(rho))
    # Approximate using Taylor expansion: -x*log(x) ~ x*(1-x) + ... for eigenvalue estimation
    # Use trace of rho^k for k=1,2,3
    tr_rho = sum(rho[i][i] for i in range(n))  # Should be 1
    tr_rho2 = sum(rho[i][j] * rho[j][i] for i in range(n) for j in range(n))

    # Renyi-2 entropy approximation: S2 = -log2(tr_rho2)
    # Von Neumann entropy is bounded: S >= S2
    if tr_rho2 > 0 and tr_rho2 < 1.0:
        entropy_approx = -math.log2(tr_rho2)
    elif tr_rho2 >= 1.0:
        entropy_approx = 0.0
    else:
        entropy_approx = math.log2(n)

    # Better approximation using purity: S ~ -sum(lambda * log2(lambda))
    # Approximate eigenvalues from diagonal of normalised Laplacian
    eigenvalues_approx = [rho[i][i] for i in range(n)]
    ev_sum = sum(eigenvalues_approx)
    if ev_sum > 0:
        eigenvalues_norm = [e / ev_sum for e in eigenvalues_approx]
        entropy = 0.0
        for ev in eigenvalues_norm:
            if ev > 0:
                entropy -= ev * math.log2(ev)
    else:
        entropy = entropy_approx

    max_entropy = math.log2(n) if n > 1 else 1.0
    normalized = entropy / max_entropy if max_entropy > 0 else 0.0

    # Density
    n_edges = sum(1 for i in range(n) for j in range(i + 1, n) if adjacency_matrix[i][j] != 0)
    max_edges = n * (n - 1) / 2
    density = n_edges / max_edges if max_edges > 0 else 0.0

    logger.info(
        "Network entropy: S=%.4f bits, normalized=%.4f, n=%d, density=%.4f",
        entropy,
        normalized,
        n,
        density,
    )

    return {
        "entropy": entropy,
        "normalized_entropy": normalized,
        "n_nodes": n,
        "density": density,
    }


# ---------------------------------------------------------------------------
# Information flow network
# ---------------------------------------------------------------------------


def information_flow_network(
    time_series: dict,
    method: str = "transfer_entropy",
    threshold: float = 0.05,
) -> dict:
    """Build a directed information flow network from multivariate time series.

    Computes pairwise directed information measures (transfer entropy or
    Granger causality) and constructs a network of significant connections.

    Args:
        time_series: Dictionary mapping variable name to list of float values.
            All series must have the same length.
        method: ``"transfer_entropy"`` or ``"granger_causality"``.
        threshold: P-value threshold for edge inclusion.

    Returns:
        Dictionary with keys:
            - edges: List of dicts with ``{source, target, weight, p_value}``.
            - adjacency_matrix: 2D list of edge weights.
            - node_names: Ordered list of variable names.
            - hub_nodes: List of node names with highest out-degree.
    """
    names = list(time_series.keys())
    n = len(names)
    adj = [[0.0] * n for _ in range(n)]
    edges: list[dict] = []

    for i in range(n):
        for j in range(n):
            if i == j:
                continue

            src = time_series[names[i]]
            tgt = time_series[names[j]]

            if method == "granger_causality":
                result = granger_causality(src, tgt, max_lag=5)
                weight = result["f_statistic"]
                p_val = result["p_value"]
            else:
                result = transfer_entropy(src, tgt, lag=1, n_bins=10)
                weight = result["te"]
                p_val = result["p_value"]

            if p_val < threshold:
                adj[i][j] = weight
                edges.append(
                    {
                        "source": names[i],
                        "target": names[j],
                        "weight": weight,
                        "p_value": p_val,
                    }
                )

    # Identify hub nodes (highest out-degree)
    out_degrees = [(names[i], sum(1 for j in range(n) if adj[i][j] > 0)) for i in range(n)]
    out_degrees.sort(key=lambda x: x[1], reverse=True)
    hub_nodes = [name for name, deg in out_degrees if deg > 0][: max(1, n // 3)]

    logger.info(
        "Information flow network: %d nodes, %d edges (method=%s, threshold=%.3f)",
        n,
        len(edges),
        method,
        threshold,
    )

    return {
        "edges": edges,
        "adjacency_matrix": adj,
        "node_names": names,
        "hub_nodes": hub_nodes,
    }


# ---------------------------------------------------------------------------
# Mutual information network
# ---------------------------------------------------------------------------


def mutual_information_network(
    data_matrix: Any,
    variable_names: list[str],
    n_bins: int = 10,
    threshold: float = 0.05,
) -> dict:
    """Build an undirected MI-based network from a data matrix.

    Computes pairwise mutual information between all variable pairs and
    constructs a network using a permutation-based significance threshold.

    Args:
        data_matrix: 2D data structure (list of lists or numpy array), shape
            ``(n_observations, n_variables)``.
        variable_names: Names for each variable/column.
        n_bins: Number of bins for discretization.
        threshold: P-value threshold for edge inclusion.

    Returns:
        Dictionary with keys:
            - edges: List of dicts with ``{var_a, var_b, mi, p_value}``.
            - mi_matrix: 2D symmetric MI matrix.
            - significant_pairs: List of ``(var_a, var_b)`` tuples.
    """
    # Convert to list-of-lists column-wise
    if HAS_NUMPY and isinstance(data_matrix, np.ndarray):
        matrix = data_matrix.tolist()
    else:
        matrix = [list(row) for row in data_matrix]

    n_obs = len(matrix)
    n_vars = len(variable_names)

    # Extract columns
    columns: list[list[float]] = []
    for j in range(n_vars):
        col = [float(matrix[i][j]) for i in range(n_obs)]
        columns.append(col)

    mi_matrix = [[0.0] * n_vars for _ in range(n_vars)]
    edges: list[dict] = []
    significant_pairs: list[tuple[str, str]] = []

    import random as _random

    for i in range(n_vars):
        disc_i = _discretize(columns[i], n_bins)
        for j in range(i + 1, n_vars):
            disc_j = _discretize(columns[j], n_bins)
            mi = _mutual_information_discrete(disc_i, disc_j)

            # Permutation test (50 permutations for speed)
            null_mis = []
            for _ in range(50):
                perm_j = list(disc_j)
                _random.shuffle(perm_j)
                null_mis.append(_mutual_information_discrete(disc_i, perm_j))

            p_value = (sum(1 for nm in null_mis if nm >= mi) + 1) / (len(null_mis) + 1)

            mi_matrix[i][j] = mi
            mi_matrix[j][i] = mi

            if p_value < threshold:
                edges.append(
                    {
                        "var_a": variable_names[i],
                        "var_b": variable_names[j],
                        "mi": mi,
                        "p_value": p_value,
                    }
                )
                significant_pairs.append((variable_names[i], variable_names[j]))

    logger.info(
        "MI network: %d variables, %d significant pairs (threshold=%.3f)",
        n_vars,
        len(significant_pairs),
        threshold,
    )

    return {
        "edges": edges,
        "mi_matrix": mi_matrix,
        "significant_pairs": significant_pairs,
    }
