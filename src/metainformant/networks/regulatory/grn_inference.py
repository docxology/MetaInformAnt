"""Gene regulatory network inference from expression data.

This module provides multiple methods for inferring gene regulatory networks
(GRNs) from gene expression matrices, including correlation-based, mutual
information-based (ARACNE-like), and regression-based approaches. It also
provides tools for scoring transcription factor activity, detecting network
motifs, and validating inferred networks against known interactions.
"""

from __future__ import annotations

import math
import random
from collections import defaultdict
from typing import Any, cast

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependencies
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = cast("Any", None)

try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = cast("Any", None)


def _pearson_correlation_pure(x: list[float], y: list[float]) -> tuple[float, float]:
    """Compute Pearson correlation coefficient using pure Python.

    Args:
        x: First variable values.
        y: Second variable values.

    Returns:
        Tuple of (correlation, p_value).
    """
    n = len(x)
    if n < 3:
        return 0.0, 1.0

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    cov_xy = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    var_x = sum((xi - mean_x) ** 2 for xi in x)
    var_y = sum((yi - mean_y) ** 2 for yi in y)

    denom = math.sqrt(var_x * var_y)
    if denom < 1e-15:
        return 0.0, 1.0

    r = cov_xy / denom

    # t-statistic for p-value approximation
    if abs(r) >= 1.0:
        return float(max(-1.0, min(1.0, r))), 0.0

    t_stat = r * math.sqrt((n - 2) / (1.0 - r * r))
    # Rough two-sided p-value approximation using normal for large n
    p_value = 2.0 * math.exp(-0.5 * t_stat * t_stat) / math.sqrt(2.0 * math.pi)
    p_value = max(0.0, min(1.0, p_value))

    return float(r), float(p_value)


def _spearman_correlation_pure(x: list[float], y: list[float]) -> tuple[float, float]:
    """Compute Spearman rank correlation using pure Python.

    Args:
        x: First variable values.
        y: Second variable values.

    Returns:
        Tuple of (correlation, p_value).
    """
    n = len(x)
    if n < 3:
        return 0.0, 1.0

    def _rank(vals: list[float]) -> list[float]:
        indexed = sorted(enumerate(vals), key=lambda t: t[1])
        ranks = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j < n - 1 and indexed[j + 1][1] == indexed[j][1]:
                j += 1
            avg_rank = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                ranks[indexed[k][0]] = avg_rank
            i = j + 1
        return ranks

    rank_x = _rank(x)
    rank_y = _rank(y)
    return _pearson_correlation_pure(rank_x, rank_y)


def _mutual_information_pure(x: list[float], y: list[float], n_bins: int = 10) -> float:
    """Compute mutual information between two variables using pure Python.

    Uses histogram-based estimation with uniform binning.

    Args:
        x: First variable values.
        y: Second variable values.
        n_bins: Number of bins for discretization.

    Returns:
        Mutual information in nats.
    """
    n = len(x)
    if n == 0:
        return 0.0

    min_x, max_x = min(x), max(x)
    min_y, max_y = min(y), max(y)

    range_x = max_x - min_x if max_x > min_x else 1.0
    range_y = max_y - min_y if max_y > min_y else 1.0

    # Bin assignments
    def _bin_idx(val: float, vmin: float, vrange: float) -> int:
        idx = int((val - vmin) / vrange * n_bins)
        return max(0, min(n_bins - 1, idx))

    # Joint and marginal counts
    joint: dict[tuple[int, int], int] = defaultdict(int)
    margin_x: dict[int, int] = defaultdict(int)
    margin_y: dict[int, int] = defaultdict(int)

    for xi, yi in zip(x, y):
        bx = _bin_idx(xi, min_x, range_x)
        by = _bin_idx(yi, min_y, range_y)
        joint[(bx, by)] += 1
        margin_x[bx] += 1
        margin_y[by] += 1

    # MI = sum p(x,y) * log(p(x,y) / (p(x)*p(y)))
    mi = 0.0
    for (bx, by), count in joint.items():
        p_xy = count / n
        p_x = margin_x[bx] / n
        p_y = margin_y[by] / n
        if p_xy > 0 and p_x > 0 and p_y > 0:
            mi += p_xy * math.log(p_xy / (p_x * p_y))

    return max(0.0, mi)


def _to_list_of_lists(matrix: Any, n_rows: int, n_cols: int) -> list[list[float]]:
    """Convert a matrix (numpy array or list of lists) to list of lists.

    Args:
        matrix: Input matrix.
        n_rows: Expected number of rows.
        n_cols: Expected number of columns.

    Returns:
        Matrix as list of lists of floats.
    """
    if HAS_NUMPY and isinstance(matrix, np.ndarray):
        return [[float(matrix[i, j]) for j in range(n_cols)] for i in range(n_rows)]
    return [[float(matrix[i][j]) for j in range(n_cols)] for i in range(n_rows)]


def _get_matrix_shape(matrix: Any) -> tuple[int, int]:
    """Get the shape of a matrix (numpy array or list of lists).

    Args:
        matrix: Input matrix.

    Returns:
        Tuple of (n_rows, n_cols).
    """
    if HAS_NUMPY and isinstance(matrix, np.ndarray):
        return int(matrix.shape[0]), int(matrix.shape[1])
    n_rows = len(matrix)
    n_cols = len(matrix[0]) if n_rows > 0 else 0
    return n_rows, n_cols


def _resolve_tf_indices(tf_list: list[str] | None, gene_names: list[str]) -> list[int]:
    """Row indices of the regulators: every gene when ``tf_list`` is ``None``, otherwise the listed genes.

    Entries of ``tf_list`` that are not in ``gene_names`` are ignored (a genome-wide TF catalogue is routinely longer
    than the measured genes) and repeated entries count once. An earlier version silently replaced a ``tf_list`` that
    matched nothing (a misspelt or wrongly-cased catalogue, an empty list) with *all* genes, returning the full
    all-against-all network where the caller had asked for a regulator-restricted one; repeated entries produced
    duplicate edges.

    Raises:
        ValueError: If ``tf_list`` is given but none of its entries is a gene in ``gene_names``.
    """
    if tf_list is None:
        return list(range(len(gene_names)))
    gene_to_idx = {name: i for i, name in enumerate(gene_names)}
    indices = list(dict.fromkeys(gene_to_idx[tf] for tf in tf_list if tf in gene_to_idx))
    if not indices:
        preview = ", ".join(repr(tf) for tf in list(tf_list)[:5]) or "an empty list"
        raise ValueError(
            f"none of the regulators in tf_list ({preview}) is in gene_names; pass tf_list=None to use every gene "
            "as a potential regulator"
        )
    return indices


def infer_grn_correlation(
    expression_matrix: Any,
    gene_names: list[str],
    tf_list: list[str] | None = None,
    method: str = "pearson",
    threshold: float = 0.3,
) -> dict:
    """Infer gene regulatory network from expression correlation.

    Computes pairwise correlation between transcription factors (TFs) and
    target genes, filtering edges by absolute correlation above threshold.

    Args:
        expression_matrix: Gene expression matrix (genes x samples). Can be
            a numpy array or list of lists.
        gene_names: List of gene names corresponding to matrix rows.
        tf_list: List of transcription factor gene names. If None, all genes
            are considered potential regulators. Names that are not in
            ``gene_names`` are ignored and repeated names count once; a list
            that matches no gene raises ``ValueError`` instead of silently
            widening to all genes.
        method: Correlation method, one of 'pearson' or 'spearman'.
        threshold: Minimum absolute correlation to include an edge.

    Returns:
        Dictionary containing:
            - edges: List of dicts with source, target, weight, p_value.
            - n_edges: Number of edges in the network.
            - n_tfs: Number of transcription factors with edges.
            - n_targets: Number of target genes with edges.
            - adjacency_matrix: Correlation matrix as list of lists.

    Raises:
        ValueError: If gene_names length does not match matrix rows, method
            is not supported, or ``tf_list`` matches no gene.
    """
    n_genes, n_samples = _get_matrix_shape(expression_matrix)

    if len(gene_names) != n_genes:
        raise ValueError(f"gene_names length ({len(gene_names)}) must match matrix rows ({n_genes})")

    if method not in ("pearson", "spearman"):
        raise ValueError(f"Unsupported method: {method}. Use 'pearson' or 'spearman'.")

    logger.info(
        "Inferring GRN via %s correlation: %d genes, %d samples, threshold=%.3f",
        method,
        n_genes,
        n_samples,
        threshold,
    )

    data = _to_list_of_lists(expression_matrix, n_genes, n_samples)
    tf_indices = _resolve_tf_indices(tf_list, gene_names)

    # Select correlation function
    if method == "pearson":
        if HAS_SCIPY:

            def _corr(x: list[float], y: list[float]) -> tuple[float, float]:
                r, p = scipy_stats.pearsonr(x, y)
                return float(r), float(p)

        else:
            _corr = _pearson_correlation_pure
    else:
        if HAS_SCIPY:

            def _corr(x: list[float], y: list[float]) -> tuple[float, float]:
                r, p = scipy_stats.spearmanr(x, y)
                return float(r), float(p)

        else:
            _corr = _spearman_correlation_pure

    # Build adjacency and collect edges
    adjacency: list[list[float]] = [[0.0] * n_genes for _ in range(n_genes)]
    edges: list[dict] = []
    active_tfs: set[str] = set()
    active_targets: set[str] = set()

    for tf_idx in tf_indices:
        for target_idx in range(n_genes):
            if tf_idx == target_idx:
                continue

            r, p = _corr(data[tf_idx], data[target_idx])
            adjacency[tf_idx][target_idx] = r

            if abs(r) >= threshold:
                edges.append(
                    {
                        "source": gene_names[tf_idx],
                        "target": gene_names[target_idx],
                        "weight": round(r, 6),
                        "p_value": round(p, 8),
                    }
                )
                active_tfs.add(gene_names[tf_idx])
                active_targets.add(gene_names[target_idx])

    # Sort edges by absolute weight descending
    edges.sort(key=lambda e: abs(e["weight"]), reverse=True)

    logger.info(
        "GRN inference complete: %d edges, %d TFs, %d targets",
        len(edges),
        len(active_tfs),
        len(active_targets),
    )

    return {
        "edges": edges,
        "n_edges": len(edges),
        "n_tfs": len(active_tfs),
        "n_targets": len(active_targets),
        "adjacency_matrix": adjacency,
    }


def infer_grn_mutual_info(
    expression_matrix: Any,
    gene_names: list[str],
    tf_list: list[str] | None = None,
    n_bins: int = 10,
) -> dict:
    """Mutual information-based GRN inference (ARACNE-like).

    Computes pairwise mutual information between all gene pairs, then applies
    the Data Processing Inequality (DPI) to prune indirect edges. For each
    triplet of genes where all three pairs have edges, the weakest edge is
    removed, as it likely represents an indirect interaction.

    Args:
        expression_matrix: Gene expression matrix (genes x samples). Can be
            a numpy array or list of lists.
        gene_names: List of gene names corresponding to matrix rows.
        tf_list: List of transcription factor gene names. If None, all genes
            are considered potential regulators. Names that are not in
            ``gene_names`` are ignored and repeated names count once; a list
            that matches no gene raises ``ValueError`` instead of silently
            widening to all genes.
        n_bins: Number of bins for mutual information discretization.

    Returns:
        Dictionary containing:
            - edges: List of dicts with source, target, weight (MI value).
            - n_edges: Number of edges after DPI pruning.
            - mi_matrix: Full mutual information matrix as list of lists.

    Raises:
        ValueError: If gene_names length does not match matrix rows or
            ``tf_list`` matches no gene.
    """
    n_genes, n_samples = _get_matrix_shape(expression_matrix)

    if len(gene_names) != n_genes:
        raise ValueError(f"gene_names length ({len(gene_names)}) must match matrix rows ({n_genes})")

    logger.info(
        "Inferring GRN via mutual information: %d genes, %d samples, n_bins=%d",
        n_genes,
        n_samples,
        n_bins,
    )

    data = _to_list_of_lists(expression_matrix, n_genes, n_samples)
    tf_indices = _resolve_tf_indices(tf_list, gene_names)

    # Compute MI matrix
    mi_matrix: list[list[float]] = [[0.0] * n_genes for _ in range(n_genes)]

    for i in tf_indices:
        for j in range(n_genes):
            if i == j:
                continue
            if HAS_NUMPY:
                mi = _mutual_information_pure(
                    list(map(float, data[i])),
                    list(map(float, data[j])),
                    n_bins=n_bins,
                )
            else:
                mi = _mutual_information_pure(data[i], data[j], n_bins=n_bins)
            mi_matrix[i][j] = round(mi, 8)

    # Apply Data Processing Inequality (DPI) to remove indirect edges
    # For each triplet (i, j, k) where all three MI values are nonzero,
    # remove the edge with the smallest MI
    edges_to_remove: set[tuple[int, int]] = set()
    all_edge_indices = set()

    for i in tf_indices:
        for j in range(n_genes):
            if i != j and mi_matrix[i][j] > 0:
                all_edge_indices.add((i, j))

    for i, j in list(all_edge_indices):
        for k in range(n_genes):
            if k == i or k == j:
                continue
            mi_ij = mi_matrix[i][j]
            mi_ik = mi_matrix[i][k] if (i, k) in all_edge_indices else 0.0
            mi_jk = mi_matrix[j][k] if (j, k) in all_edge_indices else 0.0

            if mi_ij > 0 and mi_ik > 0 and mi_jk > 0:
                min_mi = min(mi_ij, mi_ik, mi_jk)
                if min_mi == mi_ij:
                    edges_to_remove.add((i, j))
                elif min_mi == mi_ik:
                    edges_to_remove.add((i, k))
                else:
                    edges_to_remove.add((j, k))

    # Build final edge list
    edges: list[dict] = []
    for i in tf_indices:
        for j in range(n_genes):
            if i == j:
                continue
            if mi_matrix[i][j] > 0 and (i, j) not in edges_to_remove:
                edges.append(
                    {
                        "source": gene_names[i],
                        "target": gene_names[j],
                        "weight": mi_matrix[i][j],
                    }
                )

    edges.sort(key=lambda e: e["weight"], reverse=True)

    logger.info(
        "MI-based GRN: %d edges (%d removed by DPI)",
        len(edges),
        len(edges_to_remove),
    )

    return {
        "edges": edges,
        "n_edges": len(edges),
        "mi_matrix": mi_matrix,
    }


def infer_grn_regression(
    expression_matrix: Any,
    gene_names: list[str],
    tf_list: list[str] | None = None,
    method: str = "lasso",
    alpha: float = 0.1,
) -> dict:
    """Regression-based GRN inference.

    For each target gene, fits a regression model using TF expression levels
    as predictors. Uses L1 (lasso) or L2 (ridge) regularization to identify
    the most important regulators. Non-zero coefficients indicate predicted
    regulatory relationships.

    Args:
        expression_matrix: Gene expression matrix (genes x samples). Can be
            a numpy array or list of lists.
        gene_names: List of gene names corresponding to matrix rows.
        tf_list: List of transcription factor gene names. If None, all genes
            are considered potential regulators. Names that are not in
            ``gene_names`` are ignored and repeated names count once; a list
            that matches no gene raises ``ValueError`` instead of silently
            widening to all genes.
        method: Regression method, one of 'lasso' or 'ridge'.
        alpha: Regularization strength. Higher values produce sparser networks
            for lasso.

    Returns:
        Dictionary containing:
            - edges: List of dicts with source, target, weight (coefficient).
            - coefficients: Full coefficient matrix as list of lists
              (targets x TFs).
            - r_squared_per_gene: Dict mapping gene name to R-squared value.

    Raises:
        ValueError: If gene_names length does not match matrix rows, method
            is not supported, or ``tf_list`` matches no gene.
    """
    n_genes, n_samples = _get_matrix_shape(expression_matrix)

    if len(gene_names) != n_genes:
        raise ValueError(f"gene_names length ({len(gene_names)}) must match matrix rows ({n_genes})")

    if method not in ("lasso", "ridge"):
        raise ValueError(f"Unsupported method: {method}. Use 'lasso' or 'ridge'.")

    logger.info(
        "Inferring GRN via %s regression: %d genes, %d samples, alpha=%.4f",
        method,
        n_genes,
        n_samples,
        alpha,
    )

    data = _to_list_of_lists(expression_matrix, n_genes, n_samples)
    tf_indices = _resolve_tf_indices(tf_list, gene_names)

    n_tfs = len(tf_indices)

    # Prepare TF expression matrix (samples x n_tfs)
    tf_data: list[list[float]] = []
    for s in range(n_samples):
        tf_data.append([data[tf_idx][s] for tf_idx in tf_indices])

    edges: list[dict] = []
    coefficients: list[list[float]] = []
    r_squared_per_gene: dict[str, float] = {}

    for target_idx in range(n_genes):
        target_values = data[target_idx]

        # Exclude target from predictors if it is a TF
        predictor_indices = [i for i, tf_idx in enumerate(tf_indices) if tf_idx != target_idx]

        if not predictor_indices:
            coefficients.append([0.0] * n_tfs)
            r_squared_per_gene[gene_names[target_idx]] = 0.0
            continue

        # Build predictor matrix for this target
        x_matrix = [[tf_data[s][pi] for pi in predictor_indices] for s in range(n_samples)]

        # Fit regression using coordinate descent (pure Python lasso/ridge)
        coefs, r2 = _coordinate_descent_regression(x_matrix, target_values, method=method, alpha=alpha, max_iter=1000)

        # Map coefficients back to full TF list
        full_coefs = [0.0] * n_tfs
        for ci, pi in enumerate(predictor_indices):
            full_coefs[pi] = coefs[ci]

        coefficients.append(full_coefs)
        r_squared_per_gene[gene_names[target_idx]] = round(r2, 6)

        # Collect non-zero edges
        for ci, pi in enumerate(predictor_indices):
            if abs(coefs[ci]) > 1e-8:
                edges.append(
                    {
                        "source": gene_names[tf_indices[pi]],
                        "target": gene_names[target_idx],
                        "weight": round(coefs[ci], 6),
                    }
                )

    edges.sort(key=lambda e: abs(e["weight"]), reverse=True)

    logger.info("Regression GRN: %d edges from %s regression", len(edges), method)

    return {
        "edges": edges,
        "n_edges": len(edges),
        "coefficients": coefficients,
        "r_squared_per_gene": r_squared_per_gene,
    }


def _coordinate_descent_regression(
    x: list[list[float]],
    y: list[float],
    method: str = "lasso",
    alpha: float = 0.1,
    max_iter: int = 1000,
    tol: float = 1e-6,
) -> tuple[list[float], float]:
    """Fit lasso or ridge regression via coordinate descent.

    Args:
        x: Predictor matrix (n_samples x n_features).
        y: Target values (n_samples,).
        method: 'lasso' for L1 or 'ridge' for L2 penalty.
        alpha: Regularization strength.
        max_iter: Maximum iterations.
        tol: Convergence tolerance.

    Returns:
        Tuple of (coefficients, r_squared).
    """
    n = len(y)
    if n == 0 or not x or not x[0]:
        return [], 0.0

    p = len(x[0])
    coefs = [0.0] * p

    # Precompute column-wise statistics
    col_sums_sq = [sum(x[i][j] ** 2 for i in range(n)) for j in range(p)]

    for _iteration in range(max_iter):
        max_change = 0.0

        for j in range(p):
            # Compute residual excluding feature j
            residual = [y[i] - sum(coefs[k] * x[i][k] for k in range(p) if k != j) for i in range(n)]

            # Compute raw update
            rho = sum(x[i][j] * residual[i] for i in range(n))

            denom = col_sums_sq[j]
            if denom < 1e-15:
                new_coef = 0.0
            elif method == "lasso":
                # Soft thresholding
                if rho > alpha * n:
                    new_coef = (rho - alpha * n) / denom
                elif rho < -alpha * n:
                    new_coef = (rho + alpha * n) / denom
                else:
                    new_coef = 0.0
            else:
                # Ridge regression
                new_coef = rho / (denom + alpha * n)

            change = abs(new_coef - coefs[j])
            if change > max_change:
                max_change = change
            coefs[j] = new_coef

        if max_change < tol:
            break

    # Compute R-squared
    y_mean = sum(y) / n
    ss_total = sum((yi - y_mean) ** 2 for yi in y)
    predictions = [sum(coefs[j] * x[i][j] for j in range(p)) for i in range(n)]
    ss_residual = sum((y[i] - predictions[i]) ** 2 for i in range(n))

    r2 = 1.0 - (ss_residual / ss_total) if ss_total > 1e-15 else 0.0
    r2 = max(0.0, r2)

    return coefs, r2


def score_regulators(
    grn: dict,
    gene_set: list[str],
) -> list[dict]:
    """Score transcription factors for enrichment of their targets in a gene set.

    For each TF in the GRN, computes how many of its targets overlap with the
    given gene set, then tests the overlap with the exact one-sided (enrichment)
    Fisher test, i.e. the hypergeometric upper tail.

    The population of that test is the set of distinct targets in the GRN, and its successes are the gene-set
    members that are in that population (matched case-insensitively). Gene-set members that the GRN never targets
    cannot be drawn by any TF, so they are left out; counting them (and duplicate entries of ``gene_set``) made
    every p-value too large.

    Args:
        grn: GRN dictionary as returned by infer_grn_* functions. Must contain
            an 'edges' key with list of dicts having 'source' and 'target'.
        gene_set: List of genes of interest (e.g., differentially expressed
            genes).

    Returns:
        List of dicts sorted by ``enrichment_p``, each containing:
            - tf: Transcription factor name.
            - n_targets: Total number of targets for this TF.
            - n_in_set: Number of targets found in gene_set.
            - enrichment_p: Exact one-sided Fisher p-value of observing at least ``n_in_set`` gene-set targets
              among the ``n_targets``, in ``[0, 1]``. It is not rounded, so a strongly enriched TF keeps its tiny
              p-value (rounding to 8 decimals reported 0.0 for anything below 5e-9).
            - regulation_score: Weighted regulation score.

    Raises:
        ValueError: If grn does not contain 'edges' key.
    """
    if "edges" not in grn:
        raise ValueError("GRN dictionary must contain 'edges' key")

    edges = grn["edges"]
    gene_set_lower = {g.lower() for g in gene_set}

    # Build TF -> targets mapping and per-TF weighted scores in one pass
    tf_targets: dict[str, set[str]] = defaultdict(set)
    tf_scores: dict[str, float] = defaultdict(float)
    all_targets: set[str] = set()

    for edge in edges:
        source = edge["source"]
        target = edge["target"]
        tf_targets[source].add(target)
        all_targets.add(target)
        if target.lower() in gene_set_lower:
            tf_scores[source] += abs(edge.get("weight", 1.0))

    total_targets = len(all_targets)
    # Successes in the population: the gene-set members that some TF targets
    set_size = sum(1 for target in all_targets if target.lower() in gene_set_lower)

    results: list[dict] = []

    for tf, targets in tf_targets.items():
        n_targets = len(targets)
        targets_in_set = {t for t in targets if t.lower() in gene_set_lower}
        n_in_set = len(targets_in_set)

        # Exact hypergeometric upper tail (one-sided Fisher test)
        enrichment_p = _hypergeometric_pvalue(n_in_set, n_targets, set_size, total_targets)

        # Regulation score: weighted by edge weights toward the gene set
        regulation_score = tf_scores.get(tf, 0.0)

        results.append(
            {
                "tf": tf,
                "n_targets": n_targets,
                "n_in_set": n_in_set,
                "enrichment_p": enrichment_p,
                "regulation_score": round(regulation_score, 6),
            }
        )

    results.sort(key=lambda r: r["enrichment_p"])

    logger.info(
        "Scored %d regulators, top TF: %s (p=%.2e)",
        len(results),
        results[0]["tf"] if results else "none",
        results[0]["enrichment_p"] if results else 1.0,
    )

    return results


def _log_binomial(n: int, k: int) -> float:
    """``log(C(n, k))`` for ``0 <= k <= n``, through ``lgamma`` so that huge binomials never overflow."""
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _hypergeometric_pvalue(k: int, n: int, K: int, N: int) -> float:
    """Exact upper tail ``P(X >= k)`` of ``X ~ Hypergeometric(N, K, n)``.

    This is the p-value of the one-sided (``alternative="greater"``) Fisher exact test of the 2x2 table
    ``[[k, n - k], [K - k, N - K - n + k]]`` and equals ``scipy.stats.hypergeom.sf(k - 1, N, K, n)``. Every
    probability mass ``C(K, x) C(N - K, n - x) / C(N, n)`` of the tail is formed in log space (``lgamma``) and the
    tail is summed after factoring out its largest term, so no binomial overflows and a tiny tail keeps its relative
    precision (about 1e-10 for a population of 20 000). SciPy is not needed and the result does not depend on whether
    it is installed. The result only underflows to 0.0 when the true p-value is below about 5e-324.

    An earlier version fell back to a normal approximation without SciPy. That approximation is poor exactly where
    enrichment is tested: small counts and a small gene set in a large population.

    Args:
        k: Number of successes observed.
        n: Number of draws (TF targets).
        K: Total successes in population (gene-set members among the targets).
        N: Population size (total targets).

    Returns:
        P-value for enrichment, in [0, 1]. ``k <= 0`` gives 1.0, and a ``k`` larger than ``min(n, K)`` (more
        successes than there are draws or successes) is an impossible event with probability 0.0.

    Raises:
        ValueError: If the arguments do not describe a draw from a population (``K`` or ``n`` negative or larger
            than ``N``). A call used to return 1.0 for what is an invalid population (for example a gene set larger
            than the target universe, where SciPy returns NaN that was then clamped to 1.0).
    """
    if N < 0 or not 0 <= K <= N or not 0 <= n <= N:
        raise ValueError(
            f"inconsistent hypergeometric parameters: population N={N}, successes K={K}, draws n={n} "
            "(need 0 <= K <= N and 0 <= n <= N)"
        )
    lowest_possible = max(0, n + K - N)
    highest_possible = min(K, n)
    if k <= lowest_possible:
        return 1.0  # every possible outcome is at least k
    if k > highest_possible:
        return 0.0  # more successes than a draw can hold

    log_total = _log_binomial(N, n)
    log_terms = [_log_binomial(K, x) + _log_binomial(N - K, n - x) - log_total for x in range(k, highest_possible + 1)]
    peak = max(log_terms)
    log_tail = peak + math.log(math.fsum(math.exp(term - peak) for term in log_terms))
    return min(1.0, math.exp(log_tail))


def compute_network_motifs(
    edges: list[dict],
    motif_size: int = 3,
    random_state: int | random.Random | None = None,
) -> dict:
    """Count network motifs in a directed regulatory network.

    Identifies and counts recurring small subgraph patterns (motifs) such as
    feed-forward loops (FFLs), mutual regulation pairs, and regulatory
    cascades. Compares observed counts against a randomized ensemble to
    compute z-scores. Motifs span two or three distinct nodes, so
    autoregulation (self-loops) and repeated edges are not counted; a z-score
    is ``+inf`` / ``-inf`` when the randomized counts never vary and the
    observed count differs from them.

    Args:
        edges: List of edge dicts with 'source' and 'target' keys.
        motif_size: Size of motifs to search for (2 or 3).
        random_state: Seed (or a ``random.Random`` generator) for the randomized ensemble behind the z-scores.
            The same seed gives the same z-scores; the global ``random`` state is neither seeded nor consumed.
            ``None`` draws fresh entropy.

    Returns:
        Dictionary containing:
            - motif_counts: Dict mapping motif type to observed count.
            - z_scores: Dict mapping motif type to z-score vs random.
            - significant_motifs: List of motif types with |z| > 2.

    Raises:
        ValueError: If motif_size is not 2 or 3.
    """
    if motif_size not in (2, 3):
        raise ValueError(f"motif_size must be 2 or 3, got {motif_size}")

    logger.info("Computing network motifs of size %d from %d edges", motif_size, len(edges))

    nodes: set[str] = set()
    for edge in edges:
        nodes.add(edge["source"])
        nodes.add(edge["target"])
    node_list = sorted(nodes)

    # The observed network and every randomized one are counted by the same function, so a motif means the same
    # thing on both sides of the z-score.
    motif_counts = _count_motifs_fast(edges, motif_size, node_list)

    # Compute z-scores via edge randomization
    n_random = 100
    random_counts: dict[str, list[int]] = {k: [] for k in motif_counts}

    edge_list = [(e["source"], e["target"]) for e in edges]
    rng = random_state if isinstance(random_state, random.Random) else random.Random(random_state)

    for _ in range(n_random):
        # Randomize by shuffling target nodes
        shuffled_targets = [t for _, t in edge_list]
        rng.shuffle(shuffled_targets)
        random_edges = [{"source": s, "target": t} for (s, _), t in zip(edge_list, shuffled_targets)]

        random_result = _count_motifs_fast(random_edges, motif_size, node_list)
        for k in motif_counts:
            random_counts[k].append(random_result.get(k, 0))

    z_scores = _motif_z_scores(motif_counts, random_counts)

    significant_motifs = [k for k, z in z_scores.items() if abs(z) > 2.0]

    logger.info(
        "Motif analysis complete: %d types, %d significant",
        len(motif_counts),
        len(significant_motifs),
    )

    return {
        "motif_counts": motif_counts,
        "z_scores": z_scores,
        "significant_motifs": significant_motifs,
    }


def _motif_z_scores(motif_counts: dict[str, int], random_counts: dict[str, list[int]]) -> dict[str, float]:
    """z-score of each observed motif count against its randomized counts.

    When the randomized counts do not vary at all the z-score is undefined: it is 0.0 if the observed count equals
    them, and ``+inf`` / ``-inf`` if it is above / below them (infinitely many standard deviations away). It used to
    be 0.0 whatever the observed count was, which hid a motif that never occurs in the random ensemble.
    """
    z_scores: dict[str, float] = {}
    for motif, observed in motif_counts.items():
        values = random_counts[motif]
        mean_val = sum(values) / len(values) if values else 0.0
        std_val = math.sqrt(sum((v - mean_val) ** 2 for v in values) / max(len(values) - 1, 1))
        if std_val > 1e-15:
            z_scores[motif] = round((observed - mean_val) / std_val, 4)
        elif observed == mean_val:
            z_scores[motif] = 0.0
        else:
            z_scores[motif] = math.copysign(math.inf, observed - mean_val)
    return z_scores


def _count_motifs_fast(edges: list[dict], motif_size: int, node_list: list[str]) -> dict[str, int]:
    """Count 2- or 3-node regulatory motifs of a directed edge list.

    Autoregulation (a self-loop) belongs to no motif of two or three distinct nodes and is ignored: counting it made
    ``A->A, A->B`` contain a feed-forward loop on two nodes, and a network with every TF autoregulated report 7
    feed-forward loops where there is 1. Repeated edges count once, and co-regulation counts the pairs of
    *distinct* regulators of each target.

    Args:
        edges: List of edge dicts with 'source' and 'target'.
        motif_size: Size of motifs to count (2 or 3).
        node_list: Sorted list of node names.

    Returns:
        Dict mapping motif type to count.
    """
    adj: dict[str, set[str]] = defaultdict(set)
    for edge in edges:
        if edge["source"] != edge["target"]:
            adj[edge["source"]].add(edge["target"])

    if motif_size == 2:
        mutual = 0
        for a in node_list:
            for b in adj.get(a, set()):
                if a < b and a in adj.get(b, set()):
                    mutual += 1
        total_directed = sum(len(t) for t in adj.values())
        return {
            "mutual_regulation": mutual,
            "single_regulation": total_directed - 2 * mutual,
        }

    # Size 3: a, b, c are distinct because there are no self-loops and c == a is skipped
    ffl = 0  # A->B, B->C, A->C
    cascade = 0  # A->B, B->C, no A->C
    for a in node_list:
        for b in adj.get(a, set()):
            for c in adj.get(b, set()):
                if c == a:
                    continue
                if c in adj[a]:
                    ffl += 1
                else:
                    cascade += 1

    regulators: dict[str, set[str]] = defaultdict(set)
    for source, targets in adj.items():
        for target in targets:
            regulators[target].add(source)
    co_reg = sum(len(r) * (len(r) - 1) // 2 for r in regulators.values())

    return {
        "feed_forward_loop": ffl,
        "cascade": cascade,
        "co_regulation": co_reg,
    }


def validate_grn(
    predicted_edges: list[dict],
    known_edges: list[dict],
) -> dict:
    """Validate a predicted GRN against known regulatory interactions.

    Computes standard classification metrics by treating each predicted edge
    as a positive prediction and each known edge as a true positive.

    Args:
        predicted_edges: List of predicted edge dicts with 'source' and
            'target' keys.
        known_edges: List of known/validated edge dicts with 'source' and
            'target' keys.

    Returns:
        Dictionary containing:
            - precision: Fraction of predicted edges that are known.
            - recall: Fraction of known edges that are predicted.
            - f1: Harmonic mean of precision and recall.
            - auroc: Area under ROC curve (based on edge weights if
              available, otherwise binary).
            - auprc: Area under precision-recall curve.

    Raises:
        ValueError: If either edge list is empty.
    """
    if not predicted_edges:
        raise ValueError("predicted_edges must not be empty")
    if not known_edges:
        raise ValueError("known_edges must not be empty")

    logger.info(
        "Validating GRN: %d predicted vs %d known edges",
        len(predicted_edges),
        len(known_edges),
    )

    predicted_set = {(e["source"], e["target"]) for e in predicted_edges}
    known_set = {(e["source"], e["target"]) for e in known_edges}

    true_positives = predicted_set & known_set
    tp = len(true_positives)
    fp = len(predicted_set - known_set)
    fn = len(known_set - predicted_set)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    # Compute AUROC and AUPRC using edge weights as scores
    all_nodes = set()
    for e in predicted_edges + known_edges:
        all_nodes.add(e["source"])
        all_nodes.add(e["target"])

    # Build score dict from predicted edges
    pred_scores: dict[tuple[str, str], float] = {}
    for e in predicted_edges:
        key = (e["source"], e["target"])
        pred_scores[key] = abs(e.get("weight", 1.0))

    # All possible edges (for ROC/PRC computation)
    all_possible = set()
    node_list = sorted(all_nodes)
    for s in node_list:
        for t in node_list:
            if s != t:
                all_possible.add((s, t))

    # Score all possible edges
    scored: list[tuple[float, bool]] = []
    for edge in all_possible:
        score = pred_scores.get(edge, 0.0)
        is_known = edge in known_set
        scored.append((score, is_known))

    # Sort by score descending
    scored.sort(key=lambda x: x[0], reverse=True)

    # Compute AUROC and AUPRC via trapezoidal rule
    auroc = _compute_auroc(scored)
    auprc = _compute_auprc(scored)

    result = {
        "precision": round(precision, 6),
        "recall": round(recall, 6),
        "f1": round(f1, 6),
        "auroc": round(auroc, 6),
        "auprc": round(auprc, 6),
    }

    logger.info(
        "GRN validation: precision=%.4f, recall=%.4f, F1=%.4f, AUROC=%.4f",
        precision,
        recall,
        f1,
        auroc,
    )

    return result


def _score_groups(scored: list[tuple[float, bool]]) -> list[tuple[int, int]]:
    """``(positives, negatives)`` for each distinct score, from the highest score to the lowest.

    Items with the same score cannot be ranked against each other, so a curve may only move from one distinct score
    to the next. Walking the items one at a time instead let the order of tied items decide the area: the candidate
    edges come out of a ``set`` of string pairs, whose iteration order changes with ``PYTHONHASHSEED``, so an
    unweighted prediction (every predicted edge scoring 1.0, every other edge 0.0) got a different AUROC and AUPRC
    from one run to the next.
    """
    groups: list[tuple[int, int]] = []
    previous_score: float | None = None
    for score, is_positive in sorted(scored, key=lambda item: item[0], reverse=True):
        if previous_score is None or score != previous_score:
            groups.append((0, 0))
            previous_score = score
        positives, negatives = groups[-1]
        groups[-1] = (positives + int(is_positive), negatives + int(not is_positive))
    return groups


def _compute_auroc(scored: list[tuple[float, bool]]) -> float:
    """Compute AUROC from scored predictions, ties handled exactly.

    The ROC curve moves from one distinct score to the next along a straight segment, so tied positives and
    negatives count half, as in the Mann-Whitney statistic.

    Args:
        scored: List of (score, is_positive) tuples.

    Returns:
        Area under ROC curve.
    """
    n_pos = sum(1 for _, is_pos in scored if is_pos)
    n_neg = len(scored) - n_pos

    if n_pos == 0 or n_neg == 0:
        return 0.5

    tpr_prev = 0.0
    fpr_prev = 0.0
    tp = 0
    fp = 0
    auc = 0.0

    for group_pos, group_neg in _score_groups(scored):
        tp += group_pos
        fp += group_neg

        tpr = tp / n_pos
        fpr = fp / n_neg

        # Trapezoidal area
        auc += (fpr - fpr_prev) * (tpr + tpr_prev) / 2.0
        tpr_prev = tpr
        fpr_prev = fpr

    return max(0.0, min(1.0, auc))


def _compute_auprc(scored: list[tuple[float, bool]]) -> float:
    """Compute AUPRC (average precision) from scored predictions, ties handled exactly.

    Precision and recall are evaluated once per distinct score, after all items with that score are counted, so the
    result is the step-wise average precision ``sum((R_n - R_{n-1}) * P_n)`` of scikit-learn's
    ``average_precision_score`` and does not depend on the order of tied items.

    Args:
        scored: List of (score, is_positive) tuples.

    Returns:
        Area under precision-recall curve.
    """
    n_pos = sum(1 for _, is_pos in scored if is_pos)
    if n_pos == 0:
        return 0.0

    tp = 0
    fp = 0
    auprc = 0.0
    prev_recall = 0.0

    for group_pos, group_neg in _score_groups(scored):
        tp += group_pos
        fp += group_neg

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / n_pos

        auprc += (recall - prev_recall) * precision
        prev_recall = recall

    return max(0.0, min(1.0, auprc))
