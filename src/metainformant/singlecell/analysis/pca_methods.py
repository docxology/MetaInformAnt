"""Linear dimensionality reduction methods for single-cell data.

This module provides PCA, ICA, and Factor Analysis methods for reducing
the dimensionality of single-cell expression data, along with highly
variable gene (HVG) selection.
"""

from __future__ import annotations

import warnings
from typing import Any, Dict

import numpy as np
import pandas as pd

from metainformant.core.data import validation
from metainformant.core.utils import logging
from metainformant.singlecell.data.preprocessing import SingleCellData

# Try to import optional dependencies
try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    PCA = None
    StandardScaler = None

logger = logging.get_logger(__name__)


def _add_embedding_to_obs(result: SingleCellData, coords: pd.DataFrame) -> None:
    """Attach embedding coordinates to ``result.obs``, creating obs if needed."""
    if result.obs is None:
        result.obs = coords
    else:
        for col in coords.columns:
            result.obs[col] = coords[col]


def _highly_variable_mask(data: SingleCellData) -> np.ndarray:
    """Return the boolean gene mask stored in ``data.var["highly_variable"]``.

    This is the annotation written by :func:`select_hvgs` and by
    ``identify_highly_variable_genes``.

    Raises:
        ValueError: If the annotation is missing, is not boolean (or 0/1) without
            missing values, or flags no gene.
    """
    var = getattr(data, "var", None)
    if var is None or "highly_variable" not in var.columns:
        raise ValueError(
            "use_hvgs=True needs a 'highly_variable' column in data.var: run select_hvgs() "
            "(or identify_highly_variable_genes()) first, or pass use_hvgs=False"
        )

    column = var["highly_variable"]
    values = column.to_numpy()
    if column.isna().any():
        raise ValueError("data.var['highly_variable'] must not contain missing values")
    if values.dtype == object and all(isinstance(value, (bool, np.bool_)) for value in values):
        values = values.astype(bool)
    if values.dtype != bool:
        if not (pd.api.types.is_numeric_dtype(values.dtype) and np.isin(values, (0, 1)).all()):
            raise ValueError("data.var['highly_variable'] must be a boolean column (True/False or 1/0)")
        values = values.astype(bool)

    mask = np.asarray(values, dtype=bool)
    if not mask.any():
        raise ValueError("use_hvgs=True but no genes are flagged in data.var['highly_variable']")
    return mask


def pca_reduction(
    data: SingleCellData,
    n_components: int = 50,
    random_state: int | None = None,
    scale_data: bool = True,
    *,
    use_hvgs: bool = False,
) -> SingleCellData:
    """Perform PCA dimensionality reduction on single-cell data.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of principal components to compute. Clamped to
            ``min(n_obs - 1, n_genes_used)``.
        random_state: Random seed for reproducibility
        scale_data: Whether to scale data before PCA
        use_hvgs: Restrict the analysis to the genes flagged in
            ``data.var["highly_variable"]`` (written by :func:`select_hvgs`).
            Scaling, the variance ratios and the component clamp then refer to
            those genes only, and ``varm["PCs"]`` keeps one row per gene with
            exact zeros for the genes that were left out.

    Returns:
        SingleCellData with PCA coordinates added

    Raises:
        ValidationError: If data is not SingleCellData
        ValueError: If n_components is invalid, or ``use_hvgs`` is set but ``data.var["highly_variable"]`` is missing,
            not boolean, or flags no gene
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn is required for PCA. " "Install with: uv pip install scikit-learn")

    validation.validate_type(data, SingleCellData, "data")

    gene_mask = _highly_variable_mask(data) if use_hvgs else None
    n_genes_used = int(gene_mask.sum()) if gene_mask is not None else data.n_vars

    # Auto-clamp n_components to the maximum possible value
    max_components = min(data.n_obs - 1, n_genes_used)
    if n_components > max_components:
        logger.info(f"Clamping n_components from {n_components} to {max_components} (max possible)")
        n_components = max_components
    if n_components < 1:
        raise ValueError("n_components must be at least 1")

    logger.info(f"Performing PCA with {n_components} components on {n_genes_used} of {data.n_vars} genes")

    # Create copy
    result = data.copy()

    # Get expression matrix (restricted to the selected genes before densifying)
    expression = data.X if gene_mask is None else data.X[:, gene_mask]
    X = expression.toarray() if hasattr(expression, "toarray") else expression

    # Scale data if requested
    if scale_data:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
    else:
        X_scaled = X

    # Perform PCA
    pca = PCA(n_components=n_components, random_state=random_state)
    X_pca = pca.fit_transform(X_scaled)

    # Store PCA results in obsm (standard location)
    result.obsm = result.obsm if hasattr(result, "obsm") and result.obsm is not None else {}
    result.obsm["X_pca"] = X_pca

    # Also add PCA coordinates to obs for convenience
    pca_coords = pd.DataFrame(
        X_pca,
        index=result.obs.index if result.obs is not None else None,
        columns=[f"PC{i+1}" for i in range(n_components)],
    )
    _add_embedding_to_obs(result, pca_coords)

    # Gene loadings: one row per gene (the AnnData varm convention). Genes outside the analysis
    # carry exactly zero loading.
    if gene_mask is None:
        loadings = pca.components_.T
    else:
        loadings = np.zeros((data.n_vars, n_components), dtype=pca.components_.dtype)
        loadings[gene_mask] = pca.components_.T
    result.varm = result.varm if hasattr(result, "varm") and result.varm is not None else {}
    result.varm["PCs"] = loadings

    # Store PCA metadata
    explained_variance = pca.explained_variance_ratio_
    cumulative_variance = np.cumsum(explained_variance)

    result.uns["pca"] = {
        "n_components": n_components,
        "explained_variance_ratio": explained_variance.tolist(),
        "variance_ratio": explained_variance.tolist(),  # Alias
        "variance": explained_variance.tolist(),  # Alias
        "cumulative_explained_variance": cumulative_variance.tolist(),
        "cumulative_variance": cumulative_variance.tolist(),  # Alias
        "singular_values": pca.singular_values_.tolist(),
        "components_shape": loadings.T.shape,
        "random_state": random_state,
        "scaled": scale_data,
        "use_hvgs": use_hvgs,
        "n_genes_used": n_genes_used,
    }

    logger.info(f"PCA completed: {n_components} components explain {cumulative_variance[-1]:.1%} of variance")
    return result


def ica_reduction(
    data: SingleCellData, n_components: int = 10, random_state: int | None = None, max_iter: int = 1000
) -> SingleCellData:
    """Perform Independent Component Analysis (ICA) on single-cell data.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of independent components
        random_state: Random seed for reproducibility
        max_iter: Maximum number of iterations

    Returns:
        SingleCellData with ICA coordinates added

    Raises:
        ValidationError: If data is not SingleCellData or n_components is outside 2..min(n_obs, n_vars)
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn is required for ICA. Install with: uv pip install scikit-learn")
    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=2, max_val=min(data.n_obs, data.n_vars), name="n_components")

    logger.info(f"Performing ICA with {n_components} components")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Scale data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Perform ICA
    from sklearn.decomposition import FastICA

    ica = FastICA(n_components=n_components, random_state=random_state, max_iter=max_iter, tol=1e-4)

    X_ica = ica.fit_transform(X_scaled)

    # Store ICA results
    ica_coords = pd.DataFrame(
        X_ica,
        index=result.obs.index if result.obs is not None else None,
        columns=[f"IC{i+1}" for i in range(n_components)],
    )

    _add_embedding_to_obs(result, ica_coords)

    # Store ICA metadata and components
    result.uns["ica"] = {
        "n_components": n_components,
        "random_state": random_state,
        "max_iter": max_iter,
        "n_iter": ica.n_iter_,
    }

    # Store mixing matrix for gene loadings
    result.varm = result.varm if hasattr(result, "varm") and result.varm is not None else {}
    result.varm["ICs"] = ica.mixing_.T

    logger.info(f"ICA completed: {n_components} components extracted")
    return result


def factor_analysis_reduction(
    data: SingleCellData, n_components: int = 10, random_state: int | None = None
) -> SingleCellData:
    """Perform Factor Analysis on single-cell data.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of factors
        random_state: Random seed for reproducibility

    Returns:
        SingleCellData with factor analysis coordinates added

    Raises:
        ValidationError: If data is not SingleCellData or n_components is outside 2..min(n_obs, n_vars)
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn is required for Factor Analysis. Install with: uv pip install scikit-learn")
    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=2, max_val=min(data.n_obs, data.n_vars), name="n_components")

    logger.info(f"Performing Factor Analysis with {n_components} components")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Scale data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Perform Factor Analysis
    from sklearn.decomposition import FactorAnalysis

    fa = FactorAnalysis(n_components=n_components, random_state=random_state, max_iter=1000)

    X_fa = fa.fit_transform(X_scaled)

    # Store FA results
    fa_coords = pd.DataFrame(
        X_fa,
        index=result.obs.index if result.obs is not None else None,
        columns=[f"FA{i+1}" for i in range(n_components)],
    )

    _add_embedding_to_obs(result, fa_coords)

    # Store FA metadata and components
    result.uns["factor_analysis"] = {
        "n_components": n_components,
        "random_state": random_state,
        "log_likelihood": float(fa.loglike_[-1]) if hasattr(fa, "loglike_") else None,
        "noise_variance": fa.noise_variance_.tolist(),
    }

    # Store factor loadings
    result.varm = result.varm if hasattr(result, "varm") and result.varm is not None else {}
    result.varm["FA_loadings"] = fa.components_.T

    logger.info(f"Factor Analysis completed: {n_components} factors extracted")
    return result


def run_pca(
    data: SingleCellData | None = None,
    n_components: int = 50,
    random_state: int | None = None,
    scale_data: bool = True,
    *,
    expression_matrix: np.ndarray | None = None,
) -> SingleCellData | Dict[str, Any]:
    """Run PCA dimensionality reduction.

    Args:
        data: SingleCellData object (preferred)
        n_components: Number of components
        random_state: Random seed
        scale_data: Whether to scale data
        expression_matrix: Raw expression matrix (alternative input, returns dict)

    Returns:
        SingleCellData with PCA results if data provided, or dict if expression_matrix provided
    """
    # Handle expression_matrix input (returns dict format for test compatibility)
    if expression_matrix is not None:
        from sklearn.decomposition import PCA as SklearnPCA

        X = expression_matrix
        n_comp = min(n_components, X.shape[0], X.shape[1])

        if scale_data:
            X_scaled = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-10)
        else:
            X_scaled = X

        pca = SklearnPCA(n_components=n_comp, random_state=random_state)
        embedding = pca.fit_transform(X_scaled)

        return {
            "embedding": embedding,
            "components": pca.components_.T,  # (n_features, n_components)
            "explained_variance": pca.explained_variance_ratio_,
        }

    # Handle SingleCellData input (original behavior)
    if data is None:
        raise ValueError("Either data or expression_matrix must be provided")
    return pca_reduction(data, n_components, random_state, scale_data)


def compute_pca(
    data: SingleCellData,
    n_components: int = 50,
    random_state: int | None = None,
    scale_data: bool = True,
    *,
    use_hvgs: bool = False,
) -> SingleCellData:
    """Compute PCA dimensionality reduction (alias for pca_reduction).

    Args:
        data: SingleCellData object
        n_components: Number of components
        random_state: Random seed
        scale_data: Whether to scale data
        use_hvgs: Restrict PCA to the genes flagged in ``data.var["highly_variable"]``
            (run :func:`select_hvgs` first)

    Returns:
        SingleCellData with PCA results

    Raises:
        ValueError: If ``use_hvgs`` is set but ``data.var`` has no usable
            ``highly_variable`` annotation
    """
    return pca_reduction(data, n_components, random_state, scale_data, use_hvgs=use_hvgs)


# Fraction of the genes that every local fit of the seurat_v3 mean-variance trend uses (scanpy's default).
_SEURAT_V3_SPAN = 0.3


def _local_polynomial_fit(x: np.ndarray, y: np.ndarray, span: float, degree: int) -> np.ndarray:
    """Tricube-weighted local polynomial regression of ``y`` on ``x``, evaluated at every ``x``.

    This is loess without the kd-tree interpolation: the ``q = int(len(x) * span + 1e-5)`` points nearest to
    each ``x[i]`` are fitted with a polynomial of the given degree, weighted by ``(1 - (d / h)**3)**3`` where
    ``d`` is the distance to ``x[i]`` and ``h`` the distance of the farthest of the ``q`` points, and the
    polynomial's value at ``x[i]`` is the fit. Because ``x`` is one-dimensional the nearest points form a
    window of the sorted values, found here by a vectorised binary search. Cost: O(len(x) * q).

    Raises:
        ValueError: If ``span`` leaves fewer than ``degree + 3`` points per fit, or if ``q`` or more of the
            points share one ``x`` so that a neighbourhood has no width (the same condition that makes
            R's and scikit-misc's loess refuse to fit).
    """
    n = len(x)
    q = int(n * span + 1e-5)
    if q < degree + 3:
        raise ValueError(
            f"a span of {span} over {n} points leaves {q} per local fit, fewer than the {degree + 3} that a "
            f"degree-{degree} fit needs: at least {int(np.ceil((degree + 3) / span))} points are required"
        )
    order = np.argsort(x, kind="stable")
    xs, ys = x[order], y[order]

    # First index of the window of q sorted points that is nearest to each point: moving the window one to
    # the right pays exactly when its left end is farther from the point than the new right end.
    position = np.arange(n)
    low = np.maximum(0, position - q + 1)
    high = np.minimum(position, n - q)
    while True:
        active = low < high
        if not active.any():
            break
        mid = (low + high) // 2
        move_right = (xs - xs[mid]) > (xs[np.minimum(mid + q, n - 1)] - xs)
        low = np.where(active & move_right, mid + 1, low)
        high = np.where(active & ~move_right, mid, high)

    offsets = np.arange(q)
    block = max(1, 2_000_000 // q)  # targets per block: bounds the (block, q) temporaries
    fitted_sorted = np.empty(n, dtype=np.float64)
    for begin in range(0, n, block):
        selection = slice(begin, min(n, begin + block))
        first = low[selection]
        window = first[:, None] + offsets[None, :]
        centre = xs[selection, None]
        width = np.maximum(centre - xs[first][:, None], xs[first + q - 1][:, None] - centre)
        if np.any(width <= 0):
            raise ValueError(
                f"{q} or more values share exactly the same x, so a local fit has a neighbourhood of zero "
                "width; remove the repeated values (for genes: filter out rarely detected ones) or use a "
                "larger span"
            )
        u = (xs[window] - centre) / width
        weight = np.clip(1.0 - np.abs(u) ** 3, 0.0, None) ** 3

        powers = [np.ones_like(u)]
        for _ in range(2 * degree):
            powers.append(powers[-1] * u)
        moments = [(weight * power).sum(axis=1) for power in powers]
        responses = [(weight * powers[k] * ys[window]).sum(axis=1) for k in range(degree + 1)]
        normal = np.empty((len(first), degree + 1, degree + 1))
        for row in range(degree + 1):
            for column in range(degree + 1):
                normal[:, row, column] = moments[row + column]
        coefficients = np.linalg.pinv(normal, rcond=1e-10) @ np.stack(responses, axis=1)[..., None]
        fitted_sorted[selection] = coefficients[:, 0, 0]  # the basis is centred on the target: intercept = fit

    fitted = np.empty(n, dtype=np.float64)
    fitted[order] = fitted_sorted
    return fitted


def _seurat_v3_normalized_variance(counts: np.ndarray, span: float = _SEURAT_V3_SPAN) -> np.ndarray:
    """Variance of the standardised, clipped counts of every gene (the ``seurat_v3`` ranking statistic).

    Follows the definition scanpy uses for ``flavor="seurat_v3"`` with a single batch: the unbiased variance
    of each gene is regressed on its mean (both log10) by a degree-2 local regression over the genes with
    non-zero variance; the fitted value gives the expected standard deviation ``s`` of a gene of that mean;
    counts are clipped at ``s * sqrt(n_cells) + mean``; and the statistic is the variance of the clipped
    counts about the mean, in units of ``s``:
    ``(n * mean**2 + sum(clipped**2) - 2 * mean * sum(clipped)) / ((n - 1) * s**2)``.
    Genes with no variance get 0.

    Raises:
        ValueError: If ``counts`` is not a 2-D array of at least two cells, holds negative or non-finite
            values, or cannot support the local regression (see :func:`_local_polynomial_fit`).
    """
    if counts.ndim != 2 or counts.shape[0] < 2:
        raise ValueError(f"seurat_v3 needs a cells x genes matrix with at least two cells, got shape {counts.shape}")
    if not np.all(np.isfinite(counts)) or np.any(counts < 0):
        raise ValueError("seurat_v3 expects raw counts: the matrix holds negative or non-finite values")
    if not np.allclose(counts, np.round(counts)):
        warnings.warn(
            "flavor='seurat_v3' expects raw count data, but non-integer values were found", UserWarning, stacklevel=3
        )

    n_cells, n_genes = counts.shape
    mean = counts.mean(axis=0)
    variance = counts.var(axis=0, ddof=1)
    varying = variance > 0

    fitted_log_variance = np.zeros(n_genes, dtype=np.float64)
    fitted_log_variance[varying] = _local_polynomial_fit(
        np.log10(mean[varying]), np.log10(variance[varying]), span, degree=2
    )
    expected_std = np.sqrt(10.0**fitted_log_variance)

    clipped = np.minimum(counts, expected_std * np.sqrt(n_cells) + mean)
    clipped_sum = clipped.sum(axis=0)
    clipped_squares = np.square(clipped).sum(axis=0)
    return np.asarray(
        (n_cells * np.square(mean) + clipped_squares - 2.0 * clipped_sum * mean)
        / ((n_cells - 1) * np.square(expected_std)),
        dtype=np.float64,
    )


def select_hvgs(
    data: SingleCellData,
    n_top_genes: int = 2000,
    flavor: str | None = None,
    *,
    method: str | None = None,
    inplace: bool = False,
) -> SingleCellData:
    """Select highly variable genes from single-cell data.

    Like the other functions of this module, this returns a modified copy and leaves ``data``
    untouched unless ``inplace=True``. The genes are only flagged (``var["highly_variable"]``); the
    expression matrix is not subset.

    Args:
        data: SingleCellData object
        n_top_genes: Number of top variable genes to select
        flavor: Method for HVG selection: ``'seurat'`` (the default), ``'cell_ranger'``, ``'seurat_v3'`` or
            ``'variance'``. ``'seurat_v3'`` is the variance-stabilising selection of Stuart et al. (2019)
            as scanpy defines it, and expects raw counts in ``X``: the log10 variance of the genes is
            regressed on their log10 mean (local regression of degree 2, span 0.3), the counts are
            standardised with the fitted standard deviation and clipped, and the genes are ranked by the
            variance of the standardised counts, which is stored in ``var["variances_norm"]`` (for the
            other flavors that column is a z-score of the variance). The local regression is evaluated
            exactly at every gene and takes O(genes squared) time, about two seconds for 20,000 genes;
            scanpy gets it from scikit-misc, which interpolates between kd-tree vertices, so its values
            differ slightly from these. That difference has not been measured here (scikit-misc is not
            installed), and genes right at the ``n_top_genes`` cut-off can rank differently.
        method: Alias for flavor
        inplace: Write the annotations into ``data`` and return ``data`` itself instead of a copy.

    Returns:
        SingleCellData with highly variable genes marked (a copy of ``data`` unless ``inplace=True``)

    Raises:
        ValueError: If ``data`` has no expression matrix ``X``, if ``flavor`` is not one of the names above
            or, for ``'seurat_v3'``, if ``X`` has fewer than two cells, holds negative or non-finite
            values, has too few genes with non-zero variance for the local regression (at least 17), or has
            ``int(0.3 * n)`` or more of those ``n`` genes with exactly the same mean. This function raises the
            built-in ``ValueError`` (not the package's ``ValidationError``).
    """
    # Handle parameter alias
    if method is not None and flavor is None:
        flavor = method
    if flavor is None:
        flavor = "seurat"

    if not hasattr(data, "X") or data.X is None:
        raise ValueError("Data must contain expression matrix X")

    X = data.X
    if hasattr(X, "toarray"):  # sparse matrix
        X = X.toarray()

    n_cells, n_genes = X.shape
    seurat_v3_variances: np.ndarray | None = None

    if flavor == "seurat":
        # Seurat-style HVG selection
        # Calculate mean and variance for each gene
        gene_means = np.mean(X, axis=0)
        gene_vars = np.var(X, axis=0)

        # Fit curve: variance = a * mean^b
        # Use only genes with mean > 0 for fitting
        nonzero_mask = gene_means > 0
        if np.sum(nonzero_mask) < 10:
            # Fallback: just select by variance
            top_indices = np.argsort(gene_vars)[::-1][:n_top_genes]
        else:
            log_means = np.log10(gene_means[nonzero_mask])
            log_vars = np.log10(gene_vars[nonzero_mask])

            # Simple linear fit
            coeffs = np.polyfit(log_means, log_vars, 1)

            # Only apply to genes with non-zero mean
            predicted_vars = np.zeros(n_genes)
            predicted_vars[nonzero_mask] = 10 ** (coeffs[0] * np.log10(gene_means[nonzero_mask]) + coeffs[1])

            # For zero-mean genes, set predicted variance to small value
            predicted_vars[~nonzero_mask] = 1e-10

            # Calculate standardized variance
            standardized_vars = gene_vars / np.maximum(predicted_vars, 1e-10)

            # Select top genes by standardized variance (only consider valid genes)
            valid_mask = gene_means > 0
            valid_std_vars = standardized_vars[valid_mask]
            sorted_indices = np.argsort(valid_std_vars)[::-1]

            # Map back to original gene indices
            valid_gene_indices = np.where(valid_mask)[0]
            n_valid = min(n_top_genes, len(valid_gene_indices))
            top_indices = valid_gene_indices[sorted_indices[:n_valid]]

            # If we need more, add zero-mean genes by variance
            if n_valid < n_top_genes:
                zero_mean_indices = np.where(~valid_mask)[0]
                zero_mean_vars = gene_vars[~valid_mask]
                remaining = n_top_genes - n_valid
                additional = zero_mean_indices[np.argsort(zero_mean_vars)[::-1][:remaining]]
                top_indices = np.concatenate([top_indices, additional])

    elif flavor == "cell_ranger":
        # Cell Ranger-style HVG selection
        # Use coefficient of variation squared
        gene_means = np.mean(X, axis=0)
        gene_vars = np.var(X, axis=0)

        # Avoid division by zero
        gene_means = np.maximum(gene_means, 1e-10)

        cv_squared = gene_vars / (gene_means**2)

        # Select top genes by CV^2
        top_indices = np.argsort(cv_squared)[::-1][:n_top_genes]

    elif flavor == "variance":
        # Simple variance-based selection
        gene_vars = np.var(X, axis=0)
        top_indices = np.argsort(gene_vars)[::-1][:n_top_genes]

    elif flavor == "seurat_v3":
        # Variance of the standardised, clipped raw counts; ties (equal statistic) keep gene order
        seurat_v3_variances = _seurat_v3_normalized_variance(X)
        top_indices = np.argsort(-seurat_v3_variances, kind="stable")[:n_top_genes]

    else:
        raise ValueError(
            f"Unknown HVG selection method: {flavor}. Use 'seurat', 'cell_ranger', 'seurat_v3' or 'variance'"
        )

    # Mark highly variable genes
    highly_variable = np.zeros(n_genes, dtype=bool)
    highly_variable[top_indices] = True

    # Calculate gene statistics for all methods
    gene_means_all = np.mean(X, axis=0)
    gene_vars_all = np.var(X, axis=0)
    # Dispersions as coefficient of variation
    gene_means_safe = np.maximum(gene_means_all, 1e-10)
    dispersions = gene_vars_all / (gene_means_safe**2)
    # Normalized variance: the variance-stabilised statistic for seurat_v3, otherwise a z-score of the variance
    if seurat_v3_variances is not None:
        variances_norm = seurat_v3_variances
    else:
        var_mean = np.mean(gene_vars_all)
        var_std = np.std(gene_vars_all)
        if var_std > 0:
            variances_norm = (gene_vars_all - var_mean) / var_std
        else:
            variances_norm = np.zeros_like(gene_vars_all)

    result = data if inplace else data.copy()

    if hasattr(result, "var") and result.var is not None:
        result.var["highly_variable"] = highly_variable
        result.var["means"] = gene_means_all
        result.var["variances"] = gene_vars_all
        result.var["variances_norm"] = variances_norm
        result.var["dispersions"] = dispersions
    else:
        # Create var dataframe if it doesn't exist
        var_df = pd.DataFrame(index=range(n_genes))
        var_df["highly_variable"] = highly_variable
        var_df["means"] = gene_means_all
        var_df["variances"] = gene_vars_all
        var_df["variances_norm"] = variances_norm
        var_df["dispersions"] = dispersions
        result.var = var_df

    result.uns["hvg"] = {"flavor": flavor, "n_top_genes": n_top_genes, "n_selected": len(top_indices)}

    logger.info(f"Selected {len(top_indices)} highly variable genes using {flavor} method")
    return result
