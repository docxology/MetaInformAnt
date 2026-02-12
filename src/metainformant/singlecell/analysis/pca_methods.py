"""Linear dimensionality reduction methods for single-cell data.

This module provides PCA, ICA, and Factor Analysis methods for reducing
the dimensionality of single-cell expression data, along with highly
variable gene (HVG) selection.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from metainformant.core.data import validation
from metainformant.core.utils import errors
from metainformant.core.utils import logging

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

# Import our SingleCellData
from metainformant.singlecell.data.preprocessing import SingleCellData


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
        n_components: Number of principal components to compute
        random_state: Random seed for reproducibility
        scale_data: Whether to scale data before PCA
        use_hvgs: Whether using only highly variable genes

    Returns:
        SingleCellData with PCA coordinates added

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If n_components is invalid
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn is required for PCA. " "Install with: uv pip install scikit-learn")

    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=2, max_val=min(data.n_obs, data.n_vars), name="n_components")

    logger.info(f"Performing PCA with {n_components} components")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

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
    if result.obs is None:
        result.obs = pca_coords
    else:
        for col in pca_coords.columns:
            result.obs[col] = pca_coords[col]

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
        "components_shape": pca.components_.shape,
        "random_state": random_state,
        "scaled": scale_data,
        "use_hvgs": use_hvgs,
    }

    # Store components for gene loadings
    result.varm = result.varm if hasattr(result, "varm") and result.varm is not None else {}
    result.varm["PCs"] = pca.components_.T

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
        TypeError: If data is not SingleCellData
    """
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

    # Add to obs
    if result.obs is None:
        result.obs = ica_coords
    else:
        for col in ica_coords.columns:
            result.obs[col] = ica_coords[col]

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
        TypeError: If data is not SingleCellData
    """
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

    # Add to obs
    if result.obs is None:
        result.obs = fa_coords
    else:
        for col in fa_coords.columns:
            result.obs[col] = fa_coords[col]

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
        use_hvgs: Whether to use only highly variable genes

    Returns:
        SingleCellData with PCA results
    """
    # Note: use_hvgs is just a flag - the data should already have HVGs marked
    # This parameter exists for API compatibility
    return pca_reduction(data, n_components, random_state, scale_data, use_hvgs=use_hvgs)


def select_hvgs(
    data: SingleCellData,
    n_top_genes: int = 2000,
    flavor: str | None = None,
    *,
    method: str | None = None,
) -> SingleCellData:
    """Select highly variable genes from single-cell data.

    Args:
        data: SingleCellData object
        n_top_genes: Number of top variable genes to select
        flavor: Method for HVG selection ('seurat', 'cell_ranger', 'seurat_v3', 'variance')
        method: Alias for flavor

    Returns:
        SingleCellData with highly variable genes marked
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

    else:
        raise ValueError(f"Unknown HVG selection method: {flavor}. Use 'seurat', 'cell_ranger', or 'variance'")

    # Mark highly variable genes
    highly_variable = np.zeros(n_genes, dtype=bool)
    highly_variable[top_indices] = True

    # Calculate gene statistics for all methods
    gene_means_all = np.mean(X, axis=0)
    gene_vars_all = np.var(X, axis=0)
    # Dispersions as coefficient of variation
    gene_means_safe = np.maximum(gene_means_all, 1e-10)
    dispersions = gene_vars_all / (gene_means_safe**2)
    # Normalized variance (z-score of variance)
    var_mean = np.mean(gene_vars_all)
    var_std = np.std(gene_vars_all)
    if var_std > 0:
        variances_norm = (gene_vars_all - var_mean) / var_std
    else:
        variances_norm = np.zeros_like(gene_vars_all)

    if hasattr(data, "var") and data.var is not None:
        data.var["highly_variable"] = highly_variable
        data.var["means"] = gene_means_all
        data.var["variances"] = gene_vars_all
        data.var["variances_norm"] = variances_norm
        data.var["dispersions"] = dispersions
    else:
        # Create var dataframe if it doesn't exist
        import pandas as pd

        var_df = pd.DataFrame(index=range(n_genes))
        var_df["highly_variable"] = highly_variable
        var_df["means"] = gene_means_all
        var_df["variances"] = gene_vars_all
        var_df["variances_norm"] = variances_norm
        var_df["dispersions"] = dispersions
        data.var = var_df

    data.uns["hvg"] = {"flavor": flavor, "n_top_genes": n_top_genes, "n_selected": len(top_indices)}

    logger.info(f"Selected {len(top_indices)} highly variable genes using {flavor} method")
    return data
