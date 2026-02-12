"""Single-cell RNA-seq data preprocessing and quality control.

This module provides functions for loading single-cell data, performing quality
control, filtering cells and genes, and basic preprocessing steps. It is designed
to be compatible with AnnData structures while providing graceful fallbacks.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from metainformant.core.data import validation
from metainformant.core.utils import errors
from metainformant.core.io import paths
from metainformant.core.utils import logging
from metainformant.core import io

# Import scipy.sparse for sparse matrix support
try:
    from scipy import sparse as sp_sparse

    HAS_SCIPY_SPARSE = True
except ImportError:
    sp_sparse = None
    HAS_SCIPY_SPARSE = False

logger = logging.get_logger(__name__)

# Try to import anndata, but provide fallback if not available
try:
    import anndata as ad

    HAS_ANNDATA = True
except ImportError:
    ad = None
    HAS_ANNDATA = False


# Custom SingleCellData class as fallback when anndata is not available
class SingleCellData:
    """Fallback single-cell data structure when anndata is not available."""

    def __init__(
        self,
        X: np.ndarray,
        obs: pd.DataFrame | None = None,
        var: pd.DataFrame | None = None,
        uns: Dict[str, Any] | None = None,
        obsm: Dict[str, np.ndarray] | None = None,
        varm: Dict[str, np.ndarray] | None = None,
        obsp: Dict[str, np.ndarray] | None = None,
        varp: Dict[str, np.ndarray] | None = None,
        layers: Dict[str, np.ndarray] | None = None,
    ):
        # Validate dimensions if obs/var provided
        if obs is not None and len(obs) != X.shape[0]:
            raise ValueError(f"obs length ({len(obs)}) must match number of rows in X ({X.shape[0]})")
        if var is not None and len(var) != X.shape[1]:
            raise ValueError(f"var length ({len(var)}) must match number of columns in X ({X.shape[1]})")

        self.X = X  # Expression matrix (cells x genes)
        self.obs = obs if obs is not None else pd.DataFrame(index=range(X.shape[0]))  # Cell metadata
        self.var = var if var is not None else pd.DataFrame(index=range(X.shape[1]))  # Gene metadata
        self.uns = uns if uns is not None else {}  # Unstructured metadata
        self.obsm = obsm if obsm is not None else {}  # Multi-dimensional obs annotations (e.g., embeddings)
        self.varm = varm if varm is not None else {}  # Multi-dimensional var annotations
        self.obsp = obsp if obsp is not None else {}  # Pairwise obs annotations
        self.varp = varp if varp is not None else {}  # Pairwise var annotations
        self.layers = layers if layers is not None else {}  # Different matrix representations

    @property
    def n_obs(self) -> int:
        return self.X.shape[0]

    @property
    def n_vars(self) -> int:
        return self.X.shape[1]

    @property
    def shape(self) -> Tuple[int, int]:
        return self.X.shape

    def copy(self) -> SingleCellData:
        """Create a copy of the data."""
        # Handle sparse matrix copy properly
        if hasattr(self.X, "copy"):
            X_copy = self.X.copy()
        else:
            X_copy = self.X
        return SingleCellData(
            X=X_copy,
            obs=self.obs.copy() if self.obs is not None else None,
            var=self.var.copy() if self.var is not None else None,
            uns=self.uns.copy(),
            obsm={k: v.copy() for k, v in self.obsm.items()} if self.obsm else None,
            varm={k: v.copy() for k, v in self.varm.items()} if self.varm else None,
            obsp={k: v.copy() for k, v in self.obsp.items()} if self.obsp else None,
            varp={k: v.copy() for k, v in self.varp.items()} if self.varp else None,
            layers={k: v.copy() for k, v in self.layers.items()} if self.layers else None,
        )

    def to_df(self) -> pd.DataFrame:
        """Convert to pandas DataFrame."""
        X_dense = self.X.toarray() if hasattr(self.X, "toarray") else self.X
        df = pd.DataFrame(X_dense, index=self.obs.index, columns=self.var.index)
        return df


def load_count_matrix(filepath: str | Path, format: str = "h5ad", **kwargs) -> SingleCellData:
    """Load single-cell count matrix from file.

    Args:
        filepath: Path to the data file
        format: File format ("h5ad", "csv", "tsv", "mtx")
        **kwargs: Additional arguments passed to the loader

    Returns:
        SingleCellData object

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If format is unsupported
    """
    filepath = Path(filepath)
    validation.validate_path_is_file(filepath, "count matrix file")

    logger.info(f"Loading single-cell data from {filepath} (format: {format})")

    if format == "h5ad":
        if not HAS_ANNDATA:
            from metainformant.core.utils.optional_deps import warn_optional_dependency

            warn_optional_dependency("anndata", "H5AD file format support")
            raise errors.ConfigurationError("h5ad format requires anndata package")

        try:
            adata = ad.read_h5ad(filepath, **kwargs)
            # Convert to our SingleCellData format
            return _annadata_to_singlecelldata(adata)
        except Exception as e:
            logger.error(f"Failed to load h5ad file: {e}")
            raise errors.FileIOError(f"Could not load h5ad file: {e}") from e

    elif format in ["csv", "tsv"]:
        separator = "\t" if format == "tsv" else ","
        transpose = kwargs.pop("transpose", False)
        try:
            df = pd.read_csv(filepath, sep=separator, index_col=0, **kwargs)
            if transpose:
                df = df.T  # Transpose so rows are cells, columns are genes
            X = df.values
            # Convert to sparse matrix if scipy is available
            if HAS_SCIPY_SPARSE:
                X = sp_sparse.csr_matrix(X)
            obs = pd.DataFrame(index=df.index)
            var = pd.DataFrame(index=df.columns)
            return SingleCellData(X=X, obs=obs, var=var)
        except Exception as e:
            logger.error(f"Failed to load {format} file: {e}")
            raise errors.FileIOError(f"Could not load {format} file: {e}") from e

    elif format == "mtx":
        # Matrix Market format - requires genes and barcodes files
        try:
            matrix_file = filepath
            genes_file = kwargs.get("genes_file")
            barcodes_file = kwargs.get("barcodes_file")

            if not genes_file or not barcodes_file:
                raise ValueError("Matrix Market format requires genes_file and barcodes_file parameters")

            # Load matrix
            from scipy.io import mmread

            X = mmread(matrix_file).T  # Transpose to cells x genes, keep sparse
            # Convert to csr_matrix if not already
            if HAS_SCIPY_SPARSE and not sp_sparse.issparse(X):
                X = sp_sparse.csr_matrix(X)
            elif HAS_SCIPY_SPARSE:
                X = sp_sparse.csr_matrix(X)

            # Load gene names
            genes_df = pd.read_csv(genes_file, header=None)
            gene_names = genes_df.iloc[:, 0].values

            # Load cell barcodes
            barcodes_df = pd.read_csv(barcodes_file, header=None)
            cell_barcodes = barcodes_df.iloc[:, 0].values

            obs = pd.DataFrame(index=cell_barcodes)
            var = pd.DataFrame(index=gene_names)

            return SingleCellData(X=X, obs=obs, var=var)

        except Exception as e:
            logger.error(f"Failed to load Matrix Market files: {e}")
            raise errors.FileIOError(f"Could not load Matrix Market files: {e}") from e

    else:
        raise errors.ValidationError(f"Unsupported format: {format}")


def _annadata_to_singlecelldata(adata: Any) -> SingleCellData:
    """Convert AnnData object to SingleCellData."""
    return SingleCellData(
        X=adata.X,
        obs=adata.obs,
        var=adata.var,
        uns=adata.uns,
        obsm=dict(adata.obsm) if hasattr(adata, "obsm") else None,
        varm=dict(adata.varm) if hasattr(adata, "varm") else None,
        obsp=dict(adata.obsp) if hasattr(adata, "obsp") else None,
        varp=dict(adata.varp) if hasattr(adata, "varp") else None,
        layers=dict(adata.layers) if hasattr(adata, "layers") else None,
    )


def calculate_qc_metrics(data: SingleCellData) -> SingleCellData:
    """Calculate quality control metrics for single-cell data.

    Args:
        data: SingleCellData object

    Returns:
        SingleCellData with QC metrics added to obs

    Raises:
        TypeError: If data is not SingleCellData
    """
    validation.validate_type(data, SingleCellData, "data")

    logger.info(f"Calculating QC metrics for {data.n_obs} cells and {data.n_vars} genes")

    # Create a copy to avoid modifying original
    result = data.copy()

    # Basic QC metrics
    X_dense = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Per-cell metrics
    total_counts = np.sum(X_dense, axis=1)  # Total counts per cell
    n_genes = np.sum(X_dense > 0, axis=1)  # Number of genes expressed per cell
    pct_mt = np.zeros(data.n_obs)  # Mitochondrial percentage
    pct_ribo = np.zeros(data.n_obs)  # Ribosomal percentage

    # Detect mitochondrial and ribosomal genes (simple heuristic)
    if data.var is not None and hasattr(data.var, "index"):
        # Convert index to string for comparison
        gene_names = [str(gene) for gene in data.var.index]

        # Mitochondrial genes
        mt_genes = [i for i, gene in enumerate(gene_names) if gene.upper().startswith(("MT-", "MT.", "MT_"))]
        if mt_genes:
            mt_counts = np.sum(X_dense[:, mt_genes], axis=1)
            pct_mt = 100 * mt_counts / np.where(total_counts > 0, total_counts, 1)
            pct_mt = np.nan_to_num(pct_mt, nan=0.0)

        # Ribosomal genes (RPS and RPL)
        ribo_genes = [i for i, gene in enumerate(gene_names) if gene.upper().startswith(("RPS", "RPL", "MRPS", "MRPL"))]
        if ribo_genes:
            ribo_counts = np.sum(X_dense[:, ribo_genes], axis=1)
            pct_ribo = 100 * ribo_counts / np.where(total_counts > 0, total_counts, 1)
            pct_ribo = np.nan_to_num(pct_ribo, nan=0.0)

    # Per-gene metrics
    n_cells = np.sum(X_dense > 0, axis=0)  # Number of cells expressing each gene
    gene_total_counts = np.sum(X_dense, axis=0)  # Total counts per gene
    mean_expression = np.mean(X_dense, axis=0)
    pct_dropout = 100 * (1 - n_cells / data.n_obs)

    # Add to obs (cell metrics)
    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    result.obs["total_counts"] = total_counts
    result.obs["n_genes"] = n_genes
    result.obs["pct_mt"] = pct_mt
    result.obs["pct_ribo"] = pct_ribo

    # Add to var (gene metrics)
    result.var = result.var.copy() if result.var is not None else pd.DataFrame()
    result.var["n_cells"] = n_cells
    result.var["total_counts"] = gene_total_counts
    result.var["mean_expression"] = mean_expression
    result.var["pct_dropout"] = pct_dropout

    # Store summary statistics in uns
    result.uns["qc_summary"] = {
        "total_cells": data.n_obs,
        "total_genes": data.n_vars,
        "mean_counts_per_cell": float(np.mean(total_counts)),
        "mean_genes_per_cell": float(np.mean(n_genes)),
        "mean_mito_pct": float(np.mean(pct_mt)),
        "median_counts_per_cell": float(np.median(total_counts)),
        "median_genes_per_cell": float(np.median(n_genes)),
    }

    logger.info("QC metrics calculated successfully")
    return result


def filter_cells(
    data: SingleCellData,
    min_counts: int | None = None,
    max_counts: int | None = None,
    min_genes: int | None = None,
    max_genes: int | None = None,
    max_pct_mt: float | None = None,
    max_mito_percent: float | None = None,  # Alias for max_pct_mt
) -> SingleCellData:
    """Filter cells based on QC metrics.

    Args:
        data: SingleCellData object with QC metrics
        min_counts: Minimum total counts per cell
        max_counts: Maximum total counts per cell
        min_genes: Minimum number of genes per cell
        max_genes: Maximum number of genes per cell
        max_pct_mt: Maximum mitochondrial percentage
        max_mito_percent: Alias for max_pct_mt (deprecated)

    Returns:
        Filtered SingleCellData

    Raises:
        TypeError: If data is not SingleCellData
    """
    validation.validate_type(data, SingleCellData, "data")

    logger.info("Filtering cells based on QC metrics")

    # Handle alias
    if max_mito_percent is not None and max_pct_mt is None:
        max_pct_mt = max_mito_percent

    # Check if QC metrics are available
    if "total_counts" not in data.obs.columns:
        logger.warning("QC metrics not found in data.obs. Run calculate_qc_metrics first.")
        data = calculate_qc_metrics(data)

    # Create filter mask
    keep_cells = np.ones(data.n_obs, dtype=bool)

    if min_counts is not None:
        keep_cells &= np.asarray(data.obs["total_counts"] >= min_counts)
    if max_counts is not None:
        keep_cells &= np.asarray(data.obs["total_counts"] <= max_counts)
    if min_genes is not None:
        keep_cells &= np.asarray(data.obs["n_genes"] >= min_genes)
    if max_genes is not None:
        keep_cells &= np.asarray(data.obs["n_genes"] <= max_genes)
    if max_pct_mt is not None:
        keep_cells &= np.asarray(data.obs["pct_mt"] <= max_pct_mt)

    n_kept = np.sum(keep_cells)
    n_filtered = data.n_obs - n_kept

    logger.info(f"Filtered {n_filtered} cells, keeping {n_kept} cells")

    # Apply filter - handle sparse matrices
    X_filtered = data.X[keep_cells, :]
    if HAS_SCIPY_SPARSE and sp_sparse.issparse(X_filtered):
        X_filtered = sp_sparse.csr_matrix(X_filtered)

    result = SingleCellData(
        X=X_filtered,
        obs=data.obs.iloc[keep_cells].copy() if data.obs is not None else None,
        var=data.var.copy() if data.var is not None else None,
        uns=data.uns.copy(),
    )

    # Update QC summary
    if "qc_summary" in result.uns:
        result.uns["qc_summary"]["cells_after_filtering"] = n_kept
        result.uns["qc_summary"]["cells_filtered"] = n_filtered

    return result


def filter_genes(
    data: SingleCellData,
    min_cells: int | None = None,
    max_cells: int | None = None,
    min_counts: int | None = None,
    max_counts: int | None = None,
) -> SingleCellData:
    """Filter genes based on expression criteria.

    Args:
        data: SingleCellData object
        min_cells: Minimum number of cells expressing the gene
        max_cells: Maximum number of cells expressing the gene
        min_counts: Minimum total counts for the gene
        max_counts: Maximum total counts for the gene

    Returns:
        Filtered SingleCellData

    Raises:
        TypeError: If data is not SingleCellData
    """
    validation.validate_type(data, SingleCellData, "data")

    logger.info("Filtering genes based on expression criteria")

    # Check if QC metrics are available
    if "n_cells" not in data.var.columns:
        logger.warning("Gene QC metrics not found in data.var. Run calculate_qc_metrics first.")
        data = calculate_qc_metrics(data)

    # Create filter mask
    keep_genes = np.ones(data.n_vars, dtype=bool)

    if min_cells is not None:
        keep_genes &= np.asarray(data.var["n_cells"] >= min_cells)
    if max_cells is not None:
        keep_genes &= np.asarray(data.var["n_cells"] <= max_cells)
    if min_counts is not None:
        keep_genes &= np.asarray(data.var["total_counts"] >= min_counts)
    if max_counts is not None:
        keep_genes &= np.asarray(data.var["total_counts"] <= max_counts)

    n_kept = np.sum(keep_genes)
    n_filtered = data.n_vars - n_kept

    logger.info(f"Filtered {n_filtered} genes, keeping {n_kept} genes")

    # Apply filter - handle sparse matrices
    X_filtered = data.X[:, keep_genes]
    if HAS_SCIPY_SPARSE and sp_sparse.issparse(X_filtered):
        X_filtered = sp_sparse.csr_matrix(X_filtered)

    result = SingleCellData(
        X=X_filtered,
        obs=data.obs.copy() if data.obs is not None else None,
        var=data.var.iloc[keep_genes].copy() if data.var is not None else None,
        uns=data.uns.copy(),
    )

    # Update QC summary
    if "qc_summary" in result.uns:
        result.uns["qc_summary"]["genes_after_filtering"] = n_kept
        result.uns["qc_summary"]["genes_filtered"] = n_filtered

    return result


def normalize_counts(
    data: SingleCellData,
    target_sum: float | None = None,
    method: str = "total_count",
    normalize_method: str | None = None,  # Deprecated alias
) -> SingleCellData:
    """Normalize gene expression counts.

    Args:
        data: SingleCellData object
        target_sum: Target sum for normalization (default: median total counts)
        method: Normalization method ("total_count", "median", "size_factors")
        normalize_method: Deprecated alias for method

    Returns:
        Normalized SingleCellData

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If method is invalid
    """
    validation.validate_type(data, SingleCellData, "data")

    # Handle deprecated alias
    if normalize_method is not None:
        method = normalize_method

    valid_methods = ["total_count", "total", "median", "size_factors"]
    if method not in valid_methods:
        raise ValueError(f"Unknown normalization method: {method}. Valid methods: {valid_methods}")

    logger.info(f"Normalizing counts using {method} method")

    # Create copy
    result = data.copy()

    # Convert to dense if needed
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X
    X_normalized = X.astype(float)

    if method in ["total_count", "total"]:
        # Library size normalization
        lib_sizes = np.sum(X, axis=1)

        if target_sum is None:
            target_sum = np.median(lib_sizes[lib_sizes > 0]) if np.any(lib_sizes > 0) else 1.0

        # Normalize to target sum (avoid division by zero)
        with np.errstate(divide="ignore", invalid="ignore"):
            size_factors = target_sum / lib_sizes
            size_factors = np.where(np.isfinite(size_factors), size_factors, 0)
        X_normalized = X * size_factors[:, np.newaxis]

    elif method == "median":
        # Median normalization
        lib_sizes = np.sum(X, axis=1)
        median_lib_size = np.median(lib_sizes[lib_sizes > 0]) if np.any(lib_sizes > 0) else 1.0

        with np.errstate(divide="ignore", invalid="ignore"):
            size_factors = median_lib_size / lib_sizes
            size_factors = np.where(np.isfinite(size_factors), size_factors, 0)
        X_normalized = X * size_factors[:, np.newaxis]
        target_sum = median_lib_size

    elif method == "size_factors":
        # Estimate size factors using geometric mean of gene expression
        # Simplified implementation
        log_expr = np.log(X + 1)
        geometric_means = np.exp(np.mean(log_expr, axis=0))

        # Size factors based on median ratio
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios = X / geometric_means[np.newaxis, :]
            ratios = np.where(np.isfinite(ratios), ratios, 0)
        size_factors = np.median(ratios, axis=1)

        if target_sum is None:
            target_sum = np.median(size_factors[size_factors > 0]) if np.any(size_factors > 0) else 1.0

        with np.errstate(divide="ignore", invalid="ignore"):
            scale = target_sum / size_factors
            scale = np.where(np.isfinite(scale), scale, 0)
        X_normalized = X * scale[:, np.newaxis]

    # Store normalization info
    result.uns["normalization"] = {
        "method": method,
        "target_sum": float(target_sum) if target_sum is not None else None,
        "size_factors": size_factors.tolist() if "size_factors" in locals() else None,
    }

    result.X = X_normalized
    logger.info("Normalization completed")
    return result


def log_transform(data: SingleCellData, base: float = np.e) -> SingleCellData:
    """Log-transform expression values.

    Args:
        data: SingleCellData object
        base: Logarithm base (default: natural log)

    Returns:
        Log-transformed SingleCellData

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If base is invalid
    """
    validation.validate_type(data, SingleCellData, "data")

    if base <= 0 or base == 1:
        raise errors.ValidationError("Log base must be positive and not equal to 1")

    logger.info(f"Log-transforming data with base {base}")

    # Create copy
    result = data.copy()

    # Convert to dense if needed
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    if base == np.e:
        X_log = np.log(X + 1)  # Add pseudocount
        transform_type = "natural_log"
    elif base == 2:
        X_log = np.log2(X + 1)
        transform_type = "log2"
    elif base == 10:
        X_log = np.log10(X + 1)
        transform_type = "log10"
    else:
        X_log = np.log(X + 1) / np.log(base)
        transform_type = f"log{base}"

    result.X = X_log

    # Store transformation info - tests expect these specific keys
    result.uns["log_transformed"] = True
    result.uns["log_base"] = base
    result.uns["transformation"] = {
        "type": transform_type,
        "base": base,
        "pseudocount": 1,
    }

    logger.info("Log transformation completed")
    return result


def scale_data(data: SingleCellData, zero_center: bool = True, max_value: float | None = None) -> SingleCellData:
    """Scale gene expression data.

    Args:
        data: SingleCellData object
        zero_center: Whether to center data around zero
        max_value: Clip values above this threshold

    Returns:
        Scaled SingleCellData

    Raises:
        TypeError: If data is not SingleCellData
    """
    validation.validate_type(data, SingleCellData, "data")

    logger.info("Scaling data")

    # Create copy
    result = data.copy()

    # Convert to dense if needed
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X
    X_scaled = X.astype(float)

    # Zero-center if requested
    if zero_center:
        gene_means = np.mean(X_scaled, axis=0)
        X_scaled = X_scaled - gene_means[np.newaxis, :]

    # Scale by gene standard deviation
    gene_stds = np.std(X_scaled, axis=0)
    gene_stds = np.where(gene_stds == 0, 1, gene_stds)  # Avoid division by zero
    X_scaled = X_scaled / gene_stds[np.newaxis, :]

    # Clip extreme values
    if max_value is not None:
        X_scaled = np.clip(X_scaled, -max_value, max_value)

    result.X = X_scaled

    # Store scaling info - tests expect "scaled" key
    result.uns["scaled"] = True
    result.uns["scaling"] = {
        "zero_centered": zero_center,
        "max_value": max_value,
        "gene_means": gene_means.tolist() if zero_center else None,
        "gene_stds": gene_stds.tolist(),
    }

    logger.info("Data scaling completed")
    return result


def identify_highly_variable_genes(
    data: SingleCellData, n_top_genes: int = 2000, flavor: str = "seurat"
) -> SingleCellData:
    """Identify highly variable genes.

    Args:
        data: SingleCellData object (should be normalized and log-transformed)
        n_top_genes: Number of top variable genes to select
        flavor: Method for calculating variance ("seurat", "cell_ranger")

    Returns:
        SingleCellData with highly variable genes marked

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If flavor is invalid
    """
    validation.validate_type(data, SingleCellData, "data")

    if flavor not in ["seurat", "cell_ranger"]:
        raise errors.ValidationError(f"Unsupported flavor: {flavor}")

    logger.info(f"Identifying highly variable genes using {flavor} method")

    # Create copy
    result = data.copy()

    # Convert to dense if needed
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Calculate gene statistics
    gene_means = np.mean(X, axis=0)
    gene_vars = np.var(X, axis=0)
    gene_cv2 = gene_vars / (gene_means**2 + 1e-10)  # Coefficient of variation squared

    if flavor == "seurat":
        # Seurat method: fit mean-variance relationship
        # Simplified implementation
        log_means = np.log10(gene_means + 1e-10)

        # Fit polynomial to mean-variance relationship
        coeffs = np.polyfit(log_means, gene_cv2, 2)
        fitted_cv2 = np.polyval(coeffs, log_means)

        # Calculate standardized variance
        standardized_var = (gene_cv2 - fitted_cv2) / fitted_cv2

    elif flavor == "cell_ranger":
        # Cell Ranger method: simple dispersion-based selection
        standardized_var = gene_cv2

    # Rank genes by variability
    valid_genes = ~np.isnan(standardized_var) & ~np.isinf(standardized_var)
    sorted_indices = np.argsort(standardized_var[valid_genes])[::-1]

    # Mark top variable genes
    highly_variable = np.zeros(data.n_vars, dtype=bool)
    top_gene_indices = np.where(valid_genes)[0][sorted_indices[:n_top_genes]]
    highly_variable[top_gene_indices] = True

    # Add to var DataFrame
    result.var = result.var.copy() if result.var is not None else pd.DataFrame()
    result.var["highly_variable"] = highly_variable
    result.var["gene_mean"] = gene_means
    result.var["gene_variance"] = gene_vars
    result.var["gene_cv2"] = gene_cv2

    if flavor == "seurat":
        result.var["fitted_cv2"] = fitted_cv2
        result.var["standardized_variance"] = standardized_var

    # Store HVG info
    result.uns["highly_variable_genes"] = {
        "method": flavor,
        "n_top_genes": n_top_genes,
        "n_selected": np.sum(highly_variable),
    }

    logger.info(f"Identified {np.sum(highly_variable)} highly variable genes")
    return result


def remove_batch_effects(data: SingleCellData, batch_key: str, method: str = "regress_out") -> SingleCellData:
    """Remove batch effects from expression data.

    Args:
        data: SingleCellData object
        batch_key: Column in obs containing batch information
        method: Method for batch correction ("regress_out", "combat")

    Returns:
        Batch-corrected SingleCellData

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If batch_key not found or method invalid
    """
    validation.validate_type(data, SingleCellData, "data")

    if method not in ["regress_out", "combat"]:
        raise errors.ValidationError(f"Unsupported batch correction method: {method}")

    if data.obs is None or batch_key not in data.obs.columns:
        raise errors.ValidationError(f"Batch key '{batch_key}' not found in data.obs")

    logger.info(f"Removing batch effects using {method} method")

    # Create copy
    result = data.copy()

    # Convert to dense if needed
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    if method == "regress_out":
        # Simple linear regression of batch effects
        # This is a simplified implementation
        batch_labels = data.obs[batch_key].values
        unique_batches = np.unique(batch_labels)

        X_corrected = X.copy()

        for batch in unique_batches:
            batch_mask = batch_labels == batch
            if np.sum(batch_mask) > 0:
                # Subtract batch-specific mean
                batch_mean = np.mean(X[batch_mask, :], axis=0)
                X_corrected[batch_mask, :] -= batch_mean[np.newaxis, :]

    elif method == "combat":
        # Simplified ComBat-like method
        # This is a highly simplified implementation
        batch_labels = data.obs[batch_key].values
        unique_batches = np.unique(batch_labels)

        X_corrected = X.copy()

        # Estimate batch effects
        for gene_idx in range(X.shape[1]):
            gene_expr = X[:, gene_idx]

            # Simple standardization per batch
            for batch in unique_batches:
                batch_mask = batch_labels == batch
                if np.sum(batch_mask) > 1:  # Need at least 2 samples
                    batch_expr = gene_expr[batch_mask]
                    batch_mean = np.mean(batch_expr)
                    batch_std = np.std(batch_expr)

                    if batch_std > 0:
                        # Standardize within batch
                        X_corrected[batch_mask, gene_idx] = (batch_expr - batch_mean) / batch_std

    result.X = X_corrected

    # Store batch correction info
    result.uns["batch_correction"] = {
        "method": method,
        "batch_key": batch_key,
        "n_batches": len(np.unique(data.obs[batch_key])),
    }

    logger.info("Batch effect correction completed")
    return result
