"""Single-cell data preprocessing and quality control.

This module provides essential preprocessing functions for single-cell genomics data,
including loading, quality control metrics, filtering, and normalization.
All implementations use real computational methods without mocking.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union, Any
import warnings

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import zscore

from ..core.io import open_text_auto


class SingleCellData:
    """Container for single-cell expression data and metadata.
    
    Provides a simplified interface similar to AnnData for single-cell analysis,
    designed to work with standard Python libraries.
    """
    
    def __init__(
        self, 
        X: Union[np.ndarray, sparse.spmatrix],
        obs: Optional[pd.DataFrame] = None,
        var: Optional[pd.DataFrame] = None,
        obsm: Optional[Dict[str, np.ndarray]] = None,
        varm: Optional[Dict[str, np.ndarray]] = None,
        uns: Optional[Dict[str, Any]] = None,
    ):
        """Initialize SingleCellData object.
        
        Args:
            X: Gene expression matrix (cells × genes)
            obs: Cell metadata (observations)
            var: Gene metadata (variables)  
            obsm: Multi-dimensional cell annotations (e.g., PCA, UMAP coordinates)
            varm: Multi-dimensional gene annotations
            uns: Unstructured metadata
        """
        self.X = X
        self.n_obs, self.n_vars = X.shape
        
        # Initialize metadata
        self.obs = obs if obs is not None else pd.DataFrame(index=range(self.n_obs))
        self.var = var if var is not None else pd.DataFrame(index=range(self.n_vars))
        self.obsm = obsm if obsm is not None else {}
        self.varm = varm if varm is not None else {}
        self.uns = uns if uns is not None else {}
        
        # Validate dimensions
        if len(self.obs) != self.n_obs:
            raise ValueError(f"obs length ({len(self.obs)}) doesn't match X rows ({self.n_obs})")
        if len(self.var) != self.n_vars:
            raise ValueError(f"var length ({len(self.var)}) doesn't match X columns ({self.n_vars})")
    
    def copy(self) -> 'SingleCellData':
        """Create a deep copy of the data."""
        return SingleCellData(
            X=self.X.copy(),
            obs=self.obs.copy(),
            var=self.var.copy(),
            obsm={k: v.copy() for k, v in self.obsm.items()},
            varm={k: v.copy() for k, v in self.varm.items()},
            uns=self.uns.copy(),
        )


def load_count_matrix(
    path: Union[str, Path],
    format: str = "mtx",
    delimiter: str = "\t",
    transpose: bool = False,
) -> SingleCellData:
    """Load single-cell count matrix from various formats.
    
    Args:
        path: Path to count matrix file
        format: File format ('mtx', 'csv', 'tsv', 'h5', 'h5ad')
        delimiter: Delimiter for text formats
        transpose: Whether to transpose the matrix (genes × cells -> cells × genes)
        
    Returns:
        SingleCellData object with loaded expression matrix
    """
    path = Path(path)
    
    if format == "mtx":
        # Matrix Market format (common for 10X Genomics)
        try:
            from scipy.io import mmread
            X = mmread(path)
            if sparse.issparse(X):
                X = X.tocsr()
            else:
                X = sparse.csr_matrix(X)
        except ImportError:
            raise ImportError("scipy required for Matrix Market format")
            
    elif format in ["csv", "tsv"]:
        # Comma or tab-separated values
        sep = "," if format == "csv" else delimiter
        df = pd.read_csv(path, sep=sep, index_col=0)
        X = sparse.csr_matrix(df.values)
        
        # Create metadata from index/columns
        var = pd.DataFrame(index=df.columns, data={"gene_name": df.columns})
        obs = pd.DataFrame(index=df.index, data={"cell_id": df.index})
        
        # Apply transpose if requested
        if transpose:
            X = X.T
            # Swap obs and var metadata
            var, obs = obs, var
            var.columns = ["gene_name"] if "gene_name" in var.columns else ["cell_id"]
            obs.columns = ["cell_id"] if "cell_id" in obs.columns else ["gene_name"]
        
        return SingleCellData(X=X, obs=obs, var=var)
        
    elif format == "h5":
        # HDF5 format
        try:
            import h5py
            with h5py.File(path, 'r') as f:
                # Assume standard 10X structure
                if 'matrix' in f:
                    group = f['matrix']
                    X = sparse.csr_matrix((
                        group['data'][:],
                        group['indices'][:],
                        group['indptr'][:]
                    ), shape=group.attrs.get('shape', (0, 0)))
                else:
                    raise ValueError("HDF5 file doesn't contain expected 'matrix' group")
        except ImportError:
            raise ImportError("h5py required for HDF5 format")
            
    elif format == "h5ad":
        # AnnData HDF5 format
        try:
            import anndata
            adata = anndata.read_h5ad(path)
            return SingleCellData(
                X=adata.X,
                obs=adata.obs,
                var=adata.var,
                obsm=dict(adata.obsm),
                varm=dict(adata.varm),
                uns=dict(adata.uns),
            )
        except ImportError:
            warnings.warn("anndata not available, using fallback loader")
            # Fallback to HDF5 loader
            return load_count_matrix(path, format="h5", transpose=transpose)
    else:
        raise ValueError(f"Unsupported format: {format}")
    
    # Create SingleCellData object first
    sc_data = SingleCellData(X=X)
    
    if transpose:
        # Transpose both X and metadata
        sc_data.X = X.T
        sc_data.n_obs, sc_data.n_vars = sc_data.X.shape
        # Recreate obs and var with correct dimensions
        sc_data.obs = pd.DataFrame(index=range(sc_data.n_obs))
        sc_data.var = pd.DataFrame(index=range(sc_data.n_vars))
    
    return sc_data


def calculate_qc_metrics(data: SingleCellData) -> SingleCellData:
    """Calculate quality control metrics for single-cell data.
    
    Computes standard QC metrics including:
    - Number of genes per cell (n_genes)
    - Total UMI/read counts per cell (total_counts)
    - Percentage of mitochondrial gene expression (pct_mt)
    - Percentage of ribosomal gene expression (pct_ribo)
    - Number of cells expressing each gene (n_cells)
    - Mean expression per gene (mean_expression)
    
    Args:
        data: SingleCellData object
        
    Returns:
        Updated SingleCellData with QC metrics in obs and var
    """
    data = data.copy()
    X = data.X
    
    if sparse.issparse(X):
        X_dense = X.toarray() if X.nnz < 1e6 else X  # Convert small matrices for efficiency
    else:
        X_dense = X
    
    # Cell-level metrics (observations)
    if sparse.issparse(X_dense):
        # For large sparse matrices
        data.obs['total_counts'] = np.array(X.sum(axis=1)).flatten()
        data.obs['n_genes'] = np.array((X > 0).sum(axis=1)).flatten()
    else:
        # For dense matrices
        data.obs['total_counts'] = X_dense.sum(axis=1)
        data.obs['n_genes'] = (X_dense > 0).sum(axis=1)
    
    # Mitochondrial gene metrics
    if 'gene_name' in data.var.columns:
        mt_genes = data.var['gene_name'].str.startswith(('MT-', 'mt-', 'Mt-'))
        if mt_genes.any():
            mt_indices = np.where(mt_genes.values)[0]  # Convert boolean mask to indices
            if sparse.issparse(X):
                mt_counts = np.array(X[:, mt_indices].sum(axis=1)).flatten()
            else:
                mt_counts = X[:, mt_indices].sum(axis=1)
            data.obs['pct_mt'] = 100 * mt_counts / data.obs['total_counts']
        else:
            data.obs['pct_mt'] = 0.0
    else:
        data.obs['pct_mt'] = 0.0
    
    # Ribosomal gene metrics  
    if 'gene_name' in data.var.columns:
        ribo_genes = data.var['gene_name'].str.startswith(('RPS', 'RPL', 'Rps', 'Rpl'))
        if ribo_genes.any():
            ribo_indices = np.where(ribo_genes.values)[0]  # Convert boolean mask to indices
            if sparse.issparse(X):
                ribo_counts = np.array(X[:, ribo_indices].sum(axis=1)).flatten()
            else:
                ribo_counts = X[:, ribo_indices].sum(axis=1)
            data.obs['pct_ribo'] = 100 * ribo_counts / data.obs['total_counts']
        else:
            data.obs['pct_ribo'] = 0.0
    else:
        data.obs['pct_ribo'] = 0.0
    
    # Gene-level metrics (variables)
    if sparse.issparse(X):
        data.var['n_cells'] = np.array((X > 0).sum(axis=0)).flatten()
        data.var['total_counts'] = np.array(X.sum(axis=0)).flatten()
    else:
        data.var['n_cells'] = (X > 0).sum(axis=0)
        data.var['total_counts'] = X.sum(axis=0)
    
    data.var['mean_expression'] = data.var['total_counts'] / data.n_obs
    
    return data


def filter_cells(
    data: SingleCellData,
    min_genes: int = 200,
    max_genes: Optional[int] = None,
    min_counts: int = 1000,
    max_counts: Optional[int] = None,
    max_pct_mt: float = 20.0,
) -> SingleCellData:
    """Filter cells based on quality control metrics.
    
    Args:
        data: SingleCellData object with QC metrics
        min_genes: Minimum number of genes per cell
        max_genes: Maximum number of genes per cell (None for no limit)
        min_counts: Minimum total counts per cell
        max_counts: Maximum total counts per cell (None for no limit)  
        max_pct_mt: Maximum percentage of mitochondrial expression
        
    Returns:
        Filtered SingleCellData object
    """
    # Ensure QC metrics are calculated
    if 'n_genes' not in data.obs.columns:
        data = calculate_qc_metrics(data)
    
    # Create filter mask
    cell_filter = (
        (data.obs['n_genes'] >= min_genes) &
        (data.obs['total_counts'] >= min_counts) &
        (data.obs['pct_mt'] <= max_pct_mt)
    )
    
    if max_genes is not None:
        cell_filter &= (data.obs['n_genes'] <= max_genes)
    
    if max_counts is not None:
        cell_filter &= (data.obs['total_counts'] <= max_counts)
    
    print(f"Filtering cells: {cell_filter.sum()} / {len(cell_filter)} cells pass filters")
    
    # Apply filter (convert pandas Series to numpy array for sparse matrix compatibility)
    cell_filter_array = cell_filter.values
    
    return SingleCellData(
        X=data.X[cell_filter_array, :],
        obs=data.obs[cell_filter].copy(),
        var=data.var.copy(),
        obsm={k: v[cell_filter_array] for k, v in data.obsm.items()},
        varm=data.varm.copy(),
        uns=data.uns.copy(),
    )


def filter_genes(
    data: SingleCellData,
    min_cells: int = 3,
    min_counts: int = 1,
) -> SingleCellData:
    """Filter genes based on expression criteria.
    
    Args:
        data: SingleCellData object
        min_cells: Minimum number of cells expressing the gene
        min_counts: Minimum total counts for the gene
        
    Returns:
        Filtered SingleCellData object
    """
    # Ensure QC metrics are calculated
    if 'n_cells' not in data.var.columns:
        data = calculate_qc_metrics(data)
    
    # Create filter mask
    gene_filter = (
        (data.var['n_cells'] >= min_cells) &
        (data.var['total_counts'] >= min_counts)
    )
    
    print(f"Filtering genes: {gene_filter.sum()} / {len(gene_filter)} genes pass filters")
    
    # Apply filter (convert pandas Series to numpy array for sparse matrix compatibility)
    gene_filter_array = gene_filter.values
    
    return SingleCellData(
        X=data.X[:, gene_filter_array],
        obs=data.obs.copy(),
        var=data.var[gene_filter].copy(),
        obsm=data.obsm.copy(),
        varm={k: v[gene_filter_array] for k, v in data.varm.items()},
        uns=data.uns.copy(),
    )


def normalize_counts(
    data: SingleCellData,
    target_sum: float = 1e4,
    method: str = "total_count",
) -> SingleCellData:
    """Normalize count data.
    
    Args:
        data: SingleCellData object
        target_sum: Target sum for normalization (e.g., 10,000 for CPM)
        method: Normalization method ('total_count', 'median', 'size_factor')
        
    Returns:
        Normalized SingleCellData object
    """
    data = data.copy()
    X = data.X.copy()
    
    if method == "total_count":
        # Total count normalization (CPM-like)
        if sparse.issparse(X):
            counts_per_cell = np.array(X.sum(axis=1)).flatten()
        else:
            counts_per_cell = X.sum(axis=1)
        
        # Avoid division by zero
        counts_per_cell[counts_per_cell == 0] = 1
        
        # Normalize to target sum
        scale_factors = target_sum / counts_per_cell
        
        if sparse.issparse(X):
            X = X.multiply(scale_factors[:, np.newaxis])
        else:
            X = X * scale_factors[:, np.newaxis]
            
    elif method == "median":
        # Median normalization
        if sparse.issparse(X):
            counts_per_cell = np.array(X.sum(axis=1)).flatten()
        else:
            counts_per_cell = X.sum(axis=1)
        
        median_counts = np.median(counts_per_cell[counts_per_cell > 0])
        scale_factors = median_counts / counts_per_cell
        scale_factors[counts_per_cell == 0] = 0
        
        if sparse.issparse(X):
            X = X.multiply(scale_factors[:, np.newaxis])
        else:
            X = X * scale_factors[:, np.newaxis]
    
    elif method == "size_factor":
        # DESeq2-style size factor normalization
        if sparse.issparse(X):
            # Convert to dense for geometric mean calculation
            X_dense = X.toarray()
        else:
            X_dense = X
        
        # Calculate geometric means per gene (avoiding zeros)
        log_counts = np.log(X_dense + 1)
        geo_means = np.exp(np.mean(log_counts, axis=0))
        geo_means[geo_means == 0] = 1
        
        # Calculate size factors per cell
        ratios = X_dense / geo_means[np.newaxis, :]
        size_factors = np.median(ratios, axis=1)
        size_factors[size_factors == 0] = 1
        
        X = X_dense / size_factors[:, np.newaxis]
        
    else:
        raise ValueError(f"Unknown normalization method: {method}")
    
    data.X = X
    data.uns['normalization'] = {'method': method, 'target_sum': target_sum}
    
    return data


def log_transform(data: SingleCellData, base: float = np.e) -> SingleCellData:
    """Apply log transformation to expression data.
    
    Args:
        data: SingleCellData object (should be normalized)
        base: Logarithm base (e for natural log, 2 for log2, 10 for log10)
        
    Returns:
        Log-transformed SingleCellData object
    """
    data = data.copy()
    
    if sparse.issparse(data.X):
        # For sparse matrices, use log1p
        data.X = data.X.copy()
        data.X.data = np.log(data.X.data + 1) / np.log(base)
    else:
        data.X = np.log(data.X + 1) / np.log(base)
    
    data.uns['log_transformed'] = True
    data.uns['log_base'] = base
    
    return data


def scale_data(
    data: SingleCellData,
    zero_center: bool = True,
    max_value: Optional[float] = None,
) -> SingleCellData:
    """Scale expression data (z-score normalization).
    
    Args:
        data: SingleCellData object (should be log-transformed)
        zero_center: Whether to center data around zero
        max_value: Maximum value after scaling (clips outliers)
        
    Returns:
        Scaled SingleCellData object
    """
    data = data.copy()
    X = data.X.copy()
    
    if sparse.issparse(X):
        X = X.toarray()  # Convert to dense for scaling
    
    if zero_center:
        # Z-score normalization per gene
        X_scaled = zscore(X, axis=0, nan_policy='omit')
        # Handle genes with zero variance
        X_scaled = np.nan_to_num(X_scaled, nan=0.0, posinf=0.0, neginf=0.0)
    else:
        # Scale to unit variance without centering
        std_genes = np.std(X, axis=0)
        std_genes[std_genes == 0] = 1  # Avoid division by zero
        X_scaled = X / std_genes[np.newaxis, :]
    
    # Clip extreme values if specified
    if max_value is not None:
        X_scaled = np.clip(X_scaled, -max_value, max_value)
    
    data.X = X_scaled
    data.uns['scaled'] = True
    
    return data
