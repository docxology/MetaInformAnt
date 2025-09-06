"""Core multi-omics data integration methods."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd


@dataclass
class MultiOmicsData:
    """Container for multi-omics datasets with sample alignment."""

    genomics: Optional[pd.DataFrame] = None
    transcriptomics: Optional[pd.DataFrame] = None
    proteomics: Optional[pd.DataFrame] = None
    metabolomics: Optional[pd.DataFrame] = None
    epigenomics: Optional[pd.DataFrame] = None
    metadata: Optional[pd.DataFrame] = None

    def __post_init__(self):
        """Validate and align samples across omics layers."""
        self.omics_layers = {}

        # Store non-None omics data
        for name in ["genomics", "transcriptomics", "proteomics", "metabolomics", "epigenomics"]:
            data = getattr(self, name)
            if data is not None:
                if not isinstance(data, pd.DataFrame):
                    raise TypeError(f"{name} must be pandas DataFrame")
                self.omics_layers[name] = data

        if not self.omics_layers:
            raise ValueError("At least one omics layer must be provided")

        # Align samples
        self._align_samples()

    def _align_samples(self):
        """Align samples across all omics layers."""
        if not self.omics_layers:
            return

        # Find common samples
        sample_sets = [set(df.index) for df in self.omics_layers.values()]
        common_samples = sample_sets[0].intersection(*sample_sets[1:])

        if len(common_samples) == 0:
            raise ValueError("No common samples found across omics layers")

        if len(common_samples) < min(len(s) for s in sample_sets):
            warnings.warn(f"Only {len(common_samples)} samples are common across all omics layers")

        # Reorder all datasets to common samples
        common_samples = sorted(common_samples)
        for name, data in self.omics_layers.items():
            self.omics_layers[name] = data.loc[common_samples]
            setattr(self, name, self.omics_layers[name])

        # Align metadata if provided
        if self.metadata is not None:
            if set(self.metadata.index).intersection(set(common_samples)):
                self.metadata = self.metadata.loc[self.metadata.index.intersection(common_samples)]

    @property
    def samples(self) -> List[str]:
        """Get list of aligned samples."""
        if self.omics_layers:
            return list(next(iter(self.omics_layers.values())).index)
        return []

    @property
    def n_samples(self) -> int:
        """Number of aligned samples."""
        return len(self.samples)

    @property
    def layer_names(self) -> List[str]:
        """Names of available omics layers."""
        return list(self.omics_layers.keys())

    def get_layer(self, layer_name: str) -> pd.DataFrame:
        """Get specific omics layer."""
        if layer_name not in self.omics_layers:
            raise KeyError(f"Layer '{layer_name}' not available. Available: {self.layer_names}")
        return self.omics_layers[layer_name]

    def subset_samples(self, sample_list: List[str]) -> "MultiOmicsData":
        """Create subset with specified samples."""
        new_data = MultiOmicsData()

        for layer_name, data in self.omics_layers.items():
            available_samples = [s for s in sample_list if s in data.index]
            if available_samples:
                setattr(new_data, layer_name, data.loc[available_samples])

        if self.metadata is not None:
            available_samples = [s for s in sample_list if s in self.metadata.index]
            if available_samples:
                new_data.metadata = self.metadata.loc[available_samples]

        new_data._align_samples()
        return new_data

    def subset_features(self, feature_dict: Dict[str, List[str]]) -> "MultiOmicsData":
        """Create subset with specified features per layer."""
        new_data = MultiOmicsData()

        for layer_name, data in self.omics_layers.items():
            if layer_name in feature_dict:
                available_features = [f for f in feature_dict[layer_name] if f in data.columns]
                if available_features:
                    setattr(new_data, layer_name, data[available_features])
            else:
                setattr(new_data, layer_name, data.copy())

        if self.metadata is not None:
            new_data.metadata = self.metadata.copy()

        new_data._align_samples()
        return new_data


def integrate_omics_data(
    data_dict: Dict[str, Union[pd.DataFrame, str, Path]],
    sample_mapping: Optional[Dict[str, str]] = None,
    feature_mapping: Optional[Dict[str, Dict[str, str]]] = None,
    metadata: Optional[Union[pd.DataFrame, str, Path]] = None,
) -> MultiOmicsData:
    """Load and integrate multiple omics datasets.

    Args:
        data_dict: Dictionary mapping omics type to DataFrame or file path
        sample_mapping: Mapping of sample IDs across datasets
        feature_mapping: Mapping of feature IDs within each dataset
        metadata: Sample metadata DataFrame or file path

    Returns:
        Integrated MultiOmicsData object
    """
    omics_data = {}

    # Load each omics dataset
    for omics_type, data in data_dict.items():
        if isinstance(data, (str, Path)):
            # Load from file
            data_path = Path(data)
            if data_path.suffix.lower() in [".csv"]:
                df = pd.read_csv(data_path, index_col=0)
            elif data_path.suffix.lower() in [".tsv", ".txt"]:
                df = pd.read_csv(data_path, sep="\t", index_col=0)
            elif data_path.suffix.lower() in [".xlsx", ".xls"]:
                df = pd.read_excel(data_path, index_col=0)
            else:
                raise ValueError(f"Unsupported file format: {data_path.suffix}")
        else:
            df = data.copy()

        # Apply sample mapping if provided
        if sample_mapping and omics_type in sample_mapping:
            # This would implement sample ID mapping logic
            pass

        # Apply feature mapping if provided
        if feature_mapping and omics_type in feature_mapping:
            feature_map = feature_mapping[omics_type]
            df = df.rename(columns=feature_map)

        omics_data[omics_type] = df

    # Load metadata if provided
    metadata_df = None
    if metadata is not None:
        if isinstance(metadata, (str, Path)):
            metadata_path = Path(metadata)
            if metadata_path.suffix.lower() == ".csv":
                metadata_df = pd.read_csv(metadata_path, index_col=0)
            else:
                metadata_df = pd.read_csv(metadata_path, sep="\t", index_col=0)
        else:
            metadata_df = metadata.copy()

    # Create MultiOmicsData object
    kwargs = {
        k: v
        for k, v in omics_data.items()
        if k in ["genomics", "transcriptomics", "proteomics", "metabolomics", "epigenomics"]
    }
    kwargs["metadata"] = metadata_df

    return MultiOmicsData(**kwargs)


def joint_pca(
    omics_data: MultiOmicsData,
    n_components: int = 50,
    layer_weights: Optional[Dict[str, float]] = None,
    standardize: bool = True,
) -> Tuple[np.ndarray, Dict[str, np.ndarray], np.ndarray]:
    """Joint Principal Component Analysis across omics layers.

    Args:
        omics_data: Multi-omics data object
        n_components: Number of principal components
        layer_weights: Relative weights for each omics layer
        standardize: Whether to standardize features

    Returns:
        Tuple of (joint_embeddings, layer_loadings, explained_variance)
    """
    if layer_weights is None:
        layer_weights = {layer: 1.0 for layer in omics_data.layer_names}

    # Standardize and concatenate data
    concatenated_data = []
    layer_loadings = {}
    feature_ranges = {}

    current_idx = 0
    for layer_name in omics_data.layer_names:
        data = omics_data.get_layer(layer_name).values

        if standardize:
            # Z-score standardization
            mean = np.mean(data, axis=0, keepdims=True)
            std = np.std(data, axis=0, keepdims=True) + 1e-8
            data = (data - mean) / std

        # Apply layer weight
        weight = layer_weights.get(layer_name, 1.0)
        data = data * np.sqrt(weight)

        concatenated_data.append(data)
        n_features = data.shape[1]
        feature_ranges[layer_name] = (current_idx, current_idx + n_features)
        current_idx += n_features

    # Concatenate all layers
    X_concat = np.hstack(concatenated_data)

    # Perform PCA
    n_components = min(n_components, min(X_concat.shape) - 1)

    # Compute covariance matrix
    X_centered = X_concat - np.mean(X_concat, axis=0, keepdims=True)
    cov_matrix = np.cov(X_centered.T)

    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

    # Sort by eigenvalues (descending)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx][:n_components]
    eigenvectors = eigenvectors[:, idx][:, :n_components]

    # Project data
    joint_embeddings = X_centered @ eigenvectors

    # Extract loadings for each layer
    for layer_name, (start, end) in feature_ranges.items():
        layer_loadings[layer_name] = eigenvectors[start:end, :]

    # Explained variance ratio
    explained_variance = eigenvalues / np.sum(eigenvalues)

    return joint_embeddings, layer_loadings, explained_variance


def joint_nmf(
    omics_data: MultiOmicsData,
    n_components: int = 20,
    max_iter: int = 200,
    regularization: float = 0.01,
    random_state: Optional[int] = None,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Joint Non-negative Matrix Factorization across omics layers.

    Args:
        omics_data: Multi-omics data object
        n_components: Number of components
        max_iter: Maximum iterations
        regularization: L2 regularization parameter
        random_state: Random seed

    Returns:
        Tuple of (sample_factors, layer_feature_factors)
    """
    if random_state is not None:
        np.random.seed(random_state)

    # Ensure non-negative data
    layer_data = {}
    for layer_name in omics_data.layer_names:
        data = omics_data.get_layer(layer_name).values
        # Shift to non-negative if needed
        data = data - np.min(data) + 0.1
        layer_data[layer_name] = data

    n_samples = omics_data.n_samples

    # Initialize factors
    W = np.random.rand(n_samples, n_components)  # Sample factors (shared)
    H = {}  # Feature factors per layer

    for layer_name, data in layer_data.items():
        n_features = data.shape[1]
        H[layer_name] = np.random.rand(n_components, n_features)

    # Alternating optimization
    for iteration in range(max_iter):
        # Update H (feature factors) for each layer
        for layer_name, data in layer_data.items():
            numerator = W.T @ data + regularization * H[layer_name]
            denominator = W.T @ W @ H[layer_name] + regularization + 1e-10
            H[layer_name] *= numerator / denominator

        # Update W (sample factors) - shared across all layers
        W_numerator = np.zeros_like(W)
        W_denominator = np.zeros_like(W)

        for layer_name, data in layer_data.items():
            W_numerator += data @ H[layer_name].T
            W_denominator += W @ H[layer_name] @ H[layer_name].T

        W_numerator += regularization * W
        W_denominator += regularization + 1e-10
        W *= W_numerator / W_denominator

        # Optional: compute reconstruction error for convergence checking
        if iteration % 50 == 0:
            total_error = 0
            for layer_name, data in layer_data.items():
                reconstruction = W @ H[layer_name]
                error = np.sum((data - reconstruction) ** 2)
                total_error += error

            if iteration > 0 and total_error < 1e-6:
                break

    return W, H


def canonical_correlation(
    omics_data: MultiOmicsData, layer_pair: Tuple[str, str], n_components: int = 10, regularization: float = 0.01
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Canonical Correlation Analysis between two omics layers.

    Args:
        omics_data: Multi-omics data object
        layer_pair: Tuple of two layer names to analyze
        n_components: Number of canonical components
        regularization: Regularization parameter

    Returns:
        Tuple of (X_canonical, Y_canonical, X_weights, Y_weights, correlations)
    """
    layer1, layer2 = layer_pair

    if layer1 not in omics_data.layer_names:
        raise ValueError(f"Layer {layer1} not found")
    if layer2 not in omics_data.layer_names:
        raise ValueError(f"Layer {layer2} not found")

    X = omics_data.get_layer(layer1).values
    Y = omics_data.get_layer(layer2).values

    # Standardize data
    X = (X - np.mean(X, axis=0)) / (np.std(X, axis=0) + 1e-8)
    Y = (Y - np.mean(Y, axis=0)) / (np.std(Y, axis=0) + 1e-8)

    n_samples, n_features_x = X.shape
    n_features_y = Y.shape[1]

    # Compute covariance matrices
    C_xx = np.cov(X.T) + regularization * np.eye(n_features_x)
    C_yy = np.cov(Y.T) + regularization * np.eye(n_features_y)
    C_xy = np.cov(X.T, Y.T)[:n_features_x, n_features_x:]

    # Solve generalized eigenvalue problem
    try:
        # C_xx^{-1} C_xy C_yy^{-1} C_yx
        C_xx_inv = np.linalg.inv(C_xx)
        C_yy_inv = np.linalg.inv(C_yy)

        M = C_xx_inv @ C_xy @ C_yy_inv @ C_xy.T
        eigenvalues, eigenvectors_x = np.linalg.eigh(M)

        # Sort by eigenvalues (descending)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors_x = eigenvectors_x[:, idx]

        # Get Y eigenvectors
        eigenvectors_y = C_yy_inv @ C_xy.T @ eigenvectors_x

        # Normalize
        for i in range(eigenvectors_y.shape[1]):
            norm_y = np.sqrt(eigenvectors_y[:, i].T @ C_yy @ eigenvectors_y[:, i])
            if norm_y > 1e-8:
                eigenvectors_y[:, i] /= norm_y

        # Take top components
        n_components = min(n_components, len(eigenvalues))
        X_weights = eigenvectors_x[:, :n_components]
        Y_weights = eigenvectors_y[:, :n_components]

        # Compute canonical variables
        X_canonical = X @ X_weights
        Y_canonical = Y @ Y_weights

        # Compute correlations
        correlations = np.sqrt(eigenvalues[:n_components])

        return X_canonical, Y_canonical, X_weights, Y_weights, correlations

    except np.linalg.LinAlgError:
        # Fallback to simpler approach
        warnings.warn("CCA failed, using SVD-based approach")

        # SVD of cross-covariance
        U, s, Vt = np.linalg.svd(C_xy, full_matrices=False)

        n_components = min(n_components, len(s))
        X_weights = U[:, :n_components]
        Y_weights = Vt[:n_components, :].T

        X_canonical = X @ X_weights
        Y_canonical = Y @ Y_weights

        correlations = s[:n_components]

        return X_canonical, Y_canonical, X_weights, Y_weights, correlations
