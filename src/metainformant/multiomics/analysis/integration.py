"""Multi-omics data integration methods.

This module provides functions for integrating data from multiple omics types
(DNA, RNA, protein, epigenome) using joint dimensionality reduction, correlation
analysis, and other integrative approaches.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy import stats

from metainformant.core.data import validation
from metainformant.core.utils import errors
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Optional scientific dependencies
try:
    from sklearn.cross_decomposition import CCA
    from sklearn.decomposition import NMF, PCA
    from sklearn.preprocessing import StandardScaler

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    StandardScaler = None
    PCA = None
    NMF = None
    CCA = None


class MultiOmicsData:
    """Container for multi-omics data integration.

    This class provides a unified interface for storing and manipulating
    data from multiple omics types with proper alignment and metadata.
    """

    def __init__(
        self,
        data: Dict[str, pd.DataFrame] | None = None,
        sample_ids: List[str] | None = None,
        feature_ids: Dict[str, List[str]] | None = None,
        metadata: pd.DataFrame | None = None,
        # Legacy aliases for individual omics types
        genomics: pd.DataFrame | None = None,
        transcriptomics: pd.DataFrame | None = None,
        proteomics: pd.DataFrame | None = None,
        epigenomics: pd.DataFrame | None = None,
        metabolomics: pd.DataFrame | None = None,
        dna_data: pd.DataFrame | None = None,
        rna_data: pd.DataFrame | None = None,
        protein_data: pd.DataFrame | None = None,
    ):
        """Initialize multi-omics data container.

        Args:
            data: Dictionary mapping omics type to data matrix
            sample_ids: Common sample identifiers
            feature_ids: Feature identifiers for each omics type
            genomics: DNA/genomics data (legacy parameter)
            transcriptomics: RNA data (legacy parameter)
            proteomics: Protein data (legacy parameter)
            epigenomics: Epigenetic data (legacy parameter)
            metabolomics: Metabolomics data (legacy parameter)
            dna_data: Alias for genomics
            rna_data: Alias for transcriptomics
            protein_data: Alias for proteomics
        """
        # Handle legacy parameters
        if data is None:
            data = {}

        # Add legacy parameters to data dict
        # Use explicit is not None checks to avoid DataFrame truth value ambiguity
        legacy_mapping = {
            "genomics": genomics if genomics is not None else dna_data,
            "transcriptomics": transcriptomics if transcriptomics is not None else rna_data,
            "proteomics": proteomics if proteomics is not None else protein_data,
            "epigenomics": epigenomics,
            "metabolomics": metabolomics,
        }
        for key, value in legacy_mapping.items():
            if value is not None and key not in data:
                data[key] = value

        if not data:
            raise ValueError("No omics data provided")

        self.data = data.copy()
        self.sample_ids = sample_ids
        self.feature_ids = feature_ids or {}
        self._metadata = metadata if metadata is not None else pd.DataFrame()

        # Validate data compatibility and align samples
        self._validate_and_align_data()

    def _validate_and_align_data(self) -> None:
        """Validate data compatibility and align samples across omics types."""
        if not self.data:
            raise ValueError("No omics data provided")

        # Get common samples across all omics types (samples are in index/rows)
        common_samples = None
        for omics_type, df in self.data.items():
            samples = set(df.index)
            if common_samples is None:
                common_samples = samples
            else:
                common_samples = common_samples.intersection(samples)

        if common_samples is None or len(common_samples) == 0:
            raise ValueError("No common samples found across omics datasets")

        # Align all data to common samples
        common_samples_sorted = sorted(list(common_samples))
        for omics_type in self.data:
            self.data[omics_type] = self.data[omics_type].loc[common_samples_sorted]

        self.sample_ids = common_samples_sorted

        # Align metadata if present
        if isinstance(self._metadata, pd.DataFrame) and not self._metadata.empty:
            metadata_samples = set(self._metadata.index)
            common_metadata = metadata_samples.intersection(set(common_samples_sorted))
            if common_metadata:
                self._metadata = self._metadata.loc[sorted(list(common_metadata))]

    @property
    def n_samples(self) -> int:
        """Get number of samples common across all layers."""
        return len(self.sample_ids) if self.sample_ids else 0

    @property
    def samples(self) -> List[str]:
        """Get list of common sample IDs."""
        return list(self.sample_ids) if self.sample_ids else []

    @property
    def layer_names(self) -> List[str]:
        """Get list of omics layer names."""
        return list(self.data.keys())

    @property
    def metadata(self) -> Optional[pd.DataFrame]:
        """Get sample metadata DataFrame."""
        if isinstance(self._metadata, pd.DataFrame) and not self._metadata.empty:
            return self._metadata
        return None

    def get_layer(self, layer_name: str) -> pd.DataFrame:
        """Get data for a specific omics layer.

        Args:
            layer_name: Name of the omics layer

        Returns:
            DataFrame with samples in rows, features in columns

        Raises:
            KeyError: If layer not found
        """
        if layer_name not in self.data:
            raise KeyError(f"Layer '{layer_name}' not available. Available layers: {self.layer_names}")
        return self.data[layer_name]

    def get_common_samples(self) -> List[str]:
        """Get samples present in all omics types."""
        return self.samples

    def subset_samples(self, sample_ids: List[str]) -> "MultiOmicsData":
        """Create subset with specified samples."""
        subset_data = {}
        available_samples = [s for s in sample_ids if s in self.samples]
        for omics_type, df in self.data.items():
            subset_data[omics_type] = df.loc[available_samples]

        return MultiOmicsData(data=subset_data, sample_ids=available_samples, feature_ids=self.feature_ids)

    def subset_features(self, feature_dict: Dict[str, List[str]]) -> "MultiOmicsData":
        """Create subset with specified features per layer.

        Args:
            feature_dict: Dictionary mapping layer names to list of features to keep

        Returns:
            New MultiOmicsData with subset features
        """
        subset_data = {}
        for omics_type, df in self.data.items():
            if omics_type in feature_dict:
                # Subset to specified features
                features = [f for f in feature_dict[omics_type] if f in df.columns]
                subset_data[omics_type] = df[features]
            else:
                # Keep all features for this layer
                subset_data[omics_type] = df.copy()

        return MultiOmicsData(data=subset_data, sample_ids=self.sample_ids, feature_ids=self.feature_ids)

    def add_metadata(self, key: str, value: Any) -> None:
        """Add metadata to the dataset."""
        if not isinstance(self._metadata, pd.DataFrame) or self._metadata.empty:
            self._metadata = pd.DataFrame(index=self.samples)
        self._metadata[key] = value

    def get_metadata(self, key: str) -> Any:
        """Get metadata value."""
        if isinstance(self._metadata, pd.DataFrame) and key in self._metadata.columns:
            return self._metadata[key]
        return None


def integrate_omics_data(
    data: Optional[Dict[str, Union[pd.DataFrame, str, Path]]] = None,
    dna_data: pd.DataFrame | None = None,
    rna_data: pd.DataFrame | None = None,
    protein_data: pd.DataFrame | None = None,
    epigenome_data: pd.DataFrame | None = None,
    metabolomics_data: pd.DataFrame | None = None,
    **kwargs,
) -> "MultiOmicsData":
    """Integrate data from multiple omics types.

    Args:
        data: Dictionary mapping omics type to data (DataFrame or file path)
        dna_data: DNA-related data (variants, copy number, etc.)
        rna_data: RNA expression data
        protein_data: Protein abundance data
        epigenome_data: Epigenetic data (methylation, ChIP-seq, etc.)
        metabolomics_data: Metabolomics data
        **kwargs: Additional integration parameters

    Returns:
        MultiOmicsData object with integrated data

    Raises:
        ValueError: If no data provided or incompatible data shapes
    """
    # Handle dict input
    if data is not None:
        # Process dict - load files if needed
        processed_data = {}
        for key, value in data.items():
            if isinstance(value, (str, Path)):
                # Load from file
                path = Path(value)
                if path.suffix == ".csv":
                    processed_data[key] = pd.read_csv(path, index_col=0)
                elif path.suffix == ".parquet":
                    processed_data[key] = pd.read_parquet(path)
                else:
                    raise errors.ValidationError(f"Unsupported file format: {path.suffix}")
            else:
                processed_data[key] = value
        return MultiOmicsData(data=processed_data, **kwargs)

    # Handle legacy individual params
    omics_data = {
        "dna": dna_data,
        "rna": rna_data,
        "protein": protein_data,
        "epigenome": epigenome_data,
        "metabolomics": metabolomics_data,
    }

    # Filter out None values
    available_omics = {k: v for k, v in omics_data.items() if v is not None}

    if not available_omics:
        raise errors.ValidationError("At least one omics dataset must be provided")

    logger.info(f"Integrating {len(available_omics)} omics types: {list(available_omics.keys())}")

    return MultiOmicsData(data=available_omics, **kwargs)


def joint_pca(
    multiomics_data: Union["MultiOmicsData", Dict[str, pd.DataFrame]],
    n_components: int = 50,
    standardize: bool = True,
    layer_weights: Optional[Dict[str, float]] = None,
    **kwargs,
) -> Tuple[np.ndarray, Dict[str, np.ndarray], np.ndarray]:
    """Perform joint PCA across multiple omics datasets.

    Args:
        multiomics_data: MultiOmicsData object or dictionary of omics datasets
        n_components: Number of joint components
        standardize: Whether to standardize the data
        layer_weights: Optional weights for each layer
        **kwargs: Additional PCA parameters

    Returns:
        Tuple of (embeddings, loadings_dict, explained_variance_ratio)

    Raises:
        ImportError: If scikit-learn not available
        ValueError: If data shapes incompatible
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn is required for joint PCA. " "Install with: uv pip install scikit-learn")

    validation.validate_range(n_components, min_val=1, name="n_components")

    # Handle MultiOmicsData input
    if hasattr(multiomics_data, "data"):
        data_dict = multiomics_data.data
    else:
        data_dict = multiomics_data

    logger.info(f"Performing joint PCA with {n_components} components")

    # Apply layer weights if specified
    weighted_data = {}
    for layer_name, df in data_dict.items():
        weight = layer_weights.get(layer_name, 1.0) if layer_weights else 1.0
        weighted_data[layer_name] = df * weight

    # Concatenate all datasets horizontally (features from different omics)
    concatenated_data = pd.concat(list(weighted_data.values()), axis=1)

    # Handle missing values
    concatenated_data = concatenated_data.fillna(concatenated_data.mean())

    # Scale data if requested
    if standardize:
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(concatenated_data)
    else:
        scaled_data = concatenated_data.values

    # Perform PCA
    n_comp = min(n_components, scaled_data.shape[1], scaled_data.shape[0])
    pca = PCA(n_components=n_comp, **kwargs)
    embeddings = pca.fit_transform(scaled_data)

    # Create component loadings for each omics type
    loadings = {}
    feature_start = 0

    for omics_type, df in data_dict.items():
        n_features = df.shape[1]
        loadings[omics_type] = pca.components_[:, feature_start : feature_start + n_features].T
        feature_start += n_features

    variance = pca.explained_variance_ratio_

    logger.info(f"Joint PCA completed: {len(variance)} components explain {np.sum(variance):.1%} variance")

    return embeddings, loadings, variance


def joint_nmf(
    multiomics_data: Union["MultiOmicsData", Dict[str, pd.DataFrame]],
    n_components: int = 50,
    max_iter: int = 200,
    regularization: float = 0.0,
    random_state: Optional[int] = None,
    **kwargs,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Perform joint NMF across multiple omics datasets.

    Args:
        multiomics_data: MultiOmicsData object or dictionary of omics datasets
        n_components: Number of joint components
        max_iter: Maximum iterations for NMF
        regularization: L2 regularization strength
        random_state: Random seed for reproducibility
        **kwargs: Additional NMF parameters

    Returns:
        Tuple of (W matrix, H_dict) where H_dict maps layer names to component loadings

    Raises:
        ImportError: If scikit-learn not available
        ValueError: If data contains negative values
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn is required for joint NMF. " "Install with: uv pip install scikit-learn")

    validation.validate_range(n_components, min_val=1, name="n_components")
    validation.validate_range(max_iter, min_val=10, name="max_iter")

    # Handle MultiOmicsData input
    if hasattr(multiomics_data, "data"):
        data_dict = multiomics_data.data
    else:
        data_dict = multiomics_data

    logger.info(f"Performing joint NMF with {n_components} components")

    # Concatenate all datasets
    concatenated_data = pd.concat(list(data_dict.values()), axis=1)

    # Check for negative values
    if (concatenated_data < 0).any().any():
        raise errors.ValidationError("NMF requires non-negative data. Use joint_pca for data with negative values.")

    # Handle missing values
    concatenated_data = concatenated_data.fillna(concatenated_data.mean())

    # Perform NMF
    n_comp = min(n_components, concatenated_data.shape[1])
    nmf = NMF(
        n_components=n_comp,
        max_iter=max_iter,
        alpha_H=regularization,
        alpha_W=regularization,
        random_state=random_state,
        **kwargs,
    )
    W = nmf.fit_transform(concatenated_data.values)
    H_full = nmf.components_

    # Split H matrix by omics type
    H = {}
    feature_start = 0

    for omics_type, data in data_dict.items():
        n_features = data.shape[1]
        H[omics_type] = H_full[:, feature_start : feature_start + n_features]
        feature_start += n_features

    logger.info(f"Joint NMF completed: reconstruction error = {nmf.reconstruction_err_:.4f}")
    return W, H


def canonical_correlation(
    multiomics_data: Union["MultiOmicsData", Dict[str, pd.DataFrame]],
    layers: Optional[Tuple[str, str]] = None,
    n_components: int = 10,
    **kwargs,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Perform canonical correlation analysis between omics datasets.

    Args:
        multiomics_data: MultiOmicsData object or dictionary of omics datasets
        layers: Tuple of (layer1, layer2) to correlate. Required if more than 2 layers.
        n_components: Number of canonical components
        **kwargs: Additional CCA parameters

    Returns:
        Tuple of (X_c, Y_c, X_weights, Y_weights, correlations)

    Raises:
        ImportError: If scikit-learn not available
        ValueError: If layers not specified and not exactly 2 datasets provided
    """
    if not HAS_SKLEARN:
        raise ImportError(
            "scikit-learn is required for canonical correlation analysis. " "Install with: uv pip install scikit-learn"
        )

    # Handle MultiOmicsData input
    if hasattr(multiomics_data, "data"):
        data_dict = multiomics_data.data
    else:
        data_dict = multiomics_data

    # Determine layers to use
    if layers is None:
        if len(data_dict) != 2:
            raise errors.ValidationError("CCA requires exactly 2 omics datasets or layers tuple specified")
        omics_types = list(data_dict.keys())
    else:
        omics_types = list(layers)
        # Check layers exist
        for layer in omics_types:
            if layer not in data_dict:
                raise ValueError(f"Layer {layer} not found in data. Available: {list(data_dict.keys())}")

    validation.validate_range(n_components, min_val=1, name="n_components")

    logger.info(f"Performing CCA between {omics_types[0]} and {omics_types[1]}")

    # Get the two datasets
    X = data_dict[omics_types[0]].values
    Y = data_dict[omics_types[1]].values

    # Handle missing values
    X = np.nan_to_num(X, nan=np.nanmean(X))
    Y = np.nan_to_num(Y, nan=np.nanmean(Y))

    # Scale data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    Y_scaled = scaler.fit_transform(Y)

    # Perform CCA
    n_comp = min(n_components, X_scaled.shape[1], Y_scaled.shape[1])
    cca = CCA(n_components=n_comp, **kwargs)
    X_c, Y_c = cca.fit_transform(X_scaled, Y_scaled)

    # Compute canonical correlations for each component
    correlations = []
    for i in range(X_c.shape[1]):
        corr = np.corrcoef(X_c[:, i], Y_c[:, i])[0, 1]
        correlations.append(abs(corr))
    correlations = np.array(correlations)

    logger.info(f"CCA completed: {len(correlations)} components")
    return X_c, Y_c, cca.x_weights_, cca.y_weights_, correlations


def from_dna_variants(vcf_data: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """Convert DNA variant data for multi-omics integration.

    Args:
        vcf_data: VCF-style variant data
        **kwargs: Conversion parameters

    Returns:
        Processed DNA data suitable for integration
    """
    logger.info("Converting DNA variant data for integration")

    # Simple conversion - in practice this would handle VCF format properly
    # Convert genotypes to numeric (0, 1, 2 for homozygous ref, het, homozygous alt)
    processed_data = vcf_data.copy()

    # This is a placeholder - real implementation would parse VCF properly
    logger.warning("DNA variant conversion is simplified - implement proper VCF parsing for production use")

    return processed_data


def from_rna_expression(expression_data: pd.DataFrame, normalize: bool = True, **kwargs) -> pd.DataFrame:
    """Convert RNA expression data for multi-omics integration.

    Args:
        expression_data: RNA-seq or microarray expression data
        normalize: Whether to normalize the data
        **kwargs: Conversion parameters

    Returns:
        Processed RNA data suitable for integration
    """
    logger.info("Converting RNA expression data for integration")

    processed_data = expression_data.copy()

    if normalize:
        # Simple normalization - log transform and z-score
        if (processed_data > 0).all().all():
            processed_data = np.log1p(processed_data)

        scaler = StandardScaler()
        processed_data = pd.DataFrame(
            scaler.fit_transform(processed_data), index=processed_data.index, columns=processed_data.columns
        )

    return processed_data


def from_protein_abundance(protein_data: pd.DataFrame, normalize: bool = True, **kwargs) -> pd.DataFrame:
    """Convert protein abundance data for multi-omics integration.

    Args:
        protein_data: Protein abundance measurements
        normalize: Whether to normalize the data
        **kwargs: Conversion parameters

    Returns:
        Processed protein data suitable for integration
    """
    logger.info("Converting protein abundance data for integration")

    processed_data = protein_data.copy()

    if normalize:
        # Z-score normalization
        scaler = StandardScaler()
        processed_data = pd.DataFrame(
            scaler.fit_transform(processed_data), index=processed_data.index, columns=processed_data.columns
        )

    return processed_data


def from_epigenome_data(epigenome_data: pd.DataFrame, data_type: str = "methylation", **kwargs) -> pd.DataFrame:
    """Convert epigenome data for multi-omics integration.

    Args:
        epigenome_data: Epigenetic data (methylation, ChIP-seq, etc.)
        data_type: Type of epigenetic data
        **kwargs: Conversion parameters

    Returns:
        Processed epigenome data suitable for integration
    """
    logger.info(f"Converting {data_type} epigenome data for integration")

    processed_data = epigenome_data.copy()

    if data_type == "methylation":
        # Methylation data often needs beta-value transformation
        # Assume data is already in appropriate format
        pass
    elif data_type in ["chipseq", "atacseq"]:
        # Peak data - convert to binary or intensity values
        pass

    return processed_data


def from_metabolomics(metabolomics_data: pd.DataFrame, normalize: bool = True, **kwargs) -> pd.DataFrame:
    """Convert metabolomics data for multi-omics integration.

    Args:
        metabolomics_data: Metabolomics measurements
        normalize: Whether to normalize the data
        **kwargs: Conversion parameters

    Returns:
        Processed metabolomics data suitable for integration
    """
    logger.info("Converting metabolomics data for integration")

    processed_data = metabolomics_data.copy()

    if normalize:
        # Metabolomics data often has large dynamic range
        # Log transform if positive, then z-score
        if (processed_data > 0).all().all():
            processed_data = np.log(processed_data)

        scaler = StandardScaler()
        processed_data = pd.DataFrame(
            scaler.fit_transform(processed_data), index=processed_data.index, columns=processed_data.columns
        )

    return processed_data


def _integrate_by_correlation(aligned_data: Dict[str, pd.DataFrame], **kwargs) -> Dict[str, Any]:
    """Integrate omics data by computing cross-omics correlations."""
    logger.info("Integrating by correlation analysis")

    results = {}

    # Compute pairwise correlations between all omics types
    omics_types = list(aligned_data.keys())
    correlation_matrices = {}

    for i, omics1 in enumerate(omics_types):
        for j, omics2 in enumerate(omics_types):
            if i < j:  # Upper triangle only
                data1 = aligned_data[omics1].values
                data2 = aligned_data[omics2].values

                # Compute correlation matrix
                corr_matrix = np.corrcoef(data1.T, data2.T)[: data1.shape[1], data1.shape[1] :]

                correlation_matrices[f"{omics1}_{omics2}"] = {
                    "correlation_matrix": corr_matrix,
                    "mean_correlation": np.mean(np.abs(corr_matrix)),
                    "max_correlation": np.max(np.abs(corr_matrix)),
                    "omics1_features": aligned_data[omics1].columns.tolist(),
                    "omics2_features": aligned_data[omics2].columns.tolist(),
                }

    results["correlation_matrices"] = correlation_matrices

    # Find most correlated feature pairs
    top_correlations = []
    for pair_name, corr_data in correlation_matrices.items():
        corr_matrix = corr_data["correlation_matrix"]
        features1 = corr_data["omics1_features"]
        features2 = corr_data["omics2_features"]

        # Find top correlations
        n_top = min(100, corr_matrix.size)  # Top 100 or all if fewer
        flat_indices = np.argsort(np.abs(corr_matrix).flatten())[-n_top:]

        for idx in flat_indices:
            i, j = np.unravel_index(idx, corr_matrix.shape)
            top_correlations.append(
                {
                    "omics_pair": pair_name,
                    "feature1": features1[i],
                    "feature2": features2[j],
                    "correlation": corr_matrix[i, j],
                }
            )

    results["top_correlations"] = sorted(top_correlations, key=lambda x: abs(x["correlation"]), reverse=True)

    return results


def _validate_omics_data_compatibility(omics_data: Dict[str, pd.DataFrame]) -> None:
    """Validate that omics datasets are compatible for integration."""
    if not omics_data:
        raise errors.ValidationError("No omics data provided")

    # Check that all datasets have samples (rows)
    for omics_type, data in omics_data.items():
        if data.shape[0] == 0:
            raise errors.ValidationError(f"{omics_type} data has no samples")

        if data.shape[1] == 0:
            raise errors.ValidationError(f"{omics_type} data has no features")

    # Check for reasonable sample overlap (at least some samples should be shared)
    sample_sets = []
    for data in omics_data.values():
        if hasattr(data, "index"):
            sample_sets.append(set(data.index))
        else:
            sample_sets.append(set(range(data.shape[0])))

    intersection = set.intersection(*sample_sets)
    if len(intersection) == 0:
        raise errors.ValidationError("No common samples found across all omics datasets")

    union = set.union(*sample_sets)
    overlap_fraction = len(intersection) / len(union)

    if overlap_fraction < 0.1:  # Less than 10% overlap
        logger.warning(f"Low sample overlap across datasets: {overlap_fraction:.1%}")


def compute_multiomics_similarity(omics_data: Dict[str, pd.DataFrame], method: str = "correlation") -> np.ndarray:
    """Compute similarity matrix across all samples using multi-omics data.

    Args:
        omics_data: Dictionary of aligned omics datasets
        method: Similarity computation method

    Returns:
        Similarity matrix between samples

    Raises:
        ValueError: If method not supported
        ImportError: If scikit-learn required but not available
    """
    if method not in ["correlation", "euclidean", "cosine"]:
        raise errors.ValidationError(f"Unsupported similarity method: {method}")

    if method == "cosine" and not HAS_SKLEARN:
        raise ImportError(
            "scikit-learn is required for cosine similarity. " "Install with: uv pip install scikit-learn"
        )

    logger.info(f"Computing multi-omics similarity using {method}")

    # Concatenate all omics data
    concatenated = np.concatenate([data.values for data in omics_data.values()], axis=1)

    # Handle missing values
    concatenated = np.nan_to_num(concatenated, nan=np.nanmean(concatenated))

    if method == "correlation":
        # Correlation-based similarity
        similarity = np.corrcoef(concatenated)
    elif method == "euclidean":
        # Convert distance to similarity
        from scipy.spatial.distance import pdist, squareform

        distances = squareform(pdist(concatenated, metric="euclidean"))
        # Convert to similarity (inverse distance)
        similarity = 1 / (1 + distances)
    elif method == "cosine":
        from sklearn.metrics.pairwise import cosine_similarity

        similarity = cosine_similarity(concatenated)

    return similarity


def find_multiomics_modules(omics_data: Dict[str, pd.DataFrame], n_modules: int = 10, **kwargs) -> Dict[str, Any]:
    """Identify multi-omics modules (co-regulated features across omics types).

    Args:
        omics_data: Dictionary of aligned omics datasets
        n_modules: Number of modules to identify
        **kwargs: Additional parameters

    Returns:
        Dictionary with module information

    Raises:
        ValueError: If parameters invalid
    """
    validation.validate_range(n_modules, min_val=2, name="n_modules")

    logger.info(f"Finding multi-omics modules: {n_modules} modules")

    # Use joint NMF to find modules
    nmf_results = joint_nmf(omics_data, n_components=n_modules, **kwargs)

    # Interpret modules
    modules = {}

    for i in range(n_modules):
        module_features = {}

        for omics_type, components in nmf_results["omics_components"].items():
            # Find features with high loading in this component
            loadings = components[i, :]
            top_features = np.argsort(loadings)[-10:]  # Top 10 features

            module_features[omics_type] = {
                "features": [omics_data[omics_type].columns[j] for j in top_features],
                "loadings": loadings[top_features].tolist(),
            }

        modules[f"module_{i+1}"] = {
            "features": module_features,
            "sample_weights": nmf_results["W_matrix"][:, i].tolist(),
        }

    results = {
        "modules": modules,
        "nmf_results": nmf_results,
        "n_modules": n_modules,
    }

    logger.info(f"Multi-omics module detection completed: {n_modules} modules found")
    return results
