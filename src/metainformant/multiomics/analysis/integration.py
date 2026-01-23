"""Multi-omics data integration methods.

This module provides functions for integrating data from multiple omics types
(DNA, RNA, protein, epigenome) using joint dimensionality reduction, correlation
analysis, and other integrative approaches.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Any, Tuple, Union
import numpy as np
import pandas as pd
from scipy import stats

from metainformant.core import logging, errors, validation

logger = logging.get_logger(__name__)

# Optional scientific dependencies
try:
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA, NMF
    from sklearn.cross_decomposition import CCA

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
        data: Dict[str, pd.DataFrame],
        sample_ids: List[str] | None = None,
        feature_ids: Dict[str, List[str]] | None = None,
    ):
        """Initialize multi-omics data container.

        Args:
            data: Dictionary mapping omics type to data matrix
            sample_ids: Common sample identifiers
            feature_ids: Feature identifiers for each omics type
        """
        self.data = data.copy()
        self.sample_ids = sample_ids
        self.feature_ids = feature_ids or {}
        self.metadata = {}

        # Validate data compatibility
        self._validate_data()

    def _validate_data(self) -> None:
        """Validate data compatibility across omics types."""
        if not self.data:
            raise ValueError("No omics data provided")

        # Check sample alignment
        sample_counts = {}
        for omics_type, df in self.data.items():
            if self.sample_ids is None:
                sample_counts[omics_type] = len(df.columns)
            else:
                if len(df.columns) != len(self.sample_ids):
                    raise ValueError(f"Sample count mismatch for {omics_type}")

        if self.sample_ids is None and len(set(sample_counts.values())) > 1:
            logger.warning("Different sample counts across omics types - ensure proper alignment")

    def get_common_samples(self) -> List[str]:
        """Get samples present in all omics types."""
        if not self.data:
            return []

        common_samples = None
        for df in self.data.values():
            samples = set(df.columns)
            if common_samples is None:
                common_samples = samples
            else:
                common_samples = common_samples.intersection(samples)

        return sorted(list(common_samples)) if common_samples else []

    def subset_samples(self, sample_ids: List[str]) -> "MultiOmicsData":
        """Create subset with specified samples."""
        subset_data = {}
        for omics_type, df in self.data.items():
            available_samples = [s for s in sample_ids if s in df.columns]
            if available_samples:
                subset_data[omics_type] = df[available_samples]

        return MultiOmicsData(subset_data, sample_ids=available_samples, feature_ids=self.feature_ids)

    def add_metadata(self, key: str, value: Any) -> None:
        """Add metadata to the dataset."""
        self.metadata[key] = value

    def get_metadata(self, key: str) -> Any:
        """Get metadata value."""
        return self.metadata.get(key)


def integrate_omics_data(
    dna_data: pd.DataFrame | None = None,
    rna_data: pd.DataFrame | None = None,
    protein_data: pd.DataFrame | None = None,
    epigenome_data: pd.DataFrame | None = None,
    metabolomics_data: pd.DataFrame | None = None,
    **kwargs,
) -> Dict[str, Any]:
    """Integrate data from multiple omics types.

    Args:
        dna_data: DNA-related data (variants, copy number, etc.)
        rna_data: RNA expression data
        protein_data: Protein abundance data
        epigenome_data: Epigenetic data (methylation, ChIP-seq, etc.)
        metabolomics_data: Metabolomics data
        **kwargs: Additional integration parameters

    Returns:
        Dictionary containing integrated results and metadata

    Raises:
        ValueError: If no data provided or incompatible data shapes
    """
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

    # Validate data compatibility
    _validate_omics_data_compatibility(available_omics)

    # Store original shapes and metadata
    metadata = {
        "n_omics_types": len(available_omics),
        "omics_types": list(available_omics.keys()),
        "original_shapes": {k: v.shape for k, v in available_omics.items()},
        "n_samples": None,
        "n_features": {k: v.shape[1] for k, v in available_omics.items()},
    }

    # Determine common samples
    sample_intersection = None
    for data in available_omics.values():
        current_samples = set(data.index if hasattr(data, "index") else range(data.shape[0]))
        if sample_intersection is None:
            sample_intersection = current_samples
        else:
            sample_intersection = sample_intersection.intersection(current_samples)

    if not sample_intersection:
        raise errors.ValidationError("No common samples found across omics datasets")

    metadata["n_samples"] = len(sample_intersection)
    metadata["common_samples"] = sorted(list(sample_intersection))

    # Align data to common samples
    aligned_data = {}
    for omics_type, data in available_omics.items():
        if hasattr(data, "index"):
            # pandas DataFrame
            aligned_data[omics_type] = data.loc[list(sample_intersection)]
        else:
            # numpy array - assume samples are in same order
            sample_indices = [
                i
                for i, sample in enumerate(data.index if hasattr(data, "index") else range(data.shape[0]))
                if sample in sample_intersection
            ]
            aligned_data[omics_type] = data[sample_indices]

    # Perform multi-omics integration
    integration_method = kwargs.get("method", "correlation")
    integrated_results = {}

    if integration_method == "correlation":
        integrated_results = _integrate_by_correlation(aligned_data, **kwargs)
    elif integration_method == "joint_pca":
        integrated_results = joint_pca(aligned_data, **kwargs)
    elif integration_method == "joint_nmf":
        integrated_results = joint_nmf(aligned_data, **kwargs)
    elif integration_method == "cca":
        integrated_results = canonical_correlation(aligned_data, **kwargs)
    else:
        raise errors.ValidationError(f"Unsupported integration method: {integration_method}")

    # Add metadata to results
    integrated_results["metadata"] = metadata
    integrated_results["integration_method"] = integration_method

    logger.info(f"Multi-omics integration completed using {integration_method}")
    return integrated_results


def joint_pca(multiomics_data: Dict[str, pd.DataFrame], n_components: int = 50, **kwargs) -> Dict[str, Any]:
    """Perform joint PCA across multiple omics datasets.

    Args:
        multiomics_data: Dictionary of omics datasets
        n_components: Number of joint components
        **kwargs: Additional PCA parameters

    Returns:
        Dictionary with joint PCA results

    Raises:
        ImportError: If scikit-learn not available
        ValueError: If data shapes incompatible
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn is required for joint PCA. " "Install with: uv pip install scikit-learn")

    validation.validate_range(n_components, min_val=2, name="n_components")

    logger.info(f"Performing joint PCA with {n_components} components")

    # Concatenate all datasets horizontally (features from different omics)
    concatenated_data = pd.concat(list(multiomics_data.values()), axis=1)

    # Handle missing values
    concatenated_data = concatenated_data.fillna(concatenated_data.mean())

    # Scale data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(concatenated_data)

    # Perform PCA
    pca = PCA(n_components=min(n_components, scaled_data.shape[1]), **kwargs)
    pca_scores = pca.fit_transform(scaled_data)

    # Create component loadings for each omics type
    loadings = {}
    feature_start = 0

    for omics_type, data in multiomics_data.items():
        n_features = data.shape[1]
        loadings[omics_type] = pca.components_[:, feature_start : feature_start + n_features].T
        feature_start += n_features

    results = {
        "pca_scores": pca_scores,
        "pca_loadings": pca.components_,
        "omics_loadings": loadings,
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "cumulative_explained_variance": np.cumsum(pca.explained_variance_ratio_),
        "singular_values": pca.singular_values_,
    }

    logger.info(
        f"Joint PCA completed: {len(results['explained_variance_ratio'])} components explain {results['cumulative_explained_variance'][-1]:.1%} variance"
    )
    return results


def joint_nmf(
    multiomics_data: Dict[str, pd.DataFrame], n_components: int = 50, max_iter: int = 200, **kwargs
) -> Dict[str, Any]:
    """Perform joint NMF across multiple omics datasets.

    Args:
        multiomics_data: Dictionary of omics datasets
        n_components: Number of joint components
        max_iter: Maximum iterations for NMF
        **kwargs: Additional NMF parameters

    Returns:
        Dictionary with joint NMF results

    Raises:
        ImportError: If scikit-learn not available
        ValueError: If data contains negative values
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn is required for joint NMF. " "Install with: uv pip install scikit-learn")

    validation.validate_range(n_components, min_val=2, name="n_components")
    validation.validate_range(max_iter, min_val=10, name="max_iter")

    logger.info(f"Performing joint NMF with {n_components} components")

    # Concatenate all datasets
    concatenated_data = pd.concat(list(multiomics_data.values()), axis=1)

    # Check for negative values
    if (concatenated_data < 0).any().any():
        raise errors.ValidationError("NMF requires non-negative data. Use joint_pca for data with negative values.")

    # Handle missing values
    concatenated_data = concatenated_data.fillna(concatenated_data.mean())

    # Perform NMF
    nmf = NMF(n_components=min(n_components, concatenated_data.shape[1]), max_iter=max_iter, **kwargs)
    W = nmf.fit_transform(concatenated_data.values)
    H = nmf.components_

    # Split H matrix by omics type
    omics_components = {}
    feature_start = 0

    for omics_type, data in multiomics_data.items():
        n_features = data.shape[1]
        omics_components[omics_type] = H[:, feature_start : feature_start + n_features]
        feature_start += n_features

    results = {
        "W_matrix": W,  # Sample factors
        "H_matrix": H,  # Feature factors
        "omics_components": omics_components,
        "reconstruction_error": nmf.reconstruction_err_,
        "n_iter": nmf.n_iter_,
    }

    logger.info(f"Joint NMF completed: reconstruction error = {results['reconstruction_error']:.4f}")
    return results


def canonical_correlation(multiomics_data: Dict[str, pd.DataFrame], n_components: int = 10, **kwargs) -> Dict[str, Any]:
    """Perform canonical correlation analysis between omics datasets.

    Args:
        multiomics_data: Dictionary of omics datasets (requires exactly 2 datasets)
        n_components: Number of canonical components
        **kwargs: Additional CCA parameters

    Returns:
        Dictionary with CCA results

    Raises:
        ImportError: If scikit-learn not available
        ValueError: If not exactly 2 datasets provided
    """
    if not HAS_SKLEARN:
        raise ImportError(
            "scikit-learn is required for canonical correlation analysis. " "Install with: uv pip install scikit-learn"
        )

    if len(multiomics_data) != 2:
        raise errors.ValidationError("CCA requires exactly 2 omics datasets")

    validation.validate_range(n_components, min_val=1, name="n_components")

    omics_types = list(multiomics_data.keys())
    logger.info(f"Performing CCA between {omics_types[0]} and {omics_types[1]}")

    # Get the two datasets
    X = multiomics_data[omics_types[0]].values
    Y = multiomics_data[omics_types[1]].values

    # Handle missing values
    X = np.nan_to_num(X, nan=np.nanmean(X))
    Y = np.nan_to_num(Y, nan=np.nanmean(Y))

    # Scale data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    Y_scaled = scaler.fit_transform(Y)

    # Perform CCA
    cca = CCA(n_components=min(n_components, X_scaled.shape[1], Y_scaled.shape[1]), **kwargs)
    X_c, Y_c = cca.fit_transform(X_scaled, Y_scaled)

    results = {
        "cca_scores_X": X_c,
        "cca_scores_Y": Y_c,
        "cca_coefficients_X": cca.x_weights_,
        "cca_coefficients_Y": cca.y_weights_,
        "canonical_correlations": cca.score(X_scaled, Y_scaled),
        "omics_types": omics_types,
    }

    logger.info(f"CCA completed: {len(results['canonical_correlations'])} components")
    return results


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
