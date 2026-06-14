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

from metainformant.core.data import validation
from metainformant.core.io import dump_json, ensure_directory, load_json, open_text_auto
from metainformant.core.utils import errors, logging

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
            raise ValueError("At least one omics layer must be provided")

        # Validate data types
        for key, value in data.items():
            if not isinstance(value, pd.DataFrame):
                raise TypeError(f"Data for '{key}' must be pandas DataFrame, got {type(value).__name__}")

        # Validate non-empty features
        for key, value in data.items():
            if value.shape[1] == 0:
                raise ValueError(f"Layer '{key}' has no features")

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

        # Warn if samples don't fully overlap
        import warnings

        all_samples = set()
        for df in self.data.values():
            all_samples.update(df.index)
        if len(common_samples) < len(all_samples):
            warnings.warn(
                f"Only {len(common_samples)} samples are common across all layers "
                f"(out of {len(all_samples)} total unique samples)",
                UserWarning,
                stacklevel=3,
            )

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

    def __getattr__(self, name: str) -> Any:
        """Allow direct attribute access to layers by name (e.g., .transcriptomics)."""
        if name.startswith("_"):
            raise AttributeError(name)
        if hasattr(self, "data") and name in self.data:
            return self.data[name]
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def subset_samples(self, sample_ids: List[str]) -> "MultiOmicsData":
        """Create subset with specified samples."""
        subset_data = {}
        available_samples = [s for s in sample_ids if s in self.samples]
        for omics_type, df in self.data.items():
            subset_data[omics_type] = df.loc[available_samples]

        # Subset metadata too
        subset_metadata = None
        if isinstance(self._metadata, pd.DataFrame) and not self._metadata.empty:
            meta_available = [s for s in available_samples if s in self._metadata.index]
            if meta_available:
                subset_metadata = self._metadata.loc[meta_available]

        return MultiOmicsData(
            data=subset_data, sample_ids=available_samples, feature_ids=self.feature_ids, metadata=subset_metadata
        )

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

    def to_dict(self) -> Dict[str, Any]:
        """Convert the container to a dictionary with copy-safe DataFrames."""
        metadata = self._metadata.copy() if isinstance(self._metadata, pd.DataFrame) else pd.DataFrame()
        return {
            "data": {layer_name: layer_data.copy() for layer_name, layer_data in self.data.items()},
            "sample_ids": list(self.sample_ids or []),
            "feature_ids": {layer: list(features) for layer, features in self.feature_ids.items()},
            "metadata": metadata,
        }

    def save(self, output_path: Union[str, Path]) -> None:
        """Save all omics layers and metadata to a directory."""
        output_dir = ensure_directory(output_path)
        manifest = {
            "layers": list(self.data.keys()),
            "sample_ids": list(self.sample_ids or []),
            "feature_ids": {layer: list(features) for layer, features in self.feature_ids.items()},
            "metadata": bool(isinstance(self._metadata, pd.DataFrame) and not self._metadata.empty),
        }

        for layer_name, layer_data in self.data.items():
            layer_data.to_csv(output_dir / f"{layer_name}.csv")

        if manifest["metadata"]:
            self._metadata.to_csv(output_dir / "metadata.csv")

        dump_json(manifest, output_dir / "manifest.json")

    @classmethod
    def load(cls, input_path: Union[str, Path]) -> "MultiOmicsData":
        """Load a MultiOmicsData directory written by :meth:`save`."""
        input_dir = Path(input_path)
        manifest_path = input_dir / "manifest.json"
        if not manifest_path.exists():
            raise FileNotFoundError(f"MultiOmicsData manifest not found: {manifest_path}")

        manifest = load_json(manifest_path)
        layers = manifest.get("layers", [])
        if not layers:
            raise ValueError(f"MultiOmicsData manifest contains no layers: {manifest_path}")

        data = {}
        for layer_name in layers:
            layer_path = input_dir / f"{layer_name}.csv"
            if not layer_path.exists():
                raise FileNotFoundError(f"Layer file not found for '{layer_name}': {layer_path}")
            data[layer_name] = pd.read_csv(layer_path, index_col=0)

        metadata = None
        metadata_path = input_dir / "metadata.csv"
        if manifest.get("metadata") and metadata_path.exists():
            metadata = pd.read_csv(metadata_path, index_col=0)

        return cls(
            data=data,
            sample_ids=manifest.get("sample_ids"),
            feature_ids=manifest.get("feature_ids") or {},
            metadata=metadata,
        )


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
    sample_mapping = kwargs.pop("sample_mapping", None)
    feature_mapping = kwargs.pop("feature_mapping", None)

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
                elif path.suffix == ".tsv":
                    processed_data[key] = pd.read_csv(path, sep="\t", index_col=0)
                elif path.suffix == ".parquet":
                    processed_data[key] = pd.read_parquet(path)
                else:
                    raise errors.ValidationError(f"Unsupported file format: {path.suffix}")
            else:
                processed_data[key] = value

        if sample_mapping:
            for layer, mapping in sample_mapping.items():
                if layer in processed_data:
                    processed_data[layer] = processed_data[layer].rename(index=mapping)

        if feature_mapping:
            for layer, mapping in feature_mapping.items():
                if layer in processed_data:
                    processed_data[layer] = processed_data[layer].rename(columns=mapping)

        # Load metadata from file if it's a path
        if "metadata" in kwargs and isinstance(kwargs["metadata"], (str, Path)):
            meta_path = Path(kwargs["metadata"])
            if meta_path.suffix == ".csv":
                kwargs["metadata"] = pd.read_csv(meta_path, index_col=0)
            elif meta_path.suffix == ".tsv":
                kwargs["metadata"] = pd.read_csv(meta_path, sep="\t", index_col=0)
            elif meta_path.suffix == ".parquet":
                kwargs["metadata"] = pd.read_parquet(meta_path)

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
    if sample_mapping:
        for layer, mapping in sample_mapping.items():
            if layer in available_omics:
                available_omics[layer] = available_omics[layer].rename(index=mapping)

    if feature_mapping:
        for layer, mapping in feature_mapping.items():
            if layer in available_omics:
                available_omics[layer] = available_omics[layer].rename(columns=mapping)

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
    layer_pair: Optional[Tuple[str, str]] = None,
    n_components: int = 10,
    regularization: float = 0.0,
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

    # Determine layers to use - support both 'layers' and 'layer_pair' kwargs
    if layers is None and layer_pair is not None:
        layers = layer_pair
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
    order = np.argsort(correlations)[::-1]
    X_c = X_c[:, order]
    Y_c = Y_c[:, order]
    x_weights = cca.x_weights_[:, order]
    y_weights = cca.y_weights_[:, order]
    correlations = correlations[order]

    logger.info(f"CCA completed: {len(correlations)} components")
    return X_c, Y_c, x_weights, y_weights, correlations


def _normalize_variant_id(chrom: object, pos: object, variant_id: object, ref: object, alt: object) -> str:
    """Return a stable variant identifier from VCF fields."""
    variant_name = "" if pd.isna(variant_id) else str(variant_id)
    if variant_name and variant_name != ".":
        return variant_name
    return f"{chrom}:{pos}:{ref}:{alt}"


def _genotype_to_dosage(sample_value: object, format_keys: List[str]) -> float:
    """Convert a VCF sample genotype field to reference-alt dosage."""
    if sample_value is None or pd.isna(sample_value):
        return np.nan

    fields = str(sample_value).split(":")
    gt_index = format_keys.index("GT") if "GT" in format_keys else 0
    if gt_index >= len(fields):
        return np.nan

    genotype = fields[gt_index]
    if not genotype or genotype == "." or "." in genotype:
        return np.nan

    alleles = genotype.replace("|", "/").split("/")
    dosage = 0
    for allele in alleles:
        try:
            dosage += 0 if int(allele) == 0 else 1
        except ValueError:
            return np.nan
    return float(dosage)


def _subset_variant_matrix(
    matrix: pd.DataFrame,
    sample_ids: Optional[List[str]] = None,
    variant_ids: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Subset a sample-by-variant matrix with clear missing-id errors."""
    if sample_ids is not None:
        missing_samples = [sample for sample in sample_ids if sample not in matrix.index]
        if missing_samples:
            raise ValueError(f"Sample IDs not found in variant data: {missing_samples}")
        matrix = matrix.loc[sample_ids]

    if variant_ids is not None:
        missing_variants = [variant for variant in variant_ids if variant not in matrix.columns]
        if missing_variants:
            raise ValueError(f"Variant IDs not found in variant data: {missing_variants}")
        matrix = matrix.loc[:, variant_ids]

    return matrix


def _read_tabular_omics_data(data: Union[pd.DataFrame, str, Path], data_name: str) -> pd.DataFrame:
    """Read a DataFrame or CSV/TSV matrix with sample IDs in the first column."""
    if isinstance(data, pd.DataFrame):
        return data.copy()

    if not isinstance(data, (str, Path)):
        raise TypeError(f"{data_name} must be a path or pandas DataFrame, got {type(data).__name__}")

    path = Path(data)
    if not path.exists():
        raise FileNotFoundError(f"{data_name} file not found: {path}")

    suffix = path.suffix.lower()
    if suffix in {".tsv", ".tab", ".txt"}:
        return pd.read_csv(path, sep="\t", index_col=0)
    if suffix == ".csv":
        return pd.read_csv(path, index_col=0)
    raise ValueError(f"Unsupported {data_name} file extension: {suffix}")


def _subset_omics_matrix(
    matrix: pd.DataFrame,
    sample_ids: Optional[List[str]] = None,
    feature_ids: Optional[List[str]] = None,
    *,
    feature_label: str,
) -> pd.DataFrame:
    """Subset a sample-by-feature matrix with explicit missing-id errors."""
    if sample_ids is not None:
        missing_samples = [sample for sample in sample_ids if sample not in matrix.index]
        if missing_samples:
            raise ValueError(f"Sample IDs not found in omics data: {missing_samples}")
        matrix = matrix.loc[sample_ids]

    if feature_ids is not None:
        missing_features = [feature for feature in feature_ids if feature not in matrix.columns]
        if missing_features:
            raise ValueError(f"{feature_label} IDs not found in omics data: {missing_features}")
        matrix = matrix.loc[:, feature_ids]

    return matrix


def _zscore_dataframe(data: pd.DataFrame) -> pd.DataFrame:
    """Return column-wise z-scores, preserving DataFrame labels."""
    if HAS_SKLEARN and StandardScaler is not None:
        scaler = StandardScaler()
        values = scaler.fit_transform(data)
    else:
        std = data.std(axis=0, ddof=0).replace(0, 1.0)
        values = (data - data.mean(axis=0)) / std
    return pd.DataFrame(values, index=data.index, columns=data.columns)


def _vcf_dataframe_to_matrix(
    vcf_df: pd.DataFrame,
    sample_ids: Optional[List[str]] = None,
    variant_ids: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Convert a DataFrame containing VCF records into a genotype dosage matrix."""
    chrom_column = "#CHROM" if "#CHROM" in vcf_df.columns else "CHROM"
    required = {chrom_column, "POS", "REF", "ALT", "FORMAT"}
    missing = required - set(vcf_df.columns)
    if missing:
        raise ValueError(f"VCF DataFrame missing required columns: {sorted(missing)}")

    format_index = list(vcf_df.columns).index("FORMAT")
    sample_columns = list(vcf_df.columns)[format_index + 1 :]
    if not sample_columns:
        raise ValueError("VCF data does not contain sample genotype columns")

    variant_names: List[str] = []
    sample_values: Dict[str, List[float]] = {sample: [] for sample in sample_columns}

    for _, row in vcf_df.iterrows():
        variant_name = _normalize_variant_id(
            row[chrom_column],
            row["POS"],
            row["ID"] if "ID" in vcf_df.columns else ".",
            row["REF"],
            row["ALT"],
        )
        variant_names.append(variant_name)

        format_keys = str(row["FORMAT"]).split(":") if not pd.isna(row["FORMAT"]) else ["GT"]
        for sample in sample_columns:
            sample_values[sample].append(_genotype_to_dosage(row[sample], format_keys))

    matrix = pd.DataFrame(sample_values, index=variant_names, dtype=float).T
    return _subset_variant_matrix(matrix, sample_ids=sample_ids, variant_ids=variant_ids)


def _read_vcf_records(vcf_path: Path) -> pd.DataFrame:
    """Read VCF records into a DataFrame without requiring external VCF packages."""
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    header: Optional[List[str]] = None
    records: List[List[str]] = []

    with open_text_auto(vcf_path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.split("\t")
                continue
            if line.startswith("#"):
                continue
            if header is None:
                raise ValueError(f"VCF record encountered before #CHROM header at {vcf_path}:{line_number}")

            fields = line.split("\t")
            if len(fields) != len(header):
                raise ValueError(
                    f"VCF row has {len(fields)} fields but header has {len(header)} at {vcf_path}:{line_number}"
                )
            records.append(fields)

    if header is None:
        raise ValueError(f"VCF header not found in {vcf_path}")
    if not records:
        raise ValueError(f"No variant records found in {vcf_path}")

    return pd.DataFrame(records, columns=header)


def from_dna_variants(
    vcf_data: Union[pd.DataFrame, str, Path],
    sample_ids: Optional[List[str]] = None,
    variant_ids: Optional[List[str]] = None,
    **kwargs,
) -> pd.DataFrame:
    """Convert DNA variant data for multi-omics integration.

    Args:
        vcf_data: VCF path, VCF-style DataFrame, or sample-by-variant matrix
        sample_ids: Optional sample IDs to retain
        variant_ids: Optional variant IDs to retain
        **kwargs: Conversion parameters

    Returns:
        Sample-by-variant genotype dosage matrix suitable for integration
    """
    logger.info("Converting DNA variant data for integration")

    if isinstance(vcf_data, (str, Path)):
        vcf_records = _read_vcf_records(Path(vcf_data))
        return _vcf_dataframe_to_matrix(vcf_records, sample_ids=sample_ids, variant_ids=variant_ids)

    if not isinstance(vcf_data, pd.DataFrame):
        raise TypeError(f"vcf_data must be a VCF path or pandas DataFrame, got {type(vcf_data).__name__}")

    if "FORMAT" in vcf_data.columns and ({"#CHROM", "CHROM"} & set(vcf_data.columns)):
        return _vcf_dataframe_to_matrix(vcf_data, sample_ids=sample_ids, variant_ids=variant_ids)

    processed_data = vcf_data.copy()
    return _subset_variant_matrix(processed_data, sample_ids=sample_ids, variant_ids=variant_ids)


def from_rna_expression(
    expression_data: Union[pd.DataFrame, str, Path],
    normalize: bool = True,
    sample_ids: Optional[List[str]] = None,
    gene_ids: Optional[List[str]] = None,
    transpose: bool = False,
    **kwargs,
) -> pd.DataFrame:
    """Convert RNA expression data for multi-omics integration.

    Args:
        expression_data: RNA-seq or microarray expression matrix path/DataFrame
        normalize: Whether to normalize the data
        sample_ids: Optional sample IDs to retain
        gene_ids: Optional gene IDs to retain
        transpose: If True, transpose input from genes-by-samples to samples-by-genes
        **kwargs: Conversion parameters

    Returns:
        Processed RNA data suitable for integration
    """
    logger.info("Converting RNA expression data for integration")

    processed_data = _read_tabular_omics_data(expression_data, "expression_data")
    if transpose:
        processed_data = processed_data.T
    processed_data = processed_data.apply(pd.to_numeric, errors="raise")
    processed_data = _subset_omics_matrix(
        processed_data, sample_ids=sample_ids, feature_ids=gene_ids, feature_label="Gene"
    )

    if normalize:
        # Simple normalization - log transform and z-score
        if (processed_data > 0).all().all():
            processed_data = np.log1p(processed_data)

        processed_data = _zscore_dataframe(processed_data)

    return processed_data


def from_protein_abundance(
    protein_data: Union[pd.DataFrame, str, Path],
    normalize: bool = True,
    sample_ids: Optional[List[str]] = None,
    protein_ids: Optional[List[str]] = None,
    transpose: bool = False,
    **kwargs,
) -> pd.DataFrame:
    """Convert protein abundance data for multi-omics integration.

    Args:
        protein_data: Protein abundance matrix path/DataFrame
        normalize: Whether to normalize the data
        sample_ids: Optional sample IDs to retain
        protein_ids: Optional protein IDs to retain
        transpose: If True, transpose input from proteins-by-samples to samples-by-proteins
        **kwargs: Conversion parameters

    Returns:
        Processed protein data suitable for integration
    """
    logger.info("Converting protein abundance data for integration")

    processed_data = _read_tabular_omics_data(protein_data, "protein_data")
    if transpose:
        processed_data = processed_data.T
    processed_data = processed_data.apply(pd.to_numeric, errors="raise")
    processed_data = _subset_omics_matrix(
        processed_data, sample_ids=sample_ids, feature_ids=protein_ids, feature_label="Protein"
    )

    if normalize:
        # Z-score normalization
        processed_data = _zscore_dataframe(processed_data)

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


def find_multiomics_modules(
    omics_data: Union["MultiOmicsData", Dict[str, pd.DataFrame]], n_modules: int = 10, **kwargs
) -> Dict[str, Any]:
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
    sample_weights, feature_loadings = joint_nmf(omics_data, n_components=n_modules, **kwargs)
    data_dict = omics_data.data if hasattr(omics_data, "data") else omics_data
    actual_modules = sample_weights.shape[1]

    # Interpret modules
    modules = {}

    for i in range(actual_modules):
        module_features = {}

        for omics_type, components in feature_loadings.items():
            # Find features with high loading in this component
            loadings = components[i, :]
            top_n = min(10, len(loadings))
            top_features = np.argsort(loadings)[-top_n:][::-1]

            module_features[omics_type] = {
                "features": [data_dict[omics_type].columns[j] for j in top_features],
                "loadings": loadings[top_features].tolist(),
            }

        modules[f"module_{i+1}"] = {
            "features": module_features,
            "sample_weights": sample_weights[:, i].tolist(),
        }

    results = {
        "modules": modules,
        "nmf_results": {"W": sample_weights, "H_dict": feature_loadings},
        "n_modules": actual_modules,
    }

    logger.info(f"Multi-omics module detection completed: {actual_modules} modules found")
    return results
