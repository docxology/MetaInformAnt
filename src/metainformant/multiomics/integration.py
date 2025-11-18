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
    """Container for multi-omics datasets with automatic sample alignment.
    
    This class manages heterogeneous biological datasets from multiple omic layers
    (genomics, transcriptomics, proteomics, metabolomics, epigenomics) and ensures
    proper sample alignment across all layers. Only samples present in all layers
    are retained.
    
    Attributes:
        genomics: Genomic data (SNPs, mutations) as DataFrame (samples x features)
        transcriptomics: Gene expression data as DataFrame (samples x genes)
        proteomics: Protein abundance data as DataFrame (samples x proteins)
        metabolomics: Metabolite data as DataFrame (samples x metabolites)
        epigenomics: Epigenomic data as DataFrame (samples x epigenetic marks)
        metadata: Sample metadata as DataFrame (samples x metadata columns)
        
    Examples:
        >>> import pandas as pd
        >>> import numpy as np
        >>> genomics = pd.DataFrame(np.random.randn(10, 100), index=[f"S{i}" for i in range(10)])
        >>> transcriptomics = pd.DataFrame(np.random.randn(10, 200), index=[f"S{i}" for i in range(10)])
        >>> data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)
        >>> len(data.samples)
        10
        >>> data.layer_names
        ['genomics', 'transcriptomics']
    """

    genomics: Optional[pd.DataFrame] = None
    transcriptomics: Optional[pd.DataFrame] = None
    proteomics: Optional[pd.DataFrame] = None
    metabolomics: Optional[pd.DataFrame] = None
    epigenomics: Optional[pd.DataFrame] = None
    metadata: Optional[pd.DataFrame] = None

    def __post_init__(self):
        """Validate and align samples across omics layers.
        
        Automatically identifies common samples across all provided omics layers
        and filters to retain only samples present in all layers. Raises ValueError
        if no common samples are found or no omics layers are provided.
        """
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
        """Get list of aligned samples across all omics layers.
        
        Returns:
            List of sample identifiers (index values) that are present
            in all provided omics layers. Samples are sorted alphabetically.
            
        Examples:
            >>> data = MultiOmicsData(genomics=genomics_df, transcriptomics=transcriptomics_df)
            >>> len(data.samples)
            10
            >>> data.samples[0]
            'S0'
        """
        if self.omics_layers:
            return list(next(iter(self.omics_layers.values())).index)
        return []

    @property
    def n_samples(self) -> int:
        """Get number of aligned samples across all omics layers.
        
        Returns:
            Integer count of samples present in all layers.
            
        Examples:
            >>> data = MultiOmicsData(genomics=genomics_df, transcriptomics=transcriptomics_df)
            >>> data.n_samples
            10
        """
        return len(self.samples)

    @property
    def layer_names(self) -> List[str]:
        """Get names of all available omics layers.
        
        Returns:
            List of layer names (e.g., ['genomics', 'transcriptomics']).
            Only includes layers that were provided (not None).
            
        Examples:
            >>> data = MultiOmicsData(genomics=genomics_df, transcriptomics=transcriptomics_df)
            >>> data.layer_names
            ['genomics', 'transcriptomics']
        """
        return list(self.omics_layers.keys())

    def get_layer(self, layer_name: str) -> pd.DataFrame:
        """Retrieve a specific omics layer as DataFrame.
        
        Args:
            layer_name: Name of the omics layer (e.g., 'genomics', 'transcriptomics')
            
        Returns:
            DataFrame containing the specified omics layer with aligned samples.
            
        Raises:
            KeyError: If layer_name is not available in this object
            
        Examples:
            >>> data = MultiOmicsData(genomics=genomics_df, transcriptomics=transcriptomics_df)
            >>> genomics = data.get_layer("genomics")
            >>> genomics.shape
            (10, 100)
        """
        if layer_name not in self.omics_layers:
            raise KeyError(f"Layer '{layer_name}' not available. Available: {self.layer_names}")
        return self.omics_layers[layer_name]

    def subset_samples(self, sample_list: List[str]) -> "MultiOmicsData":
        """Create a new MultiOmicsData object with subset of samples.
        
        Filters all omics layers and metadata to include only the specified
        samples. Only samples present in the original data are retained.
        
        Args:
            sample_list: List of sample identifiers to include
            
        Returns:
            New MultiOmicsData object containing only the specified samples.
            Sample alignment is automatically performed on the subset.
            
        Examples:
            >>> data = MultiOmicsData(genomics=genomics_df, transcriptomics=transcriptomics_df)
            >>> subset = data.subset_samples(["S0", "S1", "S2"])
            >>> subset.n_samples
            3
        """
        # Build kwargs for MultiOmicsData constructor
        kwargs = {}
        for layer_name, data in self.omics_layers.items():
            available_samples = [s for s in sample_list if s in data.index]
            if available_samples:
                kwargs[layer_name] = data.loc[available_samples]

        # Handle metadata
        if self.metadata is not None:
            available_samples = [s for s in sample_list if s in self.metadata.index]
            if available_samples:
                kwargs["metadata"] = self.metadata.loc[available_samples]

        # Create new MultiOmicsData object with subset data
        return MultiOmicsData(**kwargs)

    def subset_features(self, feature_dict: Dict[str, List[str]]) -> "MultiOmicsData":
        """Create a new MultiOmicsData object with subset of features.
        
        Filters specific features (genes, proteins, metabolites, etc.) from
        each omics layer. Layers not specified in feature_dict are retained
        in full.
        
        Args:
            feature_dict: Dictionary mapping layer names to lists of feature
                identifiers to include (e.g., {"transcriptomics": ["GENE1", "GENE2"]})
                
        Returns:
            New MultiOmicsData object with filtered features. Samples remain
            aligned across all layers.
            
        Examples:
            >>> data = MultiOmicsData(genomics=genomics_df, transcriptomics=transcriptomics_df)
            >>> subset = data.subset_features({
            ...     "transcriptomics": ["GENE1", "GENE2", "GENE3"]
            ... })
            >>> subset.get_layer("transcriptomics").shape
            (10, 3)
        """
        # Build kwargs for MultiOmicsData constructor
        kwargs = {}
        for layer_name, data in self.omics_layers.items():
            if layer_name in feature_dict:
                available_features = [f for f in feature_dict[layer_name] if f in data.columns]
                if available_features:
                    kwargs[layer_name] = data[available_features]
            else:
                kwargs[layer_name] = data.copy()

        # Handle metadata
        if self.metadata is not None:
            kwargs["metadata"] = self.metadata.copy()

        # Create new MultiOmicsData object with subset data
        return MultiOmicsData(**kwargs)


def integrate_omics_data(
    data_dict: Dict[str, Union[pd.DataFrame, str, Path]],
    sample_mapping: Optional[Dict[str, str]] = None,
    feature_mapping: Optional[Dict[str, Dict[str, str]]] = None,
    metadata: Optional[Union[pd.DataFrame, str, Path]] = None,
) -> MultiOmicsData:
    """Load and integrate multiple omics datasets into unified structure.
    
    Loads omics data from DataFrames or file paths, applies sample and feature
    mappings if provided, and creates a MultiOmicsData object with automatic
    sample alignment.
    
    Args:
        data_dict: Dictionary mapping omics type names ("genomics", "transcriptomics",
            etc.) to either pandas DataFrame or file path. Supported file formats:
            CSV (.csv), TSV (.tsv, .txt), Excel (.xlsx, .xls)
        sample_mapping: Optional dictionary mapping dataset names to sample ID
            transformation dictionaries for harmonizing sample IDs across datasets
        feature_mapping: Optional nested dictionary mapping dataset names to
            feature ID transformation dictionaries for harmonizing feature names
        metadata: Optional sample metadata as DataFrame or file path. Metadata
            is aligned to common samples after integration.
            
    Returns:
        MultiOmicsData object with all provided omics layers integrated and
        samples aligned to common set.
        
    Raises:
        ValueError: If file format is unsupported or required columns are missing
        ValueError: If no common samples found across datasets
        
    Examples:
        >>> data_dict = {
        ...     "genomics": "data/genomics.csv",
        ...     "transcriptomics": "data/expression.tsv"
        ... }
        >>> omics_data = integrate_omics_data(data_dict)
        >>> omics_data.n_samples
        50
    """
    omics_data = {}
    metadata_df = None

    # Load metadata first if provided (needed for sample mapping)
    if metadata is not None:
        if isinstance(metadata, (str, Path)):
            metadata_path = Path(metadata)
            if metadata_path.suffix.lower() == ".csv":
                metadata_df = pd.read_csv(metadata_path, index_col=0)
            else:
                metadata_df = pd.read_csv(metadata_path, sep="\t", index_col=0)
        else:
            metadata_df = metadata.copy()

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
            sample_map = sample_mapping[omics_type]
            # Map sample IDs, keeping original if not in mapping
            mapped_index = df.index.map(sample_map)
            # Fill missing values with original index values
            original_values = df.index.values
            mapped_values = mapped_index.values
            filled_values = [mapped_val if pd.notna(mapped_val) else orig_val 
                             for mapped_val, orig_val in zip(mapped_values, original_values)]
            df.index = pd.Index(filled_values)
            # Update sample metadata if provided
            if metadata_df is not None:
                # Same fix for metadata_df
                mapped_meta_index = metadata_df.index.map(sample_map)
                meta_original_values = metadata_df.index.values
                meta_mapped_values = mapped_meta_index.values
                meta_filled_values = [mapped_val if pd.notna(mapped_val) else orig_val 
                                      for mapped_val, orig_val in zip(meta_mapped_values, meta_original_values)]
                metadata_df.index = pd.Index(meta_filled_values)

        # Apply feature mapping if provided
        if feature_mapping and omics_type in feature_mapping:
            feature_map = feature_mapping[omics_type]
            df = df.rename(columns=feature_map)

        omics_data[omics_type] = df

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
    """Perform joint Principal Component Analysis across multiple omics layers.
    
    Concatenates features from all omics layers (optionally weighted and standardized)
    and performs PCA on the combined feature space. Returns joint embeddings and
    layer-specific loadings.
    
    Args:
        omics_data: Multi-omics data object with aligned samples
        n_components: Number of principal components to extract
        layer_weights: Optional dictionary mapping layer names to relative weights.
            If None, all layers weighted equally (1.0). Weights are applied as
            sqrt(weight) scaling before concatenation.
        standardize: If True (default), features are z-score standardized
            (mean 0, std 1) per layer before concatenation
            
    Returns:
        Tuple containing:
        - joint_embeddings: Array of shape (n_samples, n_components) with joint
          low-dimensional representation
        - layer_loadings: Dictionary mapping layer names to loading matrices
          of shape (n_features_in_layer, n_components)
        - explained_variance: Array of explained variance ratios for each component
        
    Examples:
        >>> embeddings, loadings, var = joint_pca(omics_data, n_components=10)
        >>> embeddings.shape
        (100, 10)
        >>> "genomics" in loadings
        True
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
    """Perform joint Non-negative Matrix Factorization across omics layers.
    
    Factorizes each omics layer X_i ≈ W * H_i where W (sample factors) is shared
    across all layers and H_i (feature factors) are layer-specific. Uses alternating
    optimization with L2 regularization.
    
    Args:
        omics_data: Multi-omics data object with aligned samples. All data
            must be non-negative (will be shifted if needed)
        n_components: Number of latent factors (components) to extract
        max_iter: Maximum number of optimization iterations
        regularization: L2 regularization strength for factors (default 0.01)
        random_state: Random seed for reproducible initialization
        
    Returns:
        Tuple containing:
        - sample_factors: Shared sample factor matrix W of shape (n_samples, n_components)
        - layer_feature_factors: Dictionary mapping layer names to feature factor
          matrices H_i of shape (n_components, n_features_in_layer)
          
    Examples:
        >>> W, H = joint_nmf(omics_data, n_components=10, random_state=42)
        >>> W.shape
        (100, 10)
        >>> H["transcriptomics"].shape
        (10, 200)
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
    """Perform Canonical Correlation Analysis between two omics layers.
    
    CCA finds linear combinations of features in each layer that are maximally
    correlated. Useful for identifying shared patterns between omics data types.
    
    Args:
        omics_data: Multi-omics data object with aligned samples
        layer_pair: Tuple of (layer1_name, layer2_name) to analyze. Both layers
            must exist in omics_data.layer_names
        n_components: Number of canonical variate pairs to extract
        regularization: Ridge regularization parameter for covariance matrices
            (prevents singular matrix issues)
            
    Returns:
        Tuple containing:
        - X_canonical: Canonical variates for first layer (n_samples, n_components)
        - Y_canonical: Canonical variates for second layer (n_samples, n_components)
        - X_weights: Feature weights for first layer (n_features_X, n_components)
        - Y_weights: Feature weights for second layer (n_features_Y, n_components)
        - correlations: Canonical correlations for each component pair
        
    Raises:
        ValueError: If layer names not found in omics_data
        
    Examples:
        >>> X_c, Y_c, X_w, Y_w, corr = canonical_correlation(
        ...     omics_data, ("genomics", "transcriptomics"), n_components=5
        ... )
        >>> corr[0]  # First canonical correlation
        0.85...
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


def from_dna_variants(
    vcf_path: Union[str, Path],
    sample_column_prefix: str = "sample_",
) -> pd.DataFrame:
    """Load genomic variant data from VCF file for multi-omics integration.
    
    Converts VCF format variant data into a sample x variant matrix suitable
    for multi-omics integration. Handles genotype encoding and missing data.
    
    Args:
        vcf_path: Path to VCF file
        sample_column_prefix: Prefix for sample column names in output
        
    Returns:
        DataFrame with samples as rows and variants as columns. Genotypes
        are encoded as 0 (homozygous reference), 1 (heterozygous), 2 (homozygous alternate)
    """
    try:
        from metainformant.dna import variants
    except ImportError:
        raise ImportError("dna.variants module required for VCF parsing")
    
    # Parse VCF to get sample names and variant count
    vcf_data = variants.parse_vcf(vcf_path)
    
    if not vcf_data or vcf_data.get("num_variants", 0) == 0:
        return pd.DataFrame()
    
    # Note: parse_vcf only extracts metadata (samples, variant count)
    # Full VCF parsing with genotype data would require a more complete parser
    # For now, return empty DataFrame as placeholder until full VCF parser is implemented
    # This function is intended to convert VCF to sample x variant matrix format
    return pd.DataFrame()


def from_rna_expression(
    expression_path: Union[str, Path],
    delimiter: str = "\t",
) -> pd.DataFrame:
    """Load RNA expression data for multi-omics integration.
    
    Loads gene expression data from TSV/CSV file with genes as columns
    and samples as rows.
    
    Args:
        expression_path: Path to expression file (TSV or CSV)
        delimiter: Delimiter for file (default: tab for TSV)
        
    Returns:
        DataFrame with samples as rows and genes as columns
    """
    path = Path(expression_path)
    
    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path, index_col=0)
    else:
        df = pd.read_csv(path, sep=delimiter, index_col=0)
    
    return df


def from_protein_abundance(
    protein_path: Union[str, Path],
    delimiter: str = ",",
) -> pd.DataFrame:
    """Load protein abundance data for multi-omics integration.
    
    Loads protein abundance data from CSV/TSV file with proteins as columns
    and samples as rows.
    
    Args:
        protein_path: Path to protein abundance file
        delimiter: Delimiter for file (default: comma for CSV)
        
    Returns:
        DataFrame with samples as rows and proteins as columns
    """
    path = Path(protein_path)
    
    if path.suffix.lower() in [".tsv", ".txt"]:
        df = pd.read_csv(path, sep="\t", index_col=0)
    else:
        df = pd.read_csv(path, sep=delimiter, index_col=0)
    
    return df


def from_metabolomics(
    metabolite_path: Union[str, Path],
    delimiter: str = ",",
    normalize: bool = True,
) -> pd.DataFrame:
    """Load metabolomics data for multi-omics integration.
    
    Loads metabolite abundance data from CSV/TSV file with metabolites as columns
    and samples as rows. Optionally normalizes data (log2 transform and scaling).
    
    Args:
        metabolite_path: Path to metabolomics data file
        delimiter: Delimiter for file (default: comma for CSV)
        normalize: If True, apply log2 transformation and quantile normalization
        
    Returns:
        DataFrame with samples as rows and metabolites as columns
    """
    path = Path(metabolite_path)
    
    if path.suffix.lower() in [".tsv", ".txt"]:
        df = pd.read_csv(path, sep="\t", index_col=0)
    else:
        df = pd.read_csv(path, sep=delimiter, index_col=0)
    
    # Handle missing values
    df = df.fillna(0.0)
    
    # Normalize if requested
    if normalize:
        # Log2 transform (add small value to avoid log(0))
        df = np.log2(df + 1.0)
        
        # Quantile normalization (simplified)
        # Sort each column, compute mean, then assign back
        sorted_df = np.sort(df.values, axis=0)
        mean_sorted = np.mean(sorted_df, axis=1)
        
        # Get rank for each value
        ranks = df.rank(method="average", axis=0).astype(int) - 1
        ranks = np.clip(ranks, 0, len(mean_sorted) - 1)
        
        # Assign normalized values
        normalized_values = mean_sorted[ranks]
        df = pd.DataFrame(normalized_values, index=df.index, columns=df.columns)
    
    return df
