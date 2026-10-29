# Multi-Omics Integration API Reference

Complete API reference for the multi-omics integration module.

## Classes

### MultiOmicsData

Container for multi-omics datasets with automatic sample alignment.

```python
class MultiOmicsData:
    genomics: Optional[pd.DataFrame] = None
    transcriptomics: Optional[pd.DataFrame] = None
    proteomics: Optional[pd.DataFrame] = None
    metabolomics: Optional[pd.DataFrame] = None
    epigenomics: Optional[pd.DataFrame] = None
    metadata: Optional[pd.DataFrame] = None
```

#### Attributes

- **genomics**: Genomic data (SNPs, mutations) as DataFrame (samples × features)
- **transcriptomics**: Gene expression data as DataFrame (samples × genes)
- **proteomics**: Protein abundance data as DataFrame (samples × proteins)
- **metabolomics**: Metabolite data as DataFrame (samples × metabolites)
- **epigenomics**: Epigenomic data as DataFrame (samples × epigenetic marks)
- **metadata**: Sample metadata as DataFrame (samples × metadata columns)

#### Properties

##### `samples: List[str]`
Get list of aligned samples across all omics layers.

**Returns:**
- List of sample identifiers (index values) that are present in all provided omics layers. Samples are sorted alphabetically.

##### `n_samples: int`
Get number of aligned samples across all omics layers.

**Returns:**
- Integer count of samples present in all layers.

##### `layer_names: List[str]`
Get names of all available omics layers.

**Returns:**
- List of layer names (e.g., ['genomics', 'transcriptomics']). Only includes layers that were provided (not None).

#### Methods

##### `get_layer(layer_name: str) -> pd.DataFrame`
Retrieve a specific omics layer as DataFrame.

**Parameters:**
- `layer_name` (str): Name of the omics layer (e.g., 'genomics', 'transcriptomics')

**Returns:**
- DataFrame containing the specified omics layer with aligned samples.

**Raises:**
- `KeyError`: If layer_name is not available in this object

##### `subset_samples(sample_list: List[str]) -> MultiOmicsData`
Create a new MultiOmicsData object with subset of samples.

**Parameters:**
- `sample_list` (List[str]): List of sample identifiers to include

**Returns:**
- New MultiOmicsData object containing only the specified samples. Sample alignment is automatically performed on the subset.

**Raises:**
- `ValueError`: If sample_list is empty or no samples found

##### `subset_features(feature_dict: Dict[str, List[str]]) -> MultiOmicsData`
Create a new MultiOmicsData object with subset of features.

**Parameters:**
- `feature_dict` (Dict[str, List[str]]): Dictionary mapping layer names to lists of feature identifiers to include (e.g., {"transcriptomics": ["GENE1", "GENE2"]})

**Returns:**
- New MultiOmicsData object with filtered features. Samples remain aligned across all layers.

**Raises:**
- `ValueError`: If feature_dict is empty or no valid features found

##### `to_dict() -> Dict[str, Any]`
Convert MultiOmicsData to dictionary for serialization.

**Returns:**
- Dictionary containing all omics layers and metadata as DataFrames.

##### `save(output_path: Union[str, Path]) -> None`
Save MultiOmicsData to directory using core.io.

**Parameters:**
- `output_path` (Union[str, Path]): Directory path to save data files

**Note:**
- Saves each omics layer as a separate CSV file and metadata if present.
- Uses core.io for file operations to ensure proper path handling.

##### `load(input_path: Union[str, Path]) -> MultiOmicsData` (classmethod)
Load MultiOmicsData from directory using core.io.

**Parameters:**
- `input_path` (Union[str, Path]): Directory path containing data files

**Returns:**
- MultiOmicsData object with loaded layers

**Raises:**
- `FileNotFoundError`: If directory not found
- `ValueError`: If no data files found

## Functions

### integrate_omics_data

Load and integrate multiple omics datasets into unified structure.

```python-snippet
def integrate_omics_data(
    data_dict: Dict[str, Union[pd.DataFrame, str, Path]],
    sample_mapping: Optional[Dict[str, Dict[str, str]]] = None,
    feature_mapping: Optional[Dict[str, Dict[str, str]]] = None,
    metadata: Optional[Union[pd.DataFrame, str, Path]] = None,
) -> MultiOmicsData:
```

**Parameters:**

- `data_dict` (Dict[str, Union[pd.DataFrame, str, Path]]): Dictionary mapping omics type names ("genomics", "transcriptomics", etc.) to either pandas DataFrame or file path. Supported file formats: CSV (.csv), TSV (.tsv, .txt), Excel (.xlsx, .xls). Gzip compression is automatically handled via core.io.
- `sample_mapping` (Optional[Dict[str, Dict[str, str]]]): Optional nested dictionary mapping dataset names to sample ID transformation dictionaries for harmonizing sample IDs across datasets. Format: {omics_type: {old_id: new_id}}
- `feature_mapping` (Optional[Dict[str, Dict[str, str]]]): Optional nested dictionary mapping dataset names to feature ID transformation dictionaries for harmonizing feature names. Format: {omics_type: {old_feature: new_feature}}
- `metadata` (Optional[Union[pd.DataFrame, str, Path]]): Optional sample metadata as DataFrame or file path. Metadata is aligned to common samples after integration.

**Returns:**
- MultiOmicsData object with all provided omics layers integrated and samples aligned to common set.

**Raises:**
- `ValueError`: If file format is unsupported or required columns are missing, or if no common samples found across datasets
- `FileNotFoundError`: If file paths do not exist
- `IOError`: If file read fails

### joint_pca

Perform joint Principal Component Analysis across multiple omics layers.

```python-snippet
def joint_pca(
    omics_data: MultiOmicsData,
    n_components: int = 50,
    layer_weights: Optional[Dict[str, float]] = None,
    standardize: bool = True,
) -> Tuple[np.ndarray, Dict[str, np.ndarray], np.ndarray]:
```

**Parameters:**

- `omics_data` (MultiOmicsData): Multi-omics data object with aligned samples
- `n_components` (int): Number of principal components to extract. Must be positive and less than min(n_samples, n_features). Will be automatically capped if too large. Default: 50
- `layer_weights` (Optional[Dict[str, float]]): Optional dictionary mapping layer names to relative weights. If None, all layers weighted equally (1.0). Weights are applied as sqrt(weight) scaling before concatenation.
- `standardize` (bool): If True (default), features are z-score standardized (mean 0, std 1) per layer before concatenation

**Returns:**
- Tuple containing:
  - `joint_embeddings`: Array of shape (n_samples, n_components) with joint low-dimensional representation
  - `layer_loadings`: Dictionary mapping layer names to loading matrices of shape (n_features_in_layer, n_components)
  - `explained_variance`: Array of explained variance ratios for each component

**Raises:**
- `ValueError`: If n_components is invalid or data is insufficient
- `MemoryError`: If insufficient memory for covariance matrix computation

### joint_nmf

Perform joint Non-negative Matrix Factorization across omics layers.

```python-snippet
def joint_nmf(
    omics_data: MultiOmicsData,
    n_components: int = 20,
    max_iter: int = 200,
    regularization: float = 0.01,
    random_state: Optional[int] = None,
    tolerance: float = 1e-6,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
```

**Parameters:**

- `omics_data` (MultiOmicsData): Multi-omics data object with aligned samples. All data must be non-negative (will be shifted if needed)
- `n_components` (int): Number of latent factors (components) to extract. Must be positive and less than min(n_samples, min_features_per_layer). Default: 20
- `max_iter` (int): Maximum number of optimization iterations. Default: 200
- `regularization` (float): L2 regularization strength for factors. Default: 0.01
- `random_state` (Optional[int]): Random seed for reproducible initialization
- `tolerance` (float): Convergence tolerance for early stopping. Algorithm stops if relative change in reconstruction error is below this. Default: 1e-6

**Returns:**
- Tuple containing:
  - `sample_factors`: Shared sample factor matrix W of shape (n_samples, n_components)
  - `layer_feature_factors`: Dictionary mapping layer names to feature factor matrices H_i of shape (n_components, n_features_in_layer)

**Raises:**
- `ValueError`: If n_components is invalid or data is insufficient

### canonical_correlation

Perform ridge-regularized canonical correlation analysis (CCA) between two omics layers.

```python-snippet
def canonical_correlation(
    multiomics_data: Union[MultiOmicsData, Dict[str, pd.DataFrame]],
    layers: Optional[Tuple[str, str]] = None,
    layer_pair: Optional[Tuple[str, str]] = None,
    n_components: int = 10,
    regularization: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
```

Each layer is standardized first (columns centred and scaled to unit sample variance), so the within-layer covariance matrices `Cxx` and `Cyy` are correlation matrices. The canonical directions `a` and `b` maximise the correlation of the projected layers under the constraints `a'(Cxx + λI)a = b'(Cyy + λI)b = 1`, where `λ` is `regularization`. The canonical correlations are the singular values of `(Cxx + λI)^(-1/2) Cxy (Cyy + λI)^(-1/2)`. The solution is closed form: it comes from thin SVDs of the standardized layers, with no iteration and no feature-by-feature matrices, so layers with many more features than samples are supported.

**Parameters:**

- `multiomics_data` (MultiOmicsData or Dict[str, pd.DataFrame]): Multi-omics data object, or a dictionary that maps layer names to samples × features DataFrames. A `MultiOmicsData` object is already sample-aligned. The rows of a plain dictionary are matched by position, so list the samples in the same order in both layers.
- `layers` (Optional[Tuple[str, str]]): Names of the two layers to correlate; both must exist in the data. If omitted, the data must hold exactly two layers, which are used in their stored order.
- `layer_pair` (Optional[Tuple[str, str]]): Alias for `layers`.
- `n_components` (int): Number of canonical pairs to extract; must be at least 1. Default: 10. The number of pairs is capped at the smaller of the two layers' feature counts. Fewer pairs are also returned when the numerical rank of the layers supports fewer (the rank is at most `n_samples - 1`); that case logs a warning.
- `regularization` (float): Ridge parameter `λ >= 0` on the scale of a correlation matrix, added to the diagonal of each within-layer covariance matrix. It must be finite. Default: 0.0, which is classical, unregularized CCA. Larger values shrink every canonical correlation, and each correlation is non-increasing in `λ`. Without regularization, when the ranks of the two layers add up to more than `n_samples - 1` (for example with far more features than samples), at least that excess number of canonical correlations equal exactly 1, a perfect in-sample fit, and `ValueError` is raised. Use `λ > 0` for such data. The same refusal applies when either layer alone has numerical rank `n_samples - 1`.

The function takes no other keyword arguments.

**Returns:**
- Tuple `(X_c, Y_c, X_weights, Y_weights, correlations)`:
  - `X_c`, `Y_c`: Canonical variates for the first and second layer, shape (n_samples, n_pairs)
  - `X_weights`, `Y_weights`: Canonical coefficients on the standardized features, shape (n_features, n_pairs). They reproduce the variates: `X_c = standardized_X @ X_weights`, where `standardized_X` is the layer with each column centred and divided by its sample standard deviation (`ddof=1`); likewise for `Y_c`. With `regularization == 0` the variates have unit variance. With `regularization > 0` the weights satisfy the penalised constraint above, so the variates have variance at most 1.
  - `correlations`: Canonical correlations, shape (n_pairs,), in [0, 1], in descending order. With `regularization == 0` each equals the sample correlation of the matching pair of variates. With `regularization > 0` each is the penalised correlation, which is a lower bound on the sample correlation of the returned variates.
- The sign of each pair is fixed so that the largest-magnitude X weight is positive, and both members of a pair flip together.

**Note:**
- Missing values (NaN) are replaced by the observed mean of their own feature before standardization. A wholly unobserved feature raises `ValueError`; it has no mean to impute.
- Constant columns become all zeros after standardization and receive zero weight, so the result equals the one computed without them.

**Raises:**
- `ImportError`: If scikit-learn is not installed
- `ValidationError` (`metainformant.core.utils.errors`, also a `ValueError` subclass): If `regularization` is negative or not finite, if `n_components < 1`, or if `layers` is omitted and the data does not hold exactly two layers
- `ValueError`: If a requested layer is not in the data, if the two layers have different numbers of samples, if there are fewer than 2 samples, if a layer has no non-constant feature or has a wholly unobserved feature, or if unregularized layers have combined numerical rank above `n_samples - 1`
- `TypeError`: If an unknown keyword argument is passed

### from_dna_variants

Convert DNA variant data (VCF) to genomics layer DataFrame.

```python-snippet
def from_dna_variants(
    vcf_data: Union[str, Path, pd.DataFrame],
    sample_ids: Optional[List[str]] = None,
    variant_ids: Optional[List[str]] = None,
) -> pd.DataFrame:
```

**Parameters:**

- `vcf_data` (Union[str, Path, pd.DataFrame]): Path to a VCF file, a VCF-style DataFrame, or an existing sample-by-variant matrix.
- `sample_ids` (Optional[List[str]]): Optional list of sample IDs to include. If None, includes all samples.
- `variant_ids` (Optional[List[str]]): Optional list of variant IDs to include. If None, includes all variants.

**Returns:**
- DataFrame with samples as index, variants as columns, and genotype values (0, 1, 2 for homozygous ref, heterozygous, homozygous alt) as values.

**Raises:**
- `FileNotFoundError`: If VCF file does not exist
- `ValueError`: If VCF parsing fails or requested sample/variant IDs are absent
- `ImportError`: If DNA module not available

### from_rna_expression

Convert RNA expression data to transcriptomics layer DataFrame.

```python-snippet
def from_rna_expression(
    expression_data: Union[str, Path, pd.DataFrame],
    normalize: bool = True,
    sample_ids: Optional[List[str]] = None,
    gene_ids: Optional[List[str]] = None,
    transpose: bool = False,
) -> pd.DataFrame:
```

**Parameters:**

- `expression_data` (Union[str, Path, pd.DataFrame]): Path to expression file (CSV or TSV) or an existing expression matrix.
- `normalize` (bool): If True, apply `log1p` when every observed value is non-negative, then column-wise z-score normalization. Zeros and missing values do not disable the log transform. If any value is negative, only z-score normalization is applied.
- `sample_ids` (Optional[List[str]]): Optional list of sample IDs to include
- `gene_ids` (Optional[List[str]]): Optional list of gene IDs to include
- `transpose` (bool): If True, assumes file has genes as rows and samples as columns. If False (default), assumes samples as rows and genes as columns.

**Returns:**
- DataFrame with samples as index, genes as columns, and expression values.

**Raises:**
- `FileNotFoundError`: If expression file does not exist
- `ValueError`: If file format is invalid

### from_protein_abundance

Convert protein abundance data to proteomics layer DataFrame.

```python-snippet
def from_protein_abundance(
    protein_data: Union[str, Path, pd.DataFrame],
    normalize: bool = True,
    sample_ids: Optional[List[str]] = None,
    protein_ids: Optional[List[str]] = None,
    transpose: bool = False,
) -> pd.DataFrame:
```

**Parameters:**

- `protein_data` (Union[str, Path, pd.DataFrame]): Path to protein matrix file (CSV or TSV) or an existing abundance matrix.
- `normalize` (bool): If True, apply column-wise z-score normalization.
- `sample_ids` (Optional[List[str]]): Optional list of sample IDs to include
- `protein_ids` (Optional[List[str]]): Optional list of protein IDs to include
- `transpose` (bool): If True, assumes file has proteins as rows and samples as columns. If False (default), assumes samples as rows and proteins as columns. Not applicable for FASTA files.

**Returns:**
- DataFrame with samples as index, proteins as columns, and abundance values.

**Raises:**
- `FileNotFoundError`: If protein file does not exist
- `ValueError`: If file format is invalid
- `ImportError`: If protein module not available (for FASTA files)

## Error Handling

All functions include comprehensive error handling:

- **Validation errors**: Raised for invalid inputs (empty data, wrong types, etc.)
- **File I/O errors**: Raised for missing files or read failures
- **Computation errors**: Raised for numerical issues (singular matrices, etc.)
- **Memory errors**: Raised for insufficient memory

All errors include descriptive messages and context.

## Logging

The module uses structured logging via `metainformant.core.utils.logging`:

- **INFO**: Major operations (loading data, starting analysis)
- **DEBUG**: Detailed information (file paths, data shapes)
- **WARNING**: Non-fatal issues (missing samples, data shifts)
- **ERROR**: Fatal errors (file read failures, computation errors)

## Examples

See [Module README](README.md) and [Index Documentation](./index.md) for comprehensive usage examples.
