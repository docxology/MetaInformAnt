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

```python
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

```python
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

```python
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

Perform Canonical Correlation Analysis between two omics layers.

```python
def canonical_correlation(
    omics_data: MultiOmicsData,
    layer_pair: Tuple[str, str],
    n_components: int = 10,
    regularization: float = 0.01,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
```

**Parameters:**

- `omics_data` (MultiOmicsData): Multi-omics data object with aligned samples
- `layer_pair` (Tuple[str, str]): Tuple of (layer1_name, layer2_name) to analyze. Both layers must exist in omics_data.layer_names
- `n_components` (int): Number of canonical variate pairs to extract. Must be positive and less than min(n_samples, n_features_X, n_features_Y). Default: 10
- `regularization` (float): Ridge regularization parameter for covariance matrices (prevents singular matrix issues). Should be positive and typically small (0.001-0.1). Larger values provide more stability but may bias results. Default: 0.01

**Returns:**
- Tuple containing:
  - `X_canonical`: Canonical variates for first layer (n_samples, n_components)
  - `Y_canonical`: Canonical variates for second layer (n_samples, n_components)
  - `X_weights`: Feature weights for first layer (n_features_X, n_components)
  - `Y_weights`: Feature weights for second layer (n_features_Y, n_components)
  - `correlations`: Canonical correlations for each component pair

**Raises:**
- `ValueError`: If layer names not found in omics_data or invalid parameters

### from_dna_variants

Convert DNA variant data (VCF) to genomics layer DataFrame.

```python
def from_dna_variants(
    vcf_path: Union[str, Path],
    sample_ids: Optional[List[str]] = None,
    variant_ids: Optional[List[str]] = None,
) -> pd.DataFrame:
```

**Parameters:**

- `vcf_path` (Union[str, Path]): Path to VCF file
- `sample_ids` (Optional[List[str]]): Optional list of sample IDs to include. If None, includes all samples.
- `variant_ids` (Optional[List[str]]): Optional list of variant IDs to include. If None, includes all variants.

**Returns:**
- DataFrame with samples as index, variants as columns, and genotype values (0, 1, 2 for homozygous ref, heterozygous, homozygous alt) as values.

**Raises:**
- `FileNotFoundError`: If VCF file does not exist
- `ValueError`: If VCF parsing fails
- `ImportError`: If DNA module not available

### from_rna_expression

Convert RNA expression data to transcriptomics layer DataFrame.

```python
def from_rna_expression(
    expression_path: Union[str, Path],
    sample_ids: Optional[List[str]] = None,
    gene_ids: Optional[List[str]] = None,
    transpose: bool = False,
) -> pd.DataFrame:
```

**Parameters:**

- `expression_path` (Union[str, Path]): Path to expression file (CSV or TSV)
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

```python
def from_protein_abundance(
    protein_path: Union[str, Path],
    sample_ids: Optional[List[str]] = None,
    protein_ids: Optional[List[str]] = None,
    transpose: bool = False,
) -> pd.DataFrame:
```

**Parameters:**

- `protein_path` (Union[str, Path]): Path to protein file (CSV, TSV, or FASTA)
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

The module uses structured logging via `metainformant.core.logging`:

- **INFO**: Major operations (loading data, starting analysis)
- **DEBUG**: Detailed information (file paths, data shapes)
- **WARNING**: Non-fatal issues (missing samples, data shifts)
- **ERROR**: Fatal errors (file read failures, computation errors)

## Examples

See [Module README](../src/metainformant/multiomics/README.md) and [Index Documentation](./index.md) for comprehensive usage examples.

