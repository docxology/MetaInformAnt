### GWAS: Population Structure Analysis

Population structure analysis for GWAS, including PCA and kinship matrix computation.

**Functions**: `compute_pca`, `compute_kinship_matrix`, `estimate_population_structure`.

## Overview

Population structure analysis is essential for GWAS to:
- Detect population stratification
- Identify related samples
- Adjust for ancestry in association tests
- Quality control (outlier detection)

This module provides Principal Component Analysis (PCA) and kinship matrix
computation for genotype data.

## Functions

### `compute_pca(genotype_matrix, n_components=10)`

Compute Principal Component Analysis on genotype matrix.

**Parameters**:
- `genotype_matrix`: Genotype matrix (samples × variants), encoded as 0/1/2/-1
- `n_components`: Number of principal components to compute (default: 10)

**Returns**: Dictionary with:
- `pcs`: Principal components (samples × n_components)
- `explained_variance`: Explained variance per component
- `explained_variance_ratio`: Explained variance ratio per component
- `status`: "success" or "failed"

**Example**:
```python
from metainformant.gwas import compute_pca

# Genotype matrix: 5 samples, 10 variants
genotypes = [
    [0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
    [0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
    [2, 2, 1, 1, 2, 2, 1, 1, 2, 2],
    [2, 2, 1, 1, 2, 2, 1, 1, 2, 2],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
]

result = compute_pca(genotypes, n_components=3)
if result["status"] == "success":
    pcs = result["pcs"]
    variance_ratio = result["explained_variance_ratio"]
```

### `compute_kinship_matrix(genotype_matrix, method="vanraden")`

Compute kinship matrix from genotype data.

**Parameters**:
- `genotype_matrix`: Genotype matrix (samples × variants), encoded as 0/1/2/-1
- `method`: Kinship calculation method:
  - `"vanraden"` (default): VanRaden method (2008)
  - `"astle"`: Astle-Balding method
  - `"yang"`: Yang method (2011)

**Returns**: Dictionary with:
- `kinship_matrix`: Symmetric kinship matrix (samples × samples)
- `status`: "success" or "failed"

**Example**:
```python
from metainformant.gwas import compute_kinship_matrix

genotypes = [
    [0, 1, 2, 0, 1],
    [0, 1, 2, 0, 1],  # Identical to first
    [2, 1, 0, 2, 1],  # Different
]

result = compute_kinship_matrix(genotypes, method="vanraden")
if result["status"] == "success":
    kinship = result["kinship_matrix"]
    # High kinship between samples 0 and 1 (identical)
```

### `estimate_population_structure(vcf_path, config, output_dir=None)`

Estimate population structure from VCF file.

**Parameters**:
- `vcf_path`: Path to VCF file
- `config`: Dictionary with structure configuration:
  - `compute_pca`: Boolean (default: True)
  - `n_components`: Number of PCs (default: 10)
  - `compute_relatedness`: Boolean (default: True)
  - `kinship_method`: Method name (default: "vanraden")
- `output_dir`: Optional directory to write results

**Returns**: Dictionary with PCA and kinship results

**Example**:
```python
from metainformant.gwas import estimate_population_structure

config = {
    "compute_pca": True,
    "n_components": 3,
    "compute_relatedness": True,
    "kinship_method": "vanraden",
}

result = estimate_population_structure(
    vcf_path="variants.vcf",
    config=config,
    output_dir="output/structure"
)
```

## Usage Examples

### Basic PCA Analysis

```python
from metainformant.gwas import compute_pca
import numpy as np

# Load or create genotype matrix
genotypes = np.random.randint(0, 3, size=(100, 1000)).tolist()

# Compute PCA
result = compute_pca(genotypes, n_components=5)

if result["status"] == "success":
    pcs = np.array(result["pcs"])
    print(f"PC1 explains {result['explained_variance_ratio'][0]:.2%} of variance")
```

### Kinship Matrix for Relatedness

```python
from metainformant.gwas import compute_kinship_matrix

# Genotype matrix
genotypes = [
    [0, 1, 2, 0, 1, 2],
    [0, 1, 2, 0, 1, 2],  # Identical to sample 0
    [2, 1, 0, 2, 1, 0],  # Different
]

result = compute_kinship_matrix(genotypes, method="vanraden")
kinship = result["kinship_matrix"]

# High kinship indicates related samples
print(f"Kinship between samples 0 and 1: {kinship[0][1]:.4f}")
```

### Full Workflow from VCF

```python
from metainformant.gwas import estimate_population_structure

config = {
    "compute_pca": True,
    "n_components": 10,
    "compute_relatedness": True,
    "kinship_method": "vanraden",
}

result = estimate_population_structure(
    vcf_path="data/genotypes.vcf",
    config=config,
    output_dir="output/population_structure"
)

# Results saved to:
# - output/population_structure/pca_components.tsv
# - output/population_structure/kinship_matrix.tsv
# - output/population_structure/structure_summary.json
```

## Integration with GWAS Workflow

```python
from metainformant.gwas import execute_gwas_workflow, load_gwas_config

# Load GWAS configuration
config = load_gwas_config("config/gwas/analysis.yaml")

# Structure analysis is automatically included in workflow
results = execute_gwas_workflow(config)

# Access structure results
if "structure" in results:
    structure = results["structure"]
    if "pca" in structure:
        pcs = structure["pca"]["pcs"]
```

## Methods

### VanRaden Method (Default)

The VanRaden method (2008) computes kinship using:
- Centered genotype matrix
- Standardized by allele frequencies
- Suitable for dense SNP arrays

### Astle-Balding Method

Alternative method that may be more robust to:
- Missing data
- Rare variants
- Population structure

### Yang Method

Yang method (2011) from GCTA:
- Uses different normalization
- May be more appropriate for certain datasets

## Visualization

See `metainformant.gwas.visualization_population` for plotting functions:
- `pca_plot()`: 2D/3D PCA scatter plots
- `pca_scree_plot()`: Variance explained by each PC
- `kinship_heatmap()`: Kinship matrix visualization
- `admixture_plot()`: Ancestry proportions (requires ADMIXTURE output)
- `population_tree()`: Hierarchical clustering dendrogram

## Notes

1. **Missing Data**: Missing genotypes (encoded as -1) are imputed using mean
   imputation per variant before PCA computation.

2. **Data Requirements**: Requires genotype data (0/1/2 encoding) or VCF file.

3. **Computational Cost**: PCA scales with number of samples and variants.
   Consider LD pruning for large datasets.

## References

- Price, A. L., et al. (2006). Principal components analysis corrects for
  stratification in genome-wide association studies. *Nature Genetics*, 38(8), 904-909.
- VanRaden, P. M. (2008). Efficient methods to compute genomic predictions.
  *Journal of Dairy Science*, 91(11), 4414-4423.
- Yang, J., et al. (2011). GCTA: a tool for genome-wide complex trait analysis.
  *The American Journal of Human Genetics*, 88(1), 76-82.

