# Spatial Deconvolution

The deconvolution module estimates the proportion of different cell types within each spatial spot, which is especially relevant for Visium data where each 55um-diameter spot captures multiple cells. It implements Non-Negative Least Squares (NNLS) and NMF-based approaches, along with reference profile construction and niche identification.

## Key Concepts

### Why Deconvolution

Visium spots are larger than individual cells, so each spot contains a mixture of cell types. Deconvolution estimates the fractional composition of each cell type per spot using reference expression profiles from single-cell RNA-seq data.

### NNLS Deconvolution

Non-Negative Least Squares (NNLS) solves for cell type weights such that the observed spot expression is approximated as a non-negative linear combination of reference profiles. Weights are then normalized to sum to 1 per spot to yield cell type fractions.

### Reference Profile Construction

Reference profiles are built from scRNA-seq data by computing the mean or median expression per cell type. These profiles represent the expected gene expression signature of each cell type.

### Niche Identification

After deconvolution, spatial niches are identified by clustering spots based on their cell type composition, revealing microenvironments with characteristic cell type mixtures.

## Data Structures

### DeconvolutionResult

```python
@dataclass
class DeconvolutionResult:
    weights: Any              # Raw weights (n_spots x n_types)
    fractions: Any            # Normalized fractions (rows sum to 1)
    cell_type_names: list[str]
    residuals: Any            # Per-spot fitting residuals
    method: str               # "nnls" or "nmf"
    metadata: dict
```

Properties: `n_spots`, `n_types`.

## Function Reference

### create_reference_profiles

```python
def create_reference_profiles(
    scrna_data: Any,
    cell_type_labels: Any,
    *,
    gene_names: list[str] | None = None,
    method: str = "mean",
) -> tuple[Any, list[str], list[str]]
```

Build cell type reference expression profiles from scRNA-seq data. `method` can be "mean" or "median". Returns a tuple of (reference matrix, cell type names, gene names).

### nnls_deconvolution

```python
def nnls_deconvolution(
    spatial_expression: Any,
    reference_profiles: Any,
    cell_type_names: list[str],
) -> DeconvolutionResult
```

Perform NNLS-based deconvolution. Solves per-spot non-negative least squares against reference profiles.

### deconvolve_spots

```python
def deconvolve_spots(
    spatial_expression: Any,
    reference_profiles: Any,
    cell_type_names: list[str],
    method: str = "nnls",
) -> DeconvolutionResult
```

General deconvolution dispatcher supporting "nnls" and "nmf" methods.

### estimate_cell_fractions

```python
def estimate_cell_fractions(
    spatial_expression: Any,
    scrna_data: Any,
    cell_type_labels: Any,
) -> DeconvolutionResult
```

End-to-end pipeline: builds reference profiles from scRNA-seq data and deconvolves spatial spots in a single call.

### enrichment_score

```python
def enrichment_score(
    deconv_result: DeconvolutionResult,
    cell_type: str,
) -> Any  # numpy array
```

Compute per-spot enrichment score for a specific cell type relative to the overall mean.

### niche_identification

```python
def niche_identification(
    deconv_result: DeconvolutionResult,
    coordinates: Any,
    n_niches: int = 5,
) -> dict[str, Any]
```

Identify spatial niches by clustering spots based on cell type composition. Returns niche labels, composition profiles per niche, and spatial coherence scores.

### validate_deconvolution

```python
def validate_deconvolution(
    deconv_result: DeconvolutionResult,
    ground_truth: Any | None = None,
) -> dict[str, Any]
```

Validate deconvolution results. If ground truth is provided, computes correlation and RMSE. Otherwise computes internal quality metrics (residual distribution, fraction entropy).

## Usage Examples

```python
from metainformant.spatial import io, deconvolution

# Load spatial data
dataset = io.load_visium("path/to/spaceranger_output/")

# Build reference profiles from scRNA-seq
ref_profiles, type_names, genes = deconvolution.build_reference_profiles(
    scrna_expression, cell_type_labels, gene_names=gene_list
)

# Deconvolve spots
result = deconvolution.deconvolve_spots(
    dataset.expression, ref_profiles, type_names, method="nnls"
)
print(f"Deconvolved {result.n_spots} spots into {result.n_types} cell types")

# End-to-end pipeline
result = deconvolution.estimate_cell_fractions(
    dataset.expression, scrna_expression, cell_type_labels
)

# Identify spatial niches
niches = deconvolution.niche_identification(result, dataset.coordinates, n_niches=5)

# Validate results
metrics = deconvolution.validate_deconvolution(result)
print(f"Mean residual: {metrics['mean_residual']:.4f}")
```

## Configuration

- **Environment prefix**: `SPATIAL_`
- **Optional dependencies**: numpy, scipy (nnls), scikit-learn (NMF, KMeans, normalize)
- NNLS is parameter-free and generally preferred for its stability
- NMF requires specifying the number of components

## Related Modules

- `spatial.io` -- Data loading for spatial expression matrices
- `spatial.analysis.clustering` -- Spatial clustering (alternative to deconvolution-based niches)
- `spatial.visualization` -- `plot_deconvolution_pie`, `plot_cell_type_map`
- `spatial.integration` -- scRNA-seq integration for label transfer
