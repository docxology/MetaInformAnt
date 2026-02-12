# Spatial Autocorrelation

The spatial autocorrelation module implements classical spatial statistics for measuring spatial patterns in gene expression and other continuous variables. These statistics quantify whether observed values cluster spatially, disperse, or distribute randomly.

## Key Concepts

### Spatial Weights Matrix

All spatial statistics require a spatial weights matrix that defines neighborhood relationships between observations. The module provides `spatial_weights_matrix()` which constructs weights from coordinates using KNN, distance-band, or kernel methods.

### Global Statistics

Global statistics produce a single summary measure for the entire dataset:

- **Moran's I**: Range roughly [-1, 1]. Positive values indicate spatial clustering (similar values near each other). Near zero indicates random distribution. Negative values indicate dispersion (checkerboard pattern).
- **Geary's C**: Range [0, 2]. Values less than 1 indicate positive autocorrelation, 1 indicates no autocorrelation, values greater than 1 indicate negative autocorrelation. More sensitive to local differences than Moran's I.

### Local Statistics

Local statistics provide per-observation measures of spatial association:

- **Local Moran's I (LISA)**: Identifies local clusters and spatial outliers. Classifications: HH (high-high hot spot), LL (low-low cold spot), HL (high-low spatial outlier), LH (low-high spatial outlier), NS (not significant).
- **Getis-Ord G\***: Identifies statistically significant hot spots (high values clustered) and cold spots (low values clustered).

### Variograms

Spatial variograms characterize how spatial dissimilarity changes with distance. The semivariance at each lag distance reveals the scale of spatial dependence.

## Data Structures

### MoransIResult

```python
@dataclass
class MoransIResult:
    I: float              # Moran's I statistic
    expected_I: float     # Expected I under null
    variance_I: float     # Variance of I
    z_score: float        # Standardized Z-score
    p_value: float        # Two-sided p-value
    n: int                # Number of observations
```

### GearyCResult

```python
@dataclass
class GearyCResult:
    C: float; expected_C: float; variance_C: float
    z_score: float; p_value: float; n: int
```

### LocalMoransResult

```python
@dataclass
class LocalMoransResult:
    local_I: Any          # Per-observation local I values
    expected_I: float
    z_scores: Any         # Per-observation z-scores
    p_values: Any         # Per-observation p-values
    cluster_labels: Any   # "HH", "LL", "HL", "LH", "NS" per observation
    significance_level: float
```

### GetisOrdResult, VariogramResult

```python
@dataclass
class GetisOrdResult:
    G_star: Any; z_scores: Any; p_values: Any
    hot_spots: Any; cold_spots: Any; significance_level: float

@dataclass
class VariogramResult:
    lag_distances: Any; semivariances: Any; n_pairs: Any
    nugget: float; sill: float; range_param: float; model: str
```

## Function Reference

### spatial_weights_matrix

```python
def spatial_weights_matrix(
    coordinates: Any,
    method: str = "knn",
    k: int = 6,
    bandwidth: float | None = None,
) -> Any  # scipy sparse matrix
```

Construct a spatial weights matrix from coordinates. Methods: "knn" (k-nearest neighbors), "distance" (distance band), "kernel" (Gaussian kernel weights).

### morans_i

```python
def morans_i(values: Any, weights: Any) -> MoransIResult
```

Compute global Moran's I spatial autocorrelation statistic with significance testing via normal approximation.

### gearys_c

```python
def gearys_c(values: Any, weights: Any) -> GearyCResult
```

Compute global Geary's C spatial autocorrelation statistic.

### local_morans_i

```python
def local_morans_i(
    values: Any, weights: Any, significance_level: float = 0.05,
) -> LocalMoransResult
```

Compute Local Moran's I (LISA) for each observation with cluster classification.

### getis_ord_g

```python
def getis_ord_g(
    values: Any, weights: Any, significance_level: float = 0.05,
) -> GetisOrdResult
```

Compute Getis-Ord G\* statistic to identify hot spots and cold spots.

### spatial_variogram

```python
def spatial_variogram(
    values: Any, coordinates: Any, n_lags: int = 15,
) -> VariogramResult
```

Compute an empirical spatial variogram and fit a theoretical model. Returns semivariance at each lag distance plus fitted nugget, sill, and range parameters.

## Usage Examples

```python
from metainformant.spatial import io, analysis

dataset = io.load_visium("path/to/spaceranger_output/")

# Build spatial weights
weights = analysis.spatial_weights_matrix(dataset.coordinates, method="knn", k=6)

# Global Moran's I for a gene
gene_idx = dataset.gene_names.index("CD3E")
expression = dataset.expression[:, gene_idx].toarray().flatten()
result = analysis.morans_i(expression, weights)
print(f"Moran's I = {result.I:.4f} (p = {result.p_value:.4e})")

# Geary's C
geary = analysis.gearys_c(expression, weights)
print(f"Geary's C = {geary.C:.4f}")

# Local Moran's I (LISA)
lisa = analysis.local_morans_i(expression, weights, significance_level=0.05)
n_hotspots = sum(1 for l in lisa.cluster_labels if l == "HH")
print(f"Hot spots: {n_hotspots}")

# Getis-Ord G* for hot/cold spot detection
gstar = analysis.getis_ord_g(expression, weights)

# Spatial variogram
vario = analysis.spatial_variogram(expression, dataset.coordinates, n_lags=15)
print(f"Range: {vario.range_param:.1f}, Sill: {vario.sill:.4f}")
```

## Configuration

- **Environment prefix**: `SPATIAL_`
- **Optional dependencies**: numpy, scipy (KDTree, sparse, stats)
- P-values computed via normal approximation for global statistics
- LISA significance uses Bonferroni or FDR correction options

## Related Modules

- `spatial.analysis.clustering` -- Spatial clustering (complementary to autocorrelation analysis)
- `spatial.visualization` -- `plot_spatial_autocorrelation` for LISA maps
- `spatial.io` -- Data loading for spatial datasets
