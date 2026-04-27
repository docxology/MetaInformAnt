# Spatial Autocorrelation

Quantify how gene expression clusters across tissue coordinates.

## Moran's I (global)

Range [-1, +1]: negative (dispersion), zero (random), positive (clustering).

```python
from metainformant.spatial.analysis.autocorrelation import morans_i
res = morans_i(adata, 'EPCAM', n_permutations=999)
print(res.statistic, res.pvalue, res.z_score)
```

Permutation test: shuffle expression across spots, recompute I, build null
distribution. Two‑tailed p = fraction |I_perm| ≥ |I_obs|. Default B=999 gives
resolution 0.1 %; increase to 9 999 for publication.

### Analytic (fast) variant

```python
res = morans_i(adata, 'VIM', n_permutations=0)   # uses normal approximation
```

Fast (≈5 ms) but inaccurate for N < 200 or extremely skewed distributions.

## Local Moran (LISA)

Spot‑level contribution to global Moran: reveals hot/cold spots.

```python
from metainformant.spatial.analysis.autocorrelation import local_moran
lisa = local_moran(adata, 'MKI67', n_permutations=99)
adata.obs['lisa_quad'] = lisa.quadrant  # 0=LL, 1=LH, 2=HL, 3=HH
adata.obs['lisa_sig'] = lisa.is_sig      # FDR‑adjusted
```

Quadrants:
- HH = high expression surrounded by high  (hotspot core)
- LL = low expression surrounded by low  (coldspot core)
- HL = high surrounded by low  (potential outlier)
- LH = low surrounded by high  (edge of hotspot)

Visualise:

```python
from metainformant.spatial.visualization import plot_spatial_categorical
plot_spatial_categorical(adata, 'lisa_quad', palette='RdBu')
```

## Geary's C

Focuses on dissimilarity between neighbours; range approx 0–2:

```python
from metainformant.spatial.analysis.autocorrelation import gearys_c
gc = gearys_c(adata, 'CD3D')
print(f"C = {gc.statistic:.3f} ( <1 → positive autocorrelation)")
```

No permutation by default; set `n_permutations=999` for p‑value.

## Getis‑Ord Gi*

Hot‑spot detection; more sensitive than LISA HH quadrant.

```python
from metainformant.spatial.analysis.autocorrelation import getis_ord
gi = getis_ord(adata, 'CXCL12')
hotspots = (gi.z_scores > 1.96) & (gi.pvalues < 0.05)
print(f"Hot spots: {hotspots.sum()} / {len(hotspots)}")
```

## Gene‑set pipeline

```python
import numpy as np
from statsmodels.stats.multitest import fdrcorrection
genes = adata.var_names[:300]
pvals = np.array([morans_i(adata, g, n_permutations=0).pvalue for g in genes])
reject, padj = fdrcorrection(pvals, alpha=0.05)
print(f"Autocorrelated (FDR<0.05): {reject.sum()} genes")
```

Use analytic (`n_permutations=0`) for screening; re‑run candidates with 999 perm
for publication figures.

## Edge correction

Pass `mask=adata.obs['in_tissue']` to exclude background (non‑tissue) spots from
the weight sums. All functions accept a `mask` keyword.

## Performance notes

| Statistic | N=5 000 spots | N=20 000 spots |
|-----------|---------------|----------------|
| Moran (analytic) | 8 ms | 35 ms |
| Moran (999 perm) | 340 ms | 1.4 s |
| Local Moran (99 perm) | 200 ms | 850 ms |

Permutation count dominates runtime. Parallelise across genes with `joblib` if
your machine has ≥8 cores.
