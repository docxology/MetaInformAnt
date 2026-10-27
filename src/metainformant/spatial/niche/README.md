# niche

Sub-package of `metainformant.spatial`. See the module's
`AGENTS.md`/`README.md` for the domain overview.

## Contents

| File | Purpose |
|------|---------|
| `identification.py` | `identify_niches` and the `NicheResult` dataclass |

`identify_niches(cell_type_proportions, coordinates, cell_type_names=None, n_niches=5, n_neighbors=15, spatial_weight=0.3, random_state=42)` smooths the cell-type proportions of each spot over its spatial neighbours and clusters the smoothed vectors with k-means. `NicheResult` holds `niche_labels`, `n_niches`, `compositions`, `cell_type_names`, `niche_sizes` and `niche_diversity`. The package `__init__` exports no names, so import from `metainformant.spatial.niche.identification`. See [docs/spatial/niche.md](../../../../docs/spatial/niche.md) and, for the label-based and deconvolution-based variants, `analysis.neighborhood.niche_detection` and `deconvolution.spatial_deconvolution.niche_identification`.
