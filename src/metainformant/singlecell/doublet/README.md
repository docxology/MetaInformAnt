# Doublet Detection

Simulation-based doublet detection for single-cell RNA-seq data, inspired by
Scrublet and DoubletFinder.

## Contents

| File | Purpose |
|------|---------|
| `detection.py` | Synthetic-doublet generation, PCA embedding, KNN doublet scoring |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `detect_doublets()` | Generate synthetic doublets, score observed cells by KNN neighbor composition |
| `DoubletResult` | Dataclass with `scores`, `predicted_doublets`, `threshold`, `n_doublets`, `doublet_rate`, `synthetic_scores` |

## Usage

```python
from metainformant.singlecell.doublet.detection import detect_doublets

result = detect_doublets(
    counts,                 # raw counts, cells x genes
    expected_doublet_rate=0.06,
    n_neighbors=30,
    random_state=42,
)
predicted_mask = result.predicted_doublets
```

Scoring normalizes to log1p CPM, selects up to 2000 highly variable genes by
coefficient of variation, averages random cell pairs into synthetic doublets,
embeds observed + synthetic cells with SVD-based PCA, and thresholds scores at
the expected doublet rate.
