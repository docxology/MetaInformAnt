# Batch Effect Detection

Statistical detection of batch effects in sequencing feature matrices and
parametric empirical-Bayes ComBat correction.

## Contents

| File | Purpose |
|------|---------|
| `detection.py` | PVCA-style batch variance decomposition, silhouette scoring, empirical-Bayes ComBat correction |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `BatchEffectReport` | Dataclass: sample/batch counts, PVCA variance split, silhouette score, severity |
| `detect_batch_effects()` | Per-feature F-tests against batch labels, variance decomposition, PCA-space silhouette |
| `correct_batch_combat()` | Parametric empirical-Bayes ComBat (Johnson, Li & Rabinovic 2007, no covariates): removes batch location/scale effects |

## Usage

```python
import numpy as np
from metainformant.quality.batch.detection import detect_batch_effects, correct_batch_combat

data = np.array(...)  # samples x features
report = detect_batch_effects(data, batch_labels=["A", "A", "B", "B"])
print(report.severity, report.pvca_variance, report.silhouette_score)

corrected = correct_batch_combat(data, batch_labels=["A", "A", "B", "B"])
```

`detect_batch_effects` uses scipy for exact F-distribution p-values (the module imports scipy,
a core dependency, so both functions need numpy and scipy).

## `correct_batch_combat` behaviour

Data are standardised once, each batch's location and scale effects are estimated per feature and
shrunk by empirical Bayes hyperparameters fitted across features, and the posterior means are
iterated to convergence. The result does not depend on the order or the names of the batches or on
the order of the samples; the amount of shrinkage is estimated from the data, not fixed.

- `data` must be a 2D samples x features array with one batch label per sample, finite (impute or drop
  NaN/inf first), and every batch needs at least 2 samples; otherwise `ValueError`.
- A feature that is constant within some batch cannot have its variance estimated and is returned
  unchanged. The empirical-Bayes priors are fitted across the features that vary within every batch: if
  exactly one does, `ValueError` is raised; if none does, the data come back unchanged and a warning is
  logged.
- A single batch has nothing to correct and comes back as a copy; `data` itself is never modified.
- There are no covariates: biological groups that coincide with batches are removed with them.
