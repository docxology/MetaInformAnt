# Batch Effect Detection

Statistical detection of batch effects in sequencing feature matrices and
empirical-Bayes (ComBat-like) correction.

## Contents

| File | Purpose |
|------|---------|
| `detection.py` | PVCA-style batch variance decomposition, silhouette scoring, ComBat-like correction |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `BatchEffectReport` | Dataclass: sample/batch counts, PVCA variance split, silhouette score, severity |
| `detect_batch_effects()` | Per-feature F-tests against batch labels, variance decomposition, PCA-space silhouette |
| `correct_batch_combat()` | Adjust feature matrix to remove batch location/scale effects |

## Usage

```python
import numpy as np
from metainformant.quality.batch.detection import detect_batch_effects, correct_batch_combat

data = np.array(...)  # samples x features
report = detect_batch_effects(data, batch_labels=["A", "A", "B", "B"])
print(report.severity, report.pvca_variance, report.silhouette_score)

corrected = correct_batch_combat(data, batch_labels=["A", "A", "B", "B"])
```

`detect_batch_effects` requires scipy for exact F-distribution p-values;
`correct_batch_combat` requires numpy only.
