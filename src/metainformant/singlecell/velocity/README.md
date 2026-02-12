# Single-Cell RNA Velocity

RNA velocity estimation from spliced/unspliced count matrices for inferring future transcriptional states. Supports steady-state and dynamical models, velocity embedding projection onto UMAP/tSNE, velocity-based pseudotime, and per-gene confidence metrics.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `rna_velocity` module |
| `rna_velocity.py` | Velocity computation, embedding, pseudotime, confidence, dynamical model |

## Key Functions

| Function | Description |
|----------|-------------|
| `compute_velocity()` | Estimate RNA velocity from spliced/unspliced counts |
| `velocity_embedding()` | Project velocity vectors onto low-dimensional embedding |
| `velocity_pseudotime()` | Derive temporal ordering from velocity transition probabilities |
| `velocity_confidence()` | Evaluate velocity reliability per gene and cell |
| `fit_dynamical_model()` | Fit full kinetic parameters (alpha, beta, gamma) per gene |

## Usage

```python
from metainformant.singlecell.velocity import rna_velocity

velocity = rna_velocity.compute_velocity(spliced, unspliced, model="steady_state")
embedding = rna_velocity.velocity_embedding(velocity, umap_coords)
pseudotime = rna_velocity.velocity_pseudotime(velocity, embedding)
```
