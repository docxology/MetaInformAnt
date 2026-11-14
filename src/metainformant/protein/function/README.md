# Protein Function

Protein function prediction from sequence properties including domain-based function mapping, subcellular localization, solubility estimation, physicochemical profiling, disorder prediction, active site identification, and post-translational modification site detection.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports the `prediction_core` and `prediction_analysis` modules |
| `prediction_core.py` | Domain-function mapping, localization, solubility, physicochemical properties |
| `prediction_analysis.py` | Disordered regions, active sites, PTM site prediction |

## Key Functions

| Function | Description |
|----------|-------------|
| `predict_function_from_domains()` | Map domain composition to predicted biological function |
| `predict_localization()` | Predict subcellular localization from signal motifs |
| `predict_solubility()` | Estimate protein solubility from sequence properties |
| `compute_physicochemical()` | Compute MW, pI, instability index, GRAVY, and composition |
| `predict_disordered_regions()` | Identify intrinsically disordered regions (IUPred-like) |
| `find_active_sites()` | Detect catalytic active sites from motif patterns |
| `predict_post_translational_mods()` | Predict phosphorylation, glycosylation, and other PTM sites |

## Usage

```python
from metainformant.protein.function import prediction_analysis, prediction_core

funcs = prediction_core.predict_function_from_domains(domain_list)
loc = prediction_core.predict_localization(sequence)
props = prediction_core.compute_physicochemical(sequence)
sites = prediction_analysis.find_active_sites(sequence)
```
