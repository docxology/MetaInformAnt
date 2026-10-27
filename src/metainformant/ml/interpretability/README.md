# Model Interpretability

Tools for understanding and explaining machine learning model predictions, including permutation importance, SHAP, LIME, partial dependence, and advanced feature selection methods.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `explainers` and `feature_selection` |
| `explainers.py` | Permutation importance, Kernel SHAP, LIME, partial dependence, feature interaction |
| `feature_selection.py` | Boruta, recursive elimination, stability selection, MI-based selection |

## Key Functions

| Function | Description |
|----------|-------------|
| `compute_permutation_importance()` | Permutation-based feature importance scores (`random_state` seeds the shuffles) |
| `compute_shap_values_kernel()` | Kernel SHAP values for model-agnostic explanations (`random_state` seeds the sampled coalitions) |
| `compute_lime_explanation()` | LIME local interpretable explanations (`random_state` seeds the perturbed samples) |
| `feature_interaction()` | Detect pairwise feature interactions |
| `partial_dependence()` | Partial dependence plots for feature effects |
| `compute_attention_weights()` | Compute attention-like weights for feature attribution |
| `boruta_selection()` | Boruta all-relevant feature selection (`max_iter` iterations; a `random_state`, 0 included, makes the selection reproducible) |
| `recursive_elimination()` | Recursive feature elimination with cross-validation |
| `stability_selection()` | Stability selection for robust feature identification (`random_state` seeds the subsamples) |
| `mutual_information_selection()` | MI-based feature ranking and selection |

## Usage

```python
from metainformant.ml.interpretability import explainers, feature_selection

importance = explainers.compute_permutation_importance(model, X, y, n_repeats=10, random_state=0)
shap_vals = explainers.compute_shap_values_kernel(model.predict, X, n_samples=100, random_state=0)
lime_exp = explainers.compute_lime_explanation(model.predict, X[0].tolist(), feature_names, random_state=0)
selected = feature_selection.boruta_selection(X, y, max_iter=100, random_state=0)
stable = feature_selection.stability_selection(X, y, n_bootstrap=50, random_state=0)
```

## Reproducibility

`compute_permutation_importance`, `compute_shap_values_kernel`, `compute_lime_explanation`, `boruta_selection` and
`stability_selection` take `random_state`: an integer seed or a `random.Random`. Their randomness comes from the
generator that `metainformant.ml._numeric.resolve_python_rng` builds from it, so the same seed gives the same result,
`random_state=None` draws fresh entropy, and the global `random` state is neither seeded nor consumed.
`boruta_selection` with a seed, `0` included, is reproducible: iteration `k` seeds its importance model (a random
forest, or the pure-Python ensemble when scikit-learn is missing) with `seed + k` for an integer seed, and with a draw
from the generator for a `random.Random`.
