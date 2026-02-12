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
| `compute_permutation_importance()` | Permutation-based feature importance scores |
| `compute_shap_values_kernel()` | Kernel SHAP values for model-agnostic explanations |
| `compute_lime_explanation()` | LIME local interpretable explanations |
| `feature_interaction()` | Detect pairwise feature interactions |
| `partial_dependence()` | Partial dependence plots for feature effects |
| `compute_attention_weights()` | Compute attention-like weights for feature attribution |
| `boruta_selection()` | Boruta all-relevant feature selection |
| `recursive_elimination()` | Recursive feature elimination with cross-validation |
| `stability_selection()` | Stability selection for robust feature identification |
| `mutual_information_selection()` | MI-based feature ranking and selection |

## Usage

```python
from metainformant.ml.interpretability import explainers, feature_selection

importance = explainers.compute_permutation_importance(model, X, y)
shap_vals = explainers.compute_shap_values_kernel(model.predict, X)
lime_exp = explainers.compute_lime_explanation(model.predict, X, instance_idx=0)
selected = feature_selection.boruta_selection(X, y, n_iterations=100)
stable = feature_selection.stability_selection(X, y, n_bootstrap=50)
```
