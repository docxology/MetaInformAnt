# Model Interpretability

The interpretability module provides model explanation methods for understanding which features drive predictions in biological machine learning models. It implements permutation importance, Kernel SHAP, LIME, partial dependence, feature interactions, and attention weights, with pure Python fallbacks where possible.

## Key Concepts

### Permutation Importance

Measures feature importance by randomly shuffling each feature and measuring the decrease in model performance. Features whose shuffling causes large performance drops are important. Model-agnostic -- works with any predictor.

### Kernel SHAP

Approximation of Shapley values using a weighted linear regression approach. SHAP values explain the contribution of each feature to each individual prediction, decomposing the prediction into a sum of feature effects.

### LIME (Local Interpretable Model-agnostic Explanations)

Explains individual predictions by fitting a local interpretable model (e.g., linear regression) around the point of interest. Perturbs the input and observes how predictions change.

### Partial Dependence

Shows the marginal effect of one or two features on the predicted outcome, averaging over all other features. Reveals the functional relationship between a feature and the model's predictions.

### Feature Interaction

Measures pairwise interaction strength between features using H-statistic or similar metrics. Identifies features that have synergistic effects on predictions.

## Function Reference

### compute_permutation_importance

```python
def compute_permutation_importance(
    model: Any,
    X: Any,
    y: Any,
    n_repeats: int = 10,
    scoring: Callable | None = None,
    random_state: int | None = None,
) -> dict[str, Any]
```

Compute permutation importance for all features. Returns a dictionary with `importances` (mean importance per feature), `importances_std`, `ranking` (feature indices sorted by importance), and raw `scores_per_feature`.

### compute_shap_values_kernel

```python
def compute_shap_values_kernel(
    model: Any,
    X: Any,
    background: Any | None = None,
    n_samples: int = 100,
) -> dict[str, Any]
```

Compute Kernel SHAP values for each sample and feature. Returns `shap_values` (matrix of shape n_samples x n_features), `expected_value` (base prediction), and `feature_effects` (mean absolute SHAP per feature).

### compute_lime_explanation

```python
def compute_lime_explanation(
    model: Any,
    X: Any,
    instance_index: int,
    n_features: int = 10,
    n_samples: int = 1000,
) -> dict[str, Any]
```

Generate a LIME explanation for a single instance. Returns `feature_weights` (linear model coefficients), `intercept`, `local_prediction`, `actual_prediction`, and `r_squared` (fidelity of the local model).

### partial_dependence

```python
def partial_dependence(
    model: Any,
    X: Any,
    feature_index: int,
    grid_resolution: int = 50,
) -> dict[str, Any]
```

Compute partial dependence of a single feature. Returns `feature_values` (grid of feature values) and `mean_predictions` (average prediction at each grid point).

### feature_interaction

```python
def feature_interaction(
    model: Any,
    X: Any,
    feature_pairs: list[tuple[int, int]] | None = None,
) -> dict[str, Any]
```

Measure pairwise feature interaction strength. Returns `interaction_scores` mapping feature index pairs to interaction strength, and `top_interactions` ranked by score.

### compute_attention_weights

```python
def compute_attention_weights(
    model: Any,
    X: Any,
) -> dict[str, Any]
```

Extract or compute attention-like weights showing which features the model focuses on for each prediction. Returns `attention_matrix` (n_samples x n_features) and `mean_attention` per feature.

### Feature Selection (interpretability-based)

```python
def boruta_selection(model, X, y, max_iter=100, random_state=None) -> dict
def recursive_elimination(model, X, y, n_features=None, cv=5) -> dict
def stability_selection(model, X, y, n_iterations=100, threshold=0.6) -> dict
def mutual_information_selection(X, y, k=20) -> dict
```

Feature selection methods based on interpretability. Boruta uses shadow features for statistical validation. Stability selection assesses feature selection robustness across subsamples.

## Usage Examples

```python
from metainformant.ml import (
    compute_permutation_importance, compute_shap_values_kernel,
    compute_lime_explanation, partial_dependence, feature_interaction,
    boruta_selection, stability_selection,
)

# Permutation importance
perm_imp = compute_permutation_importance(model, X_test, y_test, n_repeats=20)
for idx in perm_imp["ranking"][:10]:
    print(f"Feature {idx}: importance={perm_imp['importances'][idx]:.4f}")

# SHAP values
shap = compute_shap_values_kernel(model, X_test, background=X_train[:50])
print(f"Top feature effect: {max(shap['feature_effects']):.4f}")

# LIME explanation for a single sample
lime = compute_lime_explanation(model, X_test, instance_index=0, n_features=10)
for feat, weight in zip(range(10), lime["feature_weights"]):
    print(f"Feature {feat}: weight={weight:.4f}")

# Partial dependence plot data
pd_result = partial_dependence(model, X_test, feature_index=5)
# Plot pd_result["feature_values"] vs pd_result["mean_predictions"]

# Feature interactions
interactions = feature_interaction(model, X_test)
for pair, score in interactions["top_interactions"][:5]:
    print(f"Features {pair}: interaction={score:.4f}")

# Boruta feature selection
boruta = boruta_selection(model, X, y, max_iter=100)
selected = boruta["selected_features"]
```

## Configuration

- **Environment prefix**: `ML_`
- **Required**: numpy
- **Optional**: scipy (for statistical tests)
- All methods work with any model exposing `predict` or `predict_proba`
- Pure Python fallbacks are available when scikit-learn is not installed

## Related Modules

- `ml.models` -- Models to interpret
- `ml.features` -- Feature selection using statistical methods
- `ml.evaluation` -- Cross-validation for evaluating feature selection
- `ml.automl` -- Automated model selection
