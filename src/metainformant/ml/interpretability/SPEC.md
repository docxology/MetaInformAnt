# Specification: interpretability

## 🎯 Scope
Model interpretability and explainability subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 3 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `explainers.py`
- `feature_selection.py`

### Behavior Contracts

#### `explainers.py`, `feature_selection.py`
- `compute_permutation_importance`, `compute_shap_values_kernel`, `compute_lime_explanation`, `boruta_selection` and
  `stability_selection` take `random_state: int | random.Random | None`. All their randomness (shuffles, sampled
  coalitions, perturbed samples, shadow features, subsamples) comes from `ml._numeric.resolve_python_rng(random_state)`:
  the same seed gives the same output, `None` draws fresh entropy, and the global `random` state is neither seeded nor
  consumed.
- `boruta_selection` with an integer seed (0 included) seeds the importance model of iteration `k` (a random forest, or
  the pure-Python ensemble when scikit-learn is missing) with `seed + k`; with a `random.Random` it draws each of those
  seeds from that generator; with `None` the seeds are fresh entropy.
- `compute_lime_explanation(predict_fn, instance, feature_names, n_samples=1000, n_features=10, kernel_width=None,
  random_state=None)` takes the instance as a list of feature values and the feature names as a required argument.
