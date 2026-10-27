# Specification: ml

## Scope
Machine learning module for METAINFORMANT.

## Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## Data Structures
- **Sub-packages**: automl, deep_learning, evaluation, features, interpretability, llm, models
- **Key Concepts**: AutoML, deep learning, LLM, feature selection, model evaluation

## API Definition
### Exports
- `__init__.py`

### Behavior Contracts

#### `_numeric.py`
- `resolve_python_rng(random_state=None) -> random.Random` is the generator behind the `random_state` argument of the
  pure-Python ML routines (`automl.optimization`, `interpretability.explainers`, `interpretability.feature_selection`).
  A `random.Random` is used as given, so one generator passed through several calls continues its stream; a NumPy
  `Generator` supplies the seed (`integers(0, 2**63 - 1)`); an integer gives a fresh `random.Random(seed)` (the same
  stream for the same seed); `None` gives a generator seeded from system entropy. The global `random` state is never
  seeded or consumed.

#### Seeds across the package
- Randomised routines take `random_state` and draw from a generator of their own. `random_state=None` means fresh
  entropy and leaves the global streams untouched (it does not reseed them).
- `evaluation.validation.compare_validation_strategies` supports `"holdout"`, `"<k>fold_cv"` and `"stratified_5fold"`;
  `models.classification.cross_validate_biological(method="ensemble", random_state=...)` seeds the folds and the
  ensemble, and its `parallel_jobs` governs the ensemble unless an explicit `n_jobs` keyword is given.
- `llm.ollama.config.OllamaConfig.max_retries` is the maximum number of attempts per request (the first included) and
  must be at least 1; the constructor raises `ValueError` otherwise.
