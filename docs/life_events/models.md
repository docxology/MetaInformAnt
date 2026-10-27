# Life Events: Prediction Models

Models for learning from life event sequences: outcome prediction, multi-task prediction, survival analysis, and next-event prediction.

## Concept Overview

The models submodule provides these classes:

- **EventSequencePredictor** -- Outcome predictor with two model types: `"embedding"` (mean event embedding fed to a logistic or linear model) and `"simple"` (mean outcome per sequence length)
- **MultiTaskPredictor** -- Several outcomes from one shared sequence representation, with a ridge head per regression task and a logistic head per binary classification task
- **EnsemblePredictor** -- Weighted average of fitted `EventSequencePredictor` models
- **SurvivalPredictor** -- Time-to-event prediction with a Cox, exponential, or Kaplan-Meier model (NumPy)
- **LSTMSequenceModel** / **GRUSequenceModel** -- Recurrent next-event models (require PyTorch; they do not predict outcomes)

Every class except `EnsemblePredictor` has a `fit` method, and every class has `predict`. Only `EventSequencePredictor` can be saved and reloaded (`save_model` / `load_model`, JSON). Sequences are passed as `EventSequence` objects or as lists of event tokens; `convert_sequences_to_tokens` (see [Events](events.md)) produces `"domain:event_type"` tokens.

## Model Reference

### EventSequencePredictor

The entry point for outcome prediction.

```python
from metainformant.life_events.models.predictor import EventSequencePredictor

predictor = EventSequencePredictor(
    model_type="embedding",   # "embedding" or "simple"
    embedding_dim=100,
    random_seed=42,
)
predictor.fit(sequences, outcomes, task="classification")
predictions = predictor.predict(test_sequences)
probabilities = predictor.predict_proba(test_sequences)
```

**Methods:** `fit(sequences, outcomes, task=None, event_embeddings=None)`, `predict(sequences)`, `predict_proba(sequences)`, `save_model(path)`, `load_model(path)` (class method). `task` is `"classification"` or `"regression"`; when omitted, the constructor's `task_type` is used, then `"classification"`.

**Model types:**

- `"embedding"` -- represents a sequence by the mean of its events' embedding vectors and fits a logistic regression (classification, 0/1 outcomes) or a linear regression (regression) on it. Embeddings come from `event_embeddings` (a token-to-vector dictionary); events without a vector, and every event when none are given, get random vectors drawn from a normal distribution with standard deviation 0.1.
- `"simple"` -- predicts the mean training outcome of the sequences with the same number of events (the nearest observed length for unseen lengths, the shorter one on a tie). Classification is binary: the score is the probability of the larger label.
- Any other `model_type`, including `"lstm"` and `"gru"`, makes `fit` raise `ValueError`. Recurrent models predict the next event rather than an outcome; see below.

**Predictions:** `predict` returns class labels for classification (the larger label when its probability is at least 0.5, otherwise the smaller label) and values for regression. `predict_proba` is classification only and returns an `(n_sequences, 2)` array whose columns are the probabilities of the smaller and the larger label.

To supply learned embeddings, pass token lists for both steps and set `embedding_dim` to the length of the vectors (a mismatch raises `ValueError` whenever an event has no supplied vector, at `fit` or at `predict`). `EventSequencePredictor` identifies the events of an `EventSequence` by `event_type` alone, whereas `convert_sequences_to_tokens` and `learn_event_embeddings` use `"domain:event_type"` tokens, so embeddings learned from tokens match token lists but not `EventSequence` objects.

```python
from metainformant.life_events.core.utils import convert_sequences_to_tokens
from metainformant.life_events.models.embeddings import learn_event_embeddings

tokens = convert_sequences_to_tokens(sequences)
embeddings = learn_event_embeddings(tokens, embedding_dim=64, random_state=0)

predictor = EventSequencePredictor(model_type="embedding", embedding_dim=64, task_type="classification")
predictor.fit(tokens, outcomes, event_embeddings=embeddings)
```

### LSTMSequenceModel / GRUSequenceModel

Recurrent next-event models. They learn to predict the last event of each sequence from the events before it, and they need PyTorch: without it `fit` raises `ImportError`.

```python
from metainformant.life_events.models.sequence_models import LSTMSequenceModel, GRUSequenceModel

model = LSTMSequenceModel(
    embedding_dim=100, hidden_dim=64, num_layers=1,
    dropout=0.1, epochs=50, batch_size=32, learning_rate=0.001,
    random_state=0,
)
model.fit(sequences)
next_event_probabilities = model.predict(test_sequences)
```

`predict` returns one list of probabilities per sequence, one entry per event token in `model.event_vocab`. The optional `targets` argument of `fit` is ignored (with a logged warning), so outcome labels have no effect; use `EventSequencePredictor` or `MultiTaskPredictor` to predict outcomes. `predict` raises `ValueError` before the model is fitted. The recurrent models cannot be saved with `save_model`.

### MultiTaskPredictor

Predicts several outcomes at once from one shared representation of each sequence.

```python
from metainformant.life_events.models.statistical_models import MultiTaskPredictor

predictor = MultiTaskPredictor(
    task_types={"employment": "classification", "income": "regression"},
    embedding_dim=100, hidden_dim=64, random_state=0,
)
predictor.fit(sequences, {"employment": employment_labels, "income": income_values})
results = predictor.predict(test_sequences)                      # {"employment": [...], "income": [...]}
income = predictor.predict(test_sequences, task_name="income")   # one task only
```

- **Shared representation:** word2vec event embeddings of size `embedding_dim` are learned from the training sequences (or passed to `fit` as `event_embeddings`). A sequence is described by its standardised event counts together with the mean of its event embeddings, the latter compressed to at most `hidden_dim` principal components. `shared_representation(sequences)` returns this matrix.
- **Task heads:** regression tasks use ridge regression; classification tasks use L2-regularised logistic regression and must have exactly two distinct labels (otherwise `fit` raises `ValueError`). The larger label is the positive class, and classification predictions are its probabilities. `alpha` (default 1.0) sets the penalty of every head.
- **Reproducibility:** the ridge heads are solved in closed form and the logistic heads by Newton-Raphson, so predictions are deterministic given the embeddings. `random_state` seeds the embedding training; without it the embeddings, and the predictions, differ slightly between fits.
- Options other than those in the signature are ignored with a logged warning.

### EnsemblePredictor

Combines fitted `EventSequencePredictor` models with optional weights.

```python
from metainformant.life_events.models.predictor import EnsemblePredictor

ensemble = EnsemblePredictor(
    models=[predictor1, predictor2, predictor3],
    weights=[0.5, 0.3, 0.2]   # Optional; uniform if None
)
predictions = ensemble.predict(test_sequences)
```

`predict` returns the weighted average of the members' `predict` outputs and raises `ValueError` if any member is not fitted. Class labels are averaged as numbers, so a classification ensemble returns values between the labels rather than labels.

### SurvivalPredictor

Time-to-event prediction for right-censored data. The covariates of a sequence are the numbers of times each `"domain:event_type"` token occurs in it.

```python
from metainformant.life_events.models.statistical_models import SurvivalPredictor

surv = SurvivalPredictor(method="cox", penalizer=0.1)
surv.fit(sequences, event_times, event_occurred)   # event_occurred: True = event, False = censored
mean_survival = surv.predict(test_sequences)        # restricted mean survival time per sequence
survival_curves = surv.predict_survival_function(test_sequences, times=[1, 5, 10])   # shape (n_sequences, 3)
log_hazard_ratios = surv.coefficients               # {"domain:event_type": coefficient}
```

| `method` | Model |
|----------|-------|
| `"cox"` (default) | Cox proportional hazards: Newton-Raphson with step halving on the Breslow partial likelihood (tied times share a risk set); baseline hazard by Breslow's estimator |
| `"exponential"` | Exponential proportional hazards model (constant hazard), fitted by maximum likelihood |
| `"baseline"` (alias `"kaplan_meier"`) | Kaplan-Meier estimate, the same for every sequence |

- `predict` returns the restricted mean survival time: the area under each sequence's survival function up to the longest follow-up time in the training data. `predict_survival_function` returns `S(t | x)` for the requested times; the Cox and Kaplan-Meier estimates are step functions that stay constant after the last event time.
- `penalizer` (default 0.1) is a ridge penalty on the Cox and exponential coefficients; 0 gives the unpenalised maximum-likelihood fit. `max_iter` (default 50) and `tol` (default 1e-9) control Newton-Raphson; a fit that does not converge emits a `RuntimeWarning`.
- `coefficients` (log hazard ratio per additional occurrence of each token) and `standard_errors` are empty for `"baseline"`.
- `embedding_dim` and `random_state` are accepted for compatibility and have no effect: the fits use event counts and are deterministic.
- `fit` raises `ValueError` for empty input, inputs of different lengths, negative or non-finite times, an `event_occurred` that is not boolean or 0/1, or data without any event. An unknown `method` raises `ValueError` at construction, and `predict` raises `ValueError` before the model is fitted.

For a stand-alone Cox fit on covariate vectors, see `cox_ph_model` in [Survival](survival.md).

## Embedding Functions

Standalone functions for learning event vector representations.

| Function | Signature | Purpose |
|----------|-----------|---------|
| `learn_event_embeddings` | `(sequences, embedding_dim=100, window_size=5, min_count=1, sg=1, epochs=5, method=None, random_state=None, verbose=False, **kwargs)` | word2vec event embeddings (skip-gram or CBOW, negative sampling), trained in NumPy |
| `biological_embedding` | `(sequences, embedding_type="event_type", random_state=None, **kwargs)` | Heuristic embeddings by `embedding_type` (`"event_type"`, `"domain"`, `"temporal"`); an unknown type raises `ValueError` |
| `domain_specific_embeddings` | `(sequences, domains, embedding_dim=50, **kwargs)` | One `learn_event_embeddings` model per life domain |

```python
from metainformant.life_events.models.embeddings import learn_event_embeddings

embeddings = learn_event_embeddings(sequences, embedding_dim=100, window_size=5, epochs=5, random_state=0)
# embeddings["career:promotion"] -> np.ndarray of shape (100,)
```

`learn_event_embeddings` returns a dictionary that maps each retained event token to its vector. Sequences may be `EventSequence` objects (tokenised as `"domain:event_type"`) or token lists; co-occurrence is counted inside each sequence only.

| Option | Default | Meaning |
|--------|---------|---------|
| `method` | from `sg` | `"skipgram"` or `"cbow"`; takes precedence over `sg` (`sg=1` is skip-gram, `sg=0` is CBOW) |
| `window_size` | `5` | Maximum distance between an event and its context events |
| `min_count` | `1` | Events with fewer occurrences are dropped before training and are not returned |
| `epochs` | `5` | Passes over the sequences |
| `negative` | `5` | Negative samples per positive example |
| `learning_rate` (alias `alpha`) | `0.025` | Initial learning rate, decaying linearly to `min_learning_rate` (alias `min_alpha`, default `1e-4`) |
| `ns_exponent` | `0.75` | Exponent of the noise distribution used for negative sampling |
| `cbow_mean` | `True` | CBOW only: average rather than sum the context vectors |
| `combine` | `"sum"` | `"sum"` returns each event's input vector plus its output (context) vector; `"input"` returns the plain input vectors |
| `random_state` (alias `seed`) | `None` | Seed; identical seeds give identical embeddings |

The keyword options (`negative`, `learning_rate`, `min_learning_rate`, `ns_exponent`, `cbow_mean`, `combine`, and the aliases) go through `**kwargs`; any other keyword is ignored with a logged warning. `ValueError` is raised for empty input, a non-positive `embedding_dim`, `window_size`, or `epochs`, `negative` below 1, a non-positive `learning_rate`, a negative `min_learning_rate`, an unknown `method` or `combine`, and when no event reaches `min_count`.

`biological_embedding` returns a dictionary with the vectors under `"embeddings"`. The `"domain"` embeddings contain a small random variation around each domain's centre that is reproducible: it is drawn from `random_state` (alias `seed`) and from a fixed default seed when none is given. The `"event_type"` embeddings involve no randomness.

The `"temporal"` embeddings work for any `EventSequence` and involve no randomness either. Sequences with at least two events are grouped by their first three event types, the intervals between consecutive events (in days; datetime timestamps are converted, numeric timestamps are subtracted as given) of a group are pooled, and the group's vector holds their mean, standard deviation, minimum, maximum, count and median, zero padded to `embedding_dim` (default 30) and scaled so that its largest absolute entry is 1. The result's `"embeddings"` are keyed by the string of the three-type tuple. Because the vector holds six statistics, an `embedding_dim` below 6 raises `ValueError`.

## Typical Workflow

```python
from pathlib import Path

from metainformant.life_events.core.utils import convert_sequences_to_tokens, generate_cohort_sequences
from metainformant.life_events.models.embeddings import learn_event_embeddings
from metainformant.life_events.models.predictor import EventSequencePredictor

sequences = generate_cohort_sequences(n_sequences=500, avg_events_per_sequence=15, random_state=0)
tokens = convert_sequences_to_tokens(sequences)
outcomes = [1 if len(s) > 15 else 0 for s in tokens]    # 1 when a person has more than 15 events

train, test = tokens[:400], tokens[400:]
embeddings = learn_event_embeddings(train, embedding_dim=64, random_state=0)

model = EventSequencePredictor(model_type="embedding", embedding_dim=64, random_seed=0)
model.fit(train, outcomes[:400], task="classification", event_embeddings=embeddings)
predictions = model.predict(test)

Path("output/life_events").mkdir(parents=True, exist_ok=True)   # save_model does not create directories
model.save_model("output/life_events/trained_model.json")
reloaded = EventSequencePredictor.load_model("output/life_events/trained_model.json")
```

## Configuration

`analyze_life_course` takes its embedding and model settings from the `embedding:` and `model:` sections of the workflow configuration (`embedding_dim`, `window_size`, `epochs`, `method`, `learning_rate`, `model_type`, `task_type`, `random_state`, and others). The keys, defaults, and the environment variables that override a loaded configuration are listed under [Workflow](workflow.md#configuration).

## See Also

- [Events](events.md) -- Input data structures
- [Survival](survival.md) -- Dedicated survival analysis
- [Workflow](workflow.md) -- End-to-end analysis pipelines
