"""Statistical and multi-task prediction models for life events.

Every model in this module is implemented in NumPy (with SciPy for the logistic function); none of them
needs PyTorch.

- MultiTaskPredictor: simultaneous prediction of several outcomes. Every task has its own head, ridge
  regression for regression tasks and L2-regularised logistic regression for binary classification tasks,
  on one representation of the sequences that all tasks share.
- SurvivalPredictor: survival analysis with a Cox proportional hazards model (Newton-Raphson on the Breslow
  partial likelihood, right-censoring supported), an exponential hazard model, or Kaplan-Meier estimates.
- attention_weights: Extract attention weights from trained models
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Self, Tuple, Union

import numpy as np
from scipy.special import expit

from metainformant.core.utils import logging
from metainformant.life_events.core.utils import _sequence_to_tokens, sequence_embeddings
from metainformant.life_events.models.embeddings import learn_event_embeddings

if TYPE_CHECKING:
    from metainformant.life_events.core.events import EventSequence

logger = logging.get_logger(__name__)

_TASK_TYPES = ("classification", "regression")
# Survival method name -> the model that implements it ("baseline" is the Kaplan-Meier estimator).
_SURVIVAL_METHODS = {"cox": "cox", "exponential": "exponential", "baseline": "baseline", "kaplan_meier": "baseline"}
_MIN_SCALE = 1e-12  # standard deviations below this are treated as zero (constant features)
_RANK_TOLERANCE = 1e-8  # singular values below this fraction of the largest are treated as zero
_SINGULARITY_TOLERANCE = 1e-12  # eigenvalues below this fraction of the largest are numerical noise
_NULL_WEIGHT_TOLERANCE = 1e-8  # a parameter with more weight than this in a null direction is not identified
_MAX_LOG_HAZARD = 500.0  # linear predictors are clipped to +-this before exponentiation


# ---------------------------------------------------------------------------
# Shared numerics
# ---------------------------------------------------------------------------


def _event_vocabulary(sequences: List[Any]) -> Dict[str, int]:
    """Map every event token in ``sequences`` to a column index (sorted, so the mapping is reproducible)."""
    tokens = sorted({token for sequence in sequences for token in _sequence_to_tokens(sequence)})
    return {token: index for index, token in enumerate(tokens)}


def _count_matrix(sequences: List[Any], vocabulary: Dict[str, int]) -> np.ndarray:
    """Number of occurrences of every vocabulary event in every sequence; unknown events are ignored."""
    counts = np.zeros((len(sequences), len(vocabulary)))
    for row, sequence in enumerate(sequences):
        for token in _sequence_to_tokens(sequence):
            column = vocabulary.get(token)
            if column is not None:
                counts[row, column] += 1.0
    return counts


def _column_scale(values: np.ndarray) -> np.ndarray:
    """Per-column standard deviation, with 1 for constant columns so they standardise to zero."""
    scale = np.asarray(values.std(axis=0), dtype=float)
    return np.where(scale > _MIN_SCALE, scale, 1.0)


def _stable_inverse(matrix: np.ndarray) -> np.ndarray:
    """Inverse of a symmetric positive semi-definite matrix, ignoring its (numerically) null directions.

    The eigen-directions whose eigenvalue is below ``_SINGULARITY_TOLERANCE`` times the largest (rounding
    noise, several orders of magnitude above machine precision) are dropped, so an exactly singular matrix
    (for example collinear covariates without a penalty) yields the minimum-norm solution instead of
    coefficients that blow up.
    """
    eigenvalues, eigenvectors = np.linalg.eigh((matrix + matrix.T) / 2.0)
    keep = eigenvalues > _SINGULARITY_TOLERANCE * max(float(eigenvalues[-1]), np.finfo(float).tiny)
    return np.asarray((eigenvectors[:, keep] / eigenvalues[keep]) @ eigenvectors[:, keep].T, dtype=float)


def _standard_errors(curvature: np.ndarray, first: int = 0) -> np.ndarray:
    """Standard errors ``sqrt(diag(curvature^-1))`` of the parameters from index ``first`` on.

    ``curvature`` is the information matrix (negative Hessian of the log-likelihood). A parameter that is not
    identified, because it has weight in a null direction of the curvature (one of two collinear covariates,
    a constant covariate, no event at all), has no finite standard error and gets ``inf``. The others come
    from the generalised inverse, which is exact for identified parameters.
    """
    eigenvalues, eigenvectors = np.linalg.eigh((curvature + curvature.T) / 2.0)
    keep = eigenvalues > _SINGULARITY_TOLERANCE * max(float(eigenvalues[-1]), np.finfo(float).tiny)
    covariance = (eigenvectors[:, keep] / eigenvalues[keep]) @ eigenvectors[:, keep].T
    identified = np.all(np.abs(eigenvectors[:, ~keep]) < _NULL_WEIGHT_TOLERANCE, axis=1)
    variances = np.clip(np.diag(covariance), 0.0, None)
    return np.asarray(np.where(identified, np.sqrt(variances), np.inf)[first:], dtype=float)


def _newton_direction(curvature: np.ndarray, gradient: np.ndarray) -> np.ndarray:
    """Newton step ``curvature^-1 @ gradient`` (minimum-norm when the curvature is singular)."""
    return np.asarray(_stable_inverse(curvature) @ gradient, dtype=float)


def _maximise(
    evaluate: Callable[[np.ndarray], Tuple[float, np.ndarray, np.ndarray]],
    start: np.ndarray,
    max_iter: int,
    tol: float,
) -> Tuple[np.ndarray, float, int, bool]:
    """Maximise a concave objective with Newton-Raphson steps that are halved until they do not hurt.

    ``evaluate(params)`` returns the objective, its gradient and the negative Hessian (positive
    semi-definite). A step is accepted when the objective stays finite and does not decrease, so a step
    that overshoots (or overflows) is shortened instead of ruining the fit. Iteration stops when a full
    Newton step improves the objective by less than ``tol`` (relative), or when no step improves it any more
    and the gradient has vanished.

    Returns:
        ``(params, objective, n_iterations, converged)``
    """
    params = np.asarray(start, dtype=float)
    value, gradient, curvature = evaluate(params)
    for iteration in range(1, max_iter + 1):
        direction = _newton_direction(curvature, gradient)
        scale = 1.0
        while True:
            candidate = params + scale * direction
            new_value, new_gradient, new_curvature = evaluate(candidate)
            if np.isfinite(new_value) and new_value >= value - 1e-12 * (1.0 + abs(value)):
                break
            scale *= 0.5
            if scale < 1e-10:
                vanished = float(np.linalg.norm(gradient)) <= 1e-6 * (1.0 + abs(value))
                return params, value, iteration, vanished
        improvement = new_value - value
        params, value, gradient, curvature = candidate, new_value, new_gradient, new_curvature
        if scale == 1.0 and improvement <= tol * (1.0 + abs(value)):
            return params, value, iteration, True
    return params, value, max_iter, False


def _fit_ridge_head(features: np.ndarray, target: np.ndarray, alpha: float) -> Tuple[np.ndarray, float]:
    """Ridge regression with an unpenalised intercept; returns ``(weights, intercept)``."""
    feature_mean = features.mean(axis=0)
    target_mean = float(target.mean())
    centred = features - feature_mean
    gram = centred.T @ centred + alpha * np.eye(features.shape[1])
    weights = np.linalg.solve(gram, centred.T @ (target - target_mean))
    return weights, target_mean - float(feature_mean @ weights)


def _fit_logistic_head(features: np.ndarray, labels: np.ndarray, alpha: float) -> Tuple[np.ndarray, float]:
    """L2-regularised logistic regression (unpenalised intercept) fitted by Newton-Raphson.

    ``labels`` holds 0/1 values. The penalty makes the objective strictly concave, so the fit is unique and
    finite even when the classes are perfectly separable.
    """
    n_samples, n_features = features.shape
    design = np.column_stack([features, np.ones(n_samples)])
    penalty = np.append(np.full(n_features, alpha), 0.0)
    prior = float(np.clip(labels.mean(), 1e-6, 1.0 - 1e-6))
    start = np.zeros(n_features + 1)
    start[-1] = np.log(prior / (1.0 - prior))

    def evaluate(params: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
        linear = design @ params
        log_likelihood = float(np.sum(labels * linear - np.logaddexp(0.0, linear)))
        probability = expit(linear)
        gradient = design.T @ (labels - probability) - penalty * params
        curvature = (design.T * (probability * (1.0 - probability))) @ design + np.diag(penalty)
        return log_likelihood - 0.5 * float(penalty @ (params * params)), gradient, curvature

    params, _, _, converged = _maximise(evaluate, start, max_iter=100, tol=1e-10)
    if not converged:
        logger.warning("Logistic regression head did not converge")
    return params[:-1], float(params[-1])


# ---------------------------------------------------------------------------
# Survival numerics
# ---------------------------------------------------------------------------


def _validate_survival_data(n_sequences: int, event_times: Any, event_occurred: Any) -> Tuple[np.ndarray, np.ndarray]:
    """Check follow-up times and event indicators; returns ``(times, occurred)`` as float and bool arrays."""
    if n_sequences == 0:
        raise ValueError("sequences cannot be empty")
    times = np.asarray(event_times, dtype=float)
    indicators = np.asarray(event_occurred)
    if times.ndim != 1 or indicators.ndim != 1:
        raise ValueError("event_times and event_occurred must be one-dimensional")
    if times.shape[0] != n_sequences or indicators.shape[0] != n_sequences:
        raise ValueError(
            f"sequences, event_times and event_occurred must have the same length "
            f"({n_sequences}, {times.shape[0]}, {indicators.shape[0]})"
        )
    if not np.all(np.isfinite(times)) or np.any(times < 0):
        raise ValueError("event_times must be finite and not negative")
    if indicators.dtype == bool:
        occurred = indicators
    else:
        if not np.all(np.isin(indicators, (0, 1))):
            raise ValueError("event_occurred must contain only booleans or 0 and 1 (1 = event, 0 = censored)")
        occurred = indicators.astype(bool)
    if not occurred.any():
        raise ValueError(
            "At least one event must have occurred; a survival model cannot be fitted to censored data only"
        )
    return times, occurred


class _CoxPartialLikelihood:
    """Breslow partial likelihood of right-censored data and its derivatives.

    For subjects with linear predictor ``eta_i = x_i . beta`` and event times ``u``, the partial
    log-likelihood is ``sum_i d_i eta_i - sum_u deaths(u) log sum_{j: t_j >= u} exp(eta_j)``. Tied event
    times share one risk set (Breslow's approximation). Everything is computed with reverse cumulative sums
    over the time-sorted subjects, so one evaluation costs O(n p^2) time and O(n p) memory.
    """

    def __init__(self, features: np.ndarray, times: np.ndarray, occurred: np.ndarray):
        order = np.argsort(times, kind="stable")
        self.features = features[order]
        self.times = times[order]
        self.occurred = occurred[order]
        died = self.times[self.occurred]
        self.event_times = np.unique(died)
        # Index of the first subject still at risk at each event time (the risk set is everything after).
        self.risk_start = np.searchsorted(self.times, self.event_times, side="left")
        self.deaths = np.bincount(np.searchsorted(self.event_times, died), minlength=self.event_times.size).astype(
            float
        )
        # Number of distinct event times up to and including each subject's own time.
        self.event_times_so_far = np.searchsorted(self.event_times, self.times, side="right")

    def _risk_sums(self, eta: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
        """Shift, and the risk-set sums of exp(eta - shift) per event time and of the covariate-weighted terms."""
        shift = float(eta.max())
        weights = np.exp(eta - shift)
        weighted = np.cumsum((self.features * weights[:, None])[::-1], axis=0)[::-1]
        totals = np.cumsum(weights[::-1])[::-1]
        return shift, totals[self.risk_start], weighted[self.risk_start]

    def evaluate(self, beta: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
        """Partial log-likelihood, its gradient and the information matrix (negative Hessian) at ``beta``.

        A ``beta`` that overflows the likelihood is reported as ``-inf`` so a line search can reject it.
        """
        n_features = self.features.shape[1]
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            eta = self.features @ beta
            shift, risk_total, risk_weighted = self._risk_sums(eta)
            log_likelihood = float(self.occurred @ eta - self.deaths @ (np.log(risk_total) + shift))
            # Breslow baseline-hazard increments and the cumulative hazard each subject has been exposed to.
            cumulative = np.concatenate(([0.0], np.cumsum(self.deaths / risk_total)))
            exposure = np.exp(eta - shift) * cumulative[self.event_times_so_far]
            gradient = self.features.T @ (self.occurred - exposure)
            means = risk_weighted / risk_total[:, None]
            information = (self.features.T * exposure) @ self.features - (means.T * self.deaths) @ means
        if not (np.isfinite(log_likelihood) and np.all(np.isfinite(gradient)) and np.all(np.isfinite(information))):
            return -np.inf, np.zeros(n_features), np.zeros((n_features, n_features))
        return log_likelihood, gradient, (information + information.T) / 2.0

    def baseline_cumulative_hazard(self, beta: np.ndarray) -> np.ndarray:
        """Breslow estimate of the cumulative baseline hazard at every distinct event time."""
        with np.errstate(over="ignore", under="ignore"):
            shift, risk_total, _ = self._risk_sums(self.features @ beta)
            return np.asarray(np.cumsum(self.deaths / risk_total) * np.exp(-shift), dtype=float)


def _fit_cox(
    features: np.ndarray, times: np.ndarray, occurred: np.ndarray, penalizer: float, max_iter: int, tol: float
) -> Dict[str, Any]:
    """Cox proportional hazards model by Newton-Raphson on the (ridge-penalised) Breslow partial likelihood.

    The covariates are centred first (the partial likelihood does not depend on the centring); the
    baseline hazard then refers to a sequence with the average event counts of the training data.
    """
    n_features = features.shape[1]
    feature_mean = features.mean(axis=0)
    data = _CoxPartialLikelihood(features - feature_mean, times, occurred)
    identity = np.eye(n_features)

    def evaluate(beta: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
        log_likelihood, gradient, information = data.evaluate(beta)
        value = log_likelihood - 0.5 * penalizer * float(beta @ beta)
        return value, gradient - penalizer * beta, information + penalizer * identity

    beta, _, n_iter, converged = _maximise(evaluate, np.zeros(n_features), max_iter, tol)
    log_likelihood, _, information = data.evaluate(beta)
    return {
        "coefficients": beta,
        "std_errors": _standard_errors(information + penalizer * identity),
        "information": information,
        "feature_means": feature_mean,
        "baseline_times": data.event_times,
        "baseline_cumulative_hazard": data.baseline_cumulative_hazard(beta),
        "log_likelihood": log_likelihood,
        "n_iter": n_iter,
        "converged": converged,
    }


def _fit_exponential(
    features: np.ndarray, times: np.ndarray, occurred: np.ndarray, penalizer: float, max_iter: int, tol: float
) -> Dict[str, Any]:
    """Exponential proportional hazards model ``rate = exp(intercept + x . beta)`` fitted by maximum likelihood.

    The log-likelihood of right-censored data is ``sum d_i log(rate_i) - sum t_i rate_i`` and is concave in
    the parameters; the covariates are centred, so ``intercept`` is the log rate of the average sequence.
    """
    n_samples, n_features = features.shape
    feature_mean = features.mean(axis=0)
    design = np.column_stack([np.ones(n_samples), features - feature_mean])
    penalty = np.append(0.0, np.full(n_features, penalizer))
    died = occurred.astype(float)
    total_time = float(times.sum())
    if total_time <= 0.0:
        raise ValueError("The exponential model needs a positive total follow-up time")
    start = np.zeros(n_features + 1)
    start[0] = np.log(died.sum() / total_time)  # exact maximum-likelihood rate without covariates

    def evaluate(params: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
        with np.errstate(over="ignore", invalid="ignore"):
            linear = design @ params
            expected = times * np.exp(linear)  # expected number of events per subject at the current fit
            value = float(died @ linear - expected.sum()) - 0.5 * float(penalty @ (params * params))
            gradient = design.T @ (died - expected) - penalty * params
            curvature = (design.T * expected) @ design + np.diag(penalty)
        if not (np.isfinite(value) and np.all(np.isfinite(gradient)) and np.all(np.isfinite(curvature))):
            return -np.inf, np.zeros(n_features + 1), np.zeros((n_features + 1, n_features + 1))
        return value, gradient, curvature

    params, objective, n_iter, converged = _maximise(evaluate, start, max_iter, tol)
    return {
        "coefficients": params[1:],
        "std_errors": _standard_errors(evaluate(params)[2], first=1),
        "intercept": float(params[0]),
        "feature_means": feature_mean,
        "log_likelihood": objective + 0.5 * float(penalty @ (params * params)),  # without the penalty term
        "n_iter": n_iter,
        "converged": converged,
    }


def _fit_kaplan_meier(times: np.ndarray, occurred: np.ndarray) -> Dict[str, Any]:
    """Kaplan-Meier estimate ``S(u) = prod_{u_k <= u} (1 - deaths_k / at_risk_k)`` with right-censoring."""
    order = np.argsort(times, kind="stable")
    sorted_times = times[order]
    died = sorted_times[occurred[order]]
    event_times = np.unique(died)
    at_risk = (times.size - np.searchsorted(sorted_times, event_times, side="left")).astype(float)
    deaths = np.bincount(np.searchsorted(event_times, died), minlength=event_times.size).astype(float)
    return {"baseline_times": event_times, "baseline_survival": np.cumprod(1.0 - deaths / at_risk)}


def _restricted_mean(knots: np.ndarray, survival: np.ndarray, horizon: float) -> np.ndarray:
    """Area under step survival curves up to ``horizon``.

    Row ``i`` of ``survival`` is the curve value on ``[knots[k], knots[k + 1])`` for every ``k`` (the last
    interval ends at ``horizon``).
    """
    return np.asarray(survival @ np.diff(np.append(knots, horizon)), dtype=float)


class _SequenceRepresentation:
    """The representation of sequences that all tasks of a :class:`MultiTaskPredictor` share.

    A sequence is described by its event counts and the mean of its event embeddings (the package's
    sequence representation). Both blocks are standardised with the statistics of the training sequences;
    the embedding block is projected onto at most ``hidden_dim`` of its principal components. The counts
    block is never compressed, so which events occurred how often is always available to the task heads.
    Instances are immutable after construction.
    """

    def __init__(
        self,
        vocabulary: Dict[str, int],
        embeddings: Dict[str, np.ndarray],
        sequences: List[Any],
        hidden_dim: int,
    ):
        self.vocabulary = vocabulary
        self.embeddings = embeddings
        self.embedding_dim = len(next(iter(embeddings.values())))
        counts, embedded = self._raw(sequences)
        self.count_mean = counts.mean(axis=0)
        self.count_scale = _column_scale(counts)
        self.embedding_mean = embedded.mean(axis=0)
        self.embedding_scale = _column_scale(embedded)
        standardised = (embedded - self.embedding_mean) / self.embedding_scale
        _, singular_values, components = np.linalg.svd(standardised, full_matrices=False)
        rank = int(np.sum(singular_values > _RANK_TOLERANCE * singular_values[0])) if singular_values[0] > 0 else 0
        self.projection = components[: min(hidden_dim, rank)].T

    def _raw(self, sequences: List[Any]) -> Tuple[np.ndarray, np.ndarray]:
        """Event counts and mean event embeddings.

        Events outside the training vocabulary are dropped before anything is computed, so they cannot dilute
        the mean embedding; a vocabulary event without an embedding counts as a zero vector.
        """
        known = [
            [token for token in _sequence_to_tokens(sequence) if token in self.vocabulary] for sequence in sequences
        ]
        counts = _count_matrix(known, self.vocabulary)
        if not known:
            return counts, np.zeros((0, self.embedding_dim))
        embedded = np.asarray(sequence_embeddings(known, self.embeddings, method="mean"), dtype=float)
        return counts, embedded

    def transform(self, sequences: List[Any]) -> np.ndarray:
        """Representation of ``sequences``: standardised counts, then the compressed standardised embeddings."""
        counts, embedded = self._raw(sequences)
        standardised_counts = (counts - self.count_mean) / self.count_scale
        standardised_embedding = (embedded - self.embedding_mean) / self.embedding_scale
        return np.hstack([standardised_counts, standardised_embedding @ self.projection])


class MultiTaskPredictor:
    """Multi-task predictor for life events with multiple outcome types.

    One model predicts several outcomes (classification and regression) from the same life event
    sequences. All tasks share one representation of a sequence; every task has its own head on it.

    **Shared representation.** The event vocabulary is taken from the training sequences and word2vec
    event embeddings of dimension ``embedding_dim`` are learned from them (or supplied to :meth:`fit`).
    A sequence is represented by its standardised event counts together with the mean of its event
    embeddings (the package's sequence representation), the latter compressed to at most ``hidden_dim``
    principal components. The counts are always kept, so the representation never loses which events
    occurred how often; the embeddings add a smooth description of the sequence's composition and
    ``hidden_dim`` bounds how many of its directions are used.

    **Task heads.** Regression tasks use ridge regression (their predictions are values on the scale of the
    training outcomes). Classification tasks are binary and use L2-regularised logistic regression; the
    larger of the two labels is the positive class and the predictions are its probabilities. Both heads
    are fitted by exact optimisation (closed form, respectively Newton-Raphson), so the model is
    deterministic given its embeddings. ``random_state`` seeds the embedding training; without it the
    embeddings, and with them the predictions, differ slightly between fits.
    """

    def __init__(
        self,
        task_types: Dict[str, str],
        embedding_dim: int = 100,
        hidden_dim: int = 64,
        random_state: Optional[int] = None,
        alpha: float = 1.0,
        **kwargs: Any,
    ):
        """Initialize multi-task predictor.

        Args:
            task_types: Dictionary mapping task names to types ('classification' or 'regression')
            embedding_dim: Dimension of event embeddings
            hidden_dim: Maximum number of principal components of the sequence embeddings kept in the
                shared representation
            random_state: Random seed for the embedding training (reproducibility)
            alpha: Strength of the L2 penalty of every task head, on standardised features; must be positive
            **kwargs: Not supported; any option given is ignored with a warning

        Raises:
            ValueError: If there are no tasks, a task type is unknown or a size/penalty is out of range.
        """
        if not task_types:
            raise ValueError("task_types must not be empty")
        invalid = {task: kind for task, kind in task_types.items() if kind not in _TASK_TYPES}
        if invalid:
            raise ValueError(f"Unknown task types {invalid}; every task type must be one of {_TASK_TYPES}")
        if embedding_dim < 1 or hidden_dim < 1:
            raise ValueError("embedding_dim and hidden_dim must be positive")
        if not alpha > 0:
            raise ValueError("alpha must be positive")
        if kwargs:
            logger.warning(f"MultiTaskPredictor ignores unsupported options: {sorted(kwargs)}")

        self.task_types = task_types
        self.embedding_dim = embedding_dim
        self.hidden_dim = hidden_dim
        self.random_state = random_state
        self.alpha = alpha

        # Set by fit
        self.is_fitted = False
        self.num_events: int | None = None
        self.event_vocab: Dict[str, int] | None = None
        self.event_embeddings: Dict[str, np.ndarray] = {}
        self.models: Dict[str, Dict[str, Any]] = {}
        self._representation: _SequenceRepresentation | None = None

    def fit(
        self,
        sequences: List[EventSequence],
        outcomes: Dict[str, np.ndarray],
        event_embeddings: Optional[Dict[str, np.ndarray]] = None,
    ) -> Self:
        """Fit multi-task predictor to sequences and outcomes.

        Args:
            sequences: List of EventSequence objects (or lists of event tokens)
            outcomes: Dictionary mapping every task name to its outcome array, one value per sequence.
                Outcomes must be finite; classification outcomes must have exactly two distinct labels.
            event_embeddings: Optional pre-computed event embeddings (token -> vector). They replace the
                embeddings that are otherwise learned from ``sequences``; events without a vector count as zeros.

        Returns:
            Self for method chaining

        Raises:
            ValueError: If the inputs are empty or inconsistent, a task has no outcomes (or an outcome has
                no task), outcomes are not finite, or a classification task is not binary. A failed fit
                leaves the model as it was before.
        """
        sequence_list = list(sequences)
        if not sequence_list:
            raise ValueError("sequences cannot be empty")
        missing = [task for task in self.task_types if task not in outcomes]
        if missing:
            raise ValueError(f"outcomes are missing for the tasks {missing}")
        unknown = [task for task in outcomes if task not in self.task_types]
        if unknown:
            raise ValueError(f"outcomes contain tasks that are not in task_types: {unknown}")
        targets = {task: self._prepare_target(task, outcomes[task], len(sequence_list)) for task in self.task_types}

        vocabulary = _event_vocabulary(sequence_list)
        if not vocabulary:
            raise ValueError("The training sequences contain no events")
        if event_embeddings is None:
            embeddings = learn_event_embeddings(
                sequence_list, embedding_dim=self.embedding_dim, random_state=self.random_state
            )
        else:
            embeddings = {token: np.asarray(vector, dtype=float) for token, vector in event_embeddings.items()}
            if len({vector.shape for vector in embeddings.values()}) != 1 or next(iter(embeddings.values())).ndim != 1:
                raise ValueError("event_embeddings must map events to one-dimensional vectors of the same length")

        representation = _SequenceRepresentation(vocabulary, embeddings, sequence_list, self.hidden_dim)
        features = representation.transform(sequence_list)
        models: Dict[str, Dict[str, Any]] = {}
        for task, (values, classes) in targets.items():
            if classes is None:
                weights, bias = _fit_ridge_head(features, values, self.alpha)
                models[task] = {"type": "ridge", "weights": weights, "bias": bias}
            else:
                weights, bias = _fit_logistic_head(features, values, self.alpha)
                models[task] = {"type": "logistic", "weights": weights, "bias": bias, "classes": classes}

        self._representation = representation
        self.event_vocab = vocabulary
        self.event_embeddings = embeddings
        self.num_events = len(vocabulary)
        self.models = models
        self.is_fitted = True
        logger.info(
            f"Fitted multi-task predictor: {len(sequence_list)} sequences, {len(vocabulary)} events, "
            f"{len(models)} tasks, shared representation of width {features.shape[1]}"
        )
        return self

    def predict(
        self, sequences: List[EventSequence], task_name: Optional[str] = None
    ) -> Union[Dict[str, List[float]], List[float]]:
        """Predict outcomes for multiple tasks.

        Args:
            sequences: List of EventSequence objects (or lists of event tokens)
            task_name: If specified, return predictions for only this task

        Returns:
            Dictionary mapping task names to prediction lists, or list of predictions if task_name specified.
            Regression tasks give predicted values; classification tasks give the probability of the
            positive (larger) label.

        Raises:
            ValueError: If the model is not fitted or ``task_name`` is not a task of this model.
        """
        if not self.is_fitted or self._representation is None:
            raise ValueError("MultiTaskPredictor must be fitted before predict()")
        if task_name is not None and task_name not in self.task_types:
            raise ValueError(f"Unknown task {task_name!r}; the tasks are {list(self.task_types)}")

        tasks = [task_name] if task_name is not None else list(self.task_types)
        features = self._representation.transform(list(sequences))
        predictions = {task: self._task_output(task, features) for task in tasks}
        if task_name is not None:
            return predictions[task_name]
        return predictions

    def shared_representation(self, sequences: List[EventSequence]) -> np.ndarray:
        """The representation all task heads work on: one row per sequence.

        The columns are the standardised event counts followed by the (at most ``hidden_dim``) principal
        components of the standardised mean event embedding.

        Raises:
            ValueError: If the model is not fitted.
        """
        if not self.is_fitted or self._representation is None:
            raise ValueError("MultiTaskPredictor must be fitted before shared_representation()")
        return self._representation.transform(list(sequences))

    def _prepare_target(self, task: str, values: Any, n_sequences: int) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """Validate one task's outcomes.

        Returns:
            ``(targets, classes)``: float values for a regression task (``classes`` is None), or the 0/1
            indicator of the larger label and the sorted pair of labels for a classification task.
        """
        array = np.asarray(values)
        if array.ndim != 1 or array.shape[0] != n_sequences:
            raise ValueError(f"Outcomes of task {task!r} must be a vector with one value per sequence ({n_sequences})")
        if array.dtype.kind in "fc" and not np.all(np.isfinite(array)):
            raise ValueError(f"Outcomes of task {task!r} must be finite")
        if self.task_types[task] == "regression":
            try:
                numeric = array.astype(float)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Regression outcomes of task {task!r} must be numeric") from exc
            return numeric, None
        classes = np.unique(array)
        if classes.size != 2:
            raise ValueError(
                f"Classification task {task!r} needs exactly two distinct labels, got {classes.size}: "
                "the model is binary"
            )
        return (array == classes[1]).astype(float), classes

    def _task_output(self, task: str, features: np.ndarray) -> List[float]:
        """Apply the head of ``task`` to shared representations."""
        head = self.models[task]
        linear = features @ head["weights"] + head["bias"]
        output = expit(linear) if head["type"] == "logistic" else linear
        return [float(value) for value in output]


class SurvivalPredictor:
    """Survival analysis predictor for life events.

    Predicts survival times and survival functions from life event sequences. The covariates of a
    sequence are the number of times each event type (``"domain:event_type"``) occurs in it. Three methods
    are available; all are NumPy implementations that handle right-censored data:

    - ``"cox"``: Cox proportional hazards model. The coefficients maximise the Breslow partial
      likelihood (ties share their risk set) by Newton-Raphson with step halving. The baseline hazard is
      Breslow's estimate, so survival functions are ``S(t | x) = exp(-H0(t) exp(x . beta))``.
    - ``"exponential"``: exponential proportional hazards model ``hazard = exp(intercept + x . beta)``
      fitted by maximum likelihood (Newton-Raphson). The hazard is constant over time.
    - ``"baseline"`` (alias ``"kaplan_meier"``): the Kaplan-Meier estimate, the same for every sequence.

    A model that is not fitted refuses to predict. Unknown methods are rejected.

    Predicted survival times are restricted mean survival times: the area under the sequence's survival
    function up to the longest follow-up time in the training data. That horizon keeps the prediction
    finite and comparable across methods; Kaplan-Meier and Cox estimates carry no information beyond it.
    """

    def __init__(
        self,
        method: str = "cox",
        embedding_dim: int = 50,
        random_state: Optional[int] = None,
        penalizer: float = 0.1,
        max_iter: int = 50,
        tol: float = 1e-9,
        **kwargs: Any,
    ):
        """Initialize survival predictor.

        Args:
            method: Survival analysis method: 'cox', 'exponential' or 'baseline' ('kaplan_meier')
            embedding_dim: Retained for backward compatibility; the survival models use event-count
                covariates, so this has no effect
            random_state: Retained for backward compatibility; the fits are deterministic
            penalizer: Strength of the ridge penalty on the Cox / exponential coefficients (not on the
                exponential intercept). It keeps the fit finite for rare or perfectly predictive events and
                for collinear counts; use 0 for the unpenalised maximum-likelihood estimate.
            max_iter: Maximum number of Newton-Raphson iterations
            tol: Newton-Raphson stops when a full step improves the log-likelihood by less than this
                (relative) amount
            **kwargs: Not supported; any option given is ignored with a warning

        Raises:
            ValueError: If ``method`` is unknown or ``penalizer``/``max_iter`` are out of range.
        """
        if method not in _SURVIVAL_METHODS:
            raise ValueError(f"Unknown survival method {method!r}; supported methods are {sorted(_SURVIVAL_METHODS)}")
        if not penalizer >= 0:
            raise ValueError("penalizer must not be negative")
        if max_iter < 1:
            raise ValueError("max_iter must be at least 1")
        if kwargs:
            logger.warning(f"SurvivalPredictor ignores unsupported options: {sorted(kwargs)}")

        self.method = method
        self.embedding_dim = embedding_dim
        self.random_state = random_state
        self.penalizer = penalizer
        self.max_iter = max_iter
        self.tol = tol

        # Set by fit
        self.is_fitted = False
        self.event_vocab: Dict[str, int] | None = None
        self.num_events: int | None = None
        self.model: Dict[str, Any] | None = None

    def fit(self, sequences: List[EventSequence], event_times: np.ndarray, event_occurred: np.ndarray) -> Self:
        """Fit survival model to sequences and survival data.

        Args:
            sequences: List of EventSequence objects (or lists of event tokens)
            event_times: Array of event/censoring times (finite, not negative)
            event_occurred: Boolean (or 0/1) array indicating if the event occurred (True) or the
                observation was censored (False)

        Returns:
            Self for method chaining

        Raises:
            ValueError: If the inputs are empty, differ in length or hold invalid values, or no event
                occurred. A Newton-Raphson fit that does not converge emits a ``RuntimeWarning``.
        """
        times, occurred = _validate_survival_data(len(sequences), event_times, event_occurred)
        kind = _SURVIVAL_METHODS[self.method]
        vocabulary = _event_vocabulary(list(sequences))
        if kind != "baseline" and not vocabulary:
            raise ValueError("The training sequences contain no events, so there are no covariates to fit")
        features = _count_matrix(list(sequences), vocabulary)

        if kind == "cox":
            model = _fit_cox(features, times, occurred, self.penalizer, self.max_iter, self.tol)
        elif kind == "exponential":
            model = _fit_exponential(features, times, occurred, self.penalizer, self.max_iter, self.tol)
        else:
            model = _fit_kaplan_meier(times, occurred)
        model["type"] = kind
        model["horizon"] = float(times.max())
        if not model.get("converged", True):
            warnings.warn(
                f"The {kind} model did not converge in {self.max_iter} iterations; consider a larger penalizer "
                "(perfectly predictive events make the unpenalised likelihood unbounded) or more iterations",
                RuntimeWarning,
                stacklevel=2,
            )

        self.event_vocab = vocabulary
        self.num_events = len(vocabulary)
        self.model = model
        self.is_fitted = True
        logger.info(f"Fitted {kind} survival model: {len(times)} sequences, {int(occurred.sum())} events")
        return self

    @property
    def coefficients(self) -> Dict[str, float]:
        """Log hazard ratio of one more occurrence of each event type (empty for the baseline method).

        Raises:
            ValueError: If the model is not fitted.
        """
        return self._per_event("coefficients")

    @property
    def standard_errors(self) -> Dict[str, float]:
        """Standard error of each coefficient (empty for the baseline method).

        They come from the inverse of the information matrix at the fitted coefficients, so with a
        non-zero ``penalizer`` they describe the penalised fit.

        Raises:
            ValueError: If the model is not fitted.
        """
        return self._per_event("std_errors")

    def _per_event(self, key: str) -> Dict[str, float]:
        """One entry of the fitted model, as ``{event token: value}``."""
        model = self._fitted_model()
        if model["type"] == "baseline" or self.event_vocab is None:
            return {}
        return {token: float(model[key][index]) for token, index in self.event_vocab.items()}

    def predict(self, sequences: List[EventSequence]) -> np.ndarray:
        """Predict survival times for sequences.

        Args:
            sequences: List of EventSequence objects (or lists of event tokens)

        Returns:
            Array of predicted survival times: the restricted mean survival time of every sequence, i.e.
            the area under its survival function up to the longest follow-up time seen in training.

        Raises:
            ValueError: If the model is not fitted.
        """
        model = self._fitted_model()
        horizon = model["horizon"]
        if model["type"] == "baseline":
            knots = np.concatenate(([0.0], model["baseline_times"]))
            curve = np.concatenate(([1.0], model["baseline_survival"]))
            return np.full(len(sequences), _restricted_mean(knots, curve[None, :], horizon)[0])
        if model["type"] == "exponential":
            rate = self._rate(sequences)
            return np.asarray(-np.expm1(-rate * horizon) / rate, dtype=float)
        knots = np.concatenate(([0.0], model["baseline_times"]))
        cumulative = np.concatenate(([0.0], model["baseline_cumulative_hazard"]))
        survival = np.exp(-np.outer(self._relative_hazard(sequences), cumulative))
        return _restricted_mean(knots, survival, horizon)

    def predict_survival_function(self, sequences: List[EventSequence], times: np.ndarray) -> np.ndarray:
        """Predict survival function at given times.

        Args:
            sequences: List of EventSequence objects (or lists of event tokens)
            times: Array of time points (finite, not negative)

        Returns:
            Array of survival probabilities (sequences x times). For the Cox model this is
            ``exp(-H0(t) exp(x . beta))`` with Breslow's baseline hazard ``H0``, for the exponential model
            ``exp(-rate t)`` and for the baseline method the Kaplan-Meier curve. The estimates are step
            functions that stay constant after the last observed event time (Cox, Kaplan-Meier).

        Raises:
            ValueError: If the model is not fitted or ``times`` holds invalid values.
        """
        model = self._fitted_model()
        grid = np.atleast_1d(np.asarray(times, dtype=float))
        if grid.ndim != 1 or not np.all(np.isfinite(grid)) or np.any(grid < 0):
            raise ValueError("times must be a one-dimensional array of finite, non-negative time points")

        if model["type"] == "exponential":
            return np.asarray(np.exp(-np.outer(self._rate(sequences), grid)), dtype=float)
        steps = np.searchsorted(model["baseline_times"], grid, side="right")
        if model["type"] == "baseline":
            curve = np.concatenate(([1.0], model["baseline_survival"]))[steps]
            return np.tile(curve, (len(sequences), 1))
        cumulative = np.concatenate(([0.0], model["baseline_cumulative_hazard"]))[steps]
        return np.asarray(np.exp(-np.outer(self._relative_hazard(sequences), cumulative)), dtype=float)

    def _fitted_model(self) -> Dict[str, Any]:
        """The fitted model dictionary, or ``ValueError`` when ``fit`` has not succeeded."""
        if not self.is_fitted or self.model is None:
            raise ValueError("SurvivalPredictor must be fitted before it can predict")
        return self.model

    def _sequences_to_features(self, sequences: List[EventSequence]) -> np.ndarray:
        """Convert sequences to event-count feature vectors (events outside the training vocabulary are ignored)."""
        return _count_matrix(list(sequences), self.event_vocab or {})

    def _log_hazard(self, sequences: List[EventSequence]) -> np.ndarray:
        """Linear predictor ``(x - mean) . beta`` of every sequence, clipped so it can be exponentiated."""
        model = self._fitted_model()
        centred = self._sequences_to_features(sequences) - model["feature_means"]
        return np.asarray(np.clip(centred @ model["coefficients"], -_MAX_LOG_HAZARD, _MAX_LOG_HAZARD), dtype=float)

    def _relative_hazard(self, sequences: List[EventSequence]) -> np.ndarray:
        """Hazard of every sequence relative to the average training sequence (Cox model)."""
        return np.asarray(np.exp(self._log_hazard(sequences)), dtype=float)

    def _rate(self, sequences: List[EventSequence]) -> np.ndarray:
        """Constant hazard rate of every sequence (exponential model)."""
        return np.asarray(np.exp(self._log_hazard(sequences) + self._fitted_model()["intercept"]), dtype=float)


def attention_weights(
    model: Any, sequences: List[EventSequence], embeddings: Optional[Dict[str, np.ndarray]] = None
) -> np.ndarray:
    """Compute attention weights for event sequences using a trained model.

    Args:
        model: Trained sequence model with attention mechanism
        sequences: List of EventSequence objects

    Returns:
        Attention matrix with one row per sequence
    """
    token_sequences = [_sequence_to_tokens(seq) for seq in sequences]
    max_len = max((len(seq) for seq in token_sequences), default=0)
    if max_len == 0:
        return np.zeros((len(sequences), 0))

    if hasattr(model, "get_attention_weights"):
        rows = []
        for seq in sequences:
            weights = model.get_attention_weights(seq)
            row = (
                np.array(weights, dtype=float).ravel()
                if weights is not None
                else np.ones(len(_sequence_to_tokens(seq)))
            )
            if row.sum() > 0:
                row = row / row.sum()
            rows.append(np.pad(row, (0, max_len - len(row))))
        return np.vstack(rows)

    rows = []
    for tokens in token_sequences:
        if not tokens:
            rows.append(np.zeros(max_len))
            continue
        if embeddings:
            magnitudes = np.array([float(np.linalg.norm(embeddings.get(token, np.zeros(1)))) for token in tokens])
            weights = magnitudes if magnitudes.sum() > 0 else np.ones(len(tokens))
        else:
            weights = np.ones(len(tokens))
        weights = weights / weights.sum()
        rows.append(np.pad(weights, (0, max_len - len(weights))))

    return np.vstack(rows)
