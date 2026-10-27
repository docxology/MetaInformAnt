"""Model interpretation and feature attribution for life events analysis.

This module provides functions for interpreting life event prediction models,
including feature attribution, attention weights, and temporal pattern analysis.
"""

from __future__ import annotations

import math
from types import SimpleNamespace
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core.utils import logging
from metainformant.life_events.core.utils import _sequence_to_tokens

logger = logging.get_logger(__name__)

# Exact enumeration evaluates the predictor on 2**n coalitions; beyond 20 events that is over a million.
_MAX_EXACT_EVENTS_LIMIT = 20


def attention_weights(
    model: Any, sequences: List[List[str]], embeddings: Optional[Dict[str, np.ndarray]] = None
) -> Dict[str, Any]:
    """Extract attention weights from a model.

    Supports models that expose attention weights (e.g., transformer models).
    For sklearn or other non-attention models, computes proxy attention via
    feature importance gradients if possible.

    Args:
        model: Trained prediction model (sklearn, torch, or custom)
        sequences: List of event sequences
        embeddings: Optional event embeddings

    Returns:
        Dictionary with attention weight information
    """
    # Check for PyTorch transformer model
    try:
        import torch

        if hasattr(model, "transformer") or hasattr(model, "attention"):
            # Extract attention from transformer model
            if hasattr(model, "get_attention_weights"):
                weights = model.get_attention_weights(sequences)
                return {
                    "attention_weights": weights,
                    "supported": True,
                    "model_type": "transformer",
                }
            elif hasattr(model, "encoder") and hasattr(model.encoder, "layers"):
                # Try to extract from encoder layers
                weights = []
                model.eval()
                with torch.no_grad():
                    for layer in model.encoder.layers:
                        if hasattr(layer, "self_attn"):
                            weights.append("attention_layer_present")
                if weights:
                    return {
                        "attention_weights": weights,
                        "supported": True,
                        "model_type": "transformer",
                        "note": (
                            "Layer attention detected - use model.forward() with "
                            "output_attentions=True for full weights"
                        ),
                    }
    except ImportError:
        pass

    # Check for sklearn models with feature_importances_
    if hasattr(model, "feature_importances_"):
        importances = model.feature_importances_.tolist()
        return {
            "attention_weights": importances,
            "supported": True,
            "model_type": "sklearn_tree",
            "method": "feature_importances",
        }

    # Check for sklearn models with coef_
    if hasattr(model, "coef_"):
        coef = model.coef_
        if hasattr(coef, "tolist"):
            coef = coef.tolist()
        return {
            "attention_weights": coef,
            "supported": True,
            "model_type": "sklearn_linear",
            "method": "coefficients",
        }

    # No attention weights available
    logger.info("Model does not support attention weight extraction")
    return {
        "attention_weights": [],
        "supported": False,
        "model_type": str(type(model).__name__),
        "message": "Attention weights not available for this model type. "
        "Use transformer models or tree/linear models for feature attribution.",
    }


def event_importance(
    predictor_or_sequences: Any,
    sequences: Optional[List[Any]] = None,
    embeddings: Optional[Dict[str, np.ndarray]] = None,
    method: str = "frequency",
    normalize: bool = True,
    **kwargs: Any,
) -> Dict[str, float]:
    """Calculate event importance scores from frequencies or fitted predictors.

    With ``method="frequency"`` (or without a predictor) an event's score is how often it occurs.
    With ``method="permutation"`` it is the total absolute change in the predictor's output when that
    event is removed from a sequence, so events whose removal never changes a prediction score
    exactly zero and remain tied.

    Args:
        predictor_or_sequences: A fitted predictor, or the sequences when no predictor is used
        sequences: The sequences to score when a predictor is given
        embeddings: Accepted for backward compatibility; importance depends on the predictor's
            outputs only, so embedding norms never break ties
        method: ``"frequency"`` or ``"permutation"``
        normalize: Scale scores so the largest is 1
        **kwargs: Ignored

    Returns:
        Dictionary mapping event tokens to importance scores
    """
    if sequences is None:
        predictor = None
        sequences = predictor_or_sequences
    else:
        predictor = predictor_or_sequences

    if not sequences:
        raise ValueError("sequences cannot be empty")
    if method not in {"frequency", "permutation"}:
        raise ValueError("Invalid method")
    if predictor is not None and not getattr(predictor, "is_fitted", False):
        raise ValueError("Predictor must be fitted before computing event importance")

    scores: dict[str, float] = {}
    if method == "frequency" or predictor is None:
        for seq in sequences:
            for event in _sequence_to_tokens(seq):
                scores[event] = scores.get(event, 0.0) + 1.0
    else:
        baseline = predictor.predict(sequences)
        for i, seq in enumerate(sequences):
            tokens = _sequence_to_tokens(seq)
            for j, event in enumerate(tokens):
                reduced = tokens[:j] + tokens[j + 1 :]
                if not reduced:
                    delta = 0.0
                else:
                    changed = predictor.predict([reduced])[0]
                    delta = abs(float(baseline[i]) - float(changed))
                scores[event] = scores.get(event, 0.0) + delta

    if normalize and scores:
        max_score = max(scores.values()) or 1.0
        scores = {event: float(score / max_score) for event, score in scores.items()}

    return scores


def feature_attribution(
    predictor: Any,
    sequences: List[List[str]],
    embeddings: Dict[str, np.ndarray],
    use_shap: bool = False,
    background_samples: int = 10,
) -> Dict[str, Any]:
    """Calculate feature attribution for event predictions.

    Args:
        predictor: Trained EventSequencePredictor
        sequences: List of event sequences to analyze
        embeddings: Dictionary mapping events to embeddings
        use_shap: Whether to use SHAP (if available)
        background_samples: Number of background samples for SHAP

    Returns:
        Dictionary with feature attribution results
    """
    if not getattr(predictor, "is_fitted", False):
        raise ValueError("Predictor must be fitted before feature attribution")

    if not sequences:
        raise ValueError("sequences cannot be empty")

    attributions = {}

    try:
        if use_shap:
            # Try to use SHAP if available
            try:
                import shap  # type: ignore[import-not-found]

                logger.info("Using SHAP for feature attribution")

                # Create background data
                background_sequences = sequences[: min(background_samples, len(sequences))]

                # Create a prediction function
                def predict_func(seq_list: list[Any]) -> list[Any]:
                    predictions = []
                    for seq in seq_list:
                        pred = predictor.predict([seq])
                        predictions.append(pred[0] if hasattr(pred, "__len__") else pred)
                    return list(np.asarray(predictions))

                # Create explainer
                explainer = shap.Explainer(predict_func, background_sequences)

                # Calculate attributions for each sequence
                for i, seq in enumerate(sequences):
                    shap_values = explainer([seq])
                    attributions[f"sequence_{i}"] = {
                        "events": seq,
                        "shap_values": shap_values.values.tolist(),
                        "base_value": float(shap_values.base_values[0]),
                    }

                return {"attributions": attributions, "method": "shap", "library": "shap"}

            except ImportError:
                logger.warning("SHAP not available, falling back to permutation importance")
                use_shap = False

        # Fallback: simple permutation importance
        logger.info("Using permutation importance for feature attribution")

        # Get baseline predictions
        baseline_predictions = predictor.predict(sequences)

        for i, seq in enumerate(sequences):
            seq_attributions = {}

            tokens = _sequence_to_tokens(seq)
            for j, event in enumerate(tokens):
                # Create permuted sequence
                permuted_seq = tokens.copy()
                # Remove the event (simple ablation)
                permuted_seq.pop(j)

                if permuted_seq:  # Only if sequence is not empty after removal
                    permuted_pred = predictor.predict([permuted_seq])

                    # Calculate attribution as change in prediction
                    if hasattr(baseline_predictions[i], "__len__"):
                        baseline = baseline_predictions[i][0]
                        permuted = permuted_pred[0][0]
                    else:
                        baseline = baseline_predictions[i]
                        permuted = permuted_pred[0]

                    attribution = baseline - permuted
                else:
                    attribution = 0.0  # Can't remove from single-event sequence

                seq_attributions[event] = float(attribution)

            attributions[f"sequence_{i}"] = {
                "events": seq,
                "attributions": seq_attributions,
                "total_attribution": sum(seq_attributions.values()),
            }

        return {"attributions": attributions, "method": "permutation", "library": None}

    except (ValueError, TypeError, AttributeError) as e:
        logger.error(f"Feature attribution failed due to data error: {e}")
        return {"attributions": {}, "method": "failed", "error": str(e)}
    except (IndexError, KeyError) as e:
        logger.error(f"Feature attribution failed due to sequence access error: {e}")
        return {"attributions": {}, "method": "failed", "error": str(e)}


def temporal_patterns(sequences: List[Any], predictions: Optional[Any] = None, time_window: int = 30) -> Dict[str, Any]:
    """Analyze temporal patterns in event sequences.

    Args:
        sequences: List of EventSequence objects
        time_window: Time window for pattern analysis (days)

    Returns:
        Dictionary with temporal pattern analysis
    """
    if not sequences:
        raise ValueError("sequences cannot be empty")

    if predictions is not None and not isinstance(predictions, (int, float)):
        predictions_arr = np.asarray(predictions)
        if len(predictions_arr) != len(sequences):
            raise ValueError("predictions length must match sequences")
        token_sequences = [_sequence_to_tokens(seq) for seq in sequences]
        max_len = max((len(seq) for seq in token_sequences), default=0)
        position_scores = np.zeros(max_len)
        position_counts = np.zeros(max_len)
        for seq, pred in zip(token_sequences, predictions_arr):
            for i, _event in enumerate(seq):
                position_scores[i] += float(pred)
                position_counts[i] += 1
        position_importance = {
            i: float(position_scores[i] / position_counts[i]) if position_counts[i] else 0.0 for i in range(max_len)
        }
        return {"position_importance": position_importance, "max_sequence_length": max_len}

    patterns = []
    event_counts_by_time: Dict[int, Any] = {}

    for seq in sequences:
        if not hasattr(seq, "events"):
            continue
        events = seq.events

        for i, event in enumerate(events):
            if hasattr(event, "timestamp") and hasattr(event, "event_type"):
                timestamp = event.timestamp
                event_type = event.event_type

                # Count events by time window
                time_key = int(timestamp.timestamp()) // (time_window * 24 * 60 * 60)  # Window index

                if time_key not in event_counts_by_time:
                    event_counts_by_time[time_key] = {}

                event_counts_by_time[time_key][event_type] = event_counts_by_time[time_key].get(event_type, 0) + 1

                # Look for temporal patterns (events within short time windows)
                if i > 0 and hasattr(events[i - 1], "timestamp"):
                    prev_timestamp = events[i - 1].timestamp
                    time_diff = timestamp - prev_timestamp

                    if time_diff.total_seconds() <= time_window * 24 * 60 * 60:  # Within time window
                        pattern = {
                            "event1": events[i - 1].event_type,
                            "event2": event_type,
                            "time_diff_days": time_diff.total_seconds() / (24 * 60 * 60),
                            "sequence_count": 1,
                        }
                        patterns.append(pattern)

    # Aggregate patterns
    pattern_counts = {}
    for pattern in patterns:
        key = (pattern["event1"], pattern["event2"])
        if key not in pattern_counts:
            pattern_counts[key] = {
                "event1": pattern["event1"],
                "event2": pattern["event2"],
                "count": 0,
                "avg_time_diff": 0,
            }
        pattern_counts[key]["count"] += 1
        pattern_counts[key]["avg_time_diff"] += pattern["time_diff_days"]

    # Calculate averages
    for pattern_data in pattern_counts.values():
        if pattern_data["count"] > 0:
            pattern_data["avg_time_diff"] /= pattern_data["count"]

    return {
        "patterns": list(pattern_counts.values()),
        "time_windows": event_counts_by_time,
        "statistics": {
            "total_patterns": len(pattern_counts),
            "total_sequences": len(sequences),
            "time_window_days": time_window,
        },
    }


def learn_event_embeddings(
    sequences: List[List[str]],
    embedding_dim: int = 50,
    window_size: int = 2,
    min_count: int = 1,
    epochs: int = 10,
    random_state: Optional[int] = None,
) -> Dict[str, np.ndarray]:
    """Learn event embeddings using word2vec-style skip-gram training.

    Thin wrapper around :func:`metainformant.life_events.models.embeddings.learn_event_embeddings`
    (skip-gram with negative sampling) that keeps this module's historical defaults and returns an
    empty dictionary, instead of raising, when no event reaches ``min_count``.

    Args:
        sequences: List of event sequences
        embedding_dim: Dimension of embeddings
        window_size: Context window size
        min_count: Minimum event frequency
        epochs: Number of training epochs
        random_state: Optional seed; identical seeds give identical embeddings

    Returns:
        Dictionary mapping events to embedding vectors
    """
    from collections import Counter

    from metainformant.life_events.models.embeddings import learn_event_embeddings as train_embeddings

    logger.info(f"Learning embeddings for {len(sequences)} sequences")

    event_counts: Counter[str] = Counter()
    for seq in sequences:
        event_counts.update(_sequence_to_tokens(seq))

    if not any(count >= min_count for count in event_counts.values()):
        logger.warning("No events meet minimum count threshold")
        return {}

    embeddings = train_embeddings(
        sequences,
        embedding_dim=embedding_dim,
        window_size=window_size,
        min_count=min_count,
        epochs=epochs,
        method="skipgram",
        random_state=random_state,
    )
    logger.info(f"Learned embeddings for {len(embeddings)} events")
    return embeddings


def _prediction_scores(predictor: Any, sequences: List[List[str]]) -> np.ndarray:
    """Score ``sequences`` with ``predictor``: the first component of each prediction, as a float.

    This is the same score the permutation-based importance functions of this module compare.
    """
    predictions = predictor.predict(sequences)
    if len(predictions) != len(sequences):
        raise ValueError(f"Predictor returned {len(predictions)} predictions for {len(sequences)} sequences")
    scores = np.empty(len(sequences), dtype=float)
    for i, prediction in enumerate(predictions):
        try:
            values = np.asarray(prediction, dtype=float).reshape(-1)
        except (TypeError, ValueError) as exc:
            raise ValueError("Predictor output must be numeric to compute event contributions") from exc
        if values.size == 0:
            raise ValueError("Predictor returned an empty prediction")
        scores[i] = values[0]
    return scores


def _exact_shapley_values(predictor: Any, tokens: List[str]) -> np.ndarray:
    """Exact Shapley values of ``tokens`` by enumerating every coalition (2**n evaluations)."""
    n_events = len(tokens)
    n_coalitions = 1 << n_events
    coalitions = [[tokens[j] for j in range(n_events) if (mask >> j) & 1] for mask in range(n_coalitions)]
    values = _prediction_scores(predictor, coalitions)

    masks = np.arange(n_coalitions)
    sizes = np.zeros(n_coalitions, dtype=np.int64)
    for j in range(n_events):
        sizes += (masks >> j) & 1
    # Probability that exactly a given coalition of size k precedes a given player in a random order.
    weights = np.array(
        [math.factorial(k) * math.factorial(n_events - k - 1) / math.factorial(n_events) for k in range(n_events)]
    )

    shapley = np.zeros(n_events)
    for j in range(n_events):
        without = masks[((masks >> j) & 1) == 0]
        shapley[j] = np.sum(weights[sizes[without]] * (values[without | (1 << j)] - values[without]))
    return shapley


def _sampled_shapley_values(predictor: Any, tokens: List[str], n_samples: int, rng: np.random.Generator) -> np.ndarray:
    """Monte Carlo Shapley values from ``n_samples`` random permutations of the events.

    Each permutation adds the events one by one and credits every event with the change in the
    predictor's output when it joins. The credits of one permutation telescope to
    ``v(all) - v(empty)``, so the averaged estimate satisfies the efficiency property exactly.
    """
    n_events = len(tokens)
    baseline = float(_prediction_scores(predictor, [[]])[0])
    full = float(_prediction_scores(predictor, [list(tokens)])[0])

    cache: Dict[bytes, float] = {}
    shapley = np.zeros(n_events)
    for _ in range(n_samples):
        order = rng.permutation(n_events)
        present = np.zeros(n_events, dtype=bool)
        keys: List[bytes] = []
        pending: Dict[bytes, List[str]] = {}
        for position in order[:-1]:
            present[position] = True
            key = present.tobytes()
            keys.append(key)
            if key not in cache and key not in pending:
                pending[key] = [tokens[j] for j in np.flatnonzero(present)]
        if pending:
            cache.update(
                zip(pending, (float(score) for score in _prediction_scores(predictor, list(pending.values()))))
            )

        previous = baseline
        for position, key in zip(order[:-1], keys):
            current = cache[key]
            shapley[position] += current - previous
            previous = current
        shapley[order[-1]] += full - previous
    return shapley / n_samples


def shapley_values(
    predictor: Any,
    sequence: Any,
    background_sequences: Optional[List[Any]] = None,
    n_samples: int = 100,
    *,
    random_state: Optional[int] = None,
    max_exact_events: int = 10,
) -> Dict[str, float]:
    """Calculate Shapley values for the events in a sequence.

    Every event occurrence in ``sequence`` is a player. The value of a coalition is the predictor's
    output (the first component of ``predictor.predict``, the score the other importance functions in
    this module use) on the sequence restricted to the coalition's events. Removing an event drops it
    from the sequence; the remaining events keep their order. The empty coalition is therefore valued
    at the prediction for an empty sequence. The Shapley value of an event is its marginal
    contribution to that value averaged over all orders in which the events can be added.

    The values satisfy the efficiency property: they sum to ``v(all events) - v(empty sequence)``.
    A "null" event that never changes the prediction receives exactly zero.

    Up to ``max_exact_events`` events the values are computed exactly by enumerating all ``2**n``
    coalitions. With more events they are estimated from ``n_samples`` random permutations of the
    events (seeded by ``random_state``); the estimate still sums exactly to ``v(all) - v(empty)``.

    For classification predictors ``predict`` returns class labels, so the contributions are
    changes of the predicted label. Use a regression predictor, or one whose ``predict`` returns a
    score, for graded contributions.

    Args:
        predictor: Fitted model with ``predict(list_of_sequences)``, called on lists of event tokens
        sequence: The sequence to explain: an EventSequence or a list of event tokens
        background_sequences: Not used. Absent events are removed rather than replaced by values
            drawn from a background dataset, so no background distribution enters the value
            function. Kept for interface compatibility.
        n_samples: Number of random permutations when the sequence has more than
            ``max_exact_events`` events
        random_state: Seed for the permutation sampling
        max_exact_events: Largest event count handled by exact enumeration (at most 20)

    Returns:
        Dictionary mapping each event token to its Shapley value. Events that occur several times
        are reported once, with the sum of the values of their occurrences.

    Raises:
        ValueError: If the predictor is not fitted, an option is out of range or the predictor's
            output is not numeric.
    """
    if n_samples < 1:
        raise ValueError("n_samples must be at least 1")
    if not 0 <= max_exact_events <= _MAX_EXACT_EVENTS_LIMIT:
        raise ValueError(f"max_exact_events must be between 0 and {_MAX_EXACT_EVENTS_LIMIT}")
    if getattr(predictor, "is_fitted", True) is False:
        raise ValueError("Predictor must be fitted before computing Shapley values")
    if not callable(getattr(predictor, "predict", None)):
        raise TypeError("predictor must provide a predict() method")

    items = list(sequence.events) if hasattr(sequence, "events") else list(sequence)
    tokens = _sequence_to_tokens(SimpleNamespace(events=items) if items and hasattr(items[0], "event_type") else items)
    n_events = len(tokens)
    if n_events == 0:
        return {}

    if background_sequences:
        logger.debug("shapley_values does not use background_sequences")

    if n_events <= max_exact_events:
        logger.info(f"Computing exact Shapley values for {n_events} events")
        contributions = _exact_shapley_values(predictor, tokens)
    else:
        logger.info(f"Estimating Shapley values for {n_events} events from {n_samples} permutations")
        contributions = _sampled_shapley_values(predictor, tokens, n_samples, np.random.default_rng(random_state))

    shapley: Dict[str, float] = {}
    for token, value in zip(tokens, contributions):
        shapley[token] = shapley.get(token, 0.0) + float(value)
    return shapley
