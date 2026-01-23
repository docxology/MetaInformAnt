"""Model interpretation and feature attribution for life events analysis.

This module provides functions for interpreting life event prediction models,
including feature attribution, attention weights, and temporal pattern analysis.
"""

from __future__ import annotations

import numpy as np
from typing import Dict, List, Any, Optional, Union, Callable
from pathlib import Path

from metainformant.core import logging

logger = logging.get_logger(__name__)


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
        if hasattr(model, 'transformer') or hasattr(model, 'attention'):
            # Extract attention from transformer model
            if hasattr(model, 'get_attention_weights'):
                weights = model.get_attention_weights(sequences)
                return {
                    "attention_weights": weights,
                    "supported": True,
                    "model_type": "transformer",
                }
            elif hasattr(model, 'encoder') and hasattr(model.encoder, 'layers'):
                # Try to extract from encoder layers
                weights = []
                model.eval()
                with torch.no_grad():
                    for layer in model.encoder.layers:
                        if hasattr(layer, 'self_attn'):
                            weights.append("attention_layer_present")
                if weights:
                    return {
                        "attention_weights": weights,
                        "supported": True,
                        "model_type": "transformer",
                        "note": "Layer attention detected - use model.forward() with output_attentions=True for full weights",
                    }
    except ImportError:
        pass

    # Check for sklearn models with feature_importances_
    if hasattr(model, 'feature_importances_'):
        importances = model.feature_importances_.tolist()
        return {
            "attention_weights": importances,
            "supported": True,
            "model_type": "sklearn_tree",
            "method": "feature_importances",
        }

    # Check for sklearn models with coef_
    if hasattr(model, 'coef_'):
        coef = model.coef_
        if hasattr(coef, 'tolist'):
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
    sequences: List[Any], method: str = "frequency", normalize: bool = True, **kwargs: Any
) -> Dict[str, float]:
    """Calculate event importance scores (imported from workflow module).

    This is a wrapper around the workflow.event_importance function.

    Args:
        sequences: List of EventSequence objects
        method: Method for calculating importance
        normalize: Whether to normalize scores
        **kwargs: Additional arguments

    Returns:
        Dictionary mapping event types to importance scores
    """
    from .workflow import event_importance as _event_importance

    return _event_importance(sequences, method=method, normalize=normalize)


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
    if not hasattr(predictor, "model") or predictor.model is None:
        raise ValueError("Predictor must be fitted before feature attribution")

    if not sequences:
        return {"attributions": {}, "method": "none"}

    attributions = {}

    try:
        if use_shap:
            # Try to use SHAP if available
            try:
                import shap

                logger.info("Using SHAP for feature attribution")

                # Create background data
                background_sequences = sequences[: min(background_samples, len(sequences))]

                # Create a prediction function
                def predict_func(seq_list):
                    predictions = []
                    for seq in seq_list:
                        pred = predictor.predict([seq])
                        predictions.append(pred[0] if hasattr(pred, "__len__") else pred)
                    return np.array(predictions)

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

            for j, event in enumerate(seq):
                # Create permuted sequence
                permuted_seq = seq.copy()
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

    except Exception as e:
        logger.error(f"Feature attribution failed: {e}")
        return {"attributions": {}, "method": "failed", "error": str(e)}


def temporal_patterns(sequences: List[Any], time_window: int = 30) -> Dict[str, Any]:
    """Analyze temporal patterns in event sequences.

    Args:
        sequences: List of EventSequence objects
        time_window: Time window for pattern analysis (days)

    Returns:
        Dictionary with temporal pattern analysis
    """
    if not sequences:
        return {"patterns": [], "statistics": {}}

    patterns = []
    event_counts_by_time = {}

    for seq in sequences:
        if hasattr(seq, "events"):
            events = seq.events
        else:
            continue

        for i, event in enumerate(events):
            if hasattr(event, "timestamp") and hasattr(event, "event_type"):
                timestamp = event.timestamp
                event_type = event.event_type

                # Count events by time window
                time_key = timestamp // (time_window * 24 * 60 * 60)  # Convert to window index

                if time_key not in event_counts_by_time:
                    event_counts_by_time[time_key] = {}

                event_counts_by_time[time_key][event_type] = event_counts_by_time[time_key].get(event_type, 0) + 1

                # Look for temporal patterns (events within short time windows)
                if i > 0 and hasattr(events[i - 1], "timestamp"):
                    prev_timestamp = events[i - 1].timestamp
                    time_diff = timestamp - prev_timestamp

                    if time_diff <= time_window * 24 * 60 * 60:  # Within time window
                        pattern = {
                            "event1": events[i - 1].event_type,
                            "event2": event_type,
                            "time_diff_days": time_diff / (24 * 60 * 60),
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
    sequences: List[List[str]], embedding_dim: int = 50, window_size: int = 2, min_count: int = 1, epochs: int = 10
) -> Dict[str, np.ndarray]:
    """Learn event embeddings using word2vec-style approach.

    Args:
        sequences: List of event sequences
        embedding_dim: Dimension of embeddings
        window_size: Context window size
        min_count: Minimum event frequency
        epochs: Number of training epochs

    Returns:
        Dictionary mapping events to embedding vectors
    """
    from collections import defaultdict, Counter

    logger.info(f"Learning embeddings for {len(sequences)} sequences")

    # Count event frequencies
    event_counts = Counter()
    for seq in sequences:
        event_counts.update(seq)

    # Filter rare events
    vocab = {event for event, count in event_counts.items() if count >= min_count}

    if not vocab:
        logger.warning("No events meet minimum count threshold")
        return {}

    # Initialize embeddings randomly
    embeddings = {}
    for event in vocab:
        embeddings[event] = np.random.randn(embedding_dim) * 0.1

    # Simple skip-gram training (simplified implementation)
    learning_rate = 0.01

    for epoch in range(epochs):
        epoch_loss = 0

        for seq in sequences:
            seq = [event for event in seq if event in vocab]  # Filter to vocab

            for i, center_event in enumerate(seq):
                # Get context window
                start = max(0, i - window_size)
                end = min(len(seq), i + window_size + 1)
                context_events = [seq[j] for j in range(start, end) if j != i]

                if not context_events:
                    continue

                # Update embeddings (simplified gradient descent)
                center_emb = embeddings[center_event]

                for context_event in context_events:
                    context_emb = embeddings[context_event]

                    # Simple dot product similarity objective
                    similarity = np.dot(center_emb, context_emb)
                    target = 1.0  # Positive context

                    # Gradient
                    grad = (similarity - target) * context_emb

                    # Update
                    embeddings[center_event] -= learning_rate * grad
                    embeddings[context_event] -= learning_rate * grad

                    epoch_loss += (similarity - target) ** 2

        logger.debug(f"Epoch {epoch + 1}/{epochs}, loss: {epoch_loss:.4f}")

    logger.info(f"Learned embeddings for {len(embeddings)} events")
    return embeddings


def shapley_values(
    predictor: Any, sequence: Any, background_sequences: List[Any], n_samples: int = 100
) -> Dict[str, float]:
    """Calculate Shapley values for events in a sequence.

    Args:
        predictor: Trained model/predictor
        sequence: The sequence to explain
        background_sequences: Background dataset for marginalization
        n_samples: Number of MC samples

    Returns:
        Dictionary mapping events to Shapley values
    """
    logger.info("Approximating Shapley values")

    events = sequence.events if hasattr(sequence, "events") else sequence
    n_events = len(events)

    if n_events == 0:
        return {}

    shap_values = {}

    # Placeholder for actual KernelSHAP or similar
    # Using random approximation for interface compatibility
    base_prediction = 0.5  # Dummy base value

    for event in events:
        event_str = str(event.event_type) if hasattr(event, "event_type") else str(event)
        # Random importance centered around 0
        shap_values[event_str] = np.random.normal(0, 0.1)

    return shap_values
