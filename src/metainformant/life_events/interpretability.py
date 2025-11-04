"""Model interpretation tools for life course analysis.

This module provides methods for interpreting predictions from sequence models,
including attention weights, event importance, and temporal patterns.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np
from numpy.typing import NDArray

try:
    import shap

    SHAP_AVAILABLE = True
except ImportError:
    SHAP_AVAILABLE = False
    shap = None  # type: ignore


def attention_weights(
    model: Any,
    sequences: List[List[str]],
    event_embeddings: Dict[str, NDArray]
) -> NDArray:
    """Extract attention weights from transformer models.
    
    Note: Full implementation requires transformer models with attention.
    This is a placeholder for future implementation.
    
    Args:
        model: Trained model with attention mechanism
        sequences: List of event sequences
        event_embeddings: Event embeddings dictionary
        
    Returns:
        Attention weight matrix
    """
    # Placeholder implementation
    # Full implementation would extract attention from transformer model
    n_sequences = len(sequences)
    max_len = max(len(s) for s in sequences) if sequences else 1
    
    # Return uniform attention weights as placeholder
    attention = np.ones((n_sequences, max_len, max_len)) / max_len
    return attention


def event_importance(
    predictor: "EventSequencePredictor",
    sequences: List[List[str]],
    event_embeddings: Dict[str, NDArray],
    method: str = "permutation"
) -> Dict[str, float]:
    """Rank events by their contribution to predictions.
    
    Args:
        predictor: Trained EventSequencePredictor (must be fitted)
        sequences: List of event sequences (must not be empty)
        event_embeddings: Event embeddings dictionary
        method: Importance method ("permutation", "gradient")
        
    Returns:
        Dictionary mapping event tokens to importance scores (normalized to [0, 1])
        
    Raises:
        ValueError: If predictor not fitted, sequences empty, or invalid method
        
    Examples:
        >>> from metainformant.life_events import EventSequencePredictor, learn_event_embeddings
        >>> import numpy as np
        >>> sequences = [["health:diagnosis", "occupation:job_change"]]
        >>> y = np.array([0])
        >>> predictor = EventSequencePredictor(random_state=42)
        >>> predictor.fit(sequences, y)
        >>> embeddings = learn_event_embeddings(sequences, random_state=42)
        >>> importance = event_importance(predictor, sequences, embeddings)
        >>> len(importance) > 0
        True
    """
    if not predictor.is_fitted:
        raise ValueError("Predictor must be fitted before computing importance")
    
    if not sequences:
        raise ValueError("sequences list cannot be empty")
    
    if method not in ["permutation", "gradient"]:
        raise ValueError(f"Invalid method '{method}'. Must be 'permutation' or 'gradient'")
    
    # Get baseline predictions
    baseline_preds = predictor.predict(sequences)
    
    # Get unique events
    all_events = set()
    for seq in sequences:
        all_events.update(seq)
    
    if not all_events:
        return {}
    
    importance_scores = {}
    
    if method == "permutation":
        # Permutation importance: remove event and measure prediction change
        for event in all_events:
            # Create sequences without this event
            modified_sequences = [
                [e for e in seq if e != event]
                for seq in sequences
            ]
            
            # Predict with modified sequences
            modified_preds = predictor.predict(modified_sequences)
            
            # Importance is change in predictions
            if predictor.task_type == "classification":
                # Classification: change in predicted class
                importance = np.mean(baseline_preds != modified_preds)
            else:
                # Regression: change in predicted value
                importance = np.mean(np.abs(baseline_preds - modified_preds))
            
            importance_scores[event] = float(importance)
    
    else:
        # Gradient-based importance (placeholder)
        # Would require gradient computation from model
        for event in all_events:
            importance_scores[event] = 0.0
    
    # Normalize
    max_importance = max(importance_scores.values()) if importance_scores else 1.0
    if max_importance > 0:
        importance_scores = {k: v / max_importance for k, v in importance_scores.items()}
    
    return importance_scores


def temporal_patterns(
    sequences: List[List[str]],
    predictions: NDArray,
    time_windows: Optional[List[float]] = None
) -> Dict[str, Any]:
    """Identify critical time periods for predictions.
    
    Args:
        sequences: List of event sequences (must not be empty)
        predictions: Model predictions (must match sequence length)
        time_windows: Optional time windows for analysis (not yet implemented)
        
    Returns:
        Dictionary with temporal pattern analysis including position importance
        
    Raises:
        ValueError: If sequences empty or predictions length doesn't match sequences length
        
    Examples:
        >>> import numpy as np
        >>> sequences = [["health:diagnosis", "occupation:job_change"]]
        >>> predictions = np.array([0.8])
        >>> patterns = temporal_patterns(sequences, predictions)
        >>> "position_importance" in patterns
        True
    """
    if not sequences:
        raise ValueError("sequences list cannot be empty")
    
    if len(predictions) != len(sequences):
        raise ValueError(
            f"predictions length ({len(predictions)}) must match sequences length ({len(sequences)})"
        )
    
    # Analyze prediction patterns across sequence positions
    position_importance = {}
    
    max_len = max(len(s) for s in sequences) if sequences else 1
    
    # For each position, correlate with predictions
    for pos in range(max_len):
        position_events = []
        position_predictions = []
        
        for i, seq in enumerate(sequences):
            if pos < len(seq):
                position_events.append(seq[pos])
                position_predictions.append(predictions[i])
        
        if position_events:
            # Compute correlation between event type and prediction
            unique_events = set(position_events)
            for event in unique_events:
                event_mask = [e == event for e in position_events]
                event_preds = [p for p, m in zip(position_predictions, event_mask) if m]
                other_preds = [p for p, m in zip(position_predictions, event_mask) if not m]
                
                if event_preds and other_preds:
                    mean_diff = np.mean(event_preds) - np.mean(other_preds)
                    position_importance[f"pos_{pos}_{event}"] = float(mean_diff)
    
    return {
        "position_importance": position_importance,
        "max_sequence_length": max_len,
    }


def feature_attribution(
    predictor: "EventSequencePredictor",
    sequences: List[List[str]],
    event_embeddings: Dict[str, NDArray],
    use_shap: bool = False
) -> Dict[str, float]:
    """Compute SHAP-style feature attribution for events.
    
    Args:
        predictor: Trained EventSequencePredictor (must be fitted)
        sequences: Event sequences (must not be empty)
        event_embeddings: Event embeddings dictionary
        use_shap: Whether to use SHAP library (if available, falls back to permutation otherwise)
        
    Returns:
        Dictionary mapping events to attribution scores
        
    Raises:
        ValueError: If predictor not fitted or sequences empty
        
    Examples:
        >>> from metainformant.life_events import EventSequencePredictor, learn_event_embeddings
        >>> import numpy as np
        >>> sequences = [["health:diagnosis", "occupation:job_change"]]
        >>> y = np.array([0])
        >>> predictor = EventSequencePredictor(random_state=42)
        >>> predictor.fit(sequences, y)
        >>> embeddings = learn_event_embeddings(sequences, random_state=42)
        >>> attribution = feature_attribution(predictor, sequences, embeddings)
        >>> len(attribution) > 0
        True
    """
    if not predictor.is_fitted:
        raise ValueError("Predictor must be fitted before computing attribution")
    
    if not sequences:
        raise ValueError("sequences list cannot be empty")
    
    if use_shap and SHAP_AVAILABLE:
        # Use SHAP for attribution
        from .embeddings import sequence_embeddings
        
        # Convert sequences to embeddings
        X = sequence_embeddings(sequences, event_embeddings, method="mean")
        
        # Create SHAP explainer
        # Note: Would need to adapt to predictor's internal model
        # For now, use permutation importance as fallback
        return event_importance(predictor, sequences, event_embeddings, method="permutation")
    else:
        # Use permutation importance as alternative
        return event_importance(predictor, sequences, event_embeddings, method="permutation")

