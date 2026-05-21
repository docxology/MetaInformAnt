"""Life events prediction models module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from . import embeddings, predictor, sequence_models, statistical_models
from .embeddings import biological_embedding, domain_specific_embeddings, learn_event_embeddings
from .predictor import EventSequencePredictor


def train_event_predictor(
    sequences: list[Any],
    outcomes: list[Any],
    *,
    embedding_results: dict[str, Any] | None = None,
    model_type: str = "embedding",
    task_type: str = "classification",
    embedding_dim: int = 100,
    random_state: int = 42,
    **kwargs: Any,
) -> dict[str, Any]:
    """Train an EventSequencePredictor and return workflow-style metadata."""
    model = EventSequencePredictor(
        model_type=model_type,
        task_type=task_type,
        embedding_dim=embedding_dim,
        random_state=random_state,
        **kwargs,
    )
    model.fit(sequences, outcomes, event_embeddings=embedding_results)
    predictions = model.predict(sequences)
    return {
        "model": model,
        "predictions": predictions.tolist() if hasattr(predictions, "tolist") else list(predictions),
        "model_type": model_type,
        "task_type": task_type,
        "embedding_dim": embedding_dim,
    }


def predict_outcomes(sequences: list[Any], model: Any, **_: Any) -> list[Any]:
    """Predict outcomes with a trained model or training-result dictionary."""
    predictor_obj = model.get("model") if isinstance(model, dict) else model
    predictions = predictor_obj.predict(sequences)
    return predictions.tolist() if hasattr(predictions, "tolist") else list(predictions)


def save_model(model: Any, path: str | Path) -> None:
    """Save a trained model or training-result dictionary to disk."""
    predictor_obj = model.get("model") if isinstance(model, dict) else model
    predictor_obj.save_model(path)


__all__ = [
    "embeddings",
    "predictor",
    "sequence_models",
    "statistical_models",
    "EventSequencePredictor",
    "learn_event_embeddings",
    "biological_embedding",
    "domain_specific_embeddings",
    "train_event_predictor",
    "predict_outcomes",
    "save_model",
]
