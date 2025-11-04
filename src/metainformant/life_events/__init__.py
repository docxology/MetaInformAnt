"""Life course and event sequence analysis module.

This module provides tools for analyzing human life courses as temporal event
sequences, enabling prediction of life outcomes from event data using NLP-inspired
methods.

Key Components:
- Event data structures: Event, EventSequence, EventDatabase
- Embeddings: Learn dense vector representations of events
- Models: Sequence prediction models for classification and regression
- Workflows: End-to-end analysis pipelines
- Visualization: Plotting functions for event sequences and embeddings
- Interpretability: Model interpretation and feature importance tools

All functions use real implementations without mocking. Optional dependencies
are handled defensively.
"""

from __future__ import annotations

from .events import Event, EventDatabase, EventSequence
from .embeddings import (
    domain_specific_embeddings,
    learn_event_embeddings,
    sequence_embeddings,
)
from .models import EventSequencePredictor, LSTMSequenceModel
from .workflow import analyze_life_course, compare_populations, intervention_analysis

# Import visualization functions defensively (matplotlib may not always be available)
try:
    from .visualization import (
        plot_attention_heatmap,
        plot_event_embeddings,
        plot_event_timeline,
        plot_prediction_importance,
    )

    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False
    plot_event_timeline = None  # type: ignore
    plot_event_embeddings = None  # type: ignore
    plot_attention_heatmap = None  # type: ignore
    plot_prediction_importance = None  # type: ignore

# Import interpretability functions (no optional deps required for basic functionality)
from .interpretability import (
    attention_weights,
    event_importance,
    feature_attribution,
    temporal_patterns,
)

__all__ = [
    # Event data structures
    "Event",
    "EventSequence",
    "EventDatabase",
    # Embeddings
    "learn_event_embeddings",
    "sequence_embeddings",
    "domain_specific_embeddings",
    # Models
    "EventSequencePredictor",
    "LSTMSequenceModel",
    # Workflows
    "analyze_life_course",
    "compare_populations",
    "intervention_analysis",
    # Visualization
    "plot_event_timeline",
    "plot_event_embeddings",
    "plot_attention_heatmap",
    "plot_prediction_importance",
    # Interpretability
    "attention_weights",
    "event_importance",
    "temporal_patterns",
    "feature_attribution",
]

