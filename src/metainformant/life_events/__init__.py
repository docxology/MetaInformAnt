"""Life course event analysis and temporal modeling module for METAINFORMANT.

This module provides tools for analyzing life course trajectories, event sequence modeling,
temporal pattern recognition, and outcome prediction from longitudinal data.
"""

from __future__ import annotations

# Import all life events analysis submodules
from . import (
    config,
    embeddings,
    events,
    interpretability,
    models,
    utils,
    visualization,
    workflow,
)

# Optional imports with graceful fallbacks
try:
    from . import interpretability
except ImportError:
    interpretability = None

try:
    from . import visualization
except ImportError:
    visualization = None

# Direct imports of commonly used classes and functions
from .events import Event, EventSequence, EventDatabase
from .config import LifeEventsWorkflowConfig
from .models import EventSequencePredictor, EnsemblePredictor, GRUSequenceModel, LSTMSequenceModel, MultiTaskPredictor, SurvivalPredictor
from .interpretability import attention_weights, feature_attribution, temporal_patterns
from .embeddings import learn_event_embeddings, biological_embedding, domain_specific_embeddings
from .visualization import plot_event_timeline, plot_event_embeddings, plot_attention_heatmap, plot_prediction_importance, plot_domain_distribution, plot_domain_timeline, plot_embedding_clusters, plot_event_cooccurrence, plot_event_frequency_heatmap, plot_intervention_effects, plot_outcome_distribution, plot_population_comparison, plot_prediction_accuracy, plot_sequence_length_distribution, plot_sequence_similarity, plot_temporal_density, plot_temporal_patterns, plot_transition_network
from .workflow import analyze_life_course, compare_populations, intervention_analysis, event_importance
from .utils import add_temporal_noise, get_event_statistics, load_sequences_from_json, convert_sequences_to_tokens, validate_sequence, generate_cohort_sequences, sequence_embeddings, generate_event_chain, generate_synthetic_life_events
from .config import load_life_events_config

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core data structures and processing
    "events",
    "utils",
    "Event",
    "EventSequence",
    "EventDatabase",

    # Modeling and analysis
    "embeddings",
    "models",
    "workflow",
    "EventSequencePredictor",
    "EnsemblePredictor",
    "GRUSequenceModel",
    "LSTMSequenceModel",
    "MultiTaskPredictor",
    "SurvivalPredictor",
    "attention_weights",
    "feature_attribution",
    "temporal_patterns",
    "learn_event_embeddings",
    "biological_embedding",
    "domain_specific_embeddings",
    "analyze_life_course",
    "compare_populations",
    "intervention_analysis",
    "event_importance",

    # Configuration and visualization
    "config",
    "visualization",
    "LifeEventsWorkflowConfig",
    "load_life_events_config",
    "plot_event_timeline",
    "plot_event_embeddings",
    "plot_attention_heatmap",
    "plot_prediction_importance",
    "plot_domain_distribution",
    "plot_domain_timeline",
    "plot_embedding_clusters",
    "plot_event_cooccurrence",
    "plot_event_frequency_heatmap",
    "plot_intervention_effects",
    "plot_outcome_distribution",
    "plot_population_comparison",
    "plot_prediction_accuracy",
    "plot_sequence_length_distribution",
    "plot_sequence_similarity",
    "plot_temporal_density",
    "plot_temporal_patterns",
    "plot_transition_network",

    # Advanced analysis
    "interpretability",

    # Utilities
    "utils",
    "add_temporal_noise",
    "get_event_statistics",
    "load_sequences_from_json",
    "convert_sequences_to_tokens",
    "validate_sequence",
    "generate_cohort_sequences",
    "sequence_embeddings",
    "generate_event_chain",
    "generate_synthetic_life_events",
]



