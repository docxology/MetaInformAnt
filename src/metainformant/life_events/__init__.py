"""Life events and trajectory analysis module for METAINFORMANT.

This module provides tools for analyzing life events, trajectories,
and temporal patterns in biological data.
"""

from __future__ import annotations

# Import subpackages
from . import analysis
from . import core
from . import models
from . import visualization
from . import workflow

# Import modules from subpackages for backward compatibility
from .models import (
    embeddings,
    models as models_module,
)
from .core import (
    config,
    events,
    utils,
)
from .core.events import Event, EventSequence, EventDatabase

# Models
from .models.models import (
    EventSequencePredictor,
    EnsemblePredictor,
    GRUSequenceModel,
    LSTMSequenceModel,
    MultiTaskPredictor,
    SurvivalPredictor,
    attention_weights,
)
from .models.embeddings import (
    learn_event_embeddings,
    biological_embedding,
    domain_specific_embeddings,
)

# Analysis
from .analysis import (
    event_importance,
    feature_attribution,
)
from .analysis.interpretability import temporal_patterns

# Workflow
from .workflow import workflow as workflow_module
from .workflow.workflow import (
    analyze_life_course,
    compare_populations,
    intervention_analysis,
)

# Utils
from .core.utils import (
    load_sequences_from_json,
    get_event_statistics,
    add_temporal_noise,
    validate_sequence,
    generate_cohort_sequences,
    generate_synthetic_life_events,
    generate_event_chain,
    sequence_embeddings,
    convert_sequences_to_tokens,
)

# Config
from .core.config import LifeEventsWorkflowConfig, load_life_events_config

# Visualization
from .visualization.visualization import (
    plot_domain_distribution,
    plot_domain_timeline,
    plot_event_timeline,
    plot_event_embeddings,
    plot_attention_heatmap,
    plot_prediction_importance,
    plot_embedding_clusters,
    plot_event_cooccurrence,
    plot_event_frequency_heatmap,
    plot_intervention_effects,
    plot_outcome_distribution,
    plot_population_comparison,
    plot_prediction_accuracy,
    plot_sequence_length_distribution,
    plot_sequence_similarity,
    plot_temporal_density,
    plot_temporal_patterns,
    plot_transition_network,
)
from .visualization import visualization as visualization_module

# Optional imports with graceful fallbacks
try:
    from .models import models
except ImportError:
    models = None

try:
    from .models import embeddings
except ImportError:
    embeddings = None

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

    # Workflow Functions
    "analyze_life_course",
    "compare_populations",
    "intervention_analysis",

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
    "plot_sequence_length_distribution",
    "plot_sequence_similarity",
    "plot_temporal_density",
    "plot_temporal_patterns",
    "plot_transition_network",

    # Utilities
    "temporal_patterns",
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
