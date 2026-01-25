"""Life events and trajectory analysis module for METAINFORMANT.

This module provides comprehensive tools for analyzing life events, trajectories,
and temporal patterns in biological and longitudinal data. It supports event sequence
modeling, prediction, and interpretation.

Capabilities:
    - **Event Representation**: Model life events with timestamps, types, and metadata
    - **Sequence Analysis**: Analyze temporal sequences of events
    - **Prediction Models**: Train models to predict future events or outcomes
    - **Embedding Learning**: Learn vector representations of events
    - **Interpretability**: Explain model predictions using attention and attribution

Key Classes:
    - Event: Single life event with timestamp and type
    - EventSequence: Ordered sequence of events with filtering/aggregation
    - EventDatabase: Collection of sequences with querying capabilities
    - EventSequencePredictor: ML model for event-based predictions
    - EnsemblePredictor: Combine multiple predictors
    - GRUSequenceModel, LSTMSequenceModel: Neural sequence models

Key Functions:
    - analyze_life_course: Analyze patterns in life course data
    - compare_populations: Compare event distributions between groups
    - intervention_analysis: Analyze effects of interventions
    - temporal_patterns: Discover temporal patterns in sequences
    - learn_event_embeddings: Train word2vec-style event embeddings

Submodules:
    - core.events: Event and EventSequence data structures
    - core.config: Configuration for life events analysis
    - core.utils: Utility functions (loading, validation, generation)
    - models.models: Prediction models (sklearn, PyTorch)
    - models.embeddings: Event embedding methods
    - analysis: Interpretability and pattern analysis
    - workflow: High-level analysis workflows
    - visualization: Timeline, attention, and trajectory plots

Typical Workflow:
    1. Create Event objects with timestamps and types
    2. Build EventSequence objects from events
    3. Learn embeddings or use domain-specific embeddings
    4. Train prediction model (e.g., EventSequencePredictor)
    5. Interpret predictions with attention_weights or feature_attribution
    6. Visualize with plot_event_timeline, plot_attention_heatmap, etc.

Example:
    >>> from metainformant.life_events import Event, EventSequence, EventSequencePredictor
    >>> # Create events
    >>> events = [
    ...     Event(timestamp=0, event_type="diagnosis"),
    ...     Event(timestamp=30, event_type="treatment"),
    ...     Event(timestamp=90, event_type="followup"),
    ... ]
    >>> sequence = EventSequence(events=events, subject_id="patient_001")
    >>>
    >>> # Train predictor
    >>> predictor = EventSequencePredictor()
    >>> predictor.fit(sequences, labels)
    >>> predictions = predictor.predict(new_sequences)

See Also:
    - docs/life_events/ for detailed documentation
    - metainformant.networks for event network analysis
    - metainformant.visualization for additional plotting
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
