"""Life events visualization integration.

This module provides unified access to life events visualization functions
from the life_events module, integrating them into the central visualization package.
"""

from __future__ import annotations

# Re-export life events visualization functions
try:
    from ..life_events.visualization import (
        plot_attention_heatmap,
        plot_event_embeddings,
        plot_event_timeline,
        plot_prediction_importance,
    )

    LIFE_EVENTS_VISUALIZATION_AVAILABLE = True
except ImportError:
    LIFE_EVENTS_VISUALIZATION_AVAILABLE = False
    # Define placeholder functions
    def _not_available(*args, **kwargs):
        raise ImportError("Life events visualization functions not available. Install life_events module dependencies.")
    
    plot_event_timeline = _not_available
    plot_event_embeddings = _not_available
    plot_attention_heatmap = _not_available
    plot_prediction_importance = _not_available

__all__ = [
    "plot_event_timeline",
    "plot_event_embeddings",
    "plot_attention_heatmap",
    "plot_prediction_importance",
    "LIFE_EVENTS_VISUALIZATION_AVAILABLE",
]

