"""Life events visualization module."""

from __future__ import annotations

from . import network, statistical, timeline
from .statistical import (  # noqa: F401  (re-exports for workflow facade)
    plot_intervention_effects,
    plot_outcome_distribution,
    plot_population_comparison,
    plot_sequence_length_distribution,
)
from .timeline import plot_domain_distribution  # noqa: F401  (re-export)

__all__ = [
    "network",
    "statistical",
    "timeline",
    "plot_intervention_effects",
    "plot_outcome_distribution",
    "plot_population_comparison",
    "plot_sequence_length_distribution",
    "plot_domain_distribution",
]

__all__ = ["network", "statistical", "timeline"]
