"""Chromatin state learning submodule for epigenome analysis.

This submodule provides tools for learning and analyzing chromatin states
from histone modification data, including ChromHMM-style state discovery
using Gaussian mixture models, state assignment, biological interpretation,
enrichment analysis, and cross-condition comparison.
"""

from __future__ import annotations

from .state_learning import (
    assign_states,
    compare_chromatin_states,
    compute_state_enrichment,
    compute_state_transition_rates,
    interpret_states,
    learn_chromatin_states,
    segment_genome,
)

__all__ = [
    "learn_chromatin_states",
    "assign_states",
    "interpret_states",
    "compute_state_enrichment",
    "segment_genome",
    "compare_chromatin_states",
    "compute_state_transition_rates",
]
