"""Regulatory network analysis and modeling.

This module provides tools for analyzing and modeling gene regulatory networks,
including network inference, motif discovery, and dynamic analysis.

This module re-exports all public symbols from :mod:`regulatory_core` and
:mod:`regulatory_analysis` for backward compatibility.
"""

from __future__ import annotations

from metainformant.networks.interaction.regulatory_core import (
    GeneRegulatoryNetwork,
    _find_feed_forward_loops,
    _find_feedback_loops,
    _find_regulatory_cascades,
    analyze_regulatory_dynamics,
    analyze_regulatory_motifs,
    calculate_regulatory_influence,
    construct_regulatory_network,
    export_regulatory_network,
    identify_regulatory_hubs,
    regulatory_network_stability_analysis,
)
from metainformant.networks.interaction.regulatory_analysis import (
    detect_regulatory_cascades,
    infer_grn,
    infer_regulatory_network_from_expression,
    pathway_regulation_analysis,
    regulatory_motifs,
    validate_regulation,
)

__all__ = [
    # Core classes
    "GeneRegulatoryNetwork",
    # Core functions
    "construct_regulatory_network",
    "analyze_regulatory_motifs",
    "calculate_regulatory_influence",
    "analyze_regulatory_dynamics",
    "identify_regulatory_hubs",
    "regulatory_network_stability_analysis",
    "export_regulatory_network",
    # Analysis functions
    "infer_regulatory_network_from_expression",
    "infer_grn",
    "regulatory_motifs",
    "detect_regulatory_cascades",
    "validate_regulation",
    "pathway_regulation_analysis",
    # Private helpers (still accessible for backward compat)
    "_find_feed_forward_loops",
    "_find_feedback_loops",
    "_find_regulatory_cascades",
]
