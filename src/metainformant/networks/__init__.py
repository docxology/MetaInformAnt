"""Biological network analysis module for METAINFORMANT.

This module provides comprehensive tools for biological network analysis,
including graph construction, community detection, centrality measures,
pathway analysis, protein-protein interactions, and regulatory networks.
"""

from __future__ import annotations

# Import all network analysis submodules
from . import (
    community,
    graph,
    pathway,
    ppi,
    regulatory,
)

# Optional imports with graceful fallbacks
# Import regulatory module components
try:
    from . import regulatory
    from .regulatory import GeneRegulatoryNetwork, infer_grn
except ImportError:
    regulatory = None
    GeneRegulatoryNetwork = None
    infer_grn = None

# Import regulatory functions (these may fail if optional dependencies not available)
try:
    from .regulatory import regulatory_motifs, pathway_regulation_analysis
except ImportError:
    regulatory_motifs = None
    pathway_regulation_analysis = None

try:
    from . import pathway
except ImportError:
    pathway = None

# Direct imports of commonly used classes and functions
from .graph import BiologicalNetwork, centrality_measures
from .pathway import PathwayNetwork, pathway_enrichment, load_pathway_database, network_enrichment_analysis
from .ppi import ProteinNetwork, predict_interactions
from .graph import BiologicalNetwork

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core network functionality
    "graph",
    "community",
    "BiologicalNetwork",
    "centrality_measures",

    # Biological network types
    "ppi",
    "pathway",
    "regulatory",
    "GeneRegulatoryNetwork",
    "infer_grn",
    "regulatory_motifs",
    "pathway_regulation_analysis",
    "PathwayNetwork",
    "ProteinNetwork",
    "pathway_enrichment",
    "load_pathway_database",
    "network_enrichment_analysis",
    "predict_interactions",
]



