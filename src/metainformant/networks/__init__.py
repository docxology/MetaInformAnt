"""Network analysis module for METAINFORMANT.

This module provides tools for biological network analysis, including
protein-protein interaction (PPI) networks, gene regulatory networks,
and pathway enrichment analysis.
"""

from __future__ import annotations

# Import subpackages
from . import analysis
from . import interaction

# Import modules from subpackages for backward compatibility
from .analysis import (
    community,
    graph,
    pathway,
)
from .interaction import (
    ppi,
    regulatory,
)

# Optional imports with graceful fallbacks
try:
    from .interaction import ppi
except ImportError:
    ppi = None

try:
    from .interaction import regulatory
except ImportError:
    regulatory = None

# Import regulatory module components (these may fail if optional dependencies not available)
try:
    from .interaction.regulatory import GeneRegulatoryNetwork, infer_grn
except ImportError:
    GeneRegulatoryNetwork = None
    infer_grn = None

# Import regulatory functions (these may fail if optional dependencies not available)
try:
    from .interaction.regulatory import regulatory_motifs, pathway_regulation_analysis
except ImportError:
    regulatory_motifs = None
    pathway_regulation_analysis = None

try:
    from . import pathway
except ImportError:
    pathway = None

# Direct imports of commonly used classes and functions
from .analysis.graph import BiologicalNetwork, centrality_measures
from .analysis.pathway import PathwayNetwork, pathway_enrichment, load_pathway_database, network_enrichment_analysis
from .interaction.ppi import ProteinNetwork, predict_interactions

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
