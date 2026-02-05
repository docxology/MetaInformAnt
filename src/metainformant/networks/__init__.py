"""Network analysis module for METAINFORMANT.

This module provides tools for biological network analysis, including
protein-protein interaction (PPI) networks, gene regulatory networks,
and pathway enrichment analysis.
"""

from __future__ import annotations

# Import subpackages
from . import analysis, interaction
from .analysis import community, graph, pathway
from .interaction import ppi, regulatory
from .config import (
    CommunityDetectionConfig,
    GRNConfig,
    NetworkConfig,
    NetworkWorkflowConfig,
    PPIConfig,
    PathwayEnrichmentConfig,
)
from .workflow import NetworkWorkflow

# Direct imports of commonly used classes and functions
from .analysis.graph import (
    BiologicalNetwork,
    centrality_measures,
    create_network,
    export_network,
    filter_network,
    get_connected_components,
    get_network_summary,
    import_network,
    load_network,
    network_intersection,
    network_metrics,
    network_similarity,
    network_union,
    save_network,
    shortest_paths,
)
from .analysis.community import (
    detect_communities,
    evaluate_communities,
    community_metrics,
    louvain_communities,
    greedy_modularity_communities,
    label_propagation_communities,
)
from .analysis.pathway import (
    PathwayNetwork,
    load_pathway_database,
    network_enrichment_analysis,
    pathway_enrichment,
    pathway_enrichment_analysis,
)
from .interaction.ppi import (
    ProteinNetwork,
    predict_interactions,
    load_ppi_network,
    ppi_network_analysis,
)

# Optional imports that may fail if dependencies missing
try:
    from .interaction.regulatory import (
        GeneRegulatoryNetwork,
        infer_grn,
        regulatory_motifs,
        pathway_regulation_analysis,
    )
except ImportError:
    GeneRegulatoryNetwork = None  # type: ignore[assignment,misc]
    infer_grn = None  # type: ignore[assignment]
    regulatory_motifs = None  # type: ignore[assignment]
    pathway_regulation_analysis = None  # type: ignore[assignment]

__all__ = [
    # Subpackages
    "analysis",
    "interaction",
    # Submodules
    "community",
    "graph",
    "pathway",
    "ppi",
    "regulatory",
    # Graph classes and functions
    "BiologicalNetwork",
    "centrality_measures",
    "create_network",
    "export_network",
    "filter_network",
    "get_connected_components",
    "get_network_summary",
    "import_network",
    "load_network",
    "network_intersection",
    "network_metrics",
    "network_similarity",
    "network_union",
    "save_network",
    "shortest_paths",
    # Community detection
    "detect_communities",
    "evaluate_communities",
    "community_metrics",
    "louvain_communities",
    "greedy_modularity_communities",
    "label_propagation_communities",
    # Pathway analysis
    "PathwayNetwork",
    "load_pathway_database",
    "network_enrichment_analysis",
    "pathway_enrichment",
    "pathway_enrichment_analysis",
    # PPI analysis
    "ProteinNetwork",
    "predict_interactions",
    "load_ppi_network",
    "ppi_network_analysis",
    # Regulatory networks
    "GeneRegulatoryNetwork",
    "infer_grn",
    "regulatory_motifs",
    "pathway_regulation_analysis",
    # Config and workflow
    "CommunityDetectionConfig",
    "GRNConfig",
    "NetworkConfig",
    "NetworkWorkflowConfig",
    "PPIConfig",
    "PathwayEnrichmentConfig",
    "NetworkWorkflow",
]
