"""Biological network analysis and graph algorithms.

This module provides tools for analyzing biological networks including:
- Network construction from various biological data
- Graph metrics and centrality measures
- Community detection and clustering
- Network visualization and analysis
- Pathway and gene regulatory networks
"""

from .community import (
    compare_communities,
    community_metrics,
    community_stability,
    detect_communities,
    hierarchical_communities,
    modularity,
    optimize_resolution,
)
from .graph import (
    BiologicalNetwork,
    add_edges_from_correlation,
    add_edges_from_interactions,
    centrality_measures,
    create_network,
    export_network,
    extract_subgraph,
    filter_network,
    get_connected_components,
    import_network,
    network_intersection,
    network_metrics,
    network_similarity,
    network_union,
    remove_edge,
    remove_node,
    shortest_paths,
)
from .pathway import (
    PathwayNetwork,
    load_pathway_database,
    network_enrichment_analysis,
    pathway_activity_score,
    pathway_enrichment,
    pathway_similarity,
)
from .ppi import (
    ProteinNetwork,
    detect_complexes,
    export_to_string_format,
    load_string_interactions,
    predict_interactions,
    protein_similarity,
)
from .regulatory import (
    GeneRegulatoryNetwork,
    detect_regulatory_cascades,
    infer_grn,
    regulatory_motifs,
    validate_regulation,
)

__all__ = [
    # Graph basics
    "BiologicalNetwork",
    "create_network",
    "add_edges_from_correlation",
    "add_edges_from_interactions",
    "network_metrics",
    "centrality_measures",
    "shortest_paths",
    # Network operations
    "export_network",
    "import_network",
    "network_similarity",
    "extract_subgraph",
    "filter_network",
    "get_connected_components",
    "network_union",
    "network_intersection",
    "remove_node",
    "remove_edge",
    # Community detection
    "detect_communities",
    "modularity",
    "community_metrics",
    "hierarchical_communities",
    "community_stability",
    "compare_communities",
    "optimize_resolution",
    # Pathways
    "PathwayNetwork",
    "load_pathway_database",
    "pathway_enrichment",
    "network_enrichment_analysis",
    "pathway_similarity",
    "pathway_activity_score",
    # Protein interactions
    "ProteinNetwork",
    "load_string_interactions",
    "predict_interactions",
    "get_protein_partners",  # Method on ProteinNetwork
    "protein_similarity",
    "detect_complexes",
    "export_to_string_format",
    # Gene regulation
    "GeneRegulatoryNetwork",
    "infer_grn",
    "regulatory_motifs",
    "detect_regulatory_cascades",
    "validate_regulation",
]
