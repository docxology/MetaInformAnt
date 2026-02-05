"""Network analysis module for METAINFORMANT.

This module provides tools for biological network analysis, including
protein-protein interaction (PPI) networks, gene regulatory networks,
and pathway enrichment analysis.
"""

from __future__ import annotations

# Import subpackages
from . import analysis, interaction, regulatory
from .analysis import community, graph, pathway
from .analysis.community import (
    community_metrics,
    detect_communities,
    evaluate_communities,
    greedy_modularity_communities,
    label_propagation_communities,
    louvain_communities,
)

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
from .analysis.pathway import (
    PathwayNetwork,
    load_pathway_database,
    network_enrichment_analysis,
    pathway_enrichment,
    pathway_enrichment_analysis,
)
from .config import (
    CommunityDetectionConfig,
    GRNConfig,
    NetworkConfig,
    NetworkWorkflowConfig,
    PathwayEnrichmentConfig,
    PPIConfig,
)
from .interaction import ppi
from .interaction import regulatory as interaction_regulatory
from .interaction.ppi import (
    ProteinNetwork,
    load_ppi_network,
    ppi_network_analysis,
    predict_interactions,
)
from .workflow import NetworkWorkflow

# Optional imports that may fail if dependencies missing
try:
    from .interaction.regulatory import (
        GeneRegulatoryNetwork,
        infer_grn,
        pathway_regulation_analysis,
        regulatory_motifs,
    )
except ImportError:
    GeneRegulatoryNetwork = None  # type: ignore[assignment,misc]
    infer_grn = None  # type: ignore[assignment]
    regulatory_motifs = None  # type: ignore[assignment]
    pathway_regulation_analysis = None  # type: ignore[assignment]

# Regulatory network inference subpackage
try:
    from .regulatory.grn_inference import (
        compute_network_motifs,
        infer_grn_correlation,
        infer_grn_mutual_info,
        infer_grn_regression,
        score_regulators,
        validate_grn,
    )
    from .regulatory.motif_analysis import (
        build_pwm,
        find_tf_binding_motifs,
        scan_sequence_for_motifs,
        score_motif_match,
    )
except ImportError:
    infer_grn_correlation = None  # type: ignore[assignment]
    infer_grn_mutual_info = None  # type: ignore[assignment]
    infer_grn_regression = None  # type: ignore[assignment]
    score_regulators = None  # type: ignore[assignment]
    compute_network_motifs = None  # type: ignore[assignment]
    validate_grn = None  # type: ignore[assignment]
    find_tf_binding_motifs = None  # type: ignore[assignment]
    score_motif_match = None  # type: ignore[assignment]
    build_pwm = None  # type: ignore[assignment]
    scan_sequence_for_motifs = None  # type: ignore[assignment]

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
    # Regulatory network inference (new subpackage)
    "infer_grn_correlation",
    "infer_grn_mutual_info",
    "infer_grn_regression",
    "score_regulators",
    "compute_network_motifs",
    "validate_grn",
    "find_tf_binding_motifs",
    "score_motif_match",
    "build_pwm",
    "scan_sequence_for_motifs",
    # Config and workflow
    "CommunityDetectionConfig",
    "GRNConfig",
    "NetworkConfig",
    "NetworkWorkflowConfig",
    "PPIConfig",
    "PathwayEnrichmentConfig",
    "NetworkWorkflow",
]
