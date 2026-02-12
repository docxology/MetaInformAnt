"""Configuration classes for network analysis.

Provides typed configuration for network construction, community detection,
pathway enrichment, and workflow orchestration.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class NetworkConfig:
    """Base configuration for network analysis."""

    directed: bool = False
    weighted: bool = True
    min_edge_weight: float = 0.0
    max_nodes: Optional[int] = None
    remove_isolates: bool = False
    random_seed: Optional[int] = None


@dataclass
class CommunityDetectionConfig:
    """Configuration for community detection."""

    method: str = "greedy"
    resolution: float = 1.0
    n_communities: Optional[int] = None
    max_levels: int = 5
    random_state: Optional[int] = None
    compare_methods: bool = False
    methods_to_compare: List[str] = field(default_factory=lambda: ["greedy", "label_propagation"])


@dataclass
class PathwayEnrichmentConfig:
    """Configuration for pathway enrichment analysis."""

    method: str = "fisher"
    correction: str = "bonferroni"
    min_overlap: int = 1
    max_p_value: float = 0.05
    min_pathway_size: int = 5
    max_pathway_size: Optional[int] = 500
    background_genes: Optional[List[str]] = None


@dataclass
class PPIConfig:
    """Configuration for PPI network analysis."""

    confidence_threshold: float = 0.4
    hub_percentile: float = 95.0
    max_nodes_for_centrality: int = 1000
    prediction_method: str = "similarity"
    prediction_threshold: float = 0.5
    max_predictions_per_protein: int = 10


@dataclass
class GRNConfig:
    """Configuration for gene regulatory network inference."""

    method: str = "correlation"
    threshold: float = 0.5
    tf_genes: Optional[List[str]] = None
    motif_types: List[str] = field(default_factory=lambda: ["feed_forward", "feedback", "cascade"])


@dataclass
class NetworkWorkflowConfig:
    """Configuration for network analysis workflow."""

    network: NetworkConfig = field(default_factory=NetworkConfig)
    community: CommunityDetectionConfig = field(default_factory=CommunityDetectionConfig)
    pathway: PathwayEnrichmentConfig = field(default_factory=PathwayEnrichmentConfig)
    ppi: PPIConfig = field(default_factory=PPIConfig)
    grn: GRNConfig = field(default_factory=GRNConfig)
    output_dir: Optional[str] = None
    export_format: str = "json"
    verbose: bool = False

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "NetworkWorkflowConfig":
        """Create config from dictionary.

        Args:
            data: Configuration dictionary

        Returns:
            NetworkWorkflowConfig instance
        """
        network = NetworkConfig(**data.get("network", {}))
        community = CommunityDetectionConfig(**data.get("community", {}))
        pathway = PathwayEnrichmentConfig(**data.get("pathway", {}))
        ppi = PPIConfig(**data.get("ppi", {}))
        grn = GRNConfig(**data.get("grn", {}))
        return cls(
            network=network,
            community=community,
            pathway=pathway,
            ppi=ppi,
            grn=grn,
            output_dir=data.get("output_dir"),
            export_format=data.get("export_format", "json"),
            verbose=data.get("verbose", False),
        )
