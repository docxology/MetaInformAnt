"""Graph construction and manipulation utilities for METAINFORMANT.

This module provides comprehensive tools for creating, manipulating, and analyzing
biological networks including protein-protein interaction networks, gene regulatory
networks, and other biological graph structures.

This module re-exports all public symbols from :mod:`graph_core` and
:mod:`graph_algorithms` for backward compatibility.
"""

from __future__ import annotations

from metainformant.networks.analysis.graph_core import (
    BiologicalNetwork,
    add_edge_attributes,
    add_edges_from_correlation,
    add_edges_from_dataframe,
    add_edges_from_interactions,
    add_node_attributes,
    add_nodes_from_dataframe,
    convert_from_adjacency_matrix,
    convert_to_adjacency_matrix,
    create_network,
    create_subgraph,
    get_edge_attributes_dataframe,
    get_network_summary,
    get_node_attributes_dataframe,
    load_network,
    remove_isolated_nodes,
    save_network,
    validate_network,
)
from metainformant.networks.analysis.graph_algorithms import (
    centrality_measures,
    export_network,
    extract_subgraph,
    filter_network,
    get_connected_components,
    import_network,
    network_intersection,
    network_metrics,
    network_similarity,
    network_union,
    shortest_paths,
    subgraph,
)

__all__ = [
    # Core class
    "BiologicalNetwork",
    # Core construction/IO
    "create_network",
    "load_network",
    "save_network",
    "add_nodes_from_dataframe",
    "add_edges_from_dataframe",
    "create_subgraph",
    "remove_isolated_nodes",
    "add_node_attributes",
    "add_edge_attributes",
    "get_node_attributes_dataframe",
    "get_edge_attributes_dataframe",
    "convert_to_adjacency_matrix",
    "convert_from_adjacency_matrix",
    "get_network_summary",
    "validate_network",
    "add_edges_from_interactions",
    "add_edges_from_correlation",
    # Algorithms
    "network_metrics",
    "subgraph",
    "export_network",
    "import_network",
    "network_similarity",
    "extract_subgraph",
    "filter_network",
    "get_connected_components",
    "network_union",
    "network_intersection",
    "centrality_measures",
    "shortest_paths",
]
