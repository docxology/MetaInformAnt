"""Information theory integration with network analysis.

This module provides information-theoretic measures for network analysis,
including information flow, network entropy, and information-based
community detection.
"""

from __future__ import annotations

from typing import Any

import numpy as np

try:
    import networkx as nx
except ImportError:
    nx = None

from metainformant.information.syntactic import (
    conditional_entropy,
    joint_entropy,
    mutual_information,
    shannon_entropy,
    shannon_entropy_from_counts,
)


def network_entropy(
    graph: Any,
    attribute: str | None = None
) -> float:
    """Calculate information entropy of a network.
    
    Measures the information content of the network structure or
    node attributes.
    
    Args:
        graph: NetworkX graph or network object
        attribute: Optional node attribute to analyze
        
    Returns:
        Network entropy in bits
    """
    if nx is None:
        raise ImportError("NetworkX not available")
    
    if not hasattr(graph, "nodes"):
        raise ValueError("Graph must be a NetworkX graph or compatible object")
    
    if attribute:
        # Calculate entropy of attribute distribution
        attr_values = [graph.nodes[node].get(attribute) for node in graph.nodes()]
        from collections import Counter
        
        counts = Counter(attr_values)
        return shannon_entropy_from_counts(counts)
    else:
        # Calculate entropy of degree distribution
        degrees = [graph.degree(node) for node in graph.nodes()]
        from collections import Counter
        
        counts = Counter(degrees)
        return shannon_entropy_from_counts(counts)


def information_flow(
    graph: Any,
    source_nodes: list[str] | None = None,
    target_nodes: list[str] | None = None
) -> dict[str, Any]:
    """Calculate information flow through network edges.
    
    Measures how information propagates through the network structure.
    
    Args:
        graph: NetworkX graph
        source_nodes: Optional list of source nodes
        target_nodes: Optional list of target nodes
        
    Returns:
        Dictionary with information flow metrics
    """
    if nx is None:
        raise ImportError("NetworkX not available")
    
    if not hasattr(graph, "nodes"):
        raise ValueError("Graph must be a NetworkX graph or compatible object")
    
    # Use degree as proxy for information flow
    if source_nodes is None:
        source_nodes = list(graph.nodes())[:10]  # Sample
    
    if target_nodes is None:
        target_nodes = list(graph.nodes())
    
    # Calculate shortest paths (information distance)
    path_lengths = []
    for source in source_nodes:
        try:
            lengths = nx.single_source_shortest_path_length(graph, source)
            for target in target_nodes:
                if target in lengths:
                    path_lengths.append(lengths[target])
        except Exception:
            continue
    
    # Calculate entropy of path length distribution
    if path_lengths:
        from collections import Counter
        
        counts = Counter(path_lengths)
        entropy = shannon_entropy_from_counts(counts)
    else:
        entropy = 0.0
    
    return {
        "path_length_entropy": entropy,
        "mean_path_length": float(np.mean(path_lengths)) if path_lengths else 0.0,
        "num_paths": len(path_lengths),
    }


def mutual_information_network(
    graph: Any,
    node_attributes: dict[str, list[float]] | None = None
) -> np.ndarray:
    """Calculate mutual information matrix from network structure.
    
    Uses network topology to estimate mutual information between nodes.
    
    Args:
        graph: NetworkX graph
        node_attributes: Optional dictionary mapping node IDs to attribute vectors
        
    Returns:
        Mutual information matrix
    """
    if nx is None:
        raise ImportError("NetworkX not available")
    
    if not hasattr(graph, "nodes"):
        raise ValueError("Graph must be a NetworkX graph or compatible object")
    
    nodes = list(graph.nodes())
    n = len(nodes)
    mi_matrix = np.zeros((n, n))
    
    if node_attributes:
        # Calculate MI from attributes
        for i, node_i in enumerate(nodes):
            if node_i not in node_attributes:
                continue
            for j, node_j in enumerate(nodes):
                if node_j not in node_attributes:
                    continue
                if i != j:
                    attr_i = node_attributes[node_i]
                    attr_j = node_attributes[node_j]
                    # Convert to discrete for MI calculation
                    attr_i_discrete = [int(x) for x in attr_i]
                    attr_j_discrete = [int(y) for y in attr_j]
                    mi = mutual_information(attr_i_discrete, attr_j_discrete)
                    mi_matrix[i, j] = mi
    else:
        # Use network distance as proxy for MI
        # Closer nodes have higher MI
        for i, node_i in enumerate(nodes):
            for j, node_j in enumerate(nodes):
                if i != j:
                    try:
                        path_length = nx.shortest_path_length(graph, node_i, node_j)
                        # Inverse relationship: shorter paths = higher MI
                        mi_matrix[i, j] = 1.0 / (1.0 + path_length)
                    except Exception:
                        mi_matrix[i, j] = 0.0
    
    return mi_matrix


def information_based_communities(
    graph: Any,
    resolution: float = 1.0
) -> dict[str, list[str]]:
    """Detect communities using information-theoretic measures.
    
    Uses mutual information to identify communities with high
    internal information flow.
    
    Args:
        graph: NetworkX graph
        resolution: Resolution parameter for community detection
        
    Returns:
        Dictionary mapping community ID to list of node IDs
    """
    if nx is None:
        raise ImportError("NetworkX not available")
    
    try:
        from networkx.algorithms import community
        
        # Use Louvain algorithm as base
        communities = community.louvain_communities(graph, resolution=resolution)
        
        # Convert to dictionary format
        result = {}
        for i, comm in enumerate(communities):
            result[f"community_{i}"] = list(comm)
        
        return result
    except Exception:
        # Fallback: simple connected components
        components = list(nx.connected_components(graph))
        result = {}
        for i, comp in enumerate(components):
            result[f"component_{i}"] = list(comp)
        return result

