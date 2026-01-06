"""Network information analysis for METAINFORMANT.

This module provides information-theoretic measures for analyzing
biological networks, including network entropy, information flow,
and information-based community detection.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Optional network analysis dependencies
try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    logger.warning("networkx not available, network information analysis disabled")


def network_entropy(
    graph: Any,
    attribute: Optional[str] = None,
    method: str = "shannon"
) -> float:
    """Calculate entropy of a biological network.

    Args:
        graph: NetworkX graph or adjacency matrix
        attribute: Node attribute to use for entropy calculation
        method: Entropy calculation method ('shannon', 'von_neumann')

    Returns:
        Network entropy value

    Raises:
        ImportError: If networkx not available
        ValueError: If graph format not supported
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network entropy calculation")

    # Convert to NetworkX graph if needed
    if isinstance(graph, np.ndarray):
        G = nx.from_numpy_array(graph)
    elif hasattr(graph, 'nodes'):  # NetworkX graph
        G = graph
    else:
        raise ValueError("Graph must be NetworkX graph or numpy adjacency matrix")

    if method == "shannon":
        return _shannon_network_entropy(G, attribute)
    elif method == "von_neumann":
        return _von_neumann_entropy(G)
    else:
        raise ValueError(f"Unknown entropy method: {method}")


def _shannon_network_entropy(G: Any, attribute: Optional[str] = None) -> float:
    """Calculate Shannon entropy of network structure or attributes."""
    if attribute is None:
        # Structural entropy based on degree distribution
        degrees = [d for n, d in G.degree()]
        from collections import Counter
        degree_counts = Counter(degrees)

        # Convert to probabilities
        total_nodes = len(G)
        probs = [count / total_nodes for count in degree_counts.values()]

        # Calculate entropy
        entropy = 0.0
        for p in probs:
            if p > 0:
                entropy -= p * math.log2(p)

        return entropy

    else:
        # Attribute-based entropy
        if not nx.get_node_attributes(G, attribute):
            raise ValueError(f"Attribute '{attribute}' not found in graph nodes")

        attr_values = list(nx.get_node_attributes(G, attribute).values())
        from collections import Counter
        attr_counts = Counter(attr_values)

        # Calculate entropy
        total = sum(attr_counts.values())
        entropy = 0.0
        for count in attr_counts.values():
            p = count / total
            entropy -= p * math.log2(p)

        return entropy


def _von_neumann_entropy(G: Any) -> float:
    """Calculate von Neumann graph entropy."""
    # Convert to adjacency matrix
    if hasattr(G, 'to_numpy_array'):
        A = G.to_numpy_array()
    else:
        A = nx.to_numpy_array(G)

    # Normalize adjacency matrix
    degree_sum = np.sum(A, axis=1)
    degree_sum = np.where(degree_sum == 0, 1, degree_sum)  # Avoid division by zero

    # Create transition matrix (row-normalized adjacency)
    P = A / degree_sum[:, np.newaxis]

    # Calculate von Neumann entropy: -Tr(P log P)
    # For computational efficiency, use eigenvalues
    try:
        eigenvals = np.linalg.eigvals(P)
        eigenvals = eigenvals[eigenvals > 0]  # Only positive eigenvalues

        entropy = 0.0
        for val in eigenvals:
            entropy -= val * math.log2(val)

        return entropy.real  # Return real part

    except np.linalg.LinAlgError:
        logger.warning("Eigenvalue calculation failed, using approximation")
        # Fallback: use trace approximation
        return -np.trace(P @ np.log2(np.maximum(P, 1e-10)))


def information_flow(
    graph: Any,
    source_nodes: Optional[List[str]] = None,
    target_nodes: Optional[List[str]] = None,
    method: str = "random_walk",
    steps: int = 100
) -> Dict[str, Any]:
    """Calculate information flow in a biological network.

    Args:
        graph: NetworkX graph or adjacency matrix
        source_nodes: Source nodes for information flow
        target_nodes: Target nodes for information flow
        method: Information flow method ('random_walk', 'diffusion')
        steps: Number of steps for flow calculation

    Returns:
        Dictionary with information flow results

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for information flow calculation")

    # Convert to NetworkX graph if needed
    if isinstance(graph, np.ndarray):
        G = nx.from_numpy_array(graph)
    elif hasattr(graph, 'nodes'):
        G = graph
    else:
        raise ValueError("Graph must be NetworkX graph or numpy adjacency matrix")

    if method == "random_walk":
        return _random_walk_information_flow(G, source_nodes, target_nodes, steps)
    elif method == "diffusion":
        return _diffusion_information_flow(G, source_nodes, target_nodes, steps)
    else:
        raise ValueError(f"Unknown flow method: {method}")


def _random_walk_information_flow(
    G: Any,
    source_nodes: Optional[List[str]] = None,
    target_nodes: Optional[List[str]] = None,
    steps: int = 100
) -> Dict[str, Any]:
    """Calculate information flow using random walk model."""
    # Create transition matrix
    A = nx.to_numpy_array(G)
    degrees = np.sum(A, axis=1)
    degrees = np.where(degrees == 0, 1, degrees)  # Avoid division by zero

    P = A / degrees[:, np.newaxis]  # Transition matrix

    n_nodes = len(G)
    flow_matrix = np.zeros((n_nodes, n_nodes))

    # Initialize with source nodes
    if source_nodes is None:
        # Use all nodes as sources
        initial_dist = np.ones(n_nodes) / n_nodes
    else:
        # Map node names to indices
        node_list = list(G.nodes())
        node_to_idx = {node: i for i, node in enumerate(node_list)}

        initial_dist = np.zeros(n_nodes)
        for source in source_nodes:
            if source in node_to_idx:
                initial_dist[node_to_idx[source]] = 1.0 / len(source_nodes)

    # Simulate random walk
    current_dist = initial_dist.copy()

    for step in range(steps):
        # Update distribution
        current_dist = P.T @ current_dist

        # Record flow to target nodes
        if target_nodes:
            node_list = list(G.nodes())
            node_to_idx = {node: i for i, node in enumerate(node_list)}

            for i, source_idx in enumerate(np.where(initial_dist > 0)[0]):
                source_node = node_list[source_idx]
                if source_node in source_nodes:
                    for target in target_nodes:
                        if target in node_to_idx:
                            target_idx = node_to_idx[target]
                            flow_matrix[source_idx, target_idx] += current_dist[target_idx]

    # Convert back to node names
    node_list = list(G.nodes())
    flow_dict = {}

    for i, source in enumerate(node_list):
        if source_nodes is None or source in source_nodes:
            flow_dict[source] = {}
            for j, target in enumerate(node_list):
                if target_nodes is None or target in target_nodes:
                    flow_dict[source][target] = flow_matrix[i, j]

    return {
        'flow_matrix': flow_dict,
        'method': 'random_walk',
        'steps': steps,
        'source_nodes': source_nodes,
        'target_nodes': target_nodes,
    }


def _diffusion_information_flow(
    G: Any,
    source_nodes: Optional[List[str]] = None,
    target_nodes: Optional[List[str]] = None,
    steps: int = 100
) -> Dict[str, Any]:
    """Calculate information flow using diffusion model."""
    A = nx.to_numpy_array(G)
    n_nodes = len(G)

    # Create diffusion matrix (normalized Laplacian)
    degrees = np.sum(A, axis=1)
    degrees = np.where(degrees == 0, 1, degrees)

    # Laplacian matrix
    L = np.diag(degrees) - A

    # Normalized Laplacian (for regular graphs, this approximates diffusion)
    D_inv_sqrt = np.diag(1.0 / np.sqrt(degrees))
    L_norm = D_inv_sqrt @ L @ D_inv_sqrt

    # Initialize with source nodes
    if source_nodes is None:
        initial_signal = np.ones(n_nodes) / n_nodes
    else:
        node_list = list(G.nodes())
        node_to_idx = {node: i for i, node in enumerate(node_list)}

        initial_signal = np.zeros(n_nodes)
        for source in source_nodes:
            if source in node_to_idx:
                initial_signal[node_to_idx[source]] = 1.0 / len(source_nodes)

    # Diffusion simulation
    signal = initial_signal.copy()
    flow_history = []

    for step in range(steps):
        # Diffuse signal
        signal = signal - 0.1 * L_norm @ signal  # Simple diffusion step
        signal = np.maximum(signal, 0)  # Ensure non-negative
        signal = signal / np.sum(signal)  # Renormalize

        flow_history.append(signal.copy())

    # Calculate final flow patterns
    final_signal = flow_history[-1] if flow_history else initial_signal

    # Convert to node-wise results
    node_list = list(G.nodes())
    flow_results = {}

    for i, node in enumerate(node_list):
        flow_results[node] = {
            'final_signal': float(final_signal[i]),
            'signal_history': [float(flow_history[t][i]) for t in range(len(flow_history))],
        }

    return {
        'node_flows': flow_results,
        'method': 'diffusion',
        'steps': steps,
        'source_nodes': source_nodes,
        'diffusion_matrix_shape': L_norm.shape,
    }


def information_community_detection(
    graph: Any,
    method: str = "infomap",
    **kwargs: Any
) -> Dict[str, Any]:
    """Detect communities using information-theoretic approaches.

    Args:
        graph: NetworkX graph
        method: Community detection method ('infomap', 'map_equation')
        **kwargs: Additional parameters for detection algorithm

    Returns:
        Dictionary with community detection results

    Raises:
        ImportError: If required packages not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for community detection")

    try:
        import community as community_louvain
        HAS_COMMUNITY = True
    except ImportError:
        HAS_COMMUNITY = False

    if method == "infomap":
        return _infomap_community_detection(graph, **kwargs)
    elif method == "map_equation":
        return _map_equation_community_detection(graph, **kwargs)
    elif method == "louvain" and HAS_COMMUNITY:
        return _louvain_community_detection(graph, **kwargs)
    else:
        raise ValueError(f"Unknown or unavailable community detection method: {method}")


def _infomap_community_detection(graph: Any, **kwargs: Any) -> Dict[str, Any]:
    """InfoMap community detection algorithm."""
    # InfoMap implementation would require external library
    # For now, return placeholder
    logger.warning("InfoMap implementation requires external library")

    return {
        'method': 'infomap',
        'status': 'not_implemented',
        'message': 'InfoMap requires infomap package',
        'communities': {},
    }


def _map_equation_community_detection(graph: Any, **kwargs: Any) -> Dict[str, Any]:
    """Map equation community detection."""
    # Map equation implementation would require external library
    logger.warning("Map equation implementation requires external library")

    return {
        'method': 'map_equation',
        'status': 'not_implemented',
        'message': 'Map equation requires map equation package',
        'communities': {},
    }


def _louvain_community_detection(graph: Any, **kwargs: Any) -> Dict[str, Any]:
    """Louvain community detection with information-theoretic interpretation."""
    try:
        import community as community_louvain

        # Convert to undirected graph if needed
        if graph.is_directed():
            G = graph.to_undirected()
        else:
            G = graph

        # Run Louvain algorithm
        partition = community_louvain.best_partition(G, **kwargs)

        # Organize communities
        communities = {}
        for node, community_id in partition.items():
            if community_id not in communities:
                communities[community_id] = []
            communities[community_id].append(node)

        # Calculate modularity
        modularity = community_louvain.modularity(partition, G)

        # Calculate community information measures
        community_info = {}
        for comm_id, nodes in communities.items():
            subgraph = G.subgraph(nodes)

            # Community size and density
            size = len(nodes)
            edges = subgraph.number_of_edges()
            max_edges = size * (size - 1) / 2 if not G.is_directed() else size * (size - 1)

            community_info[comm_id] = {
                'size': size,
                'edges': edges,
                'density': edges / max_edges if max_edges > 0 else 0,
                'nodes': nodes,
            }

        return {
            'method': 'louvain',
            'status': 'completed',
            'modularity': modularity,
            'n_communities': len(communities),
            'communities': communities,
            'community_info': community_info,
        }

    except ImportError:
        return {
            'method': 'louvain',
            'status': 'failed',
            'error': 'python-louvain package required',
        }


def network_information_centrality(
    graph: Any,
    method: str = "entropy",
    normalized: bool = True
) -> Dict[str, float]:
    """Calculate information-theoretic centrality measures for network nodes.

    Args:
        graph: NetworkX graph
        method: Centrality method ('entropy', 'information_flow')
        normalized: Whether to normalize centrality values

    Returns:
        Dictionary mapping node names to centrality values

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for centrality calculation")

    if isinstance(graph, np.ndarray):
        G = nx.from_numpy_array(graph)
    elif hasattr(graph, 'nodes'):
        G = graph
    else:
        raise ValueError("Graph must be NetworkX graph or numpy adjacency matrix")

    if method == "entropy":
        return _entropy_centrality(G, normalized)
    elif method == "information_flow":
        return _information_flow_centrality(G, normalized)
    else:
        raise ValueError(f"Unknown centrality method: {method}")


def _entropy_centrality(G: Any, normalized: bool = True) -> Dict[str, float]:
    """Calculate entropy-based node centrality."""
    centrality = {}

    for node in G.nodes():
        # Calculate entropy of node's neighborhood
        neighbors = list(G.neighbors(node))
        if not neighbors:
            centrality[node] = 0.0
            continue

        # Degree distribution of neighbors
        neighbor_degrees = [G.degree(n) for n in neighbors]
        from collections import Counter
        degree_counts = Counter(neighbor_degrees)

        # Calculate entropy
        total = len(neighbors)
        entropy = 0.0
        for count in degree_counts.values():
            p = count / total
            entropy -= p * math.log2(p) if p > 0 else 0

        centrality[node] = entropy

    # Normalize if requested
    if normalized and centrality:
        max_cent = max(centrality.values())
        if max_cent > 0:
            centrality = {node: cent / max_cent for node, cent in centrality.items()}

    return centrality


def _information_flow_centrality(G: Any, normalized: bool = True) -> Dict[str, float]:
    """Calculate information flow centrality."""
    # Use betweenness centrality as proxy for information flow
    centrality = nx.betweenness_centrality(G, normalized=normalized)

    # Convert to information-theoretic interpretation
    # Higher betweenness = more information flow through node
    return centrality


def network_motif_information(
    graph: Any,
    motif_size: int = 3,
    n_random: int = 100
) -> Dict[str, Any]:
    """Analyze network motifs using information-theoretic measures.

    Args:
        graph: NetworkX graph
        motif_size: Size of motifs to analyze (3 or 4)
        n_random: Number of random graphs for null model

    Returns:
        Dictionary with motif analysis results

    Raises:
        ImportError: If networkx not available
        ValueError: If motif size not supported
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for motif analysis")

    if motif_size not in [3, 4]:
        raise ValueError("Motif size must be 3 or 4")

    # This is a complex implementation that would require external motif detection libraries
    # For now, return placeholder structure

    logger.warning("Network motif analysis requires external motif detection libraries")

    return {
        'motif_size': motif_size,
        'status': 'not_implemented',
        'message': 'Motif analysis requires specialized libraries (e.g., network-motifs)',
        'n_random': n_random,
        'motifs_found': {},
        'z_scores': {},
    }


def information_graph_distance(
    graph1: Any,
    graph2: Any,
    method: str = "entropy"
) -> float:
    """Calculate information-theoretic distance between two graphs.

    Args:
        graph1: First graph
        graph2: Second graph
        method: Distance method ('entropy', 'jensen_shannon')

    Returns:
        Information distance between graphs

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for graph distance calculation")

    # Convert to NetworkX graphs
    if isinstance(graph1, np.ndarray):
        G1 = nx.from_numpy_array(graph1)
    else:
        G1 = graph1

    if isinstance(graph2, np.ndarray):
        G2 = nx.from_numpy_array(graph2)
    else:
        G2 = graph2

    if method == "entropy":
        # Compare network entropies
        ent1 = network_entropy(G1)
        ent2 = network_entropy(G2)
        return abs(ent1 - ent2)

    elif method == "jensen_shannon":
        # Compare degree distributions using Jensen-Shannon divergence
        degrees1 = [d for n, d in G1.degree()]
        degrees2 = [d for n, d in G2.degree()]

        from collections import Counter
        dist1 = Counter(degrees1)
        dist2 = Counter(degrees2)

        # Create common support
        all_degrees = set(dist1.keys()) | set(dist2.keys())
        total1 = sum(dist1.values())
        total2 = sum(dist2.values())

        probs1 = [dist1.get(d, 0) / total1 for d in all_degrees]
        probs2 = [dist2.get(d, 0) / total2 for d in all_degrees]

        # Calculate Jensen-Shannon divergence
        from metainformant.information import syntactic
        jsd = syntactic.jensen_shannon_divergence(probs1, probs2)

        return jsd

    else:
        raise ValueError(f"Unknown distance method: {method}")






