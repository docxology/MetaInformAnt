"""Core network graph functionality for biological networks."""

from __future__ import annotations

import math
from collections import defaultdict, deque
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np


class BiologicalNetwork:
    """Biological network representation with nodes and weighted edges."""

    def __init__(self, directed: bool = False):
        """Initialize empty biological network.

        Args:
            directed: Whether the network is directed
        """
        self.nodes: Set[str] = set()
        self.edges: Dict[Tuple[str, str], float] = {}
        self.node_attrs: Dict[str, Dict[str, Any]] = defaultdict(dict)
        self.directed = directed

    def add_node(self, node: str, **attrs) -> None:
        """Add a node with optional attributes."""
        self.nodes.add(node)
        if attrs:
            self.node_attrs[node].update(attrs)

    def add_edge(self, node1: str, node2: str, weight: float = 1.0) -> None:
        """Add an edge between two nodes with optional weight."""
        self.nodes.add(node1)
        self.nodes.add(node2)

        edge = (node1, node2) if self.directed else tuple(sorted([node1, node2]))
        self.edges[edge] = weight

    def get_neighbors(self, node: str) -> List[str]:
        """Get neighbors of a node."""
        neighbors = []
        for edge in self.edges:
            if self.directed:
                if edge[0] == node:
                    neighbors.append(edge[1])
            else:
                if edge[0] == node:
                    neighbors.append(edge[1])
                elif edge[1] == node:
                    neighbors.append(edge[0])
        return neighbors

    def get_edge_weight(self, node1: str, node2: str) -> Optional[float]:
        """Get weight of edge between two nodes."""
        edge = (node1, node2) if self.directed else tuple(sorted([node1, node2]))
        return self.edges.get(edge)

    def num_nodes(self) -> int:
        """Number of nodes in network."""
        return len(self.nodes)

    def num_edges(self) -> int:
        """Number of edges in network."""
        return len(self.edges)

    def density(self) -> float:
        """Network density (fraction of possible edges present)."""
        n = self.num_nodes()
        if n < 2:
            return 0.0
        max_edges = n * (n - 1) // (1 if self.directed else 2)
        return self.num_edges() / max_edges if max_edges > 0 else 0.0


def create_network(nodes: List[str], directed: bool = False) -> BiologicalNetwork:
    """Create a biological network with given nodes.

    Args:
        nodes: List of node identifiers
        directed: Whether network is directed

    Returns:
        Empty network with specified nodes
    """
    network = BiologicalNetwork(directed=directed)
    for node in nodes:
        network.add_node(node)
    return network


def add_edges_from_correlation(
    network: BiologicalNetwork, correlation_matrix: np.ndarray, node_names: List[str], threshold: float = 0.7
) -> None:
    """Add edges to network based on correlation matrix.

    Args:
        network: Network to add edges to
        correlation_matrix: Square correlation matrix
        node_names: Names corresponding to matrix rows/columns
        threshold: Minimum absolute correlation for edge creation
    """
    n = len(node_names)
    if correlation_matrix.shape != (n, n):
        raise ValueError("Correlation matrix dimensions don't match node names")

    for i in range(n):
        for j in range(i + 1 if not network.directed else 0, n):
            if i != j:
                corr = correlation_matrix[i, j]
                if abs(corr) >= threshold:
                    network.add_edge(node_names[i], node_names[j], weight=abs(corr))


def add_edges_from_interactions(network: BiologicalNetwork, interactions: List[Tuple[str, str, float]]) -> None:
    """Add edges from interaction list.

    Args:
        network: Network to add edges to
        interactions: List of (node1, node2, weight) tuples
    """
    for node1, node2, weight in interactions:
        network.add_edge(node1, node2, weight)


def network_metrics(network: BiologicalNetwork) -> Dict[str, float]:
    """Calculate basic network metrics.

    Args:
        network: Input network

    Returns:
        Dictionary of network metrics
    """
    metrics = {
        "num_nodes": network.num_nodes(),
        "num_edges": network.num_edges(),
        "density": network.density(),
    }

    # Average degree
    degrees = []
    for node in network.nodes:
        degrees.append(len(network.get_neighbors(node)))

    if degrees:
        metrics["avg_degree"] = sum(degrees) / len(degrees)
        metrics["max_degree"] = max(degrees)
        metrics["min_degree"] = min(degrees)
    else:
        metrics["avg_degree"] = 0.0
        metrics["max_degree"] = 0
        metrics["min_degree"] = 0

    return metrics


def centrality_measures(network: BiologicalNetwork) -> Dict[str, Dict[str, float]]:
    """Calculate centrality measures for all nodes.

    Args:
        network: Input network

    Returns:
        Dictionary mapping centrality type to node centrality values
    """
    centralities = {"degree": {}, "closeness": {}, "betweenness": {}, "eigenvector": {}}

    # Degree centrality
    for node in network.nodes:
        degree = len(network.get_neighbors(node))
        max_degree = len(network.nodes) - 1
        centralities["degree"][node] = degree / max_degree if max_degree > 0 else 0.0

    # Closeness centrality
    distances = shortest_paths(network)
    for node in network.nodes:
        if node in distances:
            path_lengths = [d for d in distances[node].values() if d != float("inf")]
            if path_lengths:
                avg_distance = sum(path_lengths) / len(path_lengths)
                centralities["closeness"][node] = 1.0 / avg_distance if avg_distance > 0 else 0.0
            else:
                centralities["closeness"][node] = 0.0
        else:
            centralities["closeness"][node] = 0.0

    # Simplified betweenness (computationally expensive for large networks)
    for node in network.nodes:
        centralities["betweenness"][node] = _betweenness_centrality(network, node)

    # Simplified eigenvector centrality
    centralities["eigenvector"] = _eigenvector_centrality(network)

    return centralities


def shortest_paths(network: BiologicalNetwork) -> Dict[str, Dict[str, float]]:
    """Calculate shortest paths between all pairs of nodes using BFS.

    Args:
        network: Input network

    Returns:
        Dictionary mapping source -> target -> distance
    """
    distances = {}

    for start_node in network.nodes:
        distances[start_node] = {}
        visited = set()
        queue = deque([(start_node, 0.0)])

        while queue:
            current_node, current_dist = queue.popleft()

            if current_node in visited:
                continue

            visited.add(current_node)
            distances[start_node][current_node] = current_dist

            # Add neighbors to queue
            for neighbor in network.get_neighbors(current_node):
                if neighbor not in visited:
                    edge_weight = network.get_edge_weight(current_node, neighbor) or 1.0
                    queue.append((neighbor, current_dist + edge_weight))

        # Set infinite distance for unreachable nodes
        for node in network.nodes:
            if node not in distances[start_node]:
                distances[start_node][node] = float("inf")

    return distances


def _betweenness_centrality(network: BiologicalNetwork, node: str) -> float:
    """Calculate betweenness centrality for a single node (simplified)."""
    # This is a simplified implementation
    # For large networks, use specialized algorithms like Brandes' algorithm
    betweenness = 0.0
    nodes = list(network.nodes)

    for i, source in enumerate(nodes):
        if source == node:
            continue
        for j, target in enumerate(nodes):
            if j <= i or target == node or target == source:
                continue

            # Find if shortest path from source to target goes through node
            paths_through_node = _paths_through_node(network, source, target, node)
            total_paths = _count_shortest_paths(network, source, target)

            if total_paths > 0:
                betweenness += paths_through_node / total_paths

    # Normalize
    n = len(nodes)
    if n > 2:
        norm = (n - 1) * (n - 2) / (1 if network.directed else 2)
        betweenness /= norm

    return betweenness


def _eigenvector_centrality(network: BiologicalNetwork, max_iter: int = 100) -> Dict[str, float]:
    """Calculate eigenvector centrality using power iteration."""
    nodes = list(network.nodes)
    n = len(nodes)

    if n == 0:
        return {}

    # Initialize centrality values
    centrality = {node: 1.0 / n for node in nodes}

    # Power iteration
    for _ in range(max_iter):
        new_centrality = {node: 0.0 for node in nodes}

        for node in nodes:
            for neighbor in network.get_neighbors(node):
                weight = network.get_edge_weight(node, neighbor) or 1.0
                new_centrality[node] += centrality[neighbor] * weight

        # Normalize
        norm = sum(new_centrality.values())
        if norm > 0:
            for node in nodes:
                new_centrality[node] /= norm

        # Check convergence
        diff = sum(abs(new_centrality[node] - centrality[node]) for node in nodes)
        centrality = new_centrality

        if diff < 1e-6:
            break

    return centrality


def _paths_through_node(network: BiologicalNetwork, source: str, target: str, intermediate: str) -> int:
    """Count shortest paths from source to target that pass through intermediate node."""
    # Simplified implementation - in practice, use more efficient algorithms
    return 1 if _has_path_through(network, source, target, intermediate) else 0


def _count_shortest_paths(network: BiologicalNetwork, source: str, target: str) -> int:
    """Count number of shortest paths between source and target."""
    # Simplified implementation - returns 1 if path exists, 0 otherwise
    distances = shortest_paths(network)
    return 1 if distances[source][target] != float("inf") else 0


def _has_path_through(network: BiologicalNetwork, source: str, target: str, intermediate: str) -> bool:
    """Check if shortest path from source to target goes through intermediate."""
    distances = shortest_paths(network)
    source_to_inter = distances[source][intermediate]
    inter_to_target = distances[intermediate][target]
    source_to_target = distances[source][target]

    if source_to_inter == float("inf") or inter_to_target == float("inf"):
        return False

    return abs((source_to_inter + inter_to_target) - source_to_target) < 1e-9
