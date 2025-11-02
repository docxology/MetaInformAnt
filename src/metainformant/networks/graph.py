"""Core network graph functionality for biological networks."""

from __future__ import annotations

import math
from collections import defaultdict, deque
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np


class BiologicalNetwork:
    """Biological network representation with nodes and weighted edges.
    
    Represents biological networks (protein-protein interactions, regulatory networks,
    etc.) as graphs with nodes and weighted edges. Supports both directed and
    undirected networks with node and edge attributes.
    
    Attributes:
        nodes: Set of node identifiers (strings)
        edges: Dictionary mapping (node1, node2) tuples to edge weights
        node_attrs: Dictionary mapping node IDs to attribute dictionaries
        directed: Boolean indicating if network is directed
        
    Examples:
        >>> network = BiologicalNetwork(directed=False)
        >>> network.add_node("Gene1", function="transcription")
        >>> network.add_edge("Gene1", "Gene2", weight=0.8)
        >>> network.num_nodes()
        2
        >>> network.num_edges()
        1
    """

    def __init__(self, directed: bool = False):
        """Initialize empty biological network.

        Args:
            directed: If True, network edges have direction (from -> to).
                If False (default), edges are undirected and (A,B) == (B,A)
        """
        self.nodes: Set[str] = set()
        self.edges: Dict[Tuple[str, str], float] = {}
        self.node_attrs: Dict[str, Dict[str, Any]] = defaultdict(dict)
        self.directed = directed

    def add_node(self, node: str, **attrs) -> None:
        """Add a node to the network with optional attributes.
        
        Args:
            node: Node identifier (string)
            **attrs: Optional keyword arguments for node attributes
                (e.g., function="transcription", compartment="nucleus")
                
        Examples:
            >>> network = BiologicalNetwork()
            >>> network.add_node("Gene1", function="transcription", compartment="nucleus")
            >>> "Gene1" in network.nodes
            True
            >>> network.node_attrs["Gene1"]["function"]
            'transcription'
        """
        self.nodes.add(node)
        if attrs:
            self.node_attrs[node].update(attrs)

    def add_edge(self, node1: str, node2: str, weight: float = 1.0) -> None:
        """Add an edge between two nodes with optional weight.
        
        Automatically adds nodes if they don't exist. For undirected networks,
        the edge (A, B) is equivalent to (B, A). For directed networks,
        edges have direction from node1 to node2.
        
        Args:
            node1: First node identifier
            node2: Second node identifier
            weight: Edge weight (default 1.0). Can represent interaction
                strength, confidence, or other quantitative measure.
                
        Examples:
            >>> network = BiologicalNetwork(directed=False)
            >>> network.add_edge("Gene1", "Gene2", weight=0.8)
            >>> network.get_edge_weight("Gene2", "Gene1")  # Undirected
            0.8
            
            >>> dir_net = BiologicalNetwork(directed=True)
            >>> dir_net.add_edge("TF1", "Gene1", weight=0.9)
            >>> dir_net.get_edge_weight("Gene1", "TF1")  # Reverse direction
            None
        """
        self.nodes.add(node1)
        self.nodes.add(node2)

        edge = (node1, node2) if self.directed else tuple(sorted([node1, node2]))
        self.edges[edge] = weight

    def get_neighbors(self, node: str) -> List[str]:
        """Get all neighbors of a node.
        
        For undirected networks, returns all nodes connected by an edge.
        For directed networks, returns only outgoing neighbors (nodes
        reached by edges from this node).
        
        Args:
            node: Node identifier
            
        Returns:
            List of neighbor node identifiers. Empty list if node has no neighbors
            or doesn't exist in the network.
            
        Examples:
            >>> network = BiologicalNetwork(directed=False)
            >>> network.add_edge("A", "B"); network.add_edge("A", "C")
            >>> neighbors = network.get_neighbors("A")
            >>> set(neighbors)
            {'B', 'C'}
        """
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
        """Get weight of edge between two nodes.
        
        Args:
            node1: First node identifier
            node2: Second node identifier
            
        Returns:
            Edge weight if edge exists, None otherwise. For undirected networks,
            returns weight regardless of order (get_edge_weight(A, B) == get_edge_weight(B, A)).
            
        Examples:
            >>> network = BiologicalNetwork(directed=False)
            >>> network.add_edge("A", "B", weight=0.75)
            >>> network.get_edge_weight("A", "B")
            0.75
            >>> network.get_edge_weight("B", "A")  # Undirected
            0.75
        """
        edge = (node1, node2) if self.directed else tuple(sorted([node1, node2]))
        return self.edges.get(edge)

    def num_nodes(self) -> int:
        """Get the number of nodes in the network.
        
        Returns:
            Integer count of unique nodes in the network.
            
        Examples:
            >>> network = BiologicalNetwork()
            >>> network.add_node("A"); network.add_node("B")
            >>> network.num_nodes()
            2
        """
        return len(self.nodes)

    def num_edges(self) -> int:
        """Get the number of edges in the network.
        
        Returns:
            Integer count of unique edges. For undirected networks,
            each edge is counted once regardless of direction.
            
        Examples:
            >>> network = BiologicalNetwork(directed=False)
            >>> network.add_edge("A", "B"); network.add_edge("B", "C")
            >>> network.num_edges()
            2
        """
        return len(self.edges)

    def density(self) -> float:
        """Calculate network density (fraction of possible edges present).
        
        Density measures how connected a network is, ranging from 0 (no edges)
        to 1 (fully connected). For directed networks, maximum edges is n(n-1).
        For undirected networks, maximum edges is n(n-1)/2.
        
        Returns:
            Density value in [0, 1]. Returns 0.0 if network has fewer than 2 nodes.
            Formula: density = actual_edges / max_possible_edges
            
        Examples:
            >>> network = BiologicalNetwork(directed=False)
            >>> network.add_node("A"); network.add_node("B"); network.add_node("C")
            >>> network.add_edge("A", "B")
            >>> network.density()  # 1 edge out of 3 possible
            0.333...
        """
        n = self.num_nodes()
        if n < 2:
            return 0.0
        max_edges = n * (n - 1) // (1 if self.directed else 2)
        return self.num_edges() / max_edges if max_edges > 0 else 0.0


def create_network(nodes: List[str], directed: bool = False) -> BiologicalNetwork:
    """Create a biological network with specified nodes.
    
    Convenience function to create a network and add nodes in one step.
    
    Args:
        nodes: List of node identifier strings
        directed: If True, network is directed. If False (default), undirected.
        
    Returns:
        BiologicalNetwork object with specified nodes but no edges
        
    Examples:
        >>> network = create_network(["A", "B", "C"], directed=False)
        >>> network.num_nodes()
        3
    """
    network = BiologicalNetwork(directed=directed)
    for node in nodes:
        network.add_node(node)
    return network


def add_edges_from_correlation(
    network: BiologicalNetwork, correlation_matrix: np.ndarray, node_names: List[str], threshold: float = 0.7
) -> None:
    """Add edges to network based on correlation matrix.
    
    Constructs network edges from pairwise correlations between nodes.
    Only correlations above the threshold are added as edges, with edge
    weight equal to the absolute correlation value.
    
    Args:
        network: BiologicalNetwork object to add edges to (nodes should exist)
        correlation_matrix: Square NumPy array of shape (n, n) containing
            pairwise correlations. Must be symmetric for undirected networks.
        node_names: List of node identifiers corresponding to matrix
            rows/columns. Length must equal matrix dimension.
        threshold: Minimum absolute correlation value (0-1) to create an edge.
            Default 0.7. Only |correlation| >= threshold creates edges.
            
    Raises:
        ValueError: If correlation matrix dimensions don't match node_names length
        
    Examples:
        >>> network = create_network(["GENE1", "GENE2", "GENE3"], directed=False)
        >>> corr_matrix = np.array([[1.0, 0.8, 0.3], [0.8, 1.0, 0.2], [0.3, 0.2, 1.0]])
        >>> add_edges_from_correlation(network, corr_matrix, ["GENE1", "GENE2", "GENE3"], threshold=0.7)
        >>> network.num_edges()
        1  # Only GENE1-GENE2 above threshold
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
    """Add edges to network from a list of interactions.
    
    Convenience function to add multiple edges at once from a structured
    interaction list. Automatically adds nodes if they don't exist.
    
    Args:
        network: BiologicalNetwork object to add edges to
        interactions: List of tuples, each containing (node1, node2, weight).
            Node identifiers are strings, weight is float.
            
    Examples:
        >>> network = create_network(["A", "B", "C"], directed=False)
        >>> interactions = [("A", "B", 0.8), ("B", "C", 0.6), ("A", "C", 0.9)]
        >>> add_edges_from_interactions(network, interactions)
        >>> network.num_edges()
        3
    """
    for node1, node2, weight in interactions:
        network.add_edge(node1, node2, weight)


def network_metrics(network: BiologicalNetwork) -> Dict[str, float]:
    """Calculate basic network topology metrics.
    
    Computes fundamental network statistics including size, connectivity,
    and degree distributions.
    
    Args:
        network: Input biological network
        
    Returns:
        Dictionary with keys:
        - num_nodes: Total number of nodes
        - num_edges: Total number of edges
        - density: Edge density (fraction of possible edges present)
        - avg_degree: Average node degree
        - max_degree: Maximum node degree
        - min_degree: Minimum node degree
        
    Examples:
        >>> network = create_network(["A", "B", "C"])
        >>> network.add_edge("A", "B")
        >>> metrics = network_metrics(network)
        >>> metrics["num_nodes"]
        3
        >>> metrics["density"]
        0.333...
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
    """Calculate multiple centrality measures for all nodes.
    
    Computes degree, closeness, betweenness, and eigenvector centralities
    for each node in the network. These measures identify important/hub nodes.
    
    Args:
        network: Input biological network
        
    Returns:
        Dictionary mapping centrality type to nested dictionary of
        node -> centrality_value. Centrality types:
        - "degree": Normalized degree centrality (0 to 1)
        - "closeness": Inverse of average shortest path distance
        - "betweenness": Fraction of shortest paths passing through node
        - "eigenvector": Importance based on connections to important nodes
        
    Examples:
        >>> network = create_network(["A", "B", "C", "D"])
        >>> network.add_edge("A", "B"); network.add_edge("A", "C")
        >>> network.add_edge("B", "D")
        >>> centralities = centrality_measures(network)
        >>> centralities["degree"]["A"]
        0.666...
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
    """Calculate shortest path distances between all pairs of nodes.
    
    Uses breadth-first search (BFS) for unweighted graphs or weighted
    BFS for weighted graphs. Finds minimum path length between every
    pair of nodes.
    
    Args:
        network: Input biological network
        
    Returns:
        Nested dictionary mapping source_node -> target_node -> distance.
        Distance is measured in edge weights (1.0 for unweighted edges).
        Unreachable nodes have distance of float("inf").
        Distance to self is always 0.0.
        
    Examples:
        >>> network = create_network(["A", "B", "C", "D"], directed=False)
        >>> network.add_edge("A", "B", weight=1.0)
        >>> network.add_edge("B", "C", weight=2.0)
        >>> paths = shortest_paths(network)
        >>> paths["A"]["C"]
        3.0  # Path A -> B -> C
        >>> paths["A"]["D"]
        inf  # Not reachable
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
