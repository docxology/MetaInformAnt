"""Core network graph functionality for biological networks."""

from __future__ import annotations

import math
from collections import defaultdict, deque
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np

try:
    from ..core.io import write_delimited, read_delimited
except ImportError:
    from metainformant.core.io import write_delimited, read_delimited


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

    def has_edge(self, node1: str, node2: str) -> bool:
        """Check if edge exists between two nodes.
        
        Args:
            node1: First node identifier
            node2: Second node identifier
            
        Returns:
            True if edge exists, False otherwise
        """
        edge = (node1, node2) if self.directed else tuple(sorted([node1, node2]))
        return edge in self.edges

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
    
    Uses breadth-first search (BFS) for unweighted graphs or Dijkstra-like
    BFS for weighted graphs. Finds minimum path length between every
    pair of nodes. For large networks, this can be computationally expensive.
    
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
        
    Note:
        For very large networks (>1000 nodes), consider using subgraph
        extraction or filtering before computing shortest paths.
    """
    if not network.nodes:
        return {}
    
    distances = {}
    
    # Use priority queue for weighted shortest paths (Dijkstra-like)
    # For efficiency, we use a simple BFS approach with edge weights
    for start_node in network.nodes:
        distances[start_node] = {}
        visited = set()
        # Use list as priority queue (simple implementation)
        # For better performance, could use heapq
        queue = [(0.0, start_node)]  # (distance, node)

        while queue:
            # Sort by distance and pop minimum
            queue.sort(key=lambda x: x[0])
            current_dist, current_node = queue.pop(0)

            if current_node in visited:
                continue

            visited.add(current_node)
            distances[start_node][current_node] = current_dist

            # Add neighbors to queue
            for neighbor in network.get_neighbors(current_node):
                if neighbor not in visited:
                    edge_weight = network.get_edge_weight(current_node, neighbor) or 1.0
                    new_dist = current_dist + edge_weight
                    # Check if we already have a shorter path to this neighbor
                    if neighbor not in [n for _, n in queue] or new_dist < min(
                        [d for d, n in queue if n == neighbor], default=float("inf")
                    ):
                        queue.append((new_dist, neighbor))

        # Set infinite distance for unreachable nodes
        for node in network.nodes:
            if node not in distances[start_node]:
                distances[start_node][node] = float("inf")

    return distances


def _betweenness_centrality(network: BiologicalNetwork, node: str) -> float:
    """Calculate betweenness centrality for a single node.
    
    Simplified implementation using shortest paths. For large networks,
    Brandes' algorithm would be more efficient but this provides reasonable
    results for moderate-sized networks.
    
    Args:
        network: Input biological network
        node: Node to calculate betweenness for
        
    Returns:
        Betweenness centrality (normalized, 0 to 1)
    """
    if node not in network.nodes:
        return 0.0
    
    # This is a simplified implementation
    # For large networks, use specialized algorithms like Brandes' algorithm
    betweenness = 0.0
    nodes = list(network.nodes)
    
    if len(nodes) <= 2:
        return 0.0

    # Get all shortest paths
    all_paths = shortest_paths(network)
    
    for source in nodes:
        if source == node:
            continue
        for target in nodes:
            if target == node or target == source:
                continue

            # Find if shortest path from source to target goes through node
            source_to_target = all_paths.get(source, {}).get(target, float("inf"))
            if source_to_target == float("inf"):
                continue
            
            source_to_node = all_paths.get(source, {}).get(node, float("inf"))
            node_to_target = all_paths.get(node, {}).get(target, float("inf"))
            
            # Check if node is on shortest path
            if source_to_node != float("inf") and node_to_target != float("inf"):
                if abs((source_to_node + node_to_target) - source_to_target) < 1e-9:
                    betweenness += 1.0

    # Normalize
    n = len(nodes)
    if n > 2:
        norm = (n - 1) * (n - 2) / (1 if network.directed else 2)
        betweenness /= norm if norm > 0 else 1.0

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


def export_network(network: BiologicalNetwork, filepath: str, format: str = "json") -> None:
    """Export network to file in various formats.
    
    Supports multiple export formats for network persistence and sharing:
    - JSON: Human-readable format with all attributes
    - CSV: Simple edge list format
    - GraphML: XML-based format for graph tools
    
    Args:
        network: Biological network to export
        filepath: Path to output file
        format: Export format ("json", "csv", "graphml"). Default "json"
        
    Examples:
        >>> network = create_network(["A", "B", "C"])
        >>> network.add_edge("A", "B", weight=0.8)
        >>> export_network(network, "network.json", format="json")
        
    Raises:
        ValueError: If format is not supported
    """
    from pathlib import Path
    from ..core.io import dump_json
    
    filepath = Path(filepath)
    
    if format.lower() == "json":
        # Export as JSON
        data = {
            "directed": network.directed,
            "nodes": list(network.nodes),
            "node_attributes": dict(network.node_attrs),
            "edges": [
                {
                    "source": edge[0],
                    "target": edge[1],
                    "weight": weight
                }
                for edge, weight in network.edges.items()
            ]
        }
        dump_json(data, filepath)
        
    elif format.lower() == "csv":
        # Export as CSV edge list
        rows = []
        for edge, weight in network.edges.items():
            rows.append({
                "source": edge[0],
                "target": edge[1],
                "weight": weight
            })
        write_delimited(rows, filepath, delimiter=",")
                
    elif format.lower() == "graphml":
        # Export as GraphML (XML format)
        with open(filepath, "w") as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write('<graphml xmlns="http://graphml.graphdrawing.org/xmlns">\n')
            f.write('  <graph id="G" edgedefault="')
            f.write("directed" if network.directed else "undirected")
            f.write('">\n')
            
            # Write nodes
            for node in network.nodes:
                f.write(f'    <node id="{node}"/>\n')
            
            # Write edges
            for edge, weight in network.edges.items():
                f.write(f'    <edge source="{edge[0]}" target="{edge[1]}">\n')
                f.write(f'      <data key="weight">{weight}</data>\n')
                f.write('    </edge>\n')
            
            f.write('  </graph>\n')
            f.write('</graphml>\n')
            
    else:
        raise ValueError(f"Unsupported export format: {format}. Supported formats: json, csv, graphml")


def import_network(filepath: str, format: str = "json") -> BiologicalNetwork:
    """Import network from file in various formats.
    
    Loads network data from files exported by export_network() or
    compatible formats.
    
    Args:
        filepath: Path to input file
        format: Import format ("json", "csv", "graphml"). Default "json"
        
    Returns:
        BiologicalNetwork object with loaded data
        
    Examples:
        >>> network = import_network("network.json", format="json")
        >>> network.num_nodes()
        10
        
    Raises:
        ValueError: If format is not supported
        FileNotFoundError: If file doesn't exist
    """
    from pathlib import Path
    from ..core.io import load_json
    
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"Network file not found: {filepath}")
    
    if format.lower() == "json":
        data = load_json(filepath)
        network = BiologicalNetwork(directed=data.get("directed", False))
        
        # Add nodes
        for node in data.get("nodes", []):
            attrs = data.get("node_attributes", {}).get(node, {})
            network.add_node(node, **attrs)
        
        # Add edges
        for edge_data in data.get("edges", []):
            network.add_edge(
                edge_data["source"],
                edge_data["target"],
                weight=edge_data.get("weight", 1.0)
            )
        
        return network
        
    elif format.lower() == "csv":
        network = BiologicalNetwork(directed=False)
        
        for row in read_delimited(filepath, delimiter=","):
            source = row["source"]
            target = row["target"]
            weight = float(row.get("weight", 1.0))
            network.add_edge(source, target, weight=weight)
        
        return network
        
    elif format.lower() == "graphml":
        # Basic GraphML parser (simplified)
        import xml.etree.ElementTree as ET
        
        tree = ET.parse(filepath)
        root = tree.getroot()
        
        # Find graph element
        graph = root.find(".//{http://graphml.graphdrawing.org/xmlns}graph")
        if graph is None:
            raise ValueError("Invalid GraphML file: no graph element found")
        
        directed = graph.get("edgedefault") == "directed"
        network = BiologicalNetwork(directed=directed)
        
        # Parse nodes
        for node in graph.findall(".//{http://graphml.graphdrawing.org/xmlns}node"):
            node_id = node.get("id")
            network.add_node(node_id)
        
        # Parse edges
        for edge in graph.findall(".//{http://graphml.graphdrawing.org/xmlns}edge"):
            source = edge.get("source")
            target = edge.get("target")
            weight = 1.0
            
            # Try to get weight from data
            weight_elem = edge.find(".//{http://graphml.graphdrawing.org/xmlns}data[@key='weight']")
            if weight_elem is not None:
                weight = float(weight_elem.text)
            
            network.add_edge(source, target, weight=weight)
        
        return network
        
    else:
        raise ValueError(f"Unsupported import format: {format}. Supported formats: json, csv, graphml")


def network_similarity(network1: BiologicalNetwork, network2: BiologicalNetwork) -> Dict[str, float]:
    """Calculate similarity metrics between two networks.
    
    Computes various similarity measures including Jaccard similarity
    of nodes/edges, edge overlap, and structural similarity.
    
    Args:
        network1: First biological network
        network2: Second biological network
        
    Returns:
        Dictionary containing similarity metrics:
        - node_jaccard: Jaccard similarity of node sets
        - edge_jaccard: Jaccard similarity of edge sets
        - edge_overlap: Fraction of edges in common
        - node_overlap: Fraction of nodes in common
        - structural_similarity: Similarity of edge weights (for common edges)
        
    Examples:
        >>> net1 = create_network(["A", "B", "C"])
        >>> net1.add_edge("A", "B")
        >>> net2 = create_network(["A", "B", "D"])
        >>> net2.add_edge("A", "B")
        >>> sim = network_similarity(net1, net2)
        >>> sim["node_jaccard"]
        0.5
    """
    # Node similarity
    nodes1 = set(network1.nodes)
    nodes2 = set(network2.nodes)
    
    node_intersection = nodes1.intersection(nodes2)
    node_union = nodes1.union(nodes2)
    
    node_jaccard = len(node_intersection) / len(node_union) if node_union else 0.0
    node_overlap = len(node_intersection) / len(nodes1) if nodes1 else 0.0
    
    # Edge similarity
    edges1 = set(network1.edges.keys())
    edges2 = set(network2.edges.keys())
    
    edge_intersection = edges1.intersection(edges2)
    edge_union = edges1.union(edges2)
    
    edge_jaccard = len(edge_intersection) / len(edge_union) if edge_union else 0.0
    edge_overlap = len(edge_intersection) / len(edges1) if edges1 else 0.0
    
    # Structural similarity (weight similarity for common edges)
    structural_similarity = 0.0
    if edge_intersection:
        weight_diffs = []
        for edge in edge_intersection:
            w1 = network1.get_edge_weight(edge[0], edge[1]) or 0.0
            w2 = network2.get_edge_weight(edge[0], edge[1]) or 0.0
            weight_diffs.append(abs(w1 - w2))
        
        # Cosine similarity of weights
        weights1 = [network1.get_edge_weight(e[0], e[1]) or 0.0 for e in edge_intersection]
        weights2 = [network2.get_edge_weight(e[0], e[1]) or 0.0 for e in edge_intersection]
        
        dot_product = sum(w1 * w2 for w1, w2 in zip(weights1, weights2))
        norm1 = sum(w * w for w in weights1) ** 0.5
        norm2 = sum(w * w for w in weights2) ** 0.5
        
        if norm1 > 0 and norm2 > 0:
            structural_similarity = dot_product / (norm1 * norm2)
    
    return {
        "node_jaccard": node_jaccard,
        "edge_jaccard": edge_jaccard,
        "edge_overlap": edge_overlap,
        "node_overlap": node_overlap,
        "structural_similarity": structural_similarity,
    }


def extract_subgraph(network: BiologicalNetwork, nodes: List[str]) -> BiologicalNetwork:
    """Extract subgraph containing only specified nodes and their edges.
    
    Creates a new network containing only the specified nodes and all
    edges between them from the original network.
    
    Args:
        network: Input biological network
        nodes: List of node identifiers to include in subgraph
        
    Returns:
        New BiologicalNetwork object containing subgraph
        
    Examples:
        >>> network = create_network(["A", "B", "C", "D"])
        >>> network.add_edge("A", "B")
        >>> network.add_edge("B", "C")
        >>> subgraph = extract_subgraph(network, ["A", "B"])
        >>> subgraph.num_nodes()
        2
        >>> subgraph.num_edges()
        1
    """
    node_set = set(nodes)
    subgraph = BiologicalNetwork(directed=network.directed)
    
    # Add nodes
    for node in nodes:
        if node in network.nodes:
            attrs = network.node_attrs.get(node, {})
            subgraph.add_node(node, **attrs)
    
    # Add edges between subgraph nodes
    for edge, weight in network.edges.items():
        node1, node2 = edge
        if node1 in node_set and node2 in node_set:
            subgraph.add_edge(node1, node2, weight)
    
    return subgraph


def filter_network(
    network: BiologicalNetwork,
    min_degree: Optional[int] = None,
    max_degree: Optional[int] = None,
    min_edge_weight: Optional[float] = None,
    max_edge_weight: Optional[float] = None,
    component_size: Optional[int] = None,
) -> BiologicalNetwork:
    """Filter network by various criteria.
    
    Creates a new network containing only nodes/edges that meet the
    specified criteria. Useful for removing low-degree nodes, weak edges,
    or isolated components.
    
    Args:
        network: Input biological network
        min_degree: Minimum node degree to include (None = no limit)
        max_degree: Maximum node degree to include (None = no limit)
        min_edge_weight: Minimum edge weight to include (None = no limit)
        max_edge_weight: Maximum edge weight to include (None = no limit)
        component_size: Minimum component size to keep (None = keep all)
        
    Returns:
        New filtered BiologicalNetwork object
        
    Examples:
        >>> network = create_network(["A", "B", "C", "D"])
        >>> network.add_edge("A", "B", weight=0.9)
        >>> network.add_edge("C", "D", weight=0.3)
        >>> filtered = filter_network(network, min_edge_weight=0.5)
        >>> filtered.num_edges()
        1
    """
    filtered = BiologicalNetwork(directed=network.directed)
    
    # Filter nodes by degree
    nodes_to_keep = set()
    for node in network.nodes:
        degree = len(network.get_neighbors(node))
        
        if min_degree is not None and degree < min_degree:
            continue
        if max_degree is not None and degree > max_degree:
            continue
        
        nodes_to_keep.add(node)
    
    # Filter edges by weight and node membership
    for edge, weight in network.edges.items():
        node1, node2 = edge
        
        if node1 not in nodes_to_keep or node2 not in nodes_to_keep:
            continue
        
        if min_edge_weight is not None and weight < min_edge_weight:
            continue
        if max_edge_weight is not None and weight > max_edge_weight:
            continue
        
        filtered.add_edge(node1, node2, weight)
    
    # Filter by component size if specified
    if component_size is not None:
        components = get_connected_components(filtered)
        nodes_in_large_components = set()
        
        for component in components:
            if len(component) >= component_size:
                nodes_in_large_components.update(component)
        
        # Rebuild with only large components
        final_network = BiologicalNetwork(directed=filtered.directed)
        for node in nodes_in_large_components:
            attrs = filtered.node_attrs.get(node, {})
            final_network.add_node(node, **attrs)
        
        for edge, weight in filtered.edges.items():
            if edge[0] in nodes_in_large_components and edge[1] in nodes_in_large_components:
                final_network.add_edge(edge[0], edge[1], weight)
        
        return final_network
    
    return filtered


def get_connected_components(network: BiologicalNetwork) -> List[Set[str]]:
    """Find all connected components in the network.
    
    Uses breadth-first search to identify all connected components
    in an undirected network. For directed networks, finds weakly
    connected components.
    
    Args:
        network: Input biological network
        
    Returns:
        List of sets, each containing node IDs in one component
        
    Examples:
        >>> network = create_network(["A", "B", "C", "D"])
        >>> network.add_edge("A", "B")
        >>> network.add_edge("C", "D")
        >>> components = get_connected_components(network)
        >>> len(components)
        2
    """
    visited = set()
    components = []
    
    for node in network.nodes:
        if node in visited:
            continue
        
        # BFS to find component
        component = set()
        queue = deque([node])
        
        while queue:
            current = queue.popleft()
            if current in visited:
                continue
            
            visited.add(current)
            component.add(current)
            
            # Add neighbors
            for neighbor in network.get_neighbors(current):
                if neighbor not in visited:
                    queue.append(neighbor)
        
        if component:
            components.append(component)
    
    return components


def network_union(network1: BiologicalNetwork, network2: BiologicalNetwork) -> BiologicalNetwork:
    """Create union of two networks.
    
    Creates a new network containing all nodes and edges from both
    networks. Edge weights are summed if an edge exists in both networks.
    
    Args:
        network1: First biological network
        network2: Second biological network
        
    Returns:
        New BiologicalNetwork object containing union
        
    Examples:
        >>> net1 = create_network(["A", "B"])
        >>> net1.add_edge("A", "B", weight=0.5)
        >>> net2 = create_network(["A", "B", "C"])
        >>> net2.add_edge("A", "B", weight=0.3)
        >>> net2.add_edge("B", "C", weight=0.8)
        >>> union = network_union(net1, net2)
        >>> union.num_edges()
        2
        >>> union.get_edge_weight("A", "B")
        0.8  # Sum of weights
    """
    if network1.directed != network2.directed:
        raise ValueError("Both networks must have same directed/undirected setting")
    
    union = BiologicalNetwork(directed=network1.directed)
    
    # Add all nodes
    all_nodes = network1.nodes.union(network2.nodes)
    for node in all_nodes:
        attrs1 = network1.node_attrs.get(node, {})
        attrs2 = network2.node_attrs.get(node, {})
        # Merge attributes (network2 takes precedence for conflicts)
        attrs = {**attrs1, **attrs2}
        union.add_node(node, **attrs)
    
    # Add all edges (sum weights if edge exists in both)
    all_edges = {}
    for edge, weight in network1.edges.items():
        all_edges[edge] = weight
    
    for edge, weight in network2.edges.items():
        if edge in all_edges:
            all_edges[edge] += weight
        else:
            all_edges[edge] = weight
    
    for edge, weight in all_edges.items():
        union.add_edge(edge[0], edge[1], weight)
    
    return union


def network_intersection(network1: BiologicalNetwork, network2: BiologicalNetwork) -> BiologicalNetwork:
    """Create intersection of two networks.
    
    Creates a new network containing only nodes and edges present in
    both networks. Edge weights are averaged if an edge exists in both.
    
    Args:
        network1: First biological network
        network2: Second biological network
        
    Returns:
        New BiologicalNetwork object containing intersection
        
    Examples:
        >>> net1 = create_network(["A", "B", "C"])
        >>> net1.add_edge("A", "B", weight=0.5)
        >>> net2 = create_network(["A", "B", "D"])
        >>> net2.add_edge("A", "B", weight=0.3)
        >>> intersection = network_intersection(net1, net2)
        >>> intersection.num_nodes()
        2
        >>> intersection.num_edges()
        1
    """
    if network1.directed != network2.directed:
        raise ValueError("Both networks must have same directed/undirected setting")
    
    intersection = BiologicalNetwork(directed=network1.directed)
    
    # Add only common nodes
    common_nodes = network1.nodes.intersection(network2.nodes)
    for node in common_nodes:
        attrs1 = network1.node_attrs.get(node, {})
        attrs2 = network2.node_attrs.get(node, {})
        # Merge attributes
        attrs = {**attrs1, **attrs2}
        intersection.add_node(node, **attrs)
    
    # Add only common edges
    edges1 = set(network1.edges.keys())
    edges2 = set(network2.edges.keys())
    common_edges = edges1.intersection(edges2)
    
    for edge in common_edges:
        weight1 = network1.get_edge_weight(edge[0], edge[1]) or 0.0
        weight2 = network2.get_edge_weight(edge[0], edge[1]) or 0.0
        avg_weight = (weight1 + weight2) / 2.0
        intersection.add_edge(edge[0], edge[1], avg_weight)
    
    return intersection


def remove_node(network: BiologicalNetwork, node: str) -> None:
    """Remove a node and all its edges from the network.
    
    Args:
        network: Biological network to modify
        node: Node identifier to remove
        
    Examples:
        >>> network = create_network(["A", "B", "C"])
        >>> network.add_edge("A", "B")
        >>> remove_node(network, "B")
        >>> network.num_nodes()
        2
        >>> network.num_edges()
        0
    """
    if node not in network.nodes:
        return
    
    # Remove all edges connected to this node
    edges_to_remove = []
    for edge in network.edges:
        if edge[0] == node or edge[1] == node:
            edges_to_remove.append(edge)
    
    for edge in edges_to_remove:
        del network.edges[edge]
    
    # Remove node
    network.nodes.discard(node)
    if node in network.node_attrs:
        del network.node_attrs[node]


def remove_edge(network: BiologicalNetwork, node1: str, node2: str) -> None:
    """Remove an edge from the network.
    
    Args:
        network: Biological network to modify
        node1: First node of edge
        node2: Second node of edge
        
    Examples:
        >>> network = create_network(["A", "B", "C"])
        >>> network.add_edge("A", "B")
        >>> network.add_edge("B", "C")
        >>> remove_edge(network, "A", "B")
        >>> network.num_edges()
        1
    """
    edge = (node1, node2) if network.directed else tuple(sorted([node1, node2]))
    if edge in network.edges:
        del network.edges[edge]
