"""Graph algorithms: metrics, similarity, filtering, centrality, paths.

This module provides algorithmic functions that operate on NetworkX graphs
or BiologicalNetwork instances, including network metrics, similarity
measures, subgraph extraction, filtering, connected components, set
operations (union/intersection), centrality measures, and shortest paths.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Optional network analysis dependencies
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False

from metainformant.networks.analysis.graph_core import BiologicalNetwork


def network_metrics(graph: Any) -> Dict[str, Any]:
    """Calculate comprehensive network metrics.

    Args:
        graph: NetworkX graph

    Returns:
        Dictionary of network metrics

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network metrics")

    try:
        metrics = {}

        # Basic properties
        metrics["num_nodes"] = graph.number_of_nodes()
        metrics["num_edges"] = graph.number_of_edges()
        metrics["is_directed"] = graph.is_directed()
        metrics["is_multigraph"] = graph.is_multigraph()

        if metrics["num_nodes"] > 0:
            # Density: n*(n-1) for directed, n*(n-1)/2 for undirected
            n = metrics["num_nodes"]
            if metrics["is_directed"]:
                max_edges = n * (n - 1)
            else:
                max_edges = n * (n - 1) // 2
            metrics["density"] = metrics["num_edges"] / max_edges if max_edges > 0 else 0

            # Degree statistics
            degrees = [d for n, d in graph.degree()]
            metrics["avg_degree"] = sum(degrees) / len(degrees)
            metrics["max_degree"] = max(degrees)
            try:
                metrics["degree_assortativity"] = nx.degree_assortativity_coefficient(graph)
            except Exception:
                metrics["degree_assortativity"] = None

            # Clustering
            if not metrics["is_directed"]:
                metrics["avg_clustering"] = nx.average_clustering(graph)
                try:
                    metrics["transitivity"] = nx.transitivity(graph)
                except (nx.NetworkXError, ZeroDivisionError):
                    metrics["transitivity"] = None

            # Connected components
            if metrics["is_directed"]:
                components = list(nx.weakly_connected_components(graph))
            else:
                components = list(nx.connected_components(graph))

            metrics["num_components"] = len(components)
            if components:
                component_sizes = [len(c) for c in components]
                metrics["largest_component_size"] = max(component_sizes)
                metrics["avg_component_size"] = sum(component_sizes) / len(component_sizes)

        return metrics

    except Exception as e:
        logger.error(f"Failed to calculate network metrics: {e}")
        return {"error": str(e)}


def subgraph(graph: Any, nodes: List[str]) -> Any:
    """Create subgraph from specified nodes.

    Args:
        graph: NetworkX graph
        nodes: List of nodes to include in subgraph

    Returns:
        NetworkX subgraph

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for subgraph operations")

    # Filter to nodes that exist in the graph
    existing_nodes = [n for n in nodes if n in graph.nodes()]
    if len(existing_nodes) != len(nodes):
        missing = set(nodes) - set(existing_nodes)
        logger.warning(f"Some nodes not found in graph: {missing}")

    return graph.subgraph(existing_nodes)


def export_network(graph: Any, filepath: str | Path, format: str = "json") -> None:
    """Export network to file.

    Args:
        graph: NetworkX graph or BiologicalNetwork
        filepath: Output file path
        format: Export format ('json', 'gml', 'graphml', 'edgelist')

    Raises:
        ImportError: If networkx not available
        ValueError: If format not supported
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network export")

    # Handle BiologicalNetwork
    if isinstance(graph, BiologicalNetwork):
        graph = graph.graph

    filepath = Path(filepath)

    if format.lower() == "json":
        import json

        # Convert to node-link format
        data = nx.node_link_data(graph)
        with open(filepath, "w") as f:
            json.dump(data, f, indent=2)

    elif format.lower() == "gml":
        nx.write_gml(graph, filepath)

    elif format.lower() == "graphml":
        nx.write_graphml(graph, filepath)

    elif format.lower() == "edgelist":
        nx.write_edgelist(graph, filepath)

    else:
        raise ValueError(f"Unsupported export format: {format}")

    logger.info(f"Exported network to {filepath} in {format} format")


def import_network(filepath: str | Path, format: str = "json") -> BiologicalNetwork:
    """Import network from file.

    Args:
        filepath: Input file path
        format: Import format ('json', 'gml', 'graphml', 'edgelist')

    Returns:
        BiologicalNetwork instance

    Raises:
        ImportError: If networkx not available
        ValueError: If format not supported
        FileNotFoundError: If file doesn't exist
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network import")

    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Network file not found: {filepath}")

    if format.lower() == "json":
        import json

        with open(filepath, "r") as f:
            data = json.load(f)
        graph = nx.node_link_graph(data)

    elif format.lower() == "gml":
        graph = nx.read_gml(filepath)

    elif format.lower() == "graphml":
        graph = nx.read_graphml(filepath)

    elif format.lower() == "edgelist":
        graph = nx.read_edgelist(filepath)

    else:
        raise ValueError(f"Unsupported import format: {format}")

    # Wrap in BiologicalNetwork
    network = BiologicalNetwork(directed=graph.is_directed())
    network.graph = graph

    logger.info(f"Imported network from {filepath} in {format} format")
    return network


def network_similarity(graph1: Any, graph2: Any, method: str = "jaccard") -> float:
    """Calculate similarity between two networks.

    Args:
        graph1: First network (NetworkX graph or BiologicalNetwork)
        graph2: Second network (NetworkX graph or BiologicalNetwork)
        method: Similarity method ('jaccard', 'dice', 'overlap')

    Returns:
        Similarity score between 0 and 1
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network similarity")

    # Handle BiologicalNetwork
    if isinstance(graph1, BiologicalNetwork):
        graph1 = graph1.graph
    if isinstance(graph2, BiologicalNetwork):
        graph2 = graph2.graph

    nodes1 = set(graph1.nodes())
    nodes2 = set(graph2.nodes())
    edges1 = set(graph1.edges())
    edges2 = set(graph2.edges())

    if method == "jaccard":
        # Jaccard similarity for nodes
        intersection = len(nodes1 & nodes2)
        union = len(nodes1 | nodes2)
        return intersection / union if union > 0 else 0.0

    elif method == "dice":
        # Dice similarity for edges
        intersection = len(edges1 & edges2)
        total = len(edges1) + len(edges2)
        return 2 * intersection / total if total > 0 else 0.0

    elif method == "overlap":
        # Overlap coefficient for nodes
        intersection = len(nodes1 & nodes2)
        min_size = min(len(nodes1), len(nodes2))
        return intersection / min_size if min_size > 0 else 0.0

    else:
        raise ValueError(f"Unsupported similarity method: {method}")


def extract_subgraph(graph: Any, nodes: List[str]) -> BiologicalNetwork:
    """Extract subgraph containing specified nodes.

    Args:
        graph: Source network (NetworkX graph or BiologicalNetwork)
        nodes: List of node IDs to include

    Returns:
        BiologicalNetwork containing subgraph
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for subgraph extraction")

    # Handle BiologicalNetwork
    if isinstance(graph, BiologicalNetwork):
        source_graph = graph.graph
        directed = graph.is_directed()
    else:
        source_graph = graph
        directed = source_graph.is_directed()

    sg = source_graph.subgraph(nodes)

    # Create new BiologicalNetwork
    result = BiologicalNetwork(directed=directed)
    result.graph = sg.copy()

    return result


def filter_network(
    graph: Any, min_degree: int = 0, max_degree: int | None = None, min_weight: float = 0.0
) -> BiologicalNetwork:
    """Filter network based on node/edge properties.

    Args:
        graph: Source network (NetworkX graph or BiologicalNetwork)
        min_degree: Minimum node degree
        max_degree: Maximum node degree (None for no limit)
        min_weight: Minimum edge weight

    Returns:
        Filtered BiologicalNetwork
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network filtering")

    # Handle BiologicalNetwork
    if isinstance(graph, BiologicalNetwork):
        source_graph = graph.graph
        directed = graph.is_directed()
        metadata = graph.metadata
    else:
        source_graph = graph
        directed = source_graph.is_directed()
        metadata = {}

    # Create copy to modify
    filtered_graph = source_graph.copy()

    # Filter nodes by degree
    nodes_to_remove = []
    for node in filtered_graph.nodes():
        degree = filtered_graph.degree(node)
        if degree < min_degree:
            nodes_to_remove.append(node)
        elif max_degree is not None and degree > max_degree:
            nodes_to_remove.append(node)

    filtered_graph.remove_nodes_from(nodes_to_remove)

    # Filter edges by weight
    if min_weight > 0:
        edges_to_remove = []
        for u, v, data in filtered_graph.edges(data=True):
            weight = data.get("weight", 1.0)
            if weight < min_weight:
                edges_to_remove.append((u, v))

        filtered_graph.remove_edges_from(edges_to_remove)

    # Create result network
    result = BiologicalNetwork(directed=directed, metadata=metadata)
    result.graph = filtered_graph

    return result


def get_connected_components(graph: Any) -> List[List[str]]:
    """Get connected components of the network.

    Args:
        graph: Network (NetworkX graph or BiologicalNetwork)

    Returns:
        List of component node lists
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for connected components")

    # Handle BiologicalNetwork
    if isinstance(graph, BiologicalNetwork):
        graph = graph.graph

    if graph.is_directed():
        components = list(nx.weakly_connected_components(graph))
    else:
        components = list(nx.connected_components(graph))

    return [list(comp) for comp in components]


def network_union(graph1: Any, graph2: Any) -> BiologicalNetwork:
    """Create union of two networks.

    Args:
        graph1: First network
        graph2: Second network

    Returns:
        BiologicalNetwork containing union
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network operations")

    # Handle BiologicalNetwork
    if isinstance(graph1, BiologicalNetwork):
        graph1 = graph1.graph
    if isinstance(graph2, BiologicalNetwork):
        graph2 = graph2.graph

    union_graph = nx.compose(graph1, graph2)

    result = BiologicalNetwork(directed=union_graph.is_directed())
    result.graph = union_graph

    return result


def network_intersection(graph1: Any, graph2: Any) -> BiologicalNetwork:
    """Create intersection of two networks.

    Args:
        graph1: First network
        graph2: Second network

    Returns:
        BiologicalNetwork containing intersection
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network operations")

    # Handle BiologicalNetwork
    if isinstance(graph1, BiologicalNetwork):
        graph1 = graph1.graph
    if isinstance(graph2, BiologicalNetwork):
        graph2 = graph2.graph

    # Get common nodes and edges
    common_nodes = set(graph1.nodes()) & set(graph2.nodes())

    if not common_nodes:
        # Return empty network
        result = BiologicalNetwork(directed=graph1.is_directed())
        return result

    # Create intersection graph
    intersection_graph = nx.Graph() if not graph1.is_directed() else nx.DiGraph()

    # Add common nodes with combined attributes
    for node in common_nodes:
        attrs1 = graph1.nodes.get(node, {})
        attrs2 = graph2.nodes.get(node, {})
        combined_attrs = {**attrs1, **attrs2}
        intersection_graph.add_node(node, **combined_attrs)

    # Add common edges
    for u, v in graph1.edges():
        if u in common_nodes and v in common_nodes and graph2.has_edge(u, v):
            attrs1 = graph1.edges[u, v]
            attrs2 = graph2.edges[u, v]
            combined_attrs = {**attrs1, **attrs2}
            intersection_graph.add_edge(u, v, **combined_attrs)

    result = BiologicalNetwork(directed=intersection_graph.is_directed())
    result.graph = intersection_graph

    return result


def centrality_measures(graph: Any) -> Dict[str, Dict[str, float]]:
    """Calculate various centrality measures for the network.

    Args:
        graph: NetworkX graph or BiologicalNetwork

    Returns:
        Dictionary mapping centrality type to node centrality values
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for centrality calculations")

    # Handle BiologicalNetwork
    if isinstance(graph, BiologicalNetwork):
        graph = graph.graph

    if graph.number_of_nodes() == 0:
        return {}

    results = {}

    try:
        # Degree centrality
        results["degree"] = dict(nx.degree_centrality(graph))
    except (nx.NetworkXError, ZeroDivisionError):
        results["degree"] = {}

    try:
        # Betweenness centrality
        results["betweenness"] = dict(nx.betweenness_centrality(graph))
    except (nx.NetworkXError, ZeroDivisionError):
        results["betweenness"] = {}

    try:
        # Closeness centrality
        results["closeness"] = dict(nx.closeness_centrality(graph))
    except (nx.NetworkXError, ZeroDivisionError):
        results["closeness"] = {}

    try:
        # Eigenvector centrality
        results["eigenvector"] = dict(nx.eigenvector_centrality(graph, max_iter=100))
    except (nx.NetworkXError, nx.PowerIterationFailedConvergence, ZeroDivisionError):
        results["eigenvector"] = {}

    return results


def shortest_paths(graph: Any, source: str | None = None, target: str | None = None) -> Dict[str, Dict[str, float]]:
    """Calculate shortest paths in the network.

    Args:
        graph: NetworkX graph or BiologicalNetwork
        source: Source node (if None, compute all pairs)
        target: Target node (if None, compute all pairs)

    Returns:
        Dictionary of shortest path lengths
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for shortest path calculations")

    # Handle BiologicalNetwork
    if isinstance(graph, BiologicalNetwork):
        graph = graph.graph

    if graph.number_of_nodes() == 0:
        return {}

    try:
        if source is not None and target is not None:
            # Single pair shortest path
            path_length = nx.shortest_path_length(graph, source, target)
            return {source: {target: path_length}}
        elif source is not None:
            # Single source shortest paths
            return {source: dict(nx.shortest_path_length(graph, source))}
        else:
            # All pairs shortest paths (expensive for large graphs)
            return dict(nx.shortest_path_length(graph))
    except nx.NetworkXError:
        return {}
