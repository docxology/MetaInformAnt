"""Graph construction and manipulation utilities for METAINFORMANT.

This module provides comprehensive tools for creating, manipulating, and analyzing
biological networks including protein-protein interaction networks, gene regulatory
networks, and other biological graph structures.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)


class BiologicalNetwork:
    """A biological network class wrapping NetworkX functionality.

    This class provides a convenient interface for working with biological networks,
    including protein-protein interactions, gene regulatory networks, and pathways.
    """

    def __init__(self, directed: bool = False, **kwargs):
        """Initialize a biological network.

        Args:
            directed: Whether the network is directed
            **kwargs: Additional arguments passed to NetworkX graph constructor
        """
        if not HAS_NETWORKX:
            raise ImportError("networkx required for BiologicalNetwork")

        if directed:
            self.graph = nx.DiGraph(**kwargs)
        else:
            self.graph = nx.Graph(**kwargs)

        self.metadata = kwargs.get("metadata", {})

    def add_node(self, node_id: str, **attributes):
        """Add a node to the network.

        Args:
            node_id: Unique identifier for the node
            **attributes: Node attributes (e.g., node_type, name, etc.)
        """
        self.graph.add_node(node_id, **attributes)

    def add_edge(self, source: str, target: str, **attributes):
        """Add an edge to the network.

        Args:
            source: Source node ID
            target: Target node ID
            **attributes: Edge attributes (e.g., weight, interaction_type, etc.)
        """
        self.graph.add_edge(source, target, **attributes)

    def remove_node(self, node_id: str):
        """Remove a node from the network.

        Args:
            node_id: Node ID to remove
        """
        self.graph.remove_node(node_id)

    def remove_edge(self, source: str, target: str):
        """Remove an edge from the network.

        Args:
            source: Source node ID
            target: Target node ID
        """
        self.graph.remove_edge(source, target)

    def get_nodes(self):
        """Get all nodes in the network."""
        return list(self.graph.nodes())

    def get_edges(self):
        """Get all edges in the network."""
        return list(self.graph.edges())

    def number_of_nodes(self) -> int:
        """Get the number of nodes in the network."""
        return self.graph.number_of_nodes()

    def number_of_edges(self) -> int:
        """Get the number of edges in the network."""
        return self.graph.number_of_edges()

    def size(self, weight: str | None = None) -> int | float:
        """Return the number of edges or total of all edge weights.

        Args:
            weight: Edge attribute to use as weight. If None, returns number of edges.

        Returns:
            Number of edges or sum of edge weights.
        """
        return self.graph.size(weight=weight)

    def is_directed(self) -> bool:
        """Check if the network is directed."""
        return self.graph.is_directed()

    # Delegate common NetworkX methods to the underlying graph
    def __iter__(self):
        """Iterate over nodes."""
        return iter(self.graph)

    def __contains__(self, node):
        """Check if node is in network."""
        return node in self.graph

    def __len__(self):
        """Return number of nodes."""
        return len(self.graph)

    def nodes(self, data=False):
        """Return a NodeView of the graph."""
        return self.graph.nodes(data=data)

    def edges(self, data=False, default=None):
        """Return an EdgeView of the graph."""
        return self.graph.edges(data=data, default=default)

    def neighbors(self, node):
        """Return neighbors of a node."""
        return self.graph.neighbors(node)

    def degree(self, node=None, weight=None):
        """Return degree of nodes."""
        if node is None:
            return self.graph.degree(weight=weight)
        return self.graph.degree(node, weight=weight)

    def adj(self):
        """Return adjacency object."""
        return self.graph.adj

    def __len__(self) -> int:
        """Get the number of nodes in the network."""
        return self.number_of_nodes()

    def __contains__(self, node_id: str) -> bool:
        """Check if a node is in the network."""
        return node_id in self.graph


# Optional network analysis dependencies
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    logger.warning("networkx not available, network functionality disabled")


def create_network(edges: List[Tuple[str, str]], directed: bool = False, **kwargs: Any) -> Any:
    """Create a network from edge list.

    Args:
        edges: List of (source, target) tuples
        directed: Whether to create a directed graph
        **kwargs: Additional arguments for graph creation

    Returns:
        NetworkX graph object

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network creation")

    if directed:
        G = nx.DiGraph(**kwargs)
    else:
        G = nx.Graph(**kwargs)

    G.add_edges_from(edges)

    logger.info(
        f"Created {'directed' if directed else 'undirected'} network with {len(G.nodes())} nodes and {len(G.edges())} edges"
    )
    return G


def load_network(path: Union[str, Path], format: str = "edgelist", **kwargs: Any) -> Any:
    """Load network from file.

    Args:
        path: Path to network file
        format: File format ('edgelist', 'adjlist', 'gml', 'graphml', 'json')
        **kwargs: Additional arguments for loading

    Returns:
        NetworkX graph object

    Raises:
        ImportError: If networkx not available
        ValueError: If format not supported
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network loading")

    path = Path(path)

    if format == "edgelist":
        G = nx.read_edgelist(path, **kwargs)
    elif format == "adjlist":
        G = nx.read_adjlist(path, **kwargs)
    elif format == "gml":
        G = nx.read_gml(path, **kwargs)
    elif format == "graphml":
        G = nx.read_graphml(path, **kwargs)
    elif format == "json":
        import json

        with open(path, "r") as f:
            data = json.load(f)

        # Assume node-link format
        G = nx.node_link_graph(data, **kwargs)
    else:
        raise ValueError(f"Unsupported format: {format}")

    logger.info(f"Loaded network from {path} with {len(G.nodes())} nodes and {len(G.edges())} edges")
    return G


def save_network(graph: Any, path: Union[str, Path], format: str = "edgelist", **kwargs: Any) -> None:
    """Save network to file.

    Args:
        graph: NetworkX graph object
        path: Path to save file
        format: File format ('edgelist', 'adjlist', 'gml', 'graphml', 'json')
        **kwargs: Additional arguments for saving

    Raises:
        ImportError: If networkx not available
        ValueError: If format not supported
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network saving")

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if format == "edgelist":
        nx.write_edgelist(graph, path, **kwargs)
    elif format == "adjlist":
        nx.write_adjlist(graph, path, **kwargs)
    elif format == "gml":
        nx.write_gml(graph, path, **kwargs)
    elif format == "graphml":
        nx.write_graphml(graph, path, **kwargs)
    elif format == "json":
        import json

        data = nx.node_link_data(graph, **kwargs)
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
    else:
        raise ValueError(f"Unsupported format: {format}")

    logger.info(f"Saved network to {path}")


def add_nodes_from_dataframe(
    graph: Any, df: Any, node_column: str, attribute_columns: Optional[List[str]] = None
) -> None:
    """Add nodes to graph from DataFrame.

    Args:
        graph: NetworkX graph
        df: pandas DataFrame
        node_column: Column containing node identifiers
        attribute_columns: Columns to use as node attributes

    Raises:
        ImportError: If networkx or pandas not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for node addition")

    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas required for DataFrame operations")

    if not isinstance(df, pd.DataFrame):
        raise ValueError("df must be a pandas DataFrame")

    nodes_to_add = []
    for _, row in df.iterrows():
        node_id = row[node_column]
        attrs = {}

        if attribute_columns:
            for col in attribute_columns:
                if col in df.columns:
                    attrs[col] = row[col]

        nodes_to_add.append((node_id, attrs))

    graph.add_nodes_from(nodes_to_add)
    logger.info(f"Added {len(nodes_to_add)} nodes from DataFrame")


def add_edges_from_dataframe(
    graph: Any,
    df: Any,
    source_column: str,
    target_column: str,
    attribute_columns: Optional[List[str]] = None,
    directed: bool = False,
) -> None:
    """Add edges to graph from DataFrame.

    Args:
        graph: NetworkX graph
        df: pandas DataFrame
        source_column: Column containing source node identifiers
        target_column: Column containing target node identifiers
        attribute_columns: Columns to use as edge attributes
        directed: Whether edges are directed

    Raises:
        ImportError: If networkx or pandas not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for edge addition")

    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas required for DataFrame operations")

    if not isinstance(df, pd.DataFrame):
        raise ValueError("df must be a pandas DataFrame")

    edges_to_add = []
    for _, row in df.iterrows():
        source = row[source_column]
        target = row[target_column]
        attrs = {}

        if attribute_columns:
            for col in attribute_columns:
                if col in df.columns:
                    attrs[col] = row[col]

        if attrs:
            edges_to_add.append((source, target, attrs))
        else:
            edges_to_add.append((source, target))

    graph.add_edges_from(edges_to_add)
    logger.info(f"Added {len(edges_to_add)} edges from DataFrame")


def create_subgraph(graph: Any, nodes: List[str], **kwargs: Any) -> Any:
    """Create subgraph from node list.

    Args:
        graph: Original NetworkX graph
        nodes: List of nodes to include in subgraph
        **kwargs: Additional arguments for subgraph creation

    Returns:
        NetworkX subgraph

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for subgraph creation")

    # Filter nodes that exist in the graph
    existing_nodes = [n for n in nodes if n in graph.nodes()]
    if len(existing_nodes) != len(nodes):
        logger.warning(f"{len(nodes) - len(existing_nodes)} nodes not found in graph")

    subgraph = graph.subgraph(existing_nodes, **kwargs).copy()

    logger.info(f"Created subgraph with {len(subgraph.nodes())} nodes and {len(subgraph.edges())} edges")
    return subgraph


def remove_isolated_nodes(graph: Any) -> None:
    """Remove nodes with no edges from graph.

    Args:
        graph: NetworkX graph

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for node removal")

    isolated_nodes = list(nx.isolates(graph))
    if isolated_nodes:
        graph.remove_nodes_from(isolated_nodes)
        logger.info(f"Removed {len(isolated_nodes)} isolated nodes")


def add_node_attributes(graph: Any, attributes: Dict[str, Dict[str, Any]]) -> None:
    """Add attributes to graph nodes.

    Args:
        attributes: Dictionary mapping node IDs to attribute dictionaries

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for attribute addition")

    nx.set_node_attributes(graph, attributes)
    logger.info(f"Added attributes to {len(attributes)} nodes")


def add_edge_attributes(graph: Any, attributes: Dict[Tuple[str, str], Dict[str, Any]]) -> None:
    """Add attributes to graph edges.

    Args:
        attributes: Dictionary mapping (source, target) tuples to attribute dictionaries

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for attribute addition")

    nx.set_edge_attributes(graph, attributes)
    logger.info(f"Added attributes to {len(attributes)} edges")


def get_node_attributes_dataframe(graph: Any) -> Any:
    """Get node attributes as pandas DataFrame.

    Args:
        graph: NetworkX graph

    Returns:
        pandas DataFrame with node attributes

    Raises:
        ImportError: If networkx or pandas not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for DataFrame conversion")

    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas required for DataFrame operations")

    # Get all node attributes
    node_attrs = {}
    for node in graph.nodes():
        node_attrs[node] = dict(graph.nodes[node])

    if not node_attrs:
        return pd.DataFrame()

    df = pd.DataFrame.from_dict(node_attrs, orient="index")
    df.index.name = "node"

    return df


def get_edge_attributes_dataframe(graph: Any) -> Any:
    """Get edge attributes as pandas DataFrame.

    Args:
        graph: NetworkX graph

    Returns:
        pandas DataFrame with edge attributes

    Raises:
        ImportError: If networkx or pandas not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for DataFrame conversion")

    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas required for DataFrame operations")

    # Get all edge attributes
    edge_data = []
    for source, target, attrs in graph.edges(data=True):
        row = {"source": source, "target": target}
        row.update(attrs)
        edge_data.append(row)

    if not edge_data:
        return pd.DataFrame(columns=["source", "target"])

    df = pd.DataFrame(edge_data)
    return df


def convert_to_adjacency_matrix(graph: Any) -> Any:
    """Convert graph to adjacency matrix.

    Args:
        graph: NetworkX graph

    Returns:
        numpy adjacency matrix

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for matrix conversion")

    import numpy as np

    return nx.to_numpy_array(graph)


def convert_from_adjacency_matrix(matrix: Any, node_labels: Optional[List[str]] = None, directed: bool = False) -> Any:
    """Convert adjacency matrix to graph.

    Args:
        matrix: numpy adjacency matrix
        node_labels: Optional node labels (defaults to 0, 1, 2, ...)
        directed: Whether to create directed graph

    Returns:
        NetworkX graph

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for matrix conversion")

    import numpy as np

    if node_labels is None:
        n_nodes = matrix.shape[0]
        node_labels = [str(i) for i in range(n_nodes)]

    if directed:
        G = nx.from_numpy_array(matrix, create_using=nx.DiGraph)
    else:
        G = nx.from_numpy_array(matrix)

    # Relabel nodes
    mapping = {i: label for i, label in enumerate(node_labels)}
    G = nx.relabel_nodes(G, mapping)

    return G


def get_network_summary(graph: Any) -> Dict[str, Any]:
    """Get comprehensive network summary statistics.

    Args:
        graph: NetworkX graph

    Returns:
        Dictionary with network statistics

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network summary")

    summary = {
        "n_nodes": len(graph.nodes()),
        "n_edges": len(graph.edges()),
        "directed": graph.is_directed(),
        "weighted": any("weight" in graph[u][v] for u, v in graph.edges()),
        "multigraph": graph.is_multigraph(),
    }

    if summary["n_nodes"] > 0:
        # Degree statistics
        degrees = [d for n, d in graph.degree()]
        summary["degree_stats"] = {
            "mean": float(np.mean(degrees)),
            "std": float(np.std(degrees)),
            "min": int(np.min(degrees)),
            "max": int(np.max(degrees)),
            "median": int(np.median(degrees)),
        }

        # Connected components
        if not graph.is_directed():
            components = list(nx.connected_components(graph))
            summary["connected_components"] = {
                "n_components": len(components),
                "largest_component_size": max(len(c) for c in components) if components else 0,
                "component_sizes": [len(c) for c in components],
            }

        # Clustering coefficient (for undirected graphs)
        if not graph.is_directed():
            try:
                clustering = nx.average_clustering(graph)
                summary["average_clustering"] = clustering
            except (nx.NetworkXError, ZeroDivisionError):
                summary["average_clustering"] = None

    return summary


def validate_network(graph: Any) -> Tuple[bool, List[str]]:
    """Validate network structure and data.

    Args:
        graph: NetworkX graph

    Returns:
        Tuple of (is_valid, error_messages)

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network validation")

    errors = []

    # Check for isolated nodes
    isolated = list(nx.isolates(graph))
    if isolated:
        errors.append(f"Graph contains {len(isolated)} isolated nodes")

    # Check for self-loops
    self_loops = list(nx.selfloop_edges(graph))
    if self_loops:
        errors.append(f"Graph contains {len(self_loops)} self-loop edges")

    # Check for multiple edges in simple graph
    if not graph.is_multigraph():
        # Count edges between each pair
        edge_counts = {}
        for u, v in graph.edges():
            key = tuple(sorted([u, v]))
            edge_counts[key] = edge_counts.get(key, 0) + 1

        multi_edges = [k for k, v in edge_counts.items() if v > 1]
        if multi_edges:
            errors.append(f"Simple graph contains {len(multi_edges)} multi-edges")

    # Check node and edge attributes
    try:
        # Try to access node attributes
        node_attrs = list(graph.nodes(data=True))
        if node_attrs:
            first_attrs = node_attrs[0][1]
            if not isinstance(first_attrs, dict):
                errors.append("Node attributes are not dictionaries")
    except Exception as e:
        errors.append(f"Error accessing node attributes: {e}")

    try:
        # Try to access edge attributes
        edge_attrs = list(graph.edges(data=True))
        if edge_attrs:
            first_attrs = edge_attrs[0][2]
            if not isinstance(first_attrs, dict):
                errors.append("Edge attributes are not dictionaries")
    except Exception as e:
        errors.append(f"Error accessing edge attributes: {e}")

    is_valid = len(errors) == 0
    return is_valid, errors


def add_edges_from_interactions(
    graph: Any, interactions: List[Tuple[str, str, float]], weight_threshold: float = 0.0
) -> None:
    """Add edges to graph based on interaction data.

    Args:
        graph: NetworkX graph
        interactions: List of (source, target, weight) tuples
        weight_threshold: Minimum interaction weight for edge creation

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for graph operations")

    # Ensure all nodes exist
    all_nodes = set()
    for source, target, weight in interactions:
        all_nodes.add(source)
        all_nodes.add(target)
    graph.add_nodes_from(all_nodes)

    # Add edges with weights above threshold
    for source, target, weight in interactions:
        if abs(weight) >= weight_threshold:
            graph.add_edge(source, target, weight=weight, interaction_type="protein_interaction")

    logger.info(f"Added {len(interactions)} interaction edges (threshold={weight_threshold})")


def add_edges_from_correlation(graph: Any, correlation_matrix: Any, threshold: float = 0.5) -> None:
    """Add edges to graph based on correlation matrix.

    Args:
        graph: NetworkX graph
        correlation_matrix: Correlation matrix (pandas DataFrame or numpy array)
        threshold: Minimum absolute correlation for edge creation

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for graph operations")

    try:
        import pandas as pd

        is_dataframe = isinstance(correlation_matrix, pd.DataFrame)

        if is_dataframe:
            nodes = list(correlation_matrix.index)
            corr_data = correlation_matrix.values
        else:
            # Assume numpy array with nodes as indices
            nodes = [f"Node_{i}" for i in range(len(correlation_matrix))]
            corr_data = correlation_matrix

        # Ensure graph has all nodes
        graph.add_nodes_from(nodes)

        # Add edges based on correlation threshold
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                corr_value = abs(corr_data[i, j])
                if corr_value >= threshold:
                    graph.add_edge(nodes[i], nodes[j], weight=corr_value)

        logger.info(f"Added edges from correlation matrix (threshold={threshold})")

    except Exception as e:
        logger.error(f"Failed to add edges from correlation: {e}")
        raise


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
            # Density
            max_edges = metrics["num_nodes"] * (metrics["num_nodes"] - 1)
            if metrics["is_directed"]:
                max_edges *= 2
            metrics["density"] = metrics["num_edges"] / max_edges if max_edges > 0 else 0

            # Degree statistics
            degrees = [d for n, d in graph.degree()]
            metrics["avg_degree"] = sum(degrees) / len(degrees)
            metrics["max_degree"] = max(degrees)
            metrics["degree_assortativity"] = nx.degree_assortativity_coefficient(graph)

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

    subgraph = source_graph.subgraph(nodes)

    # Create new BiologicalNetwork
    result = BiologicalNetwork(directed=directed)
    result.graph = subgraph.copy()

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
