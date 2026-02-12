"""Graph core: BiologicalNetwork class and IO/construction utilities.

This module provides the BiologicalNetwork class and functions for creating,
loading, saving, and manipulating biological networks including adding
nodes/edges from DataFrames, adjacency matrix conversion, and network
validation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Optional network analysis dependencies
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    logger.warning("networkx not available, network functionality disabled")

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


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

    def num_nodes(self) -> int:
        """Get the number of nodes in the network."""
        return self.graph.number_of_nodes()

    def number_of_edges(self) -> int:
        """Get the number of edges in the network."""
        return self.graph.number_of_edges()

    def num_edges(self) -> int:
        """Get the number of edges in the network."""
        return self.graph.number_of_edges()

    def has_edge(self, source: str, target: str) -> bool:
        """Check if an edge exists between two nodes."""
        return self.graph.has_edge(source, target)

    def get_edge_weight(self, source: str, target: str, weight_key: str = "weight") -> float:
        """Get the weight of an edge.

        Args:
            source: Source node ID
            target: Target node ID
            weight_key: Key for weight attribute

        Returns:
            Edge weight value, or 1.0 if not found
        """
        if self.graph.has_edge(source, target):
            return self.graph[source][target].get(weight_key, 1.0)
        return 0.0

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
        if HAS_NUMPY:
            summary["degree_stats"] = {
                "mean": float(np.mean(degrees)),
                "std": float(np.std(degrees)),
                "min": int(min(degrees)),
                "max": int(max(degrees)),
                "median": float(np.median(degrees)),
            }
        else:
            sorted_degrees = sorted(degrees)
            n = len(sorted_degrees)
            summary["degree_stats"] = {
                "mean": sum(degrees) / len(degrees),
                "std": 0.0,
                "min": min(degrees),
                "max": max(degrees),
                "median": float(sorted_degrees[n // 2]),
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
