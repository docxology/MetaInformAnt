"""Community detection algorithms for biological networks.

This module provides various community detection algorithms specifically
designed for analyzing biological networks, including protein-protein
interaction networks, gene regulatory networks, and other biological graphs.
"""

from __future__ import annotations

import random
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Optional dependencies for community detection
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    logger.warning("networkx not available, community detection disabled")

try:
    import community as community_louvain

    HAS_LOUVAIN = True
except ImportError:
    HAS_LOUVAIN = False
    community_louvain = None

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def louvain_communities(
    graph: Any, resolution: float = 1.0, randomize: bool = True, random_state: Optional[int] = None, **kwargs: Any
) -> List[List[str]]:
    """Detect communities using the Louvain method.

    The Louvain method is a greedy optimization algorithm that attempts to
    optimize the modularity of the network partitioning.

    Args:
        graph: NetworkX graph
        resolution: Resolution parameter (higher = more communities)
        randomize: Whether to randomize the algorithm
        random_state: Random state for reproducibility
        **kwargs: Additional parameters for the algorithm

    Returns:
        List of community lists (each containing node identifiers)

    Raises:
        ImportError: If required packages not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for community detection")
    if not HAS_LOUVAIN:
        from metainformant.core.utils.optional_deps import warn_optional_dependency

        warn_optional_dependency("python-louvain", "Louvain community detection")
        raise ImportError("python-louvain required for Louvain method")

    # Set random seed if provided
    if random_state is not None:
        random.seed(random_state)
        np.random.seed(random_state)

    # Run Louvain algorithm
    partition = community_louvain.best_partition(graph, resolution=resolution, randomize=randomize, **kwargs)

    # Convert partition to community lists
    communities = {}
    for node, community_id in partition.items():
        if community_id not in communities:
            communities[community_id] = []
        communities[community_id].append(node)

    # Sort communities by size (largest first)
    community_lists = list(communities.values())
    community_lists.sort(key=len, reverse=True)

    logger.info(f"Louvain detected {len(community_lists)} communities")
    return community_lists


def leiden_communities(
    graph: Any, resolution: float = 1.0, n_iterations: int = -1, random_state: Optional[int] = None, **kwargs: Any
) -> List[List[str]]:
    """Detect communities using the Leiden method.

    The Leiden method is an improvement over Louvain that guarantees
    that communities are well-connected.

    Args:
        graph: NetworkX graph
        resolution: Resolution parameter
        n_iterations: Number of iterations (-1 for convergence)
        random_state: Random state for reproducibility
        **kwargs: Additional parameters

    Returns:
        List of community lists

    Raises:
        ImportError: If leidenalg package not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for community detection")

    try:
        import igraph as ig
        import leidenalg as la

        HAS_LEIDEN = True
    except ImportError:
        HAS_LEIDEN = False
        logger.warning("leidenalg and python-igraph not available, Leiden method disabled")
        raise ImportError("leidenalg and python-igraph required for Leiden method")

    # Set random seed
    if random_state is not None:
        random.seed(random_state)
        np.random.seed(random_state)

    # Convert NetworkX to igraph
    ig_graph = _nx_to_igraph(graph)

    # Run Leiden algorithm
    partition = la.find_partition(
        ig_graph, la.ModularityVertexPartition, resolution_parameter=resolution, n_iterations=n_iterations, **kwargs
    )

    # Convert back to community lists
    communities = [list(partition.subgraph(i).vs["name"]) for i in range(len(partition))]

    # Sort by size
    communities.sort(key=len, reverse=True)

    logger.info(f"Leiden detected {len(communities)} communities")
    return communities


def _nx_to_igraph(nx_graph: Any) -> Any:
    """Convert NetworkX graph to igraph format."""
    import igraph as ig

    # Get nodes and edges
    nodes = list(nx_graph.nodes())
    edges = list(nx_graph.edges())

    # Create node mapping
    node_map = {node: i for i, node in enumerate(nodes)}

    # Create igraph
    ig_graph = ig.Graph(directed=nx_graph.is_directed())

    # Add vertices with names
    ig_graph.add_vertices(len(nodes))
    ig_graph.vs["name"] = nodes

    # Add edges
    edge_list = [(node_map[u], node_map[v]) for u, v in edges]
    ig_graph.add_edges(edge_list)

    # Add edge weights if present
    if nx.get_edge_attributes(nx_graph, "weight"):
        weights = []
        for u, v in edges:
            weight = nx_graph[u][v].get("weight", 1.0)
            weights.append(weight)
        ig_graph.es["weight"] = weights

    return ig_graph


def greedy_modularity_communities(graph: Any, **kwargs: Any) -> List[List[str]]:
    """Detect communities using greedy modularity optimization.

    This is NetworkX's built-in greedy modularity algorithm.

    Args:
        graph: NetworkX graph
        **kwargs: Additional parameters

    Returns:
        List of community lists

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for community detection")

    try:
        communities = list(nx.algorithms.community.greedy_modularity_communities(graph, **kwargs))
    except AttributeError:
        # Fallback for older NetworkX versions
        communities = list(nx.algorithms.community.modularity_max.greedy_modularity_communities(graph, **kwargs))

    # Convert frozensets to lists and sort by size
    community_lists = [list(comm) for comm in communities]
    community_lists.sort(key=len, reverse=True)

    logger.info(f"Greedy modularity detected {len(community_lists)} communities")
    return community_lists


def girvan_newman_communities(graph: Any, n_communities: Optional[int] = None, **kwargs: Any) -> List[List[str]]:
    """Detect communities using Girvan-Newman algorithm.

    This algorithm progressively removes edges with highest betweenness centrality.

    Args:
        graph: NetworkX graph
        n_communities: Target number of communities (stops when reached)
        **kwargs: Additional parameters

    Returns:
        List of community lists

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for community detection")

    # Make a copy to avoid modifying original
    G = graph.copy()

    # Use NetworkX's generator
    communities_generator = nx.algorithms.community.girvan_newman(G, **kwargs)

    if n_communities is None:
        # Get the first partitioning
        communities = next(communities_generator)
    else:
        # Generate communities until we reach target number
        for communities in communities_generator:
            if len(communities) >= n_communities:
                break

    # Convert to lists and sort by size
    community_lists = [list(comm) for comm in communities]
    community_lists.sort(key=len, reverse=True)

    logger.info(f"Girvan-Newman detected {len(community_lists)} communities")
    return community_lists


def label_propagation_communities(graph: Any, **kwargs: Any) -> List[List[str]]:
    """Detect communities using label propagation algorithm.

    This algorithm assigns labels to nodes and iteratively updates them
    based on the majority label of neighbors.

    Args:
        graph: NetworkX graph
        **kwargs: Additional parameters

    Returns:
        List of community lists

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for community detection")

    communities = list(nx.algorithms.community.label_propagation_communities(graph, **kwargs))

    # Convert frozensets to lists and sort by size
    community_lists = [list(comm) for comm in communities]
    community_lists.sort(key=len, reverse=True)

    logger.info(f"Label propagation detected {len(community_lists)} communities")
    return community_lists


def asyn_lpa_communities(graph: Any, **kwargs: Any) -> List[List[str]]:
    """Detect communities using asynchronous label propagation.

    Args:
        graph: NetworkX graph
        **kwargs: Additional parameters

    Returns:
        List of community lists

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for community detection")

    communities = list(nx.algorithms.community.asyn_lpa_communities(graph, **kwargs))

    # Convert to lists and sort by size
    community_lists = [list(comm) for comm in communities]
    community_lists.sort(key=len, reverse=True)

    logger.info(f"Asynchronous LPA detected {len(community_lists)} communities")
    return community_lists


def fluid_communities(graph: Any, k: int, **kwargs: Any) -> List[List[str]]:
    """Detect communities using the fluid communities algorithm.

    Args:
        graph: NetworkX graph
        k: Number of communities to find
        **kwargs: Additional parameters

    Returns:
        List of community lists

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for community detection")

    communities = list(nx.algorithms.community.asyn_fluidc(graph, k, **kwargs))

    # Convert to lists and sort by size
    community_lists = [list(comm) for comm in communities]
    community_lists.sort(key=len, reverse=True)

    logger.info(f"Fluid communities detected {len(community_lists)} communities")
    return community_lists


def detect_communities(graph: Any, method: str = "louvain", **kwargs: Any) -> List[List[str]]:
    """Unified interface for community detection.

    Args:
        graph: NetworkX graph
        method: Community detection method
        **kwargs: Method-specific parameters

    Returns:
        List of community lists

    Raises:
        ValueError: If method not supported
    """
    method_map = {
        "louvain": louvain_communities,
        "leiden": leiden_communities,
        "greedy": greedy_modularity_communities,
        "girvan_newman": girvan_newman_communities,
        "label_propagation": label_propagation_communities,
        "asyn_lpa": asyn_lpa_communities,
        "fluid": fluid_communities,
    }

    if method not in method_map:
        available_methods = list(method_map.keys())
        raise ValueError(f"Unknown method '{method}'. Available: {available_methods}")

    return method_map[method](graph, **kwargs)


def evaluate_communities(graph: Any, communities: List[List[str]]) -> Dict[str, Any]:
    """Evaluate quality of community partitioning.

    Args:
        graph: NetworkX graph
        communities: List of community lists

    Returns:
        Dictionary with evaluation metrics

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for community evaluation")

    # Convert communities to partition format
    partition = {}
    for i, community in enumerate(communities):
        for node in community:
            partition[node] = i

    # Calculate modularity
    try:
        if HAS_LOUVAIN:
            modularity = community_louvain.modularity(partition, graph)
        else:
            # Fallback using NetworkX
            modularity = nx.algorithms.community.modularity(graph, communities)
    except (nx.NetworkXError, ValueError, ZeroDivisionError):
        modularity = None

    # Calculate community statistics
    community_sizes = [len(comm) for comm in communities]
    n_communities = len(communities)

    evaluation = {
        "n_communities": n_communities,
        "modularity": modularity,
        "community_sizes": {
            "mean": float(np.mean(community_sizes)),
            "std": float(np.std(community_sizes)),
            "min": int(np.min(community_sizes)),
            "max": int(np.max(community_sizes)),
            "sizes": community_sizes,
        },
        "coverage": sum(community_sizes) / len(graph.nodes()) if graph.nodes() else 0,
    }

    # Calculate conductance for each community
    conductances = []
    for community in communities:
        if len(community) < len(graph.nodes()):
            try:
                conductance = nx.algorithms.cuts.conductance(graph, community)
                conductances.append(conductance)
            except (nx.NetworkXError, ZeroDivisionError):
                pass

    if conductances:
        evaluation["conductance"] = {
            "mean": float(np.mean(conductances)),
            "std": float(np.std(conductances)),
            "values": conductances,
        }

    return evaluation


def compare_community_methods(graph: Any, methods: Optional[List[str]] = None, **kwargs: Any) -> Dict[str, Any]:
    """Compare different community detection methods.

    Args:
        graph: NetworkX graph
        methods: List of methods to compare
        **kwargs: Parameters for individual methods

    Returns:
        Dictionary with comparison results
    """
    if methods is None:
        methods = ["louvain", "greedy", "label_propagation"]

    results = {}

    for method in methods:
        try:
            communities = detect_communities(graph, method=method, **kwargs)
            evaluation = evaluate_communities(graph, communities)

            results[method] = {
                "communities": communities,
                "evaluation": evaluation,
                "n_communities": len(communities),
            }

        except Exception as e:
            logger.warning(f"Method {method} failed: {e}")
            results[method] = {"error": str(e)}

    # Find best method by modularity
    valid_results = [
        (method, result.get("evaluation", {}).get("modularity", 0))
        for method, result in results.items()
        if isinstance(result, dict) and "evaluation" in result
    ]

    if valid_results:
        best_method = max(valid_results, key=lambda x: x[1] if x[1] is not None else -1)[0]
        results["best_method"] = best_method

    return results


def hierarchical_communities(
    graph: Any, method: str = "louvain", max_levels: int = 5, **kwargs: Any
) -> List[List[List[str]]]:
    """Detect hierarchical community structure.

    Args:
        graph: NetworkX graph
        method: Base community detection method
        max_levels: Maximum hierarchy levels to detect
        **kwargs: Parameters for community detection

    Returns:
        List of community hierarchies (each level is a list of communities)

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for hierarchical communities")

    hierarchies = []
    current_graph = graph.copy()

    for level in range(max_levels):
        # Detect communities at current level
        communities = detect_communities(current_graph, method=method, **kwargs)

        if len(communities) <= 1:
            # No further subdivision possible
            break

        hierarchies.append(communities)

        # Create meta-graph for next level
        # Contract each community into a single supernode
        meta_graph = nx.Graph()

        # Add supernodes
        for i, community in enumerate(communities):
            meta_graph.add_node(f"community_{i}", members=community)

        # Add edges between supernodes based on inter-community edges
        for i, comm1 in enumerate(communities):
            for j, comm2 in enumerate(communities):
                if i < j:  # Avoid duplicate edges
                    # Count edges between communities
                    inter_edges = 0
                    for u in comm1:
                        for v in comm2:
                            if current_graph.has_edge(u, v):
                                inter_edges += 1

                    if inter_edges > 0:
                        meta_graph.add_edge(f"community_{i}", f"community_{j}", weight=inter_edges)

        current_graph = meta_graph

        if len(current_graph.nodes()) <= 1:
            break

    logger.info(f"Detected {len(hierarchies)} levels of hierarchical communities")
    return hierarchies


def modularity(graph: Any, communities: List[List[str]]) -> float:
    """Calculate modularity of community partitioning.

    Args:
        graph: NetworkX graph
        communities: List of community lists

    Returns:
        Modularity score (-1 to 1)

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for modularity calculation")

    try:
        return nx.algorithms.community.modularity(graph, communities)
    except Exception:
        # Fallback: simple modularity approximation
        # This is a simplified calculation for cases where NetworkX fails
        total_edges = graph.number_of_edges()
        if total_edges == 0:
            return 0.0

        modularity = 0.0
        for community in communities:
            subgraph = graph.subgraph(community)
            community_edges = subgraph.number_of_edges()
            expected_edges = sum(graph.degree(node) for node in community) ** 2 / (4 * total_edges)

            if expected_edges > 0:
                modularity += (community_edges - expected_edges) / total_edges

        return modularity


def community_metrics(graph: Any, communities: List[List[str]]) -> Dict[str, Any]:
    """Calculate comprehensive community quality metrics.

    Args:
        graph: NetworkX graph
        communities: List of community lists

    Returns:
        Dictionary with comprehensive community metrics

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for community metrics")

    metrics = {}

    # Basic community statistics
    community_sizes = [len(comm) for comm in communities]
    metrics["n_communities"] = len(communities)
    metrics["community_sizes"] = {
        "mean": float(np.mean(community_sizes)),
        "std": float(np.std(community_sizes)),
        "min": int(np.min(community_sizes)),
        "max": int(np.max(community_sizes)),
        "sizes": community_sizes,
    }

    # Modularity
    try:
        metrics["modularity"] = modularity(graph, communities)
    except Exception:
        metrics["modularity"] = None

    # Coverage (fraction of nodes in communities)
    total_nodes_in_communities = sum(len(comm) for comm in communities)
    metrics["coverage"] = total_nodes_in_communities / len(graph.nodes()) if graph.nodes() else 0

    # Conductance for each community
    conductances = []
    for community in communities:
        if len(community) < len(graph.nodes()):
            try:
                conductance = nx.algorithms.cuts.conductance(graph, community)
                conductances.append(conductance)
            except (nx.NetworkXError, ZeroDivisionError):
                pass

    if conductances:
        metrics["conductance"] = {
            "mean": float(np.mean(conductances)),
            "std": float(np.std(conductances)),
            "min": float(np.min(conductances)) if conductances else None,
            "max": float(np.max(conductances)) if conductances else None,
            "values": conductances,
        }

    # Internal density for each community
    densities = []
    for community in communities:
        subgraph = graph.subgraph(community)
        if len(community) > 1:
            max_edges = len(community) * (len(community) - 1) / 2
            if graph.is_directed():
                max_edges *= 2
            density = subgraph.number_of_edges() / max_edges if max_edges > 0 else 0
            densities.append(density)

    if densities:
        metrics["internal_density"] = {
            "mean": float(np.mean(densities)),
            "std": float(np.std(densities)),
            "values": densities,
        }

    # Community quality score (composite metric)
    mod_score = metrics["modularity"] if metrics["modularity"] is not None else 0
    cov_score = metrics["coverage"]
    cond_score = 1 - metrics["conductance"]["mean"] if "conductance" in metrics and conductances else 0.5
    dens_score = metrics["internal_density"]["mean"] if "internal_density" in metrics and densities else 0.5

    metrics["quality_score"] = (mod_score + cov_score + cond_score + dens_score) / 4

    return metrics
