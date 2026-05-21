"""Community detection algorithms for biological networks.

This module provides various community detection algorithms specifically
designed for analyzing biological networks, including protein-protein
interaction networks, gene regulatory networks, and other biological graphs.
"""

from __future__ import annotations

import random
from typing import Any, Dict, List, Optional

from metainformant.core.utils import logging

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


def _as_networkx_graph(graph: Any) -> Any:
    """Return a NetworkX graph, accepting the project BiologicalNetwork wrapper."""
    if hasattr(graph, "to_networkx"):
        graph = graph.to_networkx()
    elif hasattr(graph, "graph") and HAS_NETWORKX and isinstance(graph.graph, (nx.Graph, nx.DiGraph)):
        graph = graph.graph

    if HAS_NETWORKX and isinstance(graph, (nx.Graph, nx.DiGraph)):
        graph = graph.copy()
        for _, _, data in graph.edges(data=True):
            weight = data.get("weight")
            if weight is not None and weight < 0:
                data["weight"] = 0.0
    return graph


def _communities_to_lists(communities: Any, graph: Any | None = None) -> list[list[str]]:
    """Normalize community mappings or iterable communities to list-of-lists."""
    valid_nodes = set(graph.nodes()) if graph is not None and hasattr(graph, "nodes") else None

    if isinstance(communities, dict):
        grouped: dict[Any, list[str]] = {}
        for node, community_id in communities.items():
            if valid_nodes is not None and node not in valid_nodes:
                continue
            grouped.setdefault(community_id, []).append(node)
        community_lists = list(grouped.values())
    else:
        community_lists = []
        for community in communities:
            if isinstance(community, str):
                members = [community]
            else:
                members = list(community)
            if valid_nodes is not None:
                members = [node for node in members if node in valid_nodes]
            if members:
                community_lists.append(members)

    community_lists.sort(key=len, reverse=True)
    return community_lists


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

    graph = _as_networkx_graph(graph)
    if graph.number_of_nodes() == 0:
        return []
    if graph.number_of_edges() == 0:
        return [[node] for node in graph.nodes()]

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

    graph = _as_networkx_graph(graph)
    if graph.number_of_nodes() == 0:
        return []
    if graph.number_of_edges() == 0:
        return [[node] for node in graph.nodes()]

    try:
        import igraph as ig  # noqa: F401
        import leidenalg as la

    except ImportError:
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


def _networkx_greedy_kwargs(kwargs: Dict[str, Any]) -> Dict[str, Any]:
    """Keep only kwargs accepted by NetworkX greedy modularity."""
    allowed = {"weight", "resolution", "cutoff", "best_n"}
    return {key: value for key, value in kwargs.items() if key in allowed}


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

    graph = _as_networkx_graph(graph)
    if graph.number_of_nodes() == 0:
        return []
    if graph.number_of_edges() == 0:
        return [[node] for node in graph.nodes()]

    greedy_kwargs = _networkx_greedy_kwargs(kwargs)
    try:
        communities = list(nx.algorithms.community.greedy_modularity_communities(graph, **greedy_kwargs))
    except AttributeError:
        # Fallback for older NetworkX versions
        communities = list(
            nx.algorithms.community.modularity_max.greedy_modularity_communities(graph, **greedy_kwargs)
        )

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

    graph = _as_networkx_graph(graph)

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

    graph = _as_networkx_graph(graph)
    if graph.number_of_nodes() == 0:
        return []
    if graph.number_of_edges() == 0:
        return [[node] for node in graph.nodes()]

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

    graph = _as_networkx_graph(graph)
    if graph.number_of_nodes() == 0:
        return []
    if graph.number_of_edges() == 0:
        return [[node] for node in graph.nodes()]

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

    graph = _as_networkx_graph(graph)
    if graph.number_of_nodes() == 0:
        return []

    communities = list(nx.algorithms.community.asyn_fluidc(graph, k, **kwargs))

    # Convert to lists and sort by size
    community_lists = [list(comm) for comm in communities]
    community_lists.sort(key=len, reverse=True)

    logger.info(f"Fluid communities detected {len(community_lists)} communities")
    return community_lists


def detect_communities(graph: Any, method: str = "louvain", **kwargs: Any) -> Dict[str, int]:
    """Unified interface for community detection.

    Args:
        graph: NetworkX graph
        method: Community detection method
        **kwargs: Method-specific parameters

    Returns:
        Dict mapping node IDs to community IDs

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
        raise ValueError(f"Unknown community detection method '{method}'. Available: {available_methods}")

    # Get community lists from underlying method
    nx_graph = _as_networkx_graph(graph)
    try:
        community_lists = method_map[method](nx_graph, **kwargs)
    except ImportError:
        fallback = louvain_communities if method == "leiden" and HAS_LOUVAIN else greedy_modularity_communities
        community_lists = fallback(nx_graph, **kwargs)
    except Exception as exc:
        if method != "greedy":
            logger.warning(f"Method {method} failed ({exc}); falling back to greedy modularity")
            community_lists = greedy_modularity_communities(nx_graph)
        else:
            raise

    # Convert to dict mapping node -> community ID
    node_to_community = {}
    for comm_id, community in enumerate(community_lists):
        for node in community:
            node_to_community[node] = comm_id

    return node_to_community


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

    graph = _as_networkx_graph(graph)
    communities = _communities_to_lists(communities, graph)

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
    if not community_sizes:
        return {
            "n_communities": 0,
            "num_communities": 0,
            "modularity": 0.0,
            "community_sizes": {"mean": 0.0, "std": 0.0, "min": 0, "max": 0, "sizes": []},
            "coverage": 0.0,
        }
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
) -> Dict[int, Dict[str, int]]:
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

    levels = int(kwargs.pop("levels", max_levels))
    seed = kwargs.pop("seed", None)
    if seed is not None and "random_state" not in kwargs:
        kwargs["random_state"] = seed

    hierarchies: dict[int, dict[str, int]] = {}
    current_graph = _as_networkx_graph(graph)

    previous: dict[str, int] = {}
    for level in range(levels):
        # Detect communities at current level
        partition = detect_communities(current_graph, method=method, **kwargs)
        hierarchies[level] = partition
        communities = _communities_to_lists(partition, current_graph)

        if len(communities) <= 1:
            previous = partition
            continue

        previous = partition

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
            current_graph = _as_networkx_graph(graph)

    while len(hierarchies) < levels:
        hierarchies[len(hierarchies)] = previous.copy()

    logger.info(f"Detected {len(hierarchies)} levels of hierarchical communities")
    return hierarchies


def compare_communities(communities1: Dict[str, int], communities2: Dict[str, int]) -> Dict[str, float]:
    """Compare two node-to-community assignments."""
    common_nodes = sorted(set(communities1) & set(communities2))
    if not common_nodes:
        return {"normalized_mutual_information": 0.0, "adjusted_rand_index": 0.0}

    labels1 = [communities1[node] for node in common_nodes]
    labels2 = [communities2[node] for node in common_nodes]

    try:
        from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

        nmi = float(normalized_mutual_info_score(labels1, labels2))
        ari = float(adjusted_rand_score(labels1, labels2))
    except Exception:
        matches = sum(a == b for a, b in zip(labels1, labels2))
        nmi = matches / len(common_nodes)
        ari = nmi

    return {"normalized_mutual_information": nmi, "adjusted_rand_index": ari}


def community_stability(
    graph: Any,
    method: str = "louvain",
    n_runs: int = 10,
    seed: int | None = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Estimate community stability across repeated runs."""
    nx_graph = _as_networkx_graph(graph)
    partitions: list[dict[str, int]] = []
    modularities: list[float] = []

    for run in range(n_runs):
        run_kwargs = dict(kwargs)
        run_kwargs["random_state"] = (seed or 0) + run
        partition = detect_communities(nx_graph, method=method, **run_kwargs)
        partitions.append(partition)
        modularities.append(modularity(nx_graph, partition))

    comparisons = []
    for i, partition1 in enumerate(partitions):
        for partition2 in partitions[i + 1 :]:
            comparisons.append(compare_communities(partition1, partition2)["normalized_mutual_information"])

    stability_score = float(np.mean(comparisons)) if comparisons else 1.0
    return {
        "stability_score": stability_score,
        "avg_modularity": float(np.mean(modularities)) if modularities else 0.0,
        "modularity_values": modularities,
    }


def optimize_resolution(
    graph: Any,
    resolution_range: tuple[float, float] = (0.5, 2.0),
    n_points: int = 10,
    method: str = "louvain",
    **kwargs: Any,
) -> Dict[str, Any]:
    """Find a community resolution with the best modularity."""
    nx_graph = _as_networkx_graph(graph)
    resolutions = np.linspace(resolution_range[0], resolution_range[1], n_points)
    results = []

    for resolution in resolutions:
        partition = detect_communities(nx_graph, method=method, resolution=float(resolution), **kwargs)
        score = modularity(nx_graph, partition)
        results.append(
            {
                "resolution": float(resolution),
                "modularity": float(score),
                "num_communities": len(set(partition.values())),
            }
        )

    best = (
        max(results, key=lambda item: item["modularity"])
        if results
        else {
            "resolution": resolution_range[0],
            "modularity": 0.0,
        }
    )
    return {
        "optimal_resolution": best["resolution"],
        "optimal_modularity": best["modularity"],
        "results": results,
    }


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

    graph = _as_networkx_graph(graph)
    community_lists = _communities_to_lists(communities, graph)
    if not community_lists or graph.number_of_edges() == 0:
        return 0.0

    try:
        return float(nx.algorithms.community.modularity(graph, community_lists))
    except Exception:
        # Fallback: simple modularity approximation
        # This is a simplified calculation for cases where NetworkX fails
        total_edges = graph.number_of_edges()
        if total_edges == 0:
            return 0.0

        modularity = 0.0
        for community in community_lists:
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

    graph = _as_networkx_graph(graph)
    community_lists = _communities_to_lists(communities, graph)
    metrics = {}

    # Basic community statistics
    community_sizes = [len(comm) for comm in community_lists]
    if not community_sizes:
        return {
            "n_communities": 0,
            "num_communities": 0,
            "community_sizes": {"mean": 0.0, "std": 0.0, "min": 0, "max": 0, "sizes": []},
            "avg_community_size": 0.0,
            "modularity": 0.0,
            "coverage": 0.0,
            "internal_edge_ratio": 0.0,
            "quality_score": 0.0,
        }

    metrics["n_communities"] = len(community_lists)
    metrics["num_communities"] = len(community_lists)
    metrics["community_sizes"] = {
        "mean": float(np.mean(community_sizes)),
        "std": float(np.std(community_sizes)),
        "min": int(np.min(community_sizes)),
        "max": int(np.max(community_sizes)),
        "sizes": community_sizes,
    }

    # Modularity
    try:
        metrics["modularity"] = modularity(graph, community_lists)
    except Exception:
        metrics["modularity"] = None

    # Coverage (fraction of nodes in communities)
    total_nodes_in_communities = sum(len(comm) for comm in community_lists)
    metrics["coverage"] = total_nodes_in_communities / len(graph.nodes()) if graph.nodes() else 0
    metrics["avg_community_size"] = float(np.mean(community_sizes))

    internal_edges = 0
    for community in community_lists:
        internal_edges += graph.subgraph(community).number_of_edges()
    total_edges = graph.number_of_edges()
    metrics["internal_edge_ratio"] = internal_edges / total_edges if total_edges else 0.0

    # Conductance for each community
    conductances = []
    for community in community_lists:
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
    for community in community_lists:
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
