"""Community detection algorithms for biological networks."""

from __future__ import annotations

import random
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .graph import BiologicalNetwork


def detect_communities(
    network: BiologicalNetwork, method: str = "louvain", resolution: float = 1.0, seed: Optional[int] = None
) -> Dict[str, int]:
    """Detect communities (modules) in a biological network.
    
    Identifies groups of nodes that are more densely connected internally
    than externally. Community detection is useful for finding functional
    modules in biological networks (protein complexes, gene modules, etc.).
    
    Args:
        network: Input biological network
        method: Community detection algorithm:
            - "louvain": Fast modularity optimization (default)
            - "greedy": Greedy modularity maximization
            - "leiden": Leiden algorithm (improved over Louvain)
        resolution: Resolution parameter for modularity optimization.
            Higher values find smaller, more numerous communities.
            Lower values find larger, fewer communities. Default 1.0.
        seed: Random seed for reproducible results
        
    Returns:
        Dictionary mapping node identifier to community ID (integer).
        Nodes in the same community share the same ID.
        
    Examples:
        >>> network = create_network(["A", "B", "C", "D", "E"])
        >>> network.add_edge("A", "B"); network.add_edge("B", "C")
        >>> network.add_edge("D", "E")
        >>> communities = detect_communities(network, method="louvain")
        >>> len(set(communities.values()))  # Number of communities
        2
        
    References:
        Blondel, V. D., et al. (2008). Fast unfolding of communities in
        large networks. Journal of Statistical Mechanics, 2008(10), P10008.
        
        Traag, V. A., et al. (2019). From Louvain to Leiden: guaranteeing
        well-connected communities. Scientific Reports, 9(1), 5233.
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    if method == "louvain":
        return _louvain_communities(network, resolution)
    elif method == "greedy":
        return _greedy_modularity_communities(network)
    elif method == "leiden":
        return _leiden_communities(network, resolution)
    else:
        raise ValueError(f"Unknown community detection method: {method}")


def modularity(network: BiologicalNetwork, communities: Dict[str, int], resolution: float = 1.0) -> float:
    """Calculate modularity of a network partition.
    
    Modularity measures how much more connected nodes within communities
    are compared to a random network with the same degree distribution.
    Higher modularity indicates better community structure.
    
    Args:
        network: Input biological network
        communities: Dictionary mapping node identifier to community ID
        resolution: Resolution parameter for modularity. Higher values
            favor smaller communities, lower values favor larger communities.
            Default 1.0.
            
    Returns:
        Modularity score. Typically ranges from -1 to 1, but can exceed 1
        for very modular networks. Positive values indicate good partitioning.
        Formula: Q = (1/2m) Σᵢⱼ [Aᵢⱼ - (γ kᵢ kⱼ)/(2m)] δ(cᵢ, cⱼ)
        where m is total edge weight, A is adjacency, k is degree,
        γ is resolution, and δ is community indicator.
        
    Examples:
        >>> network = create_network(["A", "B", "C", "D"], directed=False)
        >>> network.add_edge("A", "B"); network.add_edge("C", "D")
        >>> communities = {"A": 0, "B": 0, "C": 1, "D": 1}
        >>> mod = modularity(network, communities)
        >>> mod > 0.0  # Should have positive modularity
        True
        
    References:
        Newman, M. E. J. (2006). Modularity and community structure in networks.
        Proceedings of the National Academy of Sciences, 103(23), 8577-8582.
    """
    if not network.edges:
        return 0.0

    # Calculate degrees
    degrees = {}
    total_weight = 0.0

    for node in network.nodes:
        degree = 0.0
        for neighbor in network.get_neighbors(node):
            weight = network.get_edge_weight(node, neighbor) or 1.0
            degree += weight
            if not network.directed:
                total_weight += weight / 2  # Avoid double counting
            else:
                total_weight += weight
        degrees[node] = degree

    if total_weight == 0:
        return 0.0

    # Calculate modularity
    Q = 0.0
    for edge, weight in network.edges.items():
        node1, node2 = edge

        # Same community bonus
        if communities.get(node1) == communities.get(node2):
            Q += weight

        # Expected edge weight under null model
        expected = (degrees[node1] * degrees[node2]) / (2 * total_weight)
        if communities.get(node1) == communities.get(node2):
            Q -= resolution * expected

    return Q / total_weight


def community_metrics(network: BiologicalNetwork, communities: Dict[str, int]) -> Dict[str, any]:
    """Calculate comprehensive metrics for network community structure.
    
    Evaluates the quality and characteristics of a community partition,
    including community sizes, edge distribution, and modularity.
    
    Args:
        network: Input biological network
        communities: Dictionary mapping node identifier to community ID
        
    Returns:
        Dictionary containing:
        - num_communities: Number of distinct communities
        - avg_community_size: Average number of nodes per community
        - community_sizes: Dictionary mapping community_id -> size
        - internal_edges: Number of edges within communities
        - external_edges: Number of edges between communities
        - internal_edge_ratio: Fraction of edges that are internal
        - modularity: Modularity score of the partition
        
    Examples:
        >>> network = create_network(["A", "B", "C", "D"], directed=False)
        >>> network.add_edge("A", "B"); network.add_edge("C", "D")
        >>> communities = {"A": 0, "B": 0, "C": 1, "D": 1}
        >>> metrics = community_metrics(network, communities)
        >>> metrics["num_communities"]
        2
        >>> metrics["internal_edge_ratio"]
        1.0  # All edges are internal
    """
    # Count communities and their sizes
    community_sizes = defaultdict(int)
    for node, comm_id in communities.items():
        community_sizes[comm_id] += 1

    num_communities = len(community_sizes)
    avg_size = np.mean(list(community_sizes.values())) if community_sizes else 0

    # Calculate internal and external edges
    internal_edges = 0
    external_edges = 0

    for edge, weight in network.edges.items():
        node1, node2 = edge
        comm1 = communities.get(node1, -1)
        comm2 = communities.get(node2, -1)

        if comm1 == comm2 and comm1 != -1:
            internal_edges += 1
        else:
            external_edges += 1

    total_edges = internal_edges + external_edges
    internal_ratio = internal_edges / total_edges if total_edges > 0 else 0

    return {
        "num_communities": num_communities,
        "avg_community_size": avg_size,
        "community_sizes": dict(community_sizes),
        "internal_edges": internal_edges,
        "external_edges": external_edges,
        "internal_edge_ratio": internal_ratio,
        "modularity": modularity(network, communities),
    }


def _louvain_communities(network: BiologicalNetwork, resolution: float = 1.0, max_iter: int = 100) -> Dict[str, int]:
    """Louvain community detection algorithm.

    Simplified implementation of the Louvain method for modularity optimization.
    """
    # Initialize each node in its own community
    communities = {node: i for i, node in enumerate(network.nodes)}

    improved = True
    iteration = 0

    while improved and iteration < max_iter:
        improved = False
        iteration += 1

        # Randomize node order
        nodes = list(network.nodes)
        random.shuffle(nodes)

        for node in nodes:
            current_comm = communities[node]
            best_comm = current_comm
            best_gain = 0.0

            # Get neighboring communities
            neighbor_comms = set()
            for neighbor in network.get_neighbors(node):
                neighbor_comms.add(communities[neighbor])

            # Try moving to each neighboring community
            for target_comm in neighbor_comms:
                if target_comm == current_comm:
                    continue

                # Calculate modularity gain
                gain = _modularity_gain(network, node, current_comm, target_comm, communities, resolution)

                if gain > best_gain:
                    best_gain = gain
                    best_comm = target_comm

            # Move node if improvement found
            if best_comm != current_comm:
                communities[node] = best_comm
                improved = True

    # Renumber communities to be contiguous
    return _renumber_communities(communities)


def _greedy_modularity_communities(network: BiologicalNetwork) -> Dict[str, int]:
    """Greedy modularity optimization (simplified Clauset-Newman-Moore)."""
    # Start with each node in its own community
    communities = {node: i for i, node in enumerate(network.nodes)}

    # Build initial community network
    community_graph = _build_community_network(network, communities)

    # Greedily merge communities to maximize modularity
    while len(community_graph) > 1:
        best_merge = None
        best_gain = -float("inf")

        # Find best community pair to merge
        for comm1 in community_graph:
            for comm2 in community_graph:
                if comm1 >= comm2:
                    continue

                # Calculate potential modularity gain
                gain = _merge_gain(network, communities, comm1, comm2)
                if gain > best_gain:
                    best_gain = gain
                    best_merge = (comm1, comm2)

        if best_merge is None or best_gain <= 0:
            break

        # Perform merge
        comm1, comm2 = best_merge
        for node, comm in communities.items():
            if comm == comm2:
                communities[node] = comm1

        # Rebuild community network
        community_graph = _build_community_network(network, communities)

    return _renumber_communities(communities)


def _leiden_communities(network: BiologicalNetwork, resolution: float = 1.0) -> Dict[str, int]:
    """Simplified Leiden algorithm implementation.

    This is a basic version - the full Leiden algorithm is more complex.
    """
    # Start with Louvain
    communities = _louvain_communities(network, resolution)

    # Refinement phase (simplified)
    refined = True
    while refined:
        refined = False

        nodes = list(network.nodes)
        random.shuffle(nodes)

        for node in nodes:
            current_comm = communities[node]

            # Try creating singleton community
            singleton_gain = _modularity_gain(
                network, node, current_comm, max(communities.values()) + 1, communities, resolution
            )

            if singleton_gain > 0:
                communities[node] = max(communities.values()) + 1
                refined = True

    return _renumber_communities(communities)


def _modularity_gain(
    network: BiologicalNetwork, node: str, from_comm: int, to_comm: int, communities: Dict[str, int], resolution: float
) -> float:
    """Calculate modularity gain from moving node between communities."""
    # This is a simplified calculation
    # Full implementation would need more detailed bookkeeping

    # Current modularity
    current_mod = modularity(network, communities, resolution)

    # Test modularity after move
    test_communities = communities.copy()
    test_communities[node] = to_comm
    new_mod = modularity(network, test_communities, resolution)

    return new_mod - current_mod


def _build_community_network(network: BiologicalNetwork, communities: Dict[str, int]) -> Set[int]:
    """Build set of community IDs from current communities."""
    return set(communities.values())


def _merge_gain(network: BiologicalNetwork, communities: Dict[str, int], comm1: int, comm2: int) -> float:
    """Calculate modularity gain from merging two communities."""
    # Current modularity
    current_mod = modularity(network, communities)

    # Test modularity after merge
    test_communities = communities.copy()
    for node, comm in test_communities.items():
        if comm == comm2:
            test_communities[node] = comm1

    new_mod = modularity(network, test_communities)
    return new_mod - current_mod


def _renumber_communities(communities: Dict[str, int]) -> Dict[str, int]:
    """Renumber communities to be contiguous starting from 0."""
    unique_comms = sorted(set(communities.values()))
    comm_mapping = {old: new for new, old in enumerate(unique_comms)}

    return {node: comm_mapping[comm] for node, comm in communities.items()}
