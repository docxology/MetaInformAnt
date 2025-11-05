"""Community detection algorithms for biological networks."""

from __future__ import annotations

import random
from collections import Counter, defaultdict
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


def _leiden_communities(network: BiologicalNetwork, resolution: float = 1.0, max_iter: int = 100) -> Dict[str, int]:
    """Leiden algorithm implementation for community detection.
    
    The Leiden algorithm improves upon Louvain by guaranteeing well-connected
    communities through a local moving and refinement phase. It ensures that
    communities are not only locally optimal but also well-connected.
    
    Args:
        network: Input biological network
        resolution: Resolution parameter for modularity optimization
        max_iter: Maximum iterations for refinement phase
        
    Returns:
        Dictionary mapping node identifier to community ID
        
    References:
        Traag, V. A., et al. (2019). From Louvain to Leiden: guaranteeing
        well-connected communities. Scientific Reports, 9(1), 5233.
    """
    if not network.nodes:
        return {}
    
    # Phase 1: Local moving (Louvain-like optimization)
    communities = _louvain_communities(network, resolution, max_iter=max_iter)
    
    # Phase 2: Refinement - ensure well-connected communities
    refined = True
    iteration = 0
    
    while refined and iteration < max_iter:
        refined = False
        iteration += 1
        
        # Randomize node order for refinement
        nodes = list(network.nodes)
        random.shuffle(nodes)
        
        # Try to refine each node's community assignment
        for node in nodes:
            current_comm = communities[node]
            
            # Get neighboring communities
            neighbor_comms = set()
            for neighbor in network.get_neighbors(node):
                neighbor_comms.add(communities.get(neighbor, -1))
            
            best_comm = current_comm
            best_gain = 0.0
            
            # Try moving to each neighboring community
            for target_comm in neighbor_comms:
                if target_comm == current_comm:
                    continue
                
                gain = _modularity_gain(network, node, current_comm, target_comm, communities, resolution)
                
                if gain > best_gain:
                    best_gain = gain
                    best_comm = target_comm
            
            # Also try creating singleton community if beneficial
            max_comm_id = max(communities.values()) if communities else -1
            singleton_comm = max_comm_id + 1
            singleton_gain = _modularity_gain(
                network, node, current_comm, singleton_comm, communities, resolution
            )
            
            if singleton_gain > best_gain:
                best_gain = singleton_gain
                best_comm = singleton_comm
            
            # Move node if improvement found
            if best_comm != current_comm and best_gain > 1e-10:  # Threshold to avoid floating point issues
                communities[node] = best_comm
                refined = True
    
    return _renumber_communities(communities)


def _modularity_gain(
    network: BiologicalNetwork, node: str, from_comm: int, to_comm: int, communities: Dict[str, int], resolution: float
) -> float:
    """Calculate modularity gain from moving node between communities.
    
    Uses incremental delta-Q formula for efficiency:
    ΔQ = [Σ_in - Σ_tot * ki / (2m)] - [Σ_in_from - Σ_tot_from * ki / (2m)]
    where Σ_in is sum of edge weights within community,
    Σ_tot is sum of edge weights incident to community,
    ki is degree of node, and m is total edge weight.
    
    Args:
        network: Input biological network
        node: Node to move
        from_comm: Source community ID
        to_comm: Target community ID
        communities: Current community partition
        resolution: Resolution parameter
        
    Returns:
        Modularity gain (positive = improvement, negative = degradation)
    """
    if not network.edges:
        return 0.0
    
    # Calculate total edge weight
    total_weight = sum(network.edges.values())
    if total_weight == 0:
        return 0.0
    
    # Calculate node degree
    node_degree = 0.0
    for neighbor in network.get_neighbors(node):
        weight = network.get_edge_weight(node, neighbor) or 1.0
        node_degree += weight
    
    if node_degree == 0:
        return 0.0
    
    # Calculate community statistics using delta-Q formula
    # ΔQ = Σ_in_to - Σ_in_from - γ/(2m) * [Σ_tot_to^2 - Σ_tot_from^2]
    # where Σ_in is sum of weights of edges from node to community,
    # Σ_tot is sum of degrees in community
    
    # Calculate edges from node to each community
    edges_to_from = 0.0  # Edges from node to from_comm
    edges_to_to = 0.0  # Edges from node to to_comm
    
    # Calculate total degrees in each community
    from_comm_tot_degree = 0.0
    to_comm_tot_degree = 0.0
    
    for neighbor in network.get_neighbors(node):
        neighbor_comm = communities.get(neighbor, -1)
        edge_weight = network.get_edge_weight(node, neighbor) or 1.0
        
        if neighbor_comm == from_comm:
            edges_to_from += edge_weight
        elif neighbor_comm == to_comm:
            edges_to_to += edge_weight
    
    # Calculate total degrees in communities (excluding node)
    for n in network.nodes:
        if n == node:
            continue
        n_comm = communities.get(n, -1)
        n_degree = sum(
            network.get_edge_weight(n, neighbor) or 1.0
            for neighbor in network.get_neighbors(n)
        )
        
        if n_comm == from_comm:
            from_comm_tot_degree += n_degree
        elif n_comm == to_comm:
            to_comm_tot_degree += n_degree
    
    # Calculate delta-Q using formula
    # ΔQ = [Σ_in_to - γ/(2m) * (Σ_tot_to + ki)^2] - [Σ_in_from - γ/(2m) * (Σ_tot_from + ki)^2]
    # Simplified: ΔQ = Σ_in_to - Σ_in_from - γ/(2m) * [2*ki*(Σ_tot_to - Σ_tot_from) + ki^2]
    
    term_from = edges_to_from - resolution * (from_comm_tot_degree * node_degree) / (2 * total_weight)
    term_to = edges_to_to - resolution * (to_comm_tot_degree * node_degree) / (2 * total_weight)
    
    gain = term_to - term_from
    
    return gain / total_weight


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


def hierarchical_communities(
    network: BiologicalNetwork, levels: int = 3, resolution: float = 1.0, seed: Optional[int] = None
) -> Dict[int, Dict[str, int]]:
    """Detect hierarchical community structure at multiple levels.
    
    Performs community detection at multiple resolution levels to reveal
    hierarchical organization. Higher levels show larger communities,
    lower levels show more fine-grained structure.
    
    Args:
        network: Input biological network
        levels: Number of hierarchy levels to compute (default 3)
        resolution: Base resolution parameter (will be scaled for each level)
        seed: Random seed for reproducible results
        
    Returns:
        Dictionary mapping level (0 to levels-1) to community partition.
        Level 0 has highest resolution (smallest communities),
        level (levels-1) has lowest resolution (largest communities).
        
    Examples:
        >>> network = create_network(["A", "B", "C", "D", "E", "F"])
        >>> # Add edges to create hierarchy
        >>> hier = hierarchical_communities(network, levels=3)
        >>> len(hier[0])  # Finest level
        6
        >>> len(set(hier[2].values()))  # Coarsest level
        2
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    hierarchies = {}
    
    for level in range(levels):
        # Scale resolution: higher levels = lower resolution = larger communities
        level_resolution = resolution * (2 ** (levels - 1 - level))
        
        communities = detect_communities(network, method="louvain", resolution=level_resolution, seed=seed)
        hierarchies[level] = communities
    
    return hierarchies


def community_stability(
    network: BiologicalNetwork, method: str = "louvain", resolution: float = 1.0, n_runs: int = 10, seed: Optional[int] = None
) -> Dict[str, Any]:
    """Assess stability of community detection across multiple runs.
    
    Runs community detection multiple times with different random seeds
    to assess how stable the community structure is. High stability indicates
    robust community structure.
    
    Args:
        network: Input biological network
        method: Community detection algorithm (default "louvain")
        resolution: Resolution parameter for modularity optimization
        n_runs: Number of runs to perform (default 10)
        seed: Base random seed (will be incremented for each run)
        
    Returns:
        Dictionary containing:
        - stability_score: Average normalized mutual information across runs
        - modularity_scores: List of modularity values for each run
        - avg_modularity: Average modularity across runs
        - std_modularity: Standard deviation of modularity
        - consensus_communities: Consensus partition (most common assignment)
        - pairwise_agreement: Average pairwise agreement between runs
        
    Examples:
        >>> network = create_network(["A", "B", "C", "D", "E", "F"])
        >>> stability = community_stability(network, n_runs=5)
        >>> stability["stability_score"] > 0.5
        True
    """
    if n_runs < 2:
        raise ValueError("n_runs must be at least 2 for stability assessment")
    
    all_communities = []
    modularity_scores = []
    
    base_seed = seed if seed is not None else 42
    
    for run in range(n_runs):
        run_seed = base_seed + run if seed is not None else None
        communities = detect_communities(network, method=method, resolution=resolution, seed=run_seed)
        all_communities.append(communities)
        
        mod = modularity(network, communities, resolution)
        modularity_scores.append(mod)
    
    # Calculate pairwise agreement (normalized mutual information)
    agreements = []
    for i in range(len(all_communities)):
        for j in range(i + 1, len(all_communities)):
            agreement = _normalized_mutual_information(all_communities[i], all_communities[j])
            agreements.append(agreement)
    
    stability_score = np.mean(agreements) if agreements else 0.0
    
    # Consensus communities (most frequent assignment)
    consensus = _consensus_communities(all_communities)
    
    return {
        "stability_score": stability_score,
        "modularity_scores": modularity_scores,
        "avg_modularity": np.mean(modularity_scores),
        "std_modularity": np.std(modularity_scores),
        "consensus_communities": consensus,
        "pairwise_agreement": np.mean(agreements) if agreements else 0.0,
    }


def compare_communities(communities1: Dict[str, int], communities2: Dict[str, int]) -> Dict[str, float]:
    """Compare two community partitions.
    
    Computes similarity metrics between two community assignments
    to assess how similar they are.
    
    Args:
        communities1: First community partition (node -> community_id)
        communities2: Second community partition (node -> community_id)
        
    Returns:
        Dictionary containing:
        - normalized_mutual_information: NMI score (0 to 1, higher = more similar)
        - adjusted_rand_index: ARI score (-1 to 1, higher = more similar)
        - jaccard_similarity: Jaccard similarity of community assignments
        
    Examples:
        >>> comm1 = {"A": 0, "B": 0, "C": 1, "D": 1}
        >>> comm2 = {"A": 0, "B": 1, "C": 1, "D": 1}
        >>> comparison = compare_communities(comm1, comm2)
        >>> comparison["normalized_mutual_information"] > 0.0
        True
    """
    nmi = _normalized_mutual_information(communities1, communities2)
    ari = _adjusted_rand_index(communities1, communities2)
    jaccard = _community_jaccard(communities1, communities2)
    
    return {
        "normalized_mutual_information": nmi,
        "adjusted_rand_index": ari,
        "jaccard_similarity": jaccard,
    }


def optimize_resolution(
    network: BiologicalNetwork, resolution_range: Tuple[float, float] = (0.1, 2.0), n_points: int = 20, method: str = "louvain"
) -> Dict[str, Any]:
    """Find optimal resolution parameter for community detection.
    
    Tests multiple resolution values and selects the one that maximizes
    modularity or other quality metrics.
    
    Args:
        network: Input biological network
        resolution_range: (min, max) resolution values to test
        n_points: Number of resolution values to test
        method: Community detection algorithm
        
    Returns:
        Dictionary containing:
        - optimal_resolution: Resolution value with highest modularity
        - optimal_modularity: Modularity at optimal resolution
        - resolution_values: List of tested resolution values
        - modularity_values: List of modularity values for each resolution
        - optimal_communities: Community partition at optimal resolution
        
    Examples:
        >>> network = create_network(["A", "B", "C", "D", "E", "F"])
        >>> result = optimize_resolution(network, resolution_range=(0.5, 2.0), n_points=10)
        >>> result["optimal_resolution"] > 0.0
        True
    """
    resolution_min, resolution_max = resolution_range
    resolution_values = np.linspace(resolution_min, resolution_max, n_points)
    modularity_values = []
    all_communities = []
    
    for res in resolution_values:
        communities = detect_communities(network, method=method, resolution=res)
        mod = modularity(network, communities, resolution=res)
        modularity_values.append(mod)
        all_communities.append(communities)
    
    # Find optimal resolution
    optimal_idx = np.argmax(modularity_values)
    optimal_resolution = resolution_values[optimal_idx]
    optimal_modularity = modularity_values[optimal_idx]
    optimal_communities = all_communities[optimal_idx]
    
    return {
        "optimal_resolution": optimal_resolution,
        "optimal_modularity": optimal_modularity,
        "resolution_values": resolution_values.tolist(),
        "modularity_values": modularity_values,
        "optimal_communities": optimal_communities,
    }


def _normalized_mutual_information(communities1: Dict[str, int], communities2: Dict[str, int]) -> float:
    """Calculate normalized mutual information between two partitions."""
    # Get all nodes
    all_nodes = set(communities1.keys()).union(set(communities2.keys()))
    
    if not all_nodes:
        return 1.0
    
    # Count co-occurrence matrix
    comm1_ids = sorted(set(communities1.values()))
    comm2_ids = sorted(set(communities2.values()))
    
    contingency = np.zeros((len(comm1_ids), len(comm2_ids)))
    comm1_map = {cid: i for i, cid in enumerate(comm1_ids)}
    comm2_map = {cid: i for i, cid in enumerate(comm2_ids)}
    
    for node in all_nodes:
        c1 = communities1.get(node, -1)
        c2 = communities2.get(node, -1)
        if c1 != -1 and c2 != -1:
            i = comm1_map.get(c1, 0)
            j = comm2_map.get(c2, 0)
            contingency[i, j] += 1
    
    # Calculate NMI
    n = len(all_nodes)
    if n == 0:
        return 0.0
    
    # Marginal sums
    sum_comm1 = np.sum(contingency, axis=1)
    sum_comm2 = np.sum(contingency, axis=0)
    
    # Mutual information
    mi = 0.0
    for i in range(len(comm1_ids)):
        for j in range(len(comm2_ids)):
            if contingency[i, j] > 0:
                p_ij = contingency[i, j] / n
                p_i = sum_comm1[i] / n
                p_j = sum_comm2[j] / n
                
                if p_i > 0 and p_j > 0:
                    mi += p_ij * np.log2(p_ij / (p_i * p_j))
    
    # Entropies
    h1 = -np.sum([(s / n) * np.log2(s / n) for s in sum_comm1 if s > 0])
    h2 = -np.sum([(s / n) * np.log2(s / n) for s in sum_comm2 if s > 0])
    
    # Normalized MI
    if h1 + h2 == 0:
        return 1.0
    
    nmi = 2 * mi / (h1 + h2)
    return nmi


def _adjusted_rand_index(communities1: Dict[str, int], communities2: Dict[str, int]) -> float:
    """Calculate adjusted Rand index between two partitions."""
    all_nodes = set(communities1.keys()).union(set(communities2.keys()))
    
    if not all_nodes:
        return 1.0
    
    # Count pairs in same community
    same_in_1 = 0
    same_in_2 = 0
    same_in_both = 0
    total_pairs = 0
    
    nodes_list = list(all_nodes)
    for i, node1 in enumerate(nodes_list):
        for node2 in nodes_list[i + 1 :]:
            total_pairs += 1
            
            same1 = communities1.get(node1, -1) == communities1.get(node2, -1)
            same2 = communities2.get(node1, -1) == communities2.get(node2, -1)
            
            if same1:
                same_in_1 += 1
            if same2:
                same_in_2 += 1
            if same1 and same2:
                same_in_both += 1
    
    if total_pairs == 0:
        return 1.0
    
    # Adjusted Rand Index
    expected = (same_in_1 * same_in_2) / total_pairs
    max_possible = (same_in_1 + same_in_2) / 2
    
    if max_possible == expected:
        return 1.0
    
    ari = (same_in_both - expected) / (max_possible - expected) if max_possible != expected else 0.0
    return ari


def _community_jaccard(communities1: Dict[str, int], communities2: Dict[str, int]) -> float:
    """Calculate Jaccard similarity of community assignments."""
    all_nodes = set(communities1.keys()).union(set(communities2.keys()))
    
    if not all_nodes:
        return 1.0
    
    agreements = 0
    total = 0
    
    for node1 in all_nodes:
        for node2 in all_nodes:
            if node1 >= node2:
                continue
            
            total += 1
            same1 = communities1.get(node1, -1) == communities1.get(node2, -1)
            same2 = communities2.get(node1, -1) == communities2.get(node2, -1)
            
            if same1 == same2:
                agreements += 1
    
    return agreements / total if total > 0 else 1.0


def _consensus_communities(all_communities: List[Dict[str, int]]) -> Dict[str, int]:
    """Create consensus partition from multiple community detections."""
    if not all_communities:
        return {}
    
    # Get all nodes
    all_nodes = set()
    for communities in all_communities:
        all_nodes.update(communities.keys())
    
    # For each node, find most common community assignment
    consensus = {}
    for node in all_nodes:
        assignments = [comm.get(node, -1) for comm in all_communities]
        
        # Find most common assignment
        counts = Counter(assignments)
        most_common = counts.most_common(1)[0][0]
        consensus[node] = most_common
    
    return _renumber_communities(consensus)
