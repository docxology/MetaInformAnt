"""Network information analysis for METAINFORMANT.

This module provides information-theoretic measures for analyzing
biological networks, including network entropy, information flow,
and information-based community detection.
"""

from __future__ import annotations

import importlib.util
import math
from collections import Counter
from itertools import combinations, permutations
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core.utils import logging
from metainformant.information.metrics.core import syntactic
from metainformant.information.network_info.information_flow import network_entropy as _laplacian_network_entropy

logger = logging.get_logger(__name__)

# Optional network analysis dependencies
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    logger.warning("networkx not available, network information analysis disabled")


def network_entropy(graph: Any, attribute: Optional[str] = None, method: str = "shannon") -> float:
    """Calculate entropy of a biological network.

    Args:
        graph: NetworkX graph or adjacency matrix
        attribute: Node attribute to use for entropy calculation (``method="shannon"`` only)
        method: Entropy calculation method. ``'shannon'`` is the entropy of the degree distribution (or of the
            node ``attribute`` values), in bits. ``'von_neumann'`` is the Von Neumann entropy, in bits, of the
            Laplacian density matrix ``rho = L / trace(L)`` (Braunstein, Ghosh and Severini, 2006), the same
            value as ``metainformant.information.network_info.information_flow.network_entropy``: a complete
            graph on ``n`` nodes has ``log2(n - 1)``, and a graph without edges has 0. It is defined for
            undirected graphs only.

    Returns:
        Network entropy value

    Raises:
        ImportError: If networkx not available
        ValueError: If graph format not supported, ``method`` is unknown, or ``method="von_neumann"`` is given a
            directed graph
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network entropy calculation")

    # Convert to NetworkX graph if needed
    if isinstance(graph, np.ndarray):
        G = nx.from_numpy_array(graph)
    elif hasattr(graph, "nodes"):  # NetworkX graph
        G = graph
    else:
        raise ValueError("Graph must be NetworkX graph or numpy adjacency matrix")

    if method == "shannon":
        return _shannon_network_entropy(G, attribute)
    elif method == "von_neumann":
        return _von_neumann_entropy(G)
    else:
        raise ValueError(f"Unknown entropy method: {method}")


def _shannon_network_entropy(G: Any, attribute: Optional[str] = None) -> float:
    """Calculate Shannon entropy of network structure or attributes."""
    if attribute is None:
        # Structural entropy based on degree distribution
        total_nodes = len(G)
        if total_nodes == 0:
            return 0.0

        degrees = [d for n, d in G.degree()]

        degree_counts = Counter(degrees)

        # Convert to probabilities
        probs = [count / total_nodes for count in degree_counts.values()]

        # Calculate entropy
        entropy = 0.0
        for p in probs:
            if p > 0:
                entropy -= p * math.log2(p)

        return entropy

    else:
        # Attribute-based entropy
        if not nx.get_node_attributes(G, attribute):
            raise ValueError(f"Attribute '{attribute}' not found in graph nodes")

        attr_values = list(nx.get_node_attributes(G, attribute).values())

        attr_counts = Counter(attr_values)

        # Calculate entropy
        total = sum(attr_counts.values())
        entropy = 0.0
        for count in attr_counts.values():
            p = count / total
            entropy -= p * math.log2(p)

        return entropy


def _von_neumann_entropy(G: Any) -> float:
    """Von Neumann entropy (bits) of the density matrix ``rho = L / trace(L)`` of an undirected graph.

    ``L = D - W`` is the (weighted) combinatorial Laplacian, so ``rho`` is positive semidefinite with unit trace and
    ``S = -sum(lambda * log2(lambda))`` runs over its eigenvalues (Braunstein, Ghosh and Severini, 2006). The
    definition and the eigendecomposition are those of
    :func:`metainformant.information.network_info.information_flow.network_entropy`; this function only converts the
    NetworkX graph to its weighted adjacency matrix (a missing ``weight`` counts as 1, parallel edges add up).

    An earlier version took the eigenvalues of the row-normalised adjacency (random-walk transition) matrix and kept
    the positive ones. That matrix is not a density matrix: its spectrum sums to ``trace(P)`` rather than 1 and its
    top eigenvalue is always 1, which contributes ``-1 * log2(1) = 0``. Every complete graph and every star came out
    at 0 instead of ``log2(n - 1)`` and the closed form for a star.

    Raises:
        ValueError: If ``G`` is directed. The density matrix is defined for symmetric weights; convert explicitly
            (``G.to_undirected()``) to analyse the undirected skeleton.
    """
    if G.is_directed():
        raise ValueError(
            "von Neumann entropy is defined for undirected graphs: the Laplacian density matrix needs symmetric "
            "weights. Pass G.to_undirected() to analyse the undirected skeleton."
        )
    adjacency = nx.to_numpy_array(G)
    return float(_laplacian_network_entropy(adjacency.tolist())["entropy"])


def _flow_nodes(G: Any, names: Optional[List[Any]], role: str) -> List[Any]:
    """The graph nodes named in ``names`` (every node when ``names`` is None), without duplicates.

    Raises:
        ValueError: If ``names`` is empty or names a node that is not in the graph. Unknown names used to be skipped
            silently: with no valid source the diffusion signal was all zeros and its renormalisation produced NaN.
    """
    if names is None:
        return list(G.nodes())
    if len(names) == 0:
        raise ValueError(f"{role}_nodes must not be an empty list; pass None to use every node")
    present = set(G.nodes())
    unknown = [name for name in names if name not in present]
    if unknown:
        raise ValueError(f"{role}_nodes not in the graph: {unknown[:5]}")
    return list(dict.fromkeys(names))


def information_flow(
    graph: Any,
    source_nodes: Optional[List[str]] = None,
    target_nodes: Optional[List[str]] = None,
    method: str = "random_walk",
    steps: int = 100,
) -> Dict[str, Any]:
    """Calculate information flow in a biological network.

    Args:
        graph: NetworkX graph or adjacency matrix
        source_nodes: Source nodes for information flow (None: every node)
        target_nodes: Target nodes for information flow (None: every node)
        method: Information flow method ('random_walk', 'diffusion')
        steps: Number of steps for flow calculation

    Returns:
        Dictionary with information flow results. For ``'random_walk'``, ``flow_matrix[source][target]`` is the
        expected number of visits to ``target`` within ``steps`` steps of a random walk that starts at ``source``
        (each source is walked on its own), and ``path_length_entropy`` is the Shannon entropy of the normalised
        flows. For ``'diffusion'``, ``node_flows`` holds the signal of each target node over time.

    Raises:
        ImportError: If networkx not available
        ValueError: If a source or target node is not in the graph, or a node list is empty
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for information flow calculation")

    # Convert to NetworkX graph if needed
    if isinstance(graph, np.ndarray):
        G = nx.from_numpy_array(graph)
    elif hasattr(graph, "nodes"):
        G = graph
    else:
        raise ValueError("Graph must be NetworkX graph or numpy adjacency matrix")

    if method == "random_walk":
        return _random_walk_information_flow(G, source_nodes, target_nodes, steps)
    elif method == "diffusion":
        return _diffusion_information_flow(G, source_nodes, target_nodes, steps)
    else:
        raise ValueError(f"Unknown flow method: {method}")


def _random_walk_information_flow(
    G: Any, source_nodes: Optional[List[str]] = None, target_nodes: Optional[List[str]] = None, steps: int = 100
) -> Dict[str, Any]:
    """Calculate information flow using random walk model.

    Each source is walked separately: ``flow_matrix[s][t]`` is the sum over ``steps`` steps of the probability that a
    walk started at ``s`` is at ``t``. (An earlier version walked the mixture of all sources, so every source got the
    same flow to a target, and it recorded nothing at all unless ``target_nodes`` was given or the sources were
    listed, returning an all-zero matrix and a path-length entropy of 0 for the default arguments.)
    """
    node_list = list(G.nodes())
    node_to_idx = {node: i for i, node in enumerate(node_list)}
    sources = _flow_nodes(G, source_nodes, "source")
    targets = _flow_nodes(G, target_nodes, "target")

    # Transition matrix of the walk (an isolated node has no outgoing step)
    A = nx.to_numpy_array(G, nodelist=node_list)
    degrees = np.sum(A, axis=1)
    degrees = np.where(degrees == 0, 1, degrees)
    P = A / degrees[:, np.newaxis]

    # visits[i, j] = sum_{k=1..steps} (P^k)[i, j]: all sources are walked at once as the rows of the power
    n_nodes = len(node_list)
    power = np.eye(n_nodes)
    visits = np.zeros((n_nodes, n_nodes))
    for _ in range(steps):
        power = power @ P
        visits += power

    flow_dict: Dict[str, Any] = {
        source: {target: float(visits[node_to_idx[source], node_to_idx[target]]) for target in targets}
        for source in sources
    }

    # Calculate path length entropy from the flow distribution
    all_flows = [flow for source_flows in flow_dict.values() for flow in source_flows.values()]
    path_length_entropy = 0.0
    total_flow = sum(abs(f) for f in all_flows)
    if total_flow > 0:
        probs = [abs(f) / total_flow for f in all_flows if abs(f) > 0]
        for p in probs:
            path_length_entropy -= p * math.log2(p)

    return {
        "flow_matrix": flow_dict,
        "path_length_entropy": path_length_entropy,
        "method": "random_walk",
        "steps": steps,
        "source_nodes": source_nodes,
        "target_nodes": target_nodes,
    }


def _diffusion_information_flow(
    G: Any, source_nodes: Optional[List[str]] = None, target_nodes: Optional[List[str]] = None, steps: int = 100
) -> Dict[str, Any]:
    """Calculate information flow using diffusion model.

    A unit signal, spread evenly over the source nodes, diffuses by explicit Euler steps of size 0.1 of the heat
    equation with the symmetric normalised Laplacian, clipped to non-negative values and renormalised to unit total
    after every step. ``node_flows`` reports the signal of the target nodes (every node when ``target_nodes`` is
    None); the diffusion itself always involves the whole graph.
    """
    node_list = list(G.nodes())
    sources = _flow_nodes(G, source_nodes, "source")
    targets = _flow_nodes(G, target_nodes, "target")

    A = nx.to_numpy_array(G, nodelist=node_list)
    n_nodes = len(node_list)

    # Create diffusion matrix (normalized Laplacian)
    degrees = np.sum(A, axis=1)
    degrees = np.where(degrees == 0, 1, degrees)

    # Laplacian matrix
    L = np.diag(degrees) - A

    # Normalized Laplacian (for regular graphs, this approximates diffusion)
    D_inv_sqrt = np.diag(1.0 / np.sqrt(degrees))
    L_norm = D_inv_sqrt @ L @ D_inv_sqrt

    # Initialize with source nodes
    node_to_idx = {node: i for i, node in enumerate(node_list)}
    if source_nodes is None:
        initial_signal = np.ones(n_nodes) / n_nodes
    else:
        initial_signal = np.zeros(n_nodes)
        for source in sources:
            initial_signal[node_to_idx[source]] = 1.0 / len(sources)

    # Diffusion simulation
    signal = initial_signal.copy()
    flow_history = []

    for step in range(steps):
        # Diffuse signal
        signal = signal - 0.1 * L_norm @ signal  # Simple diffusion step
        signal = np.maximum(signal, 0)  # Ensure non-negative
        signal = signal / np.sum(signal)  # Renormalize

        flow_history.append(signal.copy())

    # Calculate final flow patterns
    final_signal = flow_history[-1] if flow_history else initial_signal

    # Convert to node-wise results
    flow_results = {}

    for node in targets:
        i = node_to_idx[node]
        flow_results[node] = {
            "final_signal": float(final_signal[i]),
            "signal_history": [float(flow_history[t][i]) for t in range(len(flow_history))],
        }

    return {
        "node_flows": flow_results,
        "method": "diffusion",
        "steps": steps,
        "source_nodes": source_nodes,
        "target_nodes": target_nodes,
        "diffusion_matrix_shape": L_norm.shape,
    }


def information_community_detection(graph: Any, method: str = "infomap", **kwargs: Any) -> Dict[str, Any]:
    """Detect communities using information-theoretic approaches.

    Args:
        graph: NetworkX graph
        method: Community detection method ('infomap', 'map_equation')
        **kwargs: Additional parameters for detection algorithm

    Returns:
        Dictionary with community detection results

    Raises:
        ImportError: If required packages not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for community detection")

    has_community = importlib.util.find_spec("community") is not None

    if method == "infomap":
        return _infomap_community_detection(graph, **kwargs)
    elif method == "map_equation":
        return _map_equation_community_detection(graph, **kwargs)
    elif method == "louvain" and has_community:
        return _louvain_community_detection(graph, **kwargs)
    else:
        raise ValueError(f"Unknown or unavailable community detection method: {method}")


def _infomap_community_detection(graph: Any, **kwargs: Any) -> Dict[str, Any]:
    """InfoMap community detection algorithm.

    Uses the infomap package to detect communities based on information flow.

    Args:
        graph: NetworkX graph
        **kwargs: Additional arguments for infomap

    Returns:
        Dictionary with community detection results

    Raises:
        ImportError: If infomap package is not installed
    """
    try:
        import infomap  # type: ignore[import-not-found]
    except ImportError:
        raise ImportError(
            "infomap package is required for InfoMap community detection. "
            "Install with: pip install infomap (or: uv pip install infomap)"
        )

    # Create infomap object
    im = infomap.Infomap(**kwargs)

    # Add edges from graph
    node_to_id = {node: i for i, node in enumerate(graph.nodes())}
    id_to_node = {i: node for node, i in node_to_id.items()}

    for u, v in graph.edges():
        weight = graph[u][v].get("weight", 1.0)
        im.add_link(node_to_id[u], node_to_id[v], weight)

    # Run infomap
    im.run()

    # Extract communities
    communities: Dict[Any, Any] = {}
    for node in im.tree:
        if node.is_leaf:
            node_name = id_to_node[node.node_id]
            community_id = node.module_id
            if community_id not in communities:
                communities[community_id] = []
            communities[community_id].append(node_name)

    return {
        "method": "infomap",
        "status": "success",
        "communities": communities,
        "n_communities": len(communities),
        "codelength": im.codelength,
    }


def _map_equation_community_detection(graph: Any, **kwargs: Any) -> Dict[str, Any]:
    """Map equation community detection.

    Uses the infomap package (which implements the map equation) for
    community detection based on random walk dynamics.

    Args:
        graph: NetworkX graph
        **kwargs: Additional arguments

    Returns:
        Dictionary with community detection results

    Raises:
        ImportError: If infomap package is not installed
    """
    # Use infomap with two-level partition (map equation default)
    return _infomap_community_detection(graph, two_level=True, **kwargs)


def _louvain_community_detection(graph: Any, **kwargs: Any) -> Dict[str, Any]:
    """Louvain community detection with information-theoretic interpretation."""
    try:
        import community as community_louvain

        # Convert to undirected graph if needed
        if graph.is_directed():
            G = graph.to_undirected()
        else:
            G = graph

        # Run Louvain algorithm
        partition = community_louvain.best_partition(G, **kwargs)

        # Organize communities
        communities: Dict[Any, Any] = {}
        for node, community_id in partition.items():
            if community_id not in communities:
                communities[community_id] = []
            communities[community_id].append(node)

        # Calculate modularity
        modularity = community_louvain.modularity(partition, G)

        # Calculate community information measures
        community_info = {}
        for comm_id, nodes in communities.items():
            subgraph = G.subgraph(nodes)

            # Community size and density
            size = len(nodes)
            edges = subgraph.number_of_edges()
            max_edges = size * (size - 1) / 2 if not G.is_directed() else size * (size - 1)

            community_info[comm_id] = {
                "size": size,
                "edges": edges,
                "density": edges / max_edges if max_edges > 0 else 0,
                "nodes": nodes,
            }

        return {
            "method": "louvain",
            "status": "completed",
            "modularity": modularity,
            "n_communities": len(communities),
            "communities": communities,
            "community_info": community_info,
        }

    except ImportError:
        return {
            "method": "louvain",
            "status": "failed",
            "error": "python-louvain package required",
        }


def network_information_centrality(graph: Any, method: str = "entropy", normalized: bool = True) -> Dict[str, float]:
    """Calculate information-theoretic centrality measures for network nodes.

    Args:
        graph: NetworkX graph
        method: Centrality method ('entropy', 'information_flow')
        normalized: Whether to normalize centrality values

    Returns:
        Dictionary mapping node names to centrality values

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for centrality calculation")

    if isinstance(graph, np.ndarray):
        G = nx.from_numpy_array(graph)
    elif hasattr(graph, "nodes"):
        G = graph
    else:
        raise ValueError("Graph must be NetworkX graph or numpy adjacency matrix")

    if method == "entropy":
        return _entropy_centrality(G, normalized)
    elif method == "information_flow":
        return _information_flow_centrality(G, normalized)
    else:
        raise ValueError(f"Unknown centrality method: {method}")


def _entropy_centrality(G: Any, normalized: bool = True) -> Dict[str, float]:
    """Calculate entropy-based node centrality."""
    centrality = {}

    for node in G.nodes():
        # Calculate entropy of node's neighborhood
        neighbors = list(G.neighbors(node))
        if not neighbors:
            centrality[node] = 0.0
            continue

        # Degree distribution of neighbors
        neighbor_degrees = [G.degree(n) for n in neighbors]

        degree_counts = Counter(neighbor_degrees)

        # Calculate entropy
        total = len(neighbors)
        entropy = 0.0
        for count in degree_counts.values():
            p = count / total
            entropy -= p * math.log2(p) if p > 0 else 0

        centrality[node] = entropy

    # Normalize if requested
    if normalized and centrality:
        max_cent = max(centrality.values())
        if max_cent > 0:
            centrality = {node: cent / max_cent for node, cent in centrality.items()}

    return dict(centrality)


def _information_flow_centrality(G: Any, normalized: bool = True) -> Dict[str, float]:
    """Calculate information flow centrality."""
    # Use betweenness centrality as proxy for information flow
    centrality = nx.betweenness_centrality(G, normalized=normalized)

    # Convert to information-theoretic interpretation
    # Higher betweenness = more information flow through node
    return dict(centrality)


def network_motif_information(
    graph: Any,
    motif_size: int = 3,
    n_random: int = 100,
    random_state: Any = None,
    max_subgraphs: int = 50000,
) -> Dict[str, Any]:
    """Analyze network motifs using information-theoretic measures.

    Enumerates the connected subgraphs of the specified size, classifies
    them by isomorphism class, computes z-scores against a degree-preserving
    random graph null model, and calculates the entropy of the motif
    distribution. A directed graph is analysed as directed: its subgraphs are
    classified by their directed structure and the null model keeps in- and
    out-degrees.

    Args:
        graph: NetworkX graph (directed or undirected)
        motif_size: Size of motifs to analyze (3 or 4)
        n_random: Number of random graphs for null model
        random_state: Seed (int) or ``numpy.random.Generator`` for the null graphs and for the choice of subgraphs
            when there are more than ``max_subgraphs``. The same seed gives the same result; the global ``random``
            and ``numpy.random`` states are neither seeded nor consumed. ``None`` draws fresh entropy.
        max_subgraphs: Largest number of node combinations examined per graph. A graph with more combinations than
            this is sampled uniformly, which is reported as ``sampled`` and makes the motif counts counts of the
            sample (the z-scores stay comparable: observed and random graphs are sampled alike).

    Returns:
        Dictionary with motif counts, z-scores, entropy, and significance. ``sampled`` says whether the counts come
        from a sample of the node combinations, and ``n_combinations_examined`` how many were looked at.

    Raises:
        ImportError: If networkx not available
        ValueError: If motif size not supported
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for motif analysis")

    if motif_size not in [3, 4]:
        raise ValueError("Motif size must be 3 or 4")

    if isinstance(graph, np.ndarray):
        G = nx.from_numpy_array(graph)
    elif hasattr(graph, "nodes"):
        G = graph
    else:
        raise ValueError("Graph must be NetworkX graph or numpy adjacency matrix")

    if len(G) < motif_size:
        return {
            "motif_size": motif_size,
            "status": "completed",
            "motifs_found": {},
            "z_scores": {},
            "motif_entropy": 0.0,
            "n_subgraphs_enumerated": 0,
            "sampled": False,
            "n_combinations_examined": 0,
        }

    rng = random_state if isinstance(random_state, np.random.Generator) else np.random.default_rng(random_state)

    # Enumerate connected subgraphs and classify them by isomorphism class
    motif_counts, n_examined, sampled = _count_motifs(G, motif_size, rng, max_subgraphs)

    # Null model: randomized graphs preserving the (in/out) degree sequence
    # (The old loop wrapped this in ``except Exception: continue``: a failure was skipped, and with every graph
    # skipped the z-scores below came out 0.0. NetworkX can realise any degree sequence of a graph that exists.)
    random_counts_list = []
    for _ in range(n_random):
        graph_seed = int(rng.integers(0, 2**31 - 1))
        if G.is_directed():
            R = nx.DiGraph(
                nx.directed_configuration_model(
                    [d for _, d in G.in_degree()],
                    [d for _, d in G.out_degree()],
                    seed=graph_seed,
                )
            )
        else:
            degree_seq = [d for _, d in G.degree()]
            if sum(degree_seq) % 2 != 0:
                degree_seq[0] += 1
            R = nx.Graph(nx.configuration_model(degree_seq, seed=graph_seed))
        R.remove_edges_from(list(nx.selfloop_edges(R)))
        random_counts_list.append(_count_motifs(R, motif_size, rng, max_subgraphs)[0])

    # Compute z-scores
    z_scores = {}
    for motif_key, observed_count in motif_counts.items():
        random_values = [rc.get(motif_key, 0) for rc in random_counts_list]
        if random_values:
            mean_random = np.mean(random_values)
            std_random = np.std(random_values)
            if std_random > 0:
                z_scores[motif_key] = float((observed_count - mean_random) / std_random)
            elif observed_count == mean_random:
                z_scores[motif_key] = 0.0
            else:
                z_scores[motif_key] = math.copysign(math.inf, observed_count - mean_random)
        else:
            z_scores[motif_key] = 0.0  # n_random == 0: no null model was requested

    # Motif distribution entropy
    total_motifs = sum(motif_counts.values())
    motif_entropy = 0.0
    if total_motifs > 0:
        for count in motif_counts.values():
            p = count / total_motifs
            if p > 0:
                motif_entropy -= p * math.log2(p)

    return {
        "motif_size": motif_size,
        "status": "completed",
        "motifs_found": {k: int(v) for k, v in motif_counts.items()},
        "z_scores": z_scores,
        "motif_entropy": motif_entropy,
        "n_subgraphs_enumerated": total_motifs,
        "n_random_graphs": len(random_counts_list),
        "n_unique_motifs": len(motif_counts),
        "sampled": sampled,
        "n_combinations_examined": n_examined,
    }


def _motif_class(subgraph: Any, nodes: tuple) -> str:
    """Label of the isomorphism class of a small connected subgraph.

    Undirected subgraphs of three or four nodes are identified by their edge count and degree sequence (a complete
    invariant at these sizes). Directed ones are not: the three-node feed-forward loop and the three-cycle share both
    their edge count and their total degrees, so a directed subgraph is labelled by its adjacency matrix, minimised
    over all relabellings of its nodes.
    """
    n_edges = subgraph.number_of_edges()
    if not subgraph.is_directed():
        degree_seq = tuple(sorted(d for _, d in subgraph.degree()))
        return f"n{len(nodes)}_e{n_edges}_d{''.join(str(d) for d in degree_seq)}"

    index = {node: i for i, node in enumerate(nodes)}
    arcs = {(index[u], index[v]) for u, v in subgraph.edges() if u != v}
    k = len(nodes)
    best = min(
        "".join("1" if (order[i], order[j]) in arcs else "0" for i in range(k) for j in range(k))
        for order in permutations(range(k))
    )
    return f"n{k}_e{n_edges}_dir{best}"


def _count_motifs(
    graph: Any, motif_size: int, rng: Optional[np.random.Generator] = None, max_subgraphs: int = 50000
) -> Tuple[Dict[str, int], int, bool]:
    """Count motifs in a graph by enumerating connected subgraphs.

    Args:
        graph: NetworkX graph (directed graphs are classified by their directed structure)
        motif_size: Size of subgraphs to enumerate
        rng: Generator for the sample drawn when the graph has more than ``max_subgraphs`` node combinations
            (a fresh default generator when omitted)
        max_subgraphs: Largest number of node combinations examined

    Returns:
        ``(counts, n_examined, sampled)``: the count of each isomorphism class, the number of node combinations
        examined and whether they are a uniform sample of all combinations rather than all of them.
    """
    motif_counts: Dict[str, int] = {}
    nodes = list(graph.nodes())
    n = len(nodes)

    if n < motif_size:
        return motif_counts, 0, False

    total_combinations = math.comb(n, motif_size)
    sampled = total_combinations > max_subgraphs
    if sampled:
        # Draw distinct random combinations without ever materialising all of them (the old code listed every
        # combination first: 166 million tuples for a 1000-node graph) and from a generator the caller controls (it
        # used a fixed seed, so every call, and every random graph, was examined at the same positions).
        rng = rng if rng is not None else np.random.default_rng()
        chosen: set = set()
        while len(chosen) < max_subgraphs:
            chosen.add(tuple(sorted(rng.choice(n, size=motif_size, replace=False).tolist())))
        node_combos = [tuple(nodes[i] for i in combo) for combo in sorted(chosen)]
    else:
        node_combos = list(combinations(nodes, motif_size))

    directed = graph.is_directed()
    for combo in node_combos:
        subgraph = graph.subgraph(combo)
        if nx.is_connected(subgraph.to_undirected() if directed else subgraph):
            key = _motif_class(subgraph, combo)
            motif_counts[key] = motif_counts.get(key, 0) + 1

    return motif_counts, len(node_combos), sampled


def information_graph_distance(graph1: Any, graph2: Any, method: str = "entropy") -> float:
    """Calculate information-theoretic distance between two graphs.

    Args:
        graph1: First graph
        graph2: Second graph
        method: Distance method ('entropy', 'jensen_shannon')

    Returns:
        Information distance between graphs

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for graph distance calculation")

    # Convert to NetworkX graphs
    if isinstance(graph1, np.ndarray):
        G1 = nx.from_numpy_array(graph1)
    else:
        G1 = graph1

    if isinstance(graph2, np.ndarray):
        G2 = nx.from_numpy_array(graph2)
    else:
        G2 = graph2

    if method == "entropy":
        # Compare network entropies
        ent1 = network_entropy(G1)
        ent2 = network_entropy(G2)
        return abs(ent1 - ent2)

    elif method == "jensen_shannon":
        # Compare degree distributions using Jensen-Shannon divergence
        degrees1 = [d for n, d in G1.degree()]
        degrees2 = [d for n, d in G2.degree()]

        dist1 = Counter(degrees1)
        dist2 = Counter(degrees2)

        # Create common support
        all_degrees = set(dist1.keys()) | set(dist2.keys())
        total1 = sum(dist1.values())
        total2 = sum(dist2.values())

        probs1 = [dist1.get(d, 0) / total1 for d in all_degrees]
        probs2 = [dist2.get(d, 0) / total2 for d in all_degrees]

        # Calculate Jensen-Shannon divergence

        jsd = syntactic.jensen_shannon_divergence(probs1, probs2)

        return jsd

    else:
        raise ValueError(f"Unknown distance method: {method}")
