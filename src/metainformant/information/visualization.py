"""Information theory visualization integration.

This module provides visualization functions for information-theoretic
analysis, integrating with the metainformant.visualization module.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

try:
    from metainformant.visualization import plots
except ImportError:
    plots = None


def plot_entropy_distribution(
    entropies: list[float] | np.ndarray,
    output_path: Path | str | None = None,
    title: str = "Entropy Distribution"
) -> dict[str, Any]:
    """Plot distribution of entropy values.
    
    Args:
        entropies: List of entropy values
        output_path: Path to save plot (if None, uses output/information/)
        title: Plot title
        
    Returns:
        Dictionary with plot metadata
    """
    if plots is None:
        raise ImportError("visualization.plots module not available")
    
    if output_path is None:
        from metainformant.core.paths import ensure_directory
        output_dir = Path("output/information")
        ensure_directory(output_dir)
        output_path = output_dir / "entropy_distribution.png"
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Create histogram
    fig, ax = plots.create_figure()
    ax.hist(entropies, bins=30, edgecolor="black", alpha=0.7)
    ax.set_xlabel("Entropy (bits)")
    ax.set_ylabel("Frequency")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    
    plots.save_figure(fig, output_path)
    
    return {
        "output_path": str(output_path),
        "num_values": len(entropies),
        "mean_entropy": float(np.mean(entropies)),
        "std_entropy": float(np.std(entropies)),
    }


def plot_mutual_information_matrix(
    mi_matrix: np.ndarray,
    labels: list[str] | None = None,
    output_path: Path | str | None = None,
    title: str = "Mutual Information Matrix"
) -> dict[str, Any]:
    """Plot mutual information matrix as heatmap.
    
    Args:
        mi_matrix: Square matrix of mutual information values
        labels: Optional labels for rows/columns
        output_path: Path to save plot
        title: Plot title
        
    Returns:
        Dictionary with plot metadata
    """
    if plots is None:
        raise ImportError("visualization.plots module not available")
    
    if output_path is None:
        from metainformant.core.paths import ensure_directory
        output_dir = Path("output/information")
        ensure_directory(output_dir)
        output_path = output_dir / "mi_matrix.png"
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Create heatmap
    fig, ax = plots.create_figure()
    im = ax.imshow(mi_matrix, cmap="viridis", aspect="auto")
    ax.set_title(title)
    
    if labels:
        ax.set_xticks(range(len(labels)))
        ax.set_yticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_yticklabels(labels)
    
    fig.colorbar(im, ax=ax, label="Mutual Information (bits)")
    
    plots.save_figure(fig, output_path)
    
    return {
        "output_path": str(output_path),
        "matrix_shape": mi_matrix.shape,
        "max_mi": float(np.max(mi_matrix)),
        "mean_mi": float(np.mean(mi_matrix)),
    }


def plot_information_profile(
    profile: dict[str, Any],
    output_path: Path | str | None = None,
    title: str = "Information Profile"
) -> dict[str, Any]:
    """Plot information profile visualization.
    
    Args:
        profile: Information profile dictionary from information_profile()
        output_path: Path to save plot
        title: Plot title
        
    Returns:
        Dictionary with plot metadata
    """
    if plots is None:
        raise ImportError("visualization.plots module not available")
    
    if output_path is None:
        from metainformant.core.paths import ensure_directory
        output_dir = Path("output/information")
        ensure_directory(output_dir)
        output_path = output_dir / "information_profile.png"
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Create multi-panel plot
    fig, axes = plots.create_figure(nrows=2, ncols=2, figsize=(12, 10))
    
    # Plot 1: K-mer frequencies (top 20)
    kmer_freqs = profile.get("kmer_frequencies", {})
    if kmer_freqs:
        sorted_kmers = sorted(kmer_freqs.items(), key=lambda x: x[1], reverse=True)[:20]
        kmers, counts = zip(*sorted_kmers) if sorted_kmers else ([], [])
        axes[0, 0].bar(range(len(kmers)), counts)
        axes[0, 0].set_xticks(range(len(kmers)))
        axes[0, 0].set_xticklabels(kmers, rotation=45, ha="right")
        axes[0, 0].set_ylabel("Frequency")
        axes[0, 0].set_title("Top 20 K-mers")
    
    # Plot 2: Entropy summary
    entropy = profile.get("entropy", 0.0)
    complexity = profile.get("sequence_complexity", 0.0)
    axes[0, 1].bar(["Entropy", "Complexity"], [entropy, complexity])
    axes[0, 1].set_ylabel("Value")
    axes[0, 1].set_title("Information Measures")
    
    # Plot 3: K-mer distribution
    if kmer_freqs:
        counts_list = list(kmer_freqs.values())
        axes[1, 0].hist(counts_list, bins=30, edgecolor="black", alpha=0.7)
        axes[1, 0].set_xlabel("K-mer Frequency")
        axes[1, 0].set_ylabel("Count")
        axes[1, 0].set_title("K-mer Frequency Distribution")
    
    # Plot 4: Summary statistics
    axes[1, 1].axis("off")
    stats_text = f"""
    Entropy: {entropy:.3f} bits
    Unique K-mers: {profile.get('unique_kmers', 0)}
    Complexity: {complexity:.3f}
    """
    axes[1, 1].text(0.1, 0.5, stats_text, fontsize=12, verticalalignment="center")
    
    fig.suptitle(title, fontsize=14)
    plots.save_figure(fig, output_path)
    
    return {
        "output_path": str(output_path),
        "profile": profile,
    }


def plot_renyi_spectrum(
    probs: list[float],
    alpha_range: list[float] | None = None,
    output_path: Path | str | None = None,
    title: str = "Rényi Entropy Spectrum"
) -> dict[str, Any]:
    """Plot Rényi entropy as a function of order α.
    
    Args:
        probs: Probability distribution
        alpha_range: Range of α values to plot (default: [0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0])
        output_path: Path to save plot
        title: Plot title
        
    Returns:
        Dictionary with plot metadata
    """
    if plots is None:
        raise ImportError("visualization.plots module not available")
    
    from metainformant.information.syntactic import renyi_entropy
    
    if alpha_range is None:
        alpha_range = [0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]
    
    if output_path is None:
        from metainformant.core.paths import ensure_directory
        output_dir = Path("output/information")
        ensure_directory(output_dir)
        output_path = output_dir / "renyi_spectrum.png"
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Calculate Rényi entropy for each α
    entropies = []
    for alpha in alpha_range:
        try:
            h_alpha = renyi_entropy(probs, alpha=alpha)
            entropies.append(h_alpha)
        except Exception:
            entropies.append(0.0)
    
    # Create plot
    fig, ax = plots.create_figure()
    ax.plot(alpha_range, entropies, marker="o", linewidth=2, markersize=6)
    ax.set_xlabel("Order α")
    ax.set_ylabel("Rényi Entropy (bits)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.axvline(x=1.0, color="r", linestyle="--", alpha=0.5, label="α=1 (Shannon)")
    ax.legend()
    
    plots.save_figure(fig, output_path)
    
    return {
        "output_path": str(output_path),
        "alpha_range": alpha_range,
        "entropies": entropies,
    }


def plot_information_network(
    network: Any,
    mi_matrix: np.ndarray | None = None,
    node_labels: list[str] | None = None,
    output_path: Path | str | None = None,
    title: str = "Information Network"
) -> dict[str, Any]:
    """Visualize information flow in a network using MI matrix.
    
    Args:
        network: NetworkX graph or network object
        mi_matrix: Optional mutual information matrix
        node_labels: Optional node labels
        output_path: Path to save plot
        title: Plot title
        
    Returns:
        Dictionary with plot metadata
    """
    if plots is None:
        raise ImportError("visualization.plots module not available")
    
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("NetworkX not available for network visualization")
    
    if output_path is None:
        from metainformant.core.paths import ensure_directory
        output_dir = Path("output/information")
        ensure_directory(output_dir)
        output_path = output_dir / "information_network.png"
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Convert network to NetworkX if needed
    if not isinstance(network, nx.Graph):
        if hasattr(network, "network"):
            G = network.network
        else:
            G = nx.Graph()
            # Try to create from edges
            if hasattr(network, "edges"):
                G.add_edges_from(network.edges())
    else:
        G = network
    
    # Create plot
    fig, ax = plots.create_figure(figsize=(12, 10))
    
    # Use MI matrix for edge weights if provided
    if mi_matrix is not None and len(mi_matrix) == len(G.nodes()):
        nodes = list(G.nodes())
        for i, node_i in enumerate(nodes):
            for j, node_j in enumerate(nodes):
                if i != j and mi_matrix[i, j] > 0:
                    if not G.has_edge(node_i, node_j):
                        G.add_edge(node_i, node_j, weight=mi_matrix[i, j])
    
    # Layout
    pos = nx.spring_layout(G, k=1, iterations=50)
    
    # Draw network
    nx.draw(
        G,
        pos,
        ax=ax,
        with_labels=node_labels is not None,
        node_color="lightblue",
        node_size=500,
        font_size=8,
        edge_color="gray",
        alpha=0.6,
    )
    
    if node_labels:
        labels = {node: node_labels[i] if i < len(node_labels) else str(node) for i, node in enumerate(G.nodes())}
        nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=8)
    
    ax.set_title(title)
    plots.save_figure(fig, output_path)
    
    return {
        "output_path": str(output_path),
        "num_nodes": G.number_of_nodes(),
        "num_edges": G.number_of_edges(),
    }


def plot_entropy_landscape(
    entropy_data: np.ndarray,
    x_labels: list[str] | None = None,
    y_labels: list[str] | None = None,
    output_path: Path | str | None = None,
    title: str = "Entropy Landscape"
) -> dict[str, Any]:
    """Plot 2D/3D visualization of entropy landscape.
    
    Args:
        entropy_data: 2D array of entropy values
        x_labels: Optional labels for x-axis
        y_labels: Optional labels for y-axis
        output_path: Path to save plot
        title: Plot title
        
    Returns:
        Dictionary with plot metadata
    """
    if plots is None:
        raise ImportError("visualization.plots module not available")
    
    if output_path is None:
        from metainformant.core.paths import ensure_directory
        output_dir = Path("output/information")
        ensure_directory(output_dir)
        output_path = output_dir / "entropy_landscape.png"
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    entropy_data = np.asarray(entropy_data)
    
    # Create 2D heatmap
    fig, ax = plots.create_figure()
    im = ax.imshow(entropy_data, cmap="viridis", aspect="auto", interpolation="nearest")
    
    if x_labels:
        ax.set_xticks(range(len(x_labels)))
        ax.set_xticklabels(x_labels, rotation=45, ha="right")
    if y_labels:
        ax.set_yticks(range(len(y_labels)))
        ax.set_yticklabels(y_labels)
    
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="Entropy (bits)")
    
    plots.save_figure(fig, output_path)
    
    return {
        "output_path": str(output_path),
        "data_shape": entropy_data.shape,
        "min_entropy": float(np.min(entropy_data)),
        "max_entropy": float(np.max(entropy_data)),
        "mean_entropy": float(np.mean(entropy_data)),
    }


def plot_mi_network(
    mi_matrix: np.ndarray,
    labels: list[str] | None = None,
    threshold: float = 0.1,
    output_path: Path | str | None = None,
    title: str = "Mutual Information Network"
) -> dict[str, Any]:
    """Plot network visualization colored by mutual information.
    
    Args:
        mi_matrix: Mutual information matrix
        labels: Optional node labels
        threshold: Minimum MI threshold for edges
        output_path: Path to save plot
        title: Plot title
        
    Returns:
        Dictionary with plot metadata
    """
    if plots is None:
        raise ImportError("visualization.plots module not available")
    
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("NetworkX not available for network visualization")
    
    if output_path is None:
        from metainformant.core.paths import ensure_directory
        output_dir = Path("output/information")
        ensure_directory(output_dir)
        output_path = output_dir / "mi_network.png"
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Create network from MI matrix
    n = len(mi_matrix)
    G = nx.Graph()
    
    # Add nodes
    if labels:
        G.add_nodes_from(labels)
    else:
        G.add_nodes_from(range(n))
    
    # Add edges above threshold
    for i in range(n):
        for j in range(i + 1, n):
            if mi_matrix[i, j] >= threshold:
                node_i = labels[i] if labels and i < len(labels) else i
                node_j = labels[j] if labels and j < len(labels) else j
                G.add_edge(node_i, node_j, weight=mi_matrix[i, j])
    
    # Create plot
    fig, ax = plots.create_figure(figsize=(12, 10))
    
    # Layout
    pos = nx.spring_layout(G, k=1, iterations=50)
    
    # Get edge weights for coloring
    edges = G.edges()
    edge_weights = [G[u][v]["weight"] for u, v in edges]
    
    # Draw network
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color="lightblue", node_size=500, alpha=0.8)
    nx.draw_networkx_edges(
        G, pos, ax=ax, width=[w * 2 for w in edge_weights], alpha=0.5, edge_color="gray"
    )
    
    if labels:
        labels_dict = {node: node for node in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels_dict, ax=ax, font_size=8)
    
    ax.set_title(title)
    ax.axis("off")
    
    plots.save_figure(fig, output_path)
    
    return {
        "output_path": str(output_path),
        "num_nodes": G.number_of_nodes(),
        "num_edges": G.number_of_edges(),
        "threshold": threshold,
    }


def plot_semantic_similarity_network(
    similarity_matrix: np.ndarray,
    terms: list[str],
    output_path: Path | str | None = None,
    title: str = "Semantic Similarity Network"
) -> dict[str, Any]:
    """Plot hierarchical network visualization of semantic similarity.
    
    Args:
        similarity_matrix: Semantic similarity matrix
        terms: List of terms
        output_path: Path to save plot
        title: Plot title
        
    Returns:
        Dictionary with plot metadata
    """
    if plots is None:
        raise ImportError("visualization.plots module not available")
    
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("NetworkX not available for network visualization")
    
    if output_path is None:
        from metainformant.core.paths import ensure_directory
        output_dir = Path("output/information")
        ensure_directory(output_dir)
        output_path = output_dir / "semantic_similarity_network.png"
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Create network from similarity matrix
    G = nx.Graph()
    G.add_nodes_from(terms)
    
    # Add edges for high similarity
    threshold = 0.3  # Adjust threshold as needed
    for i, term_i in enumerate(terms):
        for j, term_j in enumerate(terms):
            if i < j and similarity_matrix[i, j] >= threshold:
                G.add_edge(term_i, term_j, weight=similarity_matrix[i, j])
    
    # Create plot
    fig, ax = plots.create_figure(figsize=(12, 10))
    
    # Use hierarchical layout
    try:
        pos = nx.spring_layout(G, k=2, iterations=100)
    except Exception:
        pos = nx.circular_layout(G)
    
    # Draw network
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color="lightcoral", node_size=800, alpha=0.8)
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.3, edge_color="gray")
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8)
    
    ax.set_title(title)
    ax.axis("off")
    
    plots.save_figure(fig, output_path)
    
    return {
        "output_path": str(output_path),
        "num_nodes": G.number_of_nodes(),
        "num_edges": G.number_of_edges(),
    }

