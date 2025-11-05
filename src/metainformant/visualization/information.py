"""Information theory visualization functions.

This module provides visualization functions for information-theoretic analysis
including entropy plots, mutual information plots, information profiles,
Rényi entropy spectra, and information networks.
"""

from __future__ import annotations

from typing import Any, Sequence

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Use non-interactive backend by default for tests/headless
matplotlib.use("Agg", force=True)


def entropy_plot(
    entropies: Sequence[float],
    *,
    positions: Sequence[int] | None = None,
    ax: plt.Axes | None = None,
    title: str = "Entropy Plot",
    **kwargs
) -> plt.Axes:
    """Plot entropy values across positions or sequences.

    Args:
        entropies: Entropy values
        positions: Optional positions (if None, uses indices)
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import entropy_plot
        >>> import numpy as np
        >>> entropies = np.random.uniform(0, 2, 100)
        >>> ax = entropy_plot(entropies)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    if positions is None:
        positions = np.arange(len(entropies))

    ax.plot(positions, entropies, **kwargs)
    ax.set_xlabel("Position")
    ax.set_ylabel("Entropy (bits)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    return ax


def mutual_information_plot(
    mi_matrix: np.ndarray,
    labels: Sequence[str] | None = None,
    *,
    ax: plt.Axes | None = None,
    title: str = "Mutual Information Matrix",
    **kwargs
) -> plt.Axes:
    """Plot mutual information matrix as heatmap.

    Args:
        mi_matrix: Square matrix of mutual information values
        labels: Optional labels for rows/columns
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import mutual_information_plot
        >>> import numpy as np
        >>> mi_matrix = np.random.random((5, 5))
        >>> mi_matrix = (mi_matrix + mi_matrix.T) / 2  # Symmetric
        >>> np.fill_diagonal(mi_matrix, 1.0)
        >>> ax = mutual_information_plot(mi_matrix)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 8))

    mi_matrix = np.array(mi_matrix)

    im = ax.imshow(mi_matrix, cmap='viridis', aspect='auto', **kwargs)
    
    if labels:
        ax.set_xticks(range(len(labels)))
        ax.set_yticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.set_yticklabels(labels)
    
    plt.colorbar(im, ax=ax, label='Mutual Information (bits)')
    ax.set_title(title)

    return ax


def information_profile_plot(
    profile_data: dict[str, Any],
    *,
    figsize: tuple[int, int] | None = None,
    **kwargs
) -> plt.Figure:
    """Plot information profile visualization.

    Args:
        profile_data: Dictionary with information profile data
        figsize: Figure size tuple
        **kwargs: Additional arguments

    Returns:
        Matplotlib Figure object

    Example:
        >>> from metainformant.visualization import information_profile_plot
        >>> profile = {
        ...     'kmer_frequencies': {'AA': 10, 'AT': 8, 'TA': 7},
        ...     'entropy': 1.5,
        ...     'sequence_complexity': 0.8
        ... }
        >>> fig = information_profile_plot(profile)
    """
    if figsize is None:
        figsize = (12, 10)

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Plot 1: K-mer frequencies (top 20)
    kmer_freqs = profile_data.get('kmer_frequencies', {})
    if kmer_freqs:
        sorted_kmers = sorted(kmer_freqs.items(), key=lambda x: x[1], reverse=True)[:20]
        if sorted_kmers:
            kmers, counts = zip(*sorted_kmers)
            axes[0, 0].bar(range(len(kmers)), counts)
            axes[0, 0].set_xticks(range(len(kmers)))
            axes[0, 0].set_xticklabels(kmers, rotation=45, ha='right')
            axes[0, 0].set_ylabel("Frequency")
            axes[0, 0].set_title("Top 20 K-mers")

    # Plot 2: Entropy summary
    entropy = profile_data.get('entropy', 0.0)
    complexity = profile_data.get('sequence_complexity', 0.0)
    axes[0, 1].bar(['Entropy', 'Complexity'], [entropy, complexity])
    axes[0, 1].set_ylabel("Value")
    axes[0, 1].set_title("Information Measures")

    # Plot 3: K-mer distribution
    if kmer_freqs:
        counts_list = list(kmer_freqs.values())
        axes[1, 0].hist(counts_list, bins=30, edgecolor='black', alpha=0.7)
        axes[1, 0].set_xlabel("K-mer Frequency")
        axes[1, 0].set_ylabel("Count")
        axes[1, 0].set_title("K-mer Frequency Distribution")

    # Plot 4: Summary statistics
    axes[1, 1].axis('off')
    stats_text = f"""
    Entropy: {entropy:.3f} bits
    Unique K-mers: {profile_data.get('unique_kmers', 0)}
    Complexity: {complexity:.3f}
    """
    axes[1, 1].text(0.1, 0.5, stats_text, fontsize=12, verticalalignment='center')

    fig.suptitle("Information Profile", fontsize=14)
    plt.tight_layout()

    return fig


def renyi_spectrum_plot(
    alpha_range: Sequence[float],
    entropies: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    title: str = "Rényi Entropy Spectrum",
    **kwargs
) -> plt.Axes:
    """Plot Rényi entropy as a function of order α.

    Args:
        alpha_range: Range of α values
        entropies: Rényi entropy values for each α
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import renyi_spectrum_plot
        >>> import numpy as np
        >>> alpha = [0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]
        >>> entropies = [2.0, 1.8, 1.5, 1.3, 1.2, 1.0, 0.8, 0.6]
        >>> ax = renyi_spectrum_plot(alpha, entropies)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    ax.plot(alpha_range, entropies, marker='o', linewidth=2, markersize=6, **kwargs)
    ax.axvline(x=1.0, color='r', linestyle='--', alpha=0.5, label='α=1 (Shannon)')
    ax.set_xlabel("Order α")
    ax.set_ylabel("Rényi Entropy (bits)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()

    return ax


def information_network_plot(
    mi_matrix: np.ndarray,
    labels: Sequence[str] | None = None,
    *,
    threshold: float = 0.1,
    ax: plt.Axes | None = None,
    title: str = "Information Network",
    **kwargs
) -> plt.Axes:
    """Plot network visualization based on mutual information matrix.

    Args:
        mi_matrix: Mutual information matrix
        labels: Optional node labels
        threshold: Minimum MI threshold for edges
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import information_network_plot
        >>> import numpy as np
        >>> mi_matrix = np.random.random((5, 5))
        >>> mi_matrix = (mi_matrix + mi_matrix.T) / 2
        >>> np.fill_diagonal(mi_matrix, 0)
        >>> labels = ['A', 'B', 'C', 'D', 'E']
        >>> ax = information_network_plot(mi_matrix, labels, threshold=0.3)
    """
    try:
        import networkx as nx
    except ImportError:
        if ax is None:
            _, ax = plt.subplots()
        ax.text(0.5, 0.5, "NetworkX not available\nInstall with: uv pip install networkx",
               ha="center", va="center", transform=ax.transAxes)
        return ax

    if ax is None:
        _, ax = plt.subplots(figsize=(12, 10))

    mi_matrix = np.array(mi_matrix)
    n = len(mi_matrix)

    # Create network from MI matrix
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

    # Layout
    pos = nx.spring_layout(G, k=1, iterations=50)
    
    # Draw network
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color='lightblue', node_size=500, alpha=0.8)
    
    # Draw edges with weights
    edges = G.edges()
    edge_weights = [G[u][v]['weight'] for u, v in edges]
    nx.draw_networkx_edges(G, pos, ax=ax, width=[w * 2 for w in edge_weights],
                          alpha=0.5, edge_color='gray')
    
    if labels:
        labels_dict = {node: node for node in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels_dict, ax=ax, font_size=8)
    
    ax.set_title(title)
    ax.axis('off')

    return ax

