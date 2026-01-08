"""Visualization functions for life events analysis.

This module provides plotting and visualization utilities for life event sequences,
embeddings, and predictive models.
"""

from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List, Optional, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Optional dependencies
try:
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    logger.warning("matplotlib not available, life events visualization disabled")

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    sns = None


def plot_event_timeline(
    sequence: Any,
    output_path: Optional[str] = None,
    figsize: tuple[int, int] = (12, 6),
    show_domains: bool = True,
    title: str = "Life Event Timeline"
) -> Optional[Any]:
    """Plot a timeline of life events.

    Args:
        sequence: EventSequence to plot
        output_path: Optional path to save the plot
        figsize: Figure size (width, height)
        show_domains: Whether to show event domains
        title: Plot title

    Returns:
        Matplotlib figure if matplotlib available, None otherwise
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for timeline plotting")
        return None

    fig, ax = plt.subplots(figsize=figsize)

    # Extract events
    events = sequence.events
    if not events:
        logger.warning("No events to plot")
        return None

    # Sort events by date
    events = sorted(events, key=lambda x: x.timestamp)

    # Plot events
    dates = [event.timestamp for event in events]
    event_types = [event.event_type for event in events]

    # Create scatter plot with different colors for domains
    domains = [event.domain for event in events] if hasattr(events[0], 'domain') else ['default'] * len(events)
    unique_domains = list(set(domains))

    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_domains)))
    domain_colors = dict(zip(unique_domains, colors))

    for i, event in enumerate(events):
        color = domain_colors.get(event.domain if hasattr(event, 'domain') else 'default', 'blue')
        ax.scatter(event.timestamp, i, c=[color], s=100, alpha=0.7, edgecolors='black', linewidth=1)

        # Add event type label
        ax.annotate(event.event_type, (event.timestamp, i),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=8, ha='left', va='bottom',
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

    # Format x-axis as dates
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.YearLocator())
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

    ax.set_xlabel('Date')
    ax.set_ylabel('Event Order')
    ax.set_title(f'{title} - {sequence.person_id}')
    ax.grid(True, alpha=0.3)

    # Add domain legend if requested
    if show_domains and len(unique_domains) > 1:
        legend_elements = [plt.Line2D([0], [0], marker='o', color='w',
                                     markerfacecolor=color, markersize=10,
                                     label=domain)
                          for domain, color in domain_colors.items()]
        ax.legend(handles=legend_elements, title='Domains', loc='upper right')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved timeline plot to {output_path}")

    return fig


def plot_event_embeddings(
    embeddings: Dict[str, Any],
    output_path: Optional[str] = None,
    method: str = "pca",
    figsize: tuple[int, int] = (10, 8),
    title: str = "Event Sequence Embeddings"
) -> Optional[Any]:
    """Plot event sequence embeddings.

    Args:
        embeddings: Dictionary with embeddings and metadata
        output_path: Optional path to save the plot
        method: Dimensionality reduction method ('pca', 'tsne', 'umap')
        figsize: Figure size (width, height)
        title: Plot title

    Returns:
        Matplotlib figure if matplotlib available, None otherwise
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for embedding plotting")
        return None

    if 'embeddings' not in embeddings:
        logger.error("No embeddings data provided")
        return None

    embedding_matrix = np.array(embeddings['embeddings'])
    if embedding_matrix.shape[1] > 2:
        # Reduce to 2D if needed
        if method.lower() == "pca":
            from sklearn.decomposition import PCA
            reducer = PCA(n_components=2)
        elif method.lower() == "tsne":
            from sklearn.manifold import TSNE
            reducer = TSNE(n_components=2, random_state=42)
        elif method.lower() == "umap":
            try:
                import umap
                reducer = umap.UMAP(n_components=2, random_state=42)
            except ImportError:
                logger.warning("UMAP not available, falling back to PCA")
                from sklearn.decomposition import PCA
                reducer = PCA(n_components=2)
        else:
            from sklearn.decomposition import PCA
            reducer = PCA(n_components=2)

        embedding_2d = reducer.fit_transform(embedding_matrix)
    else:
        embedding_2d = embedding_matrix

    fig, ax = plt.subplots(figsize=figsize)

    # Plot embeddings
    scatter = ax.scatter(embedding_2d[:, 0], embedding_2d[:, 1],
                        c=range(len(embedding_2d)), cmap='viridis', alpha=0.7)

    ax.set_xlabel(f'{method.upper()} Component 1')
    ax.set_ylabel(f'{method.upper()} Component 2')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Sequence Index')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved embedding plot to {output_path}")

    return fig


def plot_attention_heatmap(
    attention_weights: Any,
    output_path: Optional[str] = None,
    figsize: tuple[int, int] = (12, 8),
    title: str = "Attention Weights Heatmap"
) -> Optional[Any]:
    """Plot attention weights as a heatmap.

    Args:
        attention_weights: Attention weights matrix or tensor
        output_path: Optional path to save the plot
        figsize: Figure size (width, height)
        title: Plot title

    Returns:
        Matplotlib figure if matplotlib available, None otherwise
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for attention heatmap")
        return None

    if not HAS_SEABORN:
        from metainformant.core.utils.optional_deps import warn_optional_dependency
        warn_optional_dependency("seaborn", "attention heatmap visualization")
        return None

    # Convert attention weights to numpy array
    if hasattr(attention_weights, 'detach'):
        # PyTorch tensor
        weights = attention_weights.detach().cpu().numpy()
    elif hasattr(attention_weights, 'numpy'):
        # NumPy array or similar
        weights = np.array(attention_weights)
    else:
        weights = np.array(attention_weights)

    # Handle different shapes
    if weights.ndim == 3:
        # Batch x Seq x Seq - take first batch
        weights = weights[0]
    elif weights.ndim == 2:
        # Seq x Seq
        pass
    else:
        logger.error(f"Unsupported attention weights shape: {weights.shape}")
        return None

    fig, ax = plt.subplots(figsize=figsize)

    # Plot heatmap
    sns.heatmap(weights, cmap='YlOrRd', ax=ax, square=True,
                cbar_kws={'label': 'Attention Weight'})

    ax.set_xlabel('Target Position')
    ax.set_ylabel('Source Position')
    ax.set_title(title)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved attention heatmap to {output_path}")

    return fig


def plot_prediction_importance(
    feature_importance: Dict[str, float],
    output_path: Optional[str] = None,
    figsize: tuple[int, int] = (10, 6),
    title: str = "Feature Importance",
    top_n: Optional[int] = None
) -> Optional[Any]:
    """Plot feature importance for predictions.

    Args:
        feature_importance: Dictionary mapping feature names to importance scores
        output_path: Optional path to save the plot
        figsize: Figure size (width, height)
        title: Plot title
        top_n: Number of top features to show (None for all)

    Returns:
        Matplotlib figure if matplotlib available, None otherwise
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for importance plotting")
        return None

    # Sort features by importance
    sorted_features = sorted(feature_importance.items(), key=lambda x: x[1], reverse=True)

    if top_n:
        sorted_features = sorted_features[:top_n]

    features, importances = zip(*sorted_features)

    fig, ax = plt.subplots(figsize=figsize)

    # Create horizontal bar plot
    bars = ax.barh(range(len(features)), importances, color='skyblue', alpha=0.8)

    ax.set_yticks(range(len(features)))
    ax.set_yticklabels(features)
    ax.set_xlabel('Importance Score')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    # Add value labels on bars
    for i, (bar, importance) in enumerate(zip(bars, importances)):
        ax.text(importance + 0.01, i, f'{importance:.3f}',
               va='center', fontsize=8)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved importance plot to {output_path}")

    return fig


def plot_domain_distribution(
    sequences: List[Any],
    output_path: Optional[str] = None,
    figsize: tuple[int, int] = (12, 6),
    title: str = "Event Domain Distribution"
) -> Optional[Any]:
    """Plot distribution of event domains across sequences.

    Args:
        sequences: List of EventSequence objects
        output_path: Optional path to save the plot
        figsize: Figure size (width, height)
        title: Plot title

    Returns:
        Matplotlib figure if matplotlib available, None otherwise
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for domain distribution plotting")
        return None

    # Count domains across all sequences
    domain_counts = {}
    for sequence in sequences:
        for event in sequence.events:
            if hasattr(event, 'domain'):
                domain = event.domain
                domain_counts[domain] = domain_counts.get(domain, 0) + 1

    if not domain_counts:
        logger.warning("No domain information found in sequences")
        return None

    domains, counts = zip(*domain_counts.items())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Bar plot
    bars = ax1.bar(range(len(domains)), counts, color='lightcoral', alpha=0.7)
    ax1.set_xticks(range(len(domains)))
    ax1.set_xticklabels(domains, rotation=45, ha='right')
    ax1.set_ylabel('Count')
    ax1.set_title('Domain Counts')
    ax1.grid(True, alpha=0.3)

    # Pie chart
    ax2.pie(counts, labels=domains, autopct='%1.1f%%', startangle=90)
    ax2.set_title('Domain Proportions')
    ax2.axis('equal')

    fig.suptitle(title, fontsize=14)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved domain distribution plot to {output_path}")

    return fig


def plot_intervention_effects(
    pre_sequences: List[EventSequence],
    post_sequences: List[EventSequence],
    pre_outcomes: Optional[np.ndarray] = None,
    post_outcomes: Optional[np.ndarray] = None,
    output_path: Optional[Union[str, Path]] = None,
    figsize: Tuple[int, int] = (12, 6)
) -> Any:
    """Plot intervention effects on life event sequences and outcomes.

    Args:
        pre_sequences: Sequences before intervention
        post_sequences: Sequences after intervention
        pre_outcomes: Outcomes before intervention
        post_outcomes: Outcomes after intervention
        output_path: Path to save the plot (optional)
        figsize: Figure size (width, height)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create intervention effects plot")
        return None

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Plot 1: Sequence length comparison
    pre_lengths = [len(seq.events) for seq in pre_sequences]
    post_lengths = [len(seq.events) for seq in post_sequences]

    ax1.boxplot([pre_lengths, post_lengths], labels=['Pre-intervention', 'Post-intervention'])
    ax1.set_ylabel('Number of Events')
    ax1.set_title('Sequence Length Changes')
    ax1.grid(True, alpha=0.3)

    # Plot 2: Outcome comparison (if provided)
    if pre_outcomes is not None and post_outcomes is not None:
        ax2.scatter(pre_outcomes, post_outcomes, alpha=0.7, s=50)
        ax2.plot([min(pre_outcomes), max(pre_outcomes)], [min(pre_outcomes), max(pre_outcomes)],
                'r--', alpha=0.7, label='No change')
        ax2.set_xlabel('Pre-intervention Outcomes')
        ax2.set_ylabel('Post-intervention Outcomes')
        ax2.set_title('Outcome Changes')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        # Add correlation coefficient
        if len(pre_outcomes) > 1:
            corr = np.corrcoef(pre_outcomes, post_outcomes)[0, 1]
            ax2.text(0.05, 0.95, f'Correlation: {corr:.3f}',
                    transform=ax2.transAxes, fontsize=10,
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    else:
        ax2.text(0.5, 0.5, 'No outcome data provided',
                transform=ax2.transAxes, ha='center', va='center', fontsize=12)
        ax2.set_title('Outcome Comparison (No Data)')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved intervention effects plot to {output_path}")

    return fig


def plot_outcome_distribution(
    outcomes: np.ndarray,
    output_path: Optional[Union[str, Path]] = None,
    plot_type: str = "histogram",
    figsize: Tuple[int, int] = (10, 6)
) -> Any:
    """Plot distribution of outcomes.

    Args:
        outcomes: Array of outcome values
        output_path: Path to save the plot (optional)
        plot_type: Type of plot ('histogram', 'boxplot', 'violin')
        figsize: Figure size (width, height)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create outcome distribution plot")
        return None

    fig, ax = plt.subplots(figsize=figsize)

    if plot_type == "histogram":
        ax.hist(outcomes, bins=20, alpha=0.7, edgecolor='black')
        ax.set_xlabel('Outcome Value')
        ax.set_ylabel('Frequency')
        ax.set_title('Outcome Distribution Histogram')

    elif plot_type == "boxplot":
        ax.boxplot(outcomes)
        ax.set_ylabel('Outcome Value')
        ax.set_title('Outcome Distribution Boxplot')
        ax.set_xticks([1])
        ax.set_xticklabels(['Outcomes'])

    elif plot_type == "violin":
        # Create violin plot manually since seaborn might not be available
        if len(outcomes) > 0:
            # Simple approximation of violin plot
            parts = ax.violinplot(outcomes, showmeans=True, showmedians=True)
            ax.set_ylabel('Outcome Value')
            ax.set_title('Outcome Distribution Violin Plot')
            ax.set_xticks([1])
            ax.set_xticklabels(['Outcomes'])

    else:
        logger.warning(f"Unknown plot type: {plot_type}, using histogram")
        ax.hist(outcomes, bins=20, alpha=0.7, edgecolor='black')
        ax.set_xlabel('Outcome Value')
        ax.set_ylabel('Frequency')
        ax.set_title('Outcome Distribution Histogram')

    # Add statistics
    if len(outcomes) > 0:
        mean_val = np.mean(outcomes)
        median_val = np.median(outcomes)
        std_val = np.std(outcomes)

        stats_text = f'Mean: {mean_val:.3f}\nMedian: {median_val:.3f}\nStd: {std_val:.3f}'
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved outcome distribution plot to {output_path}")

    return fig


def plot_event_frequency_heatmap(
    sequences: List[EventSequence],
    output_path: Optional[Union[str, Path]] = None,
    time_bins: int = 10,
    figsize: Tuple[int, int] = (12, 8)
) -> Any:
    """Plot event frequency heatmap over time.

    Args:
        sequences: List of EventSequence objects
        output_path: Path to save the plot (optional)
        time_bins: Number of time bins to divide the timeline into
        figsize: Figure size (width, height)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create frequency heatmap")
        return None

    # Find global time range
    all_timestamps = []
    for seq in sequences:
        for event in seq.events:
            all_timestamps.append(event.timestamp)

    if not all_timestamps:
        logger.warning("No timestamps found in sequences")
        return None

    min_time = min(all_timestamps)
    max_time = max(all_timestamps)
    time_range = max_time - min_time

    if time_range.total_seconds() == 0:
        logger.warning("All events have the same timestamp")
        return None

    # Create time bins
    bin_edges = [min_time + i * time_range / time_bins for i in range(time_bins + 1)]

    # Collect event frequencies by type and time bin
    event_types = set()
    for seq in sequences:
        for event in seq.events:
            event_types.add(event.event_type)

    event_types = sorted(list(event_types))

    # Initialize frequency matrix
    frequency_matrix = np.zeros((len(event_types), time_bins))

    # Count events in each bin
    for seq in sequences:
        for event in seq.events:
            if event.event_type in event_types:
                # Find which bin this event belongs to
                for i in range(time_bins):
                    if bin_edges[i] <= event.timestamp < bin_edges[i + 1]:
                        event_idx = event_types.index(event.event_type)
                        frequency_matrix[event_idx, i] += 1
                        break
                # Handle the last bin (include upper bound)
                if event.timestamp == max_time:
                    event_idx = event_types.index(event.event_type)
                    frequency_matrix[event_idx, time_bins - 1] += 1

    # Plot heatmap
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(frequency_matrix, cmap='YlOrRd', aspect='auto')

    # Set labels
    ax.set_xticks(range(time_bins))
    ax.set_yticks(range(len(event_types)))
    ax.set_xticklabels([f'Bin {i+1}' for i in range(time_bins)], rotation=45, ha='right')
    ax.set_yticklabels([event.split(':')[-1] for event in event_types])

    ax.set_title('Event Frequency Heatmap Over Time')
    ax.set_xlabel('Time Bin')
    ax.set_ylabel('Event Type')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Event Count')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved event frequency heatmap to {output_path}")

    return fig


def plot_event_cooccurrence(
    sequences: List[EventSequence],
    output_path: Optional[Union[str, Path]] = None,
    top_n: int = 20,
    figsize: Tuple[int, int] = (12, 10)
) -> Any:
    """Plot event co-occurrence matrix.

    Args:
        sequences: List of EventSequence objects
        output_path: Path to save the plot (optional)
        top_n: Number of top events to show
        figsize: Figure size (width, height)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create co-occurrence plot")
        return None

    # Calculate co-occurrence matrix
    event_counts = {}
    cooccurrence = {}

    for seq in sequences:
        events_in_seq = set(event.event_type for event in seq.events)

        for event1 in events_in_seq:
            if event1 not in event_counts:
                event_counts[event1] = 0
                cooccurrence[event1] = {}

            event_counts[event1] += 1

            for event2 in events_in_seq:
                if event1 != event2:
                    if event2 not in cooccurrence[event1]:
                        cooccurrence[event1][event2] = 0
                    cooccurrence[event1][event2] += 1

    # Get top N events by frequency
    sorted_events = sorted(event_counts.items(), key=lambda x: x[1], reverse=True)
    top_events = [event for event, count in sorted_events[:top_n]]

    # Create co-occurrence matrix for top events
    n_events = len(top_events)
    matrix = np.zeros((n_events, n_events))

    for i, event1 in enumerate(top_events):
        for j, event2 in enumerate(top_events):
            if i != j and event1 in cooccurrence and event2 in cooccurrence[event1]:
                # Normalize by geometric mean of individual frequencies
                freq1 = event_counts[event1]
                freq2 = event_counts[event2]
                expected = (freq1 * freq2) / len(sequences)
                if expected > 0:
                    matrix[i, j] = cooccurrence[event1][event2] / expected

    # Plot heatmap
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')

    # Set labels
    ax.set_xticks(range(n_events))
    ax.set_yticks(range(n_events))
    ax.set_xticklabels([event.split(':')[-1] for event in top_events], rotation=45, ha='right')
    ax.set_yticklabels([event.split(':')[-1] for event in top_events])

    ax.set_title('Event Co-occurrence Matrix')
    ax.set_xlabel('Event 2')
    ax.set_ylabel('Event 1')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Observed/Expected Ratio')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved event co-occurrence plot to {output_path}")

    return fig


def plot_embedding_clusters(
    embeddings: Dict[str, np.ndarray],
    clusters: Optional[Dict[str, int]] = None,
    output_path: Optional[Union[str, Path]] = None,
    figsize: Tuple[int, int] = (12, 8)
) -> Any:
    """Plot embedding clusters for life events.

    Args:
        embeddings: Dictionary mapping event types to embedding vectors
        clusters: Optional dictionary mapping event types to cluster IDs
        output_path: Path to save the plot (optional)
        figsize: Figure size (width, height)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create embedding clusters plot")
        return None

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Plot 1: Embedding distributions by cluster
    ax1 = axes[0]
    if clusters:
        cluster_data = {}
        for event_type, embedding in embeddings.items():
            cluster_id = clusters.get(event_type, -1)
            if cluster_id not in cluster_data:
                cluster_data[cluster_id] = []
            cluster_data[cluster_id].append(np.mean(embedding))  # Use mean for simplicity

        colors = plt.cm.tab10(np.linspace(0, 1, len(cluster_data)))
        for i, (cluster_id, values) in enumerate(cluster_data.items()):
            ax1.hist(values, alpha=0.7, label=f'Cluster {cluster_id}', color=colors[i])

        ax1.set_xlabel('Mean Embedding Value')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Embedding Distributions by Cluster')
        ax1.legend()

    # Plot 2: Embedding similarity heatmap
    ax2 = axes[1]
    event_types = list(embeddings.keys())
    n_events = len(event_types)

    if n_events > 1:
        # Calculate pairwise similarities (cosine similarity)
        similarity_matrix = np.zeros((n_events, n_events))

        for i in range(n_events):
            for j in range(n_events):
                emb1 = embeddings[event_types[i]]
                emb2 = embeddings[event_types[j]]
                # Cosine similarity
                dot_product = np.dot(emb1, emb2)
                norm1 = np.linalg.norm(emb1)
                norm2 = np.linalg.norm(emb2)
                if norm1 > 0 and norm2 > 0:
                    similarity_matrix[i, j] = dot_product / (norm1 * norm2)

        im = ax2.imshow(similarity_matrix, cmap='viridis', aspect='auto')
        ax2.set_xticks(range(n_events))
        ax2.set_yticks(range(n_events))
        ax2.set_xticklabels([et.split(':')[-1] for et in event_types], rotation=45, ha='right')
        ax2.set_yticklabels([et.split(':')[-1] for et in event_types])
        ax2.set_title('Embedding Similarity Matrix')

        # Add colorbar
        plt.colorbar(im, ax=ax2, shrink=0.8)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved embedding clusters plot to {output_path}")

    return fig


def plot_domain_timeline(
    sequences: List[EventSequence],
    output_path: Optional[Union[str, Path]] = None,
    max_sequences: int = 10,
    figsize: Tuple[int, int] = (12, 8)
) -> Any:
    """Plot domain timeline showing life domains over time.

    Args:
        sequences: List of EventSequence objects
        output_path: Path to save the plot (optional)
        max_sequences: Maximum number of sequences to plot
        figsize: Figure size (width, height)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create domain timeline plot")
        return None

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved domain timeline plot to {output_path}")

    return fig


def plot_population_comparison(
    population1: List[Any],
    population2: List[Any],
    labels: Tuple[str, str] = ("Population 1", "Population 2"),
    metric: str = "sequence_length",
    output_path: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 6)
) -> Any:
    """Compare two populations based on a metric.

    Args:
        population1: List of EventSequence objects for population 1
        population2: List of EventSequence objects for population 2
        labels: Labels for the two populations
        metric: Metric to compare ('sequence_length', 'event_count')
        output_path: Path to save the plot (optional)
        figsize: Figure size

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available")
        return None

    if metric == "sequence_length" or metric == "event_count":
        data1 = [len(seq.events) for seq in population1]
        data2 = [len(seq.events) for seq in population2]
    else:
        logger.warning(f"Unknown metric: {metric}")
        return None

    fig, ax = plt.subplots(figsize=figsize)
    ax.boxplot([data1, data2], labels=labels)
    ax.set_ylabel(metric.replace('_', ' ').title())
    ax.set_title(f'Population Comparison: {metric}')
    
    if output_path:
        plt.savefig(output_path)
        
    return fig

    fig, ax = plt.subplots(figsize=figsize)

    # Define colors for different domains
    domain_colors = {
        'education': '#1f77b4',  # blue
        'career': '#ff7f0e',     # orange
        'family': '#2ca02c',     # green
        'health': '#d62728',     # red
        'finance': '#9467bd',    # purple
        'housing': '#8c564b',    # brown
        'other': '#7f7f7f'       # gray
    }

    # Plot each sequence
    for i, seq in enumerate(sequences):
        y_pos = len(sequences) - i - 1  # Bottom to top

        # Group events by domain
        domain_periods = {}
        for event in seq.events:
            domain = getattr(event, 'domain', 'other')
            if domain not in domain_periods:
                domain_periods[domain] = []
            domain_periods[domain].append(event.timestamp)

        # Plot periods for each domain
        for domain, timestamps in domain_periods.items():
            if timestamps:
                # Sort timestamps
                timestamps.sort()
                color = domain_colors.get(domain, domain_colors['other'])

                # Plot events as points
                for ts in timestamps:
                    ax.scatter(ts, y_pos, color=color, s=50, alpha=0.7, marker='o')

                # Add domain label
                ax.text(timestamps[0], y_pos + 0.3, domain.capitalize(),
                       fontsize=8, ha='left', va='bottom', color=color)

    # Customize plot
    ax.set_yticks(range(len(sequences)))
    ax.set_yticklabels([f"Person {seq.person_id}" for seq in sequences])
    ax.set_xlabel('Time')
    ax.set_title('Life Domain Timeline')
    ax.grid(True, alpha=0.3)

    # Add legend
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w',
                                 markerfacecolor=color, markersize=8,
                                 label=domain.capitalize())
                      for domain, color in domain_colors.items()]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved domain timeline plot to {output_path}")

    return fig

def plot_population_comparison(
    group1_sequences: List[EventSequence],
    group2_sequences: List[EventSequence],
    group1_label: str = "Group 1",
    group2_label: str = "Group 2",
    output_path: Optional[Union[str, Path]] = None,
    figsize: Tuple[int, int] = (12, 8)
) -> Any:
    """Plot comparison of two populations based on their life event sequences.

    Args:
        group1_sequences: List of EventSequence objects for group 1
        group2_sequences: List of EventSequence objects for group 2
        group1_label: Label for group 1
        group2_label: Label for group 2
        output_path: Path to save the plot (optional)
        figsize: Figure size (width, height)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create population comparison plot")
        return None

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Extract sequence lengths
    group1_lengths = [len(seq.events) for seq in group1_sequences]
    group2_lengths = [len(seq.events) for seq in group2_sequences]

    # Plot 1: Sequence length comparison
    axes[0, 0].boxplot([group1_lengths, group2_lengths], labels=[group1_label, group2_label])
    axes[0, 0].set_ylabel('Number of Events')
    axes[0, 0].set_title('Sequence Length Comparison')
    axes[0, 0].grid(True, alpha=0.3)

    # Plot 2: Event type diversity
    group1_event_types = set()
    for seq in group1_sequences:
        for event in seq.events:
            group1_event_types.add(event.event_type)

    group2_event_types = set()
    for seq in group2_sequences:
        for event in seq.events:
            group2_event_types.add(event.event_type)

    diversity_data = [len(group1_event_types), len(group2_event_types)]
    axes[0, 1].bar([group1_label, group2_label], diversity_data, alpha=0.7)
    axes[0, 1].set_ylabel('Unique Event Types')
    axes[0, 1].set_title('Event Type Diversity')
    axes[0, 1].grid(True, alpha=0.3)

    # Plot 3: Event frequency comparison (top 5 events)
    all_events = set()
    for seq in group1_sequences + group2_sequences:
        for event in seq.events:
            all_events.add(event.event_type)

    event_counts_group1 = {event: 0 for event in all_events}
    event_counts_group2 = {event: 0 for event in all_events}

    for seq in group1_sequences:
        for event in seq.events:
            event_counts_group1[event.event_type] += 1

    for seq in group2_sequences:
        for event in seq.events:
            event_counts_group2[event.event_type] += 1

    # Get top 5 events by total count
    total_counts = {event: event_counts_group1[event] + event_counts_group2[event]
                   for event in all_events}
    top_events = sorted(total_counts.items(), key=lambda x: x[1], reverse=True)[:5]
    top_event_names = [event for event, count in top_events]

    group1_top_counts = [event_counts_group1[event] for event in top_event_names]
    group2_top_counts = [event_counts_group2[event] for event in top_event_names]

    x = np.arange(len(top_event_names))
    width = 0.35

    axes[1, 0].bar(x - width/2, group1_top_counts, width, label=group1_label, alpha=0.7)
    axes[1, 0].bar(x + width/2, group2_top_counts, width, label=group2_label, alpha=0.7)
    axes[1, 0].set_ylabel('Event Count')
    axes[1, 0].set_title('Top Event Frequencies')
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels([name.split(':')[-1] for name in top_event_names], rotation=45, ha='right')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)

    # Plot 4: Sequence similarity within groups
    def calculate_sequence_similarity(sequences):
        if len(sequences) < 2:
            return 0.5  # Default similarity

        # Simple similarity based on event type overlap
        similarities = []
        for i in range(len(sequences)):
            for j in range(i+1, len(sequences)):
                events1 = set(event.event_type for event in sequences[i].events)
                events2 = set(event.event_type for event in sequences[j].events)
                if events1 or events2:
                    similarity = len(events1 & events2) / len(events1 | events2)
                    similarities.append(similarity)

        return np.mean(similarities) if similarities else 0.5

    group1_similarity = calculate_sequence_similarity(group1_sequences)
    group2_similarity = calculate_sequence_similarity(group2_sequences)

    axes[1, 1].bar([group1_label, group2_label], [group1_similarity, group2_similarity], alpha=0.7)
    axes[1, 1].set_ylabel('Mean Sequence Similarity')
    axes[1, 1].set_title('Within-Group Similarity')
    axes[1, 1].set_ylim(0, 1)
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved population comparison plot to {output_path}")

    return fig

def plot_prediction_accuracy(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    task_type: str = "classification",
    output_path: Optional[Union[str, Path]] = None,
    figsize: Tuple[int, int] = (12, 8)
) -> Any:
    """Plot prediction accuracy and related metrics.

    Args:
        y_true: True target values
        y_pred: Predicted values
        task_type: Type of prediction task ('classification', 'regression')
        output_path: Path to save the plot (optional)
        figsize: Figure size (width, height)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create prediction accuracy plot")
        return None

    if task_type == "classification":
        # Binary classification metrics
        fig, axes = plt.subplots(2, 2, figsize=figsize)

        # Confusion matrix
        from sklearn.metrics import confusion_matrix
        cm = confusion_matrix(y_true, y_pred)
        axes[0, 0].imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
        axes[0, 0].set_title('Confusion Matrix')
        axes[0, 0].set_xlabel('Predicted')
        axes[0, 0].set_ylabel('True')

        # Add text annotations
        thresh = cm.max() / 2.
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                axes[0, 0].text(j, i, format(cm[i, j], 'd'),
                               ha="center", va="center",
                               color="white" if cm[i, j] > thresh else "black")

        # ROC curve
        from sklearn.metrics import roc_curve, auc
        fpr, tpr, _ = roc_curve(y_true, y_pred)
        roc_auc = auc(fpr, tpr)

        axes[0, 1].plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
        axes[0, 1].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        axes[0, 1].set_xlim([0.0, 1.0])
        axes[0, 1].set_ylim([0.0, 1.05])
        axes[0, 1].set_xlabel('False Positive Rate')
        axes[0, 1].set_ylabel('True Positive Rate')
        axes[0, 1].set_title('Receiver Operating Characteristic')
        axes[0, 1].legend(loc="lower right")

        # Precision-Recall curve
        from sklearn.metrics import precision_recall_curve, average_precision_score
        precision, recall, _ = precision_recall_curve(y_true, y_pred)
        average_precision = average_precision_score(y_true, y_pred)

        axes[1, 0].plot(recall, precision, color='blue', lw=2,
                       label=f'Precision-Recall curve (AP = {average_precision:.2f})')
        axes[1, 0].set_xlabel('Recall')
        axes[1, 0].set_ylabel('Precision')
        axes[1, 0].set_title('Precision-Recall Curve')
        axes[1, 0].legend(loc="lower left")

        # Accuracy over different thresholds (simplified)
        thresholds = np.linspace(0.1, 0.9, 9)
        accuracies = []
        for thresh in thresholds:
            pred_binary = (y_pred >= thresh).astype(int)
            acc = np.mean(pred_binary == y_true)
            accuracies.append(acc)

        axes[1, 1].plot(thresholds, accuracies, 'o-', color='green')
        axes[1, 1].set_xlabel('Threshold')
        axes[1, 1].set_ylabel('Accuracy')
        axes[1, 1].set_title('Accuracy vs Threshold')
        axes[1, 1].grid(True)

    elif task_type == "regression":
        # Regression metrics
        fig, axes = plt.subplots(2, 2, figsize=figsize)

        # Predicted vs Actual scatter
        axes[0, 0].scatter(y_true, y_pred, alpha=0.6)
        axes[0, 0].plot([y_true.min(), y_true.max()], [y_true.min(), y_true.max()], 'r--', lw=2)
        axes[0, 0].set_xlabel('True Values')
        axes[0, 0].set_ylabel('Predicted Values')
        axes[0, 0].set_title('Predicted vs Actual')

        # Residuals plot
        residuals = y_pred - y_true
        axes[0, 1].scatter(y_pred, residuals, alpha=0.6)
        axes[0, 1].axhline(y=0, color='r', linestyle='--')
        axes[0, 1].set_xlabel('Predicted Values')
        axes[0, 1].set_ylabel('Residuals')
        axes[0, 1].set_title('Residuals Plot')

        # Residuals distribution
        axes[1, 0].hist(residuals, bins=20, alpha=0.7, edgecolor='black')
        axes[1, 0].set_xlabel('Residual Value')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].set_title('Residuals Distribution')

        # Q-Q plot of residuals
        from scipy import stats
        probplot = stats.probplot(residuals, dist="norm", plot=axes[1, 1])
        axes[1, 1].set_title('Q-Q Plot of Residuals')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved prediction accuracy plot to {output_path}")

    return fig

def plot_sequence_length_distribution(
    sequences: List[EventSequence],
    output_path: Optional[Union[str, Path]] = None,
    figsize: Tuple[int, int] = (10, 6)
) -> Any:
    """Plot distribution of sequence lengths.

    Args:
        sequences: List of EventSequence objects
        output_path: Path to save the plot (optional)
        figsize: Figure size (width, height)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create sequence length distribution plot")
        return None

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Extract sequence lengths
    lengths = [len(seq.events) for seq in sequences]

    if not lengths:
        logger.warning("No sequences provided")
        return None

    # Plot histogram
    axes[0].hist(lengths, bins=20, alpha=0.7, edgecolor='black')
    axes[0].set_xlabel('Sequence Length')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title('Sequence Length Distribution')
    axes[0].grid(True, alpha=0.3)

    # Plot box plot
    axes[1].boxplot(lengths, vert=False)
    axes[1].set_xlabel('Sequence Length')
    axes[1].set_title('Sequence Length Box Plot')
    axes[1].grid(True, alpha=0.3)

    # Add statistics text
    mean_len = np.mean(lengths)
    median_len = np.median(lengths)
    std_len = np.std(lengths)
    min_len = np.min(lengths)
    max_len = np.max(lengths)

    stats_text = f'n = {len(lengths)}\nMean: {mean_len:.1f}\nMedian: {median_len:.1f}\nStd: {std_len:.1f}\nRange: {min_len}-{max_len}'
    axes[0].text(0.95, 0.95, stats_text, transform=axes[0].transAxes,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved sequence length distribution plot to {output_path}")

    return fig

def plot_sequence_similarity(
    sequences: List[EventSequence],
    output_path: Optional[Union[str, Path]] = None,
    figsize: Tuple[int, int] = (10, 8)
) -> Any:
    """Plot pairwise sequence similarities.

    Args:
        sequences: List of EventSequence objects
        output_path: Path to save the plot (optional)
        figsize: Figure size (width, height)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create sequence similarity plot")
        return None

    try:
        from sklearn.metrics.pairwise import cosine_similarity
        from sklearn.feature_extraction.text import TfidfVectorizer
        HAS_SKLEARN = True
    except ImportError:
        logger.warning("sklearn not available, cannot create sequence similarity plot")
        HAS_SKLEARN = False
        return None

    # Convert sequences to feature vectors
    sequence_strings = []
    for seq in sequences:
        # Create a string representation of the sequence
        event_types = [event.event_type for event in seq.events]
        sequence_strings.append(' '.join(event_types))

    if len(sequence_strings) < 2:
        logger.warning("Need at least 2 sequences for similarity analysis")
        return None

    # Use TF-IDF vectorization
    vectorizer = TfidfVectorizer()
    try:
        tfidf_matrix = vectorizer.fit_transform(sequence_strings)
        
        # Calculate cosine similarity
        similarity_matrix = cosine_similarity(tfidf_matrix)
        
        # Plot heatmap
        fig, ax = plt.subplots(figsize=figsize)
        im = ax.imshow(similarity_matrix, cmap='viridis', aspect='auto', vmin=0, vmax=1)
        
        # Set labels
        n_seq = len(sequences)
        ax.set_xticks(range(n_seq))
        ax.set_yticks(range(n_seq))
        ax.set_xticklabels([f'Seq {i+1}' for i in range(n_seq)], rotation=45, ha='right')
        ax.set_yticklabels([f'Seq {i+1}' for i in range(n_seq)])
        
        ax.set_title('Sequence Similarity Matrix (Cosine Similarity)')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Similarity Score')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved sequence similarity plot to {output_path}")
        
        return fig
        
    except Exception as e:
        logger.warning(f"Failed to compute similarities: {e}")
        return None

def plot_temporal_density(
    sequences: List[EventSequence],
    output_path: Optional[Union[str, Path]] = None,
    bins: int = 20,
    figsize: Tuple[int, int] = (12, 6)
) -> Any:
    """Plot temporal density of events across all sequences.

    Args:
        sequences: List of EventSequence objects
        output_path: Path to save the plot (optional)
        bins: Number of time bins for density estimation
        figsize: Figure size (width, height)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create temporal density plot")
        return None

    # Collect all event timestamps
    all_timestamps = []
    for seq in sequences:
        for event in seq.events:
            all_timestamps.append(event.timestamp)

    if not all_timestamps:
        logger.warning("No timestamps found in sequences")
        return None

    # Convert timestamps to numerical values (days since earliest event)
    timestamps_numeric = [(ts - min(all_timestamps)).total_seconds() / (24 * 3600) 
                         for ts in all_timestamps]

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Plot histogram
    counts, bin_edges, _ = axes[0].hist(timestamps_numeric, bins=bins, alpha=0.7, edgecolor='black')
    axes[0].set_xlabel('Time (days)')
    axes[0].set_ylabel('Event Count')
    axes[0].set_title('Event Temporal Distribution')
    axes[0].grid(True, alpha=0.3)

    # Plot density estimate
    try:
        from scipy import stats
        # Use Gaussian kernel density estimation
        kde = stats.gaussian_kde(timestamps_numeric)
        x_range = np.linspace(min(timestamps_numeric), max(timestamps_numeric), 200)
        density = kde(x_range)
        
        axes[1].plot(x_range, density, 'b-', linewidth=2)
        axes[1].fill_between(x_range, density, alpha=0.3)
        axes[1].set_xlabel('Time (days)')
        axes[1].set_ylabel('Density')
        axes[1].set_title('Event Temporal Density (KDE)')
        axes[1].grid(True, alpha=0.3)
        
    except ImportError:
        # Fallback: just show histogram again
        axes[1].hist(timestamps_numeric, bins=bins, alpha=0.7, edgecolor='black')
        axes[1].set_xlabel('Time (days)')
        axes[1].set_ylabel('Event Count')
        axes[1].set_title('Event Temporal Distribution (Histogram)')
        axes[1].grid(True, alpha=0.3)

    # Add statistics
    mean_time = np.mean(timestamps_numeric)
    median_time = np.median(timestamps_numeric)
    std_time = np.std(timestamps_numeric)
    
    stats_text = f'Events: {len(all_timestamps)}\nMean: {mean_time:.1f} days\nMedian: {median_time:.1f} days\nStd: {std_time:.1f} days'
    axes[0].text(0.02, 0.98, stats_text, transform=axes[0].transAxes,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved temporal density plot to {output_path}")

    return fig

def plot_temporal_patterns(
    sequences: List[EventSequence],
    importance_scores: Optional[Dict[int, float]] = None,
    output_path: Optional[Union[str, Path]] = None,
    figsize: Tuple[int, int] = (14, 8)
) -> Any:
    """Plot temporal patterns of life events with importance highlighting.

    Args:
        sequences: List of EventSequence objects
        importance_scores: Optional dict mapping sequence indices to importance scores
        output_path: Path to save the plot (optional)
        figsize: Figure size (width, height)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create temporal patterns plot")
        return None

    fig, axes = plt.subplots(2, 1, figsize=figsize)

    # Plot 1: Timeline view of all sequences
    ax1 = axes[0]
    colors = plt.cm.viridis(np.linspace(0, 1, len(sequences)))

    # Find global time range
    all_timestamps = []
    for seq in sequences:
        for event in seq.events:
            all_timestamps.append(event.timestamp)

    if all_timestamps:
        min_time = min(all_timestamps)
        max_time = max(all_timestamps)
        time_range = max_time - min_time

        for i, seq in enumerate(sequences):
            seq_color = colors[i]
            # Adjust color intensity based on importance
            if importance_scores and i in importance_scores:
                alpha = 0.3 + 0.7 * importance_scores[i]  # 0.3 to 1.0
                seq_color = (*seq_color[:3], alpha)

            for event in seq.events:
                # Convert to relative time
                rel_time = (event.timestamp - min_time).total_seconds() / (24 * 3600)  # days
                ax1.scatter(rel_time, i, color=seq_color, s=50, alpha=0.8)
                ax1.text(rel_time, i + 0.1, event.event_type.split(':')[-1],
                        fontsize=6, ha='center', va='bottom', rotation=45)

    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Sequence Index')
    ax1.set_title('Temporal Patterns of Life Events')
    ax1.grid(True, alpha=0.3)

    # Plot 2: Event type frequency over time
    ax2 = axes[1]

    # Bin events by time
    if all_timestamps:
        time_bins = np.linspace(0, (max_time - min_time).total_seconds() / (24 * 3600), 10)

        # Collect event frequencies by type and time bin
        event_types = set()
        for seq in sequences:
            for event in seq.events:
                event_types.add(event.event_type)

        event_types = sorted(list(event_types))

        # Create frequency matrix
        freq_matrix = np.zeros((len(event_types), len(time_bins) - 1))

        for seq in sequences:
            for event in seq.events:
                rel_time = (event.timestamp - min_time).total_seconds() / (24 * 3600)
                time_bin = np.digitize(rel_time, time_bins) - 1
                if 0 <= time_bin < len(time_bins) - 1:
                    event_idx = event_types.index(event.event_type)
                    freq_matrix[event_idx, time_bin] += 1

        # Plot as heatmap
        im = ax2.imshow(freq_matrix, aspect='auto', cmap='Blues')
        ax2.set_xticks(range(len(time_bins) - 1))
        ax2.set_yticks(range(len(event_types)))
        ax2.set_xticklabels([f'{time_bins[i]:.1f}' for i in range(len(time_bins) - 1)])
        ax2.set_yticklabels([et.split(':')[-1] for et in event_types])
        ax2.set_xlabel('Time Bin (days)')
        ax2.set_ylabel('Event Type')
        ax2.set_title('Event Frequency Over Time')

        # Add colorbar
        plt.colorbar(im, ax=ax2, shrink=0.8)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved temporal patterns plot to {output_path}")

    return fig

def plot_transition_network(sequences: List[EventSequence], output_path: Optional[Union[str, Path]] = None,
                          top_n: int = 10, figsize: Tuple[int, int] = (12, 8)) -> Any:
    """Plot transition network between event types.

    Args:
        sequences: List of EventSequence objects
        output_path: Path to save the plot (optional)
        top_n: Number of top transitions to show
        figsize: Figure size (width, height)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create transition network plot")
        return None

    try:
        import networkx as nx
        HAS_NETWORKX = True
    except ImportError:
        logger.warning("networkx not available, cannot create transition network plot")
        HAS_NETWORKX = False
        return None

    # Build transition matrix
    transitions = {}
    event_types = set()

    for seq in sequences:
        events = seq.events
        for i in range(len(events) - 1):
            from_event = events[i].event_type
            to_event = events[i + 1].event_type

            event_types.add(from_event)
            event_types.add(to_event)

            key = (from_event, to_event)
            transitions[key] = transitions.get(key, 0) + 1

    if not transitions:
        logger.warning("No transitions found in sequences")
        return None

    # Get top N transitions
    sorted_transitions = sorted(transitions.items(), key=lambda x: x[1], reverse=True)[:top_n]
    top_event_types = set()
    for (from_event, to_event), _ in sorted_transitions:
        top_event_types.add(from_event)
        top_event_types.add(to_event)

    # Create network
    G = nx.DiGraph()

    # Add nodes
    for event_type in top_event_types:
        G.add_node(event_type, label=event_type.split(':')[-1])

    # Add edges with weights
    for (from_event, to_event), weight in sorted_transitions:
        G.add_edge(from_event, to_event, weight=weight)

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)

    # Calculate positions
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    # Draw nodes
    node_sizes = [G.degree(node) * 100 + 300 for node in G.nodes()]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes,
                          node_color='lightblue', alpha=0.7)

    # Draw edges
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]
    max_weight = max(weights) if weights else 1

    # Scale edge widths
    edge_widths = [w / max_weight * 5 + 1 for w in weights]
    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths,
                          edge_color='gray', alpha=0.6, arrows=True,
                          arrowsize=20, arrowstyle='->')

    # Draw labels
    labels = {node: node.split(':')[-1] for node in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=8, font_weight='bold')

    # Add edge labels for weights
    edge_labels = {(u, v): str(G[u][v]['weight']) for u, v in edges}
    nx.draw_networkx_edge_labels(G, pos, edge_labels, ax=ax, font_size=6)

    ax.set_title(f'Event Transition Network\n(Top {top_n} transitions)')
    ax.axis('off')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved transition network plot to {output_path}")

    return fig


def plot_prediction_accuracy(
    y_true: Any,
    y_pred: Any,
    metric: str = 'accuracy',
    output_path: Optional[str] = None
) -> Any:
    """Plot prediction accuracy metrics.

    Args:
        y_true: True labels
        y_pred: Predicted labels
        metric: Metric name
        output_path: Output path

    Returns:
        matplotlib Figure
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available")
        return None
        
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Placeholder plot (scatter)
    ax.scatter(y_true, y_pred, alpha=0.5)
    
    if len(y_true) > 0 and len(y_pred) > 0:
        min_val = min(np.min(y_true), np.min(y_pred))
        max_val = max(np.max(y_true), np.max(y_pred))
        ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.7)
    
    ax.set_xlabel('True Values')
    ax.set_ylabel('Predicted Values')
    ax.set_title(f'Prediction {metric.title()}')
    
    if output_path:
        plt.savefig(output_path)
        
    return fig

