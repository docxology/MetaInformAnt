"""Visualization functions for event sequences.

This module provides plotting functions for visualizing event sequences,
embeddings, and analysis results. All functions use matplotlib and integrate
with the main visualization module for consistent styling.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from ..visualization import plots
from .events import EventSequence

if TYPE_CHECKING:
    from matplotlib.figure import Figure


def plot_event_timeline(
    sequence: EventSequence,
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (12, 6)
) -> "Figure":
    """Plot timeline visualization of individual life course.
    
    Args:
        sequence: EventSequence to visualize
        output_path: Optional path to save figure
        figsize: Figure size
        
    Returns:
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Get timestamps
    def get_timestamp(event):
        if isinstance(event.timestamp, float):
            return event.timestamp
        return event.timestamp.timestamp()
    
    timestamps = [get_timestamp(e) for e in sequence.events]
    
    # Map domains to colors
    domain_colors = {
        "health": "red",
        "education": "blue",
        "occupation": "green",
        "income": "orange",
        "address": "purple",
        "other": "gray"
    }
    
    # Plot events
    y_pos = 0
    for i, event in enumerate(sequence.events):
        timestamp = timestamps[i]
        color = domain_colors.get(event.domain, "gray")
        
        ax.scatter(timestamp, y_pos, c=color, s=100, alpha=0.7)
        ax.text(timestamp, y_pos + 0.1, event.event_type, 
                rotation=45, ha='left', fontsize=8)
        y_pos += 1
    
    ax.set_xlabel("Time")
    ax.set_ylabel("Event Index")
    ax.set_title(f"Life Course Timeline: {sequence.person_id}")
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_event_embeddings(
    embeddings: Dict[str, NDArray],
    method: str = "umap",
    n_components: int = 2,
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (10, 8)
) -> "Figure":
    """Plot 2D/3D visualization of event embeddings.
    
    Args:
        embeddings: Dictionary mapping event tokens to embedding vectors
        method: Dimensionality reduction method ("pca", "umap", "tsne")
        n_components: Number of components (2 or 3)
        output_path: Optional path to save figure
        figsize: Figure size
        
    Returns:
        Matplotlib figure
    """
    from ..ml.dimensionality import biological_embedding
    
    # Convert to matrix
    tokens = sorted(list(embeddings.keys()))
    embedding_matrix = np.array([embeddings[t] for t in tokens])
    
    # Reduce dimensionality
    reduced = biological_embedding(
        embedding_matrix,
        method=method,
        n_components=n_components
    )
    
    coords = reduced["embedding"]
    
    # Plot
    if n_components == 2:
        fig, ax = plt.subplots(figsize=figsize)
        ax.scatter(coords[:, 0], coords[:, 1], alpha=0.6)
        
        # Label some events
        for i, token in enumerate(tokens[:20]):  # Label first 20
            ax.annotate(token, (coords[i, 0], coords[i, 1]), fontsize=6)
        
        ax.set_xlabel(f"{method.upper()} Component 1")
        ax.set_ylabel(f"{method.upper()} Component 2")
        ax.set_title("Event Embedding Space")
        ax.grid(True, alpha=0.3)
    else:
        try:
            from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], alpha=0.6)
            ax.set_xlabel(f"{method.upper()} Component 1")
            ax.set_ylabel(f"{method.upper()} Component 2")
            ax.set_zlabel(f"{method.upper()} Component 3")
            ax.set_title("Event Embedding Space (3D)")
        except ImportError:
            # Fallback to 2D if 3D not available
            fig, ax = plt.subplots(figsize=figsize)
            ax.scatter(coords[:, 0], coords[:, 1], alpha=0.6)
            ax.set_xlabel(f"{method.upper()} Component 1")
            ax.set_ylabel(f"{method.upper()} Component 2")
            ax.set_title("Event Embedding Space (2D projection)")
    
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_attention_heatmap(
    attention_weights: NDArray,
    event_tokens: List[str],
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (12, 10)
) -> "Figure":
    """Plot attention weight heatmap for sequence models.
    
    Args:
        attention_weights: Attention weight matrix (events x events)
        event_tokens: List of event tokens for labels
        output_path: Optional path to save figure
        figsize: Figure size
        
    Returns:
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot heatmap
    plots.heatmap(
        attention_weights,
        ax=ax,
        cmap="viridis",
        annot=False
    )
    
    # Set labels
    if len(event_tokens) <= 20:
        ax.set_xticks(range(len(event_tokens)))
        ax.set_yticks(range(len(event_tokens)))
        ax.set_xticklabels(event_tokens, rotation=45, ha='right')
        ax.set_yticklabels(event_tokens)
    
    ax.set_xlabel("Event (Query)")
    ax.set_ylabel("Event (Key)")
    ax.set_title("Attention Weights")
    
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_prediction_importance(
    event_importance: Dict[str, float],
    top_n: int = 20,
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (10, 8)
) -> "Figure":
    """Plot event importance for predictions.
    
    Args:
        event_importance: Dictionary mapping event tokens to importance scores
        top_n: Number of top events to display
        output_path: Optional path to save figure
        figsize: Figure size
        
    Returns:
        Matplotlib figure
    """
    # Sort by importance
    sorted_items = sorted(
        event_importance.items(),
        key=lambda x: x[1],
        reverse=True
    )[:top_n]
    
    events = [item[0] for item in sorted_items]
    importance = [item[1] for item in sorted_items]
    
    fig, ax = plt.subplots(figsize=figsize)
    
    y_pos = np.arange(len(events))
    ax.barh(y_pos, importance, alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(events)
    ax.set_xlabel("Importance Score")
    ax.set_title(f"Top {top_n} Events by Prediction Importance")
    ax.invert_yaxis()  # Highest importance at top
    
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_domain_distribution(
    sequences: List[EventSequence],
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (10, 6),
    plot_type: str = "bar"
) -> "Figure":
    """Plot distribution of domains across event sequences."""
    domain_counts = Counter()
    for seq in sequences:
        for event in seq.events:
            domain_counts[event.domain] += 1
    
    domains = list(domain_counts.keys())
    counts = [domain_counts[d] for d in domains]
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if plot_type == "pie":
        ax.pie(counts, labels=domains, autopct='%1.1f%%', startangle=90)
        ax.set_title("Domain Distribution")
    else:
        ax.bar(domains, counts, alpha=0.7)
        ax.set_xlabel("Domain")
        ax.set_ylabel("Event Count")
        ax.set_title("Domain Distribution")
        plt.xticks(rotation=45, ha='right')
    
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig


def plot_temporal_density(
    sequences: List[EventSequence],
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (12, 6),
    bins: int = 50
) -> "Figure":
    """Plot temporal density of events over time."""
    all_timestamps = []
    for seq in sequences:
        for event in seq.events:
            if isinstance(event.timestamp, datetime):
                all_timestamps.append(event.timestamp.timestamp())
            else:
                all_timestamps.append(float(event.timestamp))
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(all_timestamps, bins=bins, alpha=0.7, edgecolor='black')
    ax.set_xlabel("Time (timestamp)")
    ax.set_ylabel("Event Count")
    ax.set_title("Temporal Event Density")
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig


def plot_event_cooccurrence(
    sequences: List[EventSequence],
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (12, 10),
    top_n: int = 20
) -> "Figure":
    """Plot heatmap of event co-occurrence patterns."""
    from .utils import convert_sequences_to_tokens
    
    sequences_tokens = convert_sequences_to_tokens(sequences)
    event_counts = Counter()
    for seq in sequences_tokens:
        event_counts.update(seq)
    
    top_events = [event for event, _ in event_counts.most_common(top_n)]
    cooccurrence = defaultdict(int)
    for seq in sequences_tokens:
        seq_set = set(seq)
        for event1 in top_events:
            if event1 in seq_set:
                for event2 in top_events:
                    if event2 in seq_set and event1 != event2:
                        cooccurrence[(event1, event2)] += 1
    
    matrix = np.zeros((len(top_events), len(top_events)))
    for i, event1 in enumerate(top_events):
        for j, event2 in enumerate(top_events):
            matrix[i, j] = cooccurrence.get((event1, event2), 0)
    
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')
    ax.set_xticks(range(len(top_events)))
    ax.set_yticks(range(len(top_events)))
    ax.set_xticklabels(top_events, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(top_events, fontsize=8)
    ax.set_xlabel("Event")
    ax.set_ylabel("Event")
    ax.set_title("Event Co-occurrence Heatmap")
    
    plt.colorbar(im, ax=ax)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig


def plot_outcome_distribution(
    outcomes: NDArray,
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (10, 6),
    plot_type: str = "histogram"
) -> "Figure":
    """Plot distribution of outcomes."""
    fig, ax = plt.subplots(figsize=figsize)
    
    if plot_type == "boxplot":
        ax.boxplot(outcomes)
        ax.set_ylabel("Outcome Value")
        ax.set_title("Outcome Distribution (Box Plot)")
    else:
        ax.hist(outcomes, bins=30, alpha=0.7, edgecolor='black')
        ax.set_xlabel("Outcome Value")
        ax.set_ylabel("Frequency")
        ax.set_title("Outcome Distribution")
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig


def plot_sequence_similarity(
    sequences: List[EventSequence],
    embeddings: Optional[Dict[str, NDArray]] = None,
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (12, 10)
) -> "Figure":
    """Plot heatmap of sequence similarity matrix."""
    from .utils import convert_sequences_to_tokens
    from .embeddings import learn_event_embeddings, sequence_embeddings
    
    sequences_tokens = convert_sequences_to_tokens(sequences)
    if embeddings is None:
        embeddings = learn_event_embeddings(sequences_tokens, embedding_dim=100)
    
    seq_embeddings = sequence_embeddings(sequences_tokens, embeddings, method="mean")
    
    try:
        from sklearn.metrics.pairwise import cosine_similarity
        similarity_matrix = cosine_similarity(seq_embeddings)
    except ImportError:
        similarity_matrix = np.dot(seq_embeddings, seq_embeddings.T)
        norms = np.linalg.norm(seq_embeddings, axis=1)
        similarity_matrix = similarity_matrix / np.outer(norms, norms)
    
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(similarity_matrix, cmap='viridis', aspect='auto', vmin=0, vmax=1)
    ax.set_xlabel("Sequence Index")
    ax.set_ylabel("Sequence Index")
    ax.set_title("Sequence Similarity Matrix")
    
    plt.colorbar(im, ax=ax)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig


def plot_transition_network(
    sequences: List[EventSequence],
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (14, 10),
    top_n: int = 15
) -> "Figure":
    """Plot network graph of event transitions."""
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx required for transition network visualization")
    
    from .utils import convert_sequences_to_tokens
    
    sequences_tokens = convert_sequences_to_tokens(sequences)
    transition_counts = defaultdict(int)
    event_counts = Counter()
    
    for seq in sequences_tokens:
        for event in seq:
            event_counts[event] += 1
        for i in range(len(seq) - 1):
            transition_counts[(seq[i], seq[i+1])] += 1
    
    top_events = [event for event, _ in event_counts.most_common(top_n)]
    G = nx.DiGraph()
    for (src, dst), count in transition_counts.items():
        if src in top_events and dst in top_events:
            if not G.has_edge(src, dst):
                G.add_edge(src, dst, weight=count)
            else:
                G[src][dst]['weight'] += count
    
    fig, ax = plt.subplots(figsize=figsize)
    pos = nx.spring_layout(G, k=2, iterations=50)
    edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, width=[w/10 for w in edge_weights], alpha=0.6, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=500, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)
    
    ax.set_title("Event Transition Network")
    ax.axis('off')
    
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig


def plot_domain_timeline(
    sequences: List[EventSequence],
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (14, 8),
    max_sequences: int = 10
) -> "Figure":
    """Plot multi-domain timeline (Gantt-style) for multiple sequences."""
    domain_colors = {
        "health": "red", "education": "blue", "occupation": "green",
        "income": "orange", "address": "purple", "other": "gray"
    }
    
    fig, ax = plt.subplots(figsize=figsize)
    sequences_to_plot = sequences[:max_sequences]
    
    for seq_idx, seq in enumerate(sequences_to_plot):
        y_pos = seq_idx
        for event in seq.events:
            timestamp = event.timestamp.timestamp() if isinstance(event.timestamp, datetime) else float(event.timestamp)
            color = domain_colors.get(event.domain, "gray")
            ax.scatter(timestamp, y_pos, c=color, s=100, alpha=0.7, label=event.domain if seq_idx == 0 else "")
    
    ax.set_yticks(range(len(sequences_to_plot)))
    ax.set_yticklabels([seq.person_id for seq in sequences_to_plot])
    ax.set_xlabel("Time")
    ax.set_title("Multi-Domain Timeline")
    ax.grid(True, alpha=0.3)
    
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right')
    
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig


def plot_prediction_accuracy(
    y_true: NDArray,
    y_pred: NDArray,
    task_type: str = "classification",
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (12, 5)
) -> "Figure":
    """Plot prediction accuracy metrics."""
    fig, axes = plt.subplots(1, 2 if task_type == "classification" else 1, figsize=figsize)
    if not isinstance(axes, np.ndarray):
        axes = [axes]
    
    if task_type == "classification":
        try:
            from sklearn.metrics import confusion_matrix
            cm = confusion_matrix(y_true, y_pred)
            axes[0].imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
            axes[0].set_title("Confusion Matrix")
            axes[0].set_ylabel("True Label")
            axes[0].set_xlabel("Predicted Label")
            thresh = cm.max() / 2.
            for i in range(cm.shape[0]):
                for j in range(cm.shape[1]):
                    axes[0].text(j, i, format(cm[i, j], 'd'),
                                ha="center", va="center",
                                color="white" if cm[i, j] > thresh else "black")
        except ImportError:
            axes[0].text(0.5, 0.5, 'sklearn required', ha='center', va='center', transform=axes[0].transAxes)
        
        if len(axes) > 1:
            try:
                from sklearn.metrics import roc_curve, auc
                if len(np.unique(y_true)) == 2:
                    fpr, tpr, _ = roc_curve(y_true, y_pred)
                    roc_auc = auc(fpr, tpr)
                    axes[1].plot(fpr, tpr, lw=2, label=f'ROC (AUC = {roc_auc:.2f})')
                    axes[1].plot([0, 1], [0, 1], 'k--', lw=2)
                    axes[1].set_xlabel('False Positive Rate')
                    axes[1].set_ylabel('True Positive Rate')
                    axes[1].set_title('ROC Curve')
                    axes[1].legend()
                    axes[1].grid(True, alpha=0.3)
            except Exception:
                pass
    else:
        axes[0].scatter(y_true, y_pred, alpha=0.6)
        min_val, max_val = min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())
        axes[0].plot([min_val, max_val], [min_val, max_val], 'r--', lw=2)
        axes[0].set_xlabel("True Values")
        axes[0].set_ylabel("Predicted Values")
        axes[0].set_title("Prediction Accuracy (Regression)")
        axes[0].grid(True, alpha=0.3)
    
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig


def plot_temporal_patterns(
    sequences: List[EventSequence],
    importance_scores: Optional[Dict[int, float]] = None,
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (12, 6)
) -> "Figure":
    """Plot time-based importance visualization."""
    fig, ax = plt.subplots(figsize=figsize)
    position_importance = defaultdict(list)
    
    for seq_idx, seq in enumerate(sequences):
        importance = importance_scores.get(seq_idx, 1.0) if importance_scores else 1.0
        for pos, event in enumerate(seq.events):
            position_importance[pos].append(importance)
    
    positions = sorted(position_importance.keys())
    avg_importance = [np.mean(position_importance[pos]) for pos in positions]
    
    ax.plot(positions, avg_importance, marker='o', linewidth=2, markersize=6)
    ax.set_xlabel("Event Position in Sequence")
    ax.set_ylabel("Average Importance")
    ax.set_title("Temporal Pattern Importance")
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig


def plot_population_comparison(
    sequences_group1: List[EventSequence],
    sequences_group2: List[EventSequence],
    group1_label: str = "Group 1",
    group2_label: str = "Group 2",
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (14, 6)
) -> "Figure":
    """Plot side-by-side comparison of two population groups."""
    domain_counts1 = Counter()
    domain_counts2 = Counter()
    
    for seq in sequences_group1:
        for event in seq.events:
            domain_counts1[event.domain] += 1
    
    for seq in sequences_group2:
        for event in seq.events:
            domain_counts2[event.domain] += 1
    
    domains = sorted(set(list(domain_counts1.keys()) + list(domain_counts2.keys())))
    counts1 = [domain_counts1.get(d, 0) for d in domains]
    counts2 = [domain_counts2.get(d, 0) for d in domains]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    x_pos = np.arange(len(domains))
    
    ax1.bar(x_pos, counts1, 0.35, alpha=0.7, label=group1_label)
    ax1.set_xlabel("Domain")
    ax1.set_ylabel("Event Count")
    ax1.set_title(f"Domain Distribution: {group1_label}")
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(domains, rotation=45, ha='right')
    ax1.grid(True, alpha=0.3)
    
    ax2.bar(x_pos, counts2, 0.35, alpha=0.7, label=group2_label, color='orange')
    ax2.set_xlabel("Domain")
    ax2.set_ylabel("Event Count")
    ax2.set_title(f"Domain Distribution: {group2_label}")
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(domains, rotation=45, ha='right')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig


def plot_intervention_effects(
    pre_sequences: List[EventSequence],
    post_sequences: List[EventSequence],
    pre_outcomes: Optional[NDArray] = None,
    post_outcomes: Optional[NDArray] = None,
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (14, 6)
) -> "Figure":
    """Plot before/after intervention visualization."""
    pre_domain_counts = Counter()
    post_domain_counts = Counter()
    
    for seq in pre_sequences:
        for event in seq.events:
            pre_domain_counts[event.domain] += 1
    
    for seq in post_sequences:
        for event in seq.events:
            post_domain_counts[event.domain] += 1
    
    domains = sorted(set(list(pre_domain_counts.keys()) + list(post_domain_counts.keys())))
    pre_counts = [pre_domain_counts.get(d, 0) for d in domains]
    post_counts = [post_domain_counts.get(d, 0) for d in domains]
    
    fig, axes = plt.subplots(1, 2 if pre_outcomes is not None else 1, figsize=figsize)
    if not isinstance(axes, np.ndarray):
        axes = [axes]
    
    x_pos = np.arange(len(domains))
    width = 0.35
    
    axes[0].bar(x_pos - width/2, pre_counts, width, alpha=0.7, label="Pre-Intervention")
    axes[0].bar(x_pos + width/2, post_counts, width, alpha=0.7, label="Post-Intervention")
    axes[0].set_xlabel("Domain")
    axes[0].set_ylabel("Event Count")
    axes[0].set_title("Domain Distribution: Before vs After Intervention")
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels(domains, rotation=45, ha='right')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    if pre_outcomes is not None and post_outcomes is not None and len(axes) > 1:
        axes[1].boxplot([pre_outcomes, post_outcomes], labels=["Pre", "Post"])
        axes[1].set_ylabel("Outcome Value")
        axes[1].set_title("Outcome Distribution: Before vs After")
        axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig


def plot_embedding_clusters(
    embeddings: Dict[str, NDArray],
    clusters: Optional[Dict[str, int]] = None,
    method: str = "umap",
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (10, 8)
) -> "Figure":
    """Plot clustered embedding visualization."""
    from ..ml.dimensionality import biological_embedding
    
    tokens = sorted(list(embeddings.keys()))
    embedding_matrix = np.array([embeddings[t] for t in tokens])
    reduced = biological_embedding(embedding_matrix, method=method, n_components=2)
    coords = reduced["embedding"]
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if clusters:
        cluster_ids = [clusters.get(t, 0) for t in tokens]
        unique_clusters = sorted(set(cluster_ids))
        colors = plt.cm.tab20(np.linspace(0, 1, len(unique_clusters)))
        
        for cluster_id in unique_clusters:
            mask = [c == cluster_id for c in cluster_ids]
            ax.scatter(coords[mask, 0], coords[mask, 1], 
                      c=[colors[unique_clusters.index(cluster_id)]],
                      label=f'Cluster {cluster_id}', alpha=0.6)
        ax.legend()
    else:
        ax.scatter(coords[:, 0], coords[:, 1], alpha=0.6)
    
    ax.set_xlabel(f"{method.upper()} Component 1")
    ax.set_ylabel(f"{method.upper()} Component 2")
    ax.set_title("Embedding Clusters")
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig


def plot_sequence_length_distribution(
    sequences: List[EventSequence],
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (10, 6)
) -> "Figure":
    """Plot histogram of sequence lengths."""
    lengths = [len(seq.events) for seq in sequences]
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(lengths, bins=30, alpha=0.7, edgecolor='black')
    ax.set_xlabel("Sequence Length (number of events)")
    ax.set_ylabel("Frequency")
    ax.set_title("Sequence Length Distribution")
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig


def plot_event_frequency_heatmap(
    sequences: List[EventSequence],
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (14, 8),
    time_bins: int = 10
) -> "Figure":
    """Plot temporal frequency heatmap of events."""
    all_timestamps = []
    for seq in sequences:
        for event in seq.events:
            if isinstance(event.timestamp, datetime):
                all_timestamps.append(event.timestamp.timestamp())
            else:
                all_timestamps.append(float(event.timestamp))
    
    if not all_timestamps:
        raise ValueError("No events found in sequences")
    
    min_time = min(all_timestamps)
    max_time = max(all_timestamps)
    time_bin_size = (max_time - min_time) / time_bins
    
    domain_time_counts = defaultdict(lambda: defaultdict(int))
    
    for seq in sequences:
        for event in seq.events:
            timestamp = event.timestamp.timestamp() if isinstance(event.timestamp, datetime) else float(event.timestamp)
            time_bin = int((timestamp - min_time) / time_bin_size)
            time_bin = min(time_bin, time_bins - 1)
            domain_time_counts[event.domain][time_bin] += 1
    
    domains = sorted(domain_time_counts.keys())
    matrix = np.zeros((len(domains), time_bins))
    
    for i, domain in enumerate(domains):
        for j in range(time_bins):
            matrix[i, j] = domain_time_counts[domain][j]
    
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')
    ax.set_yticks(range(len(domains)))
    ax.set_yticklabels(domains)
    ax.set_xlabel("Time Bin")
    ax.set_ylabel("Domain")
    ax.set_title("Event Frequency Heatmap (Temporal)")
    
    plt.colorbar(im, ax=ax)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
    return fig

