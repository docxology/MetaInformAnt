"""Visualization functions for event sequences.

This module provides plotting functions for visualizing event sequences,
embeddings, and analysis results. All functions use matplotlib and integrate
with the main visualization module for consistent styling.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional

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

