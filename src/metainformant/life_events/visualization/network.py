"""Network visualization functions for life events analysis.

This module provides network-related plotting utilities for life event sequences,
including transition networks between event types.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Optional dependencies
try:
    import matplotlib.pyplot as plt
    import numpy as np

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    logger.warning("matplotlib not available, life events network visualization disabled")


def plot_transition_network(
    sequences: List[Any],
    output_path: Optional[Union[str, Any]] = None,
    top_n: int = 10,
    figsize: Tuple[int, int] = (12, 8),
) -> Any:
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
        G.add_node(event_type, label=event_type.split(":")[-1])

    # Add edges with weights
    for (from_event, to_event), weight in sorted_transitions:
        G.add_edge(from_event, to_event, weight=weight)

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)

    # Calculate positions
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    # Draw nodes
    node_sizes = [G.degree(node) * 100 + 300 for node in G.nodes()]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes, node_color="lightblue", alpha=0.7)

    # Draw edges
    edges = G.edges()
    weights = [G[u][v]["weight"] for u, v in edges]
    max_weight = max(weights) if weights else 1

    # Scale edge widths
    edge_widths = [w / max_weight * 5 + 1 for w in weights]
    nx.draw_networkx_edges(
        G, pos, ax=ax, width=edge_widths, edge_color="gray", alpha=0.6, arrows=True, arrowsize=20, arrowstyle="->"
    )

    # Draw labels
    labels = {node: node.split(":")[-1] for node in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=8, font_weight="bold")

    # Add edge labels for weights
    edge_labels = {(u, v): str(G[u][v]["weight"]) for u, v in edges}
    nx.draw_networkx_edge_labels(G, pos, edge_labels, ax=ax, font_size=6)

    ax.set_title(f"Event Transition Network\n(Top {top_n} transitions)")
    ax.axis("off")

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved transition network plot to {output_path}")

    return fig
