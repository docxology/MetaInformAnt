"""Specialized visualization functions for advanced plotting types.

This module provides specialized visualization capabilities including Venn diagrams,
Sankey diagrams, chord diagrams, alluvial plots, and other advanced chart types
commonly used in bioinformatics and systems biology.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Circle, PathPatch, Wedge
from matplotlib.path import Path

from metainformant.core import logging, paths, validation

logger = logging.get_logger(__name__)

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    sns = None

try:
    import plotly.express as px
    import plotly.graph_objects as go

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    go = None
    px = None


def plot_venn_diagram(
    sets: Dict[str, set],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 8),
    **kwargs,
) -> Axes:
    """Create a Venn diagram from set data.

    Args:
        sets: Dictionary mapping set names to set objects
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(sets, dict, "sets")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    try:
        import matplotlib_venn as venn
    except ImportError:
        # Fallback implementation for 2-3 sets
        if len(sets) == 2:
            _plot_venn2(ax, sets, **kwargs)
        elif len(sets) == 3:
            _plot_venn3(ax, sets, **kwargs)
        else:
            raise ImportError("matplotlib-venn required for Venn diagrams with more than 3 sets")
    else:
        # Use matplotlib-venn
        set_names = list(sets.keys())
        set_data = [sets[name] for name in set_names]

        if len(set_data) == 2:
            venn.venn2(set_data, set_labels=set_names, ax=ax, **kwargs)
        elif len(set_data) == 3:
            venn.venn3(set_data, set_labels=set_names, ax=ax, **kwargs)
        else:
            venn.venn4(set_data, set_labels=set_names, ax=ax, **kwargs)

    ax.set_title("Venn Diagram")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Venn diagram saved to {output_path}")

    return ax


def _plot_venn2(ax: Axes, sets: Dict[str, set], **kwargs) -> None:
    """Simple fallback Venn diagram for 2 sets."""
    set_names = list(sets.keys())
    set1, set2 = [sets[name] for name in set_names]

    # Calculate set sizes
    only1 = len(set1 - set2)
    only2 = len(set2 - set1)
    both = len(set1 & set2)

    # Draw circles
    circle1 = Circle((-1, 0), 1.5, facecolor="red", alpha=0.5)
    circle2 = Circle((1, 0), 1.5, facecolor="blue", alpha=0.5)

    ax.add_patch(circle1)
    ax.add_patch(circle2)

    # Add labels
    ax.text(-1, 0, f"{both}\n({set_names[0]}∩{set_names[1]})", ha="center", va="center")
    ax.text(-2, 0, f"{only1}\n({set_names[0]} only)", ha="center", va="center")
    ax.text(2, 0, f"{only2}\n({set_names[1]} only)", ha="center", va="center")

    ax.set_xlim(-3, 3)
    ax.set_ylim(-2, 2)
    ax.axis("off")


def _plot_venn3(ax: Axes, sets: Dict[str, set], **kwargs) -> None:
    """Simple fallback Venn diagram for 3 sets."""
    set_names = list(sets.keys())
    set1, set2, set3 = [sets[name] for name in set_names]

    # Calculate intersections
    only1 = len(set1 - set2 - set3)
    only2 = len(set2 - set1 - set3)
    only3 = len(set3 - set1 - set2)
    only12 = len((set1 & set2) - set3)
    only13 = len((set1 & set3) - set2)
    only23 = len((set2 & set3) - set1)
    all3 = len(set1 & set2 & set3)

    # Draw circles
    circle1 = Circle((-1, 1), 1.2, facecolor="red", alpha=0.5)
    circle2 = Circle((1, 1), 1.2, facecolor="blue", alpha=0.5)
    circle3 = Circle((0, -0.8), 1.2, facecolor="green", alpha=0.5)

    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.add_patch(circle3)

    # Add labels (simplified)
    ax.text(-1, 1, f"{only1 + only12 + only13 + all3}", ha="center", va="center")
    ax.text(1, 1, f"{only2 + only12 + only23 + all3}", ha="center", va="center")
    ax.text(0, -0.8, f"{only3 + only13 + only23 + all3}", ha="center", va="center")
    ax.text(0, 0.3, f"{only12 + all3}", ha="center", va="center")
    ax.text(-0.5, -0.3, f"{only13 + all3}", ha="center", va="center")
    ax.text(0.5, -0.3, f"{only23 + all3}", ha="center", va="center")
    ax.text(0, 0, f"{all3}", ha="center", va="center")

    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.5)
    ax.axis("off")


def plot_sankey_diagram(
    flows: List[Tuple[str, str, float]],
    *,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Any:
    """Create a Sankey diagram for flow visualization.

    Args:
        flows: List of (source, target, value) tuples
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Figure or Plotly Figure object
    """
    validation.validate_type(flows, list, "flows")

    if HAS_PLOTLY:
        # Use Plotly for better Sankey diagrams
        sources = []
        targets = []
        values = []
        labels = set()

        for source, target, value in flows:
            labels.add(source)
            labels.add(target)
            sources.append(list(labels).index(source))
            targets.append(list(labels).index(target))
            values.append(value)

        label_list = list(labels)

        fig = go.Figure(
            data=[
                go.Sankey(
                    node=dict(pad=15, thickness=20, line=dict(color="black", width=0.5), label=label_list),
                    link=dict(source=sources, target=targets, value=values),
                )
            ]
        )

        fig.update_layout(title_text="Sankey Diagram", **kwargs)

        if output_path:
            paths.ensure_directory(Path(output_path).parent)
            html_path = Path(output_path).with_suffix(".html")
            fig.write_html(str(html_path))
            logger.info(f"Sankey diagram saved to {html_path}")

        return fig
    else:
        # Fallback matplotlib implementation
        fig, ax = plt.subplots(figsize=figsize)

        # Simple flow visualization
        sources = list(set(s for s, t, v in flows))
        targets = list(set(t for s, t, v in flows))

        source_positions = {src: i for i, src in enumerate(sources)}
        target_positions = {tgt: i for i, tgt in enumerate(targets)}

        for source, target, value in flows:
            src_pos = source_positions[source]
            tgt_pos = target_positions[target]

            # Draw flow line
            ax.plot([0, 1], [src_pos, tgt_pos], "b-", linewidth=value / 10, alpha=0.7)
            ax.scatter([0], [src_pos], s=value, c="red", alpha=0.7)
            ax.scatter([1], [tgt_pos], s=value, c="blue", alpha=0.7)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Sources", "Targets"])
        ax.set_yticks(list(source_positions.values()))
        ax.set_yticklabels(list(source_positions.keys()))
        ax.set_title("Flow Diagram")

        if output_path:
            paths.ensure_directory(Path(output_path).parent)
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
            logger.info(f"Flow diagram saved to {output_path}")

        return fig


def plot_chord_diagram(
    matrix: np.ndarray,
    labels: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 8),
    **kwargs,
) -> Axes:
    """Create a chord diagram for showing relationships between categories.

    Args:
        matrix: Square matrix of relationships
        labels: Optional labels for matrix axes
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(matrix, np.ndarray, "matrix")

    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Chord diagram matrix must be square")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": "polar"})

    n = matrix.shape[0]
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)

    # Normalize matrix for chord thickness
    max_val = matrix.max()
    if max_val > 0:
        matrix_norm = matrix / max_val
    else:
        matrix_norm = matrix

    # Draw chords
    for i in range(n):
        for j in range(i + 1, n):
            if matrix_norm[i, j] > 0:
                # Draw chord between i and j
                theta1 = angles[i]
                theta2 = angles[j]
                r = matrix_norm[i, j] * 0.3  # Chord radius

                # Create chord path
                verts = [
                    (theta1, 0.5),  # Start
                    (theta1, 0.5 + r),  # Arc start
                    ((theta1 + theta2) / 2, 0.5 + r),  # Middle
                    (theta2, 0.5 + r),  # Arc end
                    (theta2, 0.5),  # End
                ]

                codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CURVE4]
                path = Path(verts, codes)
                patch = PathPatch(path, facecolor="blue", alpha=0.3, edgecolor="blue")
                ax.add_patch(patch)

    # Draw node labels
    for i, angle in enumerate(angles):
        label = labels[i] if labels else f"Node {i+1}"
        ax.text(angle, 0.8, label, ha="center", va="center", fontsize=10, rotation=np.degrees(angle) - 90)

    ax.set_rlim(0, 1)
    ax.set_title("Chord Diagram")
    ax.axis("off")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Chord diagram saved to {output_path}")

    return ax


def plot_alluvial_diagram(
    data: pd.DataFrame,
    stages: List[str],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Create an alluvial diagram (flow diagram) for showing transitions.

    Args:
        data: DataFrame with columns for each stage
        stages: List of column names representing stages
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(data, pd.DataFrame, "data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    try:
        import alluvial
    except ImportError:
        # Fallback implementation
        _plot_simple_alluvial(ax, data, stages, **kwargs)
    else:
        # Use alluvial package
        alluvial.plot(data, stages, ax=ax, **kwargs)

    ax.set_title("Alluvial Diagram")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Alluvial diagram saved to {output_path}")

    return ax


def _plot_simple_alluvial(ax: Axes, data: pd.DataFrame, stages: List[str], **kwargs) -> None:
    """Simple fallback alluvial diagram."""
    n_stages = len(stages)
    colors = plt.cm.tab10(np.linspace(0, 1, len(data)))

    for i, (_, row) in enumerate(data.iterrows()):
        x_positions = np.arange(n_stages)
        y_positions = [row[stage] for stage in stages]

        # Draw flow lines
        ax.plot(x_positions, y_positions, "o-", color=colors[i], linewidth=2, alpha=0.7, markersize=6)

    ax.set_xticks(range(n_stages))
    ax.set_xticklabels(stages)
    ax.set_xlabel("Stages")
    ax.set_ylabel("Values")


def plot_circular_barplot(
    values: np.ndarray,
    labels: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 8),
    **kwargs,
) -> Axes:
    """Create a circular bar plot.

    Args:
        values: Values for each bar
        labels: Optional labels for each bar
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(values, np.ndarray, "values")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": "polar"})

    n_bars = len(values)
    angles = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)

    # Plot bars
    bars = ax.bar(angles, values, width=0.4, bottom=0.0, alpha=0.7, **kwargs)

    # Add labels
    if labels:
        for angle, label in zip(angles, labels):
            ax.text(
                angle, max(values) * 1.1, label, ha="center", va="center", rotation=np.degrees(angle) - 90, fontsize=8
            )

    ax.set_rlabel_position(0)
    ax.set_title("Circular Bar Plot")
    ax.grid(True, alpha=0.3)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Circular bar plot saved to {output_path}")

    return ax


def plot_network_circular_layout(
    graph: Any,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 10),
    **kwargs,
) -> Axes:
    """Plot a network in circular layout.

    Args:
        graph: NetworkX graph object
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx required for network visualization")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": "polar"})

    # Create circular layout
    n_nodes = len(graph.nodes())
    angles = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False)
    pos = {node: (angle, 1) for node, angle in zip(graph.nodes(), angles)}

    # Draw network
    nx.draw_networkx_nodes(graph, pos, node_size=300, alpha=0.8, ax=ax)
    nx.draw_networkx_edges(graph, pos, alpha=0.5, ax=ax)

    # Add labels
    for node, (angle, r) in pos.items():
        ax.text(angle, r + 0.1, str(node), ha="center", va="center", fontsize=8, rotation=np.degrees(angle) - 90)

    ax.set_title("Circular Network Layout")
    ax.grid(False)
    ax.set_rlabel_position(0)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Circular network plot saved to {output_path}")

    return ax


def plot_upset_plot(
    data: Dict[str, set], *, output_path: str | Path | None = None, figsize: Tuple[float, float] = (10, 6), **kwargs
) -> Any:
    """Create an UpSet plot for set intersections.

    Args:
        data: Dictionary mapping set names to sets
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Figure object
    """
    try:
        import upsetplot
    except ImportError:
        raise ImportError("upsetplot required for UpSet plot visualization")

    # Create UpSet data
    from upsetplot import from_contents

    upset_data = from_contents(data)

    # Create plot
    fig = plt.figure(figsize=figsize)
    upsetplot.plot(upset_data, fig=fig, **kwargs)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"UpSet plot saved to {output_path}")

    return fig
