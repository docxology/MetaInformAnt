"""Animated visualization functions for dynamic data exploration.

This module provides animated plotting functions for time series, evolutionary processes,
clustering algorithms, network dynamics, and trajectory analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, List, Tuple, cast

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import PathCollection
from matplotlib.lines import Line2D
from matplotlib.text import Text

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Optional imports with graceful fallbacks
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    np = cast("Any", None)
    HAS_NUMPY = False

try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    nx = None
    HAS_NETWORKX = False


def _padded_limits(low: float, high: float, margin: float = 0.05) -> Tuple[float, float]:
    """Axis limits that contain ``[low, high]`` with a margin of ``margin`` of the span on each side.

    Scaling the limits by 1.1 (``low * 1.1, high * 1.1``) cut off data whenever the lower bound is positive or
    the upper bound negative: the minimum of an all-positive series lay outside the axes.
    """
    span = high - low
    pad = margin * span if span > 0 else max(abs(high), 1.0) * margin
    return low - pad, high + pad


def animate_time_series(
    data: Any, *, interval: int = 200, output_path: str | Path | None = None, **kwargs: Any
) -> Tuple[plt.Figure, FuncAnimation]:
    """Create an animated time series plot.

    Args:
        data: Time series data: a numpy array or list (one row per series, one column per time point; a 1-D input
            is one series), or a pandas DataFrame (one column per series, one row per time point; the column
            names label the lines)
        interval: Animation interval in milliseconds
        output_path: Optional path to save animation (GIF)
        **kwargs: ``figsize`` sizes the figure, ``anim_kwargs`` is passed to FuncAnimation, and any other
            keyword argument is passed to matplotlib plot() for every line

    Returns:
        Tuple of (figure, animation) objects

    Raises:
        ImportError: If required dependencies are missing
    """
    if not HAS_NUMPY:
        raise ImportError("numpy required for animations")

    validation.validate_type(data, (np.ndarray, list, pd.DataFrame), "data")

    anim_kwargs = kwargs.pop("anim_kwargs", {})
    fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 6)))

    if isinstance(data, pd.DataFrame):
        series_labels = [str(column) for column in data.columns]
        data = data.to_numpy(dtype=float).T
    else:
        data = np.asarray(data)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        series_labels = [f"Series {i + 1}" for i in range(data.shape[0])]

    # Setup the plot
    lines: list[Line2D] = []
    for i in range(data.shape[0]):
        (line,) = ax.plot([], [], **{"label": series_labels[i], **kwargs})
        lines.append(line)

    ax.set_xlim(0, data.shape[1])
    ax.set_ylim(*_padded_limits(float(np.nanmin(data)), float(np.nanmax(data))))
    ax.set_xlabel("Time")
    ax.set_ylabel("Value")
    ax.set_title("Time Series Animation")
    if len(lines) > 1:
        ax.legend()

    # Animation function
    def animate(frame: int) -> list[Line2D]:
        for i, line in enumerate(lines):
            x_data = np.arange(frame + 1)
            y_data = data[i, : frame + 1]
            line.set_data(x_data, y_data)
        return lines

    anim = FuncAnimation(fig, animate, frames=data.shape[1], interval=interval, blit=True, **anim_kwargs)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        writer = PillowWriter(fps=cast("int", 1000 / interval))
        anim.save(output_path, writer=writer)
        logger.info(f"Time series animation saved to {output_path}")

    return fig, anim


def animate_evolution(
    sequences: List[str], *, interval: int = 500, output_path: str | Path | None = None, **kwargs: Any
) -> Tuple[plt.Figure, FuncAnimation]:
    """Create an animated sequence evolution visualization.

    The generations are shown one after another; each is revealed one character per frame, so there are
    ``sum(len(sequence))`` frames. A second line lists the first five positions at which the revealed part of the
    current generation differs from generation 1.

    Args:
        sequences: List of sequences showing evolutionary progression
        interval: Animation interval in milliseconds
        output_path: Optional path to save animation (GIF)
        **kwargs: ``figsize`` sizes the figure and ``anim_kwargs`` is passed to FuncAnimation

    Returns:
        Tuple of (figure, animation) objects

    Raises:
        ValueError: If ``sequences`` is empty or contains no characters
    """
    if not HAS_NUMPY:
        raise ImportError("numpy required for animations")

    validation.validate_type(sequences, list, "sequences")

    if not sequences:
        raise ValueError("Sequences list cannot be empty")
    if not any(sequences):
        raise ValueError("Sequences must contain at least one character to animate")

    anim_kwargs = kwargs.pop("anim_kwargs", {})
    fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 8)))

    # Setup the plot: a line for the revealed sequence and a line for its differences from generation 1
    generation_text = ax.text(0.1, 0.9, "", fontsize=12, transform=ax.transAxes, family="monospace")
    mutation_text = ax.text(0.1, 0.8, "", fontsize=12, transform=ax.transAxes, family="monospace")
    text_objects = [generation_text, mutation_text]

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Sequence Evolution Animation")

    # Frame k reveals ``shown`` characters of generation ``generation``: generation 1 first, then 2, and so on
    plan = [(generation, shown) for generation, seq in enumerate(sequences) for shown in range(1, len(seq) + 1)]

    # Animation function
    def animate(frame: int) -> list[Text]:
        generation, shown = plan[frame]
        displayed_seq = sequences[generation][:shown]

        generation_text.set_text(f"Generation {generation + 1}: {displayed_seq}")

        # Show mutations: positions where the revealed part differs from the first sequence
        mutations = []
        if generation > 0:
            for i, (orig, curr) in enumerate(zip(sequences[0][:shown], displayed_seq)):
                if orig != curr:
                    mutations.append(f"{i+1}:{orig}->{curr}")
        mutation_text.set_text(f'Mutations: {", ".join(mutations[:5])}' if mutations else "")

        return text_objects

    anim = FuncAnimation(
        fig,
        animate,
        frames=len(plan),
        interval=interval,
        blit=False,  # Text animation doesn't work well with blit
        **anim_kwargs,
    )

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        writer = PillowWriter(fps=cast("int", 1000 / interval))
        anim.save(output_path, writer=writer)
        logger.info(f"Sequence evolution animation saved to {output_path}")

    return fig, anim


def animate_clustering(
    data: Any,
    cluster_labels_over_time: List[List[int]],
    *,
    interval: int = 300,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Tuple[plt.Figure, FuncAnimation]:
    """Create an animated clustering process visualization.

    Every frame recolours the points by that frame's cluster labels. The colour of a label is the same in every
    frame (the colour scale spans the smallest to the largest label of the whole run).

    Args:
        data: 2D data points to cluster
        cluster_labels_over_time: One sequence of integer cluster labels per frame, one label per data point
        interval: Animation interval in milliseconds
        output_path: Optional path to save animation (GIF)
        **kwargs: ``figsize`` sizes the figure and ``anim_kwargs`` is passed to FuncAnimation; any other keyword
            argument (``cmap``, ``alpha``, ``s``, ...) is passed to matplotlib scatter()

    Returns:
        Tuple of (figure, animation) objects

    Raises:
        ValueError: If data is not 2D, ``cluster_labels_over_time`` is empty, or a frame does not have one label
            per data point
    """
    if not HAS_NUMPY:
        raise ImportError("numpy required for animations")

    validation.validate_type(data, (np.ndarray, list), "data")
    validation.validate_type(cluster_labels_over_time, list, "cluster_labels_over_time")

    if isinstance(data, list):
        data = np.array(data)

    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError("Data must be 2D for clustering animation")

    if not cluster_labels_over_time:
        raise ValueError("Cluster labels over time cannot be empty")

    frames = [np.asarray(labels) for labels in cluster_labels_over_time]
    for frame_index, labels in enumerate(frames):
        if labels.shape != (len(data),):
            raise ValueError(
                f"Frame {frame_index} has labels of shape {labels.shape}; "
                f"the animation needs one label for each of the {len(data)} data points"
            )

    anim_kwargs = kwargs.pop("anim_kwargs", {})
    fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    # Setup the plot: one colour scale for the whole run, so a label keeps its colour from frame to frame
    all_labels = np.concatenate(frames)
    scatter_kwargs: dict[str, Any] = {
        "cmap": "tab10",
        "alpha": 0.7,
        "vmin": all_labels.min(),
        "vmax": all_labels.max(),
        **kwargs,
    }
    scatter = ax.scatter(data[:, 0], data[:, 1], c=frames[0], **scatter_kwargs)

    ax.set_xlabel("Feature 1")
    ax.set_ylabel("Feature 2")
    ax.set_title("Clustering Animation - Iteration 1")

    # Animation function
    def animate(frame: int) -> list[PathCollection]:
        scatter.set_array(frames[min(frame, len(frames) - 1)])
        ax.set_title(f"Clustering Animation - Iteration {frame + 1}")
        return [scatter]

    anim = FuncAnimation(
        fig,
        animate,
        frames=len(frames),
        interval=interval,
        blit=True,
        **anim_kwargs,
    )

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        writer = PillowWriter(fps=cast("int", 1000 / interval))
        anim.save(output_path, writer=writer)
        logger.info(f"Clustering animation saved to {output_path}")

    return fig, anim


def animate_network(
    graphs_over_time: List[Any], *, interval: int = 500, output_path: str | Path | None = None, **kwargs: Any
) -> Tuple[plt.Figure, FuncAnimation]:
    """Create an animated network evolution visualization.

    Args:
        graphs_over_time: List of NetworkX graphs showing network evolution; later graphs may add nodes
        interval: Animation interval in milliseconds
        output_path: Optional path to save animation (GIF)
        **kwargs: ``figsize`` sizes the figure and ``anim_kwargs`` is passed to FuncAnimation

    Returns:
        Tuple of (figure, animation) objects

    Raises:
        ImportError: If NetworkX is not available
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for network animations")

    validation.validate_type(graphs_over_time, list, "graphs_over_time")

    if not graphs_over_time:
        raise ValueError("Graphs over time cannot be empty")

    anim_kwargs = kwargs.pop("anim_kwargs", {})
    fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    # Fixed positions for consistent layout, over every node that appears in any frame (a node added later
    # would otherwise have no position and nx.draw would raise)
    pos = nx.spring_layout(nx.compose_all(graphs_over_time), seed=42)

    # Setup initial plot
    G_initial = graphs_over_time[0]
    nx.draw_networkx_nodes(G_initial, pos, ax=ax, node_color="lightblue", node_size=300, alpha=0.8)
    nx.draw_networkx_edges(G_initial, pos, ax=ax, edge_color="gray", width=1, alpha=0.6)
    nx.draw_networkx_labels(G_initial, pos, ax=ax, font_size=8)

    ax.set_title("Network Evolution - Time 0")
    ax.axis("off")

    # Animation function
    def animate(frame: int) -> list[Any]:
        G = graphs_over_time[min(frame, len(graphs_over_time) - 1)]

        # Clear previous elements
        ax.clear()
        ax.axis("off")

        # Draw current network
        nx.draw(
            G,
            pos,
            ax=ax,
            with_labels=True,
            node_color="lightblue",
            node_size=300,
            edge_color="gray",
            width=1,
            alpha=0.8,
            font_size=8,
        )

        ax.set_title(f"Network Evolution - Time {frame + 1}")

        return []

    anim = FuncAnimation(
        fig,
        animate,
        frames=len(graphs_over_time),
        interval=interval,
        blit=False,  # Network drawing doesn't work well with blit
        **anim_kwargs,
    )

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        writer = PillowWriter(fps=cast("int", 1000 / interval))
        anim.save(output_path, writer=writer)
        logger.info(f"Network animation saved to {output_path}")

    return fig, anim


def animate_trajectory(
    trajectories: List[Any], *, interval: int = 200, output_path: str | Path | None = None, **kwargs: Any
) -> Tuple[plt.Figure, FuncAnimation]:
    """Create an animated trajectory visualization.

    Args:
        trajectories: List of trajectory data (each trajectory is a 2D array of positions)
        interval: Animation interval in milliseconds
        output_path: Optional path to save animation (GIF)
        **kwargs: Additional arguments for animation

    Returns:
        Tuple of (figure, animation) objects

    Raises:
        ValueError: If trajectory data is malformed
    """
    if not HAS_NUMPY:
        raise ImportError("numpy required for animations")

    validation.validate_type(trajectories, list, "trajectories")

    if not trajectories:
        raise ValueError("Trajectories list cannot be empty")

    # Convert trajectories to numpy arrays
    traj_arrays = []
    for traj in trajectories:
        if isinstance(traj, list):
            traj = np.array(traj)
        if traj.ndim != 2 or traj.shape[1] != 2:
            raise ValueError("Each trajectory must be a 2D array with shape (n_points, 2)")
        traj_arrays.append(traj)

    anim_kwargs = kwargs.pop("anim_kwargs", {})
    fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    # Setup the plot
    colors = plt.cm.tab10(np.linspace(0, 1, len(traj_arrays)))

    # Plot complete trajectories as faint lines
    for i, traj in enumerate(traj_arrays):
        ax.plot(traj[:, 0], traj[:, 1], color=colors[i], alpha=0.2, linewidth=1)

    # Plot animated points
    points: list[Line2D] = []
    trails: list[Line2D] = []
    for i, color in enumerate(colors):
        (point,) = ax.plot([], [], "o", color=color, markersize=8, alpha=0.9)
        (trail,) = ax.plot([], [], "-", color=color, alpha=0.6, linewidth=2)
        points.append(point)
        trails.append(trail)

    # Set axis limits
    all_x = np.concatenate([traj[:, 0] for traj in traj_arrays])
    all_y = np.concatenate([traj[:, 1] for traj in traj_arrays])
    ax.set_xlim(*_padded_limits(float(all_x.min()), float(all_x.max())))
    ax.set_ylim(*_padded_limits(float(all_y.min()), float(all_y.max())))

    ax.set_xlabel("X Position")
    ax.set_ylabel("Y Position")
    ax.set_title("Trajectory Animation")
    ax.grid(True, alpha=0.3)

    # Animation function
    max_frames = max(len(traj) for traj in traj_arrays)

    def animate(frame: int) -> list[Line2D]:
        for i, (point, trail, traj) in enumerate(zip(points, trails, traj_arrays)):
            # Show point up to current frame
            frame_idx = min(frame, len(traj) - 1)
            point.set_data([traj[frame_idx, 0]], [traj[frame_idx, 1]])

            # Show trail up to current frame
            trail_data = traj[: frame_idx + 1]
            trail.set_data(trail_data[:, 0], trail_data[:, 1])

        return points + trails

    anim = FuncAnimation(fig, animate, frames=max_frames, interval=interval, blit=True, **anim_kwargs)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        writer = PillowWriter(fps=cast("int", 1000 / interval))
        anim.save(output_path, writer=writer)
        logger.info(f"Trajectory animation saved to {output_path}")

    return fig, anim
