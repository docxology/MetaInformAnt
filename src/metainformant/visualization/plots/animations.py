"""Animated visualization functions for dynamic data exploration.

This module provides animated plotting functions for time series, evolutionary processes,
clustering algorithms, network dynamics, and trajectory analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, List, Tuple

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.axes import Axes

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Optional imports with graceful fallbacks
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    np = None
    HAS_NUMPY = False

try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    nx = None
    HAS_NETWORKX = False


def animate_time_series(
    data: Any, *, interval: int = 200, output_path: str | Path | None = None, **kwargs
) -> Tuple[plt.Figure, FuncAnimation]:
    """Create an animated time series plot.

    Args:
        data: Time series data (pandas DataFrame or numpy array)
        interval: Animation interval in milliseconds
        output_path: Optional path to save animation (GIF)
        **kwargs: Additional arguments for animation

    Returns:
        Tuple of (figure, animation) objects

    Raises:
        ImportError: If required dependencies are missing
    """
    if not HAS_NUMPY:
        raise ImportError("numpy required for animations")

    validation.validate_type(data, (np.ndarray, list), "data")

    fig, ax = plt.subplots(figsize=kwargs.get("figsize", (10, 6)))

    if isinstance(data, list):
        data = np.array(data)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    # Setup the plot
    lines = []
    for i in range(data.shape[0]):
        (line,) = ax.plot([], [], label=f"Series {i+1}", **kwargs)
        lines.append(line)

    ax.set_xlim(0, data.shape[1])
    ax.set_ylim(data.min() * 1.1, data.max() * 1.1)
    ax.set_xlabel("Time")
    ax.set_ylabel("Value")
    ax.set_title("Time Series Animation")
    if len(lines) > 1:
        ax.legend()

    # Animation function
    def animate(frame):
        for i, line in enumerate(lines):
            x_data = np.arange(frame + 1)
            y_data = data[i, : frame + 1]
            line.set_data(x_data, y_data)
        return lines

    anim = FuncAnimation(
        fig, animate, frames=data.shape[1], interval=interval, blit=True, **kwargs.get("anim_kwargs", {})
    )

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        writer = PillowWriter(fps=1000 / interval)
        anim.save(output_path, writer=writer)
        logger.info(f"Time series animation saved to {output_path}")

    return fig, anim


def animate_evolution(
    sequences: List[str], *, interval: int = 500, output_path: str | Path | None = None, **kwargs
) -> Tuple[plt.Figure, FuncAnimation]:
    """Create an animated sequence evolution visualization.

    Args:
        sequences: List of sequences showing evolutionary progression
        interval: Animation interval in milliseconds
        output_path: Optional path to save animation (GIF)
        **kwargs: Additional arguments for animation

    Returns:
        Tuple of (figure, animation) objects
    """
    if not HAS_NUMPY:
        raise ImportError("numpy required for animations")

    validation.validate_type(sequences, list, "sequences")

    if not sequences:
        raise ValueError("Sequences list cannot be empty")

    fig, ax = plt.subplots(figsize=kwargs.get("figsize", (12, 8)))

    # Setup the plot
    text_objects = []
    for i, seq in enumerate(sequences):
        text = ax.text(0.1, 0.9 - i * 0.1, "", fontsize=12, transform=ax.transAxes, family="monospace")
        text_objects.append(text)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Sequence Evolution Animation")

    # Animation function
    def animate(frame):
        current_seq_idx = min(frame, len(sequences) - 1)
        seq = sequences[current_seq_idx]

        # Display sequence with gradual revelation
        display_length = min(len(seq), frame - current_seq_idx * len(seq) + 1)
        displayed_seq = seq[:display_length] if display_length > 0 else ""

        text_objects[0].set_text(f"Generation {current_seq_idx + 1}: {displayed_seq}")

        # Show mutations (simplified - highlight differences from first sequence)
        if current_seq_idx > 0 and displayed_seq:
            original = sequences[0][:display_length]
            mutations = []
            for i, (orig, curr) in enumerate(zip(original, displayed_seq)):
                if orig != curr:
                    mutations.append(f"{i+1}:{orig}->{curr}")

            if mutations:
                mutation_text = f'Mutations: {", ".join(mutations[:5])}'  # Show first 5
                if len(text_objects) > 1:
                    text_objects[1].set_text(mutation_text)

        return text_objects

    total_frames = sum(len(seq) for seq in sequences)
    anim = FuncAnimation(
        fig,
        animate,
        frames=total_frames,
        interval=interval,
        blit=False,  # Text animation doesn't work well with blit
        **kwargs.get("anim_kwargs", {}),
    )

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        writer = PillowWriter(fps=1000 / interval)
        anim.save(output_path, writer=writer)
        logger.info(f"Sequence evolution animation saved to {output_path}")

    return fig, anim


def animate_clustering(
    data: Any,
    cluster_labels_over_time: List[List[int]],
    *,
    interval: int = 300,
    output_path: str | Path | None = None,
    **kwargs,
) -> Tuple[plt.Figure, FuncAnimation]:
    """Create an animated clustering process visualization.

    Args:
        data: 2D data points to cluster
        cluster_labels_over_time: List of cluster assignments for each frame
        interval: Animation interval in milliseconds
        output_path: Optional path to save animation (GIF)
        **kwargs: Additional arguments for animation

    Returns:
        Tuple of (figure, animation) objects

    Raises:
        ValueError: If data dimensions don't match cluster assignments
    """
    if not HAS_NUMPY:
        raise ImportError("numpy required for animations")

    validation.validate_type(data, (np.ndarray, list), "data")
    validation.validate_type(cluster_labels_over_time, list, "cluster_labels_over_time")

    if isinstance(data, list):
        data = np.array(data)

    if data.shape[1] != 2:
        raise ValueError("Data must be 2D for clustering animation")

    if not cluster_labels_over_time:
        raise ValueError("Cluster labels over time cannot be empty")

    fig, ax = plt.subplots(figsize=kwargs.get("figsize", (10, 8)))

    # Setup the plot
    scatter = ax.scatter(data[:, 0], data[:, 1], c=cluster_labels_over_time[0], cmap="tab10", alpha=0.7, **kwargs)

    ax.set_xlabel("Feature 1")
    ax.set_ylabel("Feature 2")
    ax.set_title("Clustering Animation - Iteration 0")

    # Animation function
    def animate(frame):
        labels = cluster_labels_over_time[min(frame, len(cluster_labels_over_time) - 1)]
        scatter.set_color(plt.cm.tab10(labels / max(labels) if max(labels) > 0 else 1))
        ax.set_title(f"Clustering Animation - Iteration {frame + 1}")
        return [scatter]

    anim = FuncAnimation(
        fig,
        animate,
        frames=len(cluster_labels_over_time),
        interval=interval,
        blit=True,
        **kwargs.get("anim_kwargs", {}),
    )

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        writer = PillowWriter(fps=1000 / interval)
        anim.save(output_path, writer=writer)
        logger.info(f"Clustering animation saved to {output_path}")

    return fig, anim


def animate_network(
    graphs_over_time: List[Any], *, interval: int = 500, output_path: str | Path | None = None, **kwargs
) -> Tuple[plt.Figure, FuncAnimation]:
    """Create an animated network evolution visualization.

    Args:
        graphs_over_time: List of NetworkX graphs showing network evolution
        interval: Animation interval in milliseconds
        output_path: Optional path to save animation (GIF)
        **kwargs: Additional arguments for animation

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

    fig, ax = plt.subplots(figsize=kwargs.get("figsize", (10, 8)))

    # Get all possible nodes across all time points
    all_nodes = set()
    for G in graphs_over_time:
        all_nodes.update(G.nodes())
    all_nodes = sorted(all_nodes)

    # Fixed positions for consistent layout
    pos = nx.spring_layout(graphs_over_time[0], seed=42)

    # Setup initial plot
    G_initial = graphs_over_time[0]
    nodes = nx.draw_networkx_nodes(G_initial, pos, ax=ax, node_color="lightblue", node_size=300, alpha=0.8)
    edges = nx.draw_networkx_edges(G_initial, pos, ax=ax, edge_color="gray", width=1, alpha=0.6)
    labels = nx.draw_networkx_labels(G_initial, pos, ax=ax, font_size=8)

    ax.set_title("Network Evolution - Time 0")
    ax.axis("off")

    # Animation function
    def animate(frame):
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
        **kwargs.get("anim_kwargs", {}),
    )

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        writer = PillowWriter(fps=1000 / interval)
        anim.save(output_path, writer=writer)
        logger.info(f"Network animation saved to {output_path}")

    return fig, anim


def animate_trajectory(
    trajectories: List[Any], *, interval: int = 200, output_path: str | Path | None = None, **kwargs
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

    fig, ax = plt.subplots(figsize=kwargs.get("figsize", (10, 8)))

    # Setup the plot
    colors = plt.cm.tab10(np.linspace(0, 1, len(traj_arrays)))

    # Plot complete trajectories as faint lines
    for i, traj in enumerate(traj_arrays):
        ax.plot(traj[:, 0], traj[:, 1], color=colors[i], alpha=0.2, linewidth=1)

    # Plot animated points
    points = []
    trails = []
    for i, color in enumerate(colors):
        (point,) = ax.plot([], [], "o", color=color, markersize=8, alpha=0.9)
        (trail,) = ax.plot([], [], "-", color=color, alpha=0.6, linewidth=2)
        points.append(point)
        trails.append(trail)

    # Set axis limits
    all_x = np.concatenate([traj[:, 0] for traj in traj_arrays])
    all_y = np.concatenate([traj[:, 1] for traj in traj_arrays])
    ax.set_xlim(all_x.min() * 1.1, all_x.max() * 1.1)
    ax.set_ylim(all_y.min() * 1.1, all_y.max() * 1.1)

    ax.set_xlabel("X Position")
    ax.set_ylabel("Y Position")
    ax.set_title("Trajectory Animation")
    ax.grid(True, alpha=0.3)

    # Animation function
    max_frames = max(len(traj) for traj in traj_arrays)

    def animate(frame):
        for i, (point, trail, traj) in enumerate(zip(points, trails, traj_arrays)):
            # Show point up to current frame
            frame_idx = min(frame, len(traj) - 1)
            point.set_data([traj[frame_idx, 0]], [traj[frame_idx, 1]])

            # Show trail up to current frame
            trail_data = traj[: frame_idx + 1]
            trail.set_data(trail_data[:, 0], trail_data[:, 1])

        return points + trails

    anim = FuncAnimation(fig, animate, frames=max_frames, interval=interval, blit=True, **kwargs.get("anim_kwargs", {}))

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        writer = PillowWriter(fps=1000 / interval)
        anim.save(output_path, writer=writer)
        logger.info(f"Trajectory animation saved to {output_path}")

    return fig, anim
