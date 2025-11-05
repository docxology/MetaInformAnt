"""Animation functions for dynamic visualizations.

This module provides animation functions for time series, evolutionary processes,
clustering iterations, network evolution, and trajectory inference.
"""

from __future__ import annotations

from typing import Sequence

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Non-interactive for tests/headless
matplotlib.use("Agg", force=True)
from matplotlib.animation import FuncAnimation  # noqa: E402


def animate_time_series(
    series: Sequence[float],
    *,
    interval_ms: int = 100,
    repeat: bool = False,
    init_points: int = 1,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, FuncAnimation]:
    """Animate a time series as a growing line.

    Args:
        series: Time series data
        interval_ms: Animation interval in milliseconds
        repeat: Whether to repeat animation
        init_points: Number of initial points to show
        ax: Matplotlib axes (creates new if None)

    Returns:
        Tuple of (figure, animation)

    Example:
        >>> from metainformant.visualization import animate_time_series
        >>> data = [1, 2, 3, 2, 4, 5, 3, 6]
        >>> fig, anim = animate_time_series(data, interval_ms=200)
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    x = list(range(len(series)))
    (line,) = ax.plot([], [], "-b")
    ax.set_xlim(0, max(1, len(series) - 1))
    ymin = min(series) if series else 0.0
    ymax = max(series) if series else 1.0
    if ymin == ymax:
        ymin -= 1.0
        ymax += 1.0
    ax.set_ylim(ymin, ymax)

    def init():
        line.set_data(x[:init_points], series[:init_points])
        return (line,)

    def update(frame_idx: int):
        i = min(frame_idx + 1, len(series))
        line.set_data(x[:i], series[:i])
        return (line,)

    anim = FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=len(series),
        interval=interval_ms,
        blit=True,
        repeat=repeat,
    )
    return fig, anim


def animate_evolution(
    generations: Sequence[dict[str, Sequence[float]]],
    *,
    interval_ms: int = 200,
    ax: plt.Axes | None = None,
    **kwargs
) -> tuple[plt.Figure, FuncAnimation]:
    """Animate evolutionary process over generations.

    Args:
        generations: List of dictionaries with generation data
        interval_ms: Animation interval in milliseconds
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments

    Returns:
        Tuple of (figure, animation)

    Example:
        >>> from metainformant.visualization import animate_evolution
        >>> generations = [
        ...     {'fitness': [0.5, 0.6, 0.4]},
        ...     {'fitness': [0.6, 0.7, 0.5]},
        ...     {'fitness': [0.7, 0.8, 0.6]}
        ... ]
        >>> fig, anim = animate_evolution(generations)
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    # Extract fitness values
    fitness_data = [gen.get('fitness', []) for gen in generations]
    max_fitness = max(max(fit) if fit else 0 for fit in fitness_data)
    min_fitness = min(min(fit) if fit else 0 for fit in fitness_data)

    ax.set_xlim(0, len(generations))
    ax.set_ylim(min_fitness - 0.1, max_fitness + 0.1)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Fitness")

    (line,) = ax.plot([], [], 'o-', **kwargs)

    def init():
        line.set_data([], [])
        return (line,)

    def update(frame_idx: int):
        if frame_idx < len(generations):
            gen_data = fitness_data[frame_idx]
            if gen_data:
                mean_fitness = np.mean(gen_data)
                line.set_data(list(range(frame_idx + 1)), 
                             [np.mean(fitness_data[i]) if fitness_data[i] else 0 
                              for i in range(frame_idx + 1)])
        return (line,)

    anim = FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=len(generations),
        interval=interval_ms,
        blit=True,
        repeat=False,
    )
    return fig, anim


def animate_clustering(
    iterations: Sequence[dict[str, np.ndarray]],
    *,
    interval_ms: int = 300,
    ax: plt.Axes | None = None,
    **kwargs
) -> tuple[plt.Figure, FuncAnimation]:
    """Animate clustering iterations.

    Args:
        iterations: List of dictionaries with iteration data ('points', 'centers', 'labels')
        interval_ms: Animation interval in milliseconds
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments

    Returns:
        Tuple of (figure, animation)

    Example:
        >>> from metainformant.visualization import animate_clustering
        >>> import numpy as np
        >>> iterations = [
        ...     {'points': np.random.random((50, 2)), 'centers': np.random.random((3, 2))},
        ...     {'points': np.random.random((50, 2)), 'centers': np.random.random((3, 2))}
        ... ]
        >>> fig, anim = animate_clustering(iterations)
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    # Get bounds from all iterations
    all_points = []
    for it in iterations:
        if 'points' in it:
            all_points.append(it['points'])
    if all_points:
        all_points = np.vstack(all_points)
        ax.set_xlim(all_points[:, 0].min() - 0.1, all_points[:, 0].max() + 0.1)
        ax.set_ylim(all_points[:, 1].min() - 0.1, all_points[:, 1].max() + 0.1)

    scatter = ax.scatter([], [], c=[], alpha=0.6, **kwargs)
    centers = ax.scatter([], [], c='red', marker='x', s=100, linewidths=2)

    def init():
        scatter.set_offsets(np.array([]).reshape(0, 2))
        scatter.set_array([])
        centers.set_offsets(np.array([]).reshape(0, 2))
        return scatter, centers

    def update(frame_idx: int):
        if frame_idx < len(iterations):
            it = iterations[frame_idx]
            if 'points' in it:
                scatter.set_offsets(it['points'])
                if 'labels' in it:
                    scatter.set_array(it['labels'])
            if 'centers' in it:
                centers.set_offsets(it['centers'])
        return scatter, centers

    anim = FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=len(iterations),
        interval=interval_ms,
        blit=True,
        repeat=False,
    )
    return fig, anim


def animate_network(
    network_states: Sequence[dict[str, Any]],
    *,
    interval_ms: int = 200,
    ax: plt.Axes | None = None,
    **kwargs
) -> tuple[plt.Figure, FuncAnimation]:
    """Animate network evolution over time.

    Args:
        network_states: List of network state dictionaries
        interval_ms: Animation interval in milliseconds
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments

    Returns:
        Tuple of (figure, animation)

    Example:
        >>> from metainformant.visualization import animate_network
        >>> network_states = [
        ...     {'nodes': ['A', 'B'], 'edges': [('A', 'B')]},
        ...     {'nodes': ['A', 'B', 'C'], 'edges': [('A', 'B'), ('B', 'C')]}
        ... ]
        >>> fig, anim = animate_network(network_states)
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("NetworkX required for network animation")

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    else:
        fig = ax.figure

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.axis('off')

    def init():
        return ()

    def update(frame_idx: int):
        ax.clear()
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.axis('off')

        if frame_idx < len(network_states):
            state = network_states[frame_idx]
            G = nx.Graph()
            if 'nodes' in state:
                G.add_nodes_from(state['nodes'])
            if 'edges' in state:
                G.add_edges_from(state['edges'])

            pos = nx.spring_layout(G, k=0.5, iterations=50)
            nx.draw(G, pos, ax=ax, with_labels=True, node_color='lightblue',
                   node_size=500, font_size=8, **kwargs)

        return ()

    anim = FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=len(network_states),
        interval=interval_ms,
        blit=False,
        repeat=False,
    )
    return fig, anim


def animate_trajectory(
    trajectory_data: Sequence[dict[str, np.ndarray]],
    *,
    interval_ms: int = 200,
    ax: plt.Axes | None = None,
    **kwargs
) -> tuple[plt.Figure, FuncAnimation]:
    """Animate trajectory inference over time.

    Args:
        trajectory_data: List of trajectory state dictionaries
        interval_ms: Animation interval in milliseconds
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments

    Returns:
        Tuple of (figure, animation)

    Example:
        >>> from metainformant.visualization import animate_trajectory
        >>> import numpy as np
        >>> trajectory_data = [
        ...     {'positions': np.random.random((10, 2))},
        ...     {'positions': np.random.random((10, 2))}
        ... ]
        >>> fig, anim = animate_trajectory(trajectory_data)
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    # Get bounds
    all_positions = []
    for td in trajectory_data:
        if 'positions' in td:
            all_positions.append(td['positions'])
    if all_positions:
        all_positions = np.vstack(all_positions)
        ax.set_xlim(all_positions[:, 0].min() - 0.1, all_positions[:, 0].max() + 0.1)
        ax.set_ylim(all_positions[:, 1].min() - 0.1, all_positions[:, 1].max() + 0.1)

    scatter = ax.scatter([], [], alpha=0.6, **kwargs)

    def init():
        scatter.set_offsets(np.array([]).reshape(0, 2))
        return scatter,

    def update(frame_idx: int):
        if frame_idx < len(trajectory_data):
            td = trajectory_data[frame_idx]
            if 'positions' in td:
                scatter.set_offsets(td['positions'])
        return scatter,

    anim = FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=len(trajectory_data),
        interval=interval_ms,
        blit=True,
        repeat=False,
    )
    return fig, anim
