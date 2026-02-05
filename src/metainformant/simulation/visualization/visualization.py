"""Simulation data visualization functions.

This module provides comprehensive visualization capabilities for simulation results
across different biological domains, including sequence evolution, population dynamics,
RNA-seq data generation, and agent-based modeling.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Circle, Rectangle

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


def plot_sequence_evolution(
    sequence_history: List[str],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot sequence evolution over generations.

    Args:
        sequence_history: List of sequences at each generation
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(sequence_history, list, "sequence_history")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    generations = len(sequence_history)
    seq_length = len(sequence_history[0]) if sequence_history else 0

    # Create mutation matrix (1 if position changed from previous generation)
    mutation_matrix = np.zeros((generations, seq_length))

    for gen in range(1, generations):
        for pos in range(seq_length):
            if sequence_history[gen][pos] != sequence_history[gen - 1][pos]:
                mutation_matrix[gen, pos] = 1

    # Plot mutation events
    im = ax.imshow(mutation_matrix, cmap="Reds", aspect="auto", origin="lower")
    ax.set_xlabel("Position")
    ax.set_ylabel("Generation")
    ax.set_title("Sequence Evolution - Mutation Events")

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Mutation (1) / No Change (0)")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Sequence evolution plot saved to {output_path}")

    return ax


def animate_sequence_evolution(
    sequence_history: List[str],
    *,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    interval: int = 500,
    **kwargs,
) -> Tuple[plt.Figure, animation.FuncAnimation]:
    """Create an animation of sequence evolution.

    Args:
        sequence_history: List of sequences at each generation
        output_path: Optional path to save the animation
        figsize: Figure size as (width, height)
        interval: Delay between frames in milliseconds
        **kwargs: Additional arguments for customization

    Returns:
        Tuple of (figure, animation) objects
    """
    validation.validate_type(sequence_history, list, "sequence_history")

    fig, ax = plt.subplots(figsize=figsize)

    generations = len(sequence_history)
    seq_length = len(sequence_history[0]) if sequence_history else 0

    # Create mutation matrix
    mutation_matrix = np.zeros((generations, seq_length))

    def animate(frame):
        ax.clear()

        # Update mutation matrix up to current frame
        for gen in range(1, min(frame + 1, generations)):
            for pos in range(seq_length):
                if sequence_history[gen][pos] != sequence_history[gen - 1][pos]:
                    mutation_matrix[gen, pos] = 1

        # Plot current state
        im = ax.imshow(mutation_matrix[: frame + 1], cmap="Reds", aspect="auto", origin="lower")
        ax.set_xlabel("Position")
        ax.set_ylabel("Generation")
        ax.set_title(f"Sequence Evolution - Generation {frame}")

        return [im]

    anim = animation.FuncAnimation(fig, animate, frames=generations, interval=interval, blit=False)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        if output_path.suffix.lower() == ".gif":
            anim.save(str(output_path), writer="pillow")
        else:
            anim.save(str(output_path.with_suffix(".mp4")), writer="ffmpeg")
        logger.info(f"Sequence evolution animation saved to {output_path}")

    return fig, anim


def plot_rnaseq_simulation_results(
    rnaseq_data: Dict[str, Any],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (14, 10),
    **kwargs,
) -> Axes:
    """Plot RNA-seq simulation results.

    Args:
        rnaseq_data: Dictionary containing RNA-seq simulation results
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(rnaseq_data, dict, "rnaseq_data")

    if ax is None:
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()
    else:
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()

    plot_idx = 0

    # Plot 1: Expression distribution
    if "expression_values" in rnaseq_data:
        expr_values = rnaseq_data["expression_values"]
        if HAS_SEABORN:
            sns.histplot(expr_values, ax=axes[plot_idx], kde=True, alpha=0.7)
        else:
            axes[plot_idx].hist(expr_values, bins=50, alpha=0.7)
        axes[plot_idx].set_xlabel("Expression Level")
        axes[plot_idx].set_ylabel("Frequency")
        axes[plot_idx].set_title("Gene Expression Distribution")
        axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot 2: Library size distribution
    if "library_sizes" in rnaseq_data:
        lib_sizes = rnaseq_data["library_sizes"]
        if HAS_SEABORN:
            sns.histplot(lib_sizes, ax=axes[plot_idx], kde=True, alpha=0.7, color="green")
        else:
            axes[plot_idx].hist(lib_sizes, bins=30, alpha=0.7, color="green")
        axes[plot_idx].set_xlabel("Library Size")
        axes[plot_idx].set_ylabel("Frequency")
        axes[plot_idx].set_title("Library Size Distribution")
        axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot 3: Mean-variance relationship
    if "mean_variance" in rnaseq_data:
        mv_data = rnaseq_data["mean_variance"]
        if "means" in mv_data and "variances" in mv_data:
            axes[plot_idx].scatter(mv_data["means"], mv_data["variances"], alpha=0.6, s=20)
            axes[plot_idx].set_xlabel("Mean Expression")
            axes[plot_idx].set_ylabel("Variance")
            axes[plot_idx].set_title("Mean-Variance Relationship")
            axes[plot_idx].set_xscale("log")
            axes[plot_idx].set_yscale("log")
            axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot 4: GC content vs expression
    if "gc_expression" in rnaseq_data:
        gc_expr = rnaseq_data["gc_expression"]
        if "gc_content" in gc_expr and "expression" in gc_expr:
            axes[plot_idx].scatter(gc_expr["gc_content"], gc_expr["expression"], alpha=0.6, s=20)
            axes[plot_idx].set_xlabel("GC Content")
            axes[plot_idx].set_ylabel("Expression Level")
            axes[plot_idx].set_title("GC Content vs Expression")
            axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot 5: Length vs expression
    if "length_expression" in rnaseq_data:
        len_expr = rnaseq_data["length_expression"]
        if "lengths" in len_expr and "expression" in len_expr:
            axes[plot_idx].scatter(len_expr["lengths"], len_expr["expression"], alpha=0.6, s=20, color="purple")
            axes[plot_idx].set_xlabel("Transcript Length")
            axes[plot_idx].set_ylabel("Expression Level")
            axes[plot_idx].set_title("Length vs Expression")
            axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot 6: Simulation parameters summary
    if "simulation_params" in rnaseq_data:
        params = rnaseq_data["simulation_params"]
        param_names = list(params.keys())
        param_values = [str(v) for v in params.values()]

        # Create text summary
        summary_text = "\n".join([f"{name}: {value}" for name, value in zip(param_names, param_values)])
        axes[plot_idx].text(
            0.1,
            0.9,
            summary_text,
            transform=axes[plot_idx].transAxes,
            verticalalignment="top",
            fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.8),
        )
        axes[plot_idx].set_title("Simulation Parameters")
        axes[plot_idx].set_xlim(0, 1)
        axes[plot_idx].set_ylim(0, 1)
        axes[plot_idx].axis("off")

    plt.tight_layout()

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"RNA-seq simulation results plot saved to {output_path}")

    return axes[0]


def plot_population_dynamics_simulation(
    population_history: List[np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot population dynamics simulation results.

    Args:
        population_history: List of population size arrays over time
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(population_history, list, "population_history")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    time_points = np.arange(len(population_history))

    # Plot population trajectories
    colors = plt.cm.tab10(np.linspace(0, 1, len(population_history[0]) if population_history else 1))

    for i in range(len(population_history[0])):
        trajectory = [pop[i] for pop in population_history]
        ax.plot(time_points, trajectory, color=colors[i], linewidth=2, alpha=0.8, label=f"Population {i+1}")

    ax.set_xlabel("Time Step")
    ax.set_ylabel("Population Size")
    ax.set_title("Population Dynamics Simulation")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Population dynamics simulation plot saved to {output_path}")

    return ax


def plot_agent_based_model_results(
    agent_data: Dict[str, Any],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot agent-based model simulation results.

    Args:
        agent_data: Dictionary containing agent simulation results
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(agent_data, dict, "agent_data")

    if ax is None:
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        axes = axes.flatten()
    else:
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        axes = axes.flatten()

    # Plot 1: Agent positions over time (if available)
    if "agent_positions" in agent_data:
        positions = agent_data["agent_positions"]
        if isinstance(positions, list) and len(positions) > 0:
            # Plot final positions
            final_pos = positions[-1]
            if len(final_pos) == 2:  # 2D positions
                x, y = final_pos
                axes[0].scatter(x, y, alpha=0.7, s=50)
                axes[0].set_xlabel("X Position")
                axes[0].set_ylabel("Y Position")
            else:  # 1D positions
                axes[0].hist(final_pos, bins=20, alpha=0.7)
                axes[0].set_xlabel("Position")
                axes[0].set_ylabel("Agent Count")
            axes[0].set_title("Agent Positions")
            axes[0].grid(True, alpha=0.3)

    # Plot 2: Agent states over time
    if "agent_states" in agent_data:
        states_history = agent_data["agent_states"]
        if isinstance(states_history, list):
            time_points = np.arange(len(states_history))
            # Count agents in each state over time
            state_counts = {}
            for t, states in enumerate(states_history):
                unique_states, counts = np.unique(states, return_counts=True)
                for state, count in zip(unique_states, counts):
                    if state not in state_counts:
                        state_counts[state] = []
                    state_counts[state].append(count)

            for state, counts in state_counts.items():
                axes[1].plot(time_points[: len(counts)], counts, label=f"State {state}", linewidth=2, alpha=0.8)

            axes[1].set_xlabel("Time Step")
            axes[1].set_ylabel("Agent Count")
            axes[1].set_title("Agent States Over Time")
            axes[1].legend()
            axes[1].grid(True, alpha=0.3)

    # Plot 3: Network structure (if applicable)
    if "network_edges" in agent_data:
        edges = agent_data["network_edges"]
        if len(edges) > 0:
            # Simple network visualization
            nodes = set()
            for edge in edges:
                nodes.update(edge)

            # Create adjacency matrix
            node_list = list(nodes)
            n_nodes = len(node_list)
            adj_matrix = np.zeros((n_nodes, n_nodes))

            for edge in edges:
                i = node_list.index(edge[0])
                j = node_list.index(edge[1])
                adj_matrix[i, j] = 1
                adj_matrix[j, i] = 1  # Undirected

            im = axes[2].imshow(adj_matrix, cmap="Blues", aspect="equal")
            axes[2].set_title("Agent Interaction Network")
            axes[2].set_xticks(range(n_nodes))
            axes[2].set_yticks(range(n_nodes))
            axes[2].set_xticklabels([f"A{i}" for i in range(n_nodes)])
            axes[2].set_yticklabels([f"A{i}" for i in range(n_nodes)])
            plt.colorbar(im, ax=axes[2])

    # Plot 4: Simulation summary statistics
    if "simulation_stats" in agent_data:
        stats = agent_data["simulation_stats"]
        stat_names = list(stats.keys())
        stat_values = list(stats.values())

        bars = axes[3].bar(range(len(stat_names)), stat_values, alpha=0.7, color="green")
        axes[3].set_xticks(range(len(stat_names)))
        axes[3].set_xticklabels(stat_names, rotation=45, ha="right")
        axes[3].set_ylabel("Value")
        axes[3].set_title("Simulation Statistics")
        axes[3].grid(True, alpha=0.3, axis="y")

        # Add value labels
        for bar, value in zip(bars, stat_values):
            axes[3].text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.01,
                f"{value:.2f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    plt.tight_layout()

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Agent-based model results plot saved to {output_path}")

    return axes[0]


def plot_evolutionary_simulation_summary(
    evolution_data: Dict[str, Any],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (14, 10),
    **kwargs,
) -> Axes:
    """Plot comprehensive evolutionary simulation summary.

    Args:
        evolution_data: Dictionary containing evolutionary simulation results
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(evolution_data, dict, "evolution_data")

    if ax is None:
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()
    else:
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()

    plot_idx = 0

    # Plot 1: Fitness over generations
    if "fitness_history" in evolution_data:
        fitness_data = evolution_data["fitness_history"]
        if "generations" in fitness_data and "mean_fitness" in fitness_data:
            axes[plot_idx].plot(
                fitness_data["generations"],
                fitness_data["mean_fitness"],
                "b-",
                linewidth=2,
                label="Mean Fitness",
                alpha=0.8,
            )
            if "max_fitness" in fitness_data:
                axes[plot_idx].plot(
                    fitness_data["generations"],
                    fitness_data["max_fitness"],
                    "r--",
                    linewidth=1,
                    label="Max Fitness",
                    alpha=0.7,
                )
            axes[plot_idx].set_xlabel("Generation")
            axes[plot_idx].set_ylabel("Fitness")
            axes[plot_idx].set_title("Fitness Evolution")
            axes[plot_idx].legend()
            axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot 2: Genetic diversity over time
    if "diversity_history" in evolution_data:
        diversity_data = evolution_data["diversity_history"]
        if "generations" in diversity_data and "diversity" in diversity_data:
            axes[plot_idx].plot(
                diversity_data["generations"], diversity_data["diversity"], "g-", linewidth=2, alpha=0.8
            )
            axes[plot_idx].set_xlabel("Generation")
            axes[plot_idx].set_ylabel("Genetic Diversity")
            axes[plot_idx].set_title("Genetic Diversity")
            axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot 3: Allele frequency changes
    if "allele_frequencies" in evolution_data:
        freq_data = evolution_data["allele_frequencies"]
        if isinstance(freq_data, dict) and "generations" in freq_data:
            gens = freq_data["generations"]
            for allele, freqs in freq_data.items():
                if allele != "generations":
                    axes[plot_idx].plot(gens, freqs, label=f"Allele {allele}", alpha=0.7)

            axes[plot_idx].set_xlabel("Generation")
            axes[plot_idx].set_ylabel("Allele Frequency")
            axes[plot_idx].set_title("Allele Frequency Changes")
            axes[plot_idx].legend()
            axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot 4: Population size changes
    if "population_size_history" in evolution_data:
        pop_data = evolution_data["population_size_history"]
        if "generations" in pop_data and "population_size" in pop_data:
            axes[plot_idx].plot(pop_data["generations"], pop_data["population_size"], "purple", linewidth=2, alpha=0.8)
            axes[plot_idx].set_xlabel("Generation")
            axes[plot_idx].set_ylabel("Population Size")
            axes[plot_idx].set_title("Population Dynamics")
            axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot 5: Mutation accumulation
    if "mutation_history" in evolution_data:
        mut_data = evolution_data["mutation_history"]
        if "generations" in mut_data and "mutation_count" in mut_data:
            axes[plot_idx].plot(mut_data["generations"], mut_data["mutation_count"], "orange", linewidth=2, alpha=0.8)
            axes[plot_idx].set_xlabel("Generation")
            axes[plot_idx].set_ylabel("Mutation Count")
            axes[plot_idx].set_title("Mutation Accumulation")
            axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot 6: Final simulation statistics
    if "final_stats" in evolution_data:
        stats = evolution_data["final_stats"]
        stat_names = list(stats.keys())
        stat_values = list(stats.values())

        bars = axes[plot_idx].bar(range(len(stat_names)), stat_values, alpha=0.7, color="brown")
        axes[plot_idx].set_xticks(range(len(stat_names)))
        axes[plot_idx].set_xticklabels(stat_names, rotation=45, ha="right")
        axes[plot_idx].set_ylabel("Value")
        axes[plot_idx].set_title("Final Statistics")
        axes[plot_idx].grid(True, alpha=0.3, axis="y")

        # Add value labels
        for bar, value in zip(bars, stat_values):
            axes[plot_idx].text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.01,
                f"{value:.2f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    plt.tight_layout()

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Evolutionary simulation summary plot saved to {output_path}")

    return axes[0]


def plot_simulation_parameter_sensitivity(
    sensitivity_data: Dict[str, Any],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes:
    """Plot parameter sensitivity analysis results.

    Args:
        sensitivity_data: Dictionary containing sensitivity analysis results
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(sensitivity_data, dict, "sensitivity_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if "parameter_values" in sensitivity_data and "output_values" in sensitivity_data:
        param_values = sensitivity_data["parameter_values"]
        output_values = sensitivity_data["output_values"]

        ax.scatter(param_values, output_values, alpha=0.7, s=50, color="blue")
        ax.plot(param_values, output_values, "b-", alpha=0.5, linewidth=1)

        ax.set_xlabel("Parameter Value")
        ax.set_ylabel("Output Value")
        ax.set_title("Parameter Sensitivity Analysis")
        ax.grid(True, alpha=0.3)

        # Add correlation coefficient
        corr = np.corrcoef(param_values, output_values)[0, 1]
        ax.text(
            0.02,
            0.98,
            f"Correlation: {corr:.3f}",
            transform=ax.transAxes,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
        )

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Simulation parameter sensitivity plot saved to {output_path}")

    return ax


def animate_population_dynamics(
    population_history: List[np.ndarray],
    *,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    interval: int = 200,
    **kwargs,
) -> Tuple[plt.Figure, animation.FuncAnimation]:
    """Create an animation of population dynamics.

    Args:
        population_history: List of population arrays over time
        output_path: Optional path to save the animation
        figsize: Figure size as (width, height)
        interval: Delay between frames in milliseconds
        **kwargs: Additional arguments for customization

    Returns:
        Tuple of (figure, animation) objects
    """
    validation.validate_type(population_history, list, "population_history")

    fig, ax = plt.subplots(figsize=figsize)

    time_points = np.arange(len(population_history))
    colors = plt.cm.tab10(np.linspace(0, 1, len(population_history[0]) if population_history else 1))

    lines = []
    for i in range(len(population_history[0])):
        (line,) = ax.plot([], [], color=colors[i], linewidth=2, alpha=0.8, label=f"Population {i+1}")
        lines.append(line)

    ax.set_xlim(0, len(population_history))
    ax.set_ylim(0, max(max(pop) for pop in population_history) * 1.1)
    ax.set_xlabel("Time Step")
    ax.set_ylabel("Population Size")
    ax.set_title("Population Dynamics Animation")
    ax.legend()
    ax.grid(True, alpha=0.3)

    def animate(frame):
        for i, line in enumerate(lines):
            x_data = time_points[: frame + 1]
            y_data = [pop[i] for pop in population_history[: frame + 1]]
            line.set_data(x_data, y_data)

        ax.set_title(f"Population Dynamics - Time Step {frame}")
        return lines

    anim = animation.FuncAnimation(fig, animate, frames=len(population_history), interval=interval, blit=True)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        if output_path.suffix.lower() == ".gif":
            anim.save(str(output_path), writer="pillow")
        else:
            anim.save(str(output_path.with_suffix(".mp4")), writer="ffmpeg")
        logger.info(f"Population dynamics animation saved to {output_path}")

    return fig, anim


def plot_simulation_validation_comparison(
    observed_data: np.ndarray,
    simulated_data: List[np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot comparison between observed and simulated data for validation.

    Args:
        observed_data: Observed data values
        simulated_data: List of simulated data replicates
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(observed_data, np.ndarray, "observed_data")
    validation.validate_type(simulated_data, list, "simulated_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot observed data
    ax.plot(observed_data, "k-", linewidth=3, label="Observed", alpha=0.8)

    # Plot simulated data (mean and confidence interval)
    sim_array = np.array(simulated_data)
    sim_mean = np.mean(sim_array, axis=0)
    sim_std = np.std(sim_array, axis=0)

    ax.plot(sim_mean, "r-", linewidth=2, label="Simulated Mean", alpha=0.8)
    ax.fill_between(
        range(len(sim_mean)), sim_mean - sim_std, sim_mean + sim_std, alpha=0.3, color="red", label="Simulated ±1σ"
    )

    ax.set_xlabel("Index")
    ax.set_ylabel("Value")
    ax.set_title("Simulation Validation Comparison")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Simulation validation comparison plot saved to {output_path}")

    return ax


def create_interactive_simulation_dashboard(
    simulation_results: Dict[str, Any], *, output_path: str | Path | None = None, **kwargs
) -> Any:
    """Create an interactive simulation results dashboard.

    Args:
        simulation_results: Dictionary containing simulation results
        output_path: Optional path to save the HTML file
        **kwargs: Additional arguments for Plotly customization

    Returns:
        Plotly Figure object
    """
    if not HAS_PLOTLY:
        raise ImportError("Plotly required for interactive simulation dashboard")

    validation.validate_type(simulation_results, dict, "simulation_results")

    # Create subplot figure
    from plotly.subplots import make_subplots

    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=["Time Series", "Distributions", "Correlations", "Statistics"],
        specs=[[{"type": "scatter"}, {"type": "histogram"}], [{"type": "heatmap"}, {"type": "bar"}]],
    )

    # Add time series plot
    if "time_series" in simulation_results:
        ts_data = simulation_results["time_series"]
        if "time" in ts_data and "values" in ts_data:
            fig.add_trace(go.Scatter(x=ts_data["time"], y=ts_data["values"], mode="lines+markers"), row=1, col=1)

    # Add distribution histogram
    if "distributions" in simulation_results:
        dist_data = simulation_results["distributions"]
        fig.add_trace(go.Histogram(x=dist_data, nbinsx=30), row=1, col=2)

    # Add correlation heatmap
    if "correlations" in simulation_results:
        corr_data = simulation_results["correlations"]
        fig.add_trace(go.Heatmap(z=corr_data, colorscale="RdBu", zmid=0), row=2, col=1)

    # Add statistics bar plot
    if "statistics" in simulation_results:
        stats = simulation_results["statistics"]
        fig.add_trace(go.Bar(x=list(stats.keys()), y=list(stats.values())), row=2, col=2)

    fig.update_layout(title="Interactive Simulation Results Dashboard", **kwargs)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        html_path = Path(output_path).with_suffix(".html")
        fig.write_html(str(html_path))
        logger.info(f"Interactive simulation dashboard saved to {html_path}")

    return fig
