"""Epigenome data visualization functions.

This module provides comprehensive visualization capabilities for epigenetic data,
including DNA methylation, histone modifications, chromatin accessibility,
and genome-wide epigenetic patterns.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

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


def plot_methylation_profile(
    methylation_data: np.ndarray,
    positions: np.ndarray | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes:
    """Plot DNA methylation profile along genomic positions.

    Args:
        methylation_data: Methylation levels (0-1) for each CpG site
        positions: Genomic positions of CpG sites
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(methylation_data, np.ndarray, "methylation_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    x_data = positions if positions is not None else np.arange(len(methylation_data))

    # Plot methylation levels
    ax.scatter(x_data, methylation_data, alpha=0.7, s=20, c=methylation_data, cmap="RdYlBu_r", vmin=0, vmax=1)

    # Add smoothed trend line
    if len(methylation_data) > 10:
        from scipy.ndimage import uniform_filter1d

        smoothed = uniform_filter1d(methylation_data.astype(float), size=min(20, len(methylation_data) // 5))
        ax.plot(x_data, smoothed, "k-", linewidth=2, alpha=0.8, label="Smoothed")

    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("Methylation Level")
    ax.set_title("DNA Methylation Profile")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap="RdYlBu_r", norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("Methylation Level")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Methylation profile saved to {output_path}")

    return ax


def plot_chipseq_peaks(
    peaks: List[Dict[str, Any]],
    chromosome: str,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 4),
    **kwargs,
) -> Axes:
    """Plot ChIP-seq peaks along a chromosome.

    Args:
        peaks: List of peak dictionaries with 'start', 'end', 'score' keys
        chromosome: Chromosome name for labeling
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(peaks, list, "peaks")
    validation.validate_type(chromosome, str, "chromosome")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Extract peak data
    starts = [peak["start"] for peak in peaks]
    ends = [peak["end"] for peak in peaks]
    scores = [peak.get("score", 1.0) for peak in peaks]

    # Create rectangles for peaks
    patches = []
    colors = []
    for start, end, score in zip(starts, ends, scores):
        width = end - start
        rect = Rectangle((start, 0), width, score)
        patches.append(rect)
        colors.append(score)

    # Create patch collection
    collection = PatchCollection(patches, cmap="Reds", alpha=0.8)
    collection.set_array(np.array(colors))
    ax.add_collection(collection)

    # Set axis limits
    if peaks:
        min_pos = min(starts)
        max_pos = max(ends)
        ax.set_xlim(min_pos - (max_pos - min_pos) * 0.05, max_pos + (max_pos - min_pos) * 0.05)
        ax.set_ylim(0, max(scores) * 1.1)

    ax.set_xlabel(f"Position on {chromosome}")
    ax.set_ylabel("Peak Score")
    ax.set_title(f"ChIP-seq Peaks - {chromosome}")
    ax.grid(True, alpha=0.3)

    # Add colorbar
    if peaks:
        cbar = plt.colorbar(collection, ax=ax)
        cbar.set_label("Peak Score")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"ChIP-seq peaks plot saved to {output_path}")

    return ax


def plot_atacseq_signal(
    signal_data: np.ndarray,
    positions: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes:
    """Plot ATAC-seq accessibility signal.

    Args:
        signal_data: Accessibility signal values
        positions: Genomic positions
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(signal_data, np.ndarray, "signal_data")
    validation.validate_type(positions, np.ndarray, "positions")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot accessibility signal
    ax.fill_between(positions, signal_data, alpha=0.7, color="skyblue")
    ax.plot(positions, signal_data, "b-", linewidth=1, alpha=0.8)

    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("Accessibility Signal")
    ax.set_title("ATAC-seq Chromatin Accessibility")
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"ATAC-seq signal plot saved to {output_path}")

    return ax


def plot_histone_modification_heatmap(
    histone_data: Dict[str, np.ndarray],
    positions: np.ndarray | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot heatmap of multiple histone modifications.

    Args:
        histone_data: Dictionary mapping modification names to signal arrays
        positions: Genomic positions (optional)
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(histone_data, dict, "histone_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Convert to matrix
    modifications = list(histone_data.keys())
    signals = np.array([histone_data[mod] for mod in modifications])

    # Plot heatmap
    im = ax.imshow(signals, aspect="auto", cmap="YlOrRd", origin="lower")

    ax.set_yticks(range(len(modifications)))
    ax.set_yticklabels(modifications)
    ax.set_xlabel("Genomic Position" if positions is None else "Position")
    ax.set_title("Histone Modifications")

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Signal Intensity")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Histone modification heatmap saved to {output_path}")

    return ax


def plot_differential_methylation(
    dm_results: List[Dict[str, Any]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot differential methylation results (volcano plot).

    Args:
        dm_results: List of DM result dictionaries with 'delta', 'pvalue', 'position' keys
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(dm_results, list, "dm_results")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Extract data
    deltas = np.array([result["delta"] for result in dm_results])
    pvalues = np.array([result["pvalue"] for result in dm_results])

    # Convert p-values to -log10
    neg_log_p = -np.log10(pvalues + 1e-10)  # Add small value to avoid log(0)

    # Color by significance
    significance_threshold = kwargs.get("sig_threshold", 0.05)
    significant = pvalues < significance_threshold
    colors = ["red" if sig else "gray" for sig in significant]

    ax.scatter(deltas, neg_log_p, c=colors, alpha=0.7, s=20)

    # Add significance line
    ax.axhline(
        -np.log10(significance_threshold), color="red", linestyle="--", alpha=0.7, label=f"p = {significance_threshold}"
    )

    ax.set_xlabel("Δ Methylation")
    ax.set_ylabel("-log₁₀(p-value)")
    ax.set_title("Differential Methylation Analysis")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Differential methylation plot saved to {output_path}")

    return ax


def plot_chromatin_states(
    states: np.ndarray,
    positions: np.ndarray,
    state_labels: Dict[int, str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 4),
    **kwargs,
) -> Axes:
    """Plot chromatin state segmentation.

    Args:
        states: Array of chromatin state assignments
        positions: Genomic positions
        state_labels: Optional mapping of state IDs to names
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(states, np.ndarray, "states")
    validation.validate_type(positions, np.ndarray, "positions")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Create color map for states
    unique_states = np.unique(states)
    colors = plt.cm.Set3(np.linspace(0, 1, len(unique_states)))
    state_color_map = dict(zip(unique_states, colors))

    # Plot state segments
    current_state = states[0]
    start_pos = positions[0]

    for i in range(1, len(states)):
        if states[i] != current_state or i == len(states) - 1:
            end_pos = positions[i]
            color = state_color_map[current_state]
            ax.fill_betweenx([0, 1], start_pos, end_pos, color=color, alpha=0.8)

            # Add state label
            if state_labels:
                label = state_labels.get(current_state, f"State {current_state}")
                ax.text((start_pos + end_pos) / 2, 0.5, label, ha="center", va="center", fontsize=8, rotation=90)

            current_state = states[i]
            start_pos = positions[i]

    ax.set_xlim(positions[0], positions[-1])
    ax.set_ylim(0, 1)
    ax.set_xlabel("Genomic Position")
    ax.set_title("Chromatin State Segmentation")
    ax.set_yticks([])

    # Add legend
    legend_elements = [
        plt.Rectangle(
            (0, 0),
            1,
            1,
            facecolor=color,
            label=state_labels.get(state, f"State {state}") if state_labels else f"State {state}",
        )
        for state, color in state_color_map.items()
    ]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Chromatin states plot saved to {output_path}")

    return ax


def plot_epigenetic_correlation_heatmap(
    epigenetic_data: Dict[str, np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Axes:
    """Plot correlation heatmap between epigenetic marks.

    Args:
        epigenetic_data: Dictionary mapping mark names to data arrays
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(epigenetic_data, dict, "epigenetic_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Calculate correlation matrix
    marks = list(epigenetic_data.keys())
    data_matrix = np.array([epigenetic_data[mark] for mark in marks])

    # Ensure same length
    min_len = min(len(arr) for arr in data_matrix)
    data_matrix = np.array([arr[:min_len] for arr in data_matrix])

    corr_matrix = np.corrcoef(data_matrix)

    # Plot heatmap
    if HAS_SEABORN:
        sns.heatmap(
            corr_matrix, annot=True, fmt=".2f", cmap="coolwarm", xticklabels=marks, yticklabels=marks, ax=ax, **kwargs
        )
    else:
        im = ax.imshow(corr_matrix, cmap="coolwarm", aspect="equal")
        ax.set_xticks(range(len(marks)))
        ax.set_yticks(range(len(marks)))
        ax.set_xticklabels(marks, rotation=45, ha="right")
        ax.set_yticklabels(marks)

        # Add correlation values
        for i in range(len(marks)):
            for j in range(len(marks)):
                ax.text(
                    j,
                    i,
                    f"{corr_matrix[i, j]:.2f}",
                    ha="center",
                    va="center",
                    color="white" if abs(corr_matrix[i, j]) > 0.7 else "black",
                )

        plt.colorbar(im, ax=ax)

    ax.set_title("Epigenetic Mark Correlations")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Epigenetic correlation heatmap saved to {output_path}")

    return ax


def plot_genome_browser_tracks(
    tracks: Dict[str, Dict[str, Any]],
    region_start: int,
    region_end: int,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot genome browser-style tracks for epigenetic data.

    Args:
        tracks: Dictionary of track data with 'data', 'color', 'label' keys
        region_start: Start position of genomic region
        region_end: End position of genomic region
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(tracks, dict, "tracks")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    track_names = list(tracks.keys())
    n_tracks = len(track_names)

    # Plot each track
    for i, track_name in enumerate(track_names):
        track_data = tracks[track_name]
        data = track_data.get("data", [])
        color = track_data.get("color", "blue")
        label = track_data.get("label", track_name)

        # Create subplot for this track
        ax_track = ax.inset_axes([0.1, 0.1 + i * 0.8 / n_tracks, 0.8, 0.8 / n_tracks])

        if isinstance(data, (list, np.ndarray)) and len(data) > 0:
            positions = np.linspace(region_start, region_end, len(data))
            ax_track.plot(positions, data, color=color, linewidth=1)
            ax_track.fill_between(positions, data, alpha=0.3, color=color)

        ax_track.set_xlim(region_start, region_end)
        ax_track.set_ylabel(label, fontsize=8)
        ax_track.set_xticks([])
        ax_track.tick_params(axis="y", labelsize=8)

        # Remove top and right spines
        ax_track.spines["top"].set_visible(False)
        ax_track.spines["right"].set_visible(False)

    # Set main axis properties
    ax.set_xlim(region_start, region_end)
    ax.set_xlabel("Genomic Position")
    ax.set_title("Genome Browser Tracks")
    ax.set_yticks([])
    ax.set_xticks([])

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Genome browser tracks saved to {output_path}")

    return ax


def plot_dna_methylation_clusters(
    methylation_matrix: np.ndarray,
    cluster_labels: np.ndarray | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot clustered DNA methylation data.

    Args:
        methylation_matrix: Samples x CpG sites methylation matrix
        cluster_labels: Optional cluster assignments for samples
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(methylation_matrix, np.ndarray, "methylation_matrix")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Sort by cluster if provided
    if cluster_labels is not None:
        sort_idx = np.argsort(cluster_labels)
        plot_data = methylation_matrix[sort_idx]
        cluster_boundaries = np.where(np.diff(cluster_labels[sort_idx]) != 0)[0] + 0.5
    else:
        plot_data = methylation_matrix
        cluster_boundaries = []

    # Plot heatmap
    im = ax.imshow(plot_data, aspect="auto", cmap="RdYlBu_r", vmin=0, vmax=1, origin="lower")

    # Add cluster boundaries
    for boundary in cluster_boundaries:
        ax.axhline(boundary, color="black", linewidth=2)

    ax.set_xlabel("CpG Sites")
    ax.set_ylabel("Samples")
    ax.set_title("DNA Methylation Clustering")

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Methylation Level")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Methylation clustering plot saved to {output_path}")

    return ax


def create_interactive_epigenome_browser(
    epigenetic_tracks: Dict[str, Any], *, output_path: str | Path | None = None, **kwargs
) -> Any:
    """Create an interactive epigenome browser using Plotly.

    Args:
        epigenetic_tracks: Dictionary containing track data
        output_path: Optional path to save the HTML file
        **kwargs: Additional arguments for Plotly customization

    Returns:
        Plotly Figure object
    """
    if not HAS_PLOTLY:
        raise ImportError("Plotly required for interactive epigenome browser")

    validation.validate_type(epigenetic_tracks, dict, "epigenetic_tracks")

    # Create subplot figure
    fig = go.Figure()

    track_names = list(epigenetic_tracks.keys())
    n_tracks = len(track_names)

    # Add each track as a subplot
    for i, track_name in enumerate(track_names):
        track_data = epigenetic_tracks[track_name]
        positions = track_data.get("positions", [])
        values = track_data.get("values", [])
        color = track_data.get("color", "blue")

        # Add trace for this track
        fig.add_trace(
            go.Scatter(
                x=positions,
                y=[i] * len(positions),  # Fixed y position for track
                mode="lines",
                line=dict(color=color, width=2),
                fill="tozeroy",
                name=track_name,
                hovertemplate=f"{track_name}: %{{x}}<br>Value: %{{customdata}}",
                customdata=values,
            )
        )

    # Update layout
    fig.update_layout(
        title="Interactive Epigenome Browser",
        xaxis_title="Genomic Position",
        yaxis=dict(tickmode="array", tickvals=list(range(n_tracks)), ticktext=track_names),
        **kwargs,
    )

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        html_path = Path(output_path).with_suffix(".html")
        fig.write_html(str(html_path))
        logger.info(f"Interactive epigenome browser saved to {html_path}")

    return fig
