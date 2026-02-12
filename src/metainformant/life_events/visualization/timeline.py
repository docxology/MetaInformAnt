"""Timeline visualization functions for life events analysis.

This module provides timeline-related plotting utilities for life event sequences,
including event timelines, domain timelines, and temporal pattern visualizations.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Optional dependencies
try:
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    import numpy as np

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    logger.warning("matplotlib not available, life events timeline visualization disabled")


def plot_event_timeline(
    sequence: Any,
    output_path: Optional[str] = None,
    figsize: tuple[int, int] = (12, 6),
    show_domains: bool = True,
    title: str = "Life Event Timeline",
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
    domains = [event.domain for event in events] if hasattr(events[0], "domain") else ["default"] * len(events)
    unique_domains = list(set(domains))

    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_domains)))
    domain_colors = dict(zip(unique_domains, colors))

    for i, event in enumerate(events):
        color = domain_colors.get(event.domain if hasattr(event, "domain") else "default", "blue")
        ax.scatter(event.timestamp, i, c=[color], s=100, alpha=0.7, edgecolors="black", linewidth=1)

        # Add event type label
        ax.annotate(
            event.event_type,
            (event.timestamp, i),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
            ha="left",
            va="bottom",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8),
        )

    # Format x-axis as dates
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.YearLocator())
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")

    ax.set_xlabel("Date")
    ax.set_ylabel("Event Order")
    ax.set_title(f"{title} - {sequence.person_id}")
    ax.grid(True, alpha=0.3)

    # Add domain legend if requested
    if show_domains and len(unique_domains) > 1:
        legend_elements = [
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, markersize=10, label=domain)
            for domain, color in domain_colors.items()
        ]
        ax.legend(handles=legend_elements, title="Domains", loc="upper right")

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved timeline plot to {output_path}")

    return fig


def plot_domain_distribution(
    sequences: List[Any],
    output_path: Optional[str] = None,
    figsize: tuple[int, int] = (12, 6),
    title: str = "Event Domain Distribution",
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
            if hasattr(event, "domain"):
                domain = event.domain
                domain_counts[domain] = domain_counts.get(domain, 0) + 1

    if not domain_counts:
        logger.warning("No domain information found in sequences")
        return None

    domains, counts = zip(*domain_counts.items())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Bar plot
    bars = ax1.bar(range(len(domains)), counts, color="lightcoral", alpha=0.7)
    ax1.set_xticks(range(len(domains)))
    ax1.set_xticklabels(domains, rotation=45, ha="right")
    ax1.set_ylabel("Count")
    ax1.set_title("Domain Counts")
    ax1.grid(True, alpha=0.3)

    # Pie chart
    ax2.pie(counts, labels=domains, autopct="%1.1f%%", startangle=90)
    ax2.set_title("Domain Proportions")
    ax2.axis("equal")

    fig.suptitle(title, fontsize=14)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved domain distribution plot to {output_path}")

    return fig


def plot_domain_timeline(
    sequences: List[Any],
    output_path: Optional[Union[str, Any]] = None,
    max_sequences: int = 10,
    figsize: Tuple[int, int] = (12, 8),
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

    fig, ax = plt.subplots(figsize=figsize)

    # Define colors for different domains
    domain_colors = {
        "education": "#1f77b4",  # blue
        "career": "#ff7f0e",  # orange
        "family": "#2ca02c",  # green
        "health": "#d62728",  # red
        "finance": "#9467bd",  # purple
        "housing": "#8c564b",  # brown
        "other": "#7f7f7f",  # gray
    }

    # Plot each sequence
    for i, seq in enumerate(sequences[:max_sequences]):
        y_pos = len(sequences[:max_sequences]) - i - 1  # Bottom to top

        # Group events by domain
        domain_periods = {}
        for event in seq.events:
            domain = getattr(event, "domain", "other")
            if domain not in domain_periods:
                domain_periods[domain] = []
            domain_periods[domain].append(event.timestamp)

        # Plot periods for each domain
        for domain, timestamps in domain_periods.items():
            if timestamps:
                timestamps.sort()
                color = domain_colors.get(domain, domain_colors["other"])

                for ts in timestamps:
                    ax.scatter(ts, y_pos, color=color, s=50, alpha=0.7, marker="o")

                ax.text(
                    timestamps[0], y_pos + 0.3, domain.capitalize(), fontsize=8, ha="left", va="bottom", color=color
                )

    # Customize plot
    ax.set_yticks(range(len(sequences[:max_sequences])))
    ax.set_yticklabels([f"Person {seq.person_id}" for seq in sequences[:max_sequences]])
    ax.set_xlabel("Time")
    ax.set_title("Life Domain Timeline")
    ax.grid(True, alpha=0.3)

    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, markersize=8, label=domain.capitalize())
        for domain, color in domain_colors.items()
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=8)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved domain timeline plot to {output_path}")

    return fig


def plot_temporal_density(
    sequences: List[Any],
    output_path: Optional[Union[str, Any]] = None,
    bins: int = 20,
    figsize: Tuple[int, int] = (12, 6),
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
    timestamps_numeric = [(ts - min(all_timestamps)).total_seconds() / (24 * 3600) for ts in all_timestamps]

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Plot histogram
    counts, bin_edges, _ = axes[0].hist(timestamps_numeric, bins=bins, alpha=0.7, edgecolor="black")
    axes[0].set_xlabel("Time (days)")
    axes[0].set_ylabel("Event Count")
    axes[0].set_title("Event Temporal Distribution")
    axes[0].grid(True, alpha=0.3)

    # Plot density estimate
    try:
        from scipy import stats

        kde = stats.gaussian_kde(timestamps_numeric)
        x_range = np.linspace(min(timestamps_numeric), max(timestamps_numeric), 200)
        density = kde(x_range)

        axes[1].plot(x_range, density, "b-", linewidth=2)
        axes[1].fill_between(x_range, density, alpha=0.3)
        axes[1].set_xlabel("Time (days)")
        axes[1].set_ylabel("Density")
        axes[1].set_title("Event Temporal Density (KDE)")
        axes[1].grid(True, alpha=0.3)

    except ImportError:
        axes[1].hist(timestamps_numeric, bins=bins, alpha=0.7, edgecolor="black")
        axes[1].set_xlabel("Time (days)")
        axes[1].set_ylabel("Event Count")
        axes[1].set_title("Event Temporal Distribution (Histogram)")
        axes[1].grid(True, alpha=0.3)

    # Add statistics
    mean_time = np.mean(timestamps_numeric)
    median_time = np.median(timestamps_numeric)
    std_time = np.std(timestamps_numeric)

    stats_text = f"Events: {len(all_timestamps)}\nMean: {mean_time:.1f} days\nMedian: {median_time:.1f} days\nStd: {std_time:.1f} days"
    axes[0].text(
        0.02,
        0.98,
        stats_text,
        transform=axes[0].transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved temporal density plot to {output_path}")

    return fig


def plot_temporal_patterns(
    sequences: List[Any],
    importance_scores: Optional[Dict[int, float]] = None,
    output_path: Optional[Union[str, Any]] = None,
    figsize: Tuple[int, int] = (14, 8),
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
            if importance_scores and i in importance_scores:
                alpha = 0.3 + 0.7 * importance_scores[i]
                seq_color = (*seq_color[:3], alpha)

            for event in seq.events:
                rel_time = (event.timestamp - min_time).total_seconds() / (24 * 3600)
                ax1.scatter(rel_time, i, color=seq_color, s=50, alpha=0.8)
                ax1.text(
                    rel_time,
                    i + 0.1,
                    event.event_type.split(":")[-1],
                    fontsize=6,
                    ha="center",
                    va="bottom",
                    rotation=45,
                )

    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("Sequence Index")
    ax1.set_title("Temporal Patterns of Life Events")
    ax1.grid(True, alpha=0.3)

    # Plot 2: Event type frequency over time
    ax2 = axes[1]

    if all_timestamps:
        time_bins = np.linspace(0, (max_time - min_time).total_seconds() / (24 * 3600), 10)

        event_types = set()
        for seq in sequences:
            for event in seq.events:
                event_types.add(event.event_type)

        event_types = sorted(list(event_types))

        freq_matrix = np.zeros((len(event_types), len(time_bins) - 1))

        for seq in sequences:
            for event in seq.events:
                rel_time = (event.timestamp - min_time).total_seconds() / (24 * 3600)
                time_bin = np.digitize(rel_time, time_bins) - 1
                if 0 <= time_bin < len(time_bins) - 1:
                    event_idx = event_types.index(event.event_type)
                    freq_matrix[event_idx, time_bin] += 1

        im = ax2.imshow(freq_matrix, aspect="auto", cmap="Blues")
        ax2.set_xticks(range(len(time_bins) - 1))
        ax2.set_yticks(range(len(event_types)))
        ax2.set_xticklabels([f"{time_bins[i]:.1f}" for i in range(len(time_bins) - 1)])
        ax2.set_yticklabels([et.split(":")[-1] for et in event_types])
        ax2.set_xlabel("Time Bin (days)")
        ax2.set_ylabel("Event Type")
        ax2.set_title("Event Frequency Over Time")

        plt.colorbar(im, ax=ax2, shrink=0.8)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved temporal patterns plot to {output_path}")

    return fig
