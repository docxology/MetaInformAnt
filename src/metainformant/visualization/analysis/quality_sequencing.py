"""Sequencing quality control visualization functions.

Functions for visualizing read-level and sequence-level QC metrics
including quality scores, GC content, read lengths, adapter content,
duplication, overrepresented sequences, and k-mer profiles.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from metainformant.core import logging, paths, validation

logger = logging.get_logger(__name__)

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    sns = None


def plot_quality_metrics(
    qc_data: Dict[str, Any], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a comprehensive quality metrics visualization.

    Args:
        qc_data: Dictionary containing various QC metrics.
        ax: Optional matplotlib axes to plot on.
        output_path: Optional path to save the figure.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(qc_data, dict, "qc_data")
    if not qc_data:
        raise ValueError("QC data dictionary cannot be empty")

    fig, axes = plt.subplots(2, 2, figsize=kwargs.pop("figsize", (12, 10)))
    axes = axes.flatten()

    if "per_base_quality" in qc_data:
        qual_data = qc_data["per_base_quality"]
        if "positions" in qual_data and "mean_qualities" in qual_data:
            axes[0].plot(qual_data["positions"], qual_data["mean_qualities"])
            axes[0].set_xlabel("Position in Read")
            axes[0].set_ylabel("Mean Quality Score")
            axes[0].set_title("Per-Base Quality")
            axes[0].grid(True, alpha=0.3)

    if "gc_content_distribution" in qc_data:
        gc_data = qc_data["gc_content_distribution"]
        if "bins" in gc_data and "counts" in gc_data:
            bin_centers = [(gc_data["bins"][i] + gc_data["bins"][i + 1]) / 2 for i in range(len(gc_data["bins"]) - 1)]
            axes[1].bar(bin_centers, gc_data["counts"], width=2, alpha=0.7)
            axes[1].set_xlabel("GC Content (%)")
            axes[1].set_ylabel("Frequency")
            axes[1].set_title("GC Content Distribution")
            axes[1].grid(True, alpha=0.3)

    if "sequence_length_distribution" in qc_data:
        length_data = qc_data["sequence_length_distribution"]
        if "lengths" in length_data and "counts" in length_data:
            axes[2].bar(length_data["lengths"], length_data["counts"], alpha=0.7)
            axes[2].set_xlabel("Read Length")
            axes[2].set_ylabel("Count")
            axes[2].set_title("Read Length Distribution")
            axes[2].grid(True, alpha=0.3)

    if "basic_statistics" in qc_data:
        stats = qc_data["basic_statistics"]
        stat_names = ["num_reads", "total_bases", "min_length", "max_length", "mean_length"]
        stat_values = [stats.get(name, 0) for name in stat_names]
        bars = axes[3].bar(stat_names, stat_values, alpha=0.7)
        axes[3].set_ylabel("Value")
        axes[3].set_title("Basic Statistics")
        axes[3].tick_params(axis="x", rotation=45)
        for bar, value in zip(bars, stat_values):
            axes[3].text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                         f"{value:.0f}", ha="center", va="bottom", fontsize=8)

    plt.tight_layout()

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Quality metrics plot saved to {output_path}")

    return axes[0]


def plot_adapter_content(
    adapter_data: Dict[str, List[float]], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create an adapter content visualization.

    Args:
        adapter_data: Dictionary with adapter names as keys and content percentages as values.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(adapter_data, dict, "adapter_data")
    if not adapter_data:
        raise ValueError("Adapter data dictionary cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 6)))

    adapters = list(adapter_data.keys())
    percentages = [max(values) if values else 0 for values in adapter_data.values()]

    sorted_indices = np.argsort(percentages)[::-1]
    adapters_sorted = [adapters[i] for i in sorted_indices]
    percentages_sorted = [percentages[i] for i in sorted_indices]

    bars = ax.bar(range(len(adapters_sorted)), percentages_sorted, alpha=0.7, **kwargs)
    ax.set_xlabel("Adapter Type")
    ax.set_ylabel("Maximum Content (%)")
    ax.set_title("Adapter Content Analysis")
    ax.set_xticks(range(len(adapters_sorted)))
    ax.set_xticklabels(adapters_sorted, rotation=45, ha="right")
    ax.grid(True, alpha=0.3)

    for bar, pct in zip(bars, percentages_sorted):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{pct:.1f}%", ha="center", va="bottom", fontsize=8)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Adapter content plot saved to {output_path}")

    return ax


def plot_gc_distribution(
    gc_data: List[float], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a GC content distribution plot.

    Args:
        gc_data: List of GC content percentages.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(gc_data, (list, np.ndarray), "gc_data")
    if len(gc_data) == 0:
        raise ValueError("GC data list cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    n, bins, patches = ax.hist(gc_data, bins=kwargs.pop("bins", 20), alpha=0.7, edgecolor="black", **kwargs)
    ax.set_xlabel("GC Content (%)")
    ax.set_ylabel("Frequency")
    ax.set_title("GC Content Distribution")
    ax.grid(True, alpha=0.3)

    mean_gc = np.mean(gc_data)
    median_gc = np.median(gc_data)
    std_gc = np.std(gc_data)
    stats_text = f"Mean: {mean_gc:.1f}%\nMedian: {median_gc:.1f}%\nStd: {std_gc:.1f}%"
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8))

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"GC distribution plot saved to {output_path}")

    return ax


def plot_length_distribution(
    length_data: List[int], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a read length distribution plot.

    Args:
        length_data: List of read lengths.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(length_data, (list, np.ndarray), "length_data")
    if len(length_data) == 0:
        raise ValueError("Length data list cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    n, bins, patches = ax.hist(length_data, bins=kwargs.pop("bins", "auto"), alpha=0.7, edgecolor="black", **kwargs)
    ax.set_xlabel("Read Length (bp)")
    ax.set_ylabel("Frequency")
    ax.set_title("Read Length Distribution")
    ax.grid(True, alpha=0.3)

    mean_len = np.mean(length_data)
    median_len = np.median(length_data)
    min_len = np.min(length_data)
    max_len = np.max(length_data)
    stats_text = f"Mean: {mean_len:.0f} bp\nMedian: {median_len:.0f} bp\nRange: {min_len}-{max_len} bp"
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8))

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Length distribution plot saved to {output_path}")

    return ax


def plot_per_base_quality_boxplot(
    per_base_qualities: Dict[str, List[float]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes:
    """Create a per-base quality boxplot from FastQC-like data.

    Args:
        per_base_qualities: Dictionary mapping position ranges to quality score lists.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.
        figsize: Figure size.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(per_base_qualities, dict, "per_base_qualities")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    positions = []
    qualities = []
    for pos_range, qual_list in per_base_qualities.items():
        if isinstance(pos_range, str) and "-" in pos_range:
            start, end = map(int, pos_range.split("-"))
            mid_pos = (start + end) / 2
        else:
            mid_pos = int(pos_range)
        positions.append(mid_pos)
        qualities.append(qual_list)

    bp = ax.boxplot(qualities, positions=positions, widths=2, patch_artist=True, **kwargs)
    colors = plt.cm.RdYlGn(np.linspace(0, 1, len(qualities)))
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xlabel("Position in Read")
    ax.set_ylabel("Quality Score")
    ax.set_title("Per-Base Quality Scores")
    ax.grid(True, alpha=0.3)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Per-base quality boxplot saved to {output_path}")

    return ax


def plot_sequence_duplication_levels(
    duplication_levels: Dict[str, float],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 6),
    **kwargs,
) -> Axes:
    """Plot sequence duplication levels from FastQC analysis.

    Args:
        duplication_levels: Dictionary mapping duplication levels to percentages.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.
        figsize: Figure size.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(duplication_levels, dict, "duplication_levels")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    levels = list(duplication_levels.keys())
    percentages = list(duplication_levels.values())

    bars = ax.bar(range(len(levels)), percentages, color="skyblue", alpha=0.8, **kwargs)
    ax.set_xlabel("Duplication Level")
    ax.set_ylabel("Percentage (%)")
    ax.set_title("Sequence Duplication Levels")
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels(levels, rotation=45, ha="right")
    ax.grid(True, alpha=0.3, axis="y")

    for bar, pct in zip(bars, percentages):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f"{pct:.1f}%", ha="center", va="bottom", fontsize=8)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Sequence duplication levels plot saved to {output_path}")

    return ax


def plot_overrepresented_sequences(
    overrepresented_seqs: List[Dict[str, Any]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot overrepresented sequences analysis.

    Args:
        overrepresented_seqs: List of dictionaries with sequence info.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.
        figsize: Figure size.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(overrepresented_seqs, list, "overrepresented_seqs")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if not overrepresented_seqs:
        ax.text(0.5, 0.5, "No overrepresented sequences found", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Overrepresented Sequences")
        return ax

    sequences = [seq.get("sequence", f"Seq {i+1}")[:20] + "..." for i, seq in enumerate(overrepresented_seqs)]
    counts = [seq.get("count", 0) for seq in overrepresented_seqs]
    percentages = [seq.get("percentage", 0) for seq in overrepresented_seqs]

    y_pos = np.arange(len(sequences))
    bars = ax.barh(y_pos, percentages, color="coral", alpha=0.8, **kwargs)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sequences)
    ax.set_xlabel("Percentage of Total Reads")
    ax.set_title("Overrepresented Sequences")
    ax.grid(True, alpha=0.3, axis="x")

    for i, (bar, pct, count) in enumerate(zip(bars, percentages, counts)):
        ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height() / 2,
                f"{pct:.2f}%\n({count})", ha="left", va="center", fontsize=8)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Overrepresented sequences plot saved to {output_path}")

    return ax


def plot_kmer_profiles(
    kmer_counts: Dict[str, int],
    top_n: int = 20,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes:
    """Plot k-mer frequency profiles.

    Args:
        kmer_counts: Dictionary mapping k-mers to their counts.
        top_n: Number of top k-mers to display.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.
        figsize: Figure size.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(kmer_counts, dict, "kmer_counts")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sorted_kmers = sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)[:top_n]
    kmers, counts = zip(*sorted_kmers)

    bars = ax.bar(range(len(kmers)), counts, color="purple", alpha=0.8, **kwargs)
    ax.set_xlabel("K-mer")
    ax.set_ylabel("Count")
    ax.set_title(f"Top {top_n} K-mer Frequencies")
    ax.set_xticks(range(len(kmers)))
    ax.set_xticklabels(kmers, rotation=45, ha="right")
    ax.grid(True, alpha=0.3, axis="y")

    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(counts) * 0.01,
                f"{count}", ha="center", va="bottom", fontsize=8)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"K-mer profiles plot saved to {output_path}")

    return ax
