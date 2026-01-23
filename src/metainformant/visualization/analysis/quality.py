"""Quality control data visualization functions.

This module provides comprehensive plotting functions for quality control across all
biological data types, including sequencing, genomic variants, protein structures,
single-cell data, and multi-omics quality assessment.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

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
        qc_data: Dictionary containing various QC metrics
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If qc_data is empty or malformed
    """
    validation.validate_type(qc_data, dict, "qc_data")

    if not qc_data:
        raise ValueError("QC data dictionary cannot be empty")

    if ax is None:
        fig, axes = plt.subplots(2, 2, figsize=kwargs.pop("figsize", (12, 10)))
        axes = axes.flatten()
    else:
        # If single ax provided, this is complex - just use it for a simple plot
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()

    # Plot 1: Quality scores distribution
    if "per_base_quality" in qc_data:
        qual_data = qc_data["per_base_quality"]
        if "positions" in qual_data and "mean_qualities" in qual_data:
            axes[0].plot(qual_data["positions"], qual_data["mean_qualities"])
            axes[0].set_xlabel("Position in Read")
            axes[0].set_ylabel("Mean Quality Score")
            axes[0].set_title("Per-Base Quality")
            axes[0].grid(True, alpha=0.3)

    # Plot 2: GC content distribution
    if "gc_content_distribution" in qc_data:
        gc_data = qc_data["gc_content_distribution"]
        if "bins" in gc_data and "counts" in gc_data:
            # Convert bins to centers for plotting
            bin_centers = [(gc_data["bins"][i] + gc_data["bins"][i + 1]) / 2 for i in range(len(gc_data["bins"]) - 1)]
            axes[1].bar(bin_centers, gc_data["counts"], width=2, alpha=0.7)
            axes[1].set_xlabel("GC Content (%)")
            axes[1].set_ylabel("Frequency")
            axes[1].set_title("GC Content Distribution")
            axes[1].grid(True, alpha=0.3)

    # Plot 3: Read length distribution
    if "sequence_length_distribution" in qc_data:
        length_data = qc_data["sequence_length_distribution"]
        if "lengths" in length_data and "counts" in length_data:
            axes[2].bar(length_data["lengths"], length_data["counts"], alpha=0.7)
            axes[2].set_xlabel("Read Length")
            axes[2].set_ylabel("Count")
            axes[2].set_title("Read Length Distribution")
            axes[2].grid(True, alpha=0.3)

    # Plot 4: Basic statistics summary
    if "basic_statistics" in qc_data:
        stats = qc_data["basic_statistics"]
        stat_names = ["num_reads", "total_bases", "min_length", "max_length", "mean_length"]
        stat_values = [stats.get(name, 0) for name in stat_names]

        bars = axes[3].bar(stat_names, stat_values, alpha=0.7)
        axes[3].set_ylabel("Value")
        axes[3].set_title("Basic Statistics")
        axes[3].tick_params(axis="x", rotation=45)

        # Add value labels on bars
        for bar, value in zip(bars, stat_values):
            axes[3].text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                f"{value:.0f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    plt.tight_layout()

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Quality metrics plot saved to {output_path}")

    return axes[0]  # Return first axis for consistency


def plot_adapter_content(
    adapter_data: Dict[str, List[float]], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create an adapter content visualization.

    Args:
        adapter_data: Dictionary with adapter names as keys and content percentages as values
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If adapter_data is empty
    """
    validation.validate_type(adapter_data, dict, "adapter_data")

    if not adapter_data:
        raise ValueError("Adapter data dictionary cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 6)))

    # Prepare data for plotting
    adapters = list(adapter_data.keys())
    # Assume each adapter has a list of percentages (could be by position)
    # For simplicity, take the maximum percentage for each adapter
    percentages = [max(values) if values else 0 for values in adapter_data.values()]

    # Sort by percentage
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

    # Add percentage labels
    for bar, pct in zip(bars, percentages_sorted):
        ax.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{pct:.1f}%", ha="center", va="bottom", fontsize=8
        )

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
        gc_data: List of GC content percentages
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib hist().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If gc_data is empty
    """
    validation.validate_type(gc_data, (list, np.ndarray), "gc_data")
    if len(gc_data) == 0:
        raise ValueError("GC data list cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    # Plot histogram
    n, bins, patches = ax.hist(gc_data, bins=kwargs.get("bins", 20), alpha=0.7, edgecolor="black", **kwargs)

    ax.set_xlabel("GC Content (%)")
    ax.set_ylabel("Frequency")
    ax.set_title("GC Content Distribution")
    ax.grid(True, alpha=0.3)

    # Add statistics as text
    mean_gc = np.mean(gc_data)
    median_gc = np.median(gc_data)
    std_gc = np.std(gc_data)

    stats_text = f"Mean: {mean_gc:.1f}%\nMedian: {median_gc:.1f}%\nStd: {std_gc:.1f}%"
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

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
        length_data: List of read lengths
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib hist().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If length_data is empty
    """
    validation.validate_type(length_data, (list, np.ndarray), "length_data")

    if len(length_data) == 0:
        raise ValueError("Length data list cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    # Plot histogram
    n, bins, patches = ax.hist(length_data, bins=kwargs.get("bins", "auto"), alpha=0.7, edgecolor="black", **kwargs)

    ax.set_xlabel("Read Length (bp)")
    ax.set_ylabel("Frequency")
    ax.set_title("Read Length Distribution")
    ax.grid(True, alpha=0.3)

    # Add statistics as text
    mean_len = np.mean(length_data)
    median_len = np.median(length_data)
    min_len = np.min(length_data)
    max_len = np.max(length_data)

    stats_text = f"Mean: {mean_len:.0f} bp\nMedian: {median_len:.0f} bp\nRange: {min_len}-{max_len} bp"
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

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
        per_base_qualities: Dictionary mapping position ranges to quality score lists
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

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
            # Handle range strings like "1-10"
            start, end = map(int, pos_range.split("-"))
            mid_pos = (start + end) / 2
        else:
            mid_pos = int(pos_range)

        positions.append(mid_pos)
        qualities.append(qual_list)

    # Create boxplot
    bp = ax.boxplot(qualities, positions=positions, widths=2, patch_artist=True, **kwargs)

    # Color boxes
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
        duplication_levels: Dictionary mapping duplication levels to percentages
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

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

    # Add percentage labels
    for bar, pct in zip(bars, percentages):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.5,
            f"{pct:.1f}%",
            ha="center",
            va="bottom",
            fontsize=8,
        )

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
        overrepresented_seqs: List of dictionaries with sequence info
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

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

    # Extract data
    sequences = [seq.get("sequence", f"Seq {i+1}")[:20] + "..." for i, seq in enumerate(overrepresented_seqs)]
    counts = [seq.get("count", 0) for seq in overrepresented_seqs]
    percentages = [seq.get("percentage", 0) for seq in overrepresented_seqs]

    # Plot as horizontal bar chart
    y_pos = np.arange(len(sequences))
    bars = ax.barh(y_pos, percentages, color="coral", alpha=0.8, **kwargs)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(sequences)
    ax.set_xlabel("Percentage of Total Reads")
    ax.set_title("Overrepresented Sequences")
    ax.grid(True, alpha=0.3, axis="x")

    # Add percentage labels
    for i, (bar, pct, count) in enumerate(zip(bars, percentages, counts)):
        ax.text(
            bar.get_width() + 0.1,
            bar.get_y() + bar.get_height() / 2,
            f"{pct:.2f}%\n({count})",
            ha="left",
            va="center",
            fontsize=8,
        )

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
        kmer_counts: Dictionary mapping k-mers to their counts
        top_n: Number of top k-mers to display
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(kmer_counts, dict, "kmer_counts")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Sort by count and take top N
    sorted_kmers = sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)[:top_n]
    kmers, counts = zip(*sorted_kmers)

    bars = ax.bar(range(len(kmers)), counts, color="purple", alpha=0.8, **kwargs)
    ax.set_xlabel("K-mer")
    ax.set_ylabel("Count")
    ax.set_title(f"Top {top_n} K-mer Frequencies")
    ax.set_xticks(range(len(kmers)))
    ax.set_xticklabels(kmers, rotation=45, ha="right")
    ax.grid(True, alpha=0.3, axis="y")

    # Add count labels
    for bar, count in zip(bars, counts):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(counts) * 0.01,
            f"{count}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"K-mer profiles plot saved to {output_path}")

    return ax


def plot_vcf_quality_metrics(
    vcf_qc_data: Dict[str, Any],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (14, 10),
    **kwargs,
) -> Axes:
    """Plot comprehensive VCF quality control metrics.

    Args:
        vcf_qc_data: Dictionary containing VCF QC metrics
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(vcf_qc_data, dict, "vcf_qc_data")

    if ax is None:
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()
    else:
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()

    # Plot 1: Quality score distribution
    if "qual_distribution" in vcf_qc_data:
        qual_data = vcf_qc_data["qual_distribution"]
        if "qualities" in qual_data and "counts" in qual_data:
            axes[0].bar(qual_data["qualities"], qual_data["counts"], alpha=0.7)
            axes[0].set_xlabel("Quality Score")
            axes[0].set_ylabel("Count")
            axes[0].set_title("Variant Quality Distribution")
            axes[0].grid(True, alpha=0.3)

    # Plot 2: Depth distribution
    if "depth_distribution" in vcf_qc_data:
        depth_data = vcf_qc_data["depth_distribution"]
        if "depths" in depth_data and "counts" in depth_data:
            axes[1].bar(depth_data["depths"], depth_data["counts"], alpha=0.7, color="green")
            axes[1].set_xlabel("Read Depth")
            axes[1].set_ylabel("Count")
            axes[1].set_title("Read Depth Distribution")
            axes[1].grid(True, alpha=0.3)

    # Plot 3: Allele frequency spectrum
    if "allele_frequencies" in vcf_qc_data:
        af_data = vcf_qc_data["allele_frequencies"]
        axes[2].hist(af_data, bins=20, alpha=0.7, color="orange")
        axes[2].set_xlabel("Allele Frequency")
        axes[2].set_ylabel("Count")
        axes[2].set_title("Allele Frequency Spectrum")
        axes[2].grid(True, alpha=0.3)

    # Plot 4: Variant type distribution
    if "variant_types" in vcf_qc_data:
        types_data = vcf_qc_data["variant_types"]
        types = list(types_data.keys())
        counts = list(types_data.values())
        axes[3].pie(counts, labels=types, autopct="%1.1f%%", startangle=90)
        axes[3].set_title("Variant Types")

    # Plot 5: Ti/Tv ratio by quality
    if "titv_by_qual" in vcf_qc_data:
        titv_data = vcf_qc_data["titv_by_qual"]
        if "qualities" in titv_data and "titv_ratios" in titv_data:
            axes[4].plot(titv_data["qualities"], titv_data["titv_ratios"], "o-", alpha=0.7)
            axes[4].axhline(y=2.1, color="red", linestyle="--", alpha=0.7, label="Expected")
            axes[4].set_xlabel("Quality Score")
            axes[4].set_ylabel("Ti/Tv Ratio")
            axes[4].set_title("Transition/Transversion Ratio")
            axes[4].legend()
            axes[4].grid(True, alpha=0.3)

    # Plot 6: Summary statistics
    if "summary_stats" in vcf_qc_data:
        stats = vcf_qc_data["summary_stats"]
        stat_names = list(stats.keys())
        stat_values = list(stats.values())

        bars = axes[5].bar(range(len(stat_names)), stat_values, alpha=0.7, color="purple")
        axes[5].set_xticks(range(len(stat_names)))
        axes[5].set_xticklabels(stat_names, rotation=45, ha="right")
        axes[5].set_ylabel("Value")
        axes[5].set_title("Summary Statistics")
        axes[5].grid(True, alpha=0.3, axis="y")

        # Add value labels
        for bar, value in zip(bars, stat_values):
            axes[5].text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(stat_values) * 0.01,
                f"{value:.0f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    plt.tight_layout()

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"VCF quality metrics plot saved to {output_path}")

    return axes[0]


def plot_singlecell_qc_metrics(
    qc_metrics: Dict[str, np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (14, 10),
    **kwargs,
) -> Axes:
    """Plot comprehensive single-cell QC metrics.

    Args:
        qc_metrics: Dictionary containing various QC metrics per cell
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(qc_metrics, dict, "qc_metrics")

    if ax is None:
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()
    else:
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()

    plot_idx = 0

    # Plot library size distribution
    if "n_counts" in qc_metrics:
        if HAS_SEABORN:
            sns.histplot(qc_metrics["n_counts"], ax=axes[plot_idx], kde=True, alpha=0.7)
        else:
            axes[plot_idx].hist(qc_metrics["n_counts"], bins=50, alpha=0.7)
        axes[plot_idx].set_xlabel("Library Size")
        axes[plot_idx].set_ylabel("Number of Cells")
        axes[plot_idx].set_title("Library Size Distribution")
        axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot number of genes detected
    if "n_genes" in qc_metrics:
        if HAS_SEABORN:
            sns.histplot(qc_metrics["n_genes"], ax=axes[plot_idx], kde=True, alpha=0.7, color="green")
        else:
            axes[plot_idx].hist(qc_metrics["n_genes"], bins=50, alpha=0.7, color="green")
        axes[plot_idx].set_xlabel("Number of Genes")
        axes[plot_idx].set_ylabel("Number of Cells")
        axes[plot_idx].set_title("Genes Detected per Cell")
        axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot mitochondrial content
    if "percent_mito" in qc_metrics:
        if HAS_SEABORN:
            sns.histplot(qc_metrics["percent_mito"], ax=axes[plot_idx], kde=True, alpha=0.7, color="red")
        else:
            axes[plot_idx].hist(qc_metrics["percent_mito"], bins=50, alpha=0.7, color="red")
        axes[plot_idx].set_xlabel("Mitochondrial Content (%)")
        axes[plot_idx].set_ylabel("Number of Cells")
        axes[plot_idx].set_title("Mitochondrial Content")
        axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot complexity (genes per UMI)
    if "n_counts" in qc_metrics and "n_genes" in qc_metrics:
        complexity = qc_metrics["n_genes"] / qc_metrics["n_counts"]
        if HAS_SEABORN:
            sns.histplot(complexity, ax=axes[plot_idx], kde=True, alpha=0.7, color="purple")
        else:
            axes[plot_idx].hist(complexity, bins=50, alpha=0.7, color="purple")
        axes[plot_idx].set_xlabel("Genes per UMI")
        axes[plot_idx].set_ylabel("Number of Cells")
        axes[plot_idx].set_title("Library Complexity")
        axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot doublet scores (if available)
    if "doublet_score" in qc_metrics:
        if HAS_SEABORN:
            sns.histplot(qc_metrics["doublet_score"], ax=axes[plot_idx], kde=True, alpha=0.7, color="orange")
        else:
            axes[plot_idx].hist(qc_metrics["doublet_score"], bins=50, alpha=0.7, color="orange")
        axes[plot_idx].set_xlabel("Doublet Score")
        axes[plot_idx].set_ylabel("Number of Cells")
        axes[plot_idx].set_title("Doublet Scores")
        axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot QC metric correlations
    if len(qc_metrics) >= 2:
        metrics_to_plot = ["n_counts", "n_genes", "percent_mito"]
        available_metrics = [m for m in metrics_to_plot if m in qc_metrics]

        if len(available_metrics) >= 2:
            # Create correlation plot
            corr_data = np.column_stack([qc_metrics[m] for m in available_metrics])
            corr_matrix = np.corrcoef(corr_data.T)

            if HAS_SEABORN:
                sns.heatmap(
                    corr_matrix,
                    annot=True,
                    fmt=".2f",
                    cmap="coolwarm",
                    xticklabels=available_metrics,
                    yticklabels=available_metrics,
                    ax=axes[plot_idx],
                    center=0,
                )
            else:
                im = axes[plot_idx].imshow(corr_matrix, cmap="coolwarm", aspect="equal", vmin=-1, vmax=1)
                axes[plot_idx].set_xticks(range(len(available_metrics)))
                axes[plot_idx].set_yticks(range(len(available_metrics)))
                axes[plot_idx].set_xticklabels(available_metrics, rotation=45, ha="right")
                axes[plot_idx].set_yticklabels(available_metrics)
                plt.colorbar(im, ax=axes[plot_idx])

            axes[plot_idx].set_title("QC Metric Correlations")

    plt.tight_layout()

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Single-cell QC metrics plot saved to {output_path}")

    return axes[0]


def plot_protein_structure_quality(
    structure_quality: Dict[str, Any],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot protein structure quality metrics.

    Args:
        structure_quality: Dictionary containing structure quality metrics
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(structure_quality, dict, "structure_quality")

    if ax is None:
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        axes = axes.flatten()
    else:
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        axes = axes.flatten()

    # Plot 1: B-factor distribution
    if "b_factors" in structure_quality:
        b_factors = structure_quality["b_factors"]
        axes[0].hist(b_factors, bins=50, alpha=0.7, color="blue")
        axes[0].set_xlabel("B-factor")
        axes[0].set_ylabel("Count")
        axes[0].set_title("B-factor Distribution")
        axes[0].grid(True, alpha=0.3)

    # Plot 2: Ramachandran plot quality
    if "ramachandran_stats" in structure_quality:
        rama_stats = structure_quality["ramachandran_stats"]
        categories = list(rama_stats.keys())
        values = list(rama_stats.values())
        axes[1].bar(categories, values, alpha=0.7, color="green")
        axes[1].set_ylabel("Percentage")
        axes[1].set_title("Ramachandran Plot")
        axes[1].tick_params(axis="x", rotation=45)
        axes[1].grid(True, alpha=0.3, axis="y")

    # Plot 3: Clash score
    if "clash_score" in structure_quality:
        clash_data = structure_quality["clash_score"]
        if isinstance(clash_data, dict):
            clash_types = list(clash_data.keys())
            clash_counts = list(clash_data.values())
            axes[2].bar(clash_types, clash_counts, alpha=0.7, color="red")
            axes[2].set_ylabel("Number of Clashes")
            axes[2].set_title("Steric Clashes")
            axes[2].tick_params(axis="x", rotation=45)
        else:
            axes[2].text(
                0.5, 0.5, f"Clash Score: {clash_data:.2f}", ha="center", va="center", transform=axes[2].transAxes
            )
            axes[2].set_title("Clash Score")
        axes[2].grid(True, alpha=0.3, axis="y")

    # Plot 4: Overall quality score
    if "overall_quality" in structure_quality:
        quality_data = structure_quality["overall_quality"]
        if isinstance(quality_data, dict):
            # Plot multiple quality metrics
            metrics = list(quality_data.keys())
            values = list(quality_data.values())
            bars = axes[3].bar(metrics, values, alpha=0.7, color="purple")
            axes[3].set_ylabel("Score")
            axes[3].set_title("Quality Metrics")
            axes[3].tick_params(axis="x", rotation=45)

            # Add value labels
            for bar, value in zip(bars, values):
                axes[3].text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.01,
                    f"{value:.2f}",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )
        else:
            # Single overall score
            axes[3].text(
                0.5, 0.5, f"Overall Quality: {quality_data:.2f}", ha="center", va="center", transform=axes[3].transAxes
            )
            axes[3].set_title("Overall Quality Score")
        axes[3].grid(True, alpha=0.3, axis="y")

    plt.tight_layout()

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Protein structure quality plot saved to {output_path}")

    return axes[0]


def plot_multiomics_quality_overview(
    quality_reports: Dict[str, Dict[str, Any]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (14, 8),
    **kwargs,
) -> Axes:
    """Plot quality overview across multiple omics layers.

    Args:
        quality_reports: Dictionary mapping omics types to quality reports
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(quality_reports, dict, "quality_reports")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    omics_types = list(quality_reports.keys())

    # Extract quality scores for each omics type
    quality_scores = []
    for omics_type in omics_types:
        report = quality_reports[omics_type]
        # Try different possible quality score keys
        score = (
            report.get("overall_quality")
            or report.get("quality_score")
            or report.get("mean_quality")
            or np.mean(list(report.values()))
            if report
            else 0
        )
        quality_scores.append(score)

    # Create radar chart
    angles = np.linspace(0, 2 * np.pi, len(omics_types), endpoint=False).tolist()
    angles += angles[:1]  # Close the plot
    quality_scores += quality_scores[:1]  # Close the plot

    ax.plot(angles, quality_scores, "o-", linewidth=2, label="Quality Score", alpha=0.8)
    ax.fill(angles, quality_scores, alpha=0.25)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(omics_types)
    ax.set_ylim(0, max(quality_scores) * 1.1)
    ax.set_title("Multi-Omics Quality Overview")
    ax.grid(True, alpha=0.3)

    # Add quality score labels
    for angle, score, omics in zip(angles[:-1], quality_scores[:-1], omics_types):
        ax.text(angle, score + max(quality_scores) * 0.05, f"{score:.2f}", ha="center", va="bottom", fontsize=8)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Multi-omics quality overview saved to {output_path}")

    return ax


def plot_coverage_uniformity(
    coverage_data: np.ndarray,
    positions: np.ndarray | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes:
    """Plot sequencing coverage uniformity.

    Args:
        coverage_data: Coverage values across positions
        positions: Optional genomic positions
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(coverage_data, np.ndarray, "coverage_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    x_data = positions if positions is not None else np.arange(len(coverage_data))

    # Plot coverage
    ax.plot(x_data, coverage_data, "b-", linewidth=1, alpha=0.8)
    ax.fill_between(x_data, coverage_data, alpha=0.3, color="blue")

    # Add statistics
    mean_cov = np.mean(coverage_data)
    median_cov = np.median(coverage_data)
    std_cov = np.std(coverage_data)

    ax.axhline(y=mean_cov, color="red", linestyle="--", alpha=0.7, label=f"Mean: {mean_cov:.1f}")
    ax.axhline(y=median_cov, color="green", linestyle="--", alpha=0.7, label=f"Median: {median_cov:.1f}")

    ax.set_xlabel("Position" if positions is None else "Genomic Position")
    ax.set_ylabel("Coverage Depth")
    ax.set_title("Coverage Uniformity")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add statistics text box
    stats_text = (
        f"CV: {std_cov/mean_cov:.3f}\nIQR: {np.percentile(coverage_data, 75) - np.percentile(coverage_data, 25):.1f}"
    )
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Coverage uniformity plot saved to {output_path}")

    return ax


def plot_error_profiles(
    error_profiles: Dict[str, np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes:
    """Plot error profiles across different error types.

    Args:
        error_profiles: Dictionary mapping error types to error rates
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(error_profiles, dict, "error_profiles")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    colors = plt.cm.tab10(np.linspace(0, 1, len(error_profiles)))

    for i, (error_type, profile) in enumerate(error_profiles.items()):
        positions = np.arange(len(profile))
        ax.plot(positions, profile, color=colors[i], linewidth=2, label=error_type, alpha=0.8)

    ax.set_xlabel("Position in Read")
    ax.set_ylabel("Error Rate")
    ax.set_title("Error Profiles")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_yscale("log")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Error profiles plot saved to {output_path}")

    return ax


def plot_batch_effects_qc(
    batch_qc_data: Dict[str, Any],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (14, 10),
    **kwargs,
) -> Axes:
    """Plot batch effects quality control analysis.

    Args:
        batch_qc_data: Dictionary containing batch effect analysis results
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(batch_qc_data, dict, "batch_qc_data")

    if ax is None:
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()
    else:
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()

    plot_idx = 0

    # Plot 1: Batch distribution
    if "batch_sizes" in batch_qc_data:
        batches = list(batch_qc_data["batch_sizes"].keys())
        sizes = list(batch_qc_data["batch_sizes"].values())
        axes[plot_idx].bar(batches, sizes, alpha=0.7)
        axes[plot_idx].set_xlabel("Batch")
        axes[plot_idx].set_ylabel("Sample Count")
        axes[plot_idx].set_title("Batch Size Distribution")
        axes[plot_idx].tick_params(axis="x", rotation=45)
        axes[plot_idx].grid(True, alpha=0.3, axis="y")
        plot_idx += 1

    # Plot 2: PC1 vs PC2 colored by batch
    if "batch_pca" in batch_qc_data:
        pca_data = batch_qc_data["batch_pca"]
        if "pc1" in pca_data and "pc2" in pca_data and "batches" in pca_data:
            unique_batches = list(set(pca_data["batches"]))
            colors = plt.cm.tab10(np.linspace(0, 1, len(unique_batches)))

            for i, batch in enumerate(unique_batches):
                mask = np.array(pca_data["batches"]) == batch
                axes[plot_idx].scatter(
                    np.array(pca_data["pc1"])[mask],
                    np.array(pca_data["pc2"])[mask],
                    color=colors[i],
                    label=batch,
                    alpha=0.7,
                    s=50,
                )

            axes[plot_idx].set_xlabel("PC1")
            axes[plot_idx].set_ylabel("PC2")
            axes[plot_idx].set_title("Batch Effects in PCA")
            axes[plot_idx].legend()
            axes[plot_idx].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot 3: Silhouette scores by batch
    if "silhouette_scores" in batch_qc_data:
        silhouette_data = batch_qc_data["silhouette_scores"]
        if isinstance(silhouette_data, dict):
            batches = list(silhouette_data.keys())
            scores = list(silhouette_data.values())
            bars = axes[plot_idx].bar(batches, scores, alpha=0.7, color="orange")
            axes[plot_idx].set_xlabel("Batch")
            axes[plot_idx].set_ylabel("Silhouette Score")
            axes[plot_idx].set_title("Batch Separation Quality")
            axes[plot_idx].tick_params(axis="x", rotation=45)
            axes[plot_idx].grid(True, alpha=0.3, axis="y")

            # Add score labels
            for bar, score in zip(bars, scores):
                axes[plot_idx].text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.01,
                    f"{score:.3f}",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )
        plot_idx += 1

    # Plot 4: Differential expression by batch
    if "batch_de_stats" in batch_qc_data:
        de_data = batch_qc_data["batch_de_stats"]
        if "n_de_genes" in de_data and "batches" in de_data:
            batches = de_data["batches"]
            n_de = de_data["n_de_genes"]
            axes[plot_idx].bar(batches, n_de, alpha=0.7, color="red")
            axes[plot_idx].set_xlabel("Batch Comparison")
            axes[plot_idx].set_ylabel("DE Genes")
            axes[plot_idx].set_title("Batch-Induced DE Genes")
            axes[plot_idx].tick_params(axis="x", rotation=45)
            axes[plot_idx].grid(True, alpha=0.3, axis="y")
        plot_idx += 1

    # Plot 5: Batch effect variance explained
    if "batch_variance" in batch_qc_data:
        variance_data = batch_qc_data["batch_variance"]
        if isinstance(variance_data, dict):
            components = list(variance_data.keys())
            variances = list(variance_data.values())
            axes[plot_idx].pie(variances, labels=components, autopct="%1.1f%%", startangle=90)
            axes[plot_idx].set_title("Variance Explained by Batch")
        plot_idx += 1

    # Plot 6: Batch correction effectiveness
    if "correction_metrics" in batch_qc_data:
        metrics = batch_qc_data["correction_metrics"]
        if isinstance(metrics, dict):
            metric_names = list(metrics.keys())
            values = list(metrics.values())
            bars = axes[plot_idx].bar(metric_names, values, alpha=0.7, color="purple")
            axes[plot_idx].set_ylabel("Metric Value")
            axes[plot_idx].set_title("Batch Correction Effectiveness")
            axes[plot_idx].tick_params(axis="x", rotation=45)
            axes[plot_idx].grid(True, alpha=0.3, axis="y")

    plt.tight_layout()

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Batch effects QC plot saved to {output_path}")

    return axes[0]


def plot_data_integrity_metrics(
    integrity_metrics: Dict[str, Any],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes:
    """Plot data integrity and completeness metrics.

    Args:
        integrity_metrics: Dictionary containing data integrity metrics
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(integrity_metrics, dict, "integrity_metrics")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Extract metrics
    metrics = []
    values = []
    colors = []

    for key, value in integrity_metrics.items():
        if isinstance(value, (int, float)):
            metrics.append(key)
            values.append(value)
            # Color based on value (green for good, red for bad)
            if "missing" in key.lower() or "error" in key.lower():
                colors.append("red" if value > 0.1 else "orange" if value > 0.01 else "green")
            elif "complete" in key.lower() or "valid" in key.lower():
                colors.append("green" if value > 0.9 else "orange" if value > 0.7 else "red")
            else:
                colors.append("blue")

    x_positions = np.arange(len(metrics))
    bars = ax.bar(x_positions, values, color=colors, alpha=0.8, **kwargs)

    ax.set_xlabel("Integrity Metric")
    ax.set_ylabel("Value")
    ax.set_title("Data Integrity Assessment")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(metrics, rotation=45, ha="right")
    ax.grid(True, alpha=0.3, axis="y")

    # Add value labels
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.01,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Data integrity metrics plot saved to {output_path}")

    return ax
