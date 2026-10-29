"""Sequencing quality control visualization functions.

Functions for visualizing read-level and sequence-level QC metrics
including quality scores, GC content, read lengths, adapter content,
duplication, overrepresented sequences, and k-mer profiles.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging
from metainformant.visualization.config.conventions import save_figure_deterministic
from metainformant.visualization.config.figures import panel_grid, root_figure

logger = logging.get_logger(__name__)


def _save_figure(ax: Axes, output_path: str | Path | None, message: str) -> None:
    """Persist the plotted figure (``ax.figure``) to *output_path* when requested."""
    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(message)


def _annotate_bars(ax: Axes, bars: Iterable[Any], values: Iterable[Any], fmt: str, offset: float = 0.0) -> None:
    """Write *fmt*-formatted value labels centered above each bar."""
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + offset,
            fmt.format(value),
            ha="center",
            va="bottom",
            fontsize=8,
        )


#: Sections of ``qc_data`` in drawing order, with the keys each section's payload must have.
_QC_SECTIONS: Dict[str, Tuple[str, ...]] = {
    "per_base_quality": ("positions", "mean_qualities"),
    "gc_content_distribution": ("bins", "counts"),
    "sequence_length_distribution": ("lengths", "counts"),
    "basic_statistics": (),
}

#: The statistics ``basic_statistics`` may hold, in plotting order.
_BASIC_STATISTICS = ("num_reads", "total_bases", "min_length", "max_length", "mean_length")


def plot_quality_metrics(
    qc_data: Dict[str, Any], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a comprehensive quality metrics visualization.

    One panel is drawn for each section present in ``qc_data`` (a compact grid of at most two columns, no empty
    panels), in this order:

    * ``per_base_quality``: ``{"positions": [...], "mean_qualities": [...]}``, a line
    * ``gc_content_distribution``: ``{"bins": [bin edges], "counts": [one per bin]}``, bars as wide as their bins
    * ``sequence_length_distribution``: ``{"lengths": [...], "counts": [...]}``, bars
    * ``basic_statistics``: any of ``num_reads``, ``total_bases``, ``min_length``, ``max_length``, ``mean_length``;
      a statistic that is absent is not drawn (it used to be drawn as a bar of 0)

    Args:
        qc_data: Dictionary containing various QC metrics.
        ax: Optional matplotlib axes that should hold the panels: they are laid out inside its grid slot and
            ``ax`` is hidden. Creates a new figure if None.
        output_path: Optional path to save the figure.
        **kwargs: ``figsize`` sets the size of a new figure.

    Returns:
        The axes of the first panel

    Raises:
        ValueError: If ``qc_data`` is empty, has no recognised section, or a section is not shaped as above
    """
    validation.validate_type(qc_data, dict, "qc_data")
    if not qc_data:
        raise ValueError("QC data dictionary cannot be empty")

    sections = [name for name in _QC_SECTIONS if name in qc_data]
    if not sections:
        raise ValueError(f"qc_data has none of the recognised sections: {list(_QC_SECTIONS)}")
    for name in sections:
        payload = qc_data[name]
        if not (isinstance(payload, dict) and all(key in payload for key in _QC_SECTIONS[name])):
            raise ValueError(f"qc_data[{name!r}] must be a dict with the keys {list(_QC_SECTIONS[name])}")
    if "basic_statistics" in sections and not any(stat in qc_data["basic_statistics"] for stat in _BASIC_STATISTICS):
        raise ValueError(f"qc_data['basic_statistics'] has none of the statistics {list(_BASIC_STATISTICS)}")
    if "gc_content_distribution" in sections:
        gc = qc_data["gc_content_distribution"]
        if len(gc["counts"]) != len(gc["bins"]) - 1:
            raise ValueError(
                f"gc_content_distribution needs one count per bin: {len(gc['bins'])} bin edges give "
                f"{len(gc['bins']) - 1} bins, got {len(gc['counts'])} counts"
            )

    axes = panel_grid(ax, len(sections), max_cols=2, figsize=kwargs.pop("figsize", (12, 10)))

    for panel, name in zip(axes, sections):
        payload = qc_data[name]
        if name == "per_base_quality":
            panel.plot(payload["positions"], payload["mean_qualities"])
            panel.set_xlabel("Position in Read")
            panel.set_ylabel("Mean Quality Score")
            panel.set_title("Per-Base Quality")
            panel.grid(True, alpha=0.3)
        elif name == "gc_content_distribution":
            edges = np.asarray(payload["bins"], dtype=float)
            panel.bar((edges[:-1] + edges[1:]) / 2, payload["counts"], width=np.diff(edges), alpha=0.7)
            panel.set_xlabel("GC Content (%)")
            panel.set_ylabel("Frequency")
            panel.set_title("GC Content Distribution")
            panel.grid(True, alpha=0.3)
        elif name == "sequence_length_distribution":
            panel.bar(payload["lengths"], payload["counts"], alpha=0.7)
            panel.set_xlabel("Read Length")
            panel.set_ylabel("Count")
            panel.set_title("Read Length Distribution")
            panel.grid(True, alpha=0.3)
        else:  # basic_statistics
            stat_names = [stat for stat in _BASIC_STATISTICS if stat in payload]
            stat_values = [payload[stat] for stat in stat_names]
            bars = panel.bar(stat_names, stat_values, alpha=0.7)
            panel.set_ylabel("Value")
            panel.set_title("Basic Statistics")
            panel.tick_params(axis="x", rotation=45)
            _annotate_bars(panel, bars, stat_values, "{:.0f}")

    if ax is None:
        root_figure(axes[0]).tight_layout()

    _save_figure(axes[0], output_path, f"Quality metrics plot saved to {output_path}")

    return axes[0]


def plot_adapter_content(
    adapter_data: Dict[str, List[float]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
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

    bars = ax.bar(range(len(adapters_sorted)), percentages_sorted, **{"alpha": 0.7, **kwargs})
    ax.set_xlabel("Adapter Type")
    ax.set_ylabel("Maximum Content (%)")
    ax.set_title("Adapter Content Analysis")
    ax.set_xticks(range(len(adapters_sorted)))
    ax.set_xticklabels(adapters_sorted, rotation=45, ha="right")
    ax.grid(True, alpha=0.3)

    _annotate_bars(ax, bars, percentages_sorted, "{:.1f}%")

    _save_figure(ax, output_path, f"Adapter content plot saved to {output_path}")

    return ax


def plot_gc_distribution(
    gc_data: List[float], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
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

    n, bins, patches = ax.hist(gc_data, **{"bins": 20, "alpha": 0.7, "edgecolor": "black", **kwargs})
    ax.set_xlabel("GC Content (%)")
    ax.set_ylabel("Frequency")
    ax.set_title("GC Content Distribution")
    ax.grid(True, alpha=0.3)

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

    _save_figure(ax, output_path, f"GC distribution plot saved to {output_path}")

    return ax


def plot_length_distribution(
    length_data: List[int], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
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

    n, bins, patches = ax.hist(length_data, **{"bins": "auto", "alpha": 0.7, "edgecolor": "black", **kwargs})
    ax.set_xlabel("Read Length (bp)")
    ax.set_ylabel("Frequency")
    ax.set_title("Read Length Distribution")
    ax.grid(True, alpha=0.3)

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

    _save_figure(ax, output_path, f"Length distribution plot saved to {output_path}")

    return ax


def plot_per_base_quality_boxplot(
    per_base_qualities: Dict[str, List[float]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs: Any,
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

    bp = ax.boxplot(qualities, **{"positions": positions, "widths": 2, "patch_artist": True, **kwargs})
    colors = plt.cm.RdYlGn(np.linspace(0, 1, len(qualities)))
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xlabel("Position in Read")
    ax.set_ylabel("Quality Score")
    ax.set_title("Per-Base Quality Scores")
    ax.grid(True, alpha=0.3)

    _save_figure(ax, output_path, f"Per-base quality boxplot saved to {output_path}")

    return ax


def plot_sequence_duplication_levels(
    duplication_levels: Dict[str, float],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 6),
    **kwargs: Any,
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

    bars = ax.bar(range(len(levels)), percentages, **{"color": "skyblue", "alpha": 0.8, **kwargs})
    ax.set_xlabel("Duplication Level")
    ax.set_ylabel("Percentage (%)")
    ax.set_title("Sequence Duplication Levels")
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels(levels, rotation=45, ha="right")
    ax.grid(True, alpha=0.3, axis="y")

    _annotate_bars(ax, bars, percentages, "{:.1f}%", offset=0.5)

    _save_figure(ax, output_path, f"Sequence duplication levels plot saved to {output_path}")

    return ax


def plot_overrepresented_sequences(
    overrepresented_seqs: List[Dict[str, Any]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs: Any,
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
        # the figure is the result the caller asked to have written: an empty finding is still a finding
        _save_figure(ax, output_path, f"Overrepresented sequences plot saved to {output_path}")
        return ax

    full_sequences = [str(seq.get("sequence", f"Seq {i+1}")) for i, seq in enumerate(overrepresented_seqs)]
    sequences = [seq[:20] + "..." if len(seq) > 20 else seq for seq in full_sequences]  # "..." only when cut
    counts = [seq.get("count", 0) for seq in overrepresented_seqs]
    percentages = [seq.get("percentage", 0) for seq in overrepresented_seqs]

    y_pos = np.arange(len(sequences))
    bars = ax.barh(y_pos, percentages, **{"color": "coral", "alpha": 0.8, **kwargs})
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sequences)
    ax.set_xlabel("Percentage of Total Reads")
    ax.set_title("Overrepresented Sequences")
    ax.grid(True, alpha=0.3, axis="x")

    for i, (bar, pct, count) in enumerate(zip(bars, percentages, counts)):
        ax.text(
            bar.get_width() + 0.1,
            bar.get_y() + bar.get_height() / 2,
            f"{pct:.2f}%\n({count})",
            ha="left",
            va="center",
            fontsize=8,
        )

    _save_figure(ax, output_path, f"Overrepresented sequences plot saved to {output_path}")

    return ax


def plot_kmer_profiles(
    kmer_counts: Dict[str, int],
    top_n: int = 20,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs: Any,
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

    bars = ax.bar(range(len(kmers)), counts, **{"color": "purple", "alpha": 0.8, **kwargs})
    ax.set_xlabel("K-mer")
    ax.set_ylabel("Count")
    ax.set_title(f"Top {top_n} K-mer Frequencies")
    ax.set_xticks(range(len(kmers)))
    ax.set_xticklabels(kmers, rotation=45, ha="right")
    ax.grid(True, alpha=0.3, axis="y")

    _annotate_bars(ax, bars, counts, "{}", offset=max(counts) * 0.01)

    _save_figure(ax, output_path, f"K-mer profiles plot saved to {output_path}")

    return ax
