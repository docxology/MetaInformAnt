"""Quality assessment visualization for coverage, errors, and batch effects.

Functions for visualizing coverage uniformity, error profiles,
batch effect diagnostics, and data integrity metrics.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging
from metainformant.visualization.config.conventions import add_log_zero_mark, save_figure_deterministic
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


def plot_coverage_uniformity(
    coverage_data: np.ndarray,
    positions: np.ndarray | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs: Any,
) -> Axes:
    """Plot sequencing coverage uniformity.

    Args:
        coverage_data: Coverage values across positions.
        positions: Optional genomic positions.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.
        figsize: Figure size.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(coverage_data, np.ndarray, "coverage_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    x_data = positions if positions is not None else np.arange(len(coverage_data))

    ax.plot(x_data, coverage_data, "b-", linewidth=1, alpha=0.8)
    ax.fill_between(x_data, coverage_data, alpha=0.3, color="blue")

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

    _save_figure(ax, output_path, f"Coverage uniformity plot saved to {output_path}")

    return ax


def plot_error_profiles(
    error_profiles: Dict[str, np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs: Any,
) -> Axes:
    """Plot error profiles across different error types.

    Args:
        error_profiles: Dictionary mapping error types to error rates.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.
        figsize: Figure size.

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
    if any(np.any(np.asarray(profile) == 0) for profile in error_profiles.values()):
        # a log axis cannot show zero: without a mark an error-free position would sit at the axis floor unexplained
        add_log_zero_mark(ax, "y")

    _save_figure(ax, output_path, f"Error profiles plot saved to {output_path}")

    return ax


#: Sections of ``batch_qc_data`` in drawing order, with the keys each section's payload must have (``None`` = any
#: mapping of label -> number).
_BATCH_SECTIONS: Dict[str, Tuple[str, ...] | None] = {
    "batch_sizes": None,
    "batch_pca": ("pc1", "pc2", "batches"),
    "silhouette_scores": None,
    "batch_de_stats": ("n_de_genes", "batches"),
    "batch_variance": None,
    "correction_metrics": None,
}


def plot_batch_effects_qc(
    batch_qc_data: Dict[str, Any],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (14, 10),
    **kwargs: Any,
) -> Axes:
    """Plot batch effects quality control analysis.

    One panel is drawn for each section present in ``batch_qc_data`` (a compact grid of at most three columns, no
    empty panels), in this order:

    * ``batch_sizes``: ``{batch: sample count}``, bars
    * ``batch_pca``: ``{"pc1": [...], "pc2": [...], "batches": [...]}``, a PCA scatter coloured by batch
    * ``silhouette_scores``: ``{batch: score}``, bars
    * ``batch_de_stats``: ``{"batches": [comparison labels], "n_de_genes": [counts]}``, bars
    * ``batch_variance``: ``{component: variance}``, a pie chart
    * ``correction_metrics``: ``{metric: value}``, bars

    Args:
        batch_qc_data: Dictionary containing batch effect analysis results.
        ax: Optional matplotlib axes that should hold the panels: they are laid out inside its grid slot and
            ``ax`` is hidden. Creates a new figure if None.
        output_path: Optional path to save.
        figsize: Figure size of a new figure.

    Returns:
        The axes of the first panel

    Raises:
        ValueError: If no recognised section is present, or a section is not shaped as described above
    """
    validation.validate_type(batch_qc_data, dict, "batch_qc_data")

    sections = [name for name in _BATCH_SECTIONS if name in batch_qc_data]
    if not sections:
        raise ValueError(f"batch_qc_data has none of the recognised sections: {list(_BATCH_SECTIONS)}")
    for name in sections:
        payload, required = batch_qc_data[name], _BATCH_SECTIONS[name]
        if not isinstance(payload, dict):
            raise ValueError(f"batch_qc_data[{name!r}] must be a dict, got {type(payload).__name__}")
        if required is not None and not all(key in payload for key in required):
            raise ValueError(f"batch_qc_data[{name!r}] needs the keys {list(required)}; got {sorted(payload)}")

    axes = panel_grid(ax, len(sections), max_cols=3, figsize=figsize)

    for panel, name in zip(axes, sections):
        payload = batch_qc_data[name]
        if name == "batch_sizes":
            panel.bar(list(payload.keys()), list(payload.values()), alpha=0.7)
            panel.set_xlabel("Batch")
            panel.set_ylabel("Sample Count")
            panel.set_title("Batch Size Distribution")
            panel.tick_params(axis="x", rotation=45)
            panel.grid(True, alpha=0.3, axis="y")
        elif name == "batch_pca":
            batch_labels = np.array(payload["batches"])
            pc1, pc2 = np.array(payload["pc1"]), np.array(payload["pc2"])
            unique_batches = list(dict.fromkeys(payload["batches"]))  # first-seen order: a set orders strings by hash
            colors = plt.cm.tab10(np.linspace(0, 1, len(unique_batches)))
            for i, batch in enumerate(unique_batches):
                mask = batch_labels == batch
                panel.scatter(pc1[mask], pc2[mask], color=colors[i], label=batch, alpha=0.7, s=50)
            panel.set_xlabel("PC1")
            panel.set_ylabel("PC2")
            panel.set_title("Batch Effects in PCA")
            panel.legend()
            panel.grid(True, alpha=0.3)
        elif name == "silhouette_scores":
            batches, scores = list(payload.keys()), list(payload.values())
            bars = panel.bar(batches, scores, alpha=0.7, color="orange")
            panel.set_xlabel("Batch")
            panel.set_ylabel("Silhouette Score")
            panel.set_title("Batch Separation Quality")
            panel.tick_params(axis="x", rotation=45)
            panel.grid(True, alpha=0.3, axis="y")
            _annotate_bars(panel, bars, scores, "{:.3f}", offset=0.01)
        elif name == "batch_de_stats":
            panel.bar(payload["batches"], payload["n_de_genes"], alpha=0.7, color="red")
            panel.set_xlabel("Batch Comparison")
            panel.set_ylabel("DE Genes")
            panel.set_title("Batch-Induced DE Genes")
            panel.tick_params(axis="x", rotation=45)
            panel.grid(True, alpha=0.3, axis="y")
        elif name == "batch_variance":
            panel.pie(list(payload.values()), labels=list(payload.keys()), autopct="%1.1f%%", startangle=90)
            panel.set_title("Variance Explained by Batch")
        else:  # correction_metrics
            panel.bar(list(payload.keys()), list(payload.values()), alpha=0.7, color="purple")
            panel.set_ylabel("Metric Value")
            panel.set_title("Batch Correction Effectiveness")
            panel.tick_params(axis="x", rotation=45)
            panel.grid(True, alpha=0.3, axis="y")

    if ax is None:
        root_figure(axes[0]).tight_layout()

    _save_figure(axes[0], output_path, f"Batch effects QC plot saved to {output_path}")

    return axes[0]


def plot_data_integrity_metrics(
    integrity_metrics: Dict[str, Any],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs: Any,
) -> Axes:
    """Plot data integrity and completeness metrics.

    Args:
        integrity_metrics: Dictionary containing data integrity metrics.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.
        figsize: Figure size.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(integrity_metrics, dict, "integrity_metrics")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    metrics = []
    values = []
    colors = []

    for key, value in integrity_metrics.items():
        if isinstance(value, (int, float)):
            metrics.append(key)
            values.append(value)
            if "missing" in key.lower() or "error" in key.lower():
                colors.append("red" if value > 0.1 else "orange" if value > 0.01 else "green")
            elif "complete" in key.lower() or "valid" in key.lower():
                colors.append("green" if value > 0.9 else "orange" if value > 0.7 else "red")
            else:
                colors.append("blue")

    x_positions = np.arange(len(metrics))
    # caller keyword arguments (color, alpha, ...) override the defaults instead of colliding with them
    bars = ax.bar(x_positions, values, **{"color": colors, "alpha": 0.8, **kwargs})

    ax.set_xlabel("Integrity Metric")
    ax.set_ylabel("Value")
    ax.set_title("Data Integrity Assessment")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(metrics, rotation=45, ha="right")
    ax.grid(True, alpha=0.3, axis="y")

    _annotate_bars(ax, bars, values, "{:.3f}", offset=0.01)

    _save_figure(ax, output_path, f"Data integrity metrics plot saved to {output_path}")

    return ax
