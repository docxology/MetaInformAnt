"""Composite and multi-panel figure builder for overview visualizations.

Provides utilities to compose multiple plot types into a single
publication-ready figure, including genomic overview dashboards
and QC summary panels.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from metainformant.core.io import paths
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


# ---------------------------------------------------------------------------
# Generic multi-panel builder
# ---------------------------------------------------------------------------


def multi_panel(
    panels: List[Dict[str, Any]],
    *,
    ncols: int = 2,
    figsize: Tuple[float, float] | None = None,
    title: str | None = None,
    output_path: str | Path | None = None,
) -> Figure:
    """Arrange multiple plots into a single figure.

    Each element in *panels* is a dict with:
        - ``"func"`` — a callable ``(ax, **kwargs) -> Axes``
        - ``"kwargs"`` — keyword arguments forwarded to the callable
        - ``"title"`` (optional) — subplot title

    Args:
        panels: List of panel specification dicts.
        ncols: Number of columns in the grid.
        figsize: Overall figure size (auto-calculated if *None*).
        title: Super-title for the figure.
        output_path: Path to save the figure.

    Returns:
        The :class:`~matplotlib.figure.Figure` object.
    """
    n = len(panels)
    nrows = max(1, -(-n // ncols))  # ceil division
    if figsize is None:
        figsize = (6 * ncols, 5 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    if n == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if hasattr(axes, "flatten") else [axes]

    for i, panel in enumerate(panels):
        ax = axes[i]
        func = panel["func"]
        kwargs = panel.get("kwargs", {})
        func(ax=ax, **kwargs)
        if "title" in panel:
            ax.set_title(panel["title"])

    # Hide unused axes
    for j in range(n, len(axes)):
        axes[j].set_visible(False)

    if title:
        fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout()

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Multi-panel figure saved to {output_path}")

    return fig


# ---------------------------------------------------------------------------
# Genomic overview dashboard
# ---------------------------------------------------------------------------


def genomic_overview(
    data: Dict[str, Any],
    *,
    figsize: Tuple[float, float] = (16, 12),
    output_path: str | Path | None = None,
) -> Figure:
    """Create a genomic overview dashboard with up to 6 panels.

    Expected keys in *data* (all optional — panels are skipped when absent):
        - ``"expression"`` — DataFrame of gene expression values
        - ``"pca_data"`` — 2-D array suitable for PCA scatter
        - ``"pvalues"`` — array of p-values for a QQ plot
        - ``"fold_changes"`` / ``"log2fc"`` + ``"pvalues"`` — volcano data
        - ``"gc_content"`` — list/array of GC percentages
        - ``"quality_scores"`` — list/array of quality scores

    Args:
        data: Data dictionary keyed by panel type.
        figsize: Overall figure size.
        output_path: Path to save the figure.

    Returns:
        The :class:`~matplotlib.figure.Figure` object.
    """
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)
    panel_idx = 0

    # Panel 1 — expression heatmap
    if "expression" in data:
        ax = fig.add_subplot(gs[0, 0])
        expr = data["expression"]
        if isinstance(expr, pd.DataFrame):
            im = ax.imshow(expr.values[:20, :20], cmap="RdBu_r", aspect="auto")
            ax.set_title("Expression Heatmap (top 20)")
            plt.colorbar(im, ax=ax, fraction=0.046)
        panel_idx += 1

    # Panel 2 — PCA scatter
    if "pca_data" in data:
        ax = fig.add_subplot(gs[0, 1])
        pca_arr = np.asarray(data["pca_data"])
        if pca_arr.ndim == 2 and pca_arr.shape[1] >= 2:
            ax.scatter(pca_arr[:, 0], pca_arr[:, 1], s=20, alpha=0.7)
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
            ax.set_title("PCA")
        panel_idx += 1

    # Panel 3 — QQ plot
    if "pvalues" in data:
        ax = fig.add_subplot(gs[0, 2])
        pvals = np.asarray(data["pvalues"])
        pvals = pvals[pvals > 0]
        observed = -np.log10(np.sort(pvals))
        expected = -np.log10(np.linspace(1 / len(pvals), 1, len(pvals)))
        ax.scatter(expected, observed, s=10, alpha=0.6)
        lim = max(observed.max(), expected.max()) * 1.05
        ax.plot([0, lim], [0, lim], "r--", alpha=0.7)
        ax.set_xlabel("Expected -log10(p)")
        ax.set_ylabel("Observed -log10(p)")
        ax.set_title("QQ Plot")
        panel_idx += 1

    # Panel 4 — volcano plot
    if "log2fc" in data and "pvalues" in data:
        ax = fig.add_subplot(gs[1, 0])
        log2fc = np.asarray(data["log2fc"])
        pvals = np.asarray(data["pvalues"])
        neg_log_p = -np.log10(np.clip(pvals, 1e-300, 1))
        colors = np.where(
            (np.abs(log2fc) > 1) & (pvals < 0.05),
            "#D32F2F",
            np.where(pvals < 0.05, "#FF9800", "#9E9E9E"),
        )
        ax.scatter(log2fc, neg_log_p, c=colors, s=10, alpha=0.6)
        ax.set_xlabel("log2 Fold Change")
        ax.set_ylabel("-log10(p)")
        ax.set_title("Volcano Plot")
        panel_idx += 1

    # Panel 5 — GC content
    if "gc_content" in data:
        ax = fig.add_subplot(gs[1, 1])
        gc = np.asarray(data["gc_content"])
        ax.hist(gc, bins=30, alpha=0.7, edgecolor="black")
        ax.set_xlabel("GC Content (%)")
        ax.set_ylabel("Frequency")
        ax.set_title("GC Distribution")
        panel_idx += 1

    # Panel 6 — quality scores
    if "quality_scores" in data:
        ax = fig.add_subplot(gs[1, 2])
        qs = np.asarray(data["quality_scores"])
        ax.hist(qs, bins=30, alpha=0.7, color="green", edgecolor="black")
        ax.set_xlabel("Quality Score")
        ax.set_ylabel("Frequency")
        ax.set_title("Quality Distribution")
        panel_idx += 1

    fig.suptitle("Genomic Overview Dashboard", fontsize=16, fontweight="bold")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Genomic overview dashboard saved to {output_path}")

    return fig


# ---------------------------------------------------------------------------
# QC summary dashboard
# ---------------------------------------------------------------------------


def qc_summary(
    qc_data: Dict[str, Any],
    *,
    figsize: Tuple[float, float] = (14, 10),
    output_path: str | Path | None = None,
) -> Figure:
    """Create a QC summary dashboard.

    Expected keys in *qc_data* (all optional):
        - ``"read_lengths"`` — array of read lengths
        - ``"gc_content"`` — array of GC percentages
        - ``"quality_scores"`` — array of per-read mean quality
        - ``"mapping_rates"`` — dict of sample -> mapping rate
        - ``"duplication_rates"`` — dict of sample -> duplication rate
        - ``"adapter_content"`` — dict of adapter -> max percentage

    Args:
        qc_data: Data dictionary.
        figsize: Overall figure size.
        output_path: Path to save the figure.

    Returns:
        The :class:`~matplotlib.figure.Figure` object.
    """
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes = axes.flatten()

    idx = 0

    if "read_lengths" in qc_data:
        axes[idx].hist(qc_data["read_lengths"], bins=30, alpha=0.7, edgecolor="black")
        axes[idx].set_xlabel("Read Length (bp)")
        axes[idx].set_ylabel("Count")
        axes[idx].set_title("Read Length Distribution")
        idx += 1

    if "gc_content" in qc_data:
        axes[idx].hist(qc_data["gc_content"], bins=30, alpha=0.7, color="teal", edgecolor="black")
        axes[idx].set_xlabel("GC Content (%)")
        axes[idx].set_ylabel("Count")
        axes[idx].set_title("GC Distribution")
        idx += 1

    if "quality_scores" in qc_data:
        axes[idx].hist(qc_data["quality_scores"], bins=30, alpha=0.7, color="green", edgecolor="black")
        axes[idx].set_xlabel("Mean Quality")
        axes[idx].set_ylabel("Count")
        axes[idx].set_title("Quality Score Distribution")
        idx += 1

    if "mapping_rates" in qc_data:
        samples = list(qc_data["mapping_rates"].keys())
        rates = list(qc_data["mapping_rates"].values())
        axes[idx].barh(samples, rates, alpha=0.8, color="steelblue")
        axes[idx].set_xlabel("Mapping Rate (%)")
        axes[idx].set_title("Mapping Rates")
        axes[idx].set_xlim(0, 100)
        idx += 1

    if "duplication_rates" in qc_data:
        samples = list(qc_data["duplication_rates"].keys())
        rates = list(qc_data["duplication_rates"].values())
        axes[idx].barh(samples, rates, alpha=0.8, color="coral")
        axes[idx].set_xlabel("Duplication Rate (%)")
        axes[idx].set_title("Duplication Rates")
        idx += 1

    if "adapter_content" in qc_data:
        adapters = list(qc_data["adapter_content"].keys())
        pcts = list(qc_data["adapter_content"].values())
        axes[idx].bar(adapters, pcts, alpha=0.8, color="goldenrod")
        axes[idx].set_ylabel("Max Content (%)")
        axes[idx].set_title("Adapter Content")
        axes[idx].tick_params(axis="x", rotation=45)
        idx += 1

    for j in range(idx, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("QC Summary Dashboard", fontsize=14, fontweight="bold")
    fig.tight_layout()

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"QC summary dashboard saved to {output_path}")

    return fig
