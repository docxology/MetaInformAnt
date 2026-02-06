"""eQTL visualization module for expression-variant analysis.

Provides specialized plots for eQTL results including volcano plots,
genotype-expression boxplots, and LocusZoom-style regional plots.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

try:
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    plt = None


def plot_eqtl_volcano(
    results: pd.DataFrame,
    fdr_threshold: float = 0.05,
    effect_size_threshold: float = 0.5,
    output_path: Path | str | None = None,
    title: str = "eQTL Volcano Plot",
    **kwargs: Any,
) -> Any:
    """Create volcano plot of eQTL associations.

    Args:
        results: DataFrame with beta, pvalue columns.
        fdr_threshold: FDR threshold for significance line.
        effect_size_threshold: Effect size threshold for vertical lines.
        output_path: Path to save figure (optional).
        title: Plot title.
        **kwargs: Additional matplotlib kwargs.

    Returns:
        Matplotlib figure object.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for plotting")
        return None

    fig, ax = plt.subplots(figsize=kwargs.get("figsize", (10, 8)))

    # Calculate -log10(p)
    log_pval = -np.log10(results["pvalue"].clip(lower=1e-300))
    beta = results["beta"]

    # Color by significance
    colors = np.where(
        (log_pval > -np.log10(fdr_threshold)) & (np.abs(beta) > effect_size_threshold),
        "red",
        np.where(log_pval > -np.log10(fdr_threshold), "orange", "gray"),
    )

    ax.scatter(beta, log_pval, c=colors, alpha=0.6, s=20)

    # Add threshold lines
    ax.axhline(-np.log10(fdr_threshold), color="red", linestyle="--", alpha=0.5)
    ax.axvline(-effect_size_threshold, color="blue", linestyle="--", alpha=0.3)
    ax.axvline(effect_size_threshold, color="blue", linestyle="--", alpha=0.3)

    ax.set_xlabel("Effect Size (β)")
    ax.set_ylabel("-log₁₀(P-value)")
    ax.set_title(title)

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved volcano plot to {output_path}")

    return fig


def plot_eqtl_boxplot(
    expression: np.ndarray,
    genotypes: np.ndarray,
    gene_id: str = "",
    variant_id: str = "",
    output_path: Path | str | None = None,
    **kwargs: Any,
) -> Any:
    """Create genotype vs expression boxplot.

    Args:
        expression: Gene expression values.
        genotypes: Genotype dosages (0, 1, 2).
        gene_id: Gene identifier for title.
        variant_id: Variant identifier for title.
        output_path: Path to save figure (optional).
        **kwargs: Additional matplotlib kwargs.

    Returns:
        Matplotlib figure object.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for plotting")
        return None

    fig, ax = plt.subplots(figsize=kwargs.get("figsize", (8, 6)))

    # Group by genotype
    gt_rounded = np.round(genotypes).astype(int)
    groups = {0: [], 1: [], 2: []}

    for gt, expr in zip(gt_rounded, expression):
        if gt in groups and not np.isnan(expr):
            groups[gt].append(expr)

    # Create boxplot
    data = [groups[0], groups[1], groups[2]]
    labels = ["0/0 (Ref)", "0/1 (Het)", "1/1 (Alt)"]

    # Filter empty groups
    valid_data = [(d, l) for d, l in zip(data, labels) if len(d) > 0]
    if valid_data:
        data, labels = zip(*valid_data)
        bp = ax.boxplot(data, labels=labels, patch_artist=True)

        colors = ["#4CAF50", "#FFC107", "#F44336"][: len(data)]
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)

    ax.set_xlabel("Genotype")
    ax.set_ylabel("Expression Level")
    title = f"{gene_id} ~ {variant_id}" if gene_id and variant_id else "eQTL Boxplot"
    ax.set_title(title)

    # Add n per group
    for i, d in enumerate(data):
        ax.annotate(f"n={len(d)}", (i + 1, ax.get_ylim()[0]), ha="center", fontsize=8)

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved boxplot to {output_path}")

    return fig


def plot_locuszoom_eqtl(
    results: pd.DataFrame,
    gene_id: str,
    gene_start: int,
    gene_end: int,
    ld_matrix: np.ndarray | None = None,
    output_path: Path | str | None = None,
    **kwargs: Any,
) -> Any:
    """Create LocusZoom-style regional plot for eQTL.

    Args:
        results: eQTL results with variant_id, position, pvalue.
        gene_id: Gene identifier.
        gene_start: Gene start position.
        gene_end: Gene end position.
        ld_matrix: LD matrix for coloring by r² (optional).
        output_path: Path to save figure.
        **kwargs: Additional matplotlib kwargs.

    Returns:
        Matplotlib figure object.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for plotting")
        return None

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=kwargs.get("figsize", (12, 8)), height_ratios=[3, 1], sharex=True)

    # Filter to gene region
    if "position" not in results.columns:
        logger.warning("No position column in results")
        return None

    positions = results["position"] / 1e6  # Convert to Mb
    log_pval = -np.log10(results["pvalue"].clip(lower=1e-300))

    # Find lead variant
    lead_idx = log_pval.argmax()

    # Color by LD if provided
    if ld_matrix is not None and len(ld_matrix) == len(results):
        r2 = ld_matrix[lead_idx] ** 2
        colors = plt.cm.RdYlBu_r(r2)
    else:
        colors = "steelblue"

    ax1.scatter(positions, log_pval, c=colors, s=50, alpha=0.7, edgecolors="black", linewidth=0.5)
    ax1.scatter(
        positions.iloc[lead_idx],
        log_pval.iloc[lead_idx],
        c="purple",
        s=150,
        marker="D",
        edgecolors="black",
        linewidth=2,
        zorder=10,
    )

    # Significance line
    ax1.axhline(-np.log10(5e-8), color="red", linestyle="--", alpha=0.5, label="Genome-wide")
    ax1.set_ylabel("-log₁₀(P-value)")
    ax1.set_title(f"eQTL LocusZoom: {gene_id}")
    ax1.legend()

    # Gene track
    gene_mid = (gene_start + gene_end) / 2 / 1e6
    gene_width = (gene_end - gene_start) / 1e6
    ax2.barh(0, gene_width, left=gene_start / 1e6, height=0.4, color="navy", alpha=0.7)
    ax2.text(gene_mid, 0, gene_id, ha="center", va="center", fontsize=10, color="white", fontweight="bold")
    ax2.set_ylim(-0.5, 0.5)
    ax2.set_xlabel("Position (Mb)")
    ax2.set_ylabel("Genes")
    ax2.set_yticks([])

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved LocusZoom plot to {output_path}")

    return fig


def plot_eqtl_summary(
    summary_stats: dict[str, Any],
    output_path: Path | str | None = None,
    **kwargs: Any,
) -> Any:
    """Create summary statistics visualization.

    Args:
        summary_stats: Dictionary from eqtl_summary_stats().
        output_path: Path to save figure.
        **kwargs: Additional matplotlib kwargs.

    Returns:
        Matplotlib figure object.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for plotting")
        return None

    fig, axes = plt.subplots(1, 3, figsize=kwargs.get("figsize", (12, 4)))

    # Bar chart of counts
    counts = [
        summary_stats.get("n_tests", 0),
        summary_stats.get("n_eqtls", 0),
        summary_stats.get("n_egenes", 0),
    ]
    labels = ["Tests", "Significant eQTLs", "eGenes"]
    colors = ["#2196F3", "#4CAF50", "#FF9800"]

    axes[0].bar(labels, counts, color=colors)
    axes[0].set_ylabel("Count")
    axes[0].set_title("eQTL Discovery")

    # Pie chart: significant vs not
    if summary_stats.get("n_tests", 0) > 0:
        sig_pct = summary_stats.get("n_eqtls", 0) / summary_stats["n_tests"] * 100
        axes[1].pie(
            [sig_pct, 100 - sig_pct],
            labels=[f"Significant\n({sig_pct:.1f}%)", "Not sig"],
            colors=["#4CAF50", "#E0E0E0"],
            startangle=90,
        )
        axes[1].set_title("Significance Rate")

    # Effect size distribution
    effect_sizes = summary_stats.get("effect_sizes")
    if effect_sizes is not None and len(effect_sizes) > 0:
        axes[2].hist(
            np.abs(effect_sizes),
            bins=min(30, max(5, len(effect_sizes) // 5)),
            color="#2196F3",
            alpha=0.7,
            edgecolor="black",
        )
        axes[2].axvline(
            np.mean(np.abs(effect_sizes)),
            color="red",
            linestyle="--",
            label=f"Mean |β|: {np.mean(np.abs(effect_sizes)):.3f}",
        )
        axes[2].set_xlabel("|Effect Size (β)|")
        axes[2].set_ylabel("Count")
        axes[2].legend(fontsize=8)
    else:
        axes[2].text(
            0.5,
            0.5,
            f"Mean |β|: {summary_stats.get('mean_effect_size', 0):.3f}",
            ha="center",
            va="center",
            fontsize=14,
            transform=axes[2].transAxes,
        )
        axes[2].set_xticks([])
        axes[2].set_yticks([])
    axes[2].set_title("Effect Size")

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved summary plot to {output_path}")

    return fig
