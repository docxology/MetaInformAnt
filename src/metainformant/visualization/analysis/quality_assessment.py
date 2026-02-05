"""Quality assessment visualization for coverage, errors, and batch effects.

Functions for visualizing coverage uniformity, error profiles,
batch effect diagnostics, and data integrity metrics.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from metainformant.core import logging, paths, validation

logger = logging.get_logger(__name__)


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
        batch_qc_data: Dictionary containing batch effect analysis results.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.
        figsize: Figure size.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(batch_qc_data, dict, "batch_qc_data")

    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes = axes.flatten()
    plot_idx = 0

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

    if "batch_variance" in batch_qc_data:
        variance_data = batch_qc_data["batch_variance"]
        if isinstance(variance_data, dict):
            components = list(variance_data.keys())
            variances = list(variance_data.values())
            axes[plot_idx].pie(variances, labels=components, autopct="%1.1f%%", startangle=90)
            axes[plot_idx].set_title("Variance Explained by Batch")
        plot_idx += 1

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
    bars = ax.bar(x_positions, values, color=colors, alpha=0.8, **kwargs)

    ax.set_xlabel("Integrity Metric")
    ax.set_ylabel("Value")
    ax.set_title("Data Integrity Assessment")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(metrics, rotation=45, ha="right")
    ax.grid(True, alpha=0.3, axis="y")

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
