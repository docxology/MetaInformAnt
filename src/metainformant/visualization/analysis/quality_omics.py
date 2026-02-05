"""Multi-omics quality control visualization functions.

Functions for visualizing QC metrics across VCF/variant data, single-cell
experiments, protein structures, and multi-omics quality assessments.
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

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    sns = None


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
        vcf_qc_data: Dictionary containing VCF QC metrics.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.
        figsize: Figure size.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(vcf_qc_data, dict, "vcf_qc_data")

    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes = axes.flatten()

    if "qual_distribution" in vcf_qc_data:
        qual_data = vcf_qc_data["qual_distribution"]
        if "qualities" in qual_data and "counts" in qual_data:
            axes[0].bar(qual_data["qualities"], qual_data["counts"], alpha=0.7)
            axes[0].set_xlabel("Quality Score")
            axes[0].set_ylabel("Count")
            axes[0].set_title("Variant Quality Distribution")
            axes[0].grid(True, alpha=0.3)

    if "depth_distribution" in vcf_qc_data:
        depth_data = vcf_qc_data["depth_distribution"]
        if "depths" in depth_data and "counts" in depth_data:
            axes[1].bar(depth_data["depths"], depth_data["counts"], alpha=0.7, color="green")
            axes[1].set_xlabel("Read Depth")
            axes[1].set_ylabel("Count")
            axes[1].set_title("Read Depth Distribution")
            axes[1].grid(True, alpha=0.3)

    if "allele_frequencies" in vcf_qc_data:
        af_data = vcf_qc_data["allele_frequencies"]
        axes[2].hist(af_data, bins=20, alpha=0.7, color="orange")
        axes[2].set_xlabel("Allele Frequency")
        axes[2].set_ylabel("Count")
        axes[2].set_title("Allele Frequency Spectrum")
        axes[2].grid(True, alpha=0.3)

    if "variant_types" in vcf_qc_data:
        types_data = vcf_qc_data["variant_types"]
        types = list(types_data.keys())
        counts = list(types_data.values())
        axes[3].pie(counts, labels=types, autopct="%1.1f%%", startangle=90)
        axes[3].set_title("Variant Types")

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
        qc_metrics: Dictionary containing various QC metrics per cell.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.
        figsize: Figure size.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(qc_metrics, dict, "qc_metrics")

    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes = axes.flatten()
    plot_idx = 0

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

    if len(qc_metrics) >= 2:
        metrics_to_plot = ["n_counts", "n_genes", "percent_mito"]
        available_metrics = [m for m in metrics_to_plot if m in qc_metrics]
        if len(available_metrics) >= 2:
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
        structure_quality: Dictionary containing structure quality metrics.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.
        figsize: Figure size.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(structure_quality, dict, "structure_quality")

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()

    if "b_factors" in structure_quality:
        b_factors = structure_quality["b_factors"]
        axes[0].hist(b_factors, bins=50, alpha=0.7, color="blue")
        axes[0].set_xlabel("B-factor")
        axes[0].set_ylabel("Count")
        axes[0].set_title("B-factor Distribution")
        axes[0].grid(True, alpha=0.3)

    if "ramachandran_stats" in structure_quality:
        rama_stats = structure_quality["ramachandran_stats"]
        categories = list(rama_stats.keys())
        values = list(rama_stats.values())
        axes[1].bar(categories, values, alpha=0.7, color="green")
        axes[1].set_ylabel("Percentage")
        axes[1].set_title("Ramachandran Plot")
        axes[1].tick_params(axis="x", rotation=45)
        axes[1].grid(True, alpha=0.3, axis="y")

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

    if "overall_quality" in structure_quality:
        quality_data = structure_quality["overall_quality"]
        if isinstance(quality_data, dict):
            metrics = list(quality_data.keys())
            values = list(quality_data.values())
            bars = axes[3].bar(metrics, values, alpha=0.7, color="purple")
            axes[3].set_ylabel("Score")
            axes[3].set_title("Quality Metrics")
            axes[3].tick_params(axis="x", rotation=45)
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
        quality_reports: Dictionary mapping omics types to quality reports.
        ax: Optional matplotlib axes.
        output_path: Optional path to save.
        figsize: Figure size.

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(quality_reports, dict, "quality_reports")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    omics_types = list(quality_reports.keys())
    quality_scores = []
    for omics_type in omics_types:
        report = quality_reports[omics_type]
        score = (
            report.get("overall_quality")
            or report.get("quality_score")
            or report.get("mean_quality")
            or np.mean(list(report.values()))
            if report
            else 0
        )
        quality_scores.append(score)

    angles = np.linspace(0, 2 * np.pi, len(omics_types), endpoint=False).tolist()
    angles += angles[:1]
    quality_scores += quality_scores[:1]

    ax.plot(angles, quality_scores, "o-", linewidth=2, label="Quality Score", alpha=0.8)
    ax.fill(angles, quality_scores, alpha=0.25)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(omics_types)
    ax.set_ylim(0, max(quality_scores) * 1.1)
    ax.set_title("Multi-Omics Quality Overview")
    ax.grid(True, alpha=0.3)

    for angle, score, omics in zip(angles[:-1], quality_scores[:-1], omics_types):
        ax.text(angle, score + max(quality_scores) * 0.05, f"{score:.2f}", ha="center", va="bottom", fontsize=8)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Multi-omics quality overview saved to {output_path}")

    return ax
