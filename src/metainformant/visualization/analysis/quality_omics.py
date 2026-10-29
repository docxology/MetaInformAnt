"""Multi-omics quality control visualization functions.

Functions for visualizing QC metrics across VCF/variant data, single-cell
experiments, protein structures, and multi-omics quality assessments.
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

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    sns = None


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


def _hist_with_optional_seaborn(
    ax: Axes, values: np.ndarray, color: str | None, xlabel: str, ylabel: str, title: str
) -> None:
    """Plot a histogram on *ax*, using seaborn KDE styling when available.

    ``color`` may be None for the default palette.
    """
    if HAS_SEABORN:
        if color is None:
            sns.histplot(values, ax=ax, kde=True, alpha=0.7)
        else:
            sns.histplot(values, ax=ax, kde=True, alpha=0.7, color=color)
    else:
        if color is None:
            ax.hist(values, bins=50, alpha=0.7)
        else:
            ax.hist(values, bins=50, alpha=0.7, color=color)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)


def _require_sections(data: Dict[str, Any], label: str, recognised: Dict[str, Tuple[str, ...] | None]) -> List[str]:
    """Return the recognised sections present in *data*, checking each payload is shaped as the plot needs.

    ``recognised`` maps a section name to the keys its payload must have: a tuple means "a dict with at least these
    keys" (an empty tuple, any dict), ``None`` means any payload is accepted.

    Raises:
        ValueError: If no recognised section is present, or a payload lacks the keys its panel is drawn from. Such
            a section used to be skipped and left an empty panel.
    """
    sections = [name for name in recognised if name in data]
    if not sections:
        raise ValueError(f"{label} has none of the recognised sections: {list(recognised)}")
    for name in sections:
        required = recognised[name]
        if required is not None and not (isinstance(data[name], dict) and all(key in data[name] for key in required)):
            raise ValueError(f"{label}[{name!r}] needs the keys {list(required)}")
    return sections


#: VCF QC sections in drawing order, with the keys each payload needs.
_VCF_SECTIONS: Dict[str, Tuple[str, ...] | None] = {
    "qual_distribution": ("qualities", "counts"),
    "depth_distribution": ("depths", "counts"),
    "allele_frequencies": None,
    "variant_types": (),
    "titv_by_qual": ("qualities", "titv_ratios"),
    "summary_stats": (),
}


def plot_vcf_quality_metrics(
    vcf_qc_data: Dict[str, Any],
    *,
    expected_titv: float | None = 2.1,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (14, 10),
    **kwargs: Any,
) -> Axes:
    """Plot comprehensive VCF quality control metrics.

    One panel is drawn for each section present in ``vcf_qc_data`` (a compact grid of at most three columns, no
    empty panels), in this order: ``qual_distribution`` (``{"qualities", "counts"}``), ``depth_distribution``
    (``{"depths", "counts"}``), ``allele_frequencies`` (a list of frequencies, drawn as a histogram),
    ``variant_types`` (``{type: count}``, a pie chart), ``titv_by_qual`` (``{"qualities", "titv_ratios"}``) and
    ``summary_stats`` (``{statistic: value}``, bars).

    Args:
        vcf_qc_data: Dictionary containing VCF QC metrics.
        expected_titv: Transition/transversion ratio drawn as the reference line of the ``titv_by_qual`` panel.
            The default 2.1 is the typical genome-wide value for human data; exomes and other species differ, so
            pass the ratio that applies to your organism, or None to draw no reference line.
        ax: Optional matplotlib axes that should hold the panels: they are laid out inside its grid slot and
            ``ax`` is hidden. Creates a new figure if None.
        output_path: Optional path to save.
        figsize: Figure size of a new figure.

    Returns:
        The axes of the first panel

    Raises:
        ValueError: If no recognised section is present, or a section is not shaped as described above
    """
    validation.validate_type(vcf_qc_data, dict, "vcf_qc_data")

    sections = _require_sections(vcf_qc_data, "vcf_qc_data", _VCF_SECTIONS)
    axes = panel_grid(ax, len(sections), max_cols=3, figsize=figsize)

    for panel, name in zip(axes, sections):
        data = vcf_qc_data[name]
        if name == "qual_distribution":
            panel.bar(data["qualities"], data["counts"], alpha=0.7)
            panel.set_xlabel("Quality Score")
            panel.set_ylabel("Count")
            panel.set_title("Variant Quality Distribution")
            panel.grid(True, alpha=0.3)
        elif name == "depth_distribution":
            panel.bar(data["depths"], data["counts"], alpha=0.7, color="green")
            panel.set_xlabel("Read Depth")
            panel.set_ylabel("Count")
            panel.set_title("Read Depth Distribution")
            panel.grid(True, alpha=0.3)
        elif name == "allele_frequencies":
            panel.hist(data, bins=20, alpha=0.7, color="orange")
            panel.set_xlabel("Allele Frequency")
            panel.set_ylabel("Count")
            panel.set_title("Allele Frequency Spectrum")
            panel.grid(True, alpha=0.3)
        elif name == "variant_types":
            panel.pie(list(data.values()), labels=list(data.keys()), autopct="%1.1f%%", startangle=90)
            panel.set_title("Variant Types")
        elif name == "titv_by_qual":
            panel.plot(data["qualities"], data["titv_ratios"], "o-", alpha=0.7)
            if expected_titv is not None:
                panel.axhline(
                    y=expected_titv, color="red", linestyle="--", alpha=0.7, label=f"Expected ({expected_titv:g})"
                )
                panel.legend()
            panel.set_xlabel("Quality Score")
            panel.set_ylabel("Ti/Tv Ratio")
            panel.set_title("Transition/Transversion Ratio")
            panel.grid(True, alpha=0.3)
        else:  # summary_stats
            stat_names, stat_values = list(data.keys()), list(data.values())
            bars = panel.bar(range(len(stat_names)), stat_values, alpha=0.7, color="purple")
            panel.set_xticks(range(len(stat_names)))
            panel.set_xticklabels(stat_names, rotation=45, ha="right")
            panel.set_ylabel("Value")
            panel.set_title("Summary Statistics")
            panel.grid(True, alpha=0.3, axis="y")
            _annotate_bars(panel, bars, stat_values, "{:.0f}", offset=max(stat_values) * 0.01)

    if ax is None:
        root_figure(axes[0]).tight_layout()

    _save_figure(axes[0], output_path, f"VCF quality metrics plot saved to {output_path}")

    return axes[0]


def plot_singlecell_qc_metrics(
    qc_metrics: Dict[str, np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (14, 10),
    **kwargs: Any,
) -> Axes:
    """Plot comprehensive single-cell QC metrics.

    One panel is drawn for each of these that ``qc_metrics`` supports (a compact grid of at most three columns,
    no empty panels), in this order: histograms of ``n_counts`` (library size), ``n_genes`` and ``percent_mito``;
    library complexity (genes per UMI, from ``n_genes`` and ``n_counts``; cells with no counts are left out and
    the panel title says how many); a histogram of ``doublet_score``; and a heatmap of the correlations between
    whichever two or three of ``n_counts``, ``n_genes`` and ``percent_mito`` are present.

    Args:
        qc_metrics: Dictionary containing various QC metrics per cell.
        ax: Optional matplotlib axes that should hold the panels: they are laid out inside its grid slot and
            ``ax`` is hidden. Creates a new figure if None.
        output_path: Optional path to save.
        figsize: Figure size of a new figure.

    Returns:
        The axes of the first panel

    Raises:
        ValueError: If ``qc_metrics`` has none of the metrics above
    """
    validation.validate_type(qc_metrics, dict, "qc_metrics")

    correlated = [m for m in ("n_counts", "n_genes", "percent_mito") if m in qc_metrics]
    sections: List[str] = list(correlated)
    if "n_counts" in qc_metrics and "n_genes" in qc_metrics:
        sections.append("complexity")
    if "doublet_score" in qc_metrics:
        sections.append("doublet_score")
    if len(correlated) >= 2:
        sections.append("correlations")
    if not sections:
        raise ValueError("qc_metrics has none of the supported metrics: n_counts, n_genes, percent_mito, doublet_score")

    axes = panel_grid(ax, len(sections), max_cols=3, figsize=figsize)

    histograms = {
        "n_counts": (None, "Library Size", "Library Size Distribution"),
        "n_genes": ("green", "Number of Genes", "Genes Detected per Cell"),
        "percent_mito": ("red", "Mitochondrial Content (%)", "Mitochondrial Content"),
        "doublet_score": ("orange", "Doublet Score", "Doublet Scores"),
    }
    for panel, name in zip(axes, sections):
        if name in histograms:
            color, xlabel, title = histograms[name]
            _hist_with_optional_seaborn(panel, qc_metrics[name], color, xlabel, "Number of Cells", title)
        elif name == "complexity":
            n_counts = np.asarray(qc_metrics["n_counts"], dtype=float)
            n_genes = np.asarray(qc_metrics["n_genes"], dtype=float)
            has_counts = n_counts > 0
            complexity = n_genes[has_counts] / n_counts[has_counts]
            title = "Library Complexity"
            if not has_counts.all():
                title += f" ({int((~has_counts).sum())} cells with no counts left out)"
            _hist_with_optional_seaborn(panel, complexity, "purple", "Genes per UMI", "Number of Cells", title)
        else:  # correlations
            corr_matrix = np.corrcoef(np.column_stack([qc_metrics[m] for m in correlated]).T)
            if HAS_SEABORN:
                sns.heatmap(
                    corr_matrix,
                    annot=True,
                    fmt=".2f",
                    cmap="coolwarm",
                    xticklabels=correlated,
                    yticklabels=correlated,
                    ax=panel,
                    center=0,
                )
            else:
                im = panel.imshow(corr_matrix, cmap="coolwarm", aspect="equal", vmin=-1, vmax=1)
                panel.set_xticks(range(len(correlated)))
                panel.set_yticks(range(len(correlated)))
                panel.set_xticklabels(correlated, rotation=45, ha="right")
                panel.set_yticklabels(correlated)
                root_figure(panel).colorbar(im, ax=panel)
            panel.set_title("QC Metric Correlations")

    if ax is None:
        root_figure(axes[0]).tight_layout()

    _save_figure(axes[0], output_path, f"Single-cell QC metrics plot saved to {output_path}")

    return axes[0]


#: Protein structure QC sections in drawing order (a clash or quality payload may be a dict or one number).
_PROTEIN_SECTIONS: Dict[str, Tuple[str, ...] | None] = {
    "b_factors": None,
    "ramachandran_stats": (),
    "clash_score": None,
    "overall_quality": None,
}


def plot_protein_structure_quality(
    structure_quality: Dict[str, Any],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs: Any,
) -> Axes:
    """Plot protein structure quality metrics.

    One panel is drawn for each section present in ``structure_quality`` (a compact grid of at most two columns,
    no empty panels), in this order: ``b_factors`` (a list, drawn as a histogram), ``ramachandran_stats``
    (``{region: percentage}``, bars), ``clash_score`` (``{clash type: count}`` as bars, or one number as text)
    and ``overall_quality`` (``{metric: score}`` as bars, or one number as text).

    Args:
        structure_quality: Dictionary containing structure quality metrics.
        ax: Optional matplotlib axes that should hold the panels: they are laid out inside its grid slot and
            ``ax`` is hidden. Creates a new figure if None.
        output_path: Optional path to save.
        figsize: Figure size of a new figure.

    Returns:
        The axes of the first panel

    Raises:
        ValueError: If no recognised section is present
    """
    validation.validate_type(structure_quality, dict, "structure_quality")

    sections = _require_sections(structure_quality, "structure_quality", _PROTEIN_SECTIONS)
    axes = panel_grid(ax, len(sections), max_cols=2, figsize=figsize)

    for panel, name in zip(axes, sections):
        data = structure_quality[name]
        if name == "b_factors":
            panel.hist(data, bins=50, alpha=0.7, color="blue")
            panel.set_xlabel("B-factor")
            panel.set_ylabel("Count")
            panel.set_title("B-factor Distribution")
            panel.grid(True, alpha=0.3)
        elif name == "ramachandran_stats":
            panel.bar(list(data.keys()), list(data.values()), alpha=0.7, color="green")
            panel.set_ylabel("Percentage")
            panel.set_title("Ramachandran Statistics")
            panel.tick_params(axis="x", rotation=45)
            panel.grid(True, alpha=0.3, axis="y")
        elif name == "clash_score":
            if isinstance(data, dict):
                panel.bar(list(data.keys()), list(data.values()), alpha=0.7, color="red")
                panel.set_ylabel("Number of Clashes")
                panel.set_title("Steric Clashes")
                panel.tick_params(axis="x", rotation=45)
            else:
                panel.text(0.5, 0.5, f"Clash Score: {data:.2f}", ha="center", va="center", transform=panel.transAxes)
                panel.set_title("Clash Score")
            panel.grid(True, alpha=0.3, axis="y")
        else:  # overall_quality
            if isinstance(data, dict):
                bars = panel.bar(list(data.keys()), list(data.values()), alpha=0.7, color="purple")
                panel.set_ylabel("Score")
                panel.set_title("Quality Metrics")
                panel.tick_params(axis="x", rotation=45)
                _annotate_bars(panel, bars, list(data.values()), "{:.2f}", offset=0.01)
            else:
                panel.text(
                    0.5, 0.5, f"Overall Quality: {data:.2f}", ha="center", va="center", transform=panel.transAxes
                )
                panel.set_title("Overall Quality Score")
            panel.grid(True, alpha=0.3, axis="y")

    if ax is None:
        root_figure(axes[0]).tight_layout()

    _save_figure(axes[0], output_path, f"Protein structure quality plot saved to {output_path}")

    return axes[0]


def _report_quality_score(report: Dict[str, Any]) -> float | None:
    """Quality score of one omics layer's report, or None when the report has nothing numeric to score.

    The score is the first of ``overall_quality``, ``quality_score`` and ``mean_quality`` that is present, even when
    it is 0.0 (the worst quality, not a missing value); failing those, the mean of the report's numeric values.
    """
    for key in ("overall_quality", "quality_score", "mean_quality"):
        if report.get(key) is not None:
            return float(report[key])
    numeric = [
        float(value)
        for value in report.values()
        if isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool)
    ]
    return float(np.mean(numeric)) if numeric else None


def plot_multiomics_quality_overview(
    quality_reports: Dict[str, Dict[str, Any]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (14, 8),
    **kwargs: Any,
) -> Axes:
    """Plot quality overview across multiple omics layers as a radar chart.

    Each layer's score is the first of ``overall_quality``, ``quality_score`` and ``mean_quality`` found in its
    report (a value of 0.0 counts), otherwise the mean of the report's numeric values. A layer whose report has
    nothing numeric is drawn at 0 and its label reads "n/a" instead of a score.

    Args:
        quality_reports: Dictionary mapping omics types to quality reports.
        ax: Optional polar matplotlib axes. Creates a new polar axes if None.
        output_path: Optional path to save.
        figsize: Figure size of a new figure.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If ``quality_reports`` is empty, or ``ax`` is not a polar axes
    """
    validation.validate_type(quality_reports, dict, "quality_reports")
    if not quality_reports:
        raise ValueError("quality_reports cannot be empty")

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, polar=True)
    elif ax.name != "polar":
        raise ValueError("A radar chart needs a polar axes (create it with projection='polar')")

    omics_types = list(quality_reports.keys())
    scores = [_report_quality_score(quality_reports[omics_type]) for omics_type in omics_types]
    quality_scores = [0.0 if score is None else score for score in scores]

    angles = np.linspace(0, 2 * np.pi, len(omics_types), endpoint=False).tolist()
    closed_angles = angles + angles[:1]
    closed_scores = quality_scores + quality_scores[:1]

    ax.plot(closed_angles, closed_scores, "o-", linewidth=2, label="Quality Score", alpha=0.8)
    ax.fill(closed_angles, closed_scores, alpha=0.25)
    ax.set_xticks(angles)
    ax.set_xticklabels(omics_types)
    top = max(quality_scores) * 1.1 or 1.0
    ax.set_ylim(0, top)
    ax.set_title("Multi-Omics Quality Overview")
    ax.grid(True, alpha=0.3)

    for angle, quality_score, score in zip(angles, quality_scores, scores):
        ax.text(
            angle,
            quality_score + top * 0.05,
            "n/a" if score is None else f"{score:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    _save_figure(ax, output_path, f"Multi-omics quality overview saved to {output_path}")

    return ax
