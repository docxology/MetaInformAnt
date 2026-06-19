"""GWAS visualization utilities.

This module provides functions for creating GWAS visualization plots,
including Manhattan plots, Q-Q plots, and regional association plots.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Import matplotlib with graceful fallback
try:
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    logger.warning("matplotlib not available, visualization functions will return None")

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    logger.warning("numpy not available, some visualizations may not work")

try:
    from scipy import stats as _scipy_stats

    HAS_SCIPY = True
except ImportError:  # pragma: no cover - used only in lean environments
    _scipy_stats = None
    HAS_SCIPY = False

EXPECTED_MEDIAN_CHI2_1DF = 0.454936423119572


def _coerce_float(value: Any, default: Optional[float] = None) -> Optional[float]:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    if not math.isfinite(result):
        return default
    return result


def _chrom_sort_key(chrom: Any) -> tuple:
    """Sort chromosome/contig names naturally across numeric and accession IDs."""
    text = str(chrom).strip()
    lowered = text.lower()
    if lowered.startswith("chr"):
        lowered = lowered[3:]

    if lowered.isdigit():
        return (0, int(lowered), text)

    special = {"x": 23, "y": 24, "m": 25, "mt": 25}
    if lowered in special:
        return (1, special[lowered], text)

    numbers = [int(part) for part in re.findall(r"\d+", lowered)]
    if numbers:
        return (2, numbers, lowered, text)
    return (3, lowered, text)


def _neg_log10_p(p_value: Any, cap: float = 300.0) -> float:
    p = _coerce_float(p_value)
    if p is None:
        return 0.0
    if p <= 0:
        return cap
    if p > 1:
        return 0.0
    return min(-math.log10(max(p, 10 ** (-cap))), cap)


def _normalise_gwas_results(results: Union[List[Dict[str, Any]], Dict[str, Any]]) -> List[Dict[str, Any]]:
    if isinstance(results, dict):
        source_rows: Iterable[Any] = [results]
    elif isinstance(results, list):
        source_rows = results
    else:
        source_rows = []

    rows: List[Dict[str, Any]] = []
    for original_index, result in enumerate(source_rows):
        if not isinstance(result, dict):
            continue
        chrom = str(result.get("chrom", result.get("chromosome", "1")))
        pos = _coerce_float(result.get("pos", result.get("position", 0)), 0.0)
        p_val = result.get("p_value", result.get("pval", result.get("pvalue", 1.0)))
        rows.append(
            {
                "chrom": chrom,
                "pos": pos if pos is not None else 0.0,
                "p_value": p_val,
                "neg_log_p": _neg_log10_p(p_val),
                "original_index": original_index,
                "source": result,
            }
        )
    rows.sort(key=lambda row: (_chrom_sort_key(row["chrom"]), row["pos"], row["original_index"]))
    return rows


def _compute_genome_axis(rows: List[Dict[str, Any]]) -> tuple[Dict[str, float], Dict[str, float], float]:
    """Assign cumulative x positions using observed contig extents."""
    if not rows:
        return {}, {}, 0.0

    chrom_order: List[str] = []
    chrom_bounds: Dict[str, tuple[float, float]] = {}
    for row in rows:
        chrom = row["chrom"]
        pos = float(row["pos"])
        if chrom not in chrom_bounds:
            chrom_order.append(chrom)
            chrom_bounds[chrom] = (pos, pos)
        else:
            lo, hi = chrom_bounds[chrom]
            chrom_bounds[chrom] = (min(lo, pos), max(hi, pos))

    largest_extent = max(max(hi - lo, hi, 1.0) for lo, hi in chrom_bounds.values())
    gap = max(1_000.0, min(largest_extent * 0.01, 5_000_000.0))

    offsets: Dict[str, float] = {}
    centers: Dict[str, float] = {}
    current = 0.0
    for chrom in chrom_order:
        lo, hi = chrom_bounds[chrom]
        offsets[chrom] = current
        centers[chrom] = current + (lo + hi) / 2.0
        current += max(hi, lo + 1.0) + gap

    for row in rows:
        row["global_pos"] = offsets[row["chrom"]] + float(row["pos"])

    return offsets, centers, current


def _lambda_gc_from_pvalues(p_values: Sequence[float]) -> Optional[float]:
    valid = [float(p) for p in p_values if math.isfinite(float(p)) and 0 < float(p) <= 1]
    if not valid:
        return None
    if HAS_SCIPY and _scipy_stats is not None:
        chi2 = np.asarray([float(_scipy_stats.chi2.isf(min(max(p, 1e-300), 1.0), 1)) for p in valid])
    else:
        chi2 = np.asarray([-2.0 * math.log(max(p, 1e-300)) for p in valid])
    chi2 = chi2[np.isfinite(chi2)]
    if chi2.size == 0:
        return None
    return float(np.median(chi2) / EXPECTED_MEDIAN_CHI2_1DF)


def _qq_confidence_band(n: int, alpha: float = 0.05) -> tuple[Any, Any, Any]:
    ranks = np.arange(1, n + 1)
    expected = (ranks - 0.5) / n
    if HAS_SCIPY and _scipy_stats is not None:
        lower_p = _scipy_stats.beta.ppf(alpha / 2.0, ranks, n - ranks + 1)
        upper_p = _scipy_stats.beta.ppf(1.0 - alpha / 2.0, ranks, n - ranks + 1)
    else:
        mean = ranks / (n + 1.0)
        sd = np.sqrt((ranks * (n - ranks + 1.0)) / (((n + 1.0) ** 2) * (n + 2.0)))
        lower_p = mean - 1.96 * sd
        upper_p = mean + 1.96 * sd
    lower_p = np.clip(lower_p, 1e-300, 1.0)
    upper_p = np.clip(upper_p, 1e-300, 1.0)
    return -np.log10(expected), -np.log10(upper_p), -np.log10(lower_p)


def manhattan_plot(
    results: Union[List[Dict[str, Any]], Dict[str, Any]],
    output_path: Optional[Union[str, Path]] = None,
    significance_threshold: float = 5e-8,
    gene_annotations: Optional[List[Dict[str, Any]]] = None,
    highlight_regions: Optional[List[Dict[str, Any]]] = None,
    suggestive_threshold: Optional[float] = 1e-5,
    label_top_n: int = 0,
    title: str = "Manhattan Plot",
    subtitle: Optional[str] = None,
    point_size: Optional[float] = None,
) -> Any:
    """Create a Manhattan plot from GWAS results.

    Args:
        results: GWAS results dictionary or list of result dictionaries
        output_path: Path to save the plot (optional)
        significance_threshold: P-value threshold for significance line
        gene_annotations: Optional list of dicts with keys: variant_index or
            (chrom, pos), and gene_name. When provided, annotate top hits with
            gene names using arrows.
        highlight_regions: Optional list of dicts with keys: chrom, start, end,
            color (optional), label (optional). When provided, draw colored
            rectangular highlights over those genomic regions.
        suggestive_threshold: Optional p-value threshold for a secondary dashed
            significance line. Defaults to 1e-5. Set to None to disable.
        label_top_n: Auto-label the top N hits with their variant_id or
            chr:pos. Defaults to 0 (disabled).
        title: Plot title.
        subtitle: Optional subtitle rendered under the title.
        point_size: Optional scatter point size. If omitted, a density-aware
            default is used.

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create Manhattan plot")
        return None

    logger.info("Creating Manhattan plot")

    rows = _normalise_gwas_results(results)
    chrom_offsets, chrom_centers, axis_max = _compute_genome_axis(rows)
    p_values = [float(row["neg_log_p"]) for row in rows]
    positions = [float(row.get("global_pos", 0.0)) for row in rows]

    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))

    if not rows:
        ax.text(
            0.5,
            0.5,
            "No GWAS results to plot",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=12,
            color="#555555",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    else:
        palette = [
            "#2f5d8c",
            "#c9792f",
            "#4f7d63",
            "#7a5ea7",
            "#9b4f5f",
            "#6b6f3f",
        ]
        chrom_order = list(chrom_offsets.keys())
        color_by_chrom = {chrom: palette[i % len(palette)] for i, chrom in enumerate(chrom_order)}
        if point_size is None:
            point_size = 9.0 if len(rows) <= 2_000 else 5.0 if len(rows) <= 20_000 else 2.0

        for chrom in chrom_order:
            chrom_rows = [row for row in rows if row["chrom"] == chrom]
            ax.scatter(
                [row["global_pos"] for row in chrom_rows],
                [row["neg_log_p"] for row in chrom_rows],
                c=color_by_chrom[chrom],
                s=point_size,
                alpha=0.72,
                linewidths=0,
            )

    # Add significance threshold line
    if significance_threshold > 0:
        threshold_line = -math.log10(significance_threshold)
        ax.axhline(
            y=threshold_line,
            color="#b23a48",
            linestyle="--",
            alpha=0.7,
            label=f"Significance threshold ({significance_threshold})",
        )

    # Add suggestive threshold line
    if suggestive_threshold is not None and suggestive_threshold > 0:
        suggestive_line = -math.log10(suggestive_threshold)
        ax.axhline(
            y=suggestive_line,
            color="#4b6c91",
            linestyle=":",
            alpha=0.75,
            label=f"Suggestive threshold ({suggestive_threshold})",
        )

    # Draw highlight regions
    if highlight_regions:
        for region in highlight_regions:
            region_chrom = str(region.get("chrom", ""))
            region_start = region.get("start", 0)
            region_end = region.get("end", 0)
            region_color = region.get("color", "yellow")
            region_label = region.get("label", None)

            if region_chrom in chrom_offsets:
                start_pos = _coerce_float(region_start, 0.0) or 0.0
                end_pos = _coerce_float(region_end, start_pos) or start_pos
                start_x = chrom_offsets[region_chrom] + start_pos
                end_x = chrom_offsets[region_chrom] + end_pos
                ax.axvspan(start_x, end_x, alpha=0.2, color=region_color, label=region_label)

    # Annotate genes
    if gene_annotations:
        for annotation in gene_annotations:
            gene_name = annotation.get("gene_name", "")
            ann_x = None
            ann_y = None

            # Look up by variant_index
            if "variant_index" in annotation:
                idx = annotation["variant_index"]
                matching = [row for row in rows if row["original_index"] == idx]
                if not matching and isinstance(idx, int) and 0 <= idx < len(rows):
                    matching = [rows[idx]]
                if matching:
                    ann_x = matching[0]["global_pos"]
                    ann_y = matching[0]["neg_log_p"]
            # Look up by chrom + pos
            elif "chrom" in annotation and "pos" in annotation:
                ann_chrom = str(annotation["chrom"])
                ann_pos = _coerce_float(annotation["pos"], 0.0)
                candidates = [row for row in rows if row["chrom"] == ann_chrom]
                if candidates and ann_pos is not None:
                    best = min(candidates, key=lambda row: abs(float(row["pos"]) - ann_pos))
                    ann_x = best["global_pos"]
                    ann_y = best["neg_log_p"]

            if ann_x is not None and ann_y is not None:
                offset = max(max(p_values) * 0.08, 0.4) if p_values else 1.0
                ax.annotate(
                    gene_name,
                    xy=(ann_x, ann_y),
                    xytext=(ann_x, ann_y + offset),
                    arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
                    fontsize=8,
                    ha="center",
                )

    # Label top N hits
    if label_top_n > 0 and rows:
        top_rows = sorted(rows, key=lambda row: row["neg_log_p"], reverse=True)[:label_top_n]
        offset = max(max(p_values) * 0.06, 0.35) if p_values else 1.0

        for rank, row in enumerate(top_rows):
            px = row["global_pos"]
            py = row["neg_log_p"]
            # Use variant_id if available, otherwise chr:pos
            result_entry = row["source"]
            variant_label = result_entry.get(
                "variant_id",
                f"{result_entry.get('chrom', result_entry.get('chromosome', '?'))}:"
                f"{result_entry.get('pos', result_entry.get('position', '?'))}",
            )
            # Stagger offsets slightly to reduce overlap
            y_offset = offset + rank * (offset * 0.4)
            ax.annotate(
                variant_label,
                xy=(px, py),
                xytext=(px, py + y_offset),
                arrowprops=dict(arrowstyle="->", color="gray", lw=0.6),
                fontsize=7,
                ha="center",
                color="darkred",
            )

    # Add chromosome labels
    if chrom_centers:
        tick_rotation = 0 if len(chrom_centers) <= 12 else 45
        tick_alignment = "center" if len(chrom_centers) <= 12 else "right"
        ax.set_xticks(list(chrom_centers.values()))
        ax.set_xticklabels(list(chrom_centers.keys()), rotation=tick_rotation, ha=tick_alignment)
        ax.set_xlim(0, axis_max if axis_max > 0 else max(positions) * 1.02)

    # Labels and title
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log₁₀(p-value)")
    if subtitle:
        ax.set_title(f"{title}\n{subtitle}", fontsize=13)
    else:
        ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.24)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Collect any auto-labeled artists (threshold lines, highlight regions)
    legend_elements: List[Any] = []
    auto_handles, auto_labels = ax.get_legend_handles_labels()
    for handle, lbl in zip(auto_handles, auto_labels):
        if lbl and lbl not in [getattr(e, "get_label", lambda: "")() for e in legend_elements]:
            legend_elements.append(handle)
    if legend_elements:
        ax.legend(handles=legend_elements, loc="upper right", frameon=False)

    plt.tight_layout()

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved Manhattan plot to {output_path}")

    return fig


def qq_plot(
    p_values: Union[List[float], List[int]],
    output_path: Optional[Union[str, Path]] = None,
    confidence_band: bool = True,
    show_lambda: bool = True,
    title: Optional[str] = None,
) -> Any:
    """Create a Q-Q plot from p-values.

    Args:
        p_values: List of p-values
        output_path: Path to save the plot (optional)
        confidence_band: Draw the 95% null envelope for order statistics.
        show_lambda: Annotate lambda GC on the plot.
        title: Optional custom plot title.

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create Q-Q plot")
        return None

    logger.info("Creating Q-Q plot")

    raw_values = list(p_values)
    valid_values = []
    for p_value in raw_values:
        p = _coerce_float(p_value)
        if p is not None and 0 < p <= 1:
            valid_values.append(p)

    if len(valid_values) == 0:
        logger.warning("No valid p-values for Q-Q plot")
        if raw_values:
            raise ValueError("No valid p-values in the interval (0, 1]")
        return None

    p_vals = np.array(valid_values, dtype=float)

    # Sort p-values
    p_vals = np.sort(p_vals)

    # Expected p-values under null hypothesis
    n = len(p_vals)
    expected = (np.arange(1, n + 1) - 0.5) / n
    expected_log = -np.log10(expected)

    # Observed -log10 p-values
    observed_log = -np.log10(p_vals)

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 6))

    if confidence_band and n >= 2:
        band_x, band_low, band_high = _qq_confidence_band(n)
        ax.fill_between(
            band_x,
            band_low,
            band_high,
            color="#8fb3d9",
            alpha=0.25,
            linewidth=0,
            label="95% null envelope",
        )

    # Plot Q-Q points
    ax.scatter(expected_log, observed_log, s=8 if n <= 2_000 else 3, alpha=0.62, color="#2f5d8c", linewidths=0)

    # Add diagonal line (expected under null)
    max_val = max(np.max(expected_log), np.max(observed_log))
    ax.plot([0, max_val], [0, max_val], color="#b23a48", linestyle="--", alpha=0.75, label="Expected (null)")

    lambda_gc = _lambda_gc_from_pvalues(p_vals.tolist()) if show_lambda else None
    if lambda_gc is not None:
        ax.text(
            0.03,
            0.97,
            f"lambda GC = {lambda_gc:.3f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#d0d0d0", alpha=0.9),
        )

    # Labels and title
    ax.set_xlabel("Expected -log₁₀(p-value)")
    ax.set_ylabel("Observed -log₁₀(p-value)")
    ax.set_title(title or f"Q-Q Plot (n={n} variants)")
    ax.set_xlim(left=0, right=max_val * 1.02 if max_val > 0 else 1)
    ax.set_ylim(bottom=0, top=max_val * 1.08 if max_val > 0 else 1)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.24)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False)

    plt.tight_layout()

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved Q-Q plot to {output_path}")

    return fig


def regional_plot(
    results: List[Dict[str, Any]], chrom: str, start: int, end: int, output_path: Optional[Union[str, Path]] = None
) -> Any:
    """Create a regional association plot.

    Args:
        results: GWAS results for the region
        chrom: Chromosome
        start: Start position
        end: End position
        output_path: Path to save the plot (optional)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create regional plot")
        return None

    logger.info(f"Creating regional plot for {chrom}:{start}-{end}")

    # Filter results to the specified region
    region_results = []
    for result in results:
        r_chrom = str(result.get("chrom", result.get("chromosome", "")))
        r_pos = result.get("pos", result.get("position", 0))
        if r_chrom == str(chrom) and start <= r_pos <= end:
            region_results.append(result)

    if not region_results:
        logger.warning(f"No results found in region {chrom}:{start}-{end}")
        return None

    # Extract positions and p-values
    positions = []
    p_values = []
    for result in region_results:
        pos = result.get("pos", result.get("position", 0))
        p_val = result.get("p_value", result.get("pval", 1.0))
        positions.append(pos)
        if p_val > 0:
            p_values.append(-math.log10(p_val))
        else:
            p_values.append(50)  # Cap very small p-values

    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot points
    ax.scatter(positions, p_values, s=20, alpha=0.7, color="blue")

    # Add significance threshold line
    threshold_line = -math.log10(5e-8)
    ax.axhline(y=threshold_line, color="red", linestyle="--", alpha=0.7, label="Genome-wide significance")

    # Labels and title
    ax.set_xlabel(f"Position on chromosome {chrom}")
    ax.set_ylabel("-log₁₀(p-value)")
    ax.set_title(f"Regional Association Plot: {chrom}:{start:,}-{end:,}")
    ax.set_xlim(start, end)
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved regional plot to {output_path}")

    return fig


def pca_plot(
    pca_result: tuple, output_path: Optional[Union[str, Path]] = None, explained_var: Optional[List[float]] = None
) -> Any:
    """Create PCA scatter plot.

    Args:
        pca_result: PCA results tuple (components, variance, loadings)
        output_path: Path to save the plot (optional)
        explained_var: Explained variance ratios

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create PCA plot")
        return None

    logger.info("Creating PCA plot")

    try:
        components, variance, loadings = pca_result

        if len(components) < 2:
            logger.warning("Need at least 2 PCA components for plotting")
            return None

        # Create 2D scatter plot of first two components
        fig, ax = plt.subplots(figsize=(10, 8))

        # Plot points
        scatter = ax.scatter(components[0], components[1], s=2, alpha=0.6, color="blue")

        # Labels and title
        pc1_var = explained_var[0] * 100 if explained_var and len(explained_var) > 0 else 0
        pc2_var = explained_var[1] * 100 if explained_var and len(explained_var) > 1 else 0

        ax.set_xlabel(f"PC1 ({pc1_var:.1f}% variance)")
        ax.set_ylabel(f"PC2 ({pc2_var:.1f}% variance)")
        ax.set_title("PCA Scatter Plot")
        ax.grid(True, alpha=0.3)

        # Add explained variance text if available
        if explained_var and len(explained_var) >= 2:
            var_text = f"PC1: {pc1_var:.1f}%, PC2: {pc2_var:.1f}%"
            ax.text(
                0.02,
                0.98,
                var_text,
                transform=ax.transAxes,
                verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
            )

        plt.tight_layout()

        # Save if output path provided
        if output_path:
            output_path = Path(output_path)
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            logger.info(f"Saved PCA plot to {output_path}")

        return fig

    except (ValueError, IndexError, TypeError) as e:
        logger.error(f"Error creating PCA plot: {e}")
        return None


def kinship_heatmap(
    kinship_matrix: Union[np.ndarray, List[List[float]]], output_path: Optional[Union[str, Path]] = None
) -> Any:
    """Create kinship matrix heatmap.

    Args:
        kinship_matrix: Kinship matrix
        output_path: Path to save the plot (optional)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create kinship heatmap")
        return None

    logger.info("Creating kinship heatmap")

    try:
        # Convert to numpy array if needed
        if isinstance(kinship_matrix, list):
            kinship_matrix = np.array(kinship_matrix)

        # Create heatmap
        fig, ax = plt.subplots(figsize=(10, 8))

        # Plot heatmap
        im = ax.imshow(kinship_matrix, cmap="viridis", aspect="equal")

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Kinship coefficient")

        # Labels and title
        ax.set_title("Kinship Matrix Heatmap")
        ax.set_xlabel("Sample")
        ax.set_ylabel("Sample")

        plt.tight_layout()

        # Save if output path provided
        if output_path:
            output_path = Path(output_path)
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            logger.info(f"Saved kinship heatmap to {output_path}")

        return fig

    except Exception as e:
        logger.error(f"Error creating kinship heatmap: {e}")
        return None


def effect_size_plot(results: List[Dict[str, Any]], output_path: Optional[Union[str, Path]] = None) -> Any:
    """Create effect size distribution plot.

    Args:
        results: GWAS results with 'beta' field for effect sizes
        output_path: Path to save the plot (optional)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create effect size plot")
        return None

    logger.info("Creating effect size plot")

    # Extract effect sizes (beta values)
    effect_sizes = []
    for result in results:
        beta = result.get("beta", result.get("effect_size", None))
        if beta is not None:
            effect_sizes.append(float(beta))

    if not effect_sizes:
        logger.warning("No effect sizes found in results")
        return None

    effect_sizes = np.array(effect_sizes)

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Histogram of effect sizes
    ax1.hist(effect_sizes, bins=50, edgecolor="black", alpha=0.7, color="steelblue")
    ax1.axvline(x=0, color="red", linestyle="--", alpha=0.7, label="Null effect")
    ax1.axvline(
        x=np.mean(effect_sizes), color="green", linestyle="-", alpha=0.7, label=f"Mean: {np.mean(effect_sizes):.4f}"
    )
    ax1.set_xlabel("Effect Size (Beta)")
    ax1.set_ylabel("Frequency")
    ax1.set_title("Effect Size Distribution")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Box plot
    ax2.boxplot(effect_sizes, vert=True)
    ax2.set_ylabel("Effect Size (Beta)")
    ax2.set_title("Effect Size Box Plot")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved effect size plot to {output_path}")

    return fig


def generate_all_plots(
    association_results: Union[str, Path],
    output_dir: Union[str, Path],
    pca_file: Optional[Union[str, Path]] = None,
    kinship_file: Optional[Union[str, Path]] = None,
    vcf_file: Optional[Union[str, Path]] = None,
    significance_threshold: float = 5e-8,
) -> Dict[str, Any]:
    """Generate all GWAS visualization plots.

    Args:
        association_results: Path to association results or results data
        output_dir: Output directory for plots
        pca_file: Path to PCA results file
        kinship_file: Path to kinship matrix file
        vcf_file: Path to VCF file
        significance_threshold: Significance threshold

    Returns:
        Dictionary with plot file paths and metadata
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Generating all GWAS plots in {output_dir}")

    plots_generated: Dict[str, Any] = {}

    # Load association results from file or use directly
    results_data: List[Dict[str, Any]] = []
    association_path = Path(association_results) if isinstance(association_results, (str, Path)) else None

    if association_path and association_path.exists():
        import csv
        import json

        suffix = association_path.suffix.lower()
        if suffix == ".json":
            with open(association_path) as fh:
                loaded = json.load(fh)
                results_data = loaded if isinstance(loaded, list) else loaded.get("results", [])
        elif suffix in (".tsv", ".txt"):
            with open(association_path, newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    entry: Dict[str, Any] = {}
                    for k, v in row.items():
                        try:
                            entry[k] = float(v)
                        except (ValueError, TypeError):
                            entry[k] = v
                    results_data.append(entry)
        elif suffix == ".csv":
            with open(association_path, newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    entry = {}
                    for k, v in row.items():
                        try:
                            entry[k] = float(v)
                        except (ValueError, TypeError):
                            entry[k] = v
                    results_data.append(entry)

    # Manhattan plot
    if results_data:
        try:
            manhattan_path = output_dir / "manhattan_plot.png"
            fig = manhattan_plot(
                results_data, output_path=manhattan_path, significance_threshold=significance_threshold
            )
            if fig is not None:
                plots_generated["manhattan"] = str(manhattan_path)
                plt.close(fig)
        except Exception as e:
            logger.warning(f"Manhattan plot failed: {e}")

        # Q-Q plot
        try:
            p_vals = [
                r.get("p_value", r.get("pval", r.get("pvalue")))
                for r in results_data
                if r.get("p_value", r.get("pval", r.get("pvalue"))) is not None
            ]
            if p_vals:
                qq_path = output_dir / "qq_plot.png"
                fig = qq_plot(p_vals, output_path=qq_path)
                if fig is not None:
                    plots_generated["qq"] = str(qq_path)
                    plt.close(fig)
        except Exception as e:
            logger.warning(f"Q-Q plot failed: {e}")

    # PCA plot
    if pca_file:
        try:
            import json as _json

            pca_path_obj = Path(pca_file)
            if pca_path_obj.exists():
                with open(pca_path_obj) as fh:
                    pca_data = _json.load(fh)
                components = pca_data.get("components", [])
                variance = pca_data.get("variance", [])
                loadings = pca_data.get("loadings", [])
                explained_var = pca_data.get("explained_variance", [])
                pca_output = output_dir / "pca_plot.png"
                fig = pca_plot((components, variance, loadings), output_path=pca_output, explained_var=explained_var)
                if fig is not None:
                    plots_generated["pca"] = str(pca_output)
                    plt.close(fig)
        except Exception as e:
            logger.warning(f"PCA plot failed: {e}")

    # Kinship heatmap
    if kinship_file:
        try:
            import json as _json

            kinship_path_obj = Path(kinship_file)
            if kinship_path_obj.exists():
                with open(kinship_path_obj) as fh:
                    kinship_data = _json.load(fh)
                matrix = kinship_data if isinstance(kinship_data, list) else kinship_data.get("matrix", [])
                kinship_output = output_dir / "kinship_plot.png"
                fig = kinship_heatmap(matrix, output_path=kinship_output)
                if fig is not None:
                    plots_generated["kinship"] = str(kinship_output)
                    plt.close(fig)
        except Exception as e:
            logger.warning(f"Kinship heatmap failed: {e}")

    logger.info(f"Generated {len(plots_generated)} plots: {list(plots_generated.keys())}")
    return plots_generated


def missingness_plot(vcf_data: Dict[str, Any], output_path: Optional[Union[str, Path]] = None) -> Any:
    """Create missingness visualization showing per-sample and per-variant missingness.

    Args:
        vcf_data: VCF data dictionary with 'variants' and 'genotypes' keys
        output_path: Path to save the plot (optional)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create missingness plot")
        return None

    logger.info("Creating missingness plot")

    genotypes = vcf_data.get("genotypes", [])
    if not genotypes:
        logger.warning("No genotype data found for missingness plot")
        return None

    genotypes = np.array(genotypes)
    n_variants, n_samples = genotypes.shape

    # Calculate missingness (assuming -1 or negative values indicate missing)
    sample_missingness = np.mean(genotypes < 0, axis=0) * 100
    variant_missingness = np.mean(genotypes < 0, axis=1) * 100

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Per-sample missingness
    ax1.bar(range(n_samples), sample_missingness, color="steelblue", alpha=0.7)
    ax1.axhline(y=5, color="red", linestyle="--", alpha=0.7, label="5% threshold")
    ax1.set_xlabel("Sample Index")
    ax1.set_ylabel("Missing Rate (%)")
    ax1.set_title(f"Per-Sample Missingness (n={n_samples})")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Per-variant missingness histogram
    ax2.hist(variant_missingness, bins=50, edgecolor="black", alpha=0.7, color="steelblue")
    ax2.axvline(x=5, color="red", linestyle="--", alpha=0.7, label="5% threshold")
    ax2.set_xlabel("Missing Rate (%)")
    ax2.set_ylabel("Number of Variants")
    ax2.set_title(f"Per-Variant Missingness Distribution (n={n_variants})")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved missingness plot to {output_path}")

    return fig


def functional_enrichment_plot(
    results: List[Dict[str, Any]], gff_path: Union[str, Path], output_path: Optional[Union[str, Path]] = None
) -> Any:
    """Create functional enrichment plot showing enrichment of significant variants in functional categories.

    Args:
        results: GWAS results with 'chrom', 'pos', and 'p_value' fields
        gff_path: Path to GFF annotation file
        output_path: Path to save the plot (optional)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create functional enrichment plot")
        return None

    logger.info("Creating functional enrichment plot")

    # Count significant variants
    significant = [r for r in results if r.get("p_value", r.get("pval", 1.0)) < 5e-8]

    if not significant:
        logger.warning("No significant variants found for enrichment analysis")
        # Create a simple summary plot instead
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No genome-wide significant variants\n(p < 5e-8)", ha="center", va="center", fontsize=14)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.set_title("Functional Enrichment Analysis")

        if output_path:
            output_path = Path(output_path)
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            logger.info(f"Saved functional enrichment plot to {output_path}")

        return fig

    # Create summary plot of p-value distribution by chromosome
    chrom_counts = {}
    for r in significant:
        chrom = str(r.get("chrom", r.get("chromosome", "unknown")))
        chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1

    fig, ax = plt.subplots(figsize=(10, 6))

    chroms = sorted(chrom_counts.keys(), key=lambda x: int(x) if x.isdigit() else 999)
    counts = [chrom_counts[c] for c in chroms]

    ax.bar(range(len(chroms)), counts, color="steelblue", alpha=0.7)
    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(chroms)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Number of Significant Variants")
    ax.set_title(f"Significant Variants by Chromosome (n={len(significant)})")
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved functional enrichment plot to {output_path}")

    return fig
