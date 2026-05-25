"""GWAS visualization utilities.

This module provides functions for creating GWAS visualization plots,
including Manhattan plots, Q-Q plots, and regional association plots.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Import matplotlib with graceful fallback
try:
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from matplotlib.offsetbox import AnchoredText

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


def set_accessible_style():
    """Apply accessible, high-contrast presentation style for visualizations."""
    if HAS_MATPLOTLIB:
        plt.rcParams.update(
            {
                "font.size": 14,
                "axes.titlesize": 16,
                "axes.labelsize": 14,
                "xtick.labelsize": 12,
                "ytick.labelsize": 12,
                "legend.fontsize": 12,
                "legend.framealpha": 0.85,
                "legend.edgecolor": "black",
                "figure.dpi": 300,
                "axes.grid": True,
                "axes.axisbelow": True,
                "grid.alpha": 0.4,
                "grid.color": "#cccccc",
                # Use STIX mathtext (bundled) to avoid DejaVu Sans dependency
                "mathtext.fontset": "stix",
                # Fallback fonts: prefer system-available sans-serif faces
                "font.sans-serif": [
                    "Helvetica",
                    "Arial",
                    "Lucida Grande",
                    "Verdana",
                    "DejaVu Sans",
                    "Bitstream Vera Sans",
                    "sans-serif",
                ],
            }
        )


if HAS_MATPLOTLIB:
    set_accessible_style()


def _get_visualization_style(style: Optional[Any] = None) -> Any:
    """Return the provided style or load the configured GWAS visualization style."""
    if style is not None:
        return style

    from metainformant.gwas.visualization.config import get_style

    return get_style()


def _ensure_axis(ax: Optional[Any], style: Any, figsize: Optional[Any] = None) -> tuple[Optional[Any], Any]:
    """Return ``(fig, ax)`` while preserving caller-owned axes."""
    if ax is not None:
        return None, ax
    fig, new_ax = plt.subplots(figsize=figsize or style.figsize)
    return fig, new_ax


def _save_owned_figure(
    fig: Optional[Any],
    output_path: Optional[Union[str, Path]],
    style: Any,
    plot_name: str,
    dpi: Optional[int] = None,
) -> None:
    """Save only figures created by this module, not caller-owned axes."""
    if output_path and fig is not None:
        saved_path = Path(output_path)
        fig.savefig(saved_path, dpi=dpi or style.dpi, bbox_inches="tight")
        logger.info(f"Saved {plot_name} to {saved_path}")


def _normalize_result_records(results: Union[List[Dict[str, Any]], Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Normalize GWAS result input to a list of dictionary records."""
    if isinstance(results, dict):
        return [results]
    return [result for result in results if isinstance(result, dict)]


def _neg_log10_p_value(p_value: Any, cap: float = 50.0) -> float:
    """Convert a p-value to -log10 space with a cap for invalid or zero values."""
    try:
        numeric = float(p_value)
    except (TypeError, ValueError):
        return cap
    return -math.log10(numeric) if numeric > 0 else cap


def _sorted_valid_p_values(p_values: Union[List[float], List[int]]) -> Any:
    """Return sorted finite p-values in the open-closed interval (0, 1]."""
    p_vals = np.array(p_values, dtype=float)
    p_vals = p_vals[~np.isnan(p_vals) & (p_vals > 0) & (p_vals <= 1)]
    return np.sort(p_vals)


def manhattan_plot(
    results: Union[List[Dict[str, Any]], Dict[str, Any]],
    output_path: Optional[Union[str, Path]] = None,
    significance_threshold: float = 5e-8,
    gene_annotations: Optional[List[Dict[str, Any]]] = None,
    highlight_regions: Optional[List[Dict[str, Any]]] = None,
    suggestive_threshold: Optional[float] = 1e-5,
    label_top_n: int = 0,
    ax: Optional[Any] = None,
    style: Optional[Any] = None,
    **kwargs: Any,
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

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create Manhattan plot")
        return None

    logger.info("Creating Manhattan plot")

    result_records = _normalize_result_records(results)

    # Extract data for plotting
    chromosomes = []
    positions = []
    p_values = []
    colors = []

    style = _get_visualization_style(style)

    # Create plot
    fig, ax = _ensure_axis(ax, style, kwargs.get("figsize"))

    # High-contrast aesthetic colors for Manhattan chromosomes
    chrom_colors = kwargs.get("colors", ["#0D47A1", "#64B5F6"])

    current_pos = 0
    chrom_offsets = {}

    for result in result_records:
        chrom = str(result.get("chrom", result.get("chromosome", "1")))
        pos = result.get("pos", result.get("position", 0))
        p_val = result.get("p_value", result.get("pval", 1.0))
        neg_log_p = _neg_log10_p_value(p_val)

        # Handle chromosome positioning
        if chrom not in chrom_offsets:
            chrom_offsets[chrom] = current_pos
            current_pos += 100000000  # Space chromosomes apart

        global_pos = chrom_offsets[chrom] + pos
        color_idx = (len(chrom_offsets) - 1) % len(chrom_colors)

        chromosomes.append(chrom)
        positions.append(global_pos)
        p_values.append(neg_log_p)
        colors.append(chrom_colors[color_idx])

    # Plot points with kwargs (auto-rasterize for very dense arrays to prevent 100MB+ PDFs)
    scatter_kwargs = kwargs.get("scatter_kwargs", {"s": style.point_size, "alpha": style.alpha})
    if len(positions) > 100000:
        scatter_kwargs.setdefault("rasterized", True)
        logger.info(f"Auto-rasterizing Manhattan plot due to high variant density (n={len(positions):,})")

    ax.scatter(positions, p_values, c=colors, **scatter_kwargs)

    # Add significance threshold line
    if significance_threshold > 0:
        threshold_line = -math.log10(significance_threshold)
        ax.axhline(
            y=threshold_line,
            color=style.significance_color,
            linestyle="--",
            alpha=0.7,
            label=f"Significance ({significance_threshold})",
        )

    # Add suggestive threshold line
    if suggestive_threshold is not None and suggestive_threshold > 0:
        suggestive_line = -math.log10(suggestive_threshold)
        ax.axhline(
            y=suggestive_line,
            color=style.suggestive_color,
            linestyle="--",
            alpha=0.5,
            label=f"Suggestive ({suggestive_threshold})",
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
                start_x = chrom_offsets[region_chrom] + region_start
                end_x = chrom_offsets[region_chrom] + region_end
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
                if 0 <= idx < len(positions):
                    ann_x = positions[idx]
                    ann_y = p_values[idx]
            # Look up by chrom + pos
            elif "chrom" in annotation and "pos" in annotation:
                ann_chrom = str(annotation["chrom"])
                ann_pos = annotation["pos"]
                if ann_chrom in chrom_offsets:
                    target_x = chrom_offsets[ann_chrom] + ann_pos
                    # Find the closest point
                    best_dist = float("inf")
                    for j, (px, py) in enumerate(zip(positions, p_values)):
                        dist = abs(px - target_x)
                        if dist < best_dist:
                            best_dist = dist
                            ann_x = px
                            ann_y = py

            if ann_x is not None and ann_y is not None:
                offset = max(p_values) * 0.08 if p_values else 1.0
                ax.annotate(
                    gene_name,
                    xy=(ann_x, ann_y),
                    xytext=(ann_x, ann_y + offset),
                    arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
                    fontsize=8,
                    ha="center",
                )

    # Label top N hits
    if label_top_n > 0 and positions:
        # Build index of (neg_log_p, position_idx) and sort by neg_log_p descending
        indexed = sorted(enumerate(p_values), key=lambda x: x[1], reverse=True)
        top_indices = [idx for idx, _ in indexed[:label_top_n]]
        offset = max(p_values) * 0.06 if p_values else 1.0

        for rank, idx in enumerate(top_indices):
            px = positions[idx]
            py = p_values[idx]
            # Use variant_id if available, otherwise chr:pos
            result_entry = result_records[idx] if idx < len(result_records) else {}
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
    chrom_centers = {}
    for chrom in sorted(chrom_offsets.keys(), key=lambda x: int(x) if x.isdigit() else 999):
        center = chrom_offsets[chrom] + 50000000  # Approximate center
        chrom_centers[center] = chrom
    ax.set_xticks(list(chrom_centers.keys()))
    ax.set_xticklabels(list(chrom_centers.values()))

    # Labels and title
    ax.set_xlabel("Chromosome", fontsize=14)
    ax.set_ylabel("-log\u2081\u2080(p-value)", fontsize=14)
    n_tested = len(positions)
    sig_threshold_log = -math.log10(significance_threshold) if significance_threshold > 0 else 50
    n_sig = sum(1 for pv in p_values if pv >= sig_threshold_log)
    ax.set_title(
        f"Manhattan Plot (n={n_tested:,} variants, {n_sig} GW-significant)",
        fontsize=16,
        fontweight="bold",
    )
    ax.tick_params(axis="both", labelsize=12)
    ax.grid(True, alpha=0.3)

    # Create legend for colors
    legend_elements = [
        mpatches.Patch(color=color, label=f"Chr {chrom}")
        for chrom, color in zip(sorted(chrom_offsets.keys()), chrom_colors)
    ]
    # Collect any auto-labeled artists (threshold lines, highlight regions)
    auto_handles, auto_labels = ax.get_legend_handles_labels()
    for handle, lbl in zip(auto_handles, auto_labels):
        if lbl and lbl not in [e.get_label() for e in legend_elements]:
            legend_elements.append(handle)
    if legend_elements:
        ax.legend(handles=legend_elements, loc="upper right")

    plt.tight_layout()

    _save_owned_figure(fig, output_path, style, "Manhattan plot")

    return ax if fig is None else fig


def qq_plot(
    p_values: Union[List[float], List[int]],
    output_path: Optional[Union[str, Path]] = None,
    ax: Optional[Any] = None,
    style: Optional[Any] = None,
    lambda_gc: Optional[float] = None,
    n_samples: Optional[int] = None,
    **kwargs: Any,
) -> Any:
    """Create a publication-quality Q-Q plot with comprehensive power diagnostics.

    Includes: λ_GC with N-adjusted interpretation, mean χ², KS goodness-of-fit,
    GC-corrected p-value overlay, 95% CI band, and power context annotation.

    Args:
        p_values: List of p-values from association tests.
        output_path: Path to save the plot (optional).
        lambda_gc: Pre-computed genomic inflation factor (computed if None).
        n_samples: Number of samples in the study (for power context).

    Returns:
        matplotlib Figure object.
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create Q-Q plot")
        return None

    logger.info("Creating enhanced Q-Q plot")

    p_vals = _sorted_valid_p_values(p_values)

    if len(p_vals) == 0:
        logger.warning("No valid p-values for Q-Q plot")
        return None

    n = len(p_vals)
    expected = np.arange(1, n + 1) / (n + 1)
    expected_log = -np.log10(expected)
    observed_log = -np.log10(p_vals)

    # ── Compute λ_GC ──
    try:
        from scipy import stats as scipy_stats

        _has_scipy = True
    except ImportError:
        _has_scipy = False

    if lambda_gc is None:
        try:
            from metainformant.gwas.analysis.correction import genomic_control

            gc_result = genomic_control(p_values=p_vals.tolist())
            lambda_gc = gc_result.get("lambda_gc", 1.0)
        except (ImportError, Exception):
            median_p = float(np.median(p_vals))
            if 0 < median_p < 1 and _has_scipy:
                chi2_median = float(scipy_stats.chi2.ppf(1.0 - median_p, df=1))
                lambda_gc = chi2_median / 0.4549364
            else:
                lambda_gc = 1.0

    # ── Compute additional diagnostic statistics ──
    mean_chi2 = None
    ks_p = None
    gc_corrected_log = None

    if _has_scipy:
        # Mean χ² (should be ~1.0 under null)
        chi2_vals = scipy_stats.chi2.ppf(1.0 - p_vals, df=1)
        chi2_vals = chi2_vals[np.isfinite(chi2_vals)]
        mean_chi2 = float(np.mean(chi2_vals)) if len(chi2_vals) > 0 else None

        # KS test of p-values against Uniform(0,1)
        ks_stat, ks_p = scipy_stats.kstest(p_vals, "uniform")

        # GC-corrected p-values for overlay
        if lambda_gc > 1.0:
            corrected_chi2 = chi2_vals / lambda_gc
            corrected_p = 1.0 - scipy_stats.chi2.cdf(corrected_chi2, df=1)
            corrected_p = np.sort(corrected_p)
            gc_corrected_log = -np.log10(np.clip(corrected_p, 1e-300, 1.0))

    style = _get_visualization_style(style)

    # ── Create figure ──
    fig, ax = _ensure_axis(ax, style, kwargs.get("figsize", (9, 9)))

    # 95% CI band
    if _has_scipy:
        lower_ci = -np.log10(scipy_stats.beta.ppf(0.975, np.arange(1, n + 1), np.arange(n, 0, -1)))
        upper_ci = -np.log10(scipy_stats.beta.ppf(0.025, np.arange(1, n + 1), np.arange(n, 0, -1)))
        ax.fill_between(
            expected_log,
            lower_ci,
            upper_ci,
            alpha=0.12,
            color="#90CAF9",
            label="95% CI (null)",
            zorder=1,
        )

    # GC-corrected overlay (if inflated)
    if gc_corrected_log is not None and lambda_gc > 1.05:
        ax.scatter(
            expected_log[: len(gc_corrected_log)],
            gc_corrected_log,
            s=4,
            alpha=0.4,
            color="#66BB6A",
            zorder=2,
            label=f"GC-corrected (λ={lambda_gc:.2f})",
        )

    # Observed Q-Q points
    s_size = max(4, min(20, 600 // max(1, int(n**0.5))))
    ax.scatter(
        expected_log,
        observed_log,
        s=s_size,
        alpha=0.75,
        color="#1565C0",
        edgecolors="none",
        zorder=4,
        label="Observed",
    )

    # Diagonal (expected under null)
    max_val = max(float(np.max(expected_log)), float(np.max(observed_log))) * 1.05
    ax.plot(
        [0, max_val],
        [0, max_val],
        linestyle="--",
        color="#E53935",
        alpha=0.7,
        linewidth=2,
        label="Expected (null)",
        zorder=3,
    )

    # ── Inflation status with N-aware interpretation ──
    if lambda_gc < 1.05:
        inflation_status = "✓ Nominal"
        status_color = "#E8F5E9"
        border_color = "#4CAF50"
    elif lambda_gc < 1.10:
        inflation_status = "~ Modest inflation"
        status_color = "#FFFDE7"
        border_color = "#FFC107"
    elif lambda_gc < 1.20:
        inflation_status = "⚠ Moderate inflation"
        status_color = "#FFF3E0"
        border_color = "#FF9800"
    else:
        inflation_status = "⚠ Elevated λ_GC"
        status_color = "#FFEBEE"
        border_color = "#F44336"

    # Context-aware interpretation for small N
    context_note = ""
    if n_samples is not None and n_samples < 100 and lambda_gc > 1.1:
        context_note = f"\nN={n_samples}: small-sample\nλ_GC often elevated"
        if lambda_gc < 1.5:
            inflation_status = "~ Expected (small N)"
            status_color = "#FFF3E0"
            border_color = "#FF9800"

    # Build annotation text
    gc_lines = [
        f"λ_GC = {lambda_gc:.4f}  {inflation_status}",
        f"n = {n:,} variants tested",
    ]
    if n_samples is not None:
        gc_lines.append(f"N = {n_samples} samples")
    if mean_chi2 is not None:
        gc_lines.append(f"mean χ² = {mean_chi2:.3f}  (expect ≈1.0)")
    if ks_p is not None:
        gc_lines.append(f"KS test p = {ks_p:.2e}")
    if context_note:
        gc_lines.append(context_note)

    gc_text = "\n".join(gc_lines)
    at = AnchoredText(
        gc_text,
        loc="upper left",
        prop=dict(size=10, fontweight="bold", fontfamily="monospace"),
        frameon=True,
    )
    at.patch.set_boxstyle("round,pad=0.4")
    at.patch.set_facecolor(status_color)
    at.patch.set_alpha(0.95)
    at.patch.set_edgecolor(border_color)
    at.patch.set_linewidth(2)
    ax.add_artist(at)

    # ── Axes formatting ──
    ax.set_xlabel("Expected −log₁₀(p)", fontsize=14, fontweight="medium")
    ax.set_ylabel("Observed −log₁₀(p)", fontsize=14, fontweight="medium")
    title_suffix = f" — {n_samples} samples" if n_samples else ""
    ax.set_title(
        f"GWAS Q-Q Plot  (n={n:,} variants{title_suffix})",
        fontsize=15,
        fontweight="bold",
    )
    ax.tick_params(axis="both", labelsize=12)
    ax.grid(True, alpha=0.25, linewidth=0.5)
    ax.set_xlim(left=-0.05)
    ax.set_ylim(bottom=-0.05)
    ax.legend(fontsize=9, loc="lower right", framealpha=0.95, edgecolor="#ccc", fancybox=True)

    plt.tight_layout()

    _save_owned_figure(fig, output_path, style, "Q-Q plot")

    return ax if fig is None else fig


def regional_plot(
    results: List[Dict[str, Any]],
    chrom: str,
    start: int,
    end: int,
    output_path: Optional[Union[str, Path]] = None,
    ax: Optional[Any] = None,
    style: Optional[Any] = None,
    **kwargs: Any,
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
    for result in _normalize_result_records(results):
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
        p_values.append(_neg_log10_p_value(p_val))

    style = _get_visualization_style(style)

    # Create plot
    fig, ax = _ensure_axis(ax, style, kwargs.get("figsize"))

    # Plot points
    scatter_kwargs = kwargs.get("scatter_kwargs", {"s": style.point_size, "alpha": style.alpha})

    # LocusZoom style mapping if LD values are provided
    ld_dict = kwargs.get("ld_values", None)

    if ld_dict and isinstance(ld_dict, dict):
        ld_vals = []
        for r in region_results:
            vid = r.get("variant_id", "")
            # fallback to chr:pos string matching if exact ID not found
            if vid not in ld_dict:
                vid = f"{r.get('chrom', '')}:{r.get('pos', 0)}"
            ld_vals.append(ld_dict.get(vid, 0.0))

        scatter = ax.scatter(positions, p_values, c=ld_vals, cmap="viridis_r", **scatter_kwargs)
        cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
        cbar.set_label("LD ($r^2$ to top hit)")
    else:
        scatter_kwargs.setdefault("color", style.suggestive_color)
        ax.scatter(positions, p_values, **scatter_kwargs)

    # Add significance threshold line
    threshold_line = -math.log10(5e-8)
    ax.axhline(
        y=threshold_line,
        color=style.significance_color,
        linestyle="--",
        alpha=0.7,
        label="Genome-wide significance",
    )

    # Labels and title
    ax.set_xlabel(f"Position on chromosome {chrom}")
    ax.set_ylabel("-log₁₀(p-value)")
    ax.set_title(f"Regional Association Plot: {chrom}:{start:,}-{end:,}")
    ax.set_xlim(start, end)
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    _save_owned_figure(fig, output_path, style, "regional plot")

    return ax if fig is None else fig


# ── Strain palette shared across all GWAS visualizations ────────────────
STRAIN_PALETTE = {
    "C": "#2196F3",  # Blue — Carniolan
    "I": "#FF9800",  # Orange — Italian
    "M": "#4CAF50",  # Green — A.m. mellifera
    "R": "#F44336",  # Red — Russian
}
STRAIN_NAMES = {
    "C": "Carniolan",
    "I": "Italian",
    "M": "A.m. mellifera",
    "R": "Russian",
}
STRAIN_MARKERS = {"C": "o", "I": "s", "M": "D", "R": "^"}


def _infer_strain(sample_id: str) -> str:
    """Infer strain code from sample ID prefix (e.g. 'C15ITQ' → 'C')."""
    s = sample_id.strip()
    if s and s[0] in STRAIN_PALETTE:
        return s[0]
    return "?"


def pca_plot(
    pca_result: tuple,
    output_path: Optional[Union[str, Path]] = None,
    explained_var: Optional[List[float]] = None,
    sample_ids: Optional[List[str]] = None,
    ax: Optional[Any] = None,
    style: Optional[Any] = None,
    **kwargs: Any,
) -> Any:
    """Create a publication-quality PCA plot with strain coloring and diagnostics.

    Features:
      • Strain-colored points with distinct markers and per-strain legend
      • 95 % confidence ellipses per strain
      • Sample-ID labels (toggled via show_labels kwarg)
      • Scree-plot inset showing top-10 PC variance
      • Title with sample/variant counts

    Args:
        pca_result: (components, variance, loadings) tuple.
        output_path: Save path (optional).
        explained_var: Explained variance ratios per PC.
        sample_ids: List of sample identifiers (same order as columns).
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create PCA plot")
        return None

    logger.info("Creating strain-aware PCA plot")

    try:
        components, variance, loadings = pca_result

        if len(components) < 2:
            logger.warning("Need at least 2 PCA components for plotting")
            return None

        if style is None:
            from metainformant.gwas.visualization.config import get_style

            style = get_style()

        pc1 = np.array(components[0], dtype=float)
        pc2 = np.array(components[1], dtype=float)
        n_samples = len(pc1)

        # ── Assign strains ──
        if sample_ids and len(sample_ids) >= n_samples:
            ids = sample_ids[:n_samples]
        else:
            ids = [f"S{i}" for i in range(n_samples)]
        strains = [_infer_strain(sid) for sid in ids]
        unique_strains = sorted(set(strains) - {"?"})

        # Variance percentages
        pc1_var = explained_var[0] * 100 if explained_var and len(explained_var) > 0 else 0
        pc2_var = explained_var[1] * 100 if explained_var and len(explained_var) > 1 else 0

        # ── Figure ──
        fig = None
        if ax is None:
            fig, ax = plt.subplots(figsize=kwargs.get("figsize", (11, 9)))

        show_labels = kwargs.get("show_labels", n_samples <= 60)
        show_ellipses = kwargs.get("show_ellipses", True)

        # ── Plot per strain ──
        for strain in unique_strains:
            mask = np.array([s == strain for s in strains])
            if not np.any(mask):
                continue
            x, y = pc1[mask], pc2[mask]
            color = STRAIN_PALETTE.get(strain, "#999")
            marker = STRAIN_MARKERS.get(strain, "o")
            label = f"{strain} — {STRAIN_NAMES.get(strain, strain)} (n={int(np.sum(mask))})"

            ax.scatter(
                x,
                y,
                c=color,
                marker=marker,
                s=70,
                alpha=0.85,
                edgecolors="white",
                linewidths=0.5,
                label=label,
                zorder=4,
            )

            # 95% confidence ellipse
            if show_ellipses and np.sum(mask) >= 3:
                try:
                    from matplotlib.patches import Ellipse

                    cov = np.cov(x, y)
                    if np.all(np.isfinite(cov)):
                        eigenvalues, eigenvectors = np.linalg.eigh(cov)
                        order = eigenvalues.argsort()[::-1]
                        eigenvalues = eigenvalues[order]
                        eigenvectors = eigenvectors[:, order]
                        angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
                        width, height = 2 * np.sqrt(eigenvalues * 5.991)  # chi2(2, 0.95)
                        ell = Ellipse(
                            xy=(np.mean(x), np.mean(y)),
                            width=width,
                            height=height,
                            angle=angle,
                            facecolor=color,
                            alpha=0.08,
                            edgecolor=color,
                            linewidth=1.5,
                            linestyle="--",
                            zorder=2,
                        )
                        ax.add_patch(ell)
                except Exception:
                    pass  # Skip ellipse on error

            # Sample labels
            if show_labels:
                for xi, yi, sid in zip(x, y, np.array(ids)[mask]):
                    ax.annotate(
                        sid,
                        (xi, yi),
                        fontsize=6,
                        alpha=0.7,
                        textcoords="offset points",
                        xytext=(4, 4),
                    )

        # Unknown strain (if any)
        unk_mask = np.array([s == "?" for s in strains])
        if np.any(unk_mask):
            ax.scatter(
                pc1[unk_mask],
                pc2[unk_mask],
                c="#999",
                marker="x",
                s=50,
                alpha=0.6,
                label=f"Unknown (n={int(np.sum(unk_mask))})",
                zorder=3,
            )

        # ── Axes ──
        ax.set_xlabel(
            f"PC1  ({pc1_var:.1f}% variance explained)",
            fontsize=13,
            fontweight="medium",
        )
        ax.set_ylabel(
            f"PC2  ({pc2_var:.1f}% variance explained)",
            fontsize=13,
            fontweight="medium",
        )
        ax.set_title(
            f"Population Structure — PCA  ({n_samples} samples)",
            fontsize=15,
            fontweight="bold",
        )
        ax.tick_params(axis="both", labelsize=11)
        ax.grid(True, alpha=0.2, linewidth=0.5)
        ax.axhline(0, color="#aaa", linewidth=0.4, zorder=1)
        ax.axvline(0, color="#aaa", linewidth=0.4, zorder=1)
        ax.legend(
            fontsize=9,
            loc="best",
            framealpha=0.95,
            edgecolor="#ccc",
            fancybox=True,
            title="Lineage",
            title_fontsize=10,
        )

        # ── Scree inset ──
        if explained_var and len(explained_var) >= 3:
            inset_ax = ax.inset_axes([0.72, 0.02, 0.26, 0.22])  # [x, y, w, h]
            n_pcs = min(10, len(explained_var))
            pcs_idx = np.arange(1, n_pcs + 1)
            var_pct = [v * 100 for v in explained_var[:n_pcs]]
            cum_var = np.cumsum(var_pct)
            inset_ax.bar(pcs_idx, var_pct, color="#1565C0", alpha=0.7, width=0.7)
            inset_ax.plot(pcs_idx, cum_var, "o-", color="#E53935", markersize=3, linewidth=1)
            inset_ax.set_xlabel("PC", fontsize=7)
            inset_ax.set_ylabel("%Var", fontsize=7)
            inset_ax.set_title("Scree", fontsize=8, fontweight="bold")
            inset_ax.tick_params(labelsize=6)
            inset_ax.set_xticks(pcs_idx)
            inset_ax.grid(True, alpha=0.2)

        plt.tight_layout()

        if output_path and fig is not None:
            output_path = Path(output_path)
            fig.savefig(output_path, dpi=style.dpi, bbox_inches="tight")
            logger.info(f"Saved PCA plot to {output_path}")

        return ax if fig is None else fig

    except (ValueError, IndexError, TypeError) as e:
        logger.error(f"Error creating PCA plot: {e}")
        return None


def kinship_heatmap(
    kinship_matrix: Union[np.ndarray, List[List[float]]],
    output_path: Optional[Union[str, Path]] = None,
    sample_ids: Optional[List[str]] = None,
    ax: Optional[Any] = None,
    style: Optional[Any] = None,
    **kwargs: Any,
) -> Any:
    """Create kinship matrix heatmap with hierarchical clustering dendrograms.

    Features:
      - Dendrograms on both axes (hierarchical clustering)
      - Sample reordering by cluster similarity
      - Strain-colored sidebar (C=blue, I=orange, M=green, R=red)
      - Publication-quality typography

    Args:
        kinship_matrix: Kinship matrix (n_samples x n_samples)
        output_path: Path to save the plot (optional)
        sample_ids: Optional sample ID labels for strain coloring
        ax: Ignored when dendrograms are enabled (needs custom GridSpec)
        style: PlotStyle configuration

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create kinship heatmap")
        return None

    logger.info("Creating kinship heatmap")

    try:
        from metainformant.gwas.visualization.strain_plots import (
            STRAIN_NAMES,
            STRAIN_PALETTE,
            _extract_strain,
        )
    except ImportError:
        STRAIN_PALETTE = {}
        STRAIN_NAMES = {}
        _extract_strain = lambda s: "?"  # noqa: E731

    try:
        # Convert to numpy array if needed
        if isinstance(kinship_matrix, list):
            kinship_matrix = np.array(kinship_matrix, dtype=float)

        if style is None:
            from metainformant.gwas.visualization.config import get_style

            style = get_style()

        n_samples = kinship_matrix.shape[0]

        # Try hierarchical clustering for dendrograms
        try:
            from scipy.cluster.hierarchy import dendrogram, linkage
            from scipy.spatial.distance import squareform

            has_scipy = True
        except ImportError:
            has_scipy = False

        if has_scipy and n_samples >= 3:
            # Convert kinship to distance: d = 1 - kinship (capped at 0)
            dist_matrix = 1.0 - np.clip(kinship_matrix, 0, 1)
            np.fill_diagonal(dist_matrix, 0)
            # Ensure symmetry
            dist_matrix = (dist_matrix + dist_matrix.T) / 2.0

            condensed = squareform(dist_matrix, checks=False)
            Z = linkage(condensed, method="average")

            # Layout: [row_dend | strain_bar | heatmap | colorbar]
            #          top_dend above heatmap
            fig = plt.figure(figsize=kwargs.get("figsize", (12, 10)))

            if sample_ids:
                gs = fig.add_gridspec(
                    2,
                    4,
                    width_ratios=[0.15, 0.03, 1, 0.04],
                    height_ratios=[0.15, 1],
                    hspace=0.02,
                    wspace=0.02,
                )
            else:
                gs = fig.add_gridspec(
                    2,
                    3,
                    width_ratios=[0.15, 1, 0.04],
                    height_ratios=[0.15, 1],
                    hspace=0.02,
                    wspace=0.02,
                )

            # Top dendrogram
            ax_top = fig.add_subplot(gs[0, -2] if sample_ids else gs[0, 1])
            dn_top = dendrogram(
                Z,
                ax=ax_top,
                orientation="top",
                no_labels=True,
                color_threshold=0,
                above_threshold_color="#333333",
            )
            ax_top.set_axis_off()
            order = dn_top["leaves"]

            # Left dendrogram
            ax_left = fig.add_subplot(gs[1, 0])
            dendrogram(
                Z,
                ax=ax_left,
                orientation="left",
                no_labels=True,
                color_threshold=0,
                above_threshold_color="#333333",
            )
            ax_left.set_axis_off()
            ax_left.invert_yaxis()

            # Reorder matrix
            reordered = kinship_matrix[np.ix_(order, order)]

            # Strain color sidebar
            if sample_ids:
                ax_bar = fig.add_subplot(gs[1, 1])
                strains = [_extract_strain(sample_ids[i]) for i in order]
                strain_colors = [STRAIN_PALETTE.get(s, "#888888") for s in strains]
                for i, color in enumerate(strain_colors):
                    ax_bar.add_patch(plt.Rectangle((0, i), 1, 1, color=color))
                ax_bar.set_xlim(0, 1)
                ax_bar.set_ylim(0, len(strains))
                ax_bar.set_axis_off()
                ax_heat = fig.add_subplot(gs[1, 2])
            else:
                ax_heat = fig.add_subplot(gs[1, 1])

            # Main heatmap
            im = ax_heat.imshow(reordered, cmap="viridis", aspect="equal", interpolation="nearest")

            # Tick labels
            if sample_ids and n_samples <= 60:
                labels = [sample_ids[i] for i in order]
                ax_heat.set_xticks(range(n_samples))
                ax_heat.set_xticklabels(labels, rotation=90, fontsize=max(5, 9 - n_samples // 15))
                ax_heat.set_yticks(range(n_samples))
                ax_heat.set_yticklabels(labels, fontsize=max(5, 9 - n_samples // 15))
            else:
                ax_heat.set_xticks([])
                ax_heat.set_yticks([])

            # Colorbar
            ax_cbar = fig.add_subplot(gs[1, -1])
            plt.colorbar(im, cax=ax_cbar, label="Kinship coefficient")

            fig.suptitle("Kinship Matrix (Clustered)", fontsize=18, fontweight="bold", y=0.98)

            # Strain legend
            if sample_ids and STRAIN_PALETTE:
                from matplotlib.patches import Patch

                unique_strains = sorted(set(_extract_strain(sample_ids[i]) for i in order))
                legend_patches = [
                    Patch(
                        facecolor=STRAIN_PALETTE.get(s, "#888"),
                        label=STRAIN_NAMES.get(s, s),
                    )
                    for s in unique_strains
                    if s != "?"
                ]
                if legend_patches:
                    fig.legend(
                        handles=legend_patches,
                        loc="lower right",
                        fontsize=10,
                        title="Strain",
                        title_fontsize=12,
                        bbox_to_anchor=(0.98, 0.02),
                    )

        else:
            # Fallback: simple heatmap without dendrograms
            fig, ax_heat = plt.subplots(figsize=kwargs.get("figsize", (10, 8)))
            im = ax_heat.imshow(kinship_matrix, cmap="viridis", aspect="equal")
            plt.colorbar(im, ax=ax_heat, label="Kinship coefficient")
            ax_heat.set_title("Kinship Matrix Heatmap", fontsize=16, fontweight="bold")
            ax_heat.set_xlabel("Sample")
            ax_heat.set_ylabel("Sample")

        plt.tight_layout()

        if output_path and fig is not None:
            output_path = Path(output_path)
            fig.savefig(output_path, dpi=style.dpi, bbox_inches="tight")
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
        x=np.mean(effect_sizes),
        color="green",
        linestyle="-",
        alpha=0.7,
        label=f"Mean: {np.mean(effect_sizes):.4f}",
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

    if isinstance(association_results, list):
        results_data = association_results

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

    # Extract sample IDs from VCF for strain coloring and sample count
    _vcf_sample_ids: List[str] = []
    if vcf_file:
        try:
            import subprocess as _sp

            _qr = _sp.run(
                ["bcftools", "query", "-l", str(vcf_file)],
                capture_output=True,
                text=True,
                timeout=10,
            )
            if _qr.returncode == 0:
                raw_lines = _qr.stdout.strip().splitlines()
                _vcf_sample_ids = [line for line in raw_lines if line and not line.startswith("[")]
        except Exception:
            pass  # bcftools not available

    if results_data:
        # Manhattan plot
        try:
            p_vals_manhattan = [
                r.get("p_value", r.get("pval", r.get("pvalue")))
                for r in results_data
                if r.get("p_value", r.get("pval", r.get("pvalue"))) is not None
            ]
            if p_vals_manhattan:
                manhattan_path = output_dir / "manhattan_plot.png"
                fig = manhattan_plot(
                    results_data,
                    output_path=manhattan_path,
                    significance_threshold=significance_threshold,
                )
                if fig is not None:
                    plots_generated["manhattan"] = str(manhattan_path)
                    plt.close(fig)
        except Exception as e:
            logger.warning(f"Manhattan plot failed: {e}")

        # Q-Q plot — with N-aware diagnostics
        try:
            p_vals = [
                r.get("p_value", r.get("pval", r.get("pvalue")))
                for r in results_data
                if r.get("p_value", r.get("pval", r.get("pvalue"))) is not None
            ]
            if p_vals:
                qq_path = output_dir / "qq_plot.png"
                fig = qq_plot(
                    p_vals,
                    output_path=qq_path,
                    n_samples=len(_vcf_sample_ids) if _vcf_sample_ids else None,
                )
                if fig is not None:
                    plots_generated["qq"] = str(qq_path)
                    plt.close(fig)
        except Exception as e:
            logger.warning(f"Q-Q plot failed: {e}")

    # PCA plot — uses _vcf_sample_ids extracted above for strain coloring
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
                fig = pca_plot(
                    (components, variance, loadings),
                    output_path=pca_output,
                    explained_var=explained_var,
                    sample_ids=_vcf_sample_ids or None,
                )
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

                # Handle both direct matrix list and wrapped dict
                if isinstance(kinship_data, list):
                    matrix = kinship_data
                elif isinstance(kinship_data, dict):
                    matrix = kinship_data.get("matrix", kinship_data.get("kinship_matrix", []))
                else:
                    matrix = []
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
    if len(genotypes) == 0:
        logger.warning("No genotype data found for missingness plot")
        return None

    genotypes = np.array(genotypes)
    # Detect major axis to avoid plotting millions of geometric bars
    if genotypes.shape[0] < genotypes.shape[1]:
        n_samples, n_variants = genotypes.shape
        variant_missingness = np.mean(genotypes < 0, axis=0) * 100
        sample_missingness = np.mean(genotypes < 0, axis=1) * 100
    else:
        n_variants, n_samples = genotypes.shape
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
    ax1.legend(loc="upper right")
    ax1.grid(True, alpha=0.3)

    # Per-variant missingness histogram
    ax2.hist(variant_missingness, bins=50, edgecolor="black", alpha=0.7, color="steelblue")
    ax2.axvline(x=5, color="red", linestyle="--", alpha=0.7, label="5% threshold")
    ax2.set_xlabel("Missing Rate (%)")
    ax2.set_ylabel("Number of Variants")
    ax2.set_title(f"Per-Variant Missingness Distribution (n={n_variants})")
    ax2.legend(loc="upper right")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved missingness plot to {output_path}")

    return fig


def functional_enrichment_plot(
    results: List[Dict[str, Any]],
    gff_path: Union[str, Path],
    output_path: Optional[Union[str, Path]] = None,
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
        ax.text(
            0.5,
            0.5,
            "No genome-wide significant variants\n(p < 5e-8)",
            ha="center",
            va="center",
            fontsize=14,
        )
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


# ---------------------------------------------------------------------------
# Power & Convergence Plots
# ---------------------------------------------------------------------------


def power_curve_plot(
    power_data: List[Dict[str, Any]],
    output_path: Optional[Union[str, Path]] = None,
) -> Optional[Any]:
    """Plot statistical power vs. sample size for multiple effect sizes.

    Args:
        power_data: List of dicts from power_curve(), each with
            sample_sizes, powers, beta, and maf.
        output_path: Path to save figure.

    Returns:
        matplotlib Figure or None.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available")
        return None

    logger.info("Creating power curve plot")

    fig, ax = plt.subplots(figsize=(10, 6))

    colors = (
        plt.cm.viridis(np.linspace(0.1, 0.9, len(power_data)))
        if HAS_NUMPY
        else [f"C{i}" for i in range(len(power_data))]
    )

    for idx, curve in enumerate(power_data):
        sizes = curve.get("sample_sizes", [])
        powers = curve.get("powers", [])
        beta = curve.get("beta", 0)
        maf = curve.get("maf", 0)
        ax.plot(
            sizes,
            powers,
            color=colors[idx],
            linewidth=2,
            marker="o",
            markersize=4,
            label=f"β={beta}, MAF={maf:.2f}",
        )

    ax.axhline(y=0.8, color="gray", linestyle="--", alpha=0.6, label="80% power")
    ax.set_xlabel("Sample Size (N)")
    ax.set_ylabel("Statistical Power")
    ax.set_title("GWAS Power Curves")
    ax.set_ylim(-0.02, 1.05)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved power curve plot to {output_path}")

    return fig


def convergence_plot(
    convergence_data: Dict[str, Any],
    metrics: Optional[List[str]] = None,
    output_path: Optional[Union[str, Path]] = None,
) -> Optional[Any]:
    """Plot GWAS metric convergence across data fractions.

    Creates a multi-panel figure with one subplot per metric, showing
    mean ± std as error ribbons.

    Args:
        convergence_data: Output from subsample_convergence().
        metrics: Which metrics to plot (default: all available).
        output_path: Path to save figure.

    Returns:
        matplotlib Figure or None.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available")
        return None

    logger.info("Creating convergence plot")

    fractions = convergence_data.get("fractions", [])
    all_metrics = convergence_data.get("metrics", {})
    pct_fractions = [f * 100 for f in fractions]

    if metrics is None:
        metrics = list(all_metrics.keys())

    n_panels = len(metrics)
    if n_panels == 0:
        return None

    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5), squeeze=False)

    metric_labels = {
        "lambda_gc": "Genomic Inflation (λ_GC)",
        "n_significant": "Significant Hits",
        "mean_abs_beta": "Mean |β|",
        "mean_neg_log10_p": "Mean -log₁₀(p)",
    }
    metric_colors = {
        "lambda_gc": "#E74C3C",
        "n_significant": "#3498DB",
        "mean_abs_beta": "#2ECC71",
        "mean_neg_log10_p": "#9B59B6",
    }

    for i, metric in enumerate(metrics):
        ax = axes[0, i]
        data = all_metrics.get(metric, [])
        if not data:
            continue

        means = [d["mean"] for d in data]
        stds = [d["std"] for d in data]
        lower = [m - s for m, s in zip(means, stds)]
        upper = [m + s for m, s in zip(means, stds)]

        color = metric_colors.get(metric, "steelblue")
        ax.plot(pct_fractions, means, color=color, linewidth=2, marker="o", markersize=5)
        ax.fill_between(pct_fractions, lower, upper, alpha=0.2, color=color)

        ax.set_xlabel("% of Variants Used")
        ax.set_ylabel(metric_labels.get(metric, metric))
        ax.set_title(metric_labels.get(metric, metric))
        ax.grid(True, alpha=0.3)

    fig.suptitle("Subsampling Convergence Analysis", fontsize=14, fontweight="bold")
    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved convergence plot to {output_path}")

    return fig


def saturation_plot(
    saturation_data: Dict[str, Any],
    output_path: Optional[Union[str, Path]] = None,
) -> Optional[Any]:
    """Plot observed data vs. fitted saturation curve.

    Shows the exponential saturation model K(f) = K_∞(1 - e^{-κf})
    with the K_∞ asymptote and convergence indicators.

    Args:
        saturation_data: Output from saturation_analysis().
        output_path: Path to save figure.

    Returns:
        matplotlib Figure or None.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available")
        return None

    logger.info("Creating saturation plot")

    fractions = saturation_data.get("fractions", [])
    observed = saturation_data.get("observed", [])
    predicted = saturation_data.get("predicted", [])
    k_inf = saturation_data.get("k_inf", 0)
    kappa = saturation_data.get("kappa", 0)
    r_sq = saturation_data.get("r_squared", 0)
    is_saturated = saturation_data.get("is_saturated", False)
    metric = saturation_data.get("metric", "metric")

    pct_fractions = [f * 100 for f in fractions]

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.scatter(pct_fractions, observed, color="#3498DB", s=60, zorder=3, label="Observed")
    ax.plot(
        pct_fractions,
        predicted,
        color="#E74C3C",
        linewidth=2,
        linestyle="--",
        label="Fitted curve",
    )
    ax.axhline(
        y=k_inf,
        color="#2ECC71",
        linestyle=":",
        linewidth=1.5,
        label=f"K_∞ = {k_inf:.3f}",
    )

    status_text = "✓ SATURATED" if is_saturated else "✗ NOT SATURATED"
    status_color = "#2ECC71" if is_saturated else "#E74C3C"
    ax.text(
        0.02,
        0.95,
        status_text,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        color=status_color,
        va="top",
    )

    ax.set_xlabel("% of Variants Used")
    ax.set_ylabel(metric.replace("_", " ").title())
    ax.set_title(f"Saturation Analysis — {metric}\n(κ={kappa:.2f}, R²={r_sq:.4f})")
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved saturation plot to {output_path}")

    return fig


def forest_plot(
    results: "List[Dict[str, Any]]",
    output_path: "Optional[Union[str, Path]]" = None,
    top_n: int = 10,
    significance_threshold: float = 5e-8,
    ax: "Optional[Any]" = None,
    style: "Optional[Any]" = None,
    **kwargs: "Any",
) -> "Any":
    """Forest plot of effect sizes (β ± 1.96·SE) for top GWAS hits.

    Each row shows one variant ordered by p-value (most significant at top).
    GW-significant hits shown in red; sub-threshold in blue.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for forest plot")
        return None
    if not results:
        return None

    if style is None:
        from metainformant.gwas.visualization.config import get_style

        style = get_style()

    sorted_r = sorted(results, key=lambda r: r.get("p_value", 1.0))[:top_n]
    plot_r = list(reversed(sorted_r))
    n = len(plot_r)

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get("figsize", (10, max(4, n * 0.55 + 1.5))))

    betas_all = [r.get("beta", 0.0) for r in plot_r]
    ses_all = [r.get("se", 0.0) for r in plot_r]
    ci_los = [b - 1.96 * s for b, s in zip(betas_all, ses_all)]
    ci_his = [b + 1.96 * s for b, s in zip(betas_all, ses_all)]
    x_min = min(ci_los) * 1.15
    x_max = max(ci_his) * 1.15

    for i, r in enumerate(plot_r):
        beta = r.get("beta", 0.0)
        ci_lo = ci_los[i]
        ci_hi = ci_his[i]
        p_val = r.get("p_value", 1.0)
        is_sig = p_val < significance_threshold
        color = style.significance_color if is_sig else style.suggestive_color
        marker = "D" if is_sig else "o"

        ax.plot([ci_lo, ci_hi], [i, i], "-", color=color, linewidth=1.8, alpha=0.75)
        ax.plot(beta, i, marker=marker, color=color, markersize=7 if is_sig else 5, zorder=5)

    ax.axvline(0, color="#555555", linewidth=1.5, linestyle="--")
    ax.set_yticks(range(n))
    ax.set_yticklabels(
        [f"{r.get('snp', '?')}  p={r.get('p_value', 1):.1e}" for r in plot_r],
        fontsize=8,
    )
    ax.set_xlim(x_min, x_max)
    ax.set_xlabel("Effect size (β) and 95% CI", fontsize=11)
    ax.set_title(
        f"Effect Size Forest Plot — Top {n} GWAS Hits\n"
        f"{sum(1 for r in sorted_r if r.get('p_value', 1) < significance_threshold)} "
        f"GW-significant (p<{significance_threshold:.0e}, ◆)",
        fontsize=11,
    )
    ax.grid(axis="x", alpha=style.grid_alpha if style.grid else 0.0)
    sig_patch = mpatches.Patch(color=style.significance_color, label=f"GW-sig (p<{significance_threshold:.0e})")
    nom_patch = mpatches.Patch(color=style.suggestive_color, label="Sub-threshold")
    ax.legend(handles=[sig_patch, nom_patch], fontsize=9, loc="lower right")
    plt.tight_layout()

    if output_path and fig is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=style.dpi, bbox_inches="tight")
        logger.info(f"Saved forest plot to {output_path}")
    return ax if fig is None else fig


def z_score_manhattan_plot(
    results: "List[Dict[str, Any]]",
    output_path: "Optional[Union[str, Path]]" = None,
    significance_threshold: float = 5e-8,
    top_n_label: int = 5,
) -> "Any":
    """Manhattan plot where point opacity encodes |z-score| = |β/SE|.

    This separates large-effect signals from high-precision (dense coverage)
    signals. GW-significant hits are annotated with their SNP ID.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for z-score Manhattan")
        return None
    if not results:
        return None

    chrom_order: "List[str]" = []
    seen: set = set()
    for r in results:
        c = r.get("chrom", "")
        if c not in seen:
            chrom_order.append(c)
            seen.add(c)

    chrom_max: "Dict[str, int]" = {c: max(r.get("pos", 0) for r in results if r.get("chrom") == c) for c in chrom_order}
    chrom_offset: "Dict[str, int]" = {}
    running = 0
    for c in chrom_order:
        chrom_offset[c] = running
        running += chrom_max.get(c, 0) + 5_000_000

    xs, ys, zs, chroms_list = [], [], [], []
    for r in results:
        p = max(r.get("p_value", 1.0), 1e-300)
        beta = r.get("beta", 0.0)
        se = r.get("se", 1.0)
        chrom = r.get("chrom", "")
        pos = r.get("pos", 0)
        xs.append(chrom_offset.get(chrom, 0) + pos)
        ys.append(-math.log10(p))
        zs.append(abs(beta / se) if se > 0 else 0.0)
        chroms_list.append(chrom)

    max_z = max(zs) if zs else 1.0
    sig_y = -math.log10(significance_threshold)

    palette = [
        "#2E86AB",
        "#A23B72",
        "#F18F01",
        "#C73E1D",
        "#44BBA4",
        "#E94F37",
        "#393E41",
        "#F4995C",
        "#D81E5B",
        "#5C6BC0",
        "#26A69A",
        "#AB47BC",
        "#7E57C2",
        "#42A5F5",
        "#6B4226",
        "#3B1F2B",
    ]
    fig, ax = plt.subplots(figsize=(14, 5))

    for ci, chrom in enumerate(chrom_order):
        base_color = palette[ci % len(palette)]
        idx_c = [i for i, c in enumerate(chroms_list) if c == chrom]
        for i in idx_c:
            z_norm = zs[i] / max_z
            alpha = 0.2 + 0.75 * z_norm
            sz = 18 if ys[i] >= sig_y else 6
            ax.scatter(xs[i], ys[i], c=base_color, s=sz, alpha=alpha, linewidths=0)

    ax.axhline(
        sig_y,
        color="#C44E52",
        linewidth=2,
        linestyle="--",
        label=f"GW significance (p={significance_threshold:.0e})",
    )

    # Annotate top |z| significant hits
    sig_pts = sorted(
        [(i, zs[i]) for i in range(len(xs)) if ys[i] >= sig_y],
        key=lambda t: -t[1],
    )[:top_n_label]
    for i, _ in sig_pts:
        snp_lbl = results[i].get("snp", "?") if i < len(results) else "?"
        ax.annotate(
            snp_lbl,
            (xs[i], ys[i]),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=7.5,
            color="#C44E52",
            fontweight="bold",
            arrowprops=dict(arrowstyle="-", color="#C44E52", lw=0.8),
        )

    for c in chrom_order:
        idx_c2 = [i for i, ch in enumerate(chroms_list) if ch == c]
        if idx_c2:
            mid_x = (min(xs[i] for i in idx_c2) + max(xs[i] for i in idx_c2)) / 2
            label = c.split(".")[-2][-2:] if "." in c else c
            ax.text(mid_x, -0.6, label, ha="center", va="top", fontsize=7, color="#555")

    ax.set_xlabel("Genomic position", fontsize=11)
    ax.set_ylabel("−log₁₀(p-value)", fontsize=11)
    ax.set_title(
        "Manhattan Plot — Colored by |z-score| (|β/SE|)\n" "Opacity encodes effect magnitude; brighter = larger |β/SE|",
        fontsize=11,
    )
    ax.set_ylim(bottom=-0.9)
    ax.legend(fontsize=9, loc="upper right")
    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved z-score Manhattan plot to {output_path}")
    return fig


def ld_heatmap_plot(
    ld_matrix: List[List[float]],
    association_results: List[Dict[str, Any]],
    *,
    top_idx: int = 0,
    half_window: int = 25,
    output_path: Optional[Union[str, Path]] = None,
) -> Optional[Any]:
    """Publication-grade LD structure heatmap centred on the top hit.

    Extracts a window of ±*half_window* variants around *top_idx* and
    renders |r| (absolute correlation) with genomic-position tick labels
    and a red crosshair marking the lead variant.

    Args:
        ld_matrix: Full variant-by-variant LD correlation matrix.
        association_results: List of result dicts with ``chrom`` and ``pos``.
        top_idx: Index of the lead variant in *association_results*.
        half_window: Half-size of the extraction window.
        output_path: If set, save the figure as PNG.

    Returns:
        matplotlib Figure, or ``None`` if dependencies are absent.
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("ld_heatmap_plot requires matplotlib and numpy")
        return None

    if not ld_matrix or len(ld_matrix) < 2:
        logger.info("Skipping LD heatmap (matrix too small or empty)")
        return None

    lo = max(0, top_idx - half_window)
    hi = min(len(ld_matrix), top_idx + half_window)
    ld_sub = np.array([row[lo:hi] for row in ld_matrix[lo:hi]])
    ld_sub = np.abs(ld_sub)  # Plot |r| not r²

    fig, ax = plt.subplots(figsize=(9, 7))
    im = ax.imshow(
        ld_sub,
        cmap="RdYlBu_r",
        vmin=0,
        vmax=1,
        aspect="auto",
        interpolation="nearest",
    )
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("|r| (correlation coefficient)", fontsize=10)

    # Genomic-position tick labels
    n_ticks = min(10, ld_sub.shape[0])
    tick_step = max(1, ld_sub.shape[0] // n_ticks)
    tick_positions = list(range(0, ld_sub.shape[0], tick_step))
    tick_labels = []
    for i in tick_positions:
        idx = lo + i
        if idx < len(association_results):
            r = association_results[idx]
            chrom_short = r["chrom"].split(".")[-2][-2:] if "." in r["chrom"] else r["chrom"]
            tick_labels.append(f"{chrom_short}:{r['pos'] // 1000}k")
        else:
            tick_labels.append("")
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(tick_labels, fontsize=7)

    # Title with top-hit coordinates
    if top_idx < len(association_results):
        top_r = association_results[top_idx]
        ax.set_title(
            f"LD Structure Near Top Hit ({top_r['chrom']}:{top_r['pos']:,})\n"
            f"Window: ±{half_window} variants, n={ld_sub.shape[0]}",
            fontsize=11,
        )

    # Red crosshair on top hit
    rel_top = min(half_window, top_idx - lo)
    ax.axhline(rel_top, color="red", lw=1.5, ls="--", alpha=0.7, label="Top hit")
    ax.axvline(rel_top, color="red", lw=1.5, ls="--", alpha=0.7)
    ax.legend(fontsize=8, loc="upper right")

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved LD heatmap to {output_path} ({ld_sub.shape[0]}×{ld_sub.shape[1]})")
    return fig
