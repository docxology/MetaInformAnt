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
        plt.rcParams.update({
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
        })

if HAS_MATPLOTLIB:
    set_accessible_style()


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

    # Convert results to list of dictionaries if it's a single dict
    if isinstance(results, dict):
        results = [results]

    # Extract data for plotting
    chromosomes = []
    positions = []
    p_values = []
    colors = []

    # Load style if not provided
    if style is None:
        from metainformant.gwas.visualization.config import get_style
        style = get_style()

    # Create plot
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get("figsize", style.figsize))

    # Color scheme for chromosomes from style
    chrom_colors = style.categorical_colors

    current_pos = 0
    chrom_offsets = {}

    for result in results:
        if not isinstance(result, dict):
            continue

        chrom = str(result.get("chrom", result.get("chromosome", "1")))
        pos = result.get("pos", result.get("position", 0))
        p_val = result.get("p_value", result.get("pval", 1.0))

        # Convert p-value to -log10 scale
        if isinstance(p_val, (int, float)) and p_val > 0:
            neg_log_p = -math.log10(p_val)
        else:
            neg_log_p = 50  # Cap very small p-values or invalid ones

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

    # Plot points with kwargs
    scatter_kwargs = kwargs.get("scatter_kwargs", {"s": style.point_size, "alpha": style.alpha})
    scatter = ax.scatter(positions, p_values, c=colors, **scatter_kwargs)

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
            result_entry = results[idx] if idx < len(results) else {}
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
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log₁₀(p-value)")
    ax.set_title("Manhattan Plot")
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

    # Save if output path provided and we own the figure
    if output_path and fig is not None:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=style.dpi, bbox_inches="tight")
        logger.info(f"Saved Manhattan plot to {output_path}")

    return ax if fig is None else fig

def qq_plot(
    p_values: Union[List[float], List[int]],
    output_path: Optional[Union[str, Path]] = None,
    ax: Optional[Any] = None,
    style: Optional[Any] = None,
    **kwargs: Any,
) -> Any:
    """Create a Q-Q plot from p-values.

    Args:
        p_values: List of p-values
        output_path: Path to save the plot (optional)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create Q-Q plot")
        return None

    logger.info("Creating Q-Q plot")

    # Convert to numpy array and filter out invalid values
    p_vals = np.array(p_values, dtype=float)
    p_vals = p_vals[~np.isnan(p_vals) & (p_vals > 0) & (p_vals <= 1)]

    if len(p_vals) == 0:
        logger.warning("No valid p-values for Q-Q plot")
        return None

    # Sort p-values
    p_vals = np.sort(p_vals)

    # Expected p-values under null hypothesis
    n = len(p_vals)
    expected = np.arange(1, n + 1) / (n + 1)
    expected_log = -np.log10(expected)

    # Observed -log10 p-values
    observed_log = -np.log10(p_vals)

    if style is None:
        from metainformant.gwas.visualization.config import get_style
        style = get_style()

    # Create plot
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get("figsize", (8, 6)))

    # Plot Q-Q points
    scatter_kwargs = kwargs.get("scatter_kwargs", {"s": max(2, style.point_size - 4), "alpha": style.alpha, "color": style.suggestive_color})
    ax.scatter(expected_log, observed_log, **scatter_kwargs)

    # Add diagonal line (expected under null)
    max_val = max(np.max(expected_log), np.max(observed_log))
    ax.plot([0, max_val], [0, max_val], linestyle="--", color=style.significance_color, alpha=0.7, label="Expected (null)")

    # Labels and title
    ax.set_xlabel("Expected -log₁₀(p-value)")
    ax.set_ylabel("Observed -log₁₀(p-value)")
    ax.set_title(f"Q-Q Plot (n={n} variants)")
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    if output_path and fig is not None:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=style.dpi, bbox_inches="tight")
        logger.info(f"Saved Q-Q plot to {output_path}")

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

    if style is None:
        from metainformant.gwas.visualization.config import get_style
        style = get_style()

    # Create plot
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get("figsize", style.figsize))

    # Plot points
    scatter_kwargs = kwargs.get("scatter_kwargs", {"s": style.point_size, "alpha": style.alpha, "color": style.suggestive_color})
    ax.scatter(positions, p_values, **scatter_kwargs)

    # Add significance threshold line
    threshold_line = -math.log10(5e-8)
    ax.axhline(y=threshold_line, color=style.significance_color, linestyle="--", alpha=0.7, label="Genome-wide significance")

    # Labels and title
    ax.set_xlabel(f"Position on chromosome {chrom}")
    ax.set_ylabel("-log₁₀(p-value)")
    ax.set_title(f"Regional Association Plot: {chrom}:{start:,}-{end:,}")
    ax.set_xlim(start, end)
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    if output_path and fig is not None:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=style.dpi, bbox_inches="tight")
        logger.info(f"Saved regional plot to {output_path}")

    return ax if fig is None else fig


def pca_plot(
    pca_result: tuple,
    output_path: Optional[Union[str, Path]] = None,
    explained_var: Optional[List[float]] = None,
    ax: Optional[Any] = None,
    style: Optional[Any] = None,
    **kwargs: Any,
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

        if style is None:
            from metainformant.gwas.visualization.config import get_style
            style = get_style()

        # Create 2D scatter plot of first two components
        fig = None
        if ax is None:
            fig, ax = plt.subplots(figsize=kwargs.get("figsize", (10, 8)))

        # Plot points
        scatter_kwargs = kwargs.get("scatter_kwargs", {"s": max(2, style.point_size - 4), "alpha": style.alpha, "color": style.suggestive_color})
        scatter = ax.scatter(components[0], components[1], **scatter_kwargs)

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
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.9, edgecolor="black"),
            )

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
    ax: Optional[Any] = None,
    style: Optional[Any] = None,
    **kwargs: Any,
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

        if style is None:
            from metainformant.gwas.visualization.config import get_style
            style = get_style()

        # Create heatmap
        fig = None
        if ax is None:
            fig, ax = plt.subplots(figsize=kwargs.get("figsize", (10, 8)))

        # Plot heatmap
        imshow_kwargs = kwargs.get("imshow_kwargs", {"cmap": style.colormap, "aspect": "equal"})
        im = ax.imshow(kinship_matrix, **imshow_kwargs)

        # Add colorbar only if we created the figure (or explicitly asked)
        if fig is not None or kwargs.get("show_colorbar", False):
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label("Kinship coefficient")

        # Labels and title
        ax.set_title("Kinship Matrix Heatmap")
        ax.set_xlabel("Sample")
        ax.set_ylabel("Sample")

        plt.tight_layout()

        if output_path and fig is not None:
            output_path = Path(output_path)
            fig.savefig(output_path, dpi=style.dpi, bbox_inches="tight")
            logger.info(f"Saved kinship heatmap to {output_path}")

        return ax if fig is None else fig

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

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(power_data))) if HAS_NUMPY else [
        f"C{i}" for i in range(len(power_data))
    ]

    for idx, curve in enumerate(power_data):
        sizes = curve.get("sample_sizes", [])
        powers = curve.get("powers", [])
        beta = curve.get("beta", 0)
        maf = curve.get("maf", 0)
        ax.plot(
            sizes, powers,
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
    ax.plot(pct_fractions, predicted, color="#E74C3C", linewidth=2, linestyle="--", label="Fitted curve")
    ax.axhline(y=k_inf, color="#2ECC71", linestyle=":", linewidth=1.5, label=f"K_∞ = {k_inf:.3f}")

    status_text = "✓ SATURATED" if is_saturated else "✗ NOT SATURATED"
    status_color = "#2ECC71" if is_saturated else "#E74C3C"
    ax.text(
        0.02, 0.95, status_text,
        transform=ax.transAxes, fontsize=12, fontweight="bold",
        color=status_color, va="top",
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
        [f"{r.get('snp','?')}  p={r.get('p_value',1):.1e}"
         for r in plot_r],
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

    chrom_max: "Dict[str, int]" = {
        c: max(r.get("pos", 0) for r in results if r.get("chrom") == c)
        for c in chrom_order
    }
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
        "#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#44BBA4",
        "#E94F37", "#393E41", "#F4995C", "#D81E5B", "#5C6BC0",
        "#26A69A", "#AB47BC", "#7E57C2", "#42A5F5", "#6B4226", "#3B1F2B",
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

    ax.axhline(sig_y, color="#C44E52", linewidth=2, linestyle="--",
               label=f"GW significance (p={significance_threshold:.0e})")

    # Annotate top |z| significant hits
    sig_pts = sorted(
        [(i, zs[i]) for i in range(len(xs)) if ys[i] >= sig_y],
        key=lambda t: -t[1],
    )[:top_n_label]
    for i, _ in sig_pts:
        snp_lbl = results[i].get("snp", "?") if i < len(results) else "?"
        ax.annotate(
            snp_lbl, (xs[i], ys[i]),
            xytext=(0, 9), textcoords="offset points",
            ha="center", fontsize=7.5, color="#C44E52", fontweight="bold",
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
        "Manhattan Plot — Colored by |z-score| (|β/SE|)\n"
        "Opacity encodes effect magnitude; brighter = larger |β/SE|",
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
