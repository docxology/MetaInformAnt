"""Genome-wide visualization functions for GWAS.

This module provides specialized plots for genomic data visualization.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def circular_manhattan_plot(
    results_df: Any,
    output_file: Optional[str | Path] = None,
    significance_threshold: float = 5e-8,
    title: str = "Circular Manhattan Plot",
    suggestive_threshold: Optional[float] = 1e-5,
    label_top_n: int = 0,
) -> Dict[str, Any]:
    """Create a circular Manhattan plot for GWAS results.

    Args:
        results_df: DataFrame with columns 'CHR', 'BP', 'P', or list of dicts
            with 'CHROM', 'POS', 'p_value' keys
        output_file: Optional output file path
        significance_threshold: P-value threshold for significance line
        title: Plot title
        suggestive_threshold: Optional p-value threshold for a secondary dashed
            ring. Defaults to 1e-5. Set to None to disable.
        label_top_n: Auto-label the top N hits with their chr:pos. Defaults to
            0 (disabled).

    Returns:
        Dictionary with status and metadata

    Example:
        >>> result = circular_manhattan_plot(gwas_results)
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Wedge
    except ImportError:
        logger.warning("matplotlib not available for circular Manhattan plot")
        return {"status": "failed", "error": "matplotlib not available"}

    # Convert list-of-dicts to parallel arrays
    if isinstance(results_df, list):
        if not results_df:
            return {"status": "failed", "error": "No results provided"}
        chroms_list = []
        positions_list = []
        pvalues_list = []
        for entry in results_df:
            chroms_list.append(entry.get("CHROM", entry.get("CHR", entry.get("chr", "1"))))
            positions_list.append(entry.get("POS", entry.get("BP", entry.get("pos", 0))))
            pvalues_list.append(entry.get("p_value", entry.get("P", entry.get("pval", 1.0))))
    elif hasattr(results_df, "columns"):
        # DataFrame path
        required_cols = ["CHR", "BP", "P"]
        missing_cols = [col for col in required_cols if col not in results_df.columns]
        if missing_cols:
            logger.error(f"Missing required columns: {missing_cols}")
            return {"status": "failed", "error": f"Missing required columns: {missing_cols}"}
        chroms_list = results_df["CHR"].tolist()
        positions_list = results_df["BP"].tolist()
        pvalues_list = results_df["P"].tolist()
    else:
        logger.error("Input data must be a DataFrame or list of dicts")
        return {"status": "failed", "error": "Input data must be a DataFrame or list of dicts"}

    fig, ax = plt.subplots(1, 1, figsize=(10, 10), subplot_kw=dict(projection="polar"))

    # Prepare data
    unique_chroms = sorted(set(str(c) for c in chroms_list))
    chrom_angles: Dict[str, tuple] = {}

    # Calculate chromosome angles (equal spacing)
    angle_per_chrom = 360 / len(unique_chroms)

    # Colors for chromosomes
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_chroms)))

    # Group data by chromosome
    chrom_data_map: Dict[str, List[tuple]] = {c: [] for c in unique_chroms}
    for chrom, pos, pval in zip(chroms_list, positions_list, pvalues_list):
        chrom_data_map[str(chrom)].append((pos, pval))

    # Plot each chromosome as a sector
    for i, chrom in enumerate(unique_chroms):
        data = chrom_data_map[chrom]
        if not data:
            continue

        start_angle = i * angle_per_chrom
        end_angle = (i + 1) * angle_per_chrom
        chrom_angles[chrom] = (start_angle, end_angle)

        pos_arr = np.array([d[0] for d in data], dtype=float)
        pval_arr = np.array([d[1] for d in data], dtype=float)

        # Convert p-values to radial positions
        neg_log_p = -np.log10(np.clip(pval_arr, 1e-50, None))

        # Normalize positions within chromosome
        pos_range = pos_arr.max() - pos_arr.min()
        if pos_range > 0:
            norm_pos = (pos_arr - pos_arr.min()) / pos_range
        else:
            norm_pos = np.full_like(pos_arr, 0.5)

        # Convert to polar coordinates
        angles = start_angle + norm_pos * angle_per_chrom
        radii = 1 + neg_log_p * 0.1  # Scale for visibility

        # Plot points
        ax.scatter(np.deg2rad(angles), radii, c=[colors[i]], s=1, alpha=0.6)

        # Add chromosome label
        label_angle = np.deg2rad(start_angle + angle_per_chrom / 2)
        label_radius = 1.2
        ax.text(label_angle, label_radius, str(chrom), ha="center", va="center", fontsize=10, fontweight="bold")

    # Add significance threshold ring
    threshold_radius = 1 + (-np.log10(significance_threshold)) * 0.1
    theta = np.linspace(0, 2 * np.pi, 100)
    ax.plot(theta, [threshold_radius] * 100, "r--", linewidth=2, alpha=0.8)

    # Add significance label
    ax.text(
        0, threshold_radius + 0.05, f"p = {significance_threshold}", ha="center", va="bottom", fontsize=10, color="red"
    )

    # Add suggestive threshold ring
    if suggestive_threshold is not None and suggestive_threshold > 0:
        sug_radius = 1 + (-np.log10(suggestive_threshold)) * 0.1
        ax.plot(theta, [sug_radius] * 100, "b--", linewidth=1.5, alpha=0.6)
        ax.text(
            np.pi / 4,
            sug_radius + 0.05,
            f"p = {suggestive_threshold}",
            ha="center",
            va="bottom",
            fontsize=8,
            color="blue",
        )

    # Label top N hits
    if label_top_n > 0:
        # Collect all points with their polar coordinates
        all_points: List[tuple] = []
        for i, chrom in enumerate(unique_chroms):
            data = chrom_data_map[chrom]
            if not data:
                continue
            start_angle = i * angle_per_chrom
            pos_arr = np.array([d[0] for d in data], dtype=float)
            pval_arr = np.array([d[1] for d in data], dtype=float)
            neg_log_p_arr = -np.log10(np.clip(pval_arr, 1e-50, None))
            pos_range = pos_arr.max() - pos_arr.min()
            if pos_range > 0:
                norm_pos = (pos_arr - pos_arr.min()) / pos_range
            else:
                norm_pos = np.full_like(pos_arr, 0.5)
            for j in range(len(data)):
                angle_deg = start_angle + norm_pos[j] * angle_per_chrom
                radius = 1 + neg_log_p_arr[j] * 0.1
                all_points.append((pval_arr[j], angle_deg, radius, chrom, int(pos_arr[j])))

        # Sort by p-value ascending (smallest p-value = most significant)
        all_points.sort(key=lambda x: x[0])
        for rank, (pval, angle_deg, radius, chrom_label, pos_val) in enumerate(all_points[:label_top_n]):
            label_text = f"{chrom_label}:{pos_val}"
            angle_rad = np.deg2rad(angle_deg)
            label_radius = radius + 0.15 + rank * 0.05
            ax.annotate(
                label_text,
                xy=(angle_rad, radius),
                xytext=(angle_rad, label_radius),
                arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
                fontsize=7,
                ha="center",
                color="darkred",
            )

    # Configure plot
    ax.set_title(title, fontsize=14, pad=20)
    ax.set_rlabel_position(0)
    ax.set_rticks([])
    ax.set_thetagrids([])
    ax.grid(False)
    ax.spines["polar"].set_visible(False)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved circular Manhattan plot to {output_file}")

    return {
        "status": "success",
        "n_variants": len(chroms_list),
        "n_chromosomes": len(unique_chroms),
        "output_path": str(output_file) if output_file else None,
    }


def chromosome_ideogram(
    chromosome_data: Any,
    output_file: Optional[str | Path] = None,
    highlighted_regions: Optional[List[Dict[str, Any]]] = None,
    title: str = "Chromosome Ideogram",
) -> Dict[str, Any]:
    """Create a chromosome ideogram showing chromosome structure.

    Args:
        chromosome_data: Dictionary mapping chromosome names to lengths,
            or a list of GWAS result dicts with 'CHROM' and 'POS' keys
            (chromosome lengths will be inferred from max positions).
        output_file: Optional output file path
        highlighted_regions: Optional list of regions to highlight
        title: Plot title

    Returns:
        Dictionary with status and metadata
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.patches import Rectangle
    except ImportError:
        logger.warning("matplotlib not available for chromosome ideogram")
        return {"status": "failed", "error": "matplotlib not available"}

    # Convert list-of-dicts (GWAS results) to chromosome_lengths dict
    if isinstance(chromosome_data, list):
        if not chromosome_data:
            return {"status": "failed", "error": "No chromosome data provided"}
        chromosome_lengths: Dict[str, int] = {}
        for entry in chromosome_data:
            chrom = str(entry.get("CHROM", entry.get("CHR", entry.get("chr", "1"))))
            pos = int(entry.get("POS", entry.get("BP", entry.get("pos", 0))))
            if chrom not in chromosome_lengths or pos > chromosome_lengths[chrom]:
                chromosome_lengths[chrom] = pos
    elif isinstance(chromosome_data, dict):
        chromosome_lengths = chromosome_data
    else:
        return {"status": "failed", "error": "Invalid chromosome data format"}

    if not chromosome_lengths:
        logger.error("No chromosome lengths provided")
        return {"status": "failed", "error": "No chromosome lengths provided"}

    # Sort chromosomes
    chroms = sorted(
        chromosome_lengths.keys(),
        key=lambda x: int(x.replace("chr", "")) if x.replace("chr", "").isdigit() else float("inf"),
    )

    fig, ax = plt.subplots(figsize=(12, 8))

    y_start = 0.1
    y_height = 0.8 / len(chroms) if chroms else 0.8
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

    max_length = max(chromosome_lengths.values()) if chromosome_lengths else 1

    for i, chrom in enumerate(chroms):
        length = chromosome_lengths[chrom]
        color = colors[i % len(colors)]

        # Draw chromosome body
        width = length / max_length * 0.8  # Scale to 80% of plot width
        rect = Rectangle(
            (0.1, y_start + i * y_height),
            width,
            y_height * 0.8,
            facecolor=color,
            alpha=0.7,
            edgecolor="black",
            linewidth=1,
        )
        ax.add_patch(rect)

        # Add chromosome label
        ax.text(
            0.05,
            y_start + i * y_height + y_height * 0.4,
            chrom,
            ha="right",
            va="center",
            fontsize=10,
            fontweight="bold",
        )

        # Add length label
        length_mb = length / 1e6
        ax.text(
            width + 0.12,
            y_start + i * y_height + y_height * 0.4,
            f"{length_mb:.1f} Mb",
            ha="left",
            va="center",
            fontsize=8,
        )

        # Highlight regions if provided
        if highlighted_regions:
            for region in highlighted_regions:
                if region.get("chrom") == chrom:
                    region_start = region.get("start", 0)
                    region_end = region.get("end", 0)
                    if region_end > region_start:
                        # Convert to plot coordinates
                        x_start = 0.1 + (region_start / max_length) * 0.8
                        x_width = ((region_end - region_start) / max_length) * 0.8

                        highlight_color = region.get("color", "red")
                        alpha = region.get("alpha", 0.5)

                        highlight_rect = Rectangle(
                            (x_start, y_start + i * y_height),
                            x_width,
                            y_height * 0.8,
                            facecolor=highlight_color,
                            alpha=alpha,
                            edgecolor="darkred",
                            linewidth=1,
                        )
                        ax.add_patch(highlight_rect)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=14, pad=20)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved chromosome ideogram to {output_file}")

    return {
        "status": "success",
        "n_chromosomes": len(chroms),
        "chromosomes": chroms,
        "output_path": str(output_file) if output_file else None,
    }


def genome_wide_ld_heatmap(
    ld_data: List[Dict[str, Any]],
    output_file: Optional[str | Path] = None,
    figsize: tuple[int, int] = (12, 10),
    title: str = "Genome-wide LD Heatmap",
    chromosomes: Optional[List[str]] = None,
    ld_threshold: float = 0.8,
) -> Dict[str, Any]:
    """Create a genome-wide LD heatmap visualization.

    Args:
        ld_data: List of LD data dictionaries with CHROM, POS1, POS2, r2 keys
        output_file: Optional output file path
        figsize: Figure size (width, height)
        title: Plot title
        chromosomes: List of chromosomes to include (None for all)
        ld_threshold: LD threshold for highlighting high LD regions

    Returns:
        Dictionary with status and metadata
    """
    try:
        import matplotlib.patches as patches
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for genome-wide LD heatmap")
        return {"status": "error", "message": "matplotlib not available"}

    if not ld_data:
        logger.warning("No LD data provided")
        return {"status": "error", "message": "No LD data provided"}

    # Group data by chromosome
    chrom_data = {}
    for entry in ld_data:
        chrom = entry["CHROM"]
        if chromosomes and chrom not in chromosomes:
            continue
        if chrom not in chrom_data:
            chrom_data[chrom] = []
        chrom_data[chrom].append(entry)

    if not chrom_data:
        logger.warning("No data found for specified chromosomes")
        return {"status": "error", "message": "No data for specified chromosomes"}

    # Create figure with subplots for each chromosome
    n_chroms = len(chrom_data)
    if n_chroms <= 3:
        nrows, ncols = 1, n_chroms
    elif n_chroms <= 6:
        nrows, ncols = 2, 3
    else:
        nrows = int(np.ceil(n_chroms / 4))
        ncols = 4

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False, sharex=True, sharey=True)
    axes = axes.flatten()

    # Plot each chromosome
    for i, (chrom, data) in enumerate(sorted(chrom_data.items())):
        if i >= len(axes):
            break

        ax = axes[i]

        # Extract positions and LD values
        pos1 = [entry["POS1"] for entry in data]
        pos2 = [entry["POS2"] for entry in data]
        r2_values = [entry["r2"] for entry in data]

        # Create scatter plot colored by LD
        scatter = ax.scatter(pos1, pos2, c=r2_values, cmap="RdYlBu_r", s=2, alpha=0.6, vmin=0, vmax=1)

        # Highlight high LD regions
        high_ld = [entry for entry in data if entry["r2"] >= ld_threshold]
        if high_ld:
            ax.scatter(
                [entry["POS1"] for entry in high_ld],
                [entry["POS2"] for entry in high_ld],
                c="red",
                s=5,
                alpha=0.8,
                marker="x",
                label=f"r² ≥ {ld_threshold}",
            )

        ax.set_xlabel("Position 1 (bp)")
        ax.set_ylabel("Position 2 (bp)")
        ax.set_title(f"{chrom}")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right", fontsize=8)

        # Format axis labels
        def format_bp(x, pos):
            if x >= 1e6:
                return f"{x/1e6:.1f}M"
            elif x >= 1e3:
                return f"{x/1e3:.0f}K"
            else:
                return f"{x:.0f}"

        ax.xaxis.set_major_formatter(plt.FuncFormatter(format_bp))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(format_bp))

    # Hide unused subplots
    for i in range(len(chrom_data), len(axes)):
        axes[i].set_visible(False)

    # Add colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(scatter, cax=cbar_ax)
    cbar.set_label("LD (r²)")

    fig.suptitle(title, fontsize=14)

    plt.tight_layout(rect=[0, 0, 0.9, 0.95])

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved genome-wide LD heatmap to {output_file}")

    return {
        "status": "success",
        "n_chromosomes": len(chrom_data),
        "total_ld_pairs": len(ld_data),
        "chromosomes": list(chrom_data.keys()),
        "ld_threshold": ld_threshold,
    }


def manhattan_plot(
    results: list[dict],
    output_path: Optional[str | Path] = None,
    significance_threshold: float = 5e-8,
    suggestiveness_threshold: Optional[float] = None,
    figsize: tuple[int, int] = (12, 8),
) -> dict[str, Any]:
    """Create a Manhattan plot for GWAS results.

    Args:
        results: List of dictionaries with CHROM, POS, p_value keys
        output_path: Path to save the plot
        significance_threshold: P-value threshold for significance line
        suggestiveness_threshold: P-value threshold for suggestive significance line
        figsize: Figure size

    Returns:
        Dictionary with plot status and metadata
    """
    try:
        import matplotlib.pyplot as plt

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    if not HAS_MATPLOTLIB:
        return {"status": "failed", "error": "matplotlib not available"}

    if not results:
        return {"status": "failed", "error": "No results provided"}

    try:
        # Convert to DataFrame-like structure
        chromosomes = []
        positions = []
        p_values = []

        for result in results:
            chromosomes.append(result.get("CHROM", result.get("chr", 1)))
            positions.append(result.get("POS", result.get("pos", 0)))
            p_values.append(result.get("p_value", result.get("P", 1.0)))

        # Create chromosome mapping
        unique_chroms = sorted(set(str(chrom) for chrom in chromosomes))
        chrom_to_num = {chrom: i + 1 for i, chrom in enumerate(unique_chroms)}

        # Convert chromosome names to numbers
        chrom_nums = [chrom_to_num[str(chrom)] for chrom in chromosomes]

        # Calculate x positions
        x_positions = []
        chrom_max_pos = {}
        current_pos = 0

        for i, (chrom, pos) in enumerate(zip(chrom_nums, positions)):
            if chrom not in chrom_max_pos:
                chrom_max_pos[chrom] = current_pos
            x_positions.append(current_pos + pos)
            current_pos = max(current_pos, chrom_max_pos[chrom] + pos + 1)

        # Convert p-values to -log10
        neg_log_p = [-math.log10(max(p, 1e-300)) for p in p_values]

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        # Color by chromosome
        colors = ["#1f77b4", "#ff7f0e"] * len(unique_chroms)
        num_to_chrom = {v: k for k, v in chrom_to_num.items()}
        chrom_colors = {}
        for i, chrom in enumerate(unique_chroms):
            chrom_colors[chrom] = colors[i % len(colors)]

        for i in range(len(x_positions)):
            chrom_name = num_to_chrom[chrom_nums[i]]
            ax.scatter(x_positions[i], neg_log_p[i], color=chrom_colors[chrom_name], alpha=0.8, s=20)

        # Add significance line
        sig_threshold = -math.log10(significance_threshold)
        ax.axhline(
            y=sig_threshold,
            color="red",
            linestyle="--",
            alpha=0.7,
            label=f"Significance threshold ({significance_threshold})",
        )

        # Add suggestiveness threshold line if provided
        if suggestiveness_threshold is not None:
            sug_threshold = -math.log10(suggestiveness_threshold)
            ax.axhline(
                y=sug_threshold,
                color="blue",
                linestyle=":",
                alpha=0.7,
                label=f"Suggestiveness threshold ({suggestiveness_threshold})",
            )

        # Add chromosome labels at centers
        chrom_centers = []
        for chrom in unique_chroms:
            chrom_num = chrom_to_num[chrom]
            chrom_x = [x for x, c in zip(x_positions, chrom_nums) if c == chrom_num]
            if chrom_x:
                chrom_centers.append((chrom, sum(chrom_x) / len(chrom_x)))

        for chrom, center in chrom_centers:
            ax.text(center, min(neg_log_p) - 0.5, str(chrom), ha="center", va="top", fontsize=8)

        ax.set_xlabel("Genomic Position")
        ax.set_ylabel("-log₁₀(p-value)")
        ax.set_title("Manhattan Plot")
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")

        return {
            "status": "success",
            "n_variants": len(results),
            "n_chromosomes": len(unique_chroms),
            "significance_threshold": significance_threshold,
            "chromosomes": unique_chroms,
            "output_path": str(output_path) if output_path else None,
        }

    except Exception as e:
        return {"status": "failed", "error": str(e)}
