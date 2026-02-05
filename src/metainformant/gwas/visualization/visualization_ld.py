"""Linkage disequilibrium (LD) visualization for GWAS.

This module provides functions for computing and visualizing LD decay curves
and regional LD heatmaps (Haploview-style triangular plots).
"""

from __future__ import annotations

import math
import random
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    plt = None  # type: ignore

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore


def _compute_r_squared(geno_a: List[int], geno_b: List[int]) -> float:
    """Compute R-squared (squared Pearson correlation) between two genotype vectors.

    Delegates to the shared implementation in analysis.utils to avoid duplication.

    Args:
        geno_a: Genotype values for variant A (0, 1, 2 or -1 for missing).
        geno_b: Genotype values for variant B.

    Returns:
        R-squared value in [0, 1].
    """
    from metainformant.gwas.analysis.utils import compute_r_squared

    return compute_r_squared(geno_a, geno_b)


def compute_ld_decay(
    genotypes_by_variant: List[List[int]],
    positions: List[int],
    max_distance: int = 500000,
    n_bins: int = 50,
    chromosomes: Optional[List[int]] = None,
) -> Dict[str, Any]:
    """Compute mean r-squared in distance bins across variant pairs.

    For each pair of variants within max_distance, computes r-squared and bins
    the result by physical distance. Returns binned mean r-squared values
    suitable for plotting an LD decay curve.

    Args:
        genotypes_by_variant: Genotype matrix (variants x samples), values 0/1/2/-1.
        positions: Physical positions for each variant (bp).
        max_distance: Maximum distance (bp) between variant pairs to consider.
        n_bins: Number of distance bins.
        chromosomes: Optional chromosome assignments per variant. When provided,
            only pairs on the same chromosome are evaluated.

    Returns:
        Dictionary with status, bin_centers, mean_r2, n_pairs per bin, and total_pairs.
    """
    n_variants = len(genotypes_by_variant)
    if n_variants < 2:
        logger.warning(f"Need at least 2 variants for LD decay, got {n_variants}")
        return {
            "status": "failed",
            "error": f"Need at least 2 variants, got {n_variants}",
            "bin_centers": [],
            "mean_r2": [],
            "n_pairs": [],
            "total_pairs": 0,
        }

    if len(positions) != n_variants:
        return {
            "status": "failed",
            "error": "Length of positions must match number of variants",
            "bin_centers": [],
            "mean_r2": [],
            "n_pairs": [],
            "total_pairs": 0,
        }

    logger.info(f"Computing LD decay: {n_variants} variants, max_distance={max_distance}, " f"n_bins={n_bins}")

    # Build bin edges
    bin_width = max_distance / n_bins
    bin_edges = [i * bin_width for i in range(n_bins + 1)]
    bin_sums = [0.0] * n_bins
    bin_counts = [0] * n_bins

    # Enumerate all valid pairs; subsample if too many potential pairs
    max_pairs = 100000
    # Build candidate pair list
    candidate_pairs: List[Tuple[int, int]] = []
    for i in range(n_variants):
        for j in range(i + 1, n_variants):
            # Skip pairs on different chromosomes
            if chromosomes is not None and chromosomes[i] != chromosomes[j]:
                continue
            dist = abs(positions[j] - positions[i])
            if dist <= max_distance:
                candidate_pairs.append((i, j))

    # Subsample if needed
    if len(candidate_pairs) > max_pairs:
        logger.info(f"Subsampling {max_pairs} pairs from {len(candidate_pairs)} candidates")
        random.seed(42)
        candidate_pairs = random.sample(candidate_pairs, max_pairs)

    total_pairs = 0
    for i, j in candidate_pairs:
        dist = abs(positions[j] - positions[i])
        r2 = _compute_r_squared(genotypes_by_variant[i], genotypes_by_variant[j])

        # Find the appropriate bin
        bin_idx = int(dist / bin_width)
        if bin_idx >= n_bins:
            bin_idx = n_bins - 1

        bin_sums[bin_idx] += r2
        bin_counts[bin_idx] += 1
        total_pairs += 1

    # Compute bin centers and mean r2
    bin_centers = [(bin_edges[k] + bin_edges[k + 1]) / 2.0 for k in range(n_bins)]
    mean_r2 = [bin_sums[k] / bin_counts[k] if bin_counts[k] > 0 else 0.0 for k in range(n_bins)]

    logger.info(f"LD decay computed: {total_pairs} pairs across {n_bins} bins")

    return {
        "status": "success",
        "bin_centers": bin_centers,
        "mean_r2": mean_r2,
        "n_pairs": bin_counts,
        "total_pairs": total_pairs,
    }


def ld_decay_plot(
    ld_decay_data: Dict[str, Any],
    output_file: Optional[Union[str, Path]] = None,
    fit_curve: bool = True,
    title: str = "LD Decay",
) -> Dict[str, Any]:
    """Plot an LD decay curve: mean r-squared vs physical distance.

    Optionally fits an exponential decay model r2 = a * exp(-b * d) + c
    and annotates the distance at which the fitted curve crosses r2 = 0.2.

    Args:
        ld_decay_data: Output from compute_ld_decay containing bin_centers,
            mean_r2, and n_pairs.
        output_file: Optional path to save the plot image.
        fit_curve: Whether to fit an exponential decay curve.
        title: Plot title.

    Returns:
        Dictionary with status and output_path.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for LD decay plot")
        return {"status": "skipped", "output_path": None}

    if not HAS_NUMPY:
        logger.warning("numpy not available for LD decay plot")
        return {"status": "skipped", "output_path": None}

    try:
        bin_centers = ld_decay_data.get("bin_centers", [])
        mean_r2 = ld_decay_data.get("mean_r2", [])
        n_pairs = ld_decay_data.get("n_pairs", [])

        if not bin_centers or not mean_r2:
            return {"status": "failed", "error": "No LD decay data to plot", "output_path": None}

        # Convert to kb for x-axis
        bin_centers_kb = [c / 1000.0 for c in bin_centers]

        # Filter bins with data
        valid_mask = [n_pairs[i] > 0 for i in range(len(n_pairs))]
        plot_x = [bin_centers_kb[i] for i in range(len(bin_centers_kb)) if valid_mask[i]]
        plot_y = [mean_r2[i] for i in range(len(mean_r2)) if valid_mask[i]]

        if not plot_x:
            return {"status": "failed", "error": "No bins with data", "output_path": None}

        fig, ax = plt.subplots(figsize=(10, 6))

        # Scatter plot of binned means
        ax.scatter(plot_x, plot_y, c="steelblue", alpha=0.7, s=30, label="Mean r\u00b2")

        # Fit exponential decay if requested
        ld_decay_distance = None
        if fit_curve and len(plot_x) >= 3:
            try:
                x_arr = np.array(plot_x)
                y_arr = np.array(plot_y)

                # Fit r2 = a * exp(-b * d) + c via grid search
                # Estimate initial parameters
                best_a = max(y_arr) - min(y_arr)
                best_c = min(y_arr)
                best_b = 0.01
                best_sse = float("inf")

                # Grid search over parameter space
                a_range = np.linspace(0.01, max(y_arr), 10)
                b_range = np.logspace(-4, -1, 15)
                c_range = np.linspace(0.0, max(min(y_arr), 0.01), 5)

                for a_try in a_range:
                    for b_try in b_range:
                        for c_try in c_range:
                            predicted = a_try * np.exp(-b_try * x_arr) + c_try
                            sse = float(np.sum((y_arr - predicted) ** 2))
                            if sse < best_sse:
                                best_sse = sse
                                best_a = a_try
                                best_b = b_try
                                best_c = c_try

                # Plot fitted curve
                x_fit = np.linspace(min(plot_x), max(plot_x), 200)
                y_fit = best_a * np.exp(-best_b * x_fit) + best_c
                ax.plot(x_fit, y_fit, "r-", linewidth=2, label="Fitted decay")

                # Find where curve crosses r2 = 0.2
                threshold = 0.2
                if best_a + best_c > threshold and best_c < threshold:
                    # Solve a * exp(-b * d) + c = 0.2 => d = -ln((0.2 - c) / a) / b
                    ratio = (threshold - best_c) / best_a
                    if ratio > 0:
                        ld_decay_distance = -math.log(ratio) / best_b
                        if min(plot_x) <= ld_decay_distance <= max(plot_x):
                            ax.axvline(
                                x=ld_decay_distance,
                                color="green",
                                linestyle="--",
                                alpha=0.7,
                            )
                            ax.annotate(
                                f"LD decay: {ld_decay_distance:.1f} kb",
                                xy=(ld_decay_distance, threshold),
                                xytext=(ld_decay_distance + max(plot_x) * 0.05, threshold + 0.05),
                                fontsize=9,
                                arrowprops=dict(arrowstyle="->", color="green"),
                                color="green",
                            )

                logger.info(f"Fitted decay: a={best_a:.4f}, b={best_b:.6f}, c={best_c:.4f}")
            except Exception as fit_err:
                logger.warning(f"Could not fit decay curve: {fit_err}")

        # Draw r2 = 0.2 threshold line
        ax.axhline(y=0.2, color="gray", linestyle=":", alpha=0.6, label="r\u00b2 = 0.2")

        ax.set_xlabel("Distance (kb)", fontsize=12)
        ax.set_ylabel("Mean r\u00b2", fontsize=12)
        ax.set_title(title, fontsize=14, pad=15)
        ax.set_ylim(bottom=0)
        ax.legend(loc="upper right")
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        output_path_str: Optional[str] = None
        if output_file:
            output_file = Path(output_file)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_file, dpi=300, bbox_inches="tight")
            output_path_str = str(output_file)
            logger.info(f"Saved LD decay plot to {output_file}")

        plt.close(fig)

        result: Dict[str, Any] = {
            "status": "success",
            "output_path": output_path_str,
        }
        if ld_decay_distance is not None:
            result["ld_decay_distance_kb"] = ld_decay_distance

        return result

    except Exception as e:
        logger.error(f"Error creating LD decay plot: {e}")
        return {"status": "failed", "error": str(e), "output_path": None}


def ld_heatmap_region(
    genotypes_by_variant: List[List[int]],
    positions: List[int],
    output_file: Optional[Union[str, Path]] = None,
    region_start: Optional[int] = None,
    region_end: Optional[int] = None,
    title: str = "LD Heatmap",
) -> Dict[str, Any]:
    """Create a triangular LD heatmap for a genomic region (Haploview style).

    Computes pairwise r-squared for all variants in the specified region and
    displays them as a rotated (45-degree) triangular heatmap with a white-to-red
    color scale.

    Args:
        genotypes_by_variant: Genotype matrix (variants x samples), values 0/1/2/-1.
        positions: Physical positions (bp) for each variant.
        output_file: Optional path to save the plot image.
        region_start: Start position of the region (bp). Defaults to min position.
        region_end: End position of the region (bp). Defaults to max position.
        title: Plot title.

    Returns:
        Dictionary with status and output_path.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for LD heatmap")
        return {"status": "skipped", "output_path": None}

    if not HAS_NUMPY:
        logger.warning("numpy not available for LD heatmap")
        return {"status": "skipped", "output_path": None}

    try:
        n_variants = len(genotypes_by_variant)
        if n_variants < 2:
            return {
                "status": "failed",
                "error": f"Need at least 2 variants, got {n_variants}",
                "output_path": None,
            }

        if len(positions) != n_variants:
            return {
                "status": "failed",
                "error": "Length of positions must match number of variants",
                "output_path": None,
            }

        # Filter to region if specified
        if region_start is None:
            region_start = min(positions)
        if region_end is None:
            region_end = max(positions)

        # Select variants within region
        region_indices = [i for i in range(n_variants) if region_start <= positions[i] <= region_end]

        if len(region_indices) < 2:
            return {
                "status": "failed",
                "error": f"Fewer than 2 variants in region [{region_start}, {region_end}]",
                "output_path": None,
            }

        n_region = len(region_indices)
        region_positions = [positions[i] for i in region_indices]

        logger.info(f"Computing LD heatmap for {n_region} variants in region " f"[{region_start}, {region_end}]")

        # Compute pairwise r-squared matrix
        r2_matrix = np.zeros((n_region, n_region), dtype=np.float64)
        for ii in range(n_region):
            r2_matrix[ii, ii] = 1.0
            for jj in range(ii + 1, n_region):
                idx_i = region_indices[ii]
                idx_j = region_indices[jj]
                r2 = _compute_r_squared(genotypes_by_variant[idx_i], genotypes_by_variant[idx_j])
                r2_matrix[ii, jj] = r2
                r2_matrix[jj, ii] = r2

        # Create the triangular heatmap (rotated 45 degrees)
        fig, ax = plt.subplots(figsize=(10, 6))

        # Build rotated coordinates for a triangular display
        # Each cell (i, j) where j > i is drawn as a diamond rotated 45 degrees
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Polygon

        patches = []
        colors_list = []

        # Normalize positions for plotting
        pos_arr = np.array(region_positions, dtype=np.float64)
        pos_min = pos_arr.min()
        pos_range = pos_arr.max() - pos_arr.min()
        if pos_range == 0:
            pos_range = 1.0

        # Scale positions to [0, n_region - 1] for even spacing
        scaled_pos = np.arange(n_region, dtype=np.float64)

        for ii in range(n_region):
            for jj in range(ii + 1, n_region):
                # Diamond vertices in rotated coordinate system
                x_center = (scaled_pos[ii] + scaled_pos[jj]) / 2.0
                y_center = (scaled_pos[jj] - scaled_pos[ii]) / 2.0
                half_w = 0.5

                diamond = [
                    (x_center, y_center - half_w),
                    (x_center + half_w, y_center),
                    (x_center, y_center + half_w),
                    (x_center - half_w, y_center),
                ]
                patches.append(Polygon(diamond, closed=True))
                colors_list.append(r2_matrix[ii, jj])

        if patches:
            collection = PatchCollection(patches, cmap="Reds", edgecolors="none")
            collection.set_array(np.array(colors_list))
            collection.set_clim(0.0, 1.0)
            ax.add_collection(collection)

            # Add colorbar
            cbar = plt.colorbar(collection, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label("r\u00b2", fontsize=11)

        # Set axis limits
        ax.set_xlim(-0.5, n_region - 0.5)
        ax.set_ylim(-0.5, n_region / 2.0 + 0.5)

        # Add position labels on x-axis
        if n_region <= 30:
            label_positions = list(range(n_region))
            label_texts = [f"{region_positions[i] / 1000:.1f}" for i in range(n_region)]
            ax.set_xticks(label_positions)
            ax.set_xticklabels(label_texts, rotation=45, ha="right", fontsize=8)
            ax.set_xlabel("Position (kb)", fontsize=11)
        else:
            # Show a subset of labels
            step = max(1, n_region // 10)
            label_positions = list(range(0, n_region, step))
            label_texts = [f"{region_positions[i] / 1000:.1f}" for i in label_positions]
            ax.set_xticks(label_positions)
            ax.set_xticklabels(label_texts, rotation=45, ha="right", fontsize=8)
            ax.set_xlabel("Position (kb)", fontsize=11)

        ax.set_yticks([])
        ax.set_title(title, fontsize=14, pad=15)

        plt.tight_layout()

        output_path_str: Optional[str] = None
        if output_file:
            output_file = Path(output_file)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_file, dpi=300, bbox_inches="tight")
            output_path_str = str(output_file)
            logger.info(f"Saved LD heatmap to {output_file}")

        plt.close(fig)
        return {
            "status": "success",
            "output_path": output_path_str,
            "n_variants": n_region,
        }

    except Exception as e:
        logger.error(f"Error creating LD heatmap: {e}")
        return {"status": "failed", "error": str(e), "output_path": None}
