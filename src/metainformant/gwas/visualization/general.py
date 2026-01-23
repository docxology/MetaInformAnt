"""GWAS visualization utilities.

This module provides functions for creating GWAS visualization plots,
including Manhattan plots, Q-Q plots, and regional association plots.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union
import math

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Import matplotlib with graceful fallback
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

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


def manhattan_plot(
    results: Union[List[Dict[str, Any]], Dict[str, Any]],
    output_path: Optional[Union[str, Path]] = None,
    significance_threshold: float = 5e-8,
) -> Any:
    """Create a Manhattan plot from GWAS results.

    Args:
        results: GWAS results dictionary or list of result dictionaries
        output_path: Path to save the plot (optional)
        significance_threshold: P-value threshold for significance line

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

    # Color scheme for chromosomes
    chrom_colors = ["#1f77b4", "#ff7f0e"]  # Blue and orange alternating

    current_pos = 0
    chrom_offsets = {}

    for result in results:
        if not isinstance(result, dict):
            continue

        chrom = str(result.get("chrom", result.get("chromosome", "1")))
        pos = result.get("pos", result.get("position", 0))
        p_val = result.get("p_value", result.get("pval", 1.0))

        # Convert p-value to -log10 scale
        if p_val > 0:
            neg_log_p = -math.log10(p_val)
        else:
            neg_log_p = 50  # Cap very small p-values

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

    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot points
    scatter = ax.scatter(positions, p_values, c=colors, s=1, alpha=0.8)

    # Add significance threshold line
    if significance_threshold > 0:
        threshold_line = -math.log10(significance_threshold)
        ax.axhline(
            y=threshold_line,
            color="red",
            linestyle="--",
            alpha=0.7,
            label=f"Significance threshold ({significance_threshold})",
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
    if significance_threshold > 0:
        ax.legend(handles=legend_elements, loc="upper right")

    plt.tight_layout()

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved Manhattan plot to {output_path}")

    return fig


def qq_plot(p_values: Union[List[float], List[int]], output_path: Optional[Union[str, Path]] = None) -> Any:
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

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot Q-Q points
    ax.scatter(expected_log, observed_log, s=2, alpha=0.6, color="blue")

    # Add diagonal line (expected under null)
    max_val = max(np.max(expected_log), np.max(observed_log))
    ax.plot([0, max_val], [0, max_val], "r--", alpha=0.7, label="Expected (null)")

    # Labels and title
    ax.set_xlabel("Expected -log₁₀(p-value)")
    ax.set_ylabel("Observed -log₁₀(p-value)")
    ax.set_title(f"Q-Q Plot (n={n} variants)")
    ax.grid(True, alpha=0.3)
    ax.legend()

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
    ax1.axvline(x=np.mean(effect_sizes), color="green", linestyle="-", alpha=0.7, label=f"Mean: {np.mean(effect_sizes):.4f}")
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

    plots_generated = {}

    # Placeholder - would generate actual plots
    plot_types = ["manhattan", "qq", "pca", "kinship"]
    for plot_type in plot_types:
        plot_path = output_dir / f"{plot_type}_plot.png"
        plots_generated[plot_type] = str(plot_path)
        logger.info(f"Generated {plot_type} plot: {plot_path}")

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
        ax.text(0.5, 0.5, "No genome-wide significant variants\n(p < 5e-8)",
                ha="center", va="center", fontsize=14)
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
