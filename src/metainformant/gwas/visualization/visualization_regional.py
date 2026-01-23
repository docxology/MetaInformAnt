"""Regional visualization functions for GWAS.

This module provides plots for regional association analysis and gene annotation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def gene_annotation_plot(
    results_df: Any,
    chrom: str,
    start: int,
    end: int,
    gene_data: Optional[Dict[str, Any]] = None,
    output_file: Optional[str | Path] = None,
    title: str = "Regional Association Plot",
) -> Optional[Any]:
    """Create a regional association plot with gene annotations.

    Args:
        results_df: DataFrame with GWAS results (columns: CHR, BP, P)
        chrom: Chromosome to plot
        start: Start position (bp)
        end: End position (bp)
        gene_data: Optional dictionary with gene annotations
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> plot = gene_annotation_plot(df, 'chr1', 100000, 200000)
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
    except ImportError:
        logger.warning("matplotlib not available for regional plot")
        return None

    # Filter data for region
    if not hasattr(results_df, "columns"):
        logger.error("Input data must be a DataFrame")
        return None

    region_data = results_df[
        (results_df["CHR"] == chrom) & (results_df["BP"] >= start) & (results_df["BP"] <= end)
    ].copy()

    if region_data.empty:
        logger.warning(f"No data found for region {chrom}:{start}-{end}")
        return None

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw={"height_ratios": [3, 1]})

    # Plot 1: Manhattan plot for region
    region_data["neg_log_p"] = -np.log10(region_data["P"].clip(lower=1e-50))

    ax1.scatter(region_data["BP"], region_data["neg_log_p"], c="blue", s=20, alpha=0.7)
    ax1.axhline(y=-np.log10(5e-8), color="red", linestyle="--", alpha=0.8, label="Genome-wide significance")
    ax1.set_ylabel("-log₁₀(P-value)", fontsize=12)
    ax1.set_title(title, fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Gene annotation track
    if gene_data:
        plot_gene_track(ax2, chrom, start, end, gene_data)
    else:
        # Default gene track (simplified)
        ax2.text(0.5, 0.5, "Gene annotations not provided", ha="center", va="center", transform=ax2.transAxes)
        ax2.set_xlim(start, end)

    ax2.set_xlabel(f"Position on {chrom} (bp)", fontsize=12)
    ax2.set_yticks([])

    # Format x-axis labels
    def format_bp(x, pos):
        if x >= 1e6:
            return ".1f"
        elif x >= 1e3:
            return ".0f"
        else:
            return str(int(x))

    ax1.xaxis.set_major_formatter(plt.FuncFormatter(format_bp))
    ax2.xaxis.set_major_formatter(plt.FuncFormatter(format_bp))

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved regional plot to {output_file}")

    return plt.gcf()


def plot_gene_track(ax, chrom: str, start: int, end: int, gene_data: Dict[str, Any]):
    """Plot gene annotation track."""
    # Simplified gene plotting - in practice would use genomic annotation data
    if "genes" in gene_data:
        genes = gene_data["genes"]
        y_pos = 0.5

        for gene in genes:
            gene_start = max(start, gene.get("start", start))
            gene_end = min(end, gene.get("end", end))

            if gene_start < gene_end:
                # Draw gene rectangle
                width = gene_end - gene_start
                rect = Rectangle((gene_start, y_pos - 0.1), width, 0.2, facecolor="lightblue", alpha=0.7)
                ax.add_patch(rect)

                # Add gene name if space allows
                if width > (end - start) * 0.05:  # If gene covers >5% of region
                    ax.text(
                        (gene_start + gene_end) / 2,
                        y_pos,
                        gene.get("name", ""),
                        ha="center",
                        va="center",
                        fontsize=8,
                        rotation=45,
                    )

    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Genes", fontsize=10)


def recombination_rate_plot(
    recombination_rates: Dict[str, List[float]],
    positions: Dict[str, List[int]],
    output_file: Optional[str | Path] = None,
    title: str = "Recombination Rate Profile",
) -> Optional[Any]:
    """Create a plot showing recombination rates across genomic regions.

    Args:
        recombination_rates: Dictionary mapping chromosomes to recombination rate lists
        positions: Dictionary mapping chromosomes to position lists
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for recombination rate plot")
        return None

    if not recombination_rates or not positions:
        logger.error("No recombination rate or position data provided")
        return None

    # Create figure with subplots for each chromosome
    chromosomes = list(recombination_rates.keys())
    n_chroms = len(chromosomes)

    if n_chroms <= 3:
        n_cols = n_chroms
        n_rows = 1
    else:
        n_cols = 3
        n_rows = (n_chroms + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows), squeeze=False, sharey=True)
    fig.suptitle(title, fontsize=16, y=0.95)

    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

    for i, chrom in enumerate(chromosomes):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]

        rates = recombination_rates[chrom]
        pos = positions[chrom]

        if len(rates) != len(pos):
            logger.warning(f"Mismatched data lengths for {chrom}, skipping")
            continue

        # Plot recombination rate
        ax.plot(pos, rates, color=colors[i % len(colors)], linewidth=1.5, alpha=0.8)

        # Add smoothing line
        if len(rates) > 10:
            from scipy import signal

            window_size = min(21, len(rates) // 10 * 2 + 1)  # Adaptive window size
            if window_size % 2 == 0:
                window_size += 1
            smoothed = signal.savgol_filter(rates, window_size, 3)
            ax.plot(pos, smoothed, color="red", linewidth=2, alpha=0.6, label="Smoothed")

        ax.set_title(f"Chromosome {chrom}", fontsize=12)
        ax.set_xlabel("Position (bp)", fontsize=10)
        if col == 0:  # Leftmost column
            ax.set_ylabel("Recombination Rate (cM/Mb)", fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend()

    # Hide unused subplots
    for i in range(n_chroms, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].set_visible(False)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved recombination rate plot to {output_file}")

    return fig


def effect_direction_plot(
    effect_sizes: List[float],
    positions: List[int],
    chrom: str,
    output_file: Optional[str | Path] = None,
    title: str = "Effect Direction Plot",
) -> Optional[Any]:
    """Create a plot showing effect directions (positive/negative) across a region.

    Args:
        effect_sizes: List of effect sizes (can be positive or negative)
        positions: List of genomic positions
        chrom: Chromosome name
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for effect direction plot")
        return None

    if not effect_sizes or not positions:
        logger.error("No effect size or position data provided")
        return None

    if len(effect_sizes) != len(positions):
        logger.error("Effect sizes and positions must have same length")
        return None

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    fig.suptitle(f"{title} - {chrom}", fontsize=14)

    # Plot 1: Effect sizes with color coding for direction
    colors = ["red" if x < 0 else "blue" for x in effect_sizes]
    ax1.scatter(positions, effect_sizes, c=colors, alpha=0.7, s=20)
    ax1.axhline(y=0, color="black", linestyle="--", alpha=0.5)
    ax1.set_ylabel("Effect Size", fontsize=12)
    ax1.set_title("Effect Sizes by Direction", fontsize=12)
    ax1.grid(True, alpha=0.3)

    # Add legend
    ax1.scatter([], [], c="red", alpha=0.7, s=20, label="Negative effect")
    ax1.scatter([], [], c="blue", alpha=0.7, s=20, label="Positive effect")
    ax1.legend()

    # Plot 2: Effect size distribution
    ax2.hist(effect_sizes, bins=30, alpha=0.7, color="skyblue", edgecolor="navy")
    ax2.axvline(x=0, color="red", linestyle="--", linewidth=2, label="No effect")
    ax2.axvline(
        x=np.mean(effect_sizes), color="green", linestyle="-", linewidth=2, label=f"Mean: {np.mean(effect_sizes):.3f}"
    )
    ax2.set_xlabel("Effect Size", fontsize=12)
    ax2.set_ylabel("Frequency", fontsize=12)
    ax2.set_title("Effect Size Distribution", fontsize=12)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Add statistics text
    stats_text = f"""Statistics:
Mean: {np.mean(effect_sizes):.3f}
Median: {np.median(effect_sizes):.3f}
SD: {np.std(effect_sizes):.3f}
Positive effects: {sum(1 for x in effect_sizes if x > 0)} ({sum(1 for x in effect_sizes if x > 0)/len(effect_sizes)*100:.1f}%)
Negative effects: {sum(1 for x in effect_sizes if x < 0)} ({sum(1 for x in effect_sizes if x < 0)/len(effect_sizes)*100:.1f}%)"""

    fig.text(
        0.02,
        0.02,
        stats_text,
        fontsize=10,
        verticalalignment="bottom",
        fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved effect direction plot to {output_file}")

    return fig


def regional_plot(
    results: pd.DataFrame,
    chr: str,
    start: int,
    end: int,
    output_file: Optional[str | Path] = None,
    figsize: tuple[int, int] = (12, 6),
    title: Optional[str] = None,
    significance_threshold: float = 5e-8,
    recombination_rate: Optional[pd.DataFrame] = None,
) -> Optional[Any]:
    """Create a regional association plot.

    Args:
        results: GWAS results DataFrame with CHR, BP, P columns
        chr: Chromosome to plot
        start: Start position (bp)
        end: End position (bp)
        output_file: Optional output file path
        figsize: Figure size (width, height)
        title: Plot title (auto-generated if None)
        significance_threshold: P-value threshold for significance line
        recombination_rate: Optional recombination rate data

    Returns:
        Matplotlib figure if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for regional plot")
        return None

    if not HAS_SEABORN:
        logger.warning("seaborn not available for regional plot")
        return None

    # Filter data for the specified region
    region_data = results[(results["CHR"] == chr) & (results["BP"] >= start) & (results["BP"] <= end)].copy()

    if region_data.empty:
        logger.warning(f"No data found for chr {chr}:{start}-{end}")
        return None

    # Convert p-values to -log10 scale
    region_data["NEG_LOG_P"] = -np.log10(region_data["P"].clip(lower=1e-50))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, gridspec_kw={"height_ratios": [3, 1]})

    # Main plot: -log10(p) vs position
    ax1.scatter(region_data["BP"], region_data["NEG_LOG_P"], c=region_data["NEG_LOG_P"], cmap="Reds", alpha=0.7, s=20)

    # Add significance threshold line
    ax1.axhline(
        y=-np.log10(significance_threshold),
        color="red",
        linestyle="--",
        alpha=0.7,
        label=f"p = {significance_threshold}",
    )

    # Add recombination rate if provided
    if recombination_rate is not None:
        # Plot recombination rate on second axis
        ax2_twin = ax2.twinx()
        recomb_data = recombination_rate[
            (recombination_rate["CHR"] == chr) & (recombination_rate["BP"] >= start) & (recombination_rate["BP"] <= end)
        ]
        if not recomb_data.empty:
            ax2_twin.plot(recomb_data["BP"], recomb_data["RATE"], color="blue", alpha=0.7, linewidth=2)
            ax2_twin.set_ylabel("Recombination Rate (cM/Mb)", color="blue")
            ax2_twin.tick_params(axis="y", labelcolor="blue")

    # Format axes
    ax1.set_xlim(start, end)
    ax1.set_ylabel("-log₁₀(p)")
    ax1.set_title(title or f"GWAS Regional Plot - Chr {chr}:{start:,}-{end:,}")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Add gene annotations if available (mock for now)
    # In a real implementation, this would use gene annotation data
    genes_in_region = []  # Would query gene database
    if genes_in_region:
        # Add gene labels
        for gene in genes_in_region:
            if start <= gene["start"] <= end:
                ax1.axvline(x=gene["start"], color="green", linestyle=":", alpha=0.5)
                ax1.text(
                    gene["start"], ax1.get_ylim()[1] * 0.9, gene["name"], rotation=90, ha="center", va="top", fontsize=8
                )

    # Bottom plot: LD heatmap (placeholder - would use LD data)
    # For now, just show position density
    ax2.hist(region_data["BP"], bins=50, alpha=0.7, color="lightblue", edgecolor="navy")
    ax2.set_xlabel("Position (bp)")
    ax2.set_ylabel("SNP Density")
    ax2.set_xlim(start, end)
    ax2.grid(True, alpha=0.3)

    # Format x-axis labels
    def format_bp(x, pos):
        if x >= 1e6:
            return f"{x/1e6:.1f}M"
        elif x >= 1e3:
            return f"{x/1e3:.0f}K"
        else:
            return f"{x:.0f}"

    ax2.xaxis.set_major_formatter(plt.FuncFormatter(format_bp))

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved regional plot to {output_file}")

    return fig


def regional_ld_plot(
    vcf_path: str | Path,
    output_file: Optional[str | Path] = None,
    chrom: str = "chr1",
    start: int = 1000,
    end: int = 2000,
    lead_snp_pos: Optional[int] = None,
    ld_threshold: float = 0.8,
    figsize: tuple[int, int] = (12, 8),
) -> Dict[str, Any]:
    """Create a regional LD plot showing linkage disequilibrium around a locus.

    Args:
        vcf_path: Path to VCF file
        output_file: Optional output file path
        chrom: Chromosome to plot
        start: Start position (bp)
        end: End position (bp)
        lead_snp_pos: Position of lead SNP to highlight
        ld_threshold: LD threshold for visualization
        figsize: Figure size (width, height)

    Returns:
        Dictionary with status and metadata
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for regional LD plot")
        return {"status": "skipped", "reason": "matplotlib not available"}

    # Check if VCF file exists
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        logger.warning(f"VCF file not found: {vcf_path}")
        return {"status": "skipped", "reason": "VCF file not found"}

    # Note: Real LD calculation would require tools like PLINK or vcftools
    # For now, generate mock LD data for testing
    logger.info("Generating mock LD data for regional plot (real LD calculation requires external tools)")

    # Generate mock SNP positions in the region
    n_snps = 20
    positions = np.linspace(start, end, n_snps, dtype=int)

    # Generate mock LD matrix (r² values)
    np.random.seed(42)  # For reproducible mock data
    ld_matrix = np.random.uniform(0, 1, (n_snps, n_snps))

    # Make LD matrix symmetric and add some structure
    ld_matrix = (ld_matrix + ld_matrix.T) / 2

    # Add distance-based decay (closer SNPs have higher LD)
    for i in range(n_snps):
        for j in range(n_snps):
            distance = abs(positions[i] - positions[j])
            decay_factor = np.exp(-distance / 10000)  # Exponential decay
            ld_matrix[i, j] *= decay_factor

    # Ensure diagonal is 1.0
    np.fill_diagonal(ld_matrix, 1.0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, gridspec_kw={"height_ratios": [1, 3]})

    # Top plot: SNP positions
    ax1.scatter(positions, [1] * n_snps, c="blue", s=50, alpha=0.7)
    if lead_snp_pos:
        # Highlight lead SNP
        lead_idx = np.argmin(np.abs(positions - lead_snp_pos))
        ax1.scatter(positions[lead_idx], 1, c="red", s=100, marker="*", label=f"Lead SNP (pos {lead_snp_pos})")

    ax1.set_xlim(start, end)
    ax1.set_ylim(0.5, 1.5)
    ax1.set_xlabel("Position (bp)")
    ax1.set_title(f"Regional LD Plot - {chrom}:{start:,}-{end:,}")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Format x-axis
    def format_bp(x, pos):
        if x >= 1e6:
            return f"{x/1e6:.1f}M"
        elif x >= 1e3:
            return f"{x/1e3:.0f}K"
        else:
            return f"{x:.0f}"

    ax1.xaxis.set_major_formatter(plt.FuncFormatter(format_bp))

    # Bottom plot: LD heatmap
    im = ax2.imshow(ld_matrix, cmap="RdYlBu_r", aspect="auto", extent=[start, end, end, start], vmin=0, vmax=1)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax2, shrink=0.8)
    cbar.set_label("LD (r²)")

    # Highlight regions above LD threshold
    if ld_threshold > 0:
        high_ld_mask = ld_matrix >= ld_threshold
        if np.any(high_ld_mask):
            # Add contour lines for high LD regions
            ax2.contour(ld_matrix, levels=[ld_threshold], colors="red", linewidths=2, extent=[start, end, end, start])

    ax2.set_xlabel("Position 1 (bp)")
    ax2.set_ylabel("Position 2 (bp)")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved regional LD plot to {output_file}")

    return {
        "status": "success",
        "n_snps": n_snps,
        "chromosome": chrom,
        "region_start": start,
        "region_end": end,
        "lead_snp_pos": lead_snp_pos,
        "ld_threshold": ld_threshold,
        "mock_data": True,  # Indicates this used mock LD data
    }
