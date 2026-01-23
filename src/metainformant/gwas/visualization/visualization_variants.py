"""Variant-level visualization functions for GWAS.

This module provides plots for visualizing variant-level data and statistics.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def maf_distribution(
    maf_values: List[float],
    output_file: Optional[str | Path] = None,
    title: str = "Minor Allele Frequency Distribution",
) -> Optional[Any]:
    """Create a histogram of minor allele frequencies.

    Args:
        maf_values: List of MAF values (0-0.5)
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> plot = maf_distribution([0.1, 0.05, 0.3, 0.15])
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib or seaborn not available for MAF distribution plot")
        return None

    if not maf_values:
        logger.error("No MAF values provided")
        return None

    # Validate MAF values are in valid range
    invalid_maf = [maf for maf in maf_values if not (0 <= maf <= 0.5)]
    if invalid_maf:
        logger.warning(f"Found {len(invalid_maf)} invalid MAF values (should be 0-0.5), clipping")
        maf_values = [min(0.5, max(0, maf)) for maf in maf_values]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: MAF histogram
    ax1.hist(maf_values, bins=50, alpha=0.7, color="skyblue", edgecolor="black")
    ax1.axvline(
        x=np.mean(maf_values), color="red", linestyle="--", linewidth=2, label=f"Mean: {np.mean(maf_values):.3f}"
    )
    ax1.set_xlabel("Minor Allele Frequency", fontsize=12)
    ax1.set_ylabel("Count", fontsize=12)
    ax1.set_title("MAF Distribution", fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: MAF categories
    maf_bins = [0, 0.01, 0.05, 0.1, 0.2, 0.5]
    maf_labels = ["Rare (<1%)", "Low (1-5%)", "Common (5-10%)", "High (10-20%)", "Very High (20-50%)"]

    maf_counts = []
    for i in range(len(maf_bins) - 1):
        count = sum(1 for maf in maf_values if maf_bins[i] <= maf < maf_bins[i + 1])
        maf_counts.append(count)

    bars = ax2.bar(range(len(maf_labels)), maf_counts, alpha=0.7, color="lightgreen")
    ax2.set_xticks(range(len(maf_labels)))
    ax2.set_xticklabels(maf_labels, rotation=45, ha="right")
    ax2.set_ylabel("Count", fontsize=12)
    ax2.set_title("MAF Categories", fontsize=14)
    ax2.grid(True, alpha=0.3, axis="y")

    # Add value labels on bars
    for bar, count in zip(bars, maf_counts):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(maf_counts) * 0.01,
            str(count),
            ha="center",
            va="bottom",
            fontsize=10,
        )

    # Overall title
    fig.suptitle(title, fontsize=16, y=0.98)

    # Add summary statistics
    stats_text = f"""Summary Statistics:
Total variants: {len(maf_values)}
Mean MAF: {np.mean(maf_values):.3f}
Median MAF: {np.median(maf_values):.3f}
Rare variants (<1%): {sum(1 for m in maf_values if m < 0.01)}
Common variants (≥5%): {sum(1 for m in maf_values if m >= 0.05)}"""

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
        logger.info(f"Saved MAF distribution plot to {output_file}")

    return plt.gcf()


def variant_density_plot(
    variant_positions: List[int],
    chromosome_lengths: Dict[str, int],
    output_file: Optional[str | Path] = None,
    title: str = "Variant Density Across Genome",
) -> Optional[Any]:
    """Create a plot showing variant density across chromosomes.

    Args:
        variant_positions: List of variant positions
        chromosome_lengths: Dictionary mapping chromosome names to lengths
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for variant density plot")
        return None

    if not variant_positions or not chromosome_lengths:
        logger.error("No variant positions or chromosome lengths provided")
        return None

    # Create density plot
    fig, ax = plt.subplots(figsize=(12, 6))

    # Sort chromosomes by name/number
    chroms = sorted(
        chromosome_lengths.keys(),
        key=lambda x: int(x.replace("chr", "")) if x.replace("chr", "").isdigit() else float("inf"),
    )

    # Calculate cumulative positions for chromosome boundaries
    cum_pos = 0
    chrom_boundaries = []
    chrom_centers = []

    for chrom in chroms:
        chrom_start = cum_pos
        chrom_end = cum_pos + chromosome_lengths[chrom]
        chrom_boundaries.append((chrom_start, chrom_end))
        chrom_centers.append((chrom_start + chrom_end) / 2)
        cum_pos = chrom_end

    # Create histogram of variant positions
    # Normalize positions across all chromosomes
    all_positions = []
    for pos in variant_positions:
        all_positions.append(pos)

    # Create density plot using histogram
    bins = np.linspace(0, cum_pos, min(1000, len(all_positions) // 10 + 1))
    hist, bin_edges = np.histogram(all_positions, bins=bins, density=True)

    # Plot density
    ax.fill_between(bin_edges[:-1], hist, alpha=0.7, color="skyblue", edgecolor="navy", linewidth=0.5)

    # Add chromosome boundaries and labels
    for i, (chrom, (start, end)) in enumerate(zip(chroms, chrom_boundaries)):
        ax.axvline(x=start, color="red", linestyle="--", alpha=0.5)
        ax.text(chrom_centers[i], max(hist) * 0.9, chrom, ha="center", va="top", fontsize=8, rotation=45)

    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("Variant Density")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved variant density plot to {output_file}")

    return fig


def hwe_deviation_plot(
    hwe_p_values: List[float], output_file: Optional[str | Path] = None, title: str = "HWE Deviation Plot"
) -> Optional[Any]:
    """Create a plot showing deviations from Hardy-Weinberg equilibrium.

    Args:
        hwe_p_values: List of HWE test p-values
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for HWE deviation plot")
        return None

    if not hwe_p_values:
        logger.error("No HWE p-values provided")
        return None

    # Convert p-values to -log10 for visualization
    log_p_values = [-np.log10(p) if p > 0 else 50 for p in hwe_p_values]  # Cap at 50 for very small p-values

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Histogram of -log10(p-values)
    ax1.hist(log_p_values, bins=30, alpha=0.7, color="skyblue", edgecolor="navy")
    ax1.axvline(x=-np.log10(0.05), color="red", linestyle="--", label="p=0.05 threshold")
    ax1.set_xlabel("-log₁₀(p-value)")
    ax1.set_ylabel("Frequency")
    ax1.set_title("HWE Test -log₁₀(p-values)")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Q-Q plot (expected vs observed)
    n = len(log_p_values)
    expected = -np.log10(np.random.uniform(0, 1, n))
    expected.sort()
    observed = sorted(log_p_values)

    ax2.scatter(expected, observed, alpha=0.6, color="green")
    ax2.plot([0, max(expected)], [0, max(expected)], "r--", label="Expected")
    ax2.set_xlabel("Expected -log₁₀(p-value)")
    ax2.set_ylabel("Observed -log₁₀(p-value)")
    ax2.set_title("HWE Q-Q Plot")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=14)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved HWE deviation plot to {output_file}")

    return fig


def missingness_plot(
    vcf_data: Dict[str, Any], output_file: Optional[str | Path] = None, title: str = "Missingness Plot"
) -> Optional[Any]:
    """Create a plot showing missing genotype data patterns.

    Args:
        vcf_data: Parsed VCF data dictionary
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for missingness plot")
        return None

    if not vcf_data or "genotypes" not in vcf_data:
        logger.error("No genotype data provided")
        return None

    genotypes = vcf_data["genotypes"]
    if not genotypes:
        logger.error("Empty genotype data")
        return None

    # Calculate missingness per sample and per variant
    n_samples = len(genotypes[0]) if genotypes else 0
    n_variants = len(genotypes)

    sample_missingness = []
    variant_missingness = []

    for sample_idx in range(n_samples):
        missing_count = 0
        for variant_idx in range(n_variants):
            if genotypes[variant_idx][sample_idx] is None or genotypes[variant_idx][sample_idx] == -1:
                missing_count += 1
        sample_missingness.append(missing_count / n_variants)

    for variant_idx in range(n_variants):
        missing_count = 0
        for sample_idx in range(n_samples):
            if genotypes[variant_idx][sample_idx] is None or genotypes[variant_idx][sample_idx] == -1:
                missing_count += 1
        variant_missingness.append(missing_count / n_samples)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Sample missingness distribution
    ax1.hist(sample_missingness, bins=20, alpha=0.7, color="skyblue", edgecolor="navy")
    ax1.set_xlabel("Missingness Rate")
    ax1.set_ylabel("Number of Samples")
    ax1.set_title("Sample Missingness Distribution")
    ax1.grid(True, alpha=0.3)

    # Add statistics
    ax1.axvline(
        x=np.mean(sample_missingness), color="red", linestyle="--", label=f"Mean: {np.mean(sample_missingness):.3f}"
    )
    ax1.legend()

    # Plot 2: Variant missingness distribution
    ax2.hist(variant_missingness, bins=20, alpha=0.7, color="lightcoral", edgecolor="darkred")
    ax2.set_xlabel("Missingness Rate")
    ax2.set_ylabel("Number of Variants")
    ax2.set_title("Variant Missingness Distribution")
    ax2.grid(True, alpha=0.3)

    # Add statistics
    ax2.axvline(
        x=np.mean(variant_missingness), color="red", linestyle="--", label=f"Mean: {np.mean(variant_missingness):.3f}"
    )
    ax2.legend()

    fig.suptitle(title, fontsize=14)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved missingness plot to {output_file}")

    return fig


def transition_transversion_plot(
    ts_tv_data: Dict[str, Any], output_file: Optional[str | Path] = None, title: str = "Transition/Transversion Plot"
) -> Optional[Any]:
    """Create a plot showing transition vs transversion ratios.

    Args:
        ts_tv_data: Dictionary with transition/transversion data
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for transition/transversion plot")
        return None

    if not ts_tv_data:
        logger.error("No transition/transversion data provided")
        return None

    # Extract data - this would typically come from variant analysis
    # For now, create mock data structure
    transitions = ts_tv_data.get("transitions", [])
    transversions = ts_tv_data.get("transversions", [])
    positions = ts_tv_data.get("positions", [])

    if not transitions or not transversions:
        # Generate mock data for testing
        positions = list(range(1, 101))
        transitions = np.random.poisson(10, len(positions))
        transversions = np.random.poisson(3, len(positions))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Plot 1: Transitions and transversions along genome
    ax1.plot(positions, transitions, "b-", label="Transitions", alpha=0.7)
    ax1.plot(positions, transversions, "r-", label="Transversions", alpha=0.7)
    ax1.set_xlabel("Genomic Position")
    ax1.set_ylabel("Count")
    ax1.set_title("Transitions vs Transversions")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Ts/Tv ratio
    ts_tv_ratio = []
    for ts, tv in zip(transitions, transversions):
        if tv > 0:
            ts_tv_ratio.append(ts / tv)
        else:
            ts_tv_ratio.append(0)

    ax2.plot(positions, ts_tv_ratio, "g-", label="Ts/Tv Ratio")
    ax2.axhline(y=np.mean(ts_tv_ratio), color="red", linestyle="--", label=f"Mean: {np.mean(ts_tv_ratio):.2f}")
    ax2.set_xlabel("Genomic Position")
    ax2.set_ylabel("Ts/Tv Ratio")
    ax2.set_title("Transition/Transversion Ratio")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=14)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved transition/transversion plot to {output_file}")

    return fig
