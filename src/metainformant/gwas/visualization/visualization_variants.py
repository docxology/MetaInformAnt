"""Variant-level visualization functions for GWAS.

This module provides plots for visualizing variant-level data and statistics.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def _extract_float_values(data: Any, key: str) -> List[float]:
    """Extract float values from input data.

    Accepts:
        - List of dicts with the given key (e.g. [{"MAF": 0.1}, ...])
        - List/array of float values
        - numpy array

    Returns:
        List of float values.
    """
    if isinstance(data, np.ndarray):
        return data.tolist()
    if not data:
        return []
    if isinstance(data, (list, tuple)):
        if len(data) > 0 and isinstance(data[0], dict):
            return [float(item[key]) for item in data if key in item]
        return [float(v) for v in data]
    return []


def maf_distribution(
    maf_values: Any,
    output_file: Optional[str | Path] = None,
    title: str = "Minor Allele Frequency Distribution",
) -> Dict[str, Any]:
    """Create a histogram of minor allele frequencies.

    Args:
        maf_values: List of MAF values (0-0.5), or list of dicts with "MAF" key
        output_file: Optional output file path
        title: Plot title

    Returns:
        Dict with status and metadata

    Example:
        >>> result = maf_distribution([0.1, 0.05, 0.3, 0.15])
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for MAF distribution plot")
        return {"status": "failed", "error": "matplotlib not available"}

    # Extract float values from input (supports list-of-dicts or plain list)
    values = _extract_float_values(maf_values, "MAF")

    if not values or len(values) == 0:
        logger.error("No MAF values provided")
        return {"status": "failed", "error": "No MAF values provided"}

    # Validate MAF values are in valid range
    invalid_maf = [maf for maf in values if not (0 <= maf <= 0.5)]
    if invalid_maf:
        logger.warning(f"Found {len(invalid_maf)} invalid MAF values (should be 0-0.5), clipping")
        values = [min(0.5, max(0, maf)) for maf in values]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: MAF histogram
    ax1.hist(values, bins=50, alpha=0.7, color="skyblue", edgecolor="black")
    ax1.axvline(x=np.mean(values), color="red", linestyle="--", linewidth=2, label=f"Mean: {np.mean(values):.3f}")
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
        count = sum(1 for maf in values if maf_bins[i] <= maf < maf_bins[i + 1])
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
Total variants: {len(values)}
Mean MAF: {np.mean(values):.3f}
Median MAF: {np.median(values):.3f}
Rare variants (<1%): {sum(1 for m in values if m < 0.01)}
Common variants (>=5%): {sum(1 for m in values if m >= 0.05)}"""

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

    plt.close(fig)
    return {"status": "success", "n_variants": len(values), "mean_maf": float(np.mean(values))}


def variant_density_plot(
    variants: Any,
    output_file: Optional[str | Path] = None,
    title: str = "Variant Density Across Genome",
    window_size: int = 1000000,
) -> Dict[str, Any]:
    """Create a plot showing variant density across chromosomes.

    Args:
        variants: List of variant dicts with "CHROM" and "POS" keys,
                  or list of positions, or (positions, chromosome_lengths) tuple
        output_file: Optional output file path
        title: Plot title
        window_size: Window size for density binning (default 1Mb)

    Returns:
        Dict with status and metadata
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for variant density plot")
        return {"status": "failed", "error": "matplotlib not available"}

    # Parse input: accept list-of-dicts with CHROM/POS keys
    if isinstance(variants, (list, tuple)) and len(variants) > 0 and isinstance(variants[0], dict):
        # Extract chromosome and position from list of dicts
        chrom_positions: Dict[str, List[int]] = {}
        for v in variants:
            chrom = str(v.get("CHROM", "1"))
            pos = int(v.get("POS", 0))
            if chrom not in chrom_positions:
                chrom_positions[chrom] = []
            chrom_positions[chrom].append(pos)
    else:
        if not variants:
            logger.error("No variant positions provided")
            return {"status": "failed", "error": "No variant positions provided"}
        # Treat as plain list of positions on chromosome "1"
        chrom_positions = {"1": [int(p) for p in variants]}

    if not chrom_positions:
        logger.error("No variant positions provided")
        return {"status": "failed", "error": "No variant positions provided"}

    # Compute chromosome lengths from data if not explicitly provided
    chromosome_lengths: Dict[str, int] = {}
    for chrom, positions in chrom_positions.items():
        if positions:
            chromosome_lengths[chrom] = max(positions) + window_size

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

    # Build all positions with cumulative offset
    all_positions = []
    chrom_offset = {chrom: boundary[0] for chrom, boundary in zip(chroms, chrom_boundaries)}
    for chrom in chroms:
        for pos in chrom_positions.get(chrom, []):
            all_positions.append(pos + chrom_offset[chrom])

    if not all_positions:
        logger.error("No variant positions to plot")
        plt.close(fig)
        return {"status": "failed", "error": "No variant positions to plot"}

    # Create density plot using histogram
    n_bins = max(2, min(1000, len(all_positions) // 2 + 1))
    bins = np.linspace(0, cum_pos, n_bins)
    hist, bin_edges = np.histogram(all_positions, bins=bins, density=True)

    # Plot density
    ax.fill_between(bin_edges[:-1], hist, alpha=0.7, color="skyblue", edgecolor="navy", linewidth=0.5)

    # Add chromosome boundaries and labels
    max_hist = max(hist) if len(hist) > 0 else 1
    for i, (chrom, (start, end)) in enumerate(zip(chroms, chrom_boundaries)):
        ax.axvline(x=start, color="red", linestyle="--", alpha=0.5)
        ax.text(chrom_centers[i], max_hist * 0.9, chrom, ha="center", va="top", fontsize=8, rotation=45)

    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("Variant Density")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved variant density plot to {output_file}")

    plt.close(fig)
    return {"status": "success", "n_variants": len(all_positions), "n_chromosomes": len(chroms)}


def hwe_deviation_plot(
    hwe_p_values: Any, output_file: Optional[str | Path] = None, title: str = "HWE Deviation Plot"
) -> Dict[str, Any]:
    """Create a plot showing deviations from Hardy-Weinberg equilibrium.

    Args:
        hwe_p_values: List of HWE test p-values, or list of dicts with "HWE_P" key
        output_file: Optional output file path
        title: Plot title

    Returns:
        Dict with status and metadata
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for HWE deviation plot")
        return {"status": "failed", "error": "matplotlib not available"}

    # Extract float values from input (supports list-of-dicts or plain list)
    values = _extract_float_values(hwe_p_values, "HWE_P")

    if not values or len(values) == 0:
        logger.error("No HWE p-values provided")
        return {"status": "failed", "error": "No HWE p-values provided"}

    # Convert p-values to -log10 for visualization
    log_p_values = [-np.log10(p) if p > 0 else 50 for p in values]  # Cap at 50 for very small p-values

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Histogram of -log10(p-values)
    ax1.hist(log_p_values, bins=30, alpha=0.7, color="skyblue", edgecolor="navy")
    ax1.axvline(x=-np.log10(0.05), color="red", linestyle="--", label="p=0.05 threshold")
    ax1.set_xlabel("-log10(p-value)")
    ax1.set_ylabel("Frequency")
    ax1.set_title("HWE Test -log10(p-values)")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Q-Q plot (expected vs observed)
    n = len(log_p_values)
    expected = -np.log10(np.random.uniform(0, 1, n))
    expected.sort()
    observed = sorted(log_p_values)

    ax2.scatter(expected, observed, alpha=0.6, color="green")
    ax2.plot([0, max(expected)], [0, max(expected)], "r--", label="Expected")
    ax2.set_xlabel("Expected -log10(p-value)")
    ax2.set_ylabel("Observed -log10(p-value)")
    ax2.set_title("HWE Q-Q Plot")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=14)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved HWE deviation plot to {output_file}")

    plt.close(fig)
    return {"status": "success", "n_variants": len(values)}


def missingness_plot(
    vcf_path: Any,
    output_file: Optional[str | Path] = None,
    title: str = "Missingness Plot",
    by_sample: Optional[bool] = None,
) -> Dict[str, Any]:
    """Create a plot showing missing genotype data patterns.

    Args:
        vcf_path: Path to VCF file (string or Path) or parsed VCF data dictionary
        output_file: Optional output file path
        title: Plot title
        by_sample: If True, plot by sample; if False, by variant. Must be bool or None.

    Returns:
        Dict with status and metadata

    Raises:
        ValueError: If by_sample is not a bool, or if paths are empty/None.
    """
    # Validate by_sample parameter
    if by_sample is not None and not isinstance(by_sample, bool):
        raise ValueError("by_sample must be of type bool")

    # Validate vcf_path
    if isinstance(vcf_path, str) and vcf_path == "":
        raise ValueError("vcf_path cannot be empty")

    # Validate output_file
    if output_file is None:
        raise ValueError("output_path cannot be empty")

    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for missingness plot")
        return {"status": "failed", "error": "matplotlib not available"}

    # If vcf_path is a string/Path, try to parse it as a VCF file
    vcf_data: Optional[Dict[str, Any]] = None
    if isinstance(vcf_path, (str, Path)):
        vcf_file = Path(vcf_path)
        if not vcf_file.exists():
            error_msg = f"VCF parsing failed: file not found: {vcf_path}"
            logger.error(error_msg)
            return {"status": "failed", "error": error_msg}
        # Parse VCF using basic parser
        try:
            # Attempt basic VCF reading
            vcf_data = _parse_vcf_basic(vcf_file)
        except Exception as e:
            error_msg = f"VCF parsing failed: {e}"
            logger.error(error_msg)
            return {"status": "failed", "error": error_msg}
    elif isinstance(vcf_path, dict):
        vcf_data = vcf_path
    else:
        return {"status": "failed", "error": "Invalid input type for vcf_path"}

    if not vcf_data or "genotypes" not in vcf_data:
        logger.error("No genotype data provided")
        return {"status": "failed", "error": "No genotype data in VCF"}

    genotypes = vcf_data["genotypes"]
    if not genotypes:
        logger.error("Empty genotype data")
        return {"status": "failed", "error": "Empty genotype data"}

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

    plt.close(fig)
    return {"status": "success", "n_samples": n_samples, "n_variants": n_variants}


def _parse_vcf_basic(vcf_file: Path) -> Dict[str, Any]:
    """Basic VCF file parser for missingness analysis.

    Args:
        vcf_file: Path to VCF file

    Returns:
        Dict with genotypes data
    """
    genotypes: List[List[Optional[int]]] = []
    with open(vcf_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                continue
            # Extract genotype fields (columns 9+)
            gt_fields = parts[9:]
            row: List[Optional[int]] = []
            for gt in gt_fields:
                gt_val = gt.split(":")[0]
                if gt_val in ("./.", ".|.", "."):
                    row.append(None)
                else:
                    try:
                        alleles = gt_val.replace("|", "/").split("/")
                        row.append(sum(int(a) for a in alleles))
                    except (ValueError, TypeError):
                        row.append(None)
            genotypes.append(row)
    return {"genotypes": genotypes}


def transition_transversion_plot(
    variants: Any, output_file: Optional[str | Path] = None, title: str = "Transition/Transversion Plot"
) -> Dict[str, Any]:
    """Create a plot showing transition vs transversion ratios.

    Args:
        variants: List of variant dicts with "REF" and "ALT" keys,
                  or dict with "transitions"/"transversions" arrays
        output_file: Optional output file path
        title: Plot title

    Returns:
        Dict with status and metadata including ts_tv_ratio
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for transition/transversion plot")
        return {"status": "failed", "error": "matplotlib not available"}

    if not variants:
        logger.error("No transition/transversion data provided")
        return {"status": "failed", "error": "No data provided"}

    # Determine input format
    transitions_list: List[int] = []
    transversions_list: List[int] = []
    positions: List[int] = []

    # Define transition pairs
    transition_pairs = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    if isinstance(variants, (list, tuple)) and len(variants) > 0 and isinstance(variants[0], dict):
        # Check if it's list-of-dicts with REF/ALT
        if "REF" in variants[0] and "ALT" in variants[0]:
            # Compute Ts/Tv from individual variant records
            ts_count = 0
            tv_count = 0
            for v in variants:
                ref = v.get("REF", "")
                alt = v.get("ALT", "")
                if len(ref) == 1 and len(alt) == 1:
                    if (ref, alt) in transition_pairs:
                        ts_count += 1
                    else:
                        tv_count += 1

            # Create a simple bar chart for Ts vs Tv
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

            # Plot 1: Bar chart of Ts vs Tv counts
            categories = ["Transitions", "Transversions"]
            counts = [ts_count, tv_count]
            colors = ["steelblue", "coral"]
            bars = ax1.bar(categories, counts, color=colors, alpha=0.8, edgecolor="black")
            ax1.set_ylabel("Count")
            ax1.set_title("Transition vs Transversion Counts")
            ax1.grid(True, alpha=0.3, axis="y")

            # Add count labels
            for bar, count in zip(bars, counts):
                ax1.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max(counts) * 0.01,
                    str(count),
                    ha="center",
                    va="bottom",
                    fontsize=12,
                    fontweight="bold",
                )

            # Plot 2: Substitution matrix
            bases = ["A", "C", "G", "T"]
            sub_matrix = np.zeros((4, 4))
            base_idx = {b: i for i, b in enumerate(bases)}
            for v in variants:
                ref = v.get("REF", "")
                alt = v.get("ALT", "")
                if ref in base_idx and alt in base_idx:
                    sub_matrix[base_idx[ref]][base_idx[alt]] += 1

            im = ax2.imshow(sub_matrix, cmap="YlOrRd", aspect="auto")
            ax2.set_xticks(range(4))
            ax2.set_yticks(range(4))
            ax2.set_xticklabels(bases)
            ax2.set_yticklabels(bases)
            ax2.set_xlabel("ALT")
            ax2.set_ylabel("REF")
            ax2.set_title("Substitution Matrix")
            plt.colorbar(im, ax=ax2)

            # Add text annotations
            for i in range(4):
                for j in range(4):
                    val = int(sub_matrix[i][j])
                    if val > 0:
                        ax2.text(j, i, str(val), ha="center", va="center", fontsize=10)

            ts_tv_ratio = ts_count / tv_count if tv_count > 0 else 0.0

            fig.suptitle(f"{title} (Ts/Tv = {ts_tv_ratio:.2f})", fontsize=14)
            plt.tight_layout()

            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches="tight")
                logger.info(f"Saved transition/transversion plot to {output_file}")

            plt.close(fig)
            return {
                "status": "success",
                "ts_tv_ratio": ts_tv_ratio,
                "transitions": ts_count,
                "transversions": tv_count,
                "n_variants": len(variants),
            }

    # Legacy dict input format
    if isinstance(variants, dict):
        transitions_list = variants.get("transitions", [])
        transversions_list = variants.get("transversions", [])
        positions = variants.get("positions", [])

        if not transitions_list or not transversions_list:
            logger.error("Transition/transversion data must include 'transitions' and 'transversions' arrays")
            return {"status": "failed", "error": "Missing transitions/transversions data"}

        if not positions:
            positions = list(range(1, len(transitions_list) + 1))

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        # Plot 1: Transitions and transversions along genome
        ax1.plot(positions, transitions_list, "b-", label="Transitions", alpha=0.7)
        ax1.plot(positions, transversions_list, "r-", label="Transversions", alpha=0.7)
        ax1.set_xlabel("Genomic Position")
        ax1.set_ylabel("Count")
        ax1.set_title("Transitions vs Transversions")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot 2: Ts/Tv ratio
        ts_tv_ratios = []
        for ts, tv in zip(transitions_list, transversions_list):
            if tv > 0:
                ts_tv_ratios.append(ts / tv)
            else:
                ts_tv_ratios.append(0)

        ax2.plot(positions, ts_tv_ratios, "g-", label="Ts/Tv Ratio")
        ax2.axhline(y=np.mean(ts_tv_ratios), color="red", linestyle="--", label=f"Mean: {np.mean(ts_tv_ratios):.2f}")
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

        mean_ratio = float(np.mean(ts_tv_ratios))
        plt.close(fig)
        return {"status": "success", "ts_tv_ratio": mean_ratio}

    logger.error("Invalid input format for transition/transversion plot")
    return {"status": "failed", "error": "Invalid input format"}
