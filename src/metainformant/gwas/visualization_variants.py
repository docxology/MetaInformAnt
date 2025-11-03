"""Variant property visualizations for GWAS.

This module provides plots for variant characteristics including
allele frequencies, density, quality metrics, and Hardy-Weinberg equilibrium.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg", force=True)

from ..core.io import ensure_directory

logger = logging.getLogger(__name__)

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    pd = None


def maf_distribution(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    bins: int = 50,
    title: str | None = None,
) -> dict[str, Any]:
    """Minor allele frequency distribution (allele frequency spectrum).
    
    Shows distribution of allele frequencies, useful for assessing
    population history and QC filtering effects.
    
    Args:
        results: Association results or path
        output_path: Output path
        bins: Number of histogram bins
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info("maf_distribution: Plotting allele frequency spectrum")
    
    # Load data
    if isinstance(results, (Path, str)):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(results, sep="\t")
            results_list = df.to_dict("records")
        else:
            from ..core.io import read_tsv
            data = read_tsv(results)
            header = data[0]
            results_list = [{header[i]: row[i] for i in range(len(header))} for row in data[1:]]
    else:
        results_list = results
    
    # Extract MAF
    mafs = []
    for r in results_list:
        try:
            maf = float(r.get("MAF", r.get("maf", 0.0)))
            if 0 <= maf <= 0.5:
                mafs.append(maf)
        except (ValueError, TypeError):
            continue
    
    if not mafs:
        return {"status": "failed", "error": "No valid MAF values"}
    
    mafs = np.array(mafs)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Histogram
    counts, bin_edges, patches = ax.hist(mafs, bins=bins, color='#2E86AB',
                                         alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # Color rare variants differently
    rare_threshold = 0.05
    for i, patch in enumerate(patches):
        if bin_edges[i] < rare_threshold:
            patch.set_facecolor('#C73E1D')
    
    # Vertical line at rare variant threshold
    ax.axvline(x=rare_threshold, color='red', linestyle='--', linewidth=2,
              label=f'Rare variants (MAF < {rare_threshold})')
    
    # Labels
    ax.set_xlabel("Minor Allele Frequency (MAF)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Number of Variants", fontsize=12, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        n_rare = (mafs < rare_threshold).sum()
        ax.set_title(f"Allele Frequency Spectrum ({len(mafs):,} variants, {n_rare:,} rare)",
                    fontsize=14, fontweight="bold")
    
    ax.grid(True, alpha=0.3, axis='y', linestyle=':')
    ax.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"maf_distribution: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_variants": len(mafs),
        "num_rare": int((mafs < rare_threshold).sum()),
        "mean_maf": float(mafs.mean()),
        "median_maf": float(np.median(mafs)),
    }


def variant_density_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    window_size: int = 1_000_000,
    title: str | None = None,
) -> dict[str, Any]:
    """SNP density across genome.
    
    Shows variant count per window, useful for identifying regions
    with high/low variant density (e.g., centromeres, genes).
    
    Args:
        results: Association results or path
        output_path: Output path
        window_size: Window size in bp (default 1 Mb)
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"variant_density_plot: Window size {window_size} bp")
    
    # Load data
    if isinstance(results, (Path, str)):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(results, sep="\t")
            results_list = df.to_dict("records")
        else:
            from ..core.io import read_tsv
            data = read_tsv(results)
            header = data[0]
            results_list = [{header[i]: row[i] for i in range(len(header))} for row in data[1:]]
    else:
        results_list = results
    
    # Group by chromosome and position
    chrom_data = {}
    for r in results_list:
        try:
            chrom = str(r.get("CHROM", ""))
            pos = int(r.get("POS", 0))
            
            if chrom not in chrom_data:
                chrom_data[chrom] = []
            chrom_data[chrom].append(pos)
        except (ValueError, TypeError):
            continue
    
    if not chrom_data:
        return {"status": "failed", "error": "No valid data"}
    
    # Calculate density per chromosome
    fig, ax = plt.subplots(figsize=(14, 6))
    
    unique_chroms = sorted(chrom_data.keys())
    cumulative_pos = 0
    chrom_starts = {}
    
    for chrom in unique_chroms:
        positions = sorted(chrom_data[chrom])
        if not positions:
            continue
        
        chrom_starts[chrom] = cumulative_pos
        
        # Bin positions
        max_pos = max(positions)
        bins = np.arange(0, max_pos + window_size, window_size)
        counts, _ = np.histogram(positions, bins=bins)
        
        # Plot density
        bin_centers = bins[:-1] + window_size / 2
        x_coords = cumulative_pos + bin_centers
        
        ax.plot(x_coords, counts, linewidth=1, alpha=0.8)
        
        cumulative_pos += max_pos + 10_000_000
    
    # Chromosome separators
    for chrom, start_pos in chrom_starts.items():
        ax.axvline(x=start_pos, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
    
    # Labels
    ax.set_xlabel("Genomic Position", fontsize=12, fontweight="bold")
    ax.set_ylabel(f"Variants per {window_size/1e6:.1f} Mb", fontsize=12, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        total_variants = sum(len(chrom_data[c]) for c in unique_chroms)
        ax.set_title(f"Variant Density Across Genome ({total_variants:,} variants)",
                    fontsize=14, fontweight="bold")
    
    # X-axis: chromosome labels
    chrom_centers = {c: chrom_starts[c] + max(chrom_data[c]) / 2 
                    for c in unique_chroms if chrom_data[c]}
    ax.set_xticks(list(chrom_centers.values()))
    ax.set_xticklabels([c.replace('chr', '').replace('NC_', '').split('.')[0] 
                        for c in chrom_centers.keys()], rotation=45, ha='right')
    
    ax.grid(True, alpha=0.3, axis='y', linestyle=':')
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"variant_density_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_chromosomes": len(unique_chroms),
        "total_variants": sum(len(chrom_data[c]) for c in unique_chroms),
    }


def hwe_deviation_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    threshold: float = 1e-6,
    title: str | None = None,
) -> dict[str, Any]:
    """Hardy-Weinberg equilibrium deviation plot.
    
    Shows distribution of HWE p-values. Deviations can indicate
    genotyping errors or population substructure.
    
    Args:
        results: Association results or path with HWE p-values
        output_path: Output path
        threshold: HWE significance threshold
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info("hwe_deviation_plot: Plotting HWE deviations")
    
    # Load data
    if isinstance(results, (Path, str)):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(results, sep="\t")
            results_list = df.to_dict("records")
        else:
            from ..core.io import read_tsv
            data = read_tsv(results)
            header = data[0]
            results_list = [{header[i]: row[i] for i in range(len(header))} for row in data[1:]]
    else:
        results_list = results
    
    # Extract HWE p-values (if available)
    hwe_pvals = []
    for r in results_list:
        try:
            hwe_p = float(r.get("hwe_pval", r.get("HWE_P", 1.0)))
            if 0 < hwe_p <= 1:
                hwe_pvals.append(hwe_p)
        except (ValueError, TypeError):
            continue
    
    if not hwe_pvals:
        return {
            "status": "skipped",
            "message": "No HWE p-values found in results",
            "note": "HWE calculation requires genotype data",
        }
    
    hwe_pvals = np.array(hwe_pvals)
    hwe_log = -np.log10(hwe_pvals)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Histogram
    ax.hist(hwe_log, bins=50, color='#2E86AB', alpha=0.7,
           edgecolor='black', linewidth=0.5)
    
    # Threshold line
    threshold_log = -math.log10(threshold)
    ax.axvline(x=threshold_log, color='red', linestyle='--', linewidth=2,
              label=f'HWE threshold (p={threshold:.0e})')
    
    # Labels
    ax.set_xlabel("-log₁₀(HWE p-value)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Number of Variants", fontsize=12, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        n_deviate = (hwe_pvals < threshold).sum()
        ax.set_title(f"Hardy-Weinberg Equilibrium Deviations ({n_deviate:,} deviant)",
                    fontsize=14, fontweight="bold")
    
    ax.grid(True, alpha=0.3, axis='y', linestyle=':')
    ax.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"hwe_deviation_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_variants": len(hwe_pvals),
        "num_deviant": int((hwe_pvals < threshold).sum()),
    }


def missingness_plot(
    vcf_path: Path,
    output_path: str | Path,
    *,
    by_sample: bool = True,
    title: str | None = None,
) -> dict[str, Any]:
    """Missing data patterns (by sample or by variant).
    
    Visualizes missingness rates to identify problematic samples
    or variants with poor genotyping.
    
    Args:
        vcf_path: Path to VCF file
        output_path: Output path
        by_sample: Plot by sample (True) or by variant (False)
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"missingness_plot: By {'sample' if by_sample else 'variant'}")
    
    return {
        "status": "skipped",
        "message": "Missingness calculation requires VCF parsing",
        "recommendation": "Use PLINK --missing or bcftools +missing2ref",
        "note": "Future enhancement: integrate VCF missingness parser",
    }


def transition_transversion_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    title: str | None = None,
) -> dict[str, Any]:
    """Transition/transversion ratio plot.
    
    Ts/Tv ratio is a quality metric for variant calling.
    Expected: ~2.0-2.1 for whole genome, higher for exomes.
    
    Args:
        results: Association results or path with REF/ALT alleles
        output_path: Output path
        title: Plot title
    
    Returns:
        Plot metadata with Ts/Tv ratio
    """
    logger.info("transition_transversion_plot: Calculating Ts/Tv ratio")
    
    # Load data
    if isinstance(results, (Path, str)):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(results, sep="\t")
            results_list = df.to_dict("records")
        else:
            from ..core.io import read_tsv
            data = read_tsv(results)
            header = data[0]
            results_list = [{header[i]: row[i] for i in range(len(header))} for row in data[1:]]
    else:
        results_list = results
    
    # Classify variants
    transitions = 0
    transversions = 0
    
    ts_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    
    for r in results_list:
        try:
            ref = str(r.get("REF", r.get("ref", ""))).upper()
            alt = str(r.get("ALT", r.get("alt", ""))).upper()
            
            if len(ref) == 1 and len(alt) == 1:
                if (ref, alt) in ts_pairs:
                    transitions += 1
                else:
                    transversions += 1
        except (ValueError, TypeError):
            continue
    
    if transitions + transversions == 0:
        return {
            "status": "skipped",
            "message": "No REF/ALT alleles found in results",
            "note": "Ensure results include REF and ALT columns",
        }
    
    ts_tv_ratio = transitions / transversions if transversions > 0 else 0
    
    # Create plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    categories = ['Transitions', 'Transversions']
    counts = [transitions, transversions]
    colors = ['#2E86AB', '#C73E1D']
    
    bars = ax.bar(categories, counts, color=colors, alpha=0.7,
                  edgecolor='black', linewidth=1.5)
    
    # Add count labels on bars
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{count:,}', ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    # Labels
    ax.set_ylabel("Number of Variants", fontsize=12, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"Ts/Tv Ratio = {ts_tv_ratio:.3f}", fontsize=14, fontweight="bold")
    
    # Add text box with interpretation
    interpretation = "Good" if 1.8 <= ts_tv_ratio <= 2.5 else "Check quality"
    textstr = f'Ts/Tv = {ts_tv_ratio:.3f}\n({interpretation})'
    ax.text(0.95, 0.95, textstr, transform=ax.transAxes,
           fontsize=11, verticalalignment='top', horizontalalignment='right',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax.grid(True, alpha=0.3, axis='y', linestyle=':')
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"transition_transversion_plot: Saved to {output_path_obj} (Ts/Tv={ts_tv_ratio:.3f})")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "transitions": transitions,
        "transversions": transversions,
        "ts_tv_ratio": float(ts_tv_ratio),
        "quality_status": "good" if 1.8 <= ts_tv_ratio <= 2.5 else "check",
    }


