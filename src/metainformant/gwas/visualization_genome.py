"""Genome-wide visualization for GWAS results.

This module provides chromosome-scale and genome-wide views of association results,
optimized for millions of SNPs with efficient rendering and publication-quality output.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Wedge, Rectangle

matplotlib.use("Agg", force=True)

from ..core.io import ensure_directory
from ..core.logging import get_logger

logger = get_logger(__name__)

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    pd = None


def manhattan_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    significance_threshold: float = 5e-8,
    suggestiveness_threshold: float = 1e-5,
    chrom_colors: list[str] | None = None,
    title: str | None = None,
    max_points_per_chrom: int | None = 50000,
    point_size: float = 15,
) -> dict[str, Any]:
    """Manhattan plot optimized for genome-scale data (millions of SNPs).
    
    Implements intelligent point thinning:
    - Keeps all significant/suggestive variants
    - Samples non-significant variants to reduce rendering time
    - Maintains visual appearance while improving performance
    
    Args:
        results: Association results or path to TSV file
        output_path: Output path for plot
        significance_threshold: Genome-wide significance (default 5e-8)
        suggestiveness_threshold: Suggestive threshold (default 1e-5)
        chrom_colors: List of colors for chromosomes
        title: Plot title
        max_points_per_chrom: Max points per chromosome for thinning
        point_size: Marker size
    
    Returns:
        Plot metadata with statistics
    """
    logger.info("manhattan_plot: Generating genome-wide Manhattan plot")
    
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
    
    if not results_list:
        return {"status": "failed", "error": "No results"}
    
    # Extract data
    chroms, positions, pvalues = [], [], []
    for r in results_list:
        try:
            chrom = str(r.get("CHROM", ""))
            pos = int(r.get("POS", 0))
            pval = float(r.get("p_value", 1.0))
            if pval > 0:
                chroms.append(chrom)
                positions.append(pos)
                pvalues.append(pval)
        except (ValueError, TypeError):
            continue
    
    if not pvalues:
        return {"status": "failed", "error": "No valid p-values"}
    
    logger.info(f"manhattan_plot: Processing {len(pvalues):,} variants")
    
    # Group by chromosome
    unique_chroms = sorted(set(chroms), key=lambda x: (
        int(''.join(filter(str.isdigit, x))) if any(c.isdigit() for c in x) else 999, x
    ))
    
    chrom_data = {c: {"pos": [], "pval": []} for c in unique_chroms}
    for chrom, pos, pval in zip(chroms, positions, pvalues):
        if chrom in chrom_data:
            chrom_data[chrom]["pos"].append(pos)
            chrom_data[chrom]["pval"].append(pval)
    
    # Calculate cumulative positions
    cumulative_pos, chrom_cumulative, chrom_spans = 0, {}, {}
    for chrom in unique_chroms:
        chrom_cumulative[chrom] = cumulative_pos
        if chrom_data[chrom]["pos"]:
            min_pos = min(chrom_data[chrom]["pos"])
            max_pos = max(chrom_data[chrom]["pos"])
            span = max_pos - min_pos
            chrom_spans[chrom] = (cumulative_pos, cumulative_pos + span)
            cumulative_pos += span + 10_000_000
    
    # Plot
    fig, ax = plt.subplots(figsize=(18, 7))
    
    if chrom_colors is None:
        chrom_colors = ["#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#06A77D", "#6A4C93"]
    
    total_plotted = 0
    for idx, chrom in enumerate(unique_chroms):
        if not chrom_data[chrom]["pos"]:
            continue
        
        pos_arr = np.array(chrom_data[chrom]["pos"])
        pval_arr = np.array(chrom_data[chrom]["pval"])
        
        x_pos = chrom_cumulative[chrom] + pos_arr - pos_arr.min()
        y_neg_log = -np.log10(pval_arr)
        
        # Intelligent thinning for performance
        if max_points_per_chrom and len(x_pos) > max_points_per_chrom:
            sig_mask = pval_arr < suggestiveness_threshold
            sig_idx = np.where(sig_mask)[0]
            non_sig_idx = np.where(~sig_mask)[0]
            
            if len(non_sig_idx) > (max_points_per_chrom - len(sig_idx)):
                sampled = np.random.choice(non_sig_idx, 
                                          size=max_points_per_chrom - len(sig_idx),
                                          replace=False)
                keep_idx = np.sort(np.concatenate([sig_idx, sampled]))
            else:
                keep_idx = np.arange(len(x_pos))
            
            x_pos, y_neg_log = x_pos[keep_idx], y_neg_log[keep_idx]
        
        total_plotted += len(x_pos)
        color = chrom_colors[idx % len(chrom_colors)]
        ax.scatter(x_pos, y_neg_log, c=color, s=point_size, alpha=0.6,
                  label=chrom.replace('chr', '').replace('NC_', '').split('.')[0] if idx < 10 else None,
                  rasterized=True)
    
    # Significance lines
    sig_line = -math.log10(significance_threshold)
    sugg_line = -math.log10(suggestiveness_threshold)
    
    ax.axhline(y=sig_line, color="red", linestyle="--", linewidth=1.5,
              label=f"Genome-wide ({significance_threshold:.0e})", zorder=100)
    ax.axhline(y=sugg_line, color="orange", linestyle="--", linewidth=1,
              label=f"Suggestive ({suggestiveness_threshold:.0e})", zorder=100)
    
    # Labels
    ax.set_xlabel("Chromosome", fontsize=14, fontweight="bold")
    ax.set_ylabel("-log₁₀(p-value)", fontsize=14, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=16, fontweight="bold")
    else:
        ax.set_title(f"Manhattan Plot ({len(pvalues):,} variants, {total_plotted:,} plotted)",
                    fontsize=16, fontweight="bold")
    
    # X-axis ticks
    chrom_centers = [(chrom_spans[c][0] + chrom_spans[c][1]) / 2 
                     for c in unique_chroms if c in chrom_spans]
    chrom_labels = [c.replace('chr', '').replace('NC_', '').split('.')[0] for c in unique_chroms]
    
    ax.set_xticks(chrom_centers)
    ax.set_xticklabels(chrom_labels, rotation=45, ha="right", fontsize=9)
    
    ax.grid(True, alpha=0.2, axis="y", linestyle=":")
    ax.legend(loc="upper right", fontsize=8, ncol=3, framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"manhattan_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_variants": len(pvalues),
        "num_plotted": total_plotted,
        "num_significant": sum(1 for p in pvalues if p < significance_threshold),
        "num_suggestive": sum(1 for p in pvalues if p < suggestiveness_threshold),
    }


def circular_manhattan_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    significance_threshold: float = 5e-8,
    title: str | None = None,
) -> dict[str, Any]:
    """Circular Manhattan plot for genome-wide view.
    
    Arranges chromosomes in a circle with -log10(p) as radial distance.
    Useful for visualizing genome-wide patterns and cross-chromosome effects.
    
    Args:
        results: Association results or path to file
        output_path: Output path
        significance_threshold: Significance threshold
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info("circular_manhattan_plot: Generating circular plot")
    
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
    
    # Extract and group data
    chroms, positions, pvalues = [], [], []
    for r in results_list:
        try:
            chrom = str(r.get("CHROM", ""))
            pos = int(r.get("POS", 0))
            pval = float(r.get("p_value", 1.0))
            if pval > 0:
                chroms.append(chrom)
                positions.append(pos)
                pvalues.append(-math.log10(pval))
        except (ValueError, TypeError):
            continue
    
    if not pvalues:
        return {"status": "failed", "error": "No valid p-values"}
    
    unique_chroms = sorted(set(chroms))
    chrom_data = {c: {"pos": [], "pval": []} for c in unique_chroms}
    
    for chrom, pos, pval in zip(chroms, positions, pvalues):
        chrom_data[chrom]["pos"].append(pos)
        chrom_data[chrom]["pval"].append(pval)
    
    # Create circular plot
    fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(projection='polar'))
    
    n_chroms = len(unique_chroms)
    theta_span = 2 * np.pi / n_chroms
    sig_threshold = -math.log10(significance_threshold)
    
    colors = plt.cm.tab20(np.linspace(0, 1, n_chroms))
    
    for idx, chrom in enumerate(unique_chroms):
        if not chrom_data[chrom]["pos"]:
            continue
        
        # Normalize positions within chromosome
        pos_norm = np.array(chrom_data[chrom]["pos"])
        pos_norm = (pos_norm - pos_norm.min()) / (pos_norm.max() - pos_norm.min() + 1)
        
        # Convert to theta (angle)
        theta = idx * theta_span + pos_norm * theta_span
        r = np.array(chrom_data[chrom]["pval"])
        
        # Plot points
        ax.scatter(theta, r, c=[colors[idx]], s=5, alpha=0.6, rasterized=True)
    
    # Significance circle
    theta_full = np.linspace(0, 2*np.pi, 100)
    r_sig = np.full_like(theta_full, sig_threshold)
    ax.plot(theta_full, r_sig, 'r--', linewidth=2, label='Significance', zorder=100)
    
    # Labels
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_xticks([idx * theta_span for idx in range(n_chroms)])
    ax.set_xticklabels([c.replace('chr', '').replace('NC_', '').split('.')[0] 
                        for c in unique_chroms], fontsize=8)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold", pad=20)
    else:
        ax.set_title(f"Circular Manhattan Plot ({len(positions):,} variants)",
                    fontsize=14, fontweight="bold", pad=20)
    
    ax.legend(loc="upper right", fontsize=10, bbox_to_anchor=(1.3, 1.1))
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"circular_manhattan_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_variants": len(positions),
    }


def chromosome_ideogram(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    significance_threshold: float = 5e-8,
    title: str | None = None,
) -> dict[str, Any]:
    """Chromosome ideogram with significant loci marked.
    
    Shows chromosome structure with bands indicating significant variants.
    Useful for visualizing distribution of associations across genome.
    
    Args:
        results: Association results or path
        output_path: Output path
        significance_threshold: Significance threshold
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info("chromosome_ideogram: Generating ideogram")
    
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
    
    # Extract data
    chrom_data = {}
    for r in results_list:
        try:
            chrom = str(r.get("CHROM", ""))
            pos = int(r.get("POS", 0))
            pval = float(r.get("p_value", 1.0))
            
            if chrom not in chrom_data:
                chrom_data[chrom] = {"positions": [], "pvalues": [], "sig_positions": []}
            
            chrom_data[chrom]["positions"].append(pos)
            chrom_data[chrom]["pvalues"].append(pval)
            
            if pval < significance_threshold:
                chrom_data[chrom]["sig_positions"].append(pos)
        except (ValueError, TypeError):
            continue
    
    if not chrom_data:
        return {"status": "failed", "error": "No valid data"}
    
    # Create plot
    unique_chroms = sorted(chrom_data.keys())
    n_chroms = len(unique_chroms)
    
    fig, ax = plt.subplots(figsize=(14, max(8, n_chroms * 0.5)))
    
    for idx, chrom in enumerate(unique_chroms):
        y_pos = n_chroms - idx - 1
        positions = chrom_data[chrom]["positions"]
        sig_pos = chrom_data[chrom]["sig_positions"]
        
        if not positions:
            continue
        
        chrom_len = max(positions)
        
        # Draw chromosome bar
        rect = Rectangle((0, y_pos - 0.3), chrom_len, 0.6,
                        linewidth=1, edgecolor='black',
                        facecolor='lightgray', alpha=0.7)
        ax.add_patch(rect)
        
        # Mark significant loci
        for pos in sig_pos:
            ax.plot([pos, pos], [y_pos - 0.3, y_pos + 0.3],
                   'r-', linewidth=2, alpha=0.8)
        
        # Chromosome label
        chrom_label = chrom.replace('chr', '').replace('NC_', '').split('.')[0]
        ax.text(-chrom_len * 0.05, y_pos, chrom_label,
               ha='right', va='center', fontsize=10, fontweight='bold')
        
        # Significant loci count
        if sig_pos:
            ax.text(chrom_len * 1.02, y_pos, f"n={len(sig_pos)}",
                   ha='left', va='center', fontsize=8, color='red')
    
    # Labels
    ax.set_xlim(-max([max(chrom_data[c]["positions"]) for c in unique_chroms]) * 0.1,
                max([max(chrom_data[c]["positions"]) for c in unique_chroms]) * 1.15)
    ax.set_ylim(-1, n_chroms)
    ax.set_xlabel("Position (bp)", fontsize=12, fontweight="bold")
    ax.set_yticks([])
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        total_sig = sum(len(chrom_data[c]["sig_positions"]) for c in unique_chroms)
        ax.set_title(f"Chromosome Ideogram ({total_sig} significant loci)",
                    fontsize=14, fontweight="bold")
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"chromosome_ideogram: Saved to {output_path_obj}")
    
    total_significant = sum(len(chrom_data[c]["sig_positions"]) for c in unique_chroms)
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_chromosomes": n_chroms,
        "num_significant_loci": total_significant,
    }


def genome_wide_ld_heatmap(
    vcf_path: Path,
    output_path: str | Path,
    *,
    sample_size: int = 1000,
    title: str | None = None,
) -> dict[str, Any]:
    """Genome-wide linkage disequilibrium heatmap (downsampled).
    
    Shows LD patterns across genome by sampling variants.
    Full genome-wide LD computation is prohibitive for millions of SNPs.
    
    Args:
        vcf_path: Path to VCF file
        output_path: Output path
        sample_size: Number of variants to sample
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"genome_wide_ld_heatmap: Sampling {sample_size} variants for LD")
    
    return {
        "status": "skipped",
        "message": "Genome-wide LD requires specialized tools (PLINK, LDmatrix)",
        "recommendation": "Use regional LD plots for specific loci instead",
    }






