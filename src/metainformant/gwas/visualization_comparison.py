"""Multi-trait and cross-cohort comparison visualizations for GWAS.

This module provides Miami plots, multi-trait Manhattan plots,
meta-analysis forest plots, and replication concordance plots.
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


def miami_plot(
    results1: list[dict[str, Any]] | Path,
    results2: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    trait1_name: str = "Trait 1",
    trait2_name: str = "Trait 2",
    significance_threshold: float = 5e-8,
    title: str | None = None,
) -> dict[str, Any]:
    """Miami plot (back-to-back Manhattan) for two traits.
    
    Shows GWAS results for two phenotypes in a mirrored Manhattan plot,
    useful for identifying pleiotropic loci.
    
    Args:
        results1: Association results for trait 1
        results2: Association results for trait 2
        output_path: Output path
        trait1_name: Name of trait 1
        trait2_name: Name of trait 2
        significance_threshold: P-value threshold
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"miami_plot: {trait1_name} vs {trait2_name}")
    
    # Load data for both traits
    def load_results(results):
        """Load GWAS results from file or dict.
        
        Args:
            results: Path to results file or dict of results
            
        Returns:
            List of result records
        """
        if isinstance(results, (Path, str)):
            if PANDAS_AVAILABLE:
                df = pd.read_csv(results, sep="\t")
                return df.to_dict("records")
            else:
                from ..core.io import read_tsv
                data = read_tsv(results)
                header = data[0]
                return [{header[i]: row[i] for i in range(len(header))} for row in data[1:]]
        return results
    
    results1_list = load_results(results1)
    results2_list = load_results(results2)
    
    # Extract data for trait 1
    chroms1, pos1, pvals1 = [], [], []
    for r in results1_list:
        try:
            chrom = str(r.get("CHROM", ""))
            pos = int(r.get("POS", 0))
            pval = float(r.get("p_value", 1.0))
            if pval > 0:
                chroms1.append(chrom)
                pos1.append(pos)
                pvals1.append(-math.log10(pval))
        except (ValueError, TypeError):
            continue
    
    # Extract data for trait 2
    chroms2, pos2, pvals2 = [], [], []
    for r in results2_list:
        try:
            chrom = str(r.get("CHROM", ""))
            pos = int(r.get("POS", 0))
            pval = float(r.get("p_value", 1.0))
            if pval > 0:
                chroms2.append(chrom)
                pos2.append(pos)
                pvals2.append(-math.log10(pval))
        except (ValueError, TypeError):
            continue
    
    if not pvals1 or not pvals2:
        return {"status": "failed", "error": "Insufficient data for both traits"}
    
    # Get unique chromosomes
    unique_chroms = sorted(set(chroms1 + chroms2))
    
    # Calculate cumulative positions
    chrom_cumulative = {}
    cumulative_pos = 0
    for chrom in unique_chroms:
        chrom_cumulative[chrom] = cumulative_pos
        # Estimate chromosome length
        chrom_pos1 = [p for c, p in zip(chroms1, pos1) if c == chrom]
        chrom_pos2 = [p for c, p in zip(chroms2, pos2) if c == chrom]
        if chrom_pos1 or chrom_pos2:
            max_pos = max((max(chrom_pos1) if chrom_pos1 else 0),
                         (max(chrom_pos2) if chrom_pos2 else 0))
            cumulative_pos += max_pos + 10_000_000
    
    # Create plot
    fig, ax = plt.subplots(figsize=(18, 10))
    
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']
    
    # Plot trait 1 (top, positive values)
    for idx, chrom in enumerate(unique_chroms):
        chrom_data1 = [(p, v) for c, p, v in zip(chroms1, pos1, pvals1) if c == chrom]
        if not chrom_data1:
            continue
        
        positions, values = zip(*chrom_data1)
        x_pos = chrom_cumulative[chrom] + np.array(positions)
        
        ax.scatter(x_pos, values, c=colors[idx % len(colors)],
                  s=15, alpha=0.6, rasterized=True)
    
    # Plot trait 2 (bottom, negative values)
    for idx, chrom in enumerate(unique_chroms):
        chrom_data2 = [(p, v) for c, p, v in zip(chroms2, pos2, pvals2) if c == chrom]
        if not chrom_data2:
            continue
        
        positions, values = zip(*chrom_data2)
        x_pos = chrom_cumulative[chrom] + np.array(positions)
        
        ax.scatter(x_pos, -np.array(values), c=colors[idx % len(colors)],
                  s=15, alpha=0.6, rasterized=True)
    
    # Significance lines
    sig_line = -math.log10(significance_threshold)
    ax.axhline(y=sig_line, color='red', linestyle='--', linewidth=1.5, zorder=100)
    ax.axhline(y=-sig_line, color='red', linestyle='--', linewidth=1.5, zorder=100)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1, zorder=100)
    
    # Labels
    ax.set_xlabel("Chromosome", fontsize=14, fontweight="bold")
    ax.set_ylabel(f"{trait1_name} ← -log₁₀(p) → {trait2_name}",
                 fontsize=14, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=16, fontweight="bold")
    else:
        ax.set_title(f"Miami Plot: {trait1_name} vs {trait2_name}",
                    fontsize=16, fontweight="bold")
    
    # X-axis ticks (chromosome centers)
    chrom_centers = []
    for chrom in unique_chroms:
        chrom_pos = [p for c, p in zip(chroms1 + chroms2, pos1 + pos2) if c == chrom]
        if chrom_pos:
            center = chrom_cumulative[chrom] + (max(chrom_pos) - min(chrom_pos)) / 2
            chrom_centers.append(center)
    
    ax.set_xticks(chrom_centers)
    ax.set_xticklabels([c.replace('chr', '').replace('NC_', '').split('.')[0]
                        for c in unique_chroms], rotation=45, ha='right')
    
    ax.grid(True, alpha=0.2, axis='y', linestyle=':')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"miami_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "trait1_variants": len(pvals1),
        "trait2_variants": len(pvals2),
    }


def multi_trait_manhattan(
    results_dict: dict[str, list[dict[str, Any]] | Path],
    output_path: str | Path,
    *,
    significance_threshold: float = 5e-8,
    title: str | None = None,
) -> dict[str, Any]:
    """Multi-trait Manhattan plot (stacked or faceted).
    
    Shows GWAS results for multiple phenotypes simultaneously.
    
    Args:
        results_dict: Dictionary of {trait_name: results}
        output_path: Output path
        significance_threshold: P-value threshold
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"multi_trait_manhattan: {len(results_dict)} traits")
    
    n_traits = len(results_dict)
    if n_traits < 2:
        return {"status": "failed", "error": "Need at least 2 traits"}
    
    # Create faceted plot
    fig, axes = plt.subplots(n_traits, 1, figsize=(16, 4 * n_traits),
                            sharex=True)
    
    if n_traits == 1:
        axes = [axes]
    
    sig_line = -math.log10(significance_threshold)
    
    for idx, (trait_name, results) in enumerate(results_dict.items()):
        ax = axes[idx]
        
        # Load results
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
        positions, pvalues = [], []
        for r in results_list:
            try:
                pos = int(r.get("POS", 0))
                pval = float(r.get("p_value", 1.0))
                if pval > 0:
                    positions.append(pos)
                    pvalues.append(-math.log10(pval))
            except (ValueError, TypeError):
                continue
        
        if positions:
            ax.scatter(positions, pvalues, s=10, alpha=0.6, c='#2E86AB',
                      rasterized=True)
            ax.axhline(y=sig_line, color='red', linestyle='--', linewidth=1)
        
        ax.set_ylabel(f"{trait_name}\n-log₁₀(p)", fontsize=10, fontweight="bold")
        ax.grid(True, alpha=0.3, linestyle=':')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    axes[-1].set_xlabel("Genomic Position", fontsize=12, fontweight="bold")
    
    if title:
        fig.suptitle(title, fontsize=14, fontweight="bold", y=0.995)
    else:
        fig.suptitle(f"Multi-Trait GWAS ({n_traits} traits)", fontsize=14, fontweight="bold")
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"multi_trait_manhattan: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_traits": n_traits,
    }


def cross_cohort_forest(
    cohort_results: dict[str, list[dict[str, Any]]],
    output_path: str | Path,
    *,
    variant_id: str,
    title: str | None = None,
) -> dict[str, Any]:
    """Forest plot for meta-analysis across cohorts.
    
    Shows effect estimates from multiple studies for a specific variant.
    
    Args:
        cohort_results: Dictionary of {cohort_name: results}
        output_path: Output path
        variant_id: Variant identifier (e.g., "chr1:12345")
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"cross_cohort_forest: Variant {variant_id} across cohorts")
    
    return {
        "status": "skipped",
        "message": "Meta-analysis requires multiple cohort datasets",
        "recommendation": "Use METAL or PLINK --meta for meta-analysis",
        "note": "Future enhancement: integrate meta-analysis tools",
    }


def concordance_plot(
    discovery_results: list[dict[str, Any]] | Path,
    replication_results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    significance_threshold: float = 5e-8,
    title: str | None = None,
) -> dict[str, Any]:
    """Replication concordance plot.
    
    Compares effect sizes/p-values between discovery and replication cohorts.
    Assesses consistency of associations.
    
    Args:
        discovery_results: Discovery cohort results
        replication_results: Replication cohort results
        output_path: Output path
        significance_threshold: P-value threshold
        title: Plot title
    
    Returns:
        Plot metadata with concordance statistics
    """
    logger.info("concordance_plot: Discovery vs replication")
    
    # Load data
    def load_results(results):
        """Load GWAS results from file or dict (concordance plot version).
        
        Args:
            results: Path to results file or dict of results
            
        Returns:
            List of result records
        """
        if isinstance(results, (Path, str)):
            if PANDAS_AVAILABLE:
                df = pd.read_csv(results, sep="\t")
                return df.to_dict("records")
            else:
                from ..core.io import read_tsv
                data = read_tsv(results)
                header = data[0]
                return [{header[i]: row[i] for i in range(len(header))} for row in data[1:]]
        return results
    
    disc_list = load_results(discovery_results)
    repl_list = load_results(replication_results)
    
    # Match variants by position
    disc_dict = {}
    for r in disc_list:
        try:
            variant_key = f"{r.get('CHROM', '')}:{r.get('POS', '')}"
            beta = float(r.get("beta", r.get("BETA", 0.0)))
            pval = float(r.get("p_value", 1.0))
            disc_dict[variant_key] = {"beta": beta, "pval": pval}
        except (ValueError, TypeError):
            continue
    
    # Match replication
    matched_disc_betas, matched_repl_betas = [], []
    for r in repl_list:
        try:
            variant_key = f"{r.get('CHROM', '')}:{r.get('POS', '')}"
            if variant_key in disc_dict:
                beta_repl = float(r.get("beta", r.get("BETA", 0.0)))
                beta_disc = disc_dict[variant_key]["beta"]
                matched_disc_betas.append(beta_disc)
                matched_repl_betas.append(beta_repl)
        except (ValueError, TypeError):
            continue
    
    if len(matched_disc_betas) < 2:
        return {
            "status": "failed",
            "error": "Insufficient overlapping variants between cohorts",
        }
    
    matched_disc_betas = np.array(matched_disc_betas)
    matched_repl_betas = np.array(matched_repl_betas)
    
    # Calculate concordance
    correlation = np.corrcoef(matched_disc_betas, matched_repl_betas)[0, 1]
    
    # Create plot
    fig, ax = plt.subplots(figsize=(9, 9))
    
    ax.scatter(matched_disc_betas, matched_repl_betas, s=20, alpha=0.6,
              c='#2E86AB', edgecolors='none', rasterized=True)
    
    # Identity line
    all_betas = np.concatenate([matched_disc_betas, matched_repl_betas])
    lim = max(abs(all_betas.min()), abs(all_betas.max()))
    ax.plot([-lim, lim], [-lim, lim], 'r--', linewidth=2, label='Perfect concordance')
    
    # Best fit line
    from scipy.stats import linregress
    slope, intercept, r_value, _, _ = linregress(matched_disc_betas, matched_repl_betas)
    x_fit = np.array([-lim, lim])
    y_fit = slope * x_fit + intercept
    ax.plot(x_fit, y_fit, 'g-', linewidth=2, alpha=0.7,
           label=f'Best fit (r²={r_value**2:.3f})')
    
    # Labels
    ax.set_xlabel("Discovery Effect Size (β)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Replication Effect Size (β)", fontsize=12, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"Replication Concordance (r={correlation:.3f}, n={len(matched_disc_betas):,})",
                    fontsize=14, fontweight="bold")
    
    # Statistics text box
    textstr = f'Correlation: {correlation:.4f}\nSlope: {slope:.4f}\nIntercept: {intercept:.4f}'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
           fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(loc='lower right', fontsize=10)
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"concordance_plot: Saved to {output_path_obj} (r={correlation:.4f})")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_overlapping_variants": len(matched_disc_betas),
        "correlation": float(correlation),
        "slope": float(slope),
        "intercept": float(intercept),
    }



