"""Statistical diagnostic visualizations for GWAS.

This module provides QQ plots, lambda GC calculations, power curves,
and volcano plots to assess statistical quality and identify issues.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

matplotlib.use("Agg", force=True)

from ..core.io import ensure_directory

logger = logging.getLogger(__name__)

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    pd = None


def qq_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    title: str | None = None,
    show_ci: bool = True,
    show_lambda_gc: bool = True,
) -> dict[str, Any]:
    """Q-Q plot with 95% confidence intervals and lambda GC.
    
    Quantile-quantile plot comparing observed vs expected p-values.
    Deviations indicate systematic bias or true associations.
    
    Args:
        results: Association results or path
        output_path: Output path
        title: Plot title
        show_ci: Show 95% confidence intervals
        show_lambda_gc: Show genomic inflation factor
    
    Returns:
        Plot metadata with lambda_GC
    """
    logger.info("qq_plot: Generating QQ plot with diagnostics")
    
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
    
    # Extract p-values
    pvalues = []
    for r in results_list:
        try:
            pval = float(r.get("p_value", 1.0))
            if 0 < pval <= 1:
                pvalues.append(pval)
        except (ValueError, TypeError):
            continue
    
    if len(pvalues) < 2:
        return {"status": "failed", "error": "Insufficient valid p-values"}
    
    pvalues = np.array(sorted(pvalues))
    n = len(pvalues)
    
    # Expected p-values (uniform distribution)
    expected = np.arange(1, n + 1) / (n + 1)
    
    # Convert to -log10
    obs_log = -np.log10(pvalues)
    exp_log = -np.log10(expected)
    
    # Calculate lambda_GC (genomic control inflation factor)
    chi2_obs = stats.chi2.ppf(1 - pvalues, df=1)
    lambda_gc = np.median(chi2_obs) / stats.chi2.ppf(0.5, df=1)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # 95% confidence intervals
    if show_ci:
        # Approximate CI using beta distribution
        alpha = 0.05
        upper_ci = -np.log10(stats.beta.ppf(alpha/2, expected * n, (1 - expected) * n))
        lower_ci = -np.log10(stats.beta.ppf(1 - alpha/2, expected * n, (1 - expected) * n))
        
        ax.fill_between(exp_log, lower_ci, upper_ci, color='gray', alpha=0.25,
                       label='95% CI')
    
    # Identity line (null expectation)
    max_val = max(exp_log.max(), obs_log.max())
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=1.5,
           label='Null expectation', zorder=1)
    
    # Observed points
    ax.scatter(exp_log, obs_log, s=15, alpha=0.6, c='#2E86AB',
              edgecolors='none', rasterized=True, label='Observed', zorder=2)
    
    # Labels
    ax.set_xlabel("Expected -log₁₀(p-value)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Observed -log₁₀(p-value)", fontsize=12, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        title_text = f"Q-Q Plot (n={n:,})"
        if show_lambda_gc:
            title_text += f" | λ_GC = {lambda_gc:.3f}"
        ax.set_title(title_text, fontsize=14, fontweight="bold")
    
    # Lambda GC text box
    if show_lambda_gc:
        textstr = f'λ_GC = {lambda_gc:.4f}'
        if lambda_gc > 1.1:
            textstr += '\n(Possible inflation)'
        elif lambda_gc < 0.9:
            textstr += '\n(Possible deflation)'
        else:
            textstr += '\n(Good calibration)'
        
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
    
    logger.info(f"qq_plot: Saved to {output_path_obj} (lambda_GC={lambda_gc:.4f})")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_variants": n,
        "lambda_gc": float(lambda_gc),
        "inflation_status": "inflated" if lambda_gc > 1.1 else "deflated" if lambda_gc < 0.9 else "good",
    }


def qq_plot_stratified(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    maf_bins: list[tuple[float, float]] | None = None,
    title: str | None = None,
) -> dict[str, Any]:
    """QQ plot stratified by MAF bins.
    
    Checks if genomic inflation varies by allele frequency,
    which can indicate technical artifacts or population stratification.
    
    Args:
        results: Association results or path
        output_path: Output path
        maf_bins: MAF bins [(min, max), ...], default [(0, 0.05), (0.05, 0.5)]
        title: Plot title
    
    Returns:
        Plot metadata with lambda_GC per bin
    """
    logger.info("qq_plot_stratified: Generating stratified QQ plot")
    
    if maf_bins is None:
        maf_bins = [(0.0, 0.05), (0.05, 0.5)]
    
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
    
    # Stratify by MAF
    stratified_data = {f"{low}-{high}": [] for low, high in maf_bins}
    
    for r in results_list:
        try:
            pval = float(r.get("p_value", 1.0))
            maf = float(r.get("MAF", r.get("maf", 0.0)))
            
            if not (0 < pval <= 1):
                continue
            
            for (low, high), key in zip(maf_bins, stratified_data.keys()):
                if low <= maf < high:
                    stratified_data[key].append(pval)
                    break
        except (ValueError, TypeError):
            continue
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    colors = plt.cm.Set2(np.linspace(0, 1, len(maf_bins)))
    lambda_gcs = {}
    
    for (bin_label, pvals), color in zip(stratified_data.items(), colors):
        if len(pvals) < 10:
            continue
        
        pvals = np.array(sorted(pvals))
        n = len(pvals)
        expected = np.arange(1, n + 1) / (n + 1)
        
        obs_log = -np.log10(pvals)
        exp_log = -np.log10(expected)
        
        # Lambda GC
        chi2_obs = stats.chi2.ppf(1 - pvals, df=1)
        lambda_gc = np.median(chi2_obs) / stats.chi2.ppf(0.5, df=1)
        lambda_gcs[bin_label] = lambda_gc
        
        # Plot
        ax.scatter(exp_log, obs_log, s=15, alpha=0.6, c=[color],
                  label=f'MAF {bin_label} (λ={lambda_gc:.3f}, n={n:,})',
                  edgecolors='none', rasterized=True)
    
    # Null expectation line
    max_val = max([ax.get_xlim()[1], ax.get_ylim()[1]])
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=1.5,
           label='Null expectation', zorder=1)
    
    # Labels
    ax.set_xlabel("Expected -log₁₀(p-value)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Observed -log₁₀(p-value)", fontsize=12, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title("Stratified Q-Q Plot by MAF", fontsize=14, fontweight="bold")
    
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(loc='lower right', fontsize=9)
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"qq_plot_stratified: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "lambda_gc_by_maf": lambda_gcs,
    }


def lambda_gc_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    subsets: list[str] | None = None,
    title: str | None = None,
) -> dict[str, Any]:
    """Lambda GC plot across different subsets.
    
    Visualizes genomic inflation across chromosomes or other subsets
    to identify localized issues.
    
    Args:
        results: Association results or path
        output_path: Output path
        subsets: Subset labels (e.g., chromosomes), default is by chromosome
        title: Plot title
    
    Returns:
        Plot metadata with lambda_GC per subset
    """
    logger.info("lambda_gc_plot: Calculating lambda GC by chromosome")
    
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
    
    # Group by chromosome
    chrom_pvals = {}
    for r in results_list:
        try:
            chrom = str(r.get("CHROM", ""))
            pval = float(r.get("p_value", 1.0))
            
            if not (0 < pval <= 1):
                continue
            
            if chrom not in chrom_pvals:
                chrom_pvals[chrom] = []
            chrom_pvals[chrom].append(pval)
        except (ValueError, TypeError):
            continue
    
    # Calculate lambda GC per chromosome
    lambda_gcs = {}
    for chrom, pvals in chrom_pvals.items():
        if len(pvals) < 10:
            continue
        
        pvals_arr = np.array(pvals)
        chi2_obs = stats.chi2.ppf(1 - pvals_arr, df=1)
        lambda_gc = np.median(chi2_obs) / stats.chi2.ppf(0.5, df=1)
        lambda_gcs[chrom] = lambda_gc
    
    if not lambda_gcs:
        return {"status": "failed", "error": "Insufficient data"}
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    chroms = sorted(lambda_gcs.keys())
    lambdas = [lambda_gcs[c] for c in chroms]
    
    x_pos = np.arange(len(chroms))
    colors = ['red' if l > 1.1 else 'orange' if l < 0.9 else 'green' for l in lambdas]
    
    ax.bar(x_pos, lambdas, color=colors, alpha=0.7, edgecolor='black')
    ax.axhline(y=1.0, color='blue', linestyle='--', linewidth=2, label='No inflation')
    ax.axhline(y=1.1, color='red', linestyle=':', linewidth=1, label='Inflation threshold')
    ax.axhline(y=0.9, color='orange', linestyle=':', linewidth=1, label='Deflation threshold')
    
    # Labels
    ax.set_xlabel("Chromosome", fontsize=12, fontweight="bold")
    ax.set_ylabel("λ_GC", fontsize=12, fontweight="bold")
    ax.set_xticks(x_pos)
    ax.set_xticklabels([c.replace('chr', '').replace('NC_', '').split('.')[0] for c in chroms],
                       rotation=45, ha='right')
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        overall_lambda = np.median(lambdas)
        ax.set_title(f"Genomic Inflation by Chromosome (Overall λ_GC = {overall_lambda:.3f})",
                    fontsize=14, fontweight="bold")
    
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3, axis='y', linestyle=':')
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"lambda_gc_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "lambda_gc_by_chrom": lambda_gcs,
        "overall_lambda_gc": float(np.median(lambdas)),
    }


def volcano_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    significance_threshold: float = 5e-8,
    effect_threshold: float = 0.5,
    title: str | None = None,
) -> dict[str, Any]:
    """Volcano plot: effect size vs significance.
    
    Visualizes both statistical significance and effect size magnitude.
    Highlights variants with large effects and strong significance.
    
    Args:
        results: Association results or path
        output_path: Output path
        significance_threshold: P-value threshold
        effect_threshold: Absolute effect size threshold
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info("volcano_plot: Generating volcano plot")
    
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
    betas, pvalues = [], []
    for r in results_list:
        try:
            beta = float(r.get("beta", r.get("BETA", 0.0)))
            pval = float(r.get("p_value", 1.0))
            
            if 0 < pval <= 1:
                betas.append(beta)
                pvalues.append(-math.log10(pval))
        except (ValueError, TypeError):
            continue
    
    if len(betas) < 2:
        return {"status": "failed", "error": "Insufficient data"}
    
    betas = np.array(betas)
    pvalues_log = np.array(pvalues)
    
    # Classify points
    sig_line = -math.log10(significance_threshold)
    
    significant = pvalues_log > sig_line
    large_effect = np.abs(betas) > effect_threshold
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Non-significant
    mask_ns = ~significant
    ax.scatter(betas[mask_ns], pvalues_log[mask_ns], c='gray', s=10, alpha=0.3,
              label='Not significant', rasterized=True)
    
    # Significant but small effect
    mask_sig_small = significant & ~large_effect
    if mask_sig_small.any():
        ax.scatter(betas[mask_sig_small], pvalues_log[mask_sig_small],
                  c='orange', s=15, alpha=0.6,
                  label='Significant (small effect)', rasterized=True)
    
    # Significant and large effect
    mask_sig_large = significant & large_effect
    if mask_sig_large.any():
        ax.scatter(betas[mask_sig_large], pvalues_log[mask_sig_large],
                  c='red', s=20, alpha=0.8,
                  label='Significant (large effect)', rasterized=True)
    
    # Threshold lines
    ax.axhline(y=sig_line, color='blue', linestyle='--', linewidth=1.5,
              label=f'p = {significance_threshold:.0e}')
    ax.axvline(x=effect_threshold, color='green', linestyle=':', linewidth=1)
    ax.axvline(x=-effect_threshold, color='green', linestyle=':', linewidth=1,
              label=f'|β| = {effect_threshold}')
    
    # Labels
    ax.set_xlabel("Effect Size (β)", fontsize=12, fontweight="bold")
    ax.set_ylabel("-log₁₀(p-value)", fontsize=12, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        n_sig_large = mask_sig_large.sum()
        ax.set_title(f"Volcano Plot ({n_sig_large} significant with large effect)",
                    fontsize=14, fontweight="bold")
    
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"volcano_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_variants": len(betas),
        "num_significant": int(significant.sum()),
        "num_large_effect": int(large_effect.sum()),
        "num_significant_large_effect": int(mask_sig_large.sum()),
    }


def power_plot(
    sample_sizes: list[int],
    output_path: str | Path,
    *,
    effect_sizes: list[float] | None = None,
    alpha: float = 5e-8,
    maf: float = 0.3,
    title: str | None = None,
) -> dict[str, Any]:
    """Statistical power curves for GWAS.
    
    Shows expected power to detect associations at different sample sizes
    and effect sizes. Useful for study design and interpretation.
    
    Args:
        sample_sizes: List of sample sizes to evaluate
        output_path: Output path
        effect_sizes: Effect sizes (OR or β) to plot
        alpha: Significance threshold
        maf: Minor allele frequency
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info("power_plot: Generating power curves")
    
    if effect_sizes is None:
        effect_sizes = [0.1, 0.2, 0.3, 0.5, 0.8]
    
    # Calculate power (simplified approximation for continuous trait)
    fig, ax = plt.subplots(figsize=(10, 7))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(effect_sizes)))
    
    for effect, color in zip(effect_sizes, colors):
        powers = []
        for n in sample_sizes:
            # Non-centrality parameter (NCP) approximation
            ncp = n * 2 * maf * (1 - maf) * effect**2
            # Power from non-central chi-square
            threshold = stats.chi2.ppf(1 - alpha, df=1)
            power = 1 - stats.ncx2.cdf(threshold, df=1, nc=ncp)
            powers.append(power)
        
        ax.plot(sample_sizes, powers, linewidth=2, color=color,
               label=f'β = {effect}', marker='o', markersize=4)
    
    # 80% power line
    ax.axhline(y=0.8, color='red', linestyle='--', linewidth=1.5,
              label='80% power', alpha=0.7)
    
    # Labels
    ax.set_xlabel("Sample Size", fontsize=12, fontweight="bold")
    ax.set_ylabel("Statistical Power", fontsize=12, fontweight="bold")
    ax.set_ylim(0, 1.05)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"GWAS Power Curves (α={alpha:.0e}, MAF={maf})",
                    fontsize=14, fontweight="bold")
    
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(loc='lower right', fontsize=10, ncol=2)
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"power_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "effect_sizes": effect_sizes,
        "sample_sizes": sample_sizes,
    }





