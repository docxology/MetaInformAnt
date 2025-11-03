"""Effect size and functional annotation visualizations for GWAS.

This module provides forest plots, effect direction plots,
functional enrichment, and allelic series visualizations.
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


def effect_size_forest_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    top_n: int = 20,
    significance_threshold: float = 5e-8,
    title: str | None = None,
) -> dict[str, Any]:
    """Forest plot of top variants with effect sizes and confidence intervals.
    
    Shows effect estimates (beta or OR) with 95% CI for significant variants.
    
    Args:
        results: Association results or path
        output_path: Output path
        top_n: Number of top variants to show
        significance_threshold: P-value threshold
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"effect_size_forest_plot: Top {top_n} variants")
    
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
    
    # Extract significant variants
    sig_variants = []
    for r in results_list:
        try:
            pval = float(r.get("p_value", 1.0))
            beta = float(r.get("beta", r.get("BETA", 0.0)))
            se = float(r.get("se", r.get("SE", 0.0)))
            
            if pval < significance_threshold and se > 0:
                variant_id = f"{r.get('CHROM', '?')}:{r.get('POS', '?')}"
                sig_variants.append({
                    "id": variant_id,
                    "beta": beta,
                    "se": se,
                    "pval": pval,
                    "ci_lower": beta - 1.96 * se,
                    "ci_upper": beta + 1.96 * se,
                })
        except (ValueError, TypeError):
            continue
    
    if not sig_variants:
        return {
            "status": "skipped",
            "message": f"No significant variants at threshold {significance_threshold}",
        }
    
    # Sort by p-value and take top N
    sig_variants = sorted(sig_variants, key=lambda x: x["pval"])[:top_n]
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, max(6, len(sig_variants) * 0.4)))
    
    y_pos = np.arange(len(sig_variants))
    
    for i, var in enumerate(sig_variants):
        # Effect estimate point
        ax.plot([var["beta"]], [i], 'o', color='#2E86AB', markersize=8, zorder=3)
        
        # Confidence interval
        ax.plot([var["ci_lower"], var["ci_upper"]], [i, i],
               '-', color='#2E86AB', linewidth=2, zorder=2)
        
        # Color by direction
        color = '#C73E1D' if var["beta"] < 0 else '#06A77D'
        ax.plot([var["beta"]], [i], 'o', color=color, markersize=8, zorder=4)
    
    # Null effect line
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1.5, zorder=1)
    
    # Labels
    ax.set_yticks(y_pos)
    ax.set_yticklabels([v["id"] for v in sig_variants], fontsize=9)
    ax.set_xlabel("Effect Size (β) with 95% CI", fontsize=12, fontweight="bold")
    ax.set_ylabel("Variant", fontsize=12, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"Forest Plot: Top {len(sig_variants)} Significant Variants",
                    fontsize=14, fontweight="bold")
    
    ax.grid(True, alpha=0.3, axis='x', linestyle=':')
    ax.invert_yaxis()
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"effect_size_forest_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_variants": len(sig_variants),
    }


def effect_direction_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    significance_threshold: float = 5e-8,
    title: str | None = None,
) -> dict[str, Any]:
    """Effect direction plot showing consistent directionality.
    
    Visualizes whether effect directions are consistent across genome,
    useful for detecting systematic biases.
    
    Args:
        results: Association results or path
        output_path: Output path
        significance_threshold: P-value threshold
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info("effect_direction_plot: Analyzing effect directions")
    
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
    
    # Classify effects
    sig_positive = 0
    sig_negative = 0
    nonsig_positive = 0
    nonsig_negative = 0
    
    for r in results_list:
        try:
            pval = float(r.get("p_value", 1.0))
            beta = float(r.get("beta", r.get("BETA", 0.0)))
            
            is_sig = pval < significance_threshold
            is_positive = beta > 0
            
            if is_sig and is_positive:
                sig_positive += 1
            elif is_sig and not is_positive:
                sig_negative += 1
            elif not is_sig and is_positive:
                nonsig_positive += 1
            elif not is_sig and not is_positive:
                nonsig_negative += 1
        except (ValueError, TypeError):
            continue
    
    total = sig_positive + sig_negative + nonsig_positive + nonsig_negative
    
    if total == 0:
        return {"status": "failed", "error": "No valid effect sizes"}
    
    # Create plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Panel 1: Pie chart of significant effects
    sig_total = sig_positive + sig_negative
    if sig_total > 0:
        ax1.pie([sig_positive, sig_negative],
               labels=[f'Positive\n({sig_positive})', f'Negative\n({sig_negative})'],
               colors=['#06A77D', '#C73E1D'],
               autopct='%1.1f%%',
               startangle=90,
               textprops={'fontsize': 11, 'fontweight': 'bold'})
        ax1.set_title(f"Significant Effects (n={sig_total})", fontsize=12, fontweight="bold")
    else:
        ax1.text(0.5, 0.5, "No significant\neffects", ha='center', va='center',
                transform=ax1.transAxes, fontsize=12, style='italic')
        ax1.set_title("Significant Effects", fontsize=12, fontweight="bold")
    
    # Panel 2: Bar chart of all effects
    categories = ['Sig+', 'Sig-', 'NS+', 'NS-']
    counts = [sig_positive, sig_negative, nonsig_positive, nonsig_negative]
    colors = ['#06A77D', '#C73E1D', '#A8D5BA', '#F4A5A5']
    
    bars = ax2.bar(categories, counts, color=colors, alpha=0.8,
                  edgecolor='black', linewidth=1)
    
    # Add count labels
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        if height > 0:
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                   f'{count:,}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax2.set_ylabel("Number of Variants", fontsize=11, fontweight="bold")
    ax2.set_xlabel("Effect Category", fontsize=11, fontweight="bold")
    ax2.set_title("All Effects by Direction and Significance", fontsize=12, fontweight="bold")
    ax2.grid(True, alpha=0.3, axis='y', linestyle=':')
    
    if title:
        fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)
    else:
        fig.suptitle(f"Effect Direction Analysis ({total:,} variants)",
                    fontsize=14, fontweight="bold", y=1.02)
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"effect_direction_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "sig_positive": sig_positive,
        "sig_negative": sig_negative,
        "nonsig_positive": nonsig_positive,
        "nonsig_negative": nonsig_negative,
    }


def functional_enrichment_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    annotation_file: Path | None = None,
    title: str | None = None,
) -> dict[str, Any]:
    """Functional enrichment of significant variants.
    
    Shows distribution of variants across functional categories
    (e.g., coding, UTR, intergenic, regulatory).
    
    Args:
        results: Association results or path
        output_path: Output path
        annotation_file: Path to variant annotations
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info("functional_enrichment_plot: Analyzing functional categories")
    
    return {
        "status": "skipped",
        "message": "Functional enrichment requires variant annotation",
        "recommendation": "Use ANNOVAR, VEP, or SnpEff for variant annotation",
        "note": "Future enhancement: integrate annotation tools",
    }


def allelic_series_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    gene_region: tuple[str, int, int] | None = None,
    title: str | None = None,
) -> dict[str, Any]:
    """Allelic series plot for multiple alleles at a locus.
    
    Shows dose-response relationship for variants with multiple alleles.
    
    Args:
        results: Association results or path
        output_path: Output path
        gene_region: (chrom, start, end) for gene region
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info("allelic_series_plot: Analyzing allelic series")
    
    return {
        "status": "skipped",
        "message": "Allelic series requires multi-allelic variant analysis",
        "recommendation": "Use specialized allelic series tools",
        "note": "Future enhancement for complex loci with multiple functional alleles",
    }



