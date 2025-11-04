"""Enhanced visualization for GWAS results with genome-scale data support.

This module provides comprehensive visualization methods for GWAS analysis,
integrating techniques from DNA sequence analysis and general visualization modules.

Supports genome-scale data (millions of SNPs) with efficient rendering and
multiple visualization types for publication-quality figures.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

matplotlib.use("Agg", force=True)

from ..core.io import ensure_directory

logger = logging.getLogger(__name__)

# Try importing pandas and seaborn for advanced plotting
try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    pd = None

try:
    import seaborn as sns
    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False
    sns = None


def manhattan_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    significance_threshold: float | None = None,
    suggestiveness_threshold: float | None = None,
    chrom_colors: list[str] | None = None,
    title: str | None = None,
    max_points: int | None = None,
) -> dict[str, Any]:
    """Generate Manhattan plot for GWAS results (genome-scale optimized).
    
    Handles millions of SNPs efficiently by implementing:
    - Point thinning for non-significant regions
    - Memory-efficient data structures
    - Optimized rendering
    
    Args:
        results: List of result dictionaries or path to results file (TSV)
        output_path: Path to save plot
        significance_threshold: P-value threshold for genome-wide significance (default: 5e-8)
        suggestiveness_threshold: P-value threshold for suggestiveness (default: 1e-5)
        chrom_colors: List of colors to alternate for chromosomes
        title: Plot title
        max_points: Maximum points to plot per chromosome (for thinning)
    
    Returns:
        Dictionary with plot metadata
    """
    logger.info(f"manhattan_plot: Generating Manhattan plot (genome-scale optimized)")
    
    if significance_threshold is None:
        significance_threshold = 5e-8
    if suggestiveness_threshold is None:
        suggestiveness_threshold = 1e-5
    
    # Load results
    if isinstance(results, (Path, str)):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(results, sep="\t")
            results_list = df.to_dict("records")
        else:
            from ..core.io import read_tsv
            data_list = read_tsv(results)
            if len(data_list) < 2:
                return {"status": "failed", "error": "Invalid results file"}
            header = data_list[0]
            results_list = [{header[i]: row[i] for i in range(len(header))} for row in data_list[1:]]
    else:
        results_list = results
    
    if not results_list:
        return {"status": "failed", "error": "No results to plot"}
    
    # Extract data
    chroms, positions, pvalues = [], [], []
    for result in results_list:
        try:
            chrom = str(result.get("CHROM", ""))
            pos = int(result.get("POS", 0))
            pval = float(result.get("p_value", 1.0))
            if pval > 0:
                chroms.append(chrom)
                positions.append(pos)
                pvalues.append(pval)
        except (ValueError, TypeError):
            continue
    
    if not pvalues:
        return {"status": "failed", "error": "No valid p-values found"}
    
    logger.info(f"manhattan_plot: Plotting {len(pvalues)} variants")
    
    # Group by chromosome
    unique_chroms = sorted(set(chroms), key=lambda x: (int(x.replace('chr', '').replace('NC_', '').split('.')[0]) if any(c.isdigit() for c in x) else 999, x))
    chrom_data = {chrom: {"pos": [], "pval": []} for chrom in unique_chroms}
    
    for chrom, pos, pval in zip(chroms, positions, pvalues):
        if chrom in chrom_data:
            chrom_data[chrom]["pos"].append(pos)
            chrom_data[chrom]["pval"].append(pval)
    
    # Calculate cumulative positions
    cumulative_pos = 0
    chrom_cumulative = {}
    chrom_spans = {}
    
    for chrom in unique_chroms:
        chrom_cumulative[chrom] = cumulative_pos
        if chrom_data[chrom]["pos"]:
            min_pos = min(chrom_data[chrom]["pos"])
            max_pos = max(chrom_data[chrom]["pos"])
            span = max_pos - min_pos
            chrom_spans[chrom] = (cumulative_pos, cumulative_pos + span)
            cumulative_pos += span + 10000000  # Gap between chromosomes
    
    # Create plot
    fig, ax = plt.subplots(figsize=(16, 7))
    
    # Colors
    if chrom_colors is None:
        chrom_colors = ["#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#06A77D", "#6A4C93"]
    
    # Plot points
    for chrom_idx, chrom in enumerate(unique_chroms):
        if not chrom_data[chrom]["pos"]:
            continue
        
        positions_chrom = np.array(chrom_data[chrom]["pos"])
        pvalues_chrom = np.array(chrom_data[chrom]["pval"])
        
        # Calculate plot positions
        x_pos = chrom_cumulative[chrom] + positions_chrom - positions_chrom.min()
        y_neg_log_p = -np.log10(pvalues_chrom)
        
        # Point thinning for non-significant variants (genome-scale optimization)
        if max_points and len(x_pos) > max_points:
            # Keep all significant/suggestive points
            significant_mask = pvalues_chrom < suggestiveness_threshold
            significant_indices = np.where(significant_mask)[0]
            non_significant_indices = np.where(~significant_mask)[0]
            
            # Sample non-significant points
            if len(non_significant_indices) > (max_points - len(significant_indices)):
                sampled_non_sig = np.random.choice(
                    non_significant_indices,
                    size=max_points - len(significant_indices),
                    replace=False
                )
                keep_indices = np.concatenate([significant_indices, sampled_non_sig])
                keep_indices = np.sort(keep_indices)
            else:
                keep_indices = np.arange(len(x_pos))
            
            x_pos = x_pos[keep_indices]
            y_neg_log_p = y_neg_log_p[keep_indices]
        
        color = chrom_colors[chrom_idx % len(chrom_colors)]
        ax.scatter(x_pos, y_neg_log_p, c=color, s=15, alpha=0.6, 
                   label=chrom if chrom_idx < 10 else None, rasterized=True)
    
    # Significance lines
    sig_line = -math.log10(significance_threshold)
    sugg_line = -math.log10(suggestiveness_threshold)
    
    ax.axhline(y=sig_line, color="r", linestyle="--", linewidth=1.5,
               label=f"Genome-wide ({significance_threshold:.0e})")
    ax.axhline(y=sugg_line, color="orange", linestyle="--", linewidth=1,
               label=f"Suggestive ({suggestiveness_threshold:.0e})")
    
    # Labels
    ax.set_xlabel("Chromosome", fontsize=13, fontweight="bold")
    ax.set_ylabel("-log₁₀(p-value)", fontsize=13, fontweight="bold")
    if title:
        ax.set_title(title, fontsize=15, fontweight="bold")
    else:
        ax.set_title(f"Manhattan Plot ({len(pvalues):,} variants)", fontsize=15, fontweight="bold")
    
    # X-axis ticks at chromosome centers
    chrom_centers = [(chrom_spans[chrom][0] + chrom_spans[chrom][1]) / 2 
                     for chrom in unique_chroms if chrom in chrom_spans]
    ax.set_xticks(chrom_centers)
    ax.set_xticklabels([c.replace('chr', '').replace('NC_', '').split('.')[0] for c in unique_chroms], 
                        rotation=45, ha="right", fontsize=9)
    
    ax.grid(True, alpha=0.2, axis="y", linestyle=":")
    ax.legend(loc="upper right", fontsize=8, ncol=2, framealpha=0.9)
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
        "num_significant": sum(1 for p in pvalues if p < significance_threshold),
        "num_suggestive": sum(1 for p in pvalues if p < suggestiveness_threshold),
    }


def qq_plot(
    pvalues: list[float] | Path,
    output_path: str | Path,
    *,
    title: str | None = None,
    show_ci: bool = True,
) -> dict[str, Any]:
    """Generate Q-Q plot for p-value distribution (genome-scale optimized).
    
    Args:
        pvalues: List of p-values or path to results file
        output_path: Path to save plot
        title: Plot title
        show_ci: Show 95% confidence interval
    
    Returns:
        Dictionary with plot metadata
    """
    logger.info(f"qq_plot: Generating Q-Q plot")
    
    # Load p-values
    if isinstance(pvalues, (Path, str)):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(pvalues, sep="\t")
            pvalues_list = df["p_value"].tolist()
        else:
            from ..core.io import read_tsv
            data_list = read_tsv(pvalues)
            header = data_list[0]
            pval_idx = header.index("p_value")
            pvalues_list = [float(row[pval_idx]) for row in data_list[1:] 
                           if len(row) > pval_idx and 0 < float(row[pval_idx]) <= 1]
    else:
        pvalues_list = [float(p) for p in pvalues if 0 < p <= 1]
    
    if not pvalues_list:
        return {"status": "failed", "error": "No valid p-values"}
    
    pvalues_arr = np.array(pvalues_list)
    pvalues_arr = pvalues_arr[(pvalues_arr > 0) & (pvalues_arr <= 1)]
    
    # Sort
    observed_pvals = np.sort(pvalues_arr)
    n = len(observed_pvals)
    expected_pvals = np.linspace(1.0 / n, 1.0, n)
    
    # Convert to -log10
    observed_neg_log = -np.log10(observed_pvals)
    expected_neg_log = -np.log10(expected_pvals)
    
    # Calculate lambda_GC
    chi2_obs = -2 * np.log(observed_pvals)
    median_chi2 = np.median(chi2_obs)
    lambda_gc = median_chi2 / 0.4549364  # Median of chi2(1)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # 95% CI
    if show_ci:
        ci_upper = -np.log10(expected_pvals * (1 + 1.96 / np.sqrt(n)))
        ci_lower = -np.log10(expected_pvals * (1 - 1.96 / np.sqrt(n)))
        ax.fill_between(expected_neg_log, ci_lower, ci_upper,
                         color='gray', alpha=0.2, label='95% CI')
    
    # Plot points
    ax.scatter(expected_neg_log, observed_neg_log, s=20, alpha=0.6, 
               c="#2E86AB", rasterized=True)
    
    # Diagonal line
    max_val = max(np.max(expected_neg_log), np.max(observed_neg_log))
    ax.plot([0, max_val], [0, max_val], "r--", linewidth=2, label="Expected (null)")
    
    # Labels
    ax.set_xlabel("Expected -log₁₀(p-value)", fontsize=13, fontweight="bold")
    ax.set_ylabel("Observed -log₁₀(p-value)", fontsize=13, fontweight="bold")
    if title:
        ax.set_title(f"{title}\nλ_GC = {lambda_gc:.4f} ({n:,} variants)", 
                     fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"Q-Q Plot\nλ_GC = {lambda_gc:.4f} ({n:,} variants)", 
                     fontsize=14, fontweight="bold")
    
    ax.grid(True, alpha=0.2, linestyle=":")
    ax.legend(loc="upper left", fontsize=10, framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"qq_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "lambda_gc": float(lambda_gc),
        "num_pvalues": len(pvalues_arr),
    }


def regional_plot(
    results: list[dict[str, Any]] | Path,
    region: str,
    output_path: str | Path,
    *,
    window: int = 500000,
    title: str | None = None,
    genes: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Generate regional association plot with optional gene annotations.
    
    Args:
        results: Results or path to file
        region: Format "chr:start-end" or "chr:center"
        output_path: Save path
        window: Window size if single position (bp)
        title: Plot title
        genes: Optional gene annotations [{"name": str, "start": int, "end": int}]
    
    Returns:
        Plot metadata
    """
    logger.info(f"regional_plot: {region}")
    
    # Parse region
    if ":" in region:
        parts = region.split(":")
        chrom = parts[0]
        if "-" in parts[1]:
            start, end = map(int, parts[1].split("-"))
        else:
            center = int(parts[1])
            start, end = center - window, center + window
    else:
        return {"status": "failed", "error": "Invalid region format"}
    
    # Load results
    if isinstance(results, (Path, str)):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(results, sep="\t")
            results_list = df.to_dict("records")
        else:
            from ..core.io import read_tsv
            data_list = read_tsv(results)
            header = data_list[0]
            results_list = [{header[i]: row[i] for i in range(len(header))} 
                           for row in data_list[1:]]
    else:
        results_list = results
    
    # Filter by region
    regional_results = []
    for result in results_list:
        try:
            result_chrom = str(result.get("CHROM", ""))
            result_pos = int(result.get("POS", 0))
            if result_chrom == chrom and start <= result_pos <= end:
                regional_results.append(result)
        except (ValueError, TypeError):
            continue
    
    if not regional_results:
        return {"status": "failed", "error": f"No variants in {region}"}
    
    # Extract data
    positions = []
    pvalues = []
    for result in regional_results:
        try:
            pos = int(result.get("POS", 0))
            pval = float(result.get("p_value", 1.0))
            if pval > 0:
                positions.append(pos)
                pvalues.append(-math.log10(pval))
        except (ValueError, TypeError):
            continue
    
    # Create plot
    if genes:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), 
                                       gridspec_kw={'height_ratios': [3, 1]},
                                       sharex=True)
        main_ax = ax1
    else:
        fig, main_ax = plt.subplots(figsize=(12, 6))
    
    # Main association plot
    main_ax.scatter(positions, pvalues, s=30, alpha=0.7, c="#2E86AB")
    
    # Significance line
    sig_line = -math.log10(5e-8)
    main_ax.axhline(y=sig_line, color="r", linestyle="--", linewidth=1,
                    label="Genome-wide significance")
    
    main_ax.set_ylabel("-log₁₀(p-value)", fontsize=12, fontweight="bold")
    if title:
        main_ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        main_ax.set_title(f"Regional Plot: {chrom}:{start:,}-{end:,}", 
                         fontsize=14, fontweight="bold")
    
    main_ax.grid(True, alpha=0.2, axis="y", linestyle=":")
    main_ax.legend(loc="upper right", fontsize=10)
    main_ax.spines['top'].set_visible(False)
    main_ax.spines['right'].set_visible(False)
    
    # Gene annotations
    if genes:
        for gene in genes:
            gene_start = gene.get("start", 0)
            gene_end = gene.get("end", 0)
            gene_name = gene.get("name", "Unknown")
            
            if start <= gene_start <= end or start <= gene_end <= end:
                # Draw gene as rectangle
                width = gene_end - gene_start
                rect = Rectangle((gene_start, 0), width, 1, 
                                linewidth=1, edgecolor='black',
                                facecolor='lightblue', alpha=0.7)
                ax2.add_patch(rect)
                
                # Gene name
                mid = (gene_start + gene_end) / 2
                ax2.text(mid, 0.5, gene_name, ha='center', va='center',
                        fontsize=8, fontweight='bold')
        
        ax2.set_xlim(start, end)
        ax2.set_ylim(0, 1)
        ax2.set_xlabel(f"Position on {chrom}", fontsize=12, fontweight="bold")
        ax2.set_ylabel("Genes", fontsize=10)
        ax2.set_yticks([])
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
    else:
        main_ax.set_xlabel(f"Position on {chrom}", fontsize=12, fontweight="bold")
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"regional_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "region": region,
        "num_variants": len(positions),
    }


def effect_size_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    top_n: int = 20,
    title: str | None = None,
) -> dict[str, Any]:
    """Forest plot showing effect sizes for top variants.
    
    Args:
        results: Results or path
        output_path: Save path
        top_n: Number of top variants to show
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"effect_size_plot: Top {top_n} variants")
    
    # Load
    if isinstance(results, (Path, str)):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(results, sep="\t")
            results_list = df.to_dict("records")
        else:
            from ..core.io import read_tsv
            data_list = read_tsv(results)
            header = data_list[0]
            results_list = [{header[i]: row[i] for i in range(len(header))} 
                           for row in data_list[1:]]
    else:
        results_list = results
    
    # Sort by p-value
    sorted_results = sorted(results_list, 
                           key=lambda x: float(x.get("p_value", 1.0)))[:top_n]
    
    # Extract data
    variant_ids = []
    betas = []
    ses = []
    
    for result in sorted_results:
        try:
            chrom = result.get("CHROM", "")
            pos = result.get("POS", "")
            variant_id = f"{chrom}:{pos}"
            beta = float(result.get("beta", 0.0))
            se = float(result.get("se", 0.0))
            
            variant_ids.append(variant_id)
            betas.append(beta)
            ses.append(se)
        except (ValueError, TypeError):
            continue
    
    if not betas:
        return {"status": "failed", "error": "No valid effect sizes"}
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, max(6, len(betas) * 0.4)))
    
    y_pos = np.arange(len(betas))
    
    # Plot effect sizes with error bars
    ax.errorbar(betas, y_pos, xerr=[1.96 * se for se in ses],
                fmt='o', markersize=8, capsize=5, capthick=2,
                color='#2E86AB', ecolor='gray', alpha=0.8)
    
    # Null line
    ax.axvline(x=0, color='r', linestyle='--', linewidth=1, alpha=0.7)
    
    # Labels
    ax.set_yticks(y_pos)
    ax.set_yticklabels(variant_ids, fontsize=9)
    ax.set_xlabel("Effect Size (β) ± 95% CI", fontsize=12, fontweight="bold")
    ax.set_ylabel("Variant", fontsize=12, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"Effect Sizes (Top {len(betas)} Variants)", 
                    fontsize=14, fontweight="bold")
    
    ax.grid(True, alpha=0.2, axis="x", linestyle=":")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"effect_size_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_variants": len(betas),
    }


def pca_plot(
    pca_components: Path | np.ndarray,
    output_path: str | Path,
    *,
    phenotypes: dict[str, Any] | None = None,
    pc1: int = 1,
    pc2: int = 2,
    title: str | None = None,
) -> dict[str, Any]:
    """Plot PCA results colored by phenotype.
    
    Args:
        pca_components: Path to PCA TSV or numpy array
        output_path: Save path
        phenotypes: Dict mapping sample_id to phenotype
        pc1: PC number for x-axis (1-indexed)
        pc2: PC number for y-axis (1-indexed)
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"pca_plot: PC{pc1} vs PC{pc2}")
    
    # Load PCA
    if isinstance(pca_components, (Path, str)):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(pca_components, sep="\t")
        else:
            return {"status": "failed", "error": "pandas required for PCA plot"}
    else:
        return {"status": "failed", "error": "PCA path required"}
    
    # Extract PCs
    pc1_col = f"PC{pc1}"
    pc2_col = f"PC{pc2}"
    
    if pc1_col not in df.columns or pc2_col not in df.columns:
        return {"status": "failed", "error": f"PCs {pc1}/{pc2} not found"}
    
    x = df[pc1_col].values
    y = df[pc2_col].values
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if phenotypes and "sample_id" in df.columns:
        # Color by phenotype
        colors = []
        for sample_id in df["sample_id"]:
            pheno = phenotypes.get(sample_id, "unknown")
            colors.append(pheno)
        
        if SEABORN_AVAILABLE:
            sns.scatterplot(x=x, y=y, hue=colors, s=100, alpha=0.7, ax=ax)
        else:
            # Simple coloring
            unique_phenos = list(set(colors))
            color_map = {p: i for i, p in enumerate(unique_phenos)}
            numeric_colors = [color_map[c] for c in colors]
            scatter = ax.scatter(x, y, c=numeric_colors, s=100, alpha=0.7,
                               cmap='tab10')
            ax.legend(handles=[plt.Line2D([0], [0], marker='o', color='w',
                      markerfacecolor=scatter.cmap(scatter.norm(color_map[p])),
                      markersize=10, label=p) for p in unique_phenos])
    else:
        ax.scatter(x, y, s=100, alpha=0.7, c='#2E86AB')
    
    # Labels
    ax.set_xlabel(f"PC{pc1}", fontsize=13, fontweight="bold")
    ax.set_ylabel(f"PC{pc2}", fontsize=13, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"PCA: PC{pc1} vs PC{pc2}", fontsize=14, fontweight="bold")
    
    ax.grid(True, alpha=0.2, linestyle=":")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"pca_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_samples": len(x),
    }


def kinship_heatmap(
    kinship_matrix: Path | np.ndarray,
    output_path: str | Path,
    *,
    title: str | None = None,
    sample_labels: list[str] | None = None,
) -> dict[str, Any]:
    """Heatmap visualization of kinship matrix.
    
    Args:
        kinship_matrix: Path to TSV or numpy array
        output_path: Save path
        title: Plot title
        sample_labels: Optional sample labels
    
    Returns:
        Plot metadata
    """
    logger.info(f"kinship_heatmap: Visualizing relatedness matrix")
    
    # Load matrix
    if isinstance(kinship_matrix, (Path, str)):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(kinship_matrix, sep="\t", index_col=0)
            matrix = df.values
            if sample_labels is None:
                sample_labels = df.index.tolist()
        else:
            matrix = np.loadtxt(kinship_matrix, skiprows=1, delimiter="\t")[:, 1:]
    else:
        matrix = kinship_matrix
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    if SEABORN_AVAILABLE:
        sns.heatmap(matrix, cmap="RdBu_r", center=0, 
                   xticklabels=sample_labels if sample_labels else False,
                   yticklabels=sample_labels if sample_labels else False,
                   cbar_kws={'label': 'Kinship Coefficient'},
                   ax=ax)
    else:
        im = ax.imshow(matrix, cmap="RdBu_r", aspect="auto")
        plt.colorbar(im, ax=ax, label="Kinship Coefficient")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"Kinship Matrix ({matrix.shape[0]} samples)", 
                    fontsize=14, fontweight="bold")
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"kinship_heatmap: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "matrix_size": matrix.shape[0],
    }


def maf_distribution(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    bins: int = 50,
    title: str | None = None,
) -> dict[str, Any]:
    """Plot minor allele frequency distribution.
    
    Args:
        results: Results or path
        output_path: Save path
        bins: Number of histogram bins
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"maf_distribution: Plotting MAF histogram")
    
    # Load
    if isinstance(results, (Path, str)):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(results, sep="\t")
            results_list = df.to_dict("records")
        else:
            from ..core.io import read_tsv
            data_list = read_tsv(results)
            header = data_list[0]
            results_list = [{header[i]: row[i] for i in range(len(header))} 
                           for row in data_list[1:]]
    else:
        results_list = results
    
    # Extract MAF if available, otherwise estimate from p-value distribution
    # (This is a placeholder - real MAF would come from VCF data)
    mafs = []
    for result in results_list:
        if "maf" in result:
            try:
                maf = float(result["maf"])
                if 0 < maf <= 0.5:
                    mafs.append(maf)
            except (ValueError, TypeError):
                continue
    
    if not mafs:
        return {"status": "skipped", "message": "No MAF data available"}
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.hist(mafs, bins=bins, color='#2E86AB', alpha=0.7, edgecolor='black')
    
    ax.set_xlabel("Minor Allele Frequency", fontsize=12, fontweight="bold")
    ax.set_ylabel("Number of Variants", fontsize=12, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"MAF Distribution ({len(mafs):,} variants)", 
                    fontsize=14, fontweight="bold")
    
    ax.grid(True, alpha=0.2, axis="y", linestyle=":")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
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
        "mean_maf": float(np.mean(mafs)),
        "median_maf": float(np.median(mafs)),
    }





