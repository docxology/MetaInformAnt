"""Cross-species visualization module.

Centralizes complex matplotlib and seaborn plotting logic for cross-species 
expression and ortholog analyses. Methods are designed to be domain-agnostic 
but accept taxonomic groupings for customized styling.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import seaborn as sns
from typing import Dict, Optional

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def _fmt(name: str) -> str:
    """Format species names for display."""
    return str(name).replace("_", " ")


def _get_family_color(sp: str, family_map: Optional[Dict[str, str]] = None, family_colors: Optional[Dict[str, str]] = None) -> str:
    """Get the color for a species based on its family mapping."""
    if not family_map or not family_colors:
        return '#34495e'  # Default dark gray-blue
    return family_colors.get(family_map.get(sp, ''), '#95a5a6')


def plot_divergence_heatmap(
    div_matrix: pd.DataFrame, 
    output_path: Path, 
    family_map: Optional[Dict[str, str]] = None,
    family_colors: Optional[Dict[str, str]] = None,
    title: str = "Expression Divergence Matrix",
    cbar_label: str = "Divergence",
) -> None:
    """Generate a clustered heatmap of expression divergence.
    
    Args:
        div_matrix: Symmetric divergence matrix (1 - correlation)
        output_path: Path to save the plot
        family_map: Dictionary mapping species names to taxonomic families
        family_colors: Dictionary mapping families to hex colors
        title: Plot title
        cbar_label: Colorbar label
    """
    n = len(div_matrix)
    condensed = squareform(div_matrix.values, checks=False)
    condensed = np.nan_to_num(condensed, nan=0.5)
    link = linkage(condensed, method='ward')
    order = [div_matrix.index[i] for i in leaves_list(link)]
    m = div_matrix.loc[order, order]

    fig, ax = plt.subplots(figsize=(max(10, n*0.5), max(9, n*0.45)))
    
    sns.heatmap(
        m, cmap="YlOrRd", annot=True if n <= 30 else False, fmt=".2f", ax=ax,
        xticklabels=[_fmt(s) for s in order],
        yticklabels=[_fmt(s) for s in order],
        linewidths=0.4 if n <= 30 else 0.1, vmin=0.25, vmax=0.75,
        annot_kws={"fontsize": 6.5, "fontweight": "bold"},
        cbar_kws={"label": cbar_label, "shrink": 0.8},
    )
    
    # Color diagonal
    for i in range(n):
        ax.add_patch(plt.Rectangle((i, i), 1, 1, fill=True, color='#2c3e50', zorder=3))

    # Color tick labels by family
    if family_map and family_colors:
        for label in ax.get_yticklabels():
            sp = label.get_text().replace(" ", "_")
            label.set_color(_get_family_color(sp, family_map, family_colors))
            label.set_fontsize(8)
        for label in ax.get_xticklabels():
            sp = label.get_text().replace(" ", "_")
            label.set_color(_get_family_color(sp, family_map, family_colors))
            label.set_fontsize(8)

    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_path, dpi=250, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved divergence heatmap to {output_path}")


def plot_dendrogram(
    div_matrix: pd.DataFrame, 
    output_path: Path,
    family_map: Optional[Dict[str, str]] = None,
    family_colors: Optional[Dict[str, str]] = None,
    title: str = "Hierarchical Clustering of Expression Divergence",
) -> None:
    """Generate a hierarchical clustering dendrogram annotated with taxonomy.
    
    Args:
        div_matrix: Symmetric divergence matrix
        output_path: Path to save the plot
        family_map: Dictionary mapping species names to taxonomic families
        family_colors: Dictionary mapping families to hex colors
        title: Plot title
    """
    condensed = squareform(div_matrix.values, checks=False)
    condensed = np.nan_to_num(condensed, nan=0.5)
    link = linkage(condensed, method='ward')

    fig, ax = plt.subplots(figsize=(max(12, len(div_matrix)*0.5), 8))
    labels = [_fmt(s) for s in div_matrix.index]
    dn = dendrogram(link, labels=labels, leaf_rotation=40, leaf_font_size=9,
                     ax=ax, color_threshold=0.4, above_threshold_color='#7f8c8d')

    if family_map and family_colors:
        for label in ax.get_xticklabels():
            sp = label.get_text().replace(" ", "_")
            label.set_color(_get_family_color(sp, family_map, family_colors))
            label.set_fontweight('bold')
        
        legend = [Patch(facecolor=c, label=f) for f, c in family_colors.items()]
        ax.legend(handles=legend, loc='upper right', fontsize=9, title="Family")

    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.set_ylabel("Ward Distance", fontsize=11)
    plt.tight_layout()
    plt.savefig(output_path, dpi=250, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved dendrogram to {output_path}")


def plot_coverage(
    coverage_series: pd.Series, 
    total_groups: int, 
    output_path: Path,
    family_map: Optional[Dict[str, str]] = None,
    family_colors: Optional[Dict[str, str]] = None,
) -> None:
    """Generate a coverage barplot with optional family colors.
    
    Args:
        coverage_series: Series mapping species names to orthogroup coverage counts
        total_groups: Total number of orthogroups
        output_path: Path to save the plot
        family_map: Dictionary mapping species names to taxonomic families
        family_colors: Dictionary mapping families to hex colors
    """
    cov = coverage_series.sort_values(ascending=True)
    fig, ax = plt.subplots(figsize=(11, max(6, len(cov)*0.3)))
    
    colors = [_get_family_color(sp, family_map, family_colors) for sp in cov.index]
    ax.barh(range(len(cov)), cov.values, color=colors, edgecolor='white', linewidth=0.5)
    ax.set_yticks(range(len(cov)))
    ax.set_yticklabels([_fmt(s) for s in cov.index], fontsize=9)

    for i, (sp, val) in enumerate(cov.items()):
        pct = 100 * val / total_groups if total_groups > 0 else 0
        label = f"{val:,} ({pct:.0f}%)" if val > 0 else "No data"
        ax.text(val + (cov.max() * 0.02), i, label, va='center', fontsize=8,
                fontweight='bold' if val > 0 else 'normal',
                color='#2c3e50' if val > 0 else '#e74c3c')

    ax.set_xlabel("Number of Orthogroups with Mapped Transcripts", fontsize=11)
    ax.set_title(f"Ortholog Coverage per Species\nTotal orthogroups: {total_groups:,}", fontsize=13, fontweight='bold')
    ax.set_xlim(0, cov.max() * 1.25)

    if family_map and family_colors:
        legend = [Patch(facecolor=c, label=f) for f, c in family_colors.items()]
        ax.legend(handles=legend, loc='lower right', fontsize=8, title="Family")
        
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved coverage barplot to {output_path}")


def plot_top_pairs(div_matrix: pd.DataFrame, output_path: Path) -> None:
    """Horizontal bar chart of top-10 most similar and most divergent pairs.
    
    Args:
        div_matrix: Symmetric divergence matrix
        output_path: Path to save the plot
    """
    pairs = []
    n = len(div_matrix)
    for i in range(n):
        for j in range(i+1, n):
            pairs.append((div_matrix.index[i], div_matrix.columns[j], div_matrix.iloc[i, j]))
    pairs.sort(key=lambda x: x[2])

    top_sim = pairs[:10]
    top_div = pairs[-10:][::-1]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))

    # Most similar
    labels_s = [f"{_fmt(a)}\n↔ {_fmt(b)}" for a, b, _ in top_sim]
    vals_s = [d for _, _, d in top_sim]
    colors_s = plt.cm.YlGn(np.linspace(0.3, 0.8, 10))
    ax1.barh(range(len(vals_s)), vals_s, color=colors_s[::-1], edgecolor='white')
    ax1.set_yticks(range(len(vals_s)))
    ax1.set_yticklabels(labels_s, fontsize=8)
    ax1.set_xlabel("Divergence", fontsize=10)
    ax1.set_title("Top 10 Most Similar Pairs", fontsize=12, fontweight='bold')
    ax1.invert_yaxis()
    for i, v in enumerate(vals_s):
        ax1.text(v + 0.003, i, f"{v:.3f}", va='center', fontsize=8, fontweight='bold')
    ax1.set_xlim(0, max(vals_s) * 1.15 if vals_s else 1)

    # Most divergent
    labels_d = [f"{_fmt(a)}\n↔ {_fmt(b)}" for a, b, _ in top_div]
    vals_d = [d for _, _, d in top_div]
    colors_d = plt.cm.OrRd(np.linspace(0.3, 0.8, 10))
    ax2.barh(range(len(vals_d)), vals_d, color=colors_d, edgecolor='white')
    ax2.set_yticks(range(len(vals_d)))
    ax2.set_yticklabels(labels_d, fontsize=8)
    ax2.set_xlabel("Divergence", fontsize=10)
    ax2.set_title("Top 10 Most Divergent Pairs", fontsize=12, fontweight='bold')
    ax2.invert_yaxis()
    for i, v in enumerate(vals_d):
        ax2.text(v + 0.003, i, f"{v:.3f}", va='center', fontsize=8, fontweight='bold')
    ax2.set_xlim(0, max(vals_d) * 1.08 if vals_d else 1)

    fig.suptitle("Pairwise Expression Divergence Extremes", fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved top pairs plot to {output_path}")


def plot_family_violin(
    div_matrix: pd.DataFrame, 
    output_path: Path,
    family_map: Dict[str, str],
) -> None:
    """Violin plot of divergence grouped by within-family vs between-family.
    
    Args:
        div_matrix: Symmetric divergence matrix
        output_path: Path to save the plot
        family_map: Dictionary mapping species names to taxonomic families
    """
    records = []
    n = len(div_matrix)
    for i in range(n):
        for j in range(i+1, n):
            sp_a, sp_b = div_matrix.index[i], div_matrix.columns[j]
            fam_a = family_map.get(sp_a, 'Unknown')
            fam_b = family_map.get(sp_b, 'Unknown')
            d = div_matrix.iloc[i, j]
            if fam_a == fam_b:
                records.append({"Category": f"Within {fam_a}", "Divergence": d, "Type": "Within-family"})
            else:
                records.append({"Category": f"{min(fam_a,fam_b)}–{max(fam_a,fam_b)}", "Divergence": d, "Type": "Between-family"})

    df = pd.DataFrame(records)
    if df.empty:
        logger.warning("No data for family violin plot.")
        return
        
    # Keep only categories with >=3 data points
    counts = df["Category"].value_counts()
    keep = counts[counts >= 3].index
    df = df[df["Category"].isin(keep)]

    order = df.groupby("Category")["Divergence"].median().sort_values().index.tolist()

    fig, ax = plt.subplots(figsize=(14, 7))
    palette = {cat: '#27ae60' if 'Within' in cat else '#e67e22' for cat in order}
    sns.violinplot(data=df, x="Category", y="Divergence", order=order, palette=palette,
                   inner="box", linewidth=1, ax=ax, cut=0)
    sns.stripplot(data=df, x="Category", y="Divergence", order=order, color='#2c3e50',
                  size=3, alpha=0.4, jitter=True, ax=ax)

    ax.set_title("Expression Divergence by Taxonomic Relationship", fontsize=14, fontweight='bold')
    ax.set_ylabel("Divergence", fontsize=11)
    ax.set_xlabel("")
    plt.xticks(rotation=35, ha='right', fontsize=9)
    legend = [Patch(facecolor='#27ae60', label='Within-family'),
              Patch(facecolor='#e67e22', label='Between-family')]
    ax.legend(handles=legend, fontsize=9)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved family violin plot to {output_path}")


def plot_method_comparison(
    ortholog_matrix: pd.DataFrame,
    fingerprint_matrix: pd.DataFrame,
    output_path: Path,
    family_map: Optional[Dict[str, str]] = None,
) -> None:
    """Scatter comparing ortholog-mapped vs fingerprint divergence for shared pairs.
    
    Args:
        ortholog_matrix: Divergence matrix based on ortholog mapping
        fingerprint_matrix: Divergence matrix based on fingerprint distribution
        output_path: Path to save the plot
        family_map: Optional dictionary mapping species to families for coloring
    """
    shared = sorted(set(ortholog_matrix.index) & set(fingerprint_matrix.index))
    if len(shared) < 2:
        logger.warning("Not enough shared species for method comparison.")
        return
        
    orth = ortholog_matrix.loc[shared, shared]
    fing = fingerprint_matrix.loc[shared, shared]

    n = len(shared)
    records = []
    for i in range(n):
        for j in range(i+1, n):
            o = orth.iloc[i, j]
            f = fing.iloc[i, j]
            if o < 1.0:  # Valid pair
                sp_a, sp_b = shared[i], shared[j]
                cat = "All Pairs"
                if family_map:
                    fam_a = family_map.get(sp_a, '?')
                    fam_b = family_map.get(sp_b, '?')
                    cat = "Within-family" if fam_a == fam_b else "Between-family"
                records.append({"Ortholog": o, "Fingerprint": f, "Category": cat})

    df = pd.DataFrame(records)
    if df.empty:
        logger.warning("No valid pairs for method comparison.")
        return

    fig, ax = plt.subplots(figsize=(9, 9))
    colors = {'Within-family': '#27ae60', 'Between-family': '#e67e22', 'All Pairs': '#3498db'}
    for cat, grp in df.groupby("Category"):
        ax.scatter(grp["Fingerprint"], grp["Ortholog"], c=colors.get(cat, '#3498db'), 
                   s=30, alpha=0.6, edgecolors='white', linewidth=0.5, label=cat, zorder=3)

    # Reference line
    lims = [min(df["Fingerprint"].min(), df["Ortholog"].min()) - 0.05,
            max(df["Fingerprint"].max(), df["Ortholog"].max()) + 0.05]
    ax.plot(lims, lims, '--', color='#7f8c8d', linewidth=1, alpha=0.7, label='y = x')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    from scipy.stats import pearsonr
    r, p = pearsonr(df["Fingerprint"], df["Ortholog"])
    ax.text(0.05, 0.95, f"Pearson r = {r:.3f}\np = {p:.2e}\nn = {len(df)} pairs",
            transform=ax.transAxes, fontsize=10, va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.set_xlabel("Fingerprint Divergence (distribution shape)", fontsize=11)
    ax.set_ylabel("Ortholog-Mapped Divergence (gene-level correlation)", fontsize=11)
    ax.set_title("Method Comparison: Fingerprint vs Ortholog-Mapped Divergence",
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved method comparison plot to {output_path}")


def plot_mean_divergence_rank(
    div_matrix: pd.DataFrame, 
    output_path: Path,
    family_map: Optional[Dict[str, str]] = None,
    family_colors: Optional[Dict[str, str]] = None,
) -> None:
    """Rank species by their mean divergence from all others.
    
    Args:
        div_matrix: Symmetric divergence matrix
        output_path: Path to save the plot
        family_map: Dictionary mapping species names to taxonomic families
        family_colors: Dictionary mapping families to hex colors
    """
    means = div_matrix.replace(1.0, np.nan).mean(axis=1).sort_values()
    n = len(means)

    fig, ax = plt.subplots(figsize=(10, max(6, n*0.3)))
    colors = [_get_family_color(sp, family_map, family_colors) for sp in means.index]
    ax.barh(range(n), means.values, color=colors, edgecolor='white', linewidth=0.5)
    ax.set_yticks(range(n))
    ax.set_yticklabels([_fmt(s) for s in means.index], fontsize=9)
    for i, v in enumerate(means.values):
        ax.text(v + 0.002, i, f"{v:.3f}", va='center', fontsize=8)
    ax.set_xlabel("Mean Divergence from All Other Species", fontsize=11)
    ax.set_title("Species Ranked by Mean Expression Divergence", fontsize=13, fontweight='bold')
    
    if family_map and family_colors:
        legend = [Patch(facecolor=c, label=f) for f, c in family_colors.items()]
        ax.legend(handles=legend, loc='lower right', fontsize=8)
        
    ax.set_xlim(0, means.max() * 1.15)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved mean divergence rank plot to {output_path}")


def plot_species_summary(gene_stats: pd.DataFrame, output_path: Path) -> None:
    """Generate a bar chart of gene counts and expression levels per species.
    
    Args:
        gene_stats: DataFrame from compute_gene_count_overlap
        output_path: Path to save the plot
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, max(8, len(gene_stats)*0.4)))
    
    gene_stats_sorted = gene_stats.sort_values("expressed_genes", ascending=True)
    
    # Gene counts
    y_pos = range(len(gene_stats_sorted))
    ax1.barh(y_pos, gene_stats_sorted["total_genes"], color="#90CAF9", label="Total Genes", height=0.6)
    ax1.barh(y_pos, gene_stats_sorted["expressed_genes"], color="#1565C0", label="Expressed Genes", height=0.6)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels([_fmt(s) for s in gene_stats_sorted["species"]], fontsize=9)
    ax1.set_xlabel("Number of Genes")
    ax1.set_title("Gene Counts per Species", fontsize=14, pad=15)
    ax1.legend()
    
    # Mean expression
    ax2.barh(y_pos, gene_stats_sorted["mean_expression"], color="#4CAF50", height=0.6)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([""] * len(gene_stats_sorted))
    ax2.set_xlabel("Mean Expression (TPM)")
    ax2.set_title("Mean Expression per Species", fontsize=14, pad=15)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved species gene summary plot to {output_path}")


def plot_combined_summary(
    div_matrix: pd.DataFrame, 
    output_path: Path,
    family_map: Optional[Dict[str, str]] = None,
    family_colors: Optional[Dict[str, str]] = None,
    title: str = "Cross-Species Expression Analysis",
) -> None:
    """4-panel combined summary figure showing dendrogram, heatmap, and distribution.
    
    Args:
        div_matrix: Symmetric divergence matrix
        output_path: Path to save the plot
        family_map: Dictionary mapping species names to taxonomic families
        family_colors: Dictionary mapping families to hex colors
        title: Overall figure title
    """
    n = len(div_matrix)
    condensed = squareform(div_matrix.values, checks=False)
    condensed = np.nan_to_num(condensed, nan=0.5)
    link = linkage(condensed, method='ward')
    order = [div_matrix.index[i] for i in leaves_list(link)]

    fig = plt.figure(figsize=(22, 18))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3)

    # A: Dendrogram
    ax1 = fig.add_subplot(gs[0, :])
    labels = [_fmt(s) for s in div_matrix.index]
    dn = dendrogram(link, labels=labels, leaf_rotation=40, leaf_font_size=8,
                     ax=ax1, color_threshold=0.4, above_threshold_color='#7f8c8d')
    if family_map and family_colors:
        for label in ax1.get_xticklabels():
            sp = label.get_text().replace(" ", "_")
            label.set_color(_get_family_color(sp, family_map, family_colors))
        legend = [Patch(facecolor=c, label=f) for f, c in family_colors.items()]
        ax1.legend(handles=legend, loc='upper right', fontsize=8)
        
    ax1.set_title(f"A. Expression-Based Phylogeny ({n} species, Ward Linkage)",
                  fontsize=12, fontweight='bold', loc='left')
    ax1.set_ylabel("Ward Distance")

    # B: Heatmap
    ax2 = fig.add_subplot(gs[1, 0])
    ordered_div = div_matrix.loc[order, order]
    short = [_fmt(s) for s in order]
    sns.heatmap(ordered_div, cmap="YlOrRd", annot=False, ax=ax2,
                xticklabels=short, yticklabels=short,
                vmin=0.25, vmax=0.75, linewidths=0.2,
                cbar_kws={"label": "Divergence", "shrink": 0.8})
    ax2.set_title("B. Divergence Matrix", fontsize=12, fontweight='bold', loc='left')
    plt.setp(ax2.get_xticklabels(), rotation=45, ha='right', fontsize=6)
    plt.setp(ax2.get_yticklabels(), fontsize=6)

    # C: Distribution
    ax3 = fig.add_subplot(gs[1, 1])
    upper = div_matrix.values[np.triu_indices(n, k=1)]
    ax3.hist(upper, bins=30, color='#6c5ce7', edgecolor='white', alpha=0.85)
    ax3.axvline(np.median(upper), color='#e74c3c', ls='--', lw=2,
                label=f'Median: {np.median(upper):.3f}')
    ax3.axvline(np.mean(upper), color='#f39c12', ls='--', lw=2,
                label=f'Mean: {np.mean(upper):.3f}')
    ax3.set_xlabel("Expression Divergence", fontsize=10)
    ax3.set_ylabel("Number of Species Pairs", fontsize=10)
    ax3.set_title("C. Divergence Distribution", fontsize=12, fontweight='bold', loc='left')
    ax3.legend(fontsize=9)

    fig.suptitle(title, fontsize=15, fontweight='bold', y=0.99)
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved combined summary figure to {output_path}")
