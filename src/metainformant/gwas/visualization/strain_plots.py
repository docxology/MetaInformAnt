"""Strain-aware GWAS visualization functions.

Provides publication-quality PCA scatter plots with strain coloring,
dendrograms from kinship matrices, Fst Manhattan plots, and allele
frequency heatmaps for Apis mellifera population group analysis.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

try:
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.offsetbox import AnchoredText

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    from scipy.cluster import hierarchy as sch
    from scipy.spatial.distance import squareform

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# ── Shared strain palette ───────────────────────────────────────────────
STRAIN_PALETTE = {"C": "#2196F3", "I": "#FF9800", "M": "#4CAF50", "R": "#F44336"}
STRAIN_MARKERS = {"C": "o", "I": "s", "M": "D", "R": "^"}
STRAIN_NAMES = {
    "C": "C-lineage (Carniolan)",
    "I": "I-lineage (Italian)",
    "M": "M-lineage (A.m. mellifera)",
    "R": "R-lineage (Russian)",
}


def _extract_strain(sample_id: str) -> str:
    """Extract strain from sample ID: 'C15ITQ' -> 'C'."""
    return sample_id[0].upper() if sample_id and sample_id[0].isalpha() else "?"


def _apply_publication_style(ax: plt.Axes, title: str, xlabel: str, ylabel: str) -> None:
    """Apply consistent publication-grade typography."""
    ax.set_title(title, fontsize=16, fontweight="bold", pad=12)
    ax.set_xlabel(xlabel, fontsize=14, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=14, labelpad=8)
    ax.tick_params(axis="both", which="major", labelsize=12)
    ax.grid(True, alpha=0.3, linewidth=0.5)


def strain_pca_plot(
    pca_result: dict,
    sample_ids: List[str],
    output_path: Optional[Union[str, Path]] = None,
    show_labels: bool = True,
    show_ellipses: bool = True,
) -> Any:
    """Strain-colored PCA scatter with confidence ellipses and scree inset.

    Args:
        pca_result: Dict with 'pcs' (n_samples x n_pc) and 'explained_variance_ratio'
        sample_ids: Sample IDs in order matching PCA rows
        output_path: Path to save figure
        show_labels: Label each point with sample ID
        show_ellipses: Draw 95% confidence ellipses per strain

    Returns:
        matplotlib Figure
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        return None

    pcs = np.array(pca_result.get("pcs", []))
    explained_var = pca_result.get("explained_variance_ratio", [])

    if len(pcs) < 2 or pcs.shape[1] < 2:
        logger.warning("Insufficient PCA components for strain PCA plot")
        return None

    strains = [_extract_strain(sid) for sid in sample_ids]
    n_pcs = min(pcs.shape[1], 4)

    # Create figure: main PC1/PC2 + scree inset
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    # ── Panel 1: PC1 vs PC2 (main) ──────────────────────────────────────
    ax_main = fig.add_subplot(gs[0, 0])
    _plot_strain_scatter(
        ax_main,
        pcs[:, 0],
        pcs[:, 1],
        strains,
        sample_ids,
        explained_var,
        0,
        1,
        show_labels,
        show_ellipses,
    )

    # ── Panel 2: PC1 vs PC3 ─────────────────────────────────────────────
    ax_13 = fig.add_subplot(gs[0, 1])
    if n_pcs >= 3:
        _plot_strain_scatter(
            ax_13,
            pcs[:, 0],
            pcs[:, 2],
            strains,
            sample_ids,
            explained_var,
            0,
            2,
            show_labels=False,
            show_ellipses=show_ellipses,
        )
    else:
        ax_13.text(0.5, 0.5, "PC3 not available", ha="center", va="center", fontsize=12)
        ax_13.set_axis_off()

    # ── Panel 3: PC2 vs PC3 ─────────────────────────────────────────────
    ax_23 = fig.add_subplot(gs[1, 0])
    if n_pcs >= 3:
        _plot_strain_scatter(
            ax_23,
            pcs[:, 1],
            pcs[:, 2],
            strains,
            sample_ids,
            explained_var,
            1,
            2,
            show_labels=False,
            show_ellipses=show_ellipses,
        )
    else:
        ax_23.text(0.5, 0.5, "PC3 not available", ha="center", va="center", fontsize=12)
        ax_23.set_axis_off()

    # ── Panel 4: Scree Plot ─────────────────────────────────────────────
    ax_scree = fig.add_subplot(gs[1, 1])
    if explained_var:
        n_show = min(len(explained_var), 10)
        pcts = [v * 100 for v in explained_var[:n_show]]
        cumulative = [sum(pcts[: i + 1]) for i in range(n_show)]
        ax_scree.bar(
            range(1, n_show + 1),
            pcts,
            color="#5C6BC0",
            alpha=0.8,
            edgecolor="black",
            linewidth=0.5,
        )
        ax_scree.plot(
            range(1, n_show + 1),
            cumulative,
            "o-",
            color="#E74C3C",
            linewidth=2,
            markersize=6,
            label="Cumulative",
        )
        ax_scree.axhline(y=80, color="gray", linestyle="--", alpha=0.5, label="80% threshold")
        _apply_publication_style(ax_scree, "Variance Explained", "Principal Component", "% Variance")
        ax_scree.legend(fontsize=10)
        ax_scree.set_xticks(range(1, n_show + 1))
    else:
        ax_scree.text(0.5, 0.5, "No variance data", ha="center", va="center")

    fig.suptitle(
        "Principal Component Analysis — Colored by Strain",
        fontsize=18,
        fontweight="bold",
        y=0.98,
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved strain PCA plot to {output_path}")

    return fig


def _plot_strain_scatter(
    ax,
    x_vals,
    y_vals,
    strains,
    sample_ids,
    explained_var,
    pc_x,
    pc_y,
    show_labels=True,
    show_ellipses=True,
):
    """Helper: plot strain-colored scatter on given axes."""
    unique_strains = sorted(set(strains))

    for strain in unique_strains:
        mask = [i for i, s in enumerate(strains) if s == strain]
        if not mask:
            continue
        xs = x_vals[mask]
        ys = y_vals[mask]
        color = STRAIN_PALETTE.get(strain, "#888888")
        marker = STRAIN_MARKERS.get(strain, "o")
        label = STRAIN_NAMES.get(strain, strain)

        ax.scatter(
            xs,
            ys,
            c=color,
            marker=marker,
            s=80,
            alpha=0.85,
            edgecolors="black",
            linewidths=0.5,
            label=label,
            zorder=3,
        )

        # 95% confidence ellipse
        if show_ellipses and len(xs) >= 3 and HAS_NUMPY:
            _draw_confidence_ellipse(ax, xs, ys, color, alpha=0.15)

        # Labels
        if show_labels:
            for idx in mask:
                ax.annotate(
                    sample_ids[idx],
                    (x_vals[idx], y_vals[idx]),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=6,
                    alpha=0.7,
                )

    var_x = explained_var[pc_x] * 100 if pc_x < len(explained_var) else 0
    var_y = explained_var[pc_y] * 100 if pc_y < len(explained_var) else 0
    _apply_publication_style(
        ax,
        f"PC{pc_x + 1} vs PC{pc_y + 1}",
        f"PC{pc_x + 1} ({var_x:.1f}%)",
        f"PC{pc_y + 1} ({var_y:.1f}%)",
    )
    ax.legend(fontsize=9, loc="best", framealpha=0.9, edgecolor="black")


def _draw_confidence_ellipse(ax, x, y, color, n_std=2.0, alpha=0.15):
    """Draw a 95% confidence ellipse around points."""
    import numpy as np
    from matplotlib.patches import Ellipse

    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]

    angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    width, height = 2 * n_std * np.sqrt(vals)

    ellipse = Ellipse(
        xy=(np.mean(x), np.mean(y)),
        width=width,
        height=height,
        angle=angle,
        facecolor=color,
        alpha=alpha,
        edgecolor=color,
        linewidth=1.5,
        linestyle="--",
    )
    ax.add_patch(ellipse)


def dendrogram_plot(
    kinship_matrix: list | Any,
    sample_ids: List[str],
    output_path: Optional[Union[str, Path]] = None,
    method: str = "ward",
) -> Any:
    """Hierarchical clustering dendrogram colored by strain.

    Args:
        kinship_matrix: Square kinship matrix (n_samples x n_samples)
        sample_ids: Sample IDs
        output_path: Path to save
        method: Linkage method

    Returns:
        matplotlib Figure
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY or not HAS_SCIPY:
        logger.warning("scipy required for dendrogram plot")
        return None

    K = np.array(kinship_matrix)
    n = K.shape[0]

    # Convert kinship to distance (1 - kinship, clipped to [0, inf))
    D = np.clip(1.0 - K, 0, None)
    np.fill_diagonal(D, 0)

    # Make symmetric
    D = (D + D.T) / 2

    # Condensed distance
    condensed = squareform(D, checks=False)

    # Linkage
    Z = sch.linkage(condensed, method=method)

    # Cophenetic correlation
    coph_corr, _ = sch.cophenet(Z, condensed)

    # Strain colors for leaves
    strains = [_extract_strain(sid) for sid in sample_ids]
    {i: STRAIN_PALETTE.get(strains[i], "#888888") for i in range(n)}

    fig, ax = plt.subplots(figsize=(16, 8))

    # Custom link coloring function
    def link_color_func(k):
        return "#555555"

    sch.dendrogram(
        Z,
        labels=sample_ids,
        ax=ax,
        leaf_rotation=90,
        leaf_font_size=9,
        above_threshold_color="#555555",
        color_threshold=0,
        link_color_func=link_color_func,
    )

    # Color leaf labels by strain
    xlbls = ax.get_xticklabels()
    for lbl in xlbls:
        sid = lbl.get_text()
        strain = _extract_strain(sid)
        lbl.set_color(STRAIN_PALETTE.get(strain, "#888888"))
        lbl.set_fontweight("bold")

    _apply_publication_style(
        ax,
        f"Hierarchical Clustering Dendrogram (method={method})",
        "Sample",
        "Distance (1 - Kinship)",
    )

    # Add cophenetic correlation and strain legend
    info_text = f"Cophenetic r = {coph_corr:.4f}\nn = {n} samples"
    at = AnchoredText(
        info_text,
        loc="upper right",
        prop=dict(size=11, fontweight="bold"),
        frameon=True,
    )
    at.patch.set_boxstyle("round,pad=0.3")
    at.patch.set_facecolor("white")
    at.patch.set_alpha(0.95)
    at.patch.set_edgecolor("black")
    ax.add_artist(at)

    # Strain legend
    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor=STRAIN_PALETTE.get(s, "#888"),
            markersize=12,
            label=STRAIN_NAMES.get(s, s),
        )
        for s in sorted(set(strains))
        if s in STRAIN_PALETTE
    ]
    ax.legend(
        handles=legend_elements,
        loc="upper left",
        fontsize=10,
        framealpha=0.95,
        edgecolor="black",
        title="Strain",
        title_fontsize=12,
    )

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved dendrogram to {output_path}")

    return fig


def kinship_heatmap_clustered(
    kinship_matrix: list | Any,
    sample_ids: List[str],
    output_path: Optional[Union[str, Path]] = None,
) -> Any:
    """Kinship heatmap with hierarchical clustering and strain color bar.

    Args:
        kinship_matrix: Square kinship matrix
        sample_ids: Sample IDs
        output_path: Save path

    Returns:
        matplotlib Figure
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY or not HAS_SCIPY:
        return None

    K = np.array(kinship_matrix)
    n = K.shape[0]

    # Distance and linkage for ordering
    D = np.clip(1.0 - K, 0, None)
    np.fill_diagonal(D, 0)
    D = (D + D.T) / 2
    condensed = squareform(D, checks=False)
    Z = sch.linkage(condensed, method="average")
    order = sch.leaves_list(Z)

    # Reorder matrix
    K_ordered = K[np.ix_(order, order)]
    ordered_ids = [sample_ids[i] for i in order]
    strains_ordered = [_extract_strain(sid) for sid in ordered_ids]

    fig, axes = plt.subplots(1, 2, figsize=(16, 12), gridspec_kw={"width_ratios": [0.03, 1]})

    # Strain color bar
    ax_bar = axes[0]
    strain_colors_arr = [STRAIN_PALETTE.get(s, "#888888") for s in strains_ordered]
    for i, color in enumerate(strain_colors_arr):
        ax_bar.barh(i, 1, color=color, edgecolor="none", height=1)
    ax_bar.set_ylim(-0.5, n - 0.5)
    ax_bar.set_xlim(0, 1)
    ax_bar.set_yticks([])
    ax_bar.set_xticks([])
    ax_bar.set_xlabel("Strain", fontsize=11, fontweight="bold")
    ax_bar.invert_yaxis()

    # Main heatmap
    ax_heat = axes[1]
    im = ax_heat.imshow(
        K_ordered,
        cmap="RdBu_r",
        aspect="equal",
        vmin=np.percentile(K_ordered, 2),
        vmax=np.percentile(K_ordered, 98),
    )

    ax_heat.set_xticks(range(n))
    ax_heat.set_yticks(range(n))
    ax_heat.set_xticklabels(ordered_ids, rotation=90, fontsize=7)
    ax_heat.set_yticklabels(ordered_ids, fontsize=7)

    # Color tick labels by strain
    for i, (xlbl, ylbl) in enumerate(zip(ax_heat.get_xticklabels(), ax_heat.get_yticklabels())):
        color = STRAIN_PALETTE.get(strains_ordered[i], "#888")
        xlbl.set_color(color)
        ylbl.set_color(color)

    plt.colorbar(im, ax=ax_heat, shrink=0.7, label="Kinship Coefficient")
    ax_heat.set_title(
        "Kinship Matrix — Clustered by Genetic Similarity",
        fontsize=16,
        fontweight="bold",
        pad=15,
    )

    # Strain legend
    legend_elements = [
        mpatches.Patch(facecolor=STRAIN_PALETTE.get(s), label=STRAIN_NAMES.get(s, s))
        for s in sorted(set(strains_ordered))
        if s in STRAIN_PALETTE
    ]
    ax_heat.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(1.15, 1.0),
        fontsize=9,
        title="Strain",
        title_fontsize=11,
        framealpha=0.95,
    )

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved clustered kinship heatmap to {output_path}")

    return fig


def fst_manhattan_plot(
    fst_data: Dict[str, List[float]],
    variants_info: List[Dict[str, Any]],
    output_path: Optional[Union[str, Path]] = None,
    top_n_label: int = 5,
) -> Any:
    """Manhattan-style plot with Fst on Y-axis, colored by strain pair.

    Args:
        fst_data: Dict mapping "A_vs_B" -> list of Fst values
        variants_info: Variant position info dicts
        output_path: Save path
        top_n_label: Annotate top N Fst hits

    Returns:
        matplotlib Figure
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        return None

    pair_colors = [
        "#E53935",
        "#1E88E5",
        "#43A047",
        "#FB8C00",
        "#8E24AA",
        "#00ACC1",
        "#D81B60",
        "#5E35B1",
    ]

    n_pairs = len(fst_data)
    fig, axes = plt.subplots(n_pairs, 1, figsize=(16, 4 * n_pairs), sharex=True, squeeze=False)

    for idx, (pair_key, fst_values) in enumerate(sorted(fst_data.items())):
        ax = axes[idx, 0]
        color = pair_colors[idx % len(pair_colors)]

        # Build genomic positions
        chrom_offsets = {}
        current_pos = 0
        xs, ys = [], []

        for v, fst in enumerate(fst_values):
            if v >= len(variants_info):
                break
            vinfo = variants_info[v]
            chrom = str(vinfo.get("chrom", "1"))
            pos = vinfo.get("pos", 0)

            if chrom not in chrom_offsets:
                chrom_offsets[chrom] = current_pos
                current_pos += 100_000_000

            gpos = chrom_offsets[chrom] + pos
            xs.append(gpos)
            ys.append(fst if not math.isnan(fst) else 0)

        ax.scatter(xs, ys, c=color, s=4, alpha=0.5, linewidths=0)

        # Label top hits
        if top_n_label > 0 and ys:
            sorted_idx = sorted(range(len(ys)), key=lambda i: ys[i], reverse=True)
            for rank, vi in enumerate(sorted_idx[:top_n_label]):
                if vi < len(variants_info):
                    vinfo = variants_info[vi]
                    lbl = f"{vinfo.get('chrom', '?')}:{vinfo.get('pos', '?')}"
                    ax.annotate(
                        lbl,
                        (xs[vi], ys[vi]),
                        xytext=(0, 8),
                        textcoords="offset points",
                        fontsize=6,
                        ha="center",
                        color=color,
                        fontweight="bold",
                    )

        pair_label = pair_key.replace("_vs_", " vs ")
        _apply_publication_style(ax, f"Fst: {pair_label}", "", "Fst")

        # Mean Fst annotation
        valid_fst = [f for f in fst_values if not math.isnan(f)]
        if valid_fst:
            mean_fst = sum(valid_fst) / len(valid_fst)
            at = AnchoredText(
                f"Mean Fst = {mean_fst:.4f}\nn = {len(valid_fst)} variants",
                loc="upper right",
                prop=dict(size=10, fontweight="bold"),
                frameon=True,
            )
            at.patch.set_boxstyle("round,pad=0.2")
            at.patch.set_facecolor("white")
            at.patch.set_alpha(0.95)
            at.patch.set_edgecolor("black")
            ax.add_artist(at)

    axes[-1, 0].set_xlabel("Genomic Position", fontsize=14)
    fig.suptitle("Per-Variant Fst Between Strain Pairs", fontsize=18, fontweight="bold", y=1.01)
    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved Fst Manhattan plot to {output_path}")

    return fig


def allele_frequency_heatmap(
    af_by_strain: Dict[str, List[float]],
    variants_info: List[Dict[str, Any]],
    output_path: Optional[Union[str, Path]] = None,
    top_n: int = 100,
) -> Any:
    """Heatmap of allele frequencies with dendrograms and strain coloring.

    Features:
      - Column dendrogram: variants clustered by AF pattern similarity
      - Row dendrogram: strains clustered by AF distance
      - Strain-colored row sidebar (C=blue, I=orange, M=green, R=red)
      - Top-N most differentiated variants selected by cross-strain variance

    Args:
        af_by_strain: Dict strain -> allele frequencies
        variants_info: Variant info dicts
        output_path: Save path
        top_n: Number of top differentially distributed variants to show

    Returns:
        matplotlib Figure
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        return None

    strains = sorted(af_by_strain.keys())
    n_variants = min(len(list(af_by_strain.values())[0]), len(variants_info))

    # Compute variance across strains for each variant to find most differentiated
    af_matrix = np.array([af_by_strain[s][:n_variants] for s in strains], dtype=float)
    var_across_strains = np.nanvar(af_matrix, axis=0)

    # Safely handle completely NaN slices by filling them with 0 variance
    var_across_strains = np.nan_to_num(var_across_strains, nan=0.0)

    top_indices = np.argsort(var_across_strains)[-top_n:][::-1]

    af_subset = af_matrix[:, top_indices]

    # Fill any remaining NaNs to prevent scipy pdist linkage crashes
    af_subset = np.nan_to_num(af_subset, nan=0.0)

    variant_labels = []
    for vi in top_indices:
        if vi < len(variants_info):
            vinfo = variants_info[vi]
            chrom = str(vinfo.get("chrom", "?"))
            if "." in chrom:
                chrom = chrom.split(".")[-2][-2:]
            pos = vinfo.get("pos", "?")
            variant_labels.append(f"{chrom}:{pos}")
        else:
            variant_labels.append(f"v{vi}")

    n_display = len(variant_labels)

    # Try scipy for dendrograms
    try:
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import pdist

        has_scipy = True
    except ImportError:
        has_scipy = False

    if has_scipy and len(strains) >= 2 and n_display >= 2:
        fig = plt.figure(figsize=(max(14, n_display * 0.15), 8))

        # GridSpec: [strain_bar | heatmap] top: col_dendrogram, left: row_dendrogram
        gs = fig.add_gridspec(
            2,
            4,
            width_ratios=[0.08, 0.03, 1, 0.03],
            height_ratios=[0.12, 1],
            hspace=0.02,
            wspace=0.02,
        )

        # Column dendrogram (cluster variants by AF similarity across strains)
        Z_col = linkage(pdist(af_subset.T, metric="euclidean"), method="average")
        ax_col_dend = fig.add_subplot(gs[0, 2])
        dn_col = dendrogram(
            Z_col,
            ax=ax_col_dend,
            orientation="top",
            no_labels=True,
            color_threshold=0,
            above_threshold_color="#666666",
        )
        ax_col_dend.set_axis_off()
        col_order = dn_col["leaves"]

        # Row dendrogram (cluster strains by AF distance)
        if len(strains) >= 2:
            Z_row = linkage(pdist(af_subset, metric="euclidean"), method="average")
            ax_row_dend = fig.add_subplot(gs[1, 0])
            dn_row = dendrogram(
                Z_row,
                ax=ax_row_dend,
                orientation="left",
                no_labels=True,
                color_threshold=0,
                above_threshold_color="#666666",
            )
            ax_row_dend.set_axis_off()
            ax_row_dend.invert_yaxis()
            row_order = dn_row["leaves"]
        else:
            row_order = list(range(len(strains)))

        # Strain color sidebar
        ax_bar = fig.add_subplot(gs[1, 1])
        ordered_strains = [strains[i] for i in row_order]
        for i, s in enumerate(ordered_strains):
            color = STRAIN_PALETTE.get(s, "#888888")
            ax_bar.add_patch(plt.Rectangle((0, i), 1, 1, color=color))
        ax_bar.set_xlim(0, 1)
        ax_bar.set_ylim(0, len(ordered_strains))
        ax_bar.set_axis_off()

        # Reorder matrix
        reordered = af_subset[np.ix_(row_order, col_order)]
        reordered_labels = [variant_labels[i] for i in col_order]

        # Main heatmap
        ax_heat = fig.add_subplot(gs[1, 2])
        im = ax_heat.imshow(
            reordered,
            cmap="YlOrRd",
            aspect="auto",
            vmin=0,
            vmax=1,
            interpolation="nearest",
        )
        ax_heat.set_yticks(range(len(ordered_strains)))
        ax_heat.set_yticklabels([STRAIN_NAMES.get(s, s) for s in ordered_strains], fontsize=11)

        # X labels
        if n_display <= 30:
            ax_heat.set_xticks(range(n_display))
            ax_heat.set_xticklabels(reordered_labels, rotation=90, fontsize=7)
        else:
            step = max(1, n_display // 20)
            ticks = list(range(0, n_display, step))
            ax_heat.set_xticks(ticks)
            ax_heat.set_xticklabels([reordered_labels[i] for i in ticks], rotation=90, fontsize=7)

        # Colorbar
        ax_cbar = fig.add_subplot(gs[1, 3])
        plt.colorbar(im, cax=ax_cbar, label="Allele Frequency")

        fig.suptitle(
            f"Allele Frequency Heatmap — Top {top_n} Differentiated Variants (Clustered)",
            fontsize=16,
            fontweight="bold",
            y=0.99,
        )

        # Strain legend
        from matplotlib.patches import Patch

        legend_patches = [
            Patch(facecolor=STRAIN_PALETTE.get(s, "#888"), label=STRAIN_NAMES.get(s, s))
            for s in sorted(set(ordered_strains))
        ]
        if legend_patches:
            fig.legend(
                handles=legend_patches,
                loc="lower left",
                fontsize=9,
                title="Strain",
                title_fontsize=10,
                bbox_to_anchor=(0.02, 0.02),
            )

    else:
        # Fallback: simple heatmap without dendrograms
        fig, ax = plt.subplots(figsize=(max(12, n_display * 0.15), 4 + len(strains) * 0.5))

        im = ax.imshow(af_subset, cmap="YlOrRd", aspect="auto", vmin=0, vmax=1)
        ax.set_yticks(range(len(strains)))
        ax.set_yticklabels([STRAIN_NAMES.get(s, s) for s in strains], fontsize=11)

        if n_display <= 30:
            ax.set_xticks(range(n_display))
            ax.set_xticklabels(variant_labels, rotation=90, fontsize=7)
        else:
            step = max(1, n_display // 20)
            ticks = list(range(0, n_display, step))
            ax.set_xticks(ticks)
            ax.set_xticklabels([variant_labels[i] for i in ticks], rotation=90, fontsize=7)

        plt.colorbar(im, ax=ax, shrink=0.6, label="Allele Frequency")
        ax.set_title(
            f"Allele Frequency Heatmap — Top {top_n} Differentiated Variants",
            fontsize=16,
            fontweight="bold",
            pad=12,
        )
        ax.set_xlabel("Variant", fontsize=14)

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved allele frequency heatmap to {output_path}")

    return fig
