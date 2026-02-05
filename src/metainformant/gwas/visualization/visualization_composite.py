"""Composite multi-panel visualization for GWAS.

This module provides composite panel plots that combine multiple GWAS
visualizations into unified summary figures for publication and reporting.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    plt = None

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None


def gwas_summary_panel(
    assoc_results: List[Dict[str, Any]],
    pca_data: Optional[Dict[str, Any]] = None,
    kinship_matrix: Optional[Any] = None,
    output_file: Optional[Union[str, Path]] = None,
    significance_threshold: float = 5e-8,
    title: str = "GWAS Summary",
) -> Dict[str, Any]:
    """Create a 2x2 GWAS summary panel with Manhattan, QQ, PCA, and kinship plots.

    Generates a composite figure with four panels:
      - Top-left: Manhattan plot of association results.
      - Top-right: QQ plot of observed vs expected p-values.
      - Bottom-left: PCA scatter of PC1 vs PC2 (or placeholder if no data).
      - Bottom-right: Kinship heatmap (or placeholder if no data).

    Args:
        assoc_results: List of association result dicts, each containing at
            minimum ``"p_value"`` and ``"chrom"`` keys.  Optional keys include
            ``"pos"`` (genomic position) and ``"variant_id"``.
        pca_data: Optional dictionary with ``"pcs"`` (2-D array of principal
            components) and optionally ``"explained_variance_ratio"``.
        kinship_matrix: Optional 2-D array-like kinship matrix.
        output_file: Path to save the composite figure.
        significance_threshold: Genome-wide significance threshold for the
            Manhattan plot horizontal line.
        title: Super-title for the figure.

    Returns:
        Dictionary with ``"status"``, ``"output_path"``, and
        ``"panels_generated"`` keys.
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib/numpy not available, skipping GWAS summary panel")
        return {"status": "skipped", "output_path": None, "panels_generated": 0}

    try:
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        panels_generated = 0

        # --- Top-left: Manhattan plot ---
        ax_manhattan = axes[0, 0]
        if assoc_results:
            _draw_manhattan(ax_manhattan, assoc_results, significance_threshold)
            panels_generated += 1
        else:
            ax_manhattan.text(
                0.5,
                0.5,
                "No data available",
                ha="center",
                va="center",
                fontsize=14,
                color="gray",
                transform=ax_manhattan.transAxes,
            )
            ax_manhattan.set_title("Manhattan Plot")

        # --- Top-right: QQ plot ---
        ax_qq = axes[0, 1]
        if assoc_results:
            _draw_qq(ax_qq, assoc_results)
            panels_generated += 1
        else:
            ax_qq.text(
                0.5,
                0.5,
                "No data available",
                ha="center",
                va="center",
                fontsize=14,
                color="gray",
                transform=ax_qq.transAxes,
            )
            ax_qq.set_title("QQ Plot")

        # --- Bottom-left: PCA scatter ---
        ax_pca = axes[1, 0]
        if pca_data is not None and "pcs" in pca_data:
            _draw_pca(ax_pca, pca_data)
            panels_generated += 1
        else:
            ax_pca.text(
                0.5,
                0.5,
                "No data available",
                ha="center",
                va="center",
                fontsize=14,
                color="gray",
                transform=ax_pca.transAxes,
            )
            ax_pca.set_title("PCA")

        # --- Bottom-right: Kinship heatmap ---
        ax_kinship = axes[1, 1]
        if kinship_matrix is not None:
            _draw_kinship(ax_kinship, kinship_matrix, fig)
            panels_generated += 1
        else:
            ax_kinship.text(
                0.5,
                0.5,
                "No data available",
                ha="center",
                va="center",
                fontsize=14,
                color="gray",
                transform=ax_kinship.transAxes,
            )
            ax_kinship.set_title("Kinship Heatmap")

        fig.suptitle(title, fontsize=16, y=0.98)
        fig.tight_layout(rect=[0, 0, 1, 0.96])

        output_path_str: Optional[str] = None
        if output_file:
            out = Path(output_file)
            out.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out, dpi=300, bbox_inches="tight")
            output_path_str = str(out)
            logger.info(f"Saved GWAS summary panel to {out}")

        plt.close(fig)
        return {
            "status": "success",
            "output_path": output_path_str,
            "panels_generated": panels_generated,
        }

    except Exception as e:
        logger.error(f"Error creating GWAS summary panel: {e}")
        return {"status": "failed", "output_path": None, "panels_generated": 0}


def population_structure_panel(
    pca_data: Dict[str, Any],
    kinship_matrix: Any,
    metadata: Optional[Dict[str, Dict]] = None,
    output_file: Optional[Union[str, Path]] = None,
    title: str = "Population Structure",
) -> Dict[str, Any]:
    """Create a 2x2 population structure panel.

    Generates a composite figure with four panels:
      - Top-left: PCA scatter coloured by population (from metadata).
      - Top-right: Scree plot of explained variance.
      - Bottom-left: Kinship heatmap.
      - Bottom-right: Population counts bar chart.

    Args:
        pca_data: Dictionary with ``"pcs"`` (2-D array), optionally
            ``"explained_variance_ratio"`` and ``"sample_ids"``.
        kinship_matrix: 2-D array-like kinship matrix.
        metadata: Optional mapping from sample id to annotation dict.
            Expected to have a ``"population"`` field per sample for colouring.
        output_file: Path to save the composite figure.
        title: Super-title for the figure.

    Returns:
        Dictionary with ``"status"``, ``"output_path"``, and
        ``"panels_generated"`` keys.
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib/numpy not available, skipping population structure panel")
        return {"status": "skipped", "output_path": None, "panels_generated": 0}

    try:
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        panels_generated = 0

        pcs = np.asarray(pca_data.get("pcs", []), dtype=np.float64)
        explained = pca_data.get("explained_variance_ratio", [])
        sample_ids = pca_data.get("sample_ids", [])

        # --- Top-left: PCA coloured by population ---
        ax_pca = axes[0, 0]
        if pcs.ndim == 2 and pcs.shape[0] > 0 and pcs.shape[1] >= 2:
            pop_labels: Optional[List[str]] = None
            if metadata is not None and sample_ids:
                pop_labels = [str(metadata.get(sid, {}).get("population", "unknown")) for sid in sample_ids]

            if pop_labels is not None:
                unique_pops = sorted(set(pop_labels))
                cmap_colors = plt.cm.tab10(np.linspace(0, 1, max(len(unique_pops), 1)))
                pop_to_color = {p: cmap_colors[i] for i, p in enumerate(unique_pops)}
                for pop in unique_pops:
                    mask = np.array([pl == pop for pl in pop_labels])
                    ax_pca.scatter(
                        pcs[mask, 0],
                        pcs[mask, 1],
                        label=pop,
                        alpha=0.7,
                        s=40,
                        c=[pop_to_color[pop]],
                    )
                ax_pca.legend(title="Population", fontsize=8, loc="best")
            else:
                ax_pca.scatter(pcs[:, 0], pcs[:, 1], alpha=0.7, s=40, c="steelblue")

            var_label_x = f" ({explained[0] * 100:.1f}%)" if len(explained) > 0 else ""
            var_label_y = f" ({explained[1] * 100:.1f}%)" if len(explained) > 1 else ""
            ax_pca.set_xlabel(f"PC1{var_label_x}", fontsize=10)
            ax_pca.set_ylabel(f"PC2{var_label_y}", fontsize=10)
            ax_pca.set_title("PCA by Population", fontsize=12)
            ax_pca.grid(True, alpha=0.3)
            panels_generated += 1
        else:
            ax_pca.text(
                0.5,
                0.5,
                "No data available",
                ha="center",
                va="center",
                fontsize=14,
                color="gray",
                transform=ax_pca.transAxes,
            )
            ax_pca.set_title("PCA by Population")

        # --- Top-right: Scree plot ---
        ax_scree = axes[0, 1]
        if explained and len(explained) > 0:
            explained_pct = [v * 100 if v <= 1.0 else v for v in explained]
            components = list(range(1, len(explained_pct) + 1))
            ax_scree.bar(
                components,
                explained_pct,
                color="skyblue",
                edgecolor="navy",
                alpha=0.7,
            )
            cumulative = np.cumsum(explained_pct)
            ax_scree.plot(
                components,
                cumulative,
                "r-o",
                linewidth=2,
                markersize=4,
                label="Cumulative",
            )
            ax_scree.set_xlabel("Principal Component", fontsize=10)
            ax_scree.set_ylabel("Explained Variance (%)", fontsize=10)
            ax_scree.set_title("Scree Plot", fontsize=12)
            ax_scree.set_xticks(components)
            ax_scree.grid(True, alpha=0.3, axis="y")
            ax_scree.legend(fontsize=8)
            panels_generated += 1
        else:
            ax_scree.text(
                0.5,
                0.5,
                "No data available",
                ha="center",
                va="center",
                fontsize=14,
                color="gray",
                transform=ax_scree.transAxes,
            )
            ax_scree.set_title("Scree Plot")

        # --- Bottom-left: Kinship heatmap ---
        ax_kinship = axes[1, 0]
        if kinship_matrix is not None:
            _draw_kinship(ax_kinship, kinship_matrix, fig)
            panels_generated += 1
        else:
            ax_kinship.text(
                0.5,
                0.5,
                "No data available",
                ha="center",
                va="center",
                fontsize=14,
                color="gray",
                transform=ax_kinship.transAxes,
            )
            ax_kinship.set_title("Kinship Heatmap")

        # --- Bottom-right: Population counts bar chart ---
        ax_bar = axes[1, 1]
        if metadata is not None and sample_ids:
            pop_counts: Dict[str, int] = {}
            for sid in sample_ids:
                pop = str(metadata.get(sid, {}).get("population", "unknown"))
                pop_counts[pop] = pop_counts.get(pop, 0) + 1

            sorted_pops = sorted(pop_counts.keys())
            counts = [pop_counts[p] for p in sorted_pops]
            bar_colors = plt.cm.tab10(np.linspace(0, 1, max(len(sorted_pops), 1)))

            bars = ax_bar.bar(
                sorted_pops,
                counts,
                color=bar_colors[: len(sorted_pops)],
                edgecolor="black",
                alpha=0.8,
            )
            for bar, count in zip(bars, counts):
                ax_bar.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    bar.get_height() + 0.3,
                    str(count),
                    ha="center",
                    va="bottom",
                    fontsize=9,
                )
            ax_bar.set_xlabel("Population", fontsize=10)
            ax_bar.set_ylabel("Sample Count", fontsize=10)
            ax_bar.set_title("Population Sample Counts", fontsize=12)
            ax_bar.grid(True, alpha=0.3, axis="y")
            panels_generated += 1
        else:
            ax_bar.text(
                0.5,
                0.5,
                "No data available",
                ha="center",
                va="center",
                fontsize=14,
                color="gray",
                transform=ax_bar.transAxes,
            )
            ax_bar.set_title("Population Sample Counts")

        fig.suptitle(title, fontsize=16, y=0.98)
        fig.tight_layout(rect=[0, 0, 1, 0.96])

        output_path_str: Optional[str] = None
        if output_file:
            out = Path(output_file)
            out.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out, dpi=300, bbox_inches="tight")
            output_path_str = str(out)
            logger.info(f"Saved population structure panel to {out}")

        plt.close(fig)
        return {
            "status": "success",
            "output_path": output_path_str,
            "panels_generated": panels_generated,
        }

    except Exception as e:
        logger.error(f"Error creating population structure panel: {e}")
        return {"status": "failed", "output_path": None, "panels_generated": 0}


def top_hit_detail_panel(
    assoc_results: List[Dict[str, Any]],
    genotypes: Optional[List[List[int]]] = None,
    phenotypes: Optional[List[float]] = None,
    variant_index: int = 0,
    output_file: Optional[Union[str, Path]] = None,
    title: str = "Top Hit Detail",
) -> Dict[str, Any]:
    """Create a 1x3 detail panel for a top GWAS hit.

    Generates a composite figure with three panels:
      - Left: Regional scatter of -log10(p) for variants within 500 kb of
        the top hit.
      - Centre: Boxplot of phenotype values grouped by genotype at the
        specified variant index.
      - Right: Text display of variant annotation details.

    Args:
        assoc_results: List of association result dicts.  Each dict should
            contain ``"p_value"``, ``"chrom"``, ``"pos"``, and optionally
            ``"variant_id"``, ``"beta"``, ``"se"``.
        genotypes: Optional list of genotype vectors (one per variant).
            The vector at ``variant_index`` is used for the boxplot.
        phenotypes: Optional list of phenotype values (one per sample).
        variant_index: Index into ``assoc_results`` for the top hit.
        output_file: Path to save the composite figure.
        title: Super-title for the figure.

    Returns:
        Dictionary with ``"status"``, ``"output_path"``, and
        ``"panels_generated"`` keys.
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib/numpy not available, skipping top hit detail panel")
        return {"status": "skipped", "output_path": None, "panels_generated": 0}

    try:
        if not assoc_results:
            return {"status": "failed", "output_path": None, "panels_generated": 0}

        # Clamp variant_index
        variant_index = min(variant_index, len(assoc_results) - 1)
        top_hit = assoc_results[variant_index]

        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        panels_generated = 0

        # --- Left: Regional plot ---
        ax_regional = axes[0]
        top_chrom = top_hit.get("chrom", top_hit.get("CHROM", ""))
        top_pos = top_hit.get("pos", top_hit.get("POS", 0))
        window = 500_000  # 500 kb window

        regional_variants = [
            r
            for r in assoc_results
            if r.get("chrom", r.get("CHROM", "")) == top_chrom
            and abs(r.get("pos", r.get("POS", 0)) - top_pos) <= window
        ]

        if regional_variants:
            reg_positions = [r.get("pos", r.get("POS", 0)) for r in regional_variants]
            reg_pvalues = [-np.log10(max(r.get("p_value", 1.0), 1e-300)) for r in regional_variants]
            ax_regional.scatter(
                reg_positions,
                reg_pvalues,
                alpha=0.7,
                s=30,
                c="steelblue",
            )
            # Highlight the top hit
            top_p_log = -np.log10(max(top_hit.get("p_value", 1.0), 1e-300))
            ax_regional.scatter(
                [top_pos],
                [top_p_log],
                s=100,
                c="red",
                zorder=5,
                edgecolors="black",
                linewidths=1,
                label="Top hit",
            )
            ax_regional.axhline(
                y=-np.log10(5e-8),
                color="red",
                linestyle="--",
                alpha=0.5,
                label="GW significance",
            )
            ax_regional.set_xlabel(f"Position on chr{top_chrom}", fontsize=10)
            ax_regional.set_ylabel("-log10(p-value)", fontsize=10)
            ax_regional.set_title("Regional Plot (+/- 500kb)", fontsize=12)
            ax_regional.legend(fontsize=8)
            ax_regional.grid(True, alpha=0.3)
            panels_generated += 1
        else:
            ax_regional.text(
                0.5,
                0.5,
                "No data available",
                ha="center",
                va="center",
                fontsize=14,
                color="gray",
                transform=ax_regional.transAxes,
            )
            ax_regional.set_title("Regional Plot")

        # --- Centre: Genotype-phenotype boxplot ---
        ax_boxplot = axes[1]
        if genotypes is not None and phenotypes is not None and variant_index < len(genotypes):
            geno_vec = genotypes[variant_index]
            pheno_arr = np.asarray(phenotypes, dtype=np.float64)
            geno_arr = np.asarray(geno_vec, dtype=int)

            groups: Dict[int, List[float]] = {}
            for g, p in zip(geno_arr, pheno_arr):
                groups.setdefault(int(g), []).append(float(p))

            sorted_genos = sorted(groups.keys())
            box_data = [groups[g] for g in sorted_genos]
            geno_labels = [f"Genotype {g}" for g in sorted_genos]

            bp = ax_boxplot.boxplot(
                box_data,
                labels=geno_labels,
                patch_artist=True,
            )
            colors = ["#a6cee3", "#b2df8a", "#fb9a99"]
            for i, patch in enumerate(bp["boxes"]):
                patch.set_facecolor(colors[i % len(colors)])

            ax_boxplot.set_xlabel("Genotype", fontsize=10)
            ax_boxplot.set_ylabel("Phenotype", fontsize=10)
            ax_boxplot.set_title("Genotype-Phenotype", fontsize=12)
            ax_boxplot.grid(True, alpha=0.3, axis="y")
            panels_generated += 1
        else:
            ax_boxplot.text(
                0.5,
                0.5,
                "No data available",
                ha="center",
                va="center",
                fontsize=14,
                color="gray",
                transform=ax_boxplot.transAxes,
            )
            ax_boxplot.set_title("Genotype-Phenotype")

        # --- Right: Variant info text ---
        ax_info = axes[2]
        ax_info.axis("off")
        variant_id = top_hit.get("variant_id", top_hit.get("SNP", "N/A"))
        info_lines = [
            f"Variant ID:  {variant_id}",
            f"Chromosome:  {top_chrom}",
            f"Position:    {top_pos:,}",
            (
                f"P-value:     {top_hit.get('p_value', 'N/A'):.2e}"
                if isinstance(top_hit.get("p_value"), (int, float))
                else f"P-value:     {top_hit.get('p_value', 'N/A')}"
            ),
            f"Beta:        {top_hit.get('beta', 'N/A')}",
            f"SE:          {top_hit.get('se', 'N/A')}",
        ]
        info_text = "\n".join(info_lines)
        ax_info.text(
            0.1,
            0.5,
            info_text,
            ha="left",
            va="center",
            fontsize=12,
            fontfamily="monospace",
            transform=ax_info.transAxes,
            bbox=dict(boxstyle="round,pad=0.8", facecolor="lightyellow", alpha=0.9),
        )
        ax_info.set_title("Variant Information", fontsize=12)
        panels_generated += 1

        fig.suptitle(title, fontsize=16, y=0.98)
        fig.tight_layout(rect=[0, 0, 1, 0.94])

        output_path_str: Optional[str] = None
        if output_file:
            out = Path(output_file)
            out.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out, dpi=300, bbox_inches="tight")
            output_path_str = str(out)
            logger.info(f"Saved top hit detail panel to {out}")

        plt.close(fig)
        return {
            "status": "success",
            "output_path": output_path_str,
            "panels_generated": panels_generated,
        }

    except Exception as e:
        logger.error(f"Error creating top hit detail panel: {e}")
        return {"status": "failed", "output_path": None, "panels_generated": 0}


# ---------------------------------------------------------------------------
# Internal drawing helpers
# ---------------------------------------------------------------------------


def _draw_manhattan(
    ax: Any,
    assoc_results: List[Dict[str, Any]],
    significance_threshold: float = 5e-8,
) -> None:
    """Draw a Manhattan plot on the given axes."""
    chrom_order: List[str] = []
    chrom_positions: Dict[str, List[int]] = {}
    chrom_pvalues: Dict[str, List[float]] = {}

    for r in assoc_results:
        chrom = str(r.get("chrom", r.get("CHROM", "1")))
        pos = int(r.get("pos", r.get("POS", 0)))
        pval = r.get("p_value", 1.0)
        if pval is None or pval <= 0:
            continue
        if chrom not in chrom_positions:
            chrom_order.append(chrom)
            chrom_positions[chrom] = []
            chrom_pvalues[chrom] = []
        chrom_positions[chrom].append(pos)
        chrom_pvalues[chrom].append(pval)

    # Sort chromosomes numerically where possible
    def _chrom_sort_key(c: str) -> Tuple[int, str]:
        try:
            return (0, str(int(c)).zfill(5))
        except ValueError:
            return (1, c)

    chrom_order.sort(key=_chrom_sort_key)

    colors = ["#1f77b4", "#aec7e8"]
    offset = 0
    tick_positions: List[float] = []
    tick_labels: List[str] = []

    for idx, chrom in enumerate(chrom_order):
        positions = np.array(chrom_positions[chrom])
        pvalues = np.array(chrom_pvalues[chrom])
        neg_log_p = -np.log10(np.maximum(pvalues, 1e-300))
        x_vals = positions + offset

        ax.scatter(
            x_vals,
            neg_log_p,
            c=colors[idx % 2],
            s=8,
            alpha=0.7,
            edgecolors="none",
        )

        mid = offset + (positions.max() + positions.min()) / 2 if len(positions) > 0 else offset
        tick_positions.append(mid)
        tick_labels.append(chrom)
        if len(positions) > 0:
            offset += positions.max() + 1_000_000

    # Significance line
    sig_line = -np.log10(significance_threshold)
    ax.axhline(y=sig_line, color="red", linestyle="--", alpha=0.7, linewidth=1)

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=7, rotation=45)
    ax.set_xlabel("Chromosome", fontsize=10)
    ax.set_ylabel("-log10(p-value)", fontsize=10)
    ax.set_title("Manhattan Plot", fontsize=12)
    ax.grid(True, alpha=0.2, axis="y")


def _draw_qq(ax: Any, assoc_results: List[Dict[str, Any]]) -> None:
    """Draw a QQ plot on the given axes."""
    p_values = [
        r.get("p_value", 1.0) for r in assoc_results if r.get("p_value") is not None and r.get("p_value", 1.0) > 0
    ]
    if not p_values:
        ax.text(
            0.5,
            0.5,
            "No valid p-values",
            ha="center",
            va="center",
            fontsize=14,
            color="gray",
            transform=ax.transAxes,
        )
        ax.set_title("QQ Plot")
        return

    p_values = sorted(p_values)
    n = len(p_values)
    expected = [-np.log10((i + 1) / (n + 1)) for i in range(n)]
    observed = [-np.log10(max(p, 1e-300)) for p in p_values]

    ax.scatter(expected, observed, alpha=0.6, s=12, c="steelblue")

    max_val = max(max(expected), max(observed))
    ax.plot([0, max_val], [0, max_val], "r--", alpha=0.7, linewidth=1)

    ax.set_xlabel("Expected -log10(p)", fontsize=10)
    ax.set_ylabel("Observed -log10(p)", fontsize=10)
    ax.set_title("QQ Plot", fontsize=12)
    ax.grid(True, alpha=0.3)


def _draw_pca(ax: Any, pca_data: Dict[str, Any]) -> None:
    """Draw a PCA scatter on the given axes."""
    pcs = np.asarray(pca_data["pcs"], dtype=np.float64)
    if pcs.ndim != 2 or pcs.shape[1] < 2:
        ax.text(
            0.5,
            0.5,
            "Insufficient PCA data",
            ha="center",
            va="center",
            fontsize=14,
            color="gray",
            transform=ax.transAxes,
        )
        ax.set_title("PCA")
        return

    explained = pca_data.get("explained_variance_ratio", [])
    ax.scatter(pcs[:, 0], pcs[:, 1], alpha=0.7, s=30, c="steelblue")

    var_x = f" ({explained[0] * 100:.1f}%)" if len(explained) > 0 else ""
    var_y = f" ({explained[1] * 100:.1f}%)" if len(explained) > 1 else ""
    ax.set_xlabel(f"PC1{var_x}", fontsize=10)
    ax.set_ylabel(f"PC2{var_y}", fontsize=10)
    ax.set_title("PCA", fontsize=12)
    ax.grid(True, alpha=0.3)


def _draw_kinship(ax: Any, kinship_matrix: Any, fig: Any) -> None:
    """Draw a kinship heatmap on the given axes."""
    km = np.asarray(kinship_matrix, dtype=np.float64)
    im = ax.imshow(km, cmap="RdYlBu_r", aspect="equal")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title("Kinship Heatmap", fontsize=12)
    ax.set_xlabel("Sample", fontsize=10)
    ax.set_ylabel("Sample", fontsize=10)
