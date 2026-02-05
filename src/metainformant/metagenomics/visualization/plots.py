"""Metagenomics visualization functions.

Provides specialized plots for metagenomic data analysis including
Krona-style hierarchical taxonomy charts, stacked bar composition plots,
rarefaction curves, ordination (PCoA/NMDS), alpha diversity boxplots,
and abundance heatmaps.

All functions accept data in standard Python data structures and
produce matplotlib figures. Figures are saved to disk when an output_path
is provided, or returned for interactive display.
"""

from __future__ import annotations

import math
import os
import random
from collections import defaultdict
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "META_"

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    logger.warning("matplotlib not available; visualization functions will return data only")

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False


def _ensure_matplotlib() -> None:
    """Raise ImportError if matplotlib is not available."""
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib is required for visualization. Install with: uv pip install matplotlib")


def _save_or_show(fig: Any, output_path: str | Path | None) -> Any:
    """Save figure to file or return it."""
    if output_path is not None:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(path), dpi=150, bbox_inches="tight")
        plt.close(fig)
        logger.info(f"Saved plot to {path}")
        return str(path)
    return fig


def plot_krona_chart(
    taxonomy_data: dict[str, float],
    output_path: str | Path | None = None,
    title: str = "Taxonomic Composition",
    max_depth: int = 6,
    min_fraction: float = 0.01,
) -> Any:
    """Create a hierarchical sunburst chart (Krona-style) of taxonomy data.

    Displays taxonomic composition as concentric rings, with each ring
    representing a taxonomic level (domain -> phylum -> class -> etc).

    Args:
        taxonomy_data: Dict mapping taxonomic lineage strings (semicolon-separated)
            to abundance values. Example: {"Bacteria;Proteobacteria;Gammaproteobacteria": 0.35}
        output_path: Path to save figure. If None, returns figure object.
        title: Plot title.
        max_depth: Maximum taxonomic depth to display.
        min_fraction: Minimum fraction to include (filter small taxa).

    Returns:
        Path to saved file or matplotlib figure.
    """
    _ensure_matplotlib()

    if not taxonomy_data:
        logger.warning("Empty taxonomy data provided")
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.text(0.5, 0.5, "No taxonomy data", ha="center", va="center", fontsize=14)
        ax.set_axis_off()
        return _save_or_show(fig, output_path)

    # Normalize abundances
    total = sum(taxonomy_data.values())
    if total <= 0:
        total = 1.0

    normalized = {k: v / total for k, v in taxonomy_data.items() if v / total >= min_fraction}

    # Build hierarchical structure
    tree: dict[str, Any] = {}
    for lineage, abundance in normalized.items():
        levels = [l.strip() for l in lineage.split(";") if l.strip()][:max_depth]
        node = tree
        for level in levels:
            if level not in node:
                node[level] = {"_abundance": 0.0, "_children": {}}
            node[level]["_abundance"] += abundance
            node = node[level]["_children"]

    # Generate sunburst rings
    fig, ax = plt.subplots(figsize=(12, 12), subplot_kw={"projection": "polar"})

    cmap = plt.cm.Set3
    color_idx = 0

    def _draw_ring(
        node: dict[str, Any],
        depth: int,
        start_angle: float,
        end_angle: float,
        parent_color: str | None = None,
    ) -> None:
        nonlocal color_idx
        if depth >= max_depth or not node:
            return

        inner_r = 0.3 + depth * 0.15
        outer_r = inner_r + 0.14

        items = sorted(node.items(), key=lambda x: x[1]["_abundance"], reverse=True)
        current_angle = start_angle

        for name, data in items:
            frac = data["_abundance"] / total if total > 0 else 0
            if frac < min_fraction:
                continue

            span = (end_angle - start_angle) * (frac / max(sum(d["_abundance"] for _, d in items), 1e-10) * total)
            span = min(span, end_angle - current_angle)

            if span <= 0:
                continue

            color = cmap(color_idx % 12 / 12.0)
            color_idx += 1

            # Draw wedge
            theta = [current_angle + i * span / 50 for i in range(51)]
            r_inner = [inner_r] * 51
            r_outer = [outer_r] * 51

            ax.fill_between(theta, r_inner, r_outer, alpha=0.7, color=color, edgecolor="white", linewidth=0.5)

            # Label
            mid_angle = current_angle + span / 2
            if span > 0.15:
                label_r = (inner_r + outer_r) / 2
                ax.text(
                    mid_angle,
                    label_r,
                    name[:15],
                    ha="center",
                    va="center",
                    fontsize=max(5, min(8, span * 20)),
                    rotation=0,
                )

            _draw_ring(data["_children"], depth + 1, current_angle, current_angle + span, color)
            current_angle += span

    _draw_ring(tree, 0, 0, 2 * math.pi)

    ax.set_title(title, fontsize=14, pad=20)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.grid(False)

    return _save_or_show(fig, output_path)


def plot_stacked_bar(
    abundance_data: dict[str, dict[str, float]],
    level: str = "phylum",
    output_path: str | Path | None = None,
    title: str | None = None,
    top_n: int = 15,
    sort_samples: bool = True,
) -> Any:
    """Create stacked bar chart of taxonomic composition across samples.

    Args:
        abundance_data: Dict mapping sample IDs to {taxon: abundance} dicts.
        level: Taxonomic level label for the legend title.
        output_path: Path to save figure.
        title: Plot title (default: auto-generated).
        top_n: Show top N most abundant taxa; rest grouped as "Other".
        sort_samples: Sort samples by dominant taxon.

    Returns:
        Path to saved file or matplotlib figure.
    """
    _ensure_matplotlib()

    if not abundance_data:
        logger.warning("Empty abundance data provided")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        return _save_or_show(fig, output_path)

    # Collect all taxa and compute totals
    taxa_totals: dict[str, float] = defaultdict(float)
    for sample_data in abundance_data.values():
        for taxon, abundance in sample_data.items():
            taxa_totals[taxon] += abundance

    # Top N taxa
    sorted_taxa = sorted(taxa_totals.items(), key=lambda x: x[1], reverse=True)
    top_taxa = [t[0] for t in sorted_taxa[:top_n]]

    sample_names = list(abundance_data.keys())
    if sort_samples and sample_names:
        sample_names.sort(
            key=lambda s: max(abundance_data[s].values()) if abundance_data[s] else 0,
            reverse=True,
        )

    n_samples = len(sample_names)
    fig, ax = plt.subplots(figsize=(max(8, n_samples * 0.5), 6))

    cmap = plt.cm.tab20 if len(top_taxa) > 10 else plt.cm.Set3
    colors = [cmap(i / max(len(top_taxa), 1)) for i in range(len(top_taxa) + 1)]

    x = list(range(n_samples))
    bottom = [0.0] * n_samples

    for i, taxon in enumerate(top_taxa):
        values = []
        for sample in sample_names:
            sample_total = sum(abundance_data[sample].values()) or 1.0
            values.append(abundance_data[sample].get(taxon, 0.0) / sample_total)
        ax.bar(x, values, bottom=bottom, label=taxon, color=colors[i], width=0.8)
        bottom = [b + v for b, v in zip(bottom, values)]

    # "Other" category
    other_values = [1.0 - b for b in bottom]
    if any(v > 0.001 for v in other_values):
        ax.bar(x, other_values, bottom=bottom, label="Other", color="#cccccc", width=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(sample_names, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Relative Abundance")
    ax.set_title(title or f"Taxonomic Composition ({level.capitalize()} Level)")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8, title=level.capitalize())
    ax.set_ylim(0, 1.05)

    fig.tight_layout()
    return _save_or_show(fig, output_path)


def plot_rarefaction_curves(
    otu_table: dict[str, dict[str, int]],
    step: int = 100,
    output_path: str | Path | None = None,
    title: str = "Rarefaction Curves",
    n_iterations: int = 10,
    seed: int = 42,
) -> Any:
    """Plot rarefaction curves showing species richness vs sequencing depth.

    Subsamples reads at increasing depths and counts observed species
    to assess whether sequencing depth is sufficient.

    Args:
        otu_table: Dict mapping sample IDs to {OTU_id: count} dicts.
        step: Increment for subsampling depths.
        output_path: Path to save figure.
        title: Plot title.
        n_iterations: Iterations per depth for averaging.
        seed: Random seed for reproducibility.

    Returns:
        Path to saved file or matplotlib figure.
    """
    _ensure_matplotlib()

    rng = random.Random(seed)

    if not otu_table:
        logger.warning("Empty OTU table provided")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        return _save_or_show(fig, output_path)

    fig, ax = plt.subplots(figsize=(10, 6))
    cmap = plt.cm.tab10

    for idx, (sample_id, otu_counts) in enumerate(otu_table.items()):
        # Expand OTU counts into read pool
        read_pool: list[str] = []
        for otu_id, count in otu_counts.items():
            read_pool.extend([otu_id] * count)

        total_reads = len(read_pool)
        if total_reads == 0:
            continue

        depths = list(range(step, total_reads + 1, step))
        if depths and depths[-1] != total_reads:
            depths.append(total_reads)

        mean_richness: list[float] = []
        for depth in depths:
            richness_sum = 0
            for _ in range(n_iterations):
                subsample = rng.sample(read_pool, min(depth, total_reads))
                richness_sum += len(set(subsample))
            mean_richness.append(richness_sum / n_iterations)

        color = cmap(idx % 10)
        ax.plot(depths, mean_richness, label=sample_id, color=color, linewidth=1.5)

    ax.set_xlabel("Sequencing Depth (reads)")
    ax.set_ylabel("Observed Species")
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    fig.tight_layout()

    return _save_or_show(fig, output_path)


def plot_ordination(
    distance_matrix: dict[str, dict[str, float]],
    method: str = "PCoA",
    output_path: str | Path | None = None,
    title: str | None = None,
    groups: dict[str, str] | None = None,
) -> Any:
    """Plot ordination (PCoA or classical MDS) from a distance matrix.

    Performs dimensionality reduction on a pairwise distance matrix and
    plots samples in 2D space.

    Args:
        distance_matrix: Nested dict of pairwise distances {sample_i: {sample_j: dist}}.
        method: Ordination method ('PCoA' or 'NMDS').
        output_path: Path to save figure.
        title: Plot title.
        groups: Optional dict mapping sample IDs to group labels for coloring.

    Returns:
        Path to saved file or matplotlib figure.
    """
    _ensure_matplotlib()

    samples = sorted(distance_matrix.keys())
    n = len(samples)

    if n < 3:
        logger.warning(f"Need at least 3 samples for ordination, got {n}")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, f"Need >= 3 samples (got {n})", ha="center", va="center")
        return _save_or_show(fig, output_path)

    # Build distance matrix as 2D list
    D = [[0.0] * n for _ in range(n)]
    for i, si in enumerate(samples):
        for j, sj in enumerate(samples):
            D[i][j] = distance_matrix.get(si, {}).get(sj, 0.0)

    # Classical MDS / PCoA via double centering
    # 1. Compute squared distances
    D2 = [[d**2 for d in row] for row in D]

    # 2. Double centering: B = -0.5 * J * D^2 * J, where J = I - 1/n * 11'
    row_means = [sum(row) / n for row in D2]
    grand_mean = sum(row_means) / n

    B = [[0.0] * n for _ in range(n)]
    for i in range(n):
        col_mean_cache = [sum(D2[k][j] for k in range(n)) / n for j in range(n)]
        for j in range(n):
            B[i][j] = -0.5 * (D2[i][j] - row_means[i] - col_mean_cache[j] + grand_mean)

    # 3. Power iteration for top 2 eigenvectors
    coords = _power_iteration_2d(B, n)

    fig, ax = plt.subplots(figsize=(8, 8))

    if groups:
        unique_groups = sorted(set(groups.values()))
        cmap = plt.cm.Set1
        group_colors = {g: cmap(i / max(len(unique_groups), 1)) for i, g in enumerate(unique_groups)}

        for i, sample in enumerate(samples):
            group = groups.get(sample, "Unknown")
            color = group_colors.get(group, "gray")
            ax.scatter(coords[i][0], coords[i][1], c=[color], s=80, zorder=3)
            ax.annotate(sample, (coords[i][0], coords[i][1]), fontsize=7, xytext=(5, 5), textcoords="offset points")

        handles = [mpatches.Patch(color=c, label=g) for g, c in group_colors.items()]
        ax.legend(handles=handles, title="Group", fontsize=9)
    else:
        for i, sample in enumerate(samples):
            ax.scatter(coords[i][0], coords[i][1], c="steelblue", s=80, zorder=3)
            ax.annotate(sample, (coords[i][0], coords[i][1]), fontsize=7, xytext=(5, 5), textcoords="offset points")

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(title or f"Ordination ({method})")
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="gray", linewidth=0.5, linestyle="--")

    fig.tight_layout()
    return _save_or_show(fig, output_path)


def _power_iteration_2d(matrix: list[list[float]], n: int, max_iter: int = 300) -> list[list[float]]:
    """Extract top 2 eigenvectors using power iteration.

    Simple pure-Python implementation for PCoA coordinate extraction.
    """
    rng = random.Random(42)

    def _mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
        return [sum(mat[i][j] * vec[j] for j in range(n)) for i in range(n)]

    def _dot(a: list[float], b: list[float]) -> float:
        return sum(ai * bi for ai, bi in zip(a, b))

    def _norm(v: list[float]) -> float:
        return math.sqrt(max(_dot(v, v), 1e-30))

    def _normalize(v: list[float]) -> list[float]:
        n_val = _norm(v)
        return [vi / n_val for vi in v]

    coords: list[list[float]] = [[0.0, 0.0] for _ in range(n)]
    eigenvectors: list[list[float]] = []

    for comp in range(2):
        v = [rng.gauss(0, 1) for _ in range(n)]
        v = _normalize(v)

        for _ in range(max_iter):
            v_new = _mat_vec(matrix, v)

            # Deflate: remove projections onto previous eigenvectors
            for prev in eigenvectors:
                proj = _dot(v_new, prev)
                v_new = [vi - proj * pi for vi, pi in zip(v_new, prev)]

            new_norm = _norm(v_new)
            if new_norm < 1e-15:
                break
            v_new = _normalize(v_new)

            # Check convergence
            diff = sum((a - b) ** 2 for a, b in zip(v, v_new))
            v = v_new
            if diff < 1e-12:
                break

        eigenvalue = _dot(_mat_vec(matrix, v), v)
        scale = math.sqrt(max(eigenvalue, 0))

        for i in range(n):
            coords[i][comp] = v[i] * scale

        eigenvectors.append(v)

    return coords


def plot_alpha_diversity(
    samples: dict[str, dict[str, int]],
    metrics: list[str] | None = None,
    output_path: str | Path | None = None,
    groups: dict[str, str] | None = None,
    title: str = "Alpha Diversity",
) -> Any:
    """Plot alpha diversity metrics across samples.

    Args:
        samples: Dict mapping sample IDs to {OTU_id: count} dicts.
        metrics: List of metrics to compute ('observed', 'shannon', 'simpson', 'chao1').
            Default: all four.
        output_path: Path to save figure.
        groups: Optional grouping for samples.
        title: Plot title.

    Returns:
        Path to saved file or matplotlib figure.
    """
    _ensure_matplotlib()

    if metrics is None:
        metrics = ["observed", "shannon", "simpson", "chao1"]

    # Compute diversity metrics for each sample
    diversity_data: dict[str, dict[str, float]] = {}

    for sample_id, otu_counts in samples.items():
        counts = [c for c in otu_counts.values() if c > 0]
        total = sum(counts)

        sample_metrics: dict[str, float] = {}

        if "observed" in metrics:
            sample_metrics["observed"] = float(len(counts))

        if "shannon" in metrics:
            if total > 0:
                shannon = 0.0
                for c in counts:
                    p = c / total
                    if p > 0:
                        shannon -= p * math.log(p)
                sample_metrics["shannon"] = shannon
            else:
                sample_metrics["shannon"] = 0.0

        if "simpson" in metrics:
            if total > 1:
                simpson = 1.0 - sum(c * (c - 1) for c in counts) / (total * (total - 1))
                sample_metrics["simpson"] = simpson
            else:
                sample_metrics["simpson"] = 0.0

        if "chao1" in metrics:
            singletons = sum(1 for c in counts if c == 1)
            doubletons = sum(1 for c in counts if c == 2)
            observed = len(counts)
            if doubletons > 0:
                chao1 = observed + (singletons**2) / (2 * doubletons)
            elif singletons > 0:
                chao1 = observed + singletons * (singletons - 1) / 2
            else:
                chao1 = float(observed)
            sample_metrics["chao1"] = chao1

        diversity_data[sample_id] = sample_metrics

    n_metrics = len(metrics)
    fig, axes = plt.subplots(1, n_metrics, figsize=(4 * n_metrics, 5))
    if n_metrics == 1:
        axes = [axes]

    for ax, metric in zip(axes, metrics):
        if groups:
            group_values: dict[str, list[float]] = defaultdict(list)
            for sample_id, sample_metrics in diversity_data.items():
                group = groups.get(sample_id, "Unknown")
                group_values[group].append(sample_metrics.get(metric, 0.0))

            group_names = sorted(group_values.keys())
            positions = list(range(len(group_names)))
            data_to_plot = [group_values[g] for g in group_names]

            bp = ax.boxplot(data_to_plot, positions=positions, patch_artist=True)
            cmap = plt.cm.Set2
            for i, patch in enumerate(bp["boxes"]):
                patch.set_facecolor(cmap(i / max(len(group_names), 1)))

            ax.set_xticks(positions)
            ax.set_xticklabels(group_names, rotation=45, ha="right")
        else:
            values = [diversity_data[s].get(metric, 0.0) for s in sorted(diversity_data)]
            ax.bar(range(len(values)), values, color="steelblue", alpha=0.7)
            ax.set_xticks(range(len(values)))
            ax.set_xticklabels(sorted(diversity_data), rotation=45, ha="right", fontsize=7)

        ax.set_title(metric.capitalize())
        ax.set_ylabel("Value")

    fig.suptitle(title, fontsize=14)
    fig.tight_layout()
    return _save_or_show(fig, output_path)


def plot_heatmap(
    abundance_matrix: dict[str, dict[str, float]],
    output_path: str | Path | None = None,
    cluster: bool = True,
    title: str = "Abundance Heatmap",
    cmap: str = "YlOrRd",
    top_n: int = 30,
) -> Any:
    """Plot abundance heatmap of taxa across samples.

    Args:
        abundance_matrix: Dict mapping sample IDs to {taxon: abundance} dicts.
        output_path: Path to save figure.
        cluster: Whether to cluster rows and columns.
        title: Plot title.
        cmap: Matplotlib colormap name.
        top_n: Maximum number of taxa to show.

    Returns:
        Path to saved file or matplotlib figure.
    """
    _ensure_matplotlib()

    if not abundance_matrix:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        return _save_or_show(fig, output_path)

    # Collect all taxa and select top N
    taxa_totals: dict[str, float] = defaultdict(float)
    for sample_data in abundance_matrix.values():
        for taxon, abundance in sample_data.items():
            taxa_totals[taxon] += abundance

    sorted_taxa = sorted(taxa_totals.items(), key=lambda x: x[1], reverse=True)
    selected_taxa = [t[0] for t in sorted_taxa[:top_n]]

    samples = sorted(abundance_matrix.keys())

    # Build matrix
    matrix = []
    for taxon in selected_taxa:
        row = []
        for sample in samples:
            val = abundance_matrix[sample].get(taxon, 0.0)
            row.append(val)
        matrix.append(row)

    n_taxa = len(selected_taxa)
    n_samples = len(samples)

    if cluster and n_taxa > 2 and n_samples > 2:
        # Simple hierarchical clustering by sorting by similarity
        # Row clustering: sort taxa by their abundance profiles
        row_order = _hierarchical_sort(matrix)
        matrix = [matrix[i] for i in row_order]
        selected_taxa = [selected_taxa[i] for i in row_order]

        # Column clustering: sort samples by their profiles
        transposed = [[matrix[i][j] for i in range(len(matrix))] for j in range(n_samples)]
        col_order = _hierarchical_sort(transposed)
        samples = [samples[i] for i in col_order]
        matrix = [[row[j] for j in col_order] for row in matrix]

    fig, ax = plt.subplots(figsize=(max(8, n_samples * 0.4), max(6, n_taxa * 0.3)))

    # Plot heatmap manually
    for i in range(n_taxa):
        for j in range(n_samples):
            val = matrix[i][j]
            max_val = max(max(row) for row in matrix) if matrix else 1.0
            intensity = val / max_val if max_val > 0 else 0
            color = plt.cm.get_cmap(cmap)(intensity)
            ax.add_patch(plt.Rectangle((j, n_taxa - 1 - i), 1, 1, facecolor=color, edgecolor="white", linewidth=0.5))

    ax.set_xlim(0, n_samples)
    ax.set_ylim(0, n_taxa)
    ax.set_xticks([j + 0.5 for j in range(n_samples)])
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax.set_yticks([n_taxa - 1 - i + 0.5 for i in range(n_taxa)])
    ax.set_yticklabels(selected_taxa, fontsize=7)
    ax.set_title(title)

    # Colorbar
    max_val = max(max(row) for row in matrix) if matrix and any(matrix) else 1.0
    sm = ScalarMappable(cmap=cmap, norm=Normalize(vmin=0, vmax=max_val))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Abundance", shrink=0.6)

    fig.tight_layout()
    return _save_or_show(fig, output_path)


def _hierarchical_sort(vectors: list[list[float]]) -> list[int]:
    """Sort vectors by similarity using agglomerative approach.

    Simple greedy nearest-neighbor ordering.
    """
    n = len(vectors)
    if n <= 1:
        return list(range(n))

    def _distance(a: list[float], b: list[float]) -> float:
        return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))

    remaining = set(range(n))
    order = [0]
    remaining.remove(0)

    while remaining:
        last = order[-1]
        best_idx = min(remaining, key=lambda i: _distance(vectors[last], vectors[i]))
        order.append(best_idx)
        remaining.remove(best_idx)

    return order


__all__ = [
    "plot_krona_chart",
    "plot_stacked_bar",
    "plot_rarefaction_curves",
    "plot_ordination",
    "plot_alpha_diversity",
    "plot_heatmap",
]
