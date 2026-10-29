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
import random
from collections import defaultdict
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger
from metainformant.metagenomics.diversity.metrics import ordination

logger = get_logger(__name__)

_ENV_PREFIX = "META_"

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    logger.warning("matplotlib not available; visualization functions will return data only")


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
            to abundance values, as counts or fractions (they are normalised to sum to 1).
            Example: {"Bacteria;Proteobacteria;Gammaproteobacteria": 0.35}
        output_path: Path to save figure. If None, returns figure object.
        title: Plot title.
        max_depth: Maximum taxonomic depth to display.
        min_fraction: Minimum fraction of the total abundance a taxon needs to be drawn.

    Returns:
        Path to saved file or matplotlib figure.

    Raises:
        ValueError: If ``taxonomy_data`` is not empty but its abundances do not sum to a
            positive number.
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
    if not total > 0:
        raise ValueError(f"taxonomy_data abundances must have a positive sum, got {total}")

    normalized = {k: v / total for k, v in taxonomy_data.items() if v / total >= min_fraction}

    # Build hierarchical structure
    tree: dict[str, Any] = {}
    for lineage, abundance in normalized.items():
        levels = [level.strip() for level in lineage.split(";") if level.strip()][:max_depth]
        node = tree
        for level in levels:
            if level not in node:
                node[level] = {"_abundance": 0.0, "_children": {}}
            node[level]["_abundance"] += abundance
            node = node[level]["_children"]

    # Generate sunburst rings
    fig, ax = plt.subplots(figsize=(12, 12), subplot_kw={"projection": "polar"})

    cmap = plt.get_cmap("Set3")
    color_idx = 0

    def _draw_ring(
        node: dict[str, Any],
        depth: int,
        start_angle: float,
        end_angle: float,
        node_abundance: float,
    ) -> None:
        nonlocal color_idx
        if depth >= max_depth or not node or node_abundance <= 0:
            return

        inner_r = 0.3 + depth * 0.15
        outer_r = inner_r + 0.14

        items = sorted(node.items(), key=lambda x: x[1]["_abundance"], reverse=True)
        current_angle = start_angle

        for name, data in items:
            # ``_abundance`` is already a fraction of the grand total (normalised above), so it is
            # compared with ``min_fraction`` as it is. A child takes its share of the parent's arc.
            frac = data["_abundance"]
            if frac < min_fraction:
                continue

            span = (end_angle - start_angle) * frac / node_abundance
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

            _draw_ring(data["_children"], depth + 1, current_angle, current_angle + span, frac)
            current_angle += span

    _draw_ring(tree, 0, 0, 2 * math.pi, sum(entry["_abundance"] for entry in tree.values()))

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

    cmap = plt.get_cmap("tab20") if len(top_taxa) > 10 else plt.get_cmap("Set3")
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
    cmap = plt.get_cmap("tab10")

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


def _distance_matrix_from_dict(
    distance_matrix: dict[str, dict[str, float]],
) -> tuple[list[str], list[list[float]]]:
    """Convert a nested ``{sample_i: {sample_j: dist}}`` mapping into sorted sample IDs and a dense matrix.

    A distance may be given in either triangle (``d[a][b]`` or ``d[b][a]``), and a sample's
    distance to itself defaults to 0. Any other missing distance is an error: it is never
    replaced by 0, which would silently place the two samples on top of each other.

    Args:
        distance_matrix: Nested dict of pairwise distances.

    Returns:
        ``(samples, matrix)`` with samples sorted and ``matrix[i][j]`` their distance.

    Raises:
        ValueError: If a distance between two different samples is missing.
    """
    samples = sorted(set(distance_matrix) | {other for row in distance_matrix.values() for other in row})
    matrix: list[list[float]] = []
    for sample_i in samples:
        row: list[float] = []
        for sample_j in samples:
            if sample_j in distance_matrix.get(sample_i, {}):
                row.append(float(distance_matrix[sample_i][sample_j]))
            elif sample_i in distance_matrix.get(sample_j, {}):
                row.append(float(distance_matrix[sample_j][sample_i]))
            elif sample_i == sample_j:
                row.append(0.0)
            else:
                raise ValueError(f"distance_matrix has no distance between {sample_i!r} and {sample_j!r}")
        matrix.append(row)
    return samples, matrix


def plot_ordination(
    distance_matrix: dict[str, dict[str, float]],
    method: str = "PCoA",
    output_path: str | Path | None = None,
    title: str | None = None,
    groups: dict[str, str] | None = None,
    seed: int = 42,
) -> Any:
    """Plot ordination (PCoA or NMDS) from a distance matrix.

    Performs dimensionality reduction on a pairwise distance matrix and
    plots samples in 2D space. Both methods come from
    ``metainformant.metagenomics.diversity.metrics.ordination``.

    PCoA axes are labelled with their share of the positive-eigenvalue variance. A
    distance matrix that is not Euclidean (Bray-Curtis, Jaccard, ...) has negative
    eigenvalues, which cannot be drawn; the plot says how much of the eigenvalue mass
    they hold, and the shares on the axes are relative to the positive eigenvalues only.
    NMDS axes are arbitrary, so they carry no variance; the plot shows Kruskal's stress-1
    instead.

    Args:
        distance_matrix: Nested dict of pairwise distances {sample_i: {sample_j: dist}}.
            Either triangle may be given; a distance that is missing in both is an error.
        method: Ordination method, case-insensitive: 'PCoA' (classical scaling by double
            centering and eigendecomposition, the default) or 'NMDS' (Kruskal's non-metric
            scaling by isotonic-regression SMACOF, as in
            ``metainformant.metagenomics.diversity.metrics.ordination``).
        output_path: Path to save figure.
        title: Plot title.
        groups: Optional dict mapping sample IDs to group labels for coloring.
        seed: Seed for the NMDS random starts (ignored by PCoA).

    Returns:
        Path to saved file or matplotlib figure.

    Raises:
        ValueError: If ``method`` is neither 'PCoA' nor 'NMDS', or the distance matrix
            has a missing, negative, non-finite or asymmetric distance.
    """
    _ensure_matplotlib()

    method_key = method.strip().lower() if isinstance(method, str) else ""
    if method_key not in ("pcoa", "nmds"):
        raise ValueError(f"Invalid ordination method {method!r}. Must be 'PCoA' or 'NMDS'.")
    method_label = "PCoA" if method_key == "pcoa" else "NMDS"

    n = len(set(distance_matrix) | {other for row in distance_matrix.values() for other in row})

    if n < 3:
        logger.warning(f"Need at least 3 samples for ordination, got {n}")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, f"Need >= 3 samples (got {n})", ha="center", va="center")
        return _save_or_show(fig, output_path)

    samples, D = _distance_matrix_from_dict(distance_matrix)
    result = ordination(D, method=method_key, n_components=2, seed=seed)
    coords = result["coordinates"]

    if method_key == "pcoa":
        shares = result["variance_explained"]
        x_label, y_label = f"PC1 ({shares[0]:.1%})", f"PC2 ({shares[1]:.1%})"
        note = None
        if result["n_negative_eigenvalues"]:
            note = (
                f"{result['n_negative_eigenvalues']} negative eigenvalues "
                f"({result['negative_eigenvalue_fraction']:.0%} of |eigenvalue| mass)"
            )
    else:
        x_label, y_label = "NMDS1", "NMDS2"
        note = f"stress = {result['stress']:.3f}"

    fig, ax = plt.subplots(figsize=(8, 8))

    if groups:
        unique_groups = sorted(set(groups.values()))
        cmap = plt.get_cmap("Set1")
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

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title or f"Ordination ({method_label})")
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="gray", linewidth=0.5, linestyle="--")
    if note:
        ax.text(0.99, 0.01, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="dimgray")

    fig.tight_layout()
    return _save_or_show(fig, output_path)


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
            cmap = plt.get_cmap("Set2")
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
        cluster: Whether to reorder rows and columns so that similar profiles are adjacent
            (a greedy nearest-neighbour ordering, not a hierarchical clustering; no dendrogram is drawn).
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
        # Greedy nearest-neighbour ordering (not hierarchical clustering)
        # Row ordering: sort taxa by their abundance profiles
        row_order = _hierarchical_sort(matrix)
        matrix = [matrix[i] for i in row_order]
        selected_taxa = [selected_taxa[i] for i in row_order]

        # Column ordering: sort samples by their profiles
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
            color = plt.get_cmap(cmap)(intensity)
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
    """Order vectors so that similar ones are adjacent (greedy nearest-neighbour chain).

    Starts from the first vector and repeatedly appends the nearest unused one by Euclidean
    distance. This is an ordering heuristic, not hierarchical (agglomerative) clustering.
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
