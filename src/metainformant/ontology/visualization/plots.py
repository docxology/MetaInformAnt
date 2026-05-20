"""Enrichment and ontology network visualization.

Provides matplotlib-based plots for GO/pathway enrichment results:

- :func:`enrichment_dotplot` — classic enrichment dot plot
  (gene ratio × GO term, dot size = overlap count, colour = -log10 adj_p)
- :func:`pathway_network_plot` — pathway similarity network
  (nodes = pathways, edges = Jaccard similarity)

All functions return a ``matplotlib.figure.Figure`` and optionally save
to disk.  No Seaborn dependency; pure matplotlib.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Optional

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Colour helpers
# ---------------------------------------------------------------------------


def _cmap_red_blue(t: float) -> tuple[float, float, float]:
    """Map t in [0,1] to a red-white-blue diverging colour."""
    t = max(0.0, min(1.0, t))
    if t < 0.5:
        r = 1.0
        g = 2 * t
        b = 2 * t
    else:
        r = 2 * (1 - t)
        g = 2 * (1 - t)
        b = 1.0
    return (r, g, b)


def _viridis_approx(t: float) -> tuple[float, float, float]:
    """Approximate viridis colormap (purple -> yellow)."""
    t = max(0.0, min(1.0, t))
    r = 0.267 + 0.004 * t + 0.658 * t**2 + 0.073 * t**3
    g = 0.005 + 1.096 * t - 0.900 * t**2 + 0.498 * t**3
    b = 0.329 + 1.648 * t - 3.294 * t**2 + 1.780 * t**3
    return (max(0, min(1, r)), max(0, min(1, g)), max(0, min(1, b)))


# ---------------------------------------------------------------------------
# Enrichment dot plot
# ---------------------------------------------------------------------------


def enrichment_dotplot(
    results: list[dict[str, Any]],
    output_path: str | Path | None = None,
    *,
    top_n: int = 20,
    title: str = "GO Enrichment",
    alpha_threshold: float = 0.05,
    figsize: tuple[float, float] | None = None,
    ax: Optional[Any] = None,
    style: Optional[Any] = None,
    **kwargs: Any,
) -> Any:
    """Plot a GO/pathway enrichment dot plot.

    Circles represent enriched terms; X-axis = gene ratio
    (overlap / query size), Y-axis = term label, dot area proportional to
    overlap count, dot colour = -log10(adjusted p-value).

    Args:
        results: List of enrichment result dicts from
            :func:`~metainformant.ontology.pathway_enrichment.enrichment.over_representation_analysis`
            or similar.  Must contain at minimum ``term_id`` / ``term_name``,
            ``adjusted_p``, ``n_overlap``, ``n_genes``.
        output_path: Optional file path to save PNG.
        top_n: Maximum number of terms to show (most significant first).
        title: Plot title.
        alpha_threshold: FDR threshold for dashed-line annotation.
        figsize: Figure size ``(width, height)`` in inches.

    Returns:
        ``matplotlib.figure.Figure``.

    Examples:
        >>> fig = enrichment_dotplot(ora_results, "enrichment_dot.png", top_n=15)
    """
    try:
        import matplotlib.pyplot as plt

        # Apply accessible, high-contrast presentation style
        plt.rcParams.update(
            {
                "font.size": 14,
                "axes.titlesize": 16,
                "axes.labelsize": 14,
                "xtick.labelsize": 12,
                "ytick.labelsize": 12,
                "legend.fontsize": 12,
                "legend.framealpha": 0.85,
                "legend.edgecolor": "black",
                "figure.dpi": 300,
                "axes.grid": True,
                "axes.axisbelow": True,
                "grid.alpha": 0.4,
                "grid.color": "#cccccc",
            }
        )
    except ImportError:
        logger.warning("matplotlib not available; skipping enrichment_dotplot")
        return None

    if not results:
        logger.warning("enrichment_dotplot: no results to plot")
        return None

    # Sort by adjusted_p ascending, take top N
    sorted_res = sorted(results, key=lambda r: r.get("adjusted_p", r.get("p_value", 1.0)))
    top = sorted_res[:top_n]

    if not top:
        return None

    # Prepare values
    labels = [r.get("term_name", r.get("term_id", f"term_{i}")) for i, r in enumerate(top)]
    adj_ps = [max(r.get("adjusted_p", r.get("p_value", 1.0)), 1e-300) for r in top]
    neg_log_ps = [-math.log10(p) for p in adj_ps]
    overlaps = [r.get("n_overlap", r.get("observed", 1)) for r in top]
    [r.get("n_genes", r.get("term_size", 1)) for r in top]
    query_sizes = [r.get("query_size", max(overlaps) or 1) for r in top]
    gene_ratios = [ov / qs if qs > 0 else 0.0 for ov, qs in zip(overlaps, query_sizes)]

    # Colour by -log10(adj_p)
    max_nlp = max(neg_log_ps) if neg_log_ps else 1.0
    max_nlp = max(max_nlp, 1e-6)  # guard: prevents div-by-zero when all p=1.0
    colours = [_viridis_approx(nlp / max_nlp) for nlp in neg_log_ps]

    # Dot size proportional to overlap (scaled to a nice range)
    max_ov = max(overlaps) if overlaps else 1
    dot_sizes = [40 + 200 * (ov / max_ov) for ov in overlaps]

    fsize = figsize or kwargs.get("figsize", (8, max(4, 0.4 * len(top) + 2)))
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=fsize)

    y_pos = list(range(len(top)))

    ax.scatter(
        gene_ratios,
        y_pos,
        s=dot_sizes,
        c=colours,
        alpha=0.85,
        zorder=3,
    )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Gene ratio (overlap / query)", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.invert_yaxis()

    # Colour bar legend
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0, max_nlp))
    sm.set_array([])
    try:
        cbar = fig.colorbar(sm, ax=ax, pad=0.01)
        cbar.set_label("-log10(adj p-value)", fontsize=9)
    except Exception:
        pass

    # Size legend
    mid_ov = max(1, max(overlaps) // 2)
    for ov_val, label in [(1, "1 gene"), (mid_ov, str(mid_ov)), (max(overlaps), str(max(overlaps)))]:
        ax.scatter([], [], s=40 + 200 * (ov_val / max_ov), c="grey", alpha=0.7, label=label)
    ax.legend(title="Overlap genes", bbox_to_anchor=(1.15, 0.3), loc="center left", fontsize=8)

    plt.tight_layout()

    if output_path and fig is not None:
        fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
        logger.info("Enrichment dotplot saved to %s", output_path)

    return ax if fig is None else fig


# ---------------------------------------------------------------------------
# Pathway network plot
# ---------------------------------------------------------------------------


def pathway_network_plot(
    network: dict[str, Any],
    output_path: str | Path | None = None,
    *,
    top_n: int = 30,
    title: str = "Pathway Similarity Network",
    figsize: tuple[float, float] = (9, 9),
    ax: Optional[Any] = None,
    style: Optional[Any] = None,
    **kwargs: Any,
) -> Any:
    """Plot a pathway similarity network.

    Nodes = pathway terms (size proportional to gene-set size, colour =
    -log10 adj_p), edges = Jaccard similarity weight.  Uses a simple
    spring-layout (no networkx required).

    Args:
        network: Dict from
            :func:`~metainformant.ontology.pathway_enrichment.enrichment.pathway_network`
            with keys ``nodes``, ``edges``, ``clusters``.
        output_path: Optional output file path.
        top_n: Maximum number of nodes to render.
        title: Plot title.
        figsize: Figure dimensions in inches.

    Returns:
        ``matplotlib.figure.Figure``.

    Examples:
        >>> fig = pathway_network_plot(net, "pathway_net.png")
    """
    try:
        import matplotlib.pyplot as plt

        # Apply accessible, high-contrast presentation style
        plt.rcParams.update(
            {
                "font.size": 14,
                "axes.titlesize": 16,
                "axes.labelsize": 14,
                "xtick.labelsize": 12,
                "ytick.labelsize": 12,
                "legend.fontsize": 12,
                "legend.framealpha": 0.85,
                "legend.edgecolor": "black",
                "figure.dpi": 300,
                "axes.grid": True,
                "axes.axisbelow": True,
                "grid.alpha": 0.4,
                "grid.color": "#cccccc",
            }
        )
    except ImportError:
        logger.warning("matplotlib not available; skipping pathway_network_plot")
        return None

    nodes = network.get("nodes", [])[:top_n]
    edges = network.get("edges", [])

    if not nodes:
        logger.warning("pathway_network_plot: no nodes to plot")
        return None

    node_ids = {n["term_id"] for n in nodes}

    # Filter edges to visible nodes only
    visible_edges = [e for e in edges if e["source"] in node_ids and e["target"] in node_ids]

    # Simple spring layout using force-directed positions (no scipy/networkx)
    import math as _math
    import random

    random.seed(42)
    positions: dict[str, list[float]] = {n["term_id"]: [random.uniform(-1, 1), random.uniform(-1, 1)] for n in nodes}

    def _layout(pos: dict, edges_: list, iterations: int = 80) -> dict:
        ids = list(pos.keys())
        adj: dict[str, list[str]] = {i: [] for i in ids}
        for e in edges_:
            adj[e["source"]].append(e["target"])
            adj[e["target"]].append(e["source"])
        k = 1.0 / max(1, _math.sqrt(len(ids)))
        t = 0.1  # temperature
        for _ in range(iterations):
            disp: dict[str, list[float]] = {i: [0.0, 0.0] for i in ids}
            for i in range(len(ids)):
                for j in range(i + 1, len(ids)):
                    a, b = ids[i], ids[j]
                    dx = pos[a][0] - pos[b][0]
                    dy = pos[a][1] - pos[b][1]
                    dist = max(0.001, _math.hypot(dx, dy))
                    rep = k * k / dist
                    disp[a][0] += dx / dist * rep
                    disp[a][1] += dy / dist * rep
                    disp[b][0] -= dx / dist * rep
                    disp[b][1] -= dy / dist * rep
            for e in edges_:
                a, b = e["source"], e["target"]
                dx = pos[a][0] - pos[b][0]
                dy = pos[a][1] - pos[b][1]
                dist = max(0.001, _math.hypot(dx, dy))
                attr = dist * dist / k
                disp[a][0] -= dx / dist * attr
                disp[a][1] -= dy / dist * attr
                disp[b][0] += dx / dist * attr
                disp[b][1] += dy / dist * attr
            for i in ids:
                dm = max(0.001, _math.hypot(*disp[i]))
                pos[i][0] += disp[i][0] / dm * min(dm, t)
                pos[i][1] += disp[i][1] / dm * min(dm, t)
            t *= 0.95
        return pos

    positions = _layout(positions, visible_edges)

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get("figsize", figsize))

    # Draw edges
    for e in visible_edges:
        x0, y0 = positions[e["source"]]
        x1, y1 = positions[e["target"]]
        ax.plot([x0, x1], [y0, y1], color="#cccccc", linewidth=0.5 + 2 * e.get("jaccard", 0), alpha=0.6, zorder=1)

    # Draw nodes
    max_n_genes = max((n.get("n_genes", 1) for n in nodes), default=1)
    max_nlp = max((-_math.log10(max(n.get("adjusted_p", 1e-3), 1e-300)) for n in nodes), default=1.0)
    max_nlp = max(max_nlp, 1e-6)  # guard: prevents div-by-zero when all p=1.0

    for node in nodes:
        x, y = positions[node["term_id"]]
        n_genes = node.get("n_genes", 1)
        adj_p = node.get("adjusted_p", node.get("p_value", 1.0))
        neg_lp = -_math.log10(max(adj_p, 1e-300))
        colour = _viridis_approx(neg_lp / max_nlp if max_nlp > 0 else 0)
        size = 30 + 150 * (n_genes / max_n_genes)
        ax.scatter(x, y, s=size, color=colour, alpha=0.85, edgecolors="#444444", linewidths=0.4, zorder=3)
        label = node["term_id"]
        ax.annotate(label, (x, y), fontsize=5.5, ha="center", va="bottom", xytext=(0, 4), textcoords="offset points")

    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.axis("off")
    plt.tight_layout()

    if output_path and fig is not None:
        fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
        logger.info("Pathway network plot saved to %s", output_path)

    return ax if fig is None else fig


# ---------------------------------------------------------------------------
# GSEA enrichment score trace plot
# ---------------------------------------------------------------------------


def enrichment_score_plot(
    running_es: list[float],
    hit_indices: list[int],
    *,
    term_name: str = "Gene Set",
    es: float | None = None,
    output_path: str | Path | None = None,
    figsize: tuple[float, float] = (8, 4),
    axes: Optional[Any] = None,
    style: Optional[Any] = None,
    **kwargs: Any,
) -> Any:
    """Plot GSEA running enrichment score trace.

    Args:
        running_es: Running enrichment score list (from ``compute_enrichment_score``).
        hit_indices: Positions of gene-set members in the ranked list.
        term_name: Label for the title.
        es: Observed enrichment score to label on plot.
        output_path: Optional output file path.
        figsize: Figure dimensions.

    Returns:
        ``matplotlib.figure.Figure``.

    Examples:
        >>> fig = enrichment_score_plot(res["running_es"], res["hit_indices"])
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    fig = None
    if axes is None:
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=kwargs.get("figsize", figsize), gridspec_kw={"height_ratios": [3, 1]}, sharex=True
        )
    else:
        ax1, ax2 = axes

    # Running ES trace
    x = list(range(len(running_es)))
    ax1.plot(x, running_es, color="#2166ac", linewidth=1.2)
    ax1.axhline(0, color="#888888", linewidth=0.8, linestyle="--")
    if es is not None:
        ax1.axhline(es, color="#d6604d", linewidth=0.8, linestyle=":", label=f"ES = {es:.3f}")
        ax1.legend(fontsize=8)
    ax1.set_ylabel("Running ES", fontsize=9)
    ax1.set_title(f"GSEA: {term_name}", fontsize=10, fontweight="bold")

    # Hit rug
    ax2.set_ylim(0, 1)
    ax2.set_yticks([])
    for idx in hit_indices:
        ax2.axvline(idx, color="#d6604d", linewidth=0.5, alpha=0.7)
    ax2.set_xlabel("Rank in gene list", fontsize=9)
    ax2.set_ylabel("Hits", fontsize=9)

    plt.tight_layout()

    if output_path and fig is not None:
        fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
        logger.info("ES trace plot saved to %s", output_path)

    return axes if fig is None else fig
