"""Interactive plot wrappers using Plotly with static fallback.

Provides interactive HTML-exportable versions of common plot types.
When Plotly is unavailable the functions fall back to static matplotlib
figures with a logged warning. An HTML ``output_path`` needs Plotly, so
without it the functions raise ``ImportError`` instead of writing a static
image under an HTML name; any other ``output_path`` is saved as a static image.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from metainformant.core.io import paths
from metainformant.core.utils import logging
from metainformant.visualization.config.conventions import save_figure_deterministic

logger = logging.get_logger(__name__)

try:
    import plotly.express as px
    import plotly.graph_objects as go

    HAS_PLOTLY = True
except ImportError:
    go = None
    px = None
    HAS_PLOTLY = False

_HTML_SUFFIXES = (".html", ".htm")


def _require_plotly_for_html(output_path: str | Path | None) -> None:
    """Refuse an HTML export that the static fallback cannot honour.

    matplotlib has no HTML writer, so without Plotly an ``.html`` path used to fail deep inside ``savefig`` with
    ``ValueError: Format 'html' is not supported``. Call this before any work is done: it raises a clear
    ``ImportError`` that names Plotly, and leaves no figure, file or directory behind.

    Raises:
        ImportError: If Plotly is not installed and ``output_path`` is an HTML file.
    """
    if HAS_PLOTLY or not output_path or Path(output_path).suffix.lower() not in _HTML_SUFFIXES:
        return
    raise ImportError(
        f"Plotly is required to write the interactive HTML file {str(output_path)!r}; install it with "
        "'uv pip install plotly', or pass a static image path (for example .png or .svg) to export the "
        "matplotlib fallback."
    )


# ---------------------------------------------------------------------------
# Interactive scatter
# ---------------------------------------------------------------------------


def _marker_sizes(sizes: pd.Series) -> np.ndarray:
    """Marker areas (points squared) for a size column: 20 for its smallest value up to 200 for its largest.

    A column that does not vary gets one middle size.
    """
    values = pd.to_numeric(sizes, errors="raise").to_numpy(dtype=float)
    low, high = np.nanmin(values), np.nanmax(values)
    if not high > low:
        return np.full(len(values), 60.0)
    return np.asarray(20.0 + (values - low) / (high - low) * 180.0)


def _scatter_by_category(
    ax: Axes,
    x: pd.Series,
    y: pd.Series,
    categories: pd.Series,
    *,
    legend_title: str,
    sizes: np.ndarray | None = None,
) -> None:
    """Scatter *x* against *y* with one coloured series per category and a legend.

    Mirrors the discrete colouring Plotly gives a non-numeric column. Categories
    keep their order of first appearance, and rows with a missing category are
    drawn grey as "missing" rather than dropped. ``sizes`` are per-row marker areas.
    """
    from ..config import palettes

    xs, ys = x.to_numpy(), y.to_numpy()
    values = categories.astype(object).to_numpy()
    present = categories.notna().to_numpy()
    labels = list(dict.fromkeys(values[present]))
    for label, color in zip(labels, palettes.categorical(len(labels))):
        mask = present & (values == label)
        ax.scatter(
            xs[mask], ys[mask], s=None if sizes is None else sizes[mask], color=color, alpha=0.7, label=str(label)
        )
    if not present.all():
        ax.scatter(
            xs[~present],
            ys[~present],
            s=None if sizes is None else sizes[~present],
            color="#9E9E9E",
            alpha=0.7,
            label="missing",
        )
    if ax.get_legend_handles_labels()[0]:  # an empty frame has nothing to label
        ax.legend(title=legend_title)


def interactive_scatter(
    data: pd.DataFrame,
    x: str,
    y: str,
    *,
    color: str | None = None,
    size: str | None = None,
    hover_data: List[str] | None = None,
    title: str = "Interactive Scatter",
    output_path: str | Path | None = None,
) -> Any:
    """Create an interactive scatter plot.

    Args:
        data: DataFrame with columns for axes.
        x: Column name for x-axis.
        y: Column name for y-axis.
        color: Column for color encoding. Without Plotly, a numeric column is
            drawn as a continuous colour scale and any other column as one
            coloured series per category with a legend.
        size: Column for size encoding (without Plotly, marker areas scale linearly from 20 to 200 with the column).
        hover_data: Extra columns shown on hover. A static image has no hover, so without Plotly they are not
            drawn and a warning says so.
        title: Figure title.
        output_path: Path to save. HTML with Plotly; without Plotly a static image
            (such as PNG), and an HTML path raises ``ImportError``.

    Returns:
        Plotly Figure or matplotlib Axes.

    Raises:
        ImportError: If ``output_path`` is an HTML file and Plotly is not installed.
    """
    _require_plotly_for_html(output_path)
    if HAS_PLOTLY:
        fig = px.scatter(
            data,
            x=x,
            y=y,
            color=color,
            size=size,
            hover_data=hover_data,
            title=title,
        )
        if output_path:
            out = Path(output_path)
            paths.ensure_directory(out.parent)
            if out.suffix == ".html":
                fig.write_html(str(out))
            else:
                fig.write_html(str(out.with_suffix(".html")))
            logger.info(f"Interactive scatter saved to {out}")
        return fig

    # Matplotlib fallback
    logger.warning("Plotly not available, falling back to static scatter")
    if hover_data:
        logger.warning(f"A static scatter has no hover text; hover_data {list(hover_data)} is not drawn")
    fig_m, ax = plt.subplots(figsize=(8, 6))
    c = data[color] if color and color in data.columns else None
    sizes = _marker_sizes(data[size]) if size is not None else None
    if c is not None and not pd.api.types.is_numeric_dtype(c):
        # matplotlib would read the strings as colour names; give each category its own colour instead
        _scatter_by_category(ax, data[x], data[y], c, legend_title=str(color), sizes=sizes)
    else:
        ax.scatter(data[x], data[y], c=c, s=sizes, alpha=0.7)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(title)
    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(fig_m, output_path, dpi=300, bbox_inches="tight")
    return ax


# ---------------------------------------------------------------------------
# Interactive heatmap
# ---------------------------------------------------------------------------


def _matplotlib_colormap(colorscale: str) -> str:
    """Name of the matplotlib colormap that stands in for a Plotly colorscale name (default "RdBu_r")."""
    registered = {name.lower(): name for name in plt.colormaps()}
    if colorscale.lower() in registered:
        return registered[colorscale.lower()]
    logger.warning(f"colorscale {colorscale!r} is not a matplotlib colormap; the static heatmap uses 'RdBu_r'")
    return "RdBu_r"


def interactive_heatmap(
    data: pd.DataFrame | np.ndarray,
    *,
    x_labels: List[str] | None = None,
    y_labels: List[str] | None = None,
    colorscale: str = "RdBu_r",
    title: str = "Interactive Heatmap",
    output_path: str | Path | None = None,
) -> Any:
    """Create an interactive heatmap.

    Args:
        data: 2-D array or DataFrame.
        x_labels: Column labels.
        y_labels: Row labels.
        colorscale: Plotly colorscale name. Without Plotly it is used as a matplotlib colormap name when one of that
            name exists (case-insensitively, so "Viridis" and "RdBu_r" work); otherwise a warning says so and
            "RdBu_r" is used.
        title: Figure title.
        output_path: Path to save. Without Plotly an HTML path raises ``ImportError``
            and any other path is saved as a static image.

    Returns:
        Plotly Figure or matplotlib Axes.

    Raises:
        ImportError: If ``output_path`` is an HTML file and Plotly is not installed.
    """
    _require_plotly_for_html(output_path)
    if isinstance(data, pd.DataFrame):
        z = data.values
        x_labels = x_labels or list(data.columns)
        y_labels = y_labels or list(data.index)
    else:
        z = np.asarray(data)

    if HAS_PLOTLY:
        fig = go.Figure(data=go.Heatmap(z=z, x=x_labels, y=y_labels, colorscale=colorscale))
        fig.update_layout(title=title)
        if output_path:
            out = Path(output_path)
            paths.ensure_directory(out.parent)
            fig.write_html(str(out.with_suffix(".html")))
            logger.info(f"Interactive heatmap saved to {out}")
        return fig

    logger.warning("Plotly not available, falling back to static heatmap")
    cmap = _matplotlib_colormap(colorscale)
    fig_m, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(z, cmap=cmap, aspect="auto")
    plt.colorbar(im, ax=ax)
    if x_labels:
        ax.set_xticks(range(len(x_labels)))
        ax.set_xticklabels(x_labels, rotation=45, ha="right")
    if y_labels:
        ax.set_yticks(range(len(y_labels)))
        ax.set_yticklabels(y_labels)
    ax.set_title(title)
    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(fig_m, output_path, dpi=300, bbox_inches="tight")
    return ax


# ---------------------------------------------------------------------------
# Interactive volcano
# ---------------------------------------------------------------------------


def interactive_volcano(
    data: pd.DataFrame,
    *,
    log2fc_col: str = "log2FoldChange",
    pvalue_col: str = "pvalue",
    gene_col: str | None = "gene",
    fc_threshold: float = 1.0,
    p_threshold: float = 0.05,
    title: str = "Interactive Volcano Plot",
    output_path: str | Path | None = None,
) -> Any:
    """Create an interactive volcano plot.

    Args:
        data: DataFrame with fold-change and p-value columns.
        log2fc_col: Column name for log2 fold-change.
        pvalue_col: Column name for p-values.
        gene_col: Column name for gene labels (shown on hover).
        fc_threshold: Absolute fold-change threshold for significance.
        p_threshold: P-value threshold for significance.
        title: Figure title.
        output_path: Path to save. Without Plotly an HTML path raises ``ImportError``
            and any other path is saved as a static image.

    Returns:
        Plotly Figure or matplotlib Axes.

    Raises:
        ImportError: If ``output_path`` is an HTML file and Plotly is not installed.
    """
    _require_plotly_for_html(output_path)
    df = data.copy()
    df["neg_log10_p"] = -np.log10(np.clip(df[pvalue_col].values, 1e-300, 1))
    df["significance"] = "Not Significant"
    sig_mask = (np.abs(df[log2fc_col]) > fc_threshold) & (df[pvalue_col] < p_threshold)
    df.loc[sig_mask, "significance"] = "Significant"
    up_mask = sig_mask & (df[log2fc_col] > 0)
    down_mask = sig_mask & (df[log2fc_col] < 0)
    df.loc[up_mask, "significance"] = "Up-regulated"
    df.loc[down_mask, "significance"] = "Down-regulated"

    color_map = {
        "Not Significant": "#9E9E9E",
        "Significant": "#FF9800",
        "Up-regulated": "#D32F2F",
        "Down-regulated": "#2166AC",
    }

    if HAS_PLOTLY:
        hover = [gene_col] if gene_col and gene_col in df.columns else None
        fig = px.scatter(
            df,
            x=log2fc_col,
            y="neg_log10_p",
            color="significance",
            color_discrete_map=color_map,
            hover_data=hover,
            title=title,
        )
        fig.add_hline(y=-np.log10(p_threshold), line_dash="dash", line_color="grey")
        fig.add_vline(x=fc_threshold, line_dash="dash", line_color="grey")
        fig.add_vline(x=-fc_threshold, line_dash="dash", line_color="grey")
        if output_path:
            out = Path(output_path)
            paths.ensure_directory(out.parent)
            fig.write_html(str(out.with_suffix(".html")))
            logger.info(f"Interactive volcano saved to {out}")
        return fig

    logger.warning("Plotly not available, falling back to static volcano")
    fig_m, ax = plt.subplots(figsize=(8, 6))
    for sig, color in color_map.items():
        mask = df["significance"] == sig
        if (
            not mask.any()
        ):  # a category with no genes (such as "Significant", which up/down always replace) has no entry
            continue
        ax.scatter(df.loc[mask, log2fc_col], df.loc[mask, "neg_log10_p"], c=color, label=sig, s=10, alpha=0.6)
    ax.axhline(-np.log10(p_threshold), ls="--", color="grey", alpha=0.7)
    ax.axvline(fc_threshold, ls="--", color="grey", alpha=0.7)
    ax.axvline(-fc_threshold, ls="--", color="grey", alpha=0.7)
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(title)
    ax.legend(markerscale=2)
    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(fig_m, output_path, dpi=300, bbox_inches="tight")
    return ax


# ---------------------------------------------------------------------------
# Interactive Manhattan
# ---------------------------------------------------------------------------


def interactive_manhattan(
    data: pd.DataFrame,
    *,
    chrom_col: str = "chromosome",
    pos_col: str = "position",
    pvalue_col: str = "pvalue",
    gene_col: str | None = None,
    significance_threshold: float = 5e-8,
    suggestive_threshold: float = 1e-5,
    title: str = "Interactive Manhattan Plot",
    output_path: str | Path | None = None,
) -> Any:
    """Create an interactive Manhattan plot.

    Args:
        data: DataFrame with chromosome, position, p-value columns.
        chrom_col: Column for chromosome.
        pos_col: Column for genomic position.
        pvalue_col: Column for p-values.
        gene_col: Column for gene annotations (hover).
        significance_threshold: Genome-wide significance line.
        suggestive_threshold: Suggestive significance line.
        title: Figure title.
        output_path: Path to save. Without Plotly an HTML path raises ``ImportError``
            and any other path is saved as a static image.

    Returns:
        Plotly Figure or matplotlib Axes.

    Raises:
        ImportError: If ``output_path`` is an HTML file and Plotly is not installed.
    """
    _require_plotly_for_html(output_path)
    df = data.copy()
    df["neg_log10_p"] = -np.log10(np.clip(df[pvalue_col].values, 1e-300, 1))

    # Compute cumulative position
    chroms = df[chrom_col].unique()
    chrom_order = sorted(
        chroms,
        key=lambda c: (
            int(str(c).replace("chr", ""))
            if str(c).replace("chr", "").isdigit()
            else 100 + ord(str(c).replace("chr", "")[0])
        ),
    )
    chrom_offsets: Dict[Any, int] = {}
    cumulative = 0
    for ch in chrom_order:
        chrom_offsets[ch] = cumulative
        cumulative += int(df.loc[df[chrom_col] == ch, pos_col].max())

    df["cumulative_pos"] = df.apply(lambda r: r[pos_col] + chrom_offsets.get(r[chrom_col], 0), axis=1)

    if HAS_PLOTLY:
        hover = [chrom_col, pos_col, pvalue_col]
        if gene_col and gene_col in df.columns:
            hover.append(gene_col)
        fig = px.scatter(
            df,
            x="cumulative_pos",
            y="neg_log10_p",
            color=chrom_col,
            hover_data=hover,
            title=title,
        )
        fig.add_hline(
            y=-np.log10(significance_threshold), line_dash="dash", line_color="red", annotation_text="Genome-wide"
        )
        fig.add_hline(
            y=-np.log10(suggestive_threshold), line_dash="dot", line_color="blue", annotation_text="Suggestive"
        )
        fig.update_layout(xaxis_title="Genomic Position", yaxis_title="-log10(p)")
        if output_path:
            out = Path(output_path)
            paths.ensure_directory(out.parent)
            fig.write_html(str(out.with_suffix(".html")))
            logger.info(f"Interactive Manhattan saved to {out}")
        return fig

    logger.warning("Plotly not available, falling back to static Manhattan")
    from ..config import palettes

    fig_m, ax = plt.subplots(figsize=(14, 6))
    colors = palettes.alternating_pair(len(chrom_order))
    for i, ch in enumerate(chrom_order):
        mask = df[chrom_col] == ch
        ax.scatter(
            df.loc[mask, "cumulative_pos"], df.loc[mask, "neg_log10_p"], c=colors[i % len(colors)], s=5, alpha=0.6
        )
    ax.axhline(-np.log10(significance_threshold), ls="--", color="red", alpha=0.7)
    ax.axhline(-np.log10(suggestive_threshold), ls=":", color="blue", alpha=0.7)
    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("-log10(p)")
    ax.set_title(title)
    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(fig_m, output_path, dpi=300, bbox_inches="tight")
    return ax
