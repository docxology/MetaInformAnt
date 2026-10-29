"""Specialized visualization functions for advanced plotting types.

This module provides specialized visualization capabilities including Venn diagrams,
Sankey diagrams, chord diagrams, alluvial plots, and other advanced chart types
commonly used in bioinformatics and systems biology.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple, cast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Circle, PathPatch, Rectangle
from matplotlib.path import Path as MplPath

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging
from metainformant.visualization.config import palettes
from metainformant.visualization.config.conventions import save_figure_deterministic
from metainformant.visualization.config.figures import root_figure

logger = logging.get_logger(__name__)

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    sns = None

try:
    import plotly.express as px
    import plotly.graph_objects as go

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    go = None
    px = None


def plot_venn_diagram(
    sets: Dict[str, set],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 8),
    **kwargs: Any,
) -> Axes:
    """Create a Venn diagram from set data.

    Each region is labelled with the number of elements that belong to exactly the sets whose circles overlap
    there. matplotlib-venn draws the diagram when it is installed; otherwise a built-in drawing with unscaled
    circles is used, which has the same labels but ignores ``kwargs``.

    Args:
        sets: Dictionary mapping set names to set objects (two or three sets)
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for matplotlib-venn's venn2/venn3 (for example ``set_colors``)

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If ``sets`` does not hold exactly two or three sets (whether or not matplotlib-venn is installed)
    """
    validation.validate_type(sets, dict, "sets")

    if len(sets) not in (2, 3):
        raise ValueError(f"Venn diagrams support two or three sets, got {len(sets)}")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    try:
        import matplotlib_venn as venn
    except ImportError:
        # Fallback implementation for 2-3 sets
        if kwargs:
            logger.warning(f"matplotlib-venn is not installed; the built-in Venn drawing ignores {sorted(kwargs)}")
        if len(sets) == 2:
            _plot_venn2(ax, sets)
        else:
            _plot_venn3(ax, sets)
    else:
        # Use matplotlib-venn
        set_names = list(sets.keys())
        set_data = [sets[name] for name in set_names]

        if len(set_data) == 2:
            venn.venn2(set_data, set_labels=set_names, ax=ax, **kwargs)
        else:
            venn.venn3(set_data, set_labels=set_names, ax=ax, **kwargs)

    ax.set_title("Venn Diagram")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Venn diagram saved to {output_path}")

    return ax


def _plot_venn2(ax: Axes, sets: Dict[str, set]) -> None:
    """Built-in Venn diagram for 2 sets: exclusive counts in the left, middle and right regions."""
    set_names = list(sets.keys())
    set1, set2 = [sets[name] for name in set_names]

    # Calculate set sizes
    only1 = len(set1 - set2)
    only2 = len(set2 - set1)
    both = len(set1 & set2)

    # Draw circles: they overlap between x = -0.5 and x = 0.5
    circle1 = Circle((-1, 0), 1.5, facecolor="red", alpha=0.5)
    circle2 = Circle((1, 0), 1.5, facecolor="blue", alpha=0.5)

    ax.add_patch(circle1)
    ax.add_patch(circle2)

    # Add labels: the intersection count sits in the overlap, not in the first circle
    ax.text(0, 0, f"{both}\n({set_names[0]}∩{set_names[1]})", ha="center", va="center")
    ax.text(-1.7, 0, f"{only1}\n({set_names[0]} only)", ha="center", va="center")
    ax.text(1.7, 0, f"{only2}\n({set_names[1]} only)", ha="center", va="center")

    ax.set_xlim(-3, 3)
    ax.set_ylim(-2, 2)
    ax.axis("off")


def _plot_venn3(ax: Axes, sets: Dict[str, set]) -> None:
    """Built-in Venn diagram for 3 sets: the count of each of the seven regions, in the region itself.

    The numbers are exclusive counts (an element is counted in the one region that lists exactly the sets it
    belongs to), as in matplotlib-venn. The old drawing printed whole-set and pairwise-overlap totals instead.
    """
    set_names = list(sets.keys())
    set1, set2, set3 = [sets[name] for name in set_names]

    # Draw circles: three radius-1.2 circles whose centres are 1.2 apart, so all seven regions exist
    ax.add_patch(Circle((-0.6, 0.35), 1.2, facecolor="red", alpha=0.5))
    ax.add_patch(Circle((0.6, 0.35), 1.2, facecolor="blue", alpha=0.5))
    ax.add_patch(Circle((0.0, -0.7), 1.2, facecolor="green", alpha=0.5))

    regions = {
        (-1.15, 0.75): len(set1 - set2 - set3),
        (1.15, 0.75): len(set2 - set1 - set3),
        (0.0, -1.35): len(set3 - set1 - set2),
        (0.0, 0.95): len((set1 & set2) - set3),
        (-0.65, -0.35): len((set1 & set3) - set2),
        (0.65, -0.35): len((set2 & set3) - set1),
        (0.0, 0.1): len(set1 & set2 & set3),
    }
    for (x, y), count in regions.items():
        ax.text(x, y, f"{count}", ha="center", va="center")

    # Set names beside their circles
    ax.text(-1.5, 1.65, set_names[0], ha="center", va="bottom", fontweight="bold")
    ax.text(1.5, 1.65, set_names[1], ha="center", va="bottom", fontweight="bold")
    ax.text(0.0, -2.05, set_names[2], ha="center", va="top", fontweight="bold")

    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.5)
    ax.axis("off")


def _build_sankey_links(
    flows: List[Tuple[str, str, float]],
) -> Tuple[List[str], List[int], List[int], List[float]]:
    """Index the nodes of ``flows`` for a Sankey diagram.

    Returns ``(labels, sources, targets, values)`` where ``sources[i]`` and ``targets[i]`` are the positions in
    ``labels`` of the endpoints of flow ``i``. Labels are unique and kept in first-seen order, so an index
    never changes after it is assigned (Plotly resolves the indices against the finished ``labels`` list).
    """
    index: Dict[str, int] = {}
    sources: List[int] = []
    targets: List[int] = []
    values: List[float] = []
    for source, target, value in flows:
        for node in (source, target):
            if node not in index:
                index[node] = len(index)
        sources.append(index[source])
        targets.append(index[target])
        values.append(value)
    return list(index), sources, targets, values


def plot_sankey_diagram(
    flows: List[Tuple[str, str, float]],
    *,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs: Any,
) -> Any:
    """Create a Sankey diagram for flow visualization.

    With Plotly installed this is a Plotly Sankey diagram, written as HTML (the suffix of ``output_path`` is
    replaced by ``.html``). Without Plotly it is a matplotlib flow diagram with the sources on the left and the
    targets on the right, each labelled with its own name; line widths and marker areas grow linearly with the flow
    value, from the thinnest flow to the largest.

    Args:
        flows: List of (source, target, value) tuples
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height) of the matplotlib fallback
        **kwargs: Passed to Plotly's ``update_layout``; the matplotlib fallback does not use them

    Returns:
        matplotlib Figure or Plotly Figure object
    """
    validation.validate_type(flows, list, "flows")

    if HAS_PLOTLY:
        # Use Plotly for better Sankey diagrams
        sources: list[Any]
        targets: list[Any]
        values: list[Any]
        label_list, sources, targets, values = _build_sankey_links(flows)

        fig = go.Figure(
            data=[
                go.Sankey(
                    node=dict(pad=15, thickness=20, line=dict(color="black", width=0.5), label=label_list),
                    link=dict(source=sources, target=targets, value=values),
                )
            ]
        )

        fig.update_layout(title_text="Sankey Diagram", **kwargs)

        if output_path:
            paths.ensure_directory(Path(output_path).parent)
            html_path = Path(output_path).with_suffix(".html")
            fig.write_html(str(html_path))
            logger.info(f"Sankey diagram saved to {html_path}")

        return fig
    else:
        # Fallback matplotlib implementation
        fig, ax = plt.subplots(figsize=figsize)

        # Simple flow visualization; rows keep the order in which the nodes first appear in ``flows`` (a
        # ``set`` would order them by hash, so the figure would change from run to run)
        sources = list(dict.fromkeys(s for s, t, v in flows))
        targets = list(dict.fromkeys(t for s, t, v in flows))

        source_positions = {src: i for i, src in enumerate(sources)}
        target_positions = {tgt: i for i, tgt in enumerate(targets)}

        # Widths and marker areas are proportional to the flow's share of the largest flow, so the figure looks the
        # same whatever the units (a flow of 0.3 or of 3000 was drawn with a hair-thin or a huge line)
        largest = max((abs(v) for _, _, v in flows), default=0.0) or 1.0
        for source, target, value in flows:
            src_pos = source_positions[source]
            tgt_pos = target_positions[target]
            share = abs(value) / largest

            # Draw flow line
            ax.plot([0, 1], [src_pos, tgt_pos], "b-", linewidth=0.5 + 9.5 * share, alpha=0.7)
            ax.scatter([0], [src_pos], s=20 + 280 * share, c="red", alpha=0.7)
            ax.scatter([1], [tgt_pos], s=20 + 280 * share, c="blue", alpha=0.7)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Sources", "Targets"])
        ax.set_yticks(list(source_positions.values()))
        ax.set_yticklabels(list(source_positions.keys()))
        # The targets are plotted at their own heights on the right: label those heights with the target names (the
        # left labels used to be read for both columns, so a target was named after the source at the same height)
        right = ax.twinx()
        right.set_ylim(ax.get_ylim())
        right.set_yticks(list(target_positions.values()))
        right.set_yticklabels(list(target_positions.keys()))
        ax.set_title("Flow Diagram")

        if output_path:
            paths.ensure_directory(Path(output_path).parent)
            save_figure_deterministic(fig, output_path, dpi=300, bbox_inches="tight")
            logger.info(f"Flow diagram saved to {output_path}")

        return fig


def _chord_layout(
    matrix: np.ndarray, gap: float
) -> Tuple[List[Tuple[float, float]], Dict[Tuple[int, int], Tuple[float, float]]]:
    """Angles of a chord diagram, laid out as in d3-chord.

    Node ``i`` spans an angle proportional to the sum of row ``i``, in matrix order round the circle with a ``gap``
    (radians) after each node. Inside node ``i`` the entry ``matrix[i, j]`` takes a sub-arc of its own, in column
    order, with width proportional to the entry; the chord between ``i`` and ``j`` joins the sub-arc of ``i`` for
    ``j`` to the sub-arc of ``j`` for ``i``.

    Returns:
        ``(spans, ends)``: ``spans[i]`` is node ``i``'s ``(start, end)`` angle and ``ends[(i, j)]`` the ``(start,
        end)`` angle of the sub-arc of node ``i`` that belongs to ``matrix[i, j]``.
    """
    n = matrix.shape[0]
    unit = (2 * np.pi - n * gap) / float(matrix.sum())
    spans: List[Tuple[float, float]] = []
    ends: Dict[Tuple[int, int], Tuple[float, float]] = {}
    angle = 0.0
    for i in range(n):
        start = angle
        for j in range(n):
            width = float(matrix[i, j]) * unit
            ends[(i, j)] = (angle, angle + width)
            angle += width
        spans.append((start, angle))
        angle += gap
    return spans, ends


def _arc_points(start: float, end: float, radius: float, steps: int = 24) -> List[Tuple[float, float]]:
    """Points of the circular arc of ``radius`` from angle ``start`` to ``end``."""
    return [(radius * np.cos(t), radius * np.sin(t)) for t in np.linspace(start, end, steps)]


def _curve_points(
    start: Tuple[float, float], control: Tuple[float, float], end: Tuple[float, float], steps: int = 24
) -> List[Tuple[float, float]]:
    """Points of the quadratic Bezier curve from ``start`` to ``end`` bent toward ``control``."""
    t = np.linspace(0.0, 1.0, steps)[:, None]
    points = (1 - t) ** 2 * np.array(start) + 2 * (1 - t) * t * np.array(control) + t**2 * np.array(end)
    return [(float(x), float(y)) for x, y in points]


def plot_chord_diagram(
    matrix: np.ndarray,
    labels: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 8),
    **kwargs: Any,
) -> Axes:
    """Create a chord diagram for showing relationships between categories.

    The nodes are arcs of a circle, each as long as its row sum; the entry ``matrix[i, j]`` takes a sub-arc of node
    ``i`` of proportional width, and a ribbon joins the sub-arc of ``i`` for ``j`` to the sub-arc of ``j`` for ``i``
    (so a symmetric matrix gives ribbons that are as wide at both ends, and an asymmetric one a different width
    at each end, one per direction). The diagonal widens its node's arc but has no ribbon. Both triangles of the
    matrix are used: the upper triangle alone used to be read, and the ribbons were distorted shapes (a Bezier path
    with the wrong number of control points, in angle/radius coordinates).

    Args:
        matrix: Square matrix of non-negative relationships
        labels: Optional labels for matrix axes
        ax: Optional matplotlib axes to plot on, Cartesian or polar. Creates a Cartesian axes if None.
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: ``gap`` (radians between nodes, default 0.06) and ``colors`` (one colour per node, default the
            colour-blind-safe palette); any other keyword argument (``alpha``, ``linewidth``, ...) is passed to
            matplotlib fill() for every ribbon.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If the matrix is not square, has a negative or non-finite entry or no positive one, or ``labels``
            or ``colors`` do not have one entry per row
    """
    validation.validate_type(matrix, np.ndarray, "matrix")

    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Chord diagram matrix must be square")
    n = matrix.shape[0]
    if labels is not None and len(labels) != n:
        raise ValueError(f"Number of labels ({len(labels)}) must match matrix dimension ({n})")
    values = np.asarray(matrix, dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Chord diagram matrix entries must be finite and non-negative")
    if not (values.sum() > 0):
        raise ValueError("Chord diagram matrix has no positive entry to draw")

    gap = float(kwargs.pop("gap", 0.06))
    colors: Sequence[Any] = kwargs.pop("colors", None) or palettes.categorical(n)
    if len(colors) != n:
        raise ValueError(f"Number of colors ({len(colors)}) must match matrix dimension ({n})")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    polar = ax.name == "polar"

    def to_axes(points: List[Tuple[float, float]]) -> Tuple[List[float], List[float]]:
        """Cartesian points in the coordinates ``ax`` reads: (x, y), or (angle, radius) on a polar axes."""
        if polar:
            return [float(np.arctan2(y, x)) for x, y in points], [float(np.hypot(x, y)) for x, y in points]
        return [x for x, _ in points], [y for _, y in points]

    spans, ends = _chord_layout(values, gap)
    ribbon_kwargs: Dict[str, Any] = {"alpha": 0.5, "linewidth": 0.5, **kwargs}
    radius, ring = 1.0, 0.06

    # Ribbons: arc on node i, curve through the centre to node j's arc, arc, curve back
    for i in range(n):
        for j in range(i + 1, n):
            if values[i, j] <= 0 and values[j, i] <= 0:
                continue
            (a0, a1), (b0, b1) = ends[(i, j)], ends[(j, i)]
            outline = (
                _arc_points(a0, a1, radius)
                + _curve_points(_arc_points(a1, a1, radius, 1)[0], (0.0, 0.0), _arc_points(b0, b0, radius, 1)[0])
                + _arc_points(b0, b1, radius)
                + _curve_points(_arc_points(b1, b1, radius, 1)[0], (0.0, 0.0), _arc_points(a0, a0, radius, 1)[0])
            )
            xs, ys = to_axes(outline)
            ax.fill(xs, ys, **{"facecolor": colors[i], "edgecolor": colors[i], **ribbon_kwargs})

    # Nodes: a ring segment over each node's span, and its label just outside
    for i, (start, end) in enumerate(spans):
        xs, ys = to_axes(_arc_points(start, end, radius + ring) + _arc_points(end, start, radius))
        ax.fill(xs, ys, facecolor=colors[i], edgecolor=colors[i], linewidth=0.5)
        middle = (start + end) / 2
        label_xy = [(radius + 0.14) * np.cos(middle), (radius + 0.14) * np.sin(middle)]
        label_pos = to_axes([(float(label_xy[0]), float(label_xy[1]))])
        ax.text(
            label_pos[0][0],
            label_pos[1][0],
            labels[i] if labels else f"Node {i+1}",
            ha="left" if label_xy[0] >= 0 else "right",
            va="center",
            fontsize=10,
        )

    if polar:
        cast("Any", ax).set_rlim(0, 1.5)
    else:
        ax.set_xlim(-1.6, 1.6)
        ax.set_ylim(-1.6, 1.6)
        ax.set_aspect("equal")
    ax.set_title("Chord Diagram")
    ax.axis("off")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Chord diagram saved to {output_path}")

    return ax


def plot_alluvial_diagram(
    data: pd.DataFrame,
    stages: List[str],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs: Any,
) -> Axes:
    """Create an alluvial diagram (flow diagram) for showing transitions.

    With the optional ``alluvial`` package installed, it draws the diagram. Otherwise the built-in drawing is used:
    one bar per category in each stage (as tall as its number of rows) and one ribbon per transition between
    adjacent stages (as tall as the number of rows making it). Each row of ``data`` is one observation.

    Args:
        data: DataFrame with columns for each stage
        stages: List of column names representing stages
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Passed to ``alluvial.plot``; the built-in drawing does not use them and says so in a warning

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(data, pd.DataFrame, "data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    try:
        import alluvial
    except ImportError:
        # Fallback implementation
        _plot_simple_alluvial(ax, data, stages, **kwargs)
    else:
        # Use alluvial package
        alluvial.plot(data, stages, ax=ax, **kwargs)

    ax.set_title("Alluvial Diagram")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Alluvial diagram saved to {output_path}")

    return ax


def _plot_simple_alluvial(ax: Axes, data: pd.DataFrame, stages: List[str], **kwargs: Any) -> None:
    """Built-in alluvial diagram: a bar per category in each stage and a ribbon per transition between stages.

    A bar is as tall as the number of rows that have its category at that stage (missing values form a "missing"
    category), and the ribbon between a category of one stage and a category of the next is as tall as the number of
    rows that go from one to the other, so the bars and ribbons at both ends of a ribbon account for every row.
    The old drawing was a parallel-coordinates line plot with one line per row, which showed no flow at all.
    """
    if kwargs:
        logger.warning(f"The built-in alluvial drawing ignores {sorted(kwargs)}")
    categories = [data[stage].map(lambda v: "missing" if pd.isna(v) else str(v)) for stage in stages]
    n_rows = len(data)
    bar_width, gap = 0.12, 0.02 * n_rows

    # Stack the categories of each stage (sorted by name) from the bottom up
    stacks: List[Dict[str, Tuple[float, float]]] = []  # category -> (bottom, row count)
    for column in categories:
        counts = column.value_counts().sort_index()
        bottom, stack = 0.0, {}
        for category, count in counts.items():
            stack[str(category)] = (bottom, float(count))
            bottom += float(count) + gap
        stacks.append(stack)

    for k, stack in enumerate(stacks):
        colors = dict(zip(stack, palettes.categorical(len(stack))))
        for category, (low, count) in stack.items():
            ax.add_patch(Rectangle((k - bar_width / 2, low), bar_width, count, facecolor=colors[category]))
            ax.text(
                k,
                low + count / 2,
                category,
                ha="center",
                va="center",
                fontsize=8,
                color="black",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.7, linewidth=0),
            )
        if k == len(stacks) - 1:
            break
        # Ribbons from stage k to stage k + 1, filling each bar from the bottom in partner order
        pairs = pd.crosstab(categories[k], categories[k + 1])
        used_out = {category: low for category, (low, _) in stack.items()}
        used_in = {category: low for category, (low, _) in stacks[k + 1].items()}
        for source in pairs.index:
            for target in pairs.columns:
                count = float(pairs.loc[source, target])
                if count <= 0:
                    continue
                y0, y1 = used_out[source], used_in[target]
                used_out[source] += count
                used_in[target] += count
                x0, x1 = k + bar_width / 2, k + 1 - bar_width / 2
                mid = (x0 + x1) / 2
                vertices = [
                    (x0, y0 + count), (mid, y0 + count), (mid, y1 + count), (x1, y1 + count),
                    (x1, y1),
                    (mid, y1), (mid, y0), (x0, y0),
                    (x0, y0 + count),
                ]  # fmt: skip
                codes = [MplPath.MOVETO] + [MplPath.CURVE4] * 3 + [MplPath.LINETO] + [MplPath.CURVE4] * 3
                codes.append(MplPath.CLOSEPOLY)
                color = dict(zip(stack, palettes.categorical(len(stack))))[source]
                ax.add_patch(PathPatch(MplPath(vertices, codes), facecolor=color, edgecolor="none", alpha=0.5))

    ax.set_xlim(-0.5, len(stages) - 0.5)
    ax.set_ylim(0, max(low + count for stack in stacks for low, count in stack.values()) + gap)
    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels(stages)
    ax.set_yticks([])
    ax.set_xlabel("Stages")
    ax.set_ylabel(f"Rows ({n_rows} in every stage)")


def plot_circular_barplot(
    values: np.ndarray,
    labels: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 8),
    **kwargs: Any,
) -> Axes:
    """Create a circular bar plot.

    Args:
        values: Values for each bar
        labels: Optional labels for each bar
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(values, np.ndarray, "values")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": "polar"})

    n_bars = len(values)
    angles = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)

    # Plot bars (caller keyword arguments override the defaults)
    bar_kwargs: dict[str, Any] = {"width": 0.4, "bottom": 0.0, "alpha": 0.7, **kwargs}
    ax.bar(angles, values, **bar_kwargs)

    # Add labels
    if labels:
        for angle, label in zip(angles, labels):
            ax.text(
                angle, max(values) * 1.1, label, ha="center", va="center", rotation=np.degrees(angle) - 90, fontsize=8
            )

    cast("Any", ax).set_rlabel_position(0)
    ax.set_title("Circular Bar Plot")
    ax.grid(True, alpha=0.3)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Circular bar plot saved to {output_path}")

    return ax


def plot_network_circular_layout(
    graph: Any,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 10),
    **kwargs: Any,
) -> Axes:
    """Plot a network in circular layout.

    Args:
        graph: NetworkX graph object
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx required for network visualization")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": "polar"})

    # Create circular layout
    n_nodes = len(graph.nodes())
    angles = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False)
    pos = {node: (angle, 1) for node, angle in zip(graph.nodes(), angles)}

    # Draw network
    nx.draw_networkx_nodes(graph, pos, node_size=300, alpha=0.8, ax=ax)
    nx.draw_networkx_edges(graph, pos, alpha=0.5, ax=ax)

    # Add labels
    for node, (angle, r) in pos.items():
        ax.text(angle, r + 0.1, str(node), ha="center", va="center", fontsize=8, rotation=np.degrees(angle) - 90)

    ax.set_title("Circular Network Layout")
    ax.grid(False)
    cast("Any", ax).set_rlabel_position(0)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Circular network plot saved to {output_path}")

    return ax


def plot_upset_plot(
    data: Dict[str, set],
    *,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs: Any,
) -> Any:
    """Create an UpSet plot for set intersections.

    Args:
        data: Dictionary mapping set names to sets
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Figure object
    """
    try:
        import upsetplot
    except ImportError:
        raise ImportError("upsetplot required for UpSet plot visualization")

    # Create UpSet data
    from upsetplot import from_contents

    upset_data = from_contents(data)

    # Create plot
    fig = plt.figure(figsize=figsize)
    upsetplot.plot(upset_data, fig=fig, **kwargs)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(fig, output_path, dpi=300, bbox_inches="tight")
        logger.info(f"UpSet plot saved to {output_path}")

    return fig
