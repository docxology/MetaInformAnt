"""Interactive dashboard data structures for web-based visualization.

Generates JSON-serializable plot data structures for interactive scatter
plots, heatmaps, genome browser tracks, volcano plots, and composed
dashboards. Includes standalone HTML export using embedded SVG and
vanilla JavaScript (no external dependencies).

All functions produce pure Python dictionaries ready for JSON serialization
and rendering by any web framework or standalone HTML viewer.
"""

from __future__ import annotations

import html
import json
import math
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import leaves_list, linkage, optimal_leaf_ordering

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _json_for_script(data: object) -> str:
    """JSON text that is safe to place inside an HTML ``<script>`` element.

    ``</script>`` (or ``<!--``) inside a string value would end the script block and let the rest of the value be
    parsed as markup, so ``<``, ``>`` and ``&`` are written as their JSON ``\\u`` escapes, which decode to the same
    characters.
    """
    text = json.dumps(data, default=str)
    return text.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")


def _is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _scale(value: float, low: float, high: float, start: float, end: float) -> float:
    """Map ``value`` from ``[low, high]`` to ``[start, end]``; a degenerate range maps to the middle."""
    if not high > low:
        return (start + end) / 2.0
    return start + (value - low) / (high - low) * (end - start)


def _panel_svg(plot_data: dict, width: int = 480, height: int = 300) -> str | None:
    """Draw ``plot_data`` from a ``create_interactive_*`` function as an SVG string, or None if it is not drawable.

    * scatter and volcano data (points with ``x`` and ``y``): one circle per point, filled with the point's
      ``color``;
    * heatmap data (cells with ``row_idx``, ``col_idx`` and ``value``): one rectangle per cell on a blue-white-red
      scale spanning the values;
    * genome browser data (``track_data`` features with ``start`` and ``end``): one bar per feature.

    Every label is HTML-escaped and attached as an SVG ``<title>``, which browsers show on hover.
    """
    margin = 30
    inner_w, inner_h = width - 2 * margin, height - 2 * margin
    body: list[str] = []
    points = plot_data.get("plot_data")

    if isinstance(points, list) and points and all(isinstance(p, dict) and "x" in p and "y" in p for p in points):
        xs = [p["x"] for p in points if _is_number(p["x"])]
        ys = [p["y"] for p in points if _is_number(p["y"])]
        if not xs or not ys:
            return None
        x_low, x_high, y_low, y_high = min(xs), max(xs), min(ys), max(ys)
        for index, point in enumerate(points):
            if not (_is_number(point["x"]) and _is_number(point["y"])):
                continue
            cx = _scale(point["x"], x_low, x_high, margin, margin + inner_w)
            cy = _scale(point["y"], y_low, y_high, margin + inner_h, margin)
            label = point.get("gene") or point.get("label") or f"Point {index}"
            color = html.escape(str(point.get("color", "#1f77b4")), quote=True)
            body.append(
                f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="4" fill="{color}">'
                f"<title>{html.escape(str(label))} ({point['x']:.3g}, {point['y']:.3g})</title></circle>"
            )
        body.append(f'<text x="{margin}" y="{height - 8}" font-size="10">x: {x_low:.3g} to {x_high:.3g}</text>')
        body.append(f'<text x="{margin}" y="{margin - 8}" font-size="10">y: {y_low:.3g} to {y_high:.3g}</text>')
    elif (
        isinstance(points, list)
        and points
        and all(isinstance(p, dict) and {"row_idx", "col_idx", "value"} <= set(p) for p in points)
    ):
        n_rows = max(p["row_idx"] for p in points) + 1
        n_cols = max(p["col_idx"] for p in points) + 1
        values = [p["value"] for p in points if _is_number(p["value"])]
        low, high = (min(values), max(values)) if values else (0.0, 1.0)
        cell_w, cell_h = inner_w / n_cols, inner_h / n_rows
        for cell in points:
            if not _is_number(cell["value"]):
                continue
            t = _scale(cell["value"], low, high, 0.0, 1.0)
            # blue (low) -> white (middle) -> red (high)
            red = 255 if t >= 0.5 else int(255 * 2 * t)
            blue = 255 if t <= 0.5 else int(255 * 2 * (1 - t))
            green = int(255 * (1 - abs(2 * t - 1)))
            body.append(
                f'<rect x="{margin + cell["col_idx"] * cell_w:.2f}" y="{margin + cell["row_idx"] * cell_h:.2f}" '
                f'width="{cell_w:.2f}" height="{cell_h:.2f}" fill="rgb({red},{green},{blue})">'
                f"<title>{html.escape(str(cell.get('row')))}, {html.escape(str(cell.get('col')))}: "
                f"{cell['value']:.3g}</title></rect>"
            )
    elif isinstance(plot_data.get("track_data"), list) and plot_data["track_data"]:
        features = plot_data["track_data"]
        low, high = plot_data["coordinate_range"]
        bar_h = min(24.0, inner_h)
        for feature in features:
            x0 = _scale(feature["start"], low, high, margin, margin + inner_w)
            x1 = _scale(feature["end"], low, high, margin, margin + inner_w)
            color = html.escape(str(feature.get("color", "#636363")), quote=True)
            body.append(
                f'<rect x="{x0:.2f}" y="{margin + (inner_h - bar_h) / 2:.2f}" width="{max(x1 - x0, 1.0):.2f}" '
                f'height="{bar_h:.2f}" fill="{color}"><title>{html.escape(str(feature.get("name")))} '
                f"({feature['start']}-{feature['end']})</title></rect>"
            )
        body.append(f'<text x="{margin}" y="{height - 8}" font-size="10">{low} to {high}</text>')
    else:
        return None

    return (
        f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" role="img">'
        f'<rect width="{width}" height="{height}" fill="white"/>' + "".join(body) + "</svg>"
    )


def _mean(values: list[float]) -> float:
    """Arithmetic mean."""
    return sum(values) / len(values) if values else 0.0


def _std(values: list[float]) -> float:
    """Sample standard deviation."""
    n = len(values)
    if n < 2:
        return 0.0
    mu = _mean(values)
    return math.sqrt(sum((x - mu) ** 2 for x in values) / (n - 1))


def _cluster_axis(data: list[list[float]], axis: int = 0) -> tuple[list[int], list[list[float]]]:
    """Hierarchically cluster the rows (``axis=0``) or columns (``axis=1``) of *data*.

    Average-linkage (UPGMA) clustering on Euclidean distance, with the leaves ordered so that the order is the
    smoothest one the tree allows (``optimal_leaf_ordering``). Every cluster occupies one contiguous block of the
    returned order, which the greedy nearest-neighbour chain this replaces did not guarantee: it appended a late
    neighbour of an early row at the very end, splitting their cluster.

    Returns:
        ``(order, merges)``. ``order`` lists the original indices in display order. ``merges`` is the dendrogram
        in SciPy's linkage format, one ``[left, right, height, size]`` row per merge: an id below the number of
        leaves is an original index, and id ``n + i`` is the cluster formed by merge ``i``. It is empty when there
        is nothing to merge.

    Raises:
        ValueError: If the data contains NaN or infinite values.
    """
    matrix = np.asarray(data, dtype=float)
    vectors = matrix.T if axis == 1 else matrix
    n = len(vectors)
    if n <= 1:
        return list(range(n)), []
    if not np.isfinite(vectors).all():
        raise ValueError("Cannot cluster data with NaN or infinite values; fill them or pass cluster_rows/cols=False")

    merges = optimal_leaf_ordering(linkage(vectors, method="average", metric="euclidean"), vectors)
    return [int(i) for i in leaves_list(merges)], merges.tolist()


# ---------------------------------------------------------------------------
# Interactive Scatter
# ---------------------------------------------------------------------------


def create_interactive_scatter(
    x: list[float],
    y: list[float],
    labels: list[str] | None = None,
    colors: list[str] | None = None,
    hover_data: dict | None = None,
    title: str = "Scatter Plot",
) -> dict:
    """Generate interactive scatter plot data structure.

    Args:
        x: X-axis values.
        y: Y-axis values.
        labels: Optional point labels for hover display.
        colors: Optional per-point color values (hex or named).
        hover_data: Optional dict of additional hover fields,
            each mapping field name to list of values.
        title: Plot title.

    Returns:
        Dictionary with keys:
            - plot_data: Point data with coordinates, labels, colors.
            - layout: Layout configuration (title, axes, margins).
            - config: Interaction configuration.

    Raises:
        ValueError: If x and y have different lengths.
    """
    if len(x) != len(y):
        raise ValueError(f"x ({len(x)}) and y ({len(y)}) must have same length")

    n = len(x)

    if labels is None:
        labels = [f"Point {i}" for i in range(n)]
    if colors is None:
        colors = ["#1f77b4"] * n

    points = []
    for i in range(n):
        point = {
            "x": x[i],
            "y": y[i],
            "label": labels[i],
            "color": colors[i],
        }
        if hover_data:
            for field, values in hover_data.items():
                if i < len(values):
                    point[field] = values[i]
        points.append(point)

    x_range = (min(x), max(x)) if x else (0, 1)
    y_range = (min(y), max(y)) if y else (0, 1)
    x_pad = (x_range[1] - x_range[0]) * 0.05 or 1.0
    y_pad = (y_range[1] - y_range[0]) * 0.05 or 1.0

    layout = {
        "title": title,
        "x_axis": {"label": "X", "range": (x_range[0] - x_pad, x_range[1] + x_pad)},
        "y_axis": {"label": "Y", "range": (y_range[0] - y_pad, y_range[1] + y_pad)},
        "width": 800,
        "height": 600,
        "margin": {"top": 50, "right": 30, "bottom": 50, "left": 60},
    }

    config = {
        "hover_enabled": True,
        "zoom_enabled": True,
        "pan_enabled": True,
        "point_radius": 5,
    }

    logger.info("Created interactive scatter: %d points", n)

    return {"plot_data": points, "layout": layout, "config": config}


# ---------------------------------------------------------------------------
# Interactive Heatmap
# ---------------------------------------------------------------------------


def create_interactive_heatmap(
    data: list[list[float]],
    row_names: list[str],
    col_names: list[str],
    cluster_rows: bool = True,
    cluster_cols: bool = True,
) -> dict:
    """Generate interactive heatmap with optional clustering dendrograms.

    Args:
        data: 2D matrix of values.
        row_names: Row labels.
        col_names: Column labels.
        cluster_rows: Whether to cluster rows (average-linkage hierarchical clustering, see :func:`_cluster_axis`).
        cluster_cols: Whether to cluster columns.

    Returns:
        Dictionary with keys:
            - plot_data: Heatmap cell data with row, col, value.
            - row_order: Ordered row indices.
            - col_order: Ordered column indices.
            - dendrograms: ``{"rows": ..., "cols": ...}``, each ``{"order", "clustered", "linkage"}``, where
              ``linkage`` is the dendrogram as ``[left, right, height, size]`` merge rows (SciPy format; empty when
              not clustered).
            - color_range: (min_value, max_value).

    Raises:
        ValueError: If the labels do not match the data, or clustering is requested for data with NaN or infinite
            values.
    """
    n_rows = len(data)
    n_cols = len(data[0]) if data else 0

    if len(row_names) != n_rows:
        raise ValueError(f"row_names ({len(row_names)}) != data rows ({n_rows})")
    if len(col_names) != n_cols:
        raise ValueError(f"col_names ({len(col_names)}) != data cols ({n_cols})")

    row_order, row_linkage = _cluster_axis(data, axis=0) if cluster_rows else (list(range(n_rows)), [])
    col_order, col_linkage = _cluster_axis(data, axis=1) if cluster_cols else (list(range(n_cols)), [])

    cells = []
    all_values = []
    for ri, r in enumerate(row_order):
        for ci, c in enumerate(col_order):
            val = data[r][c]
            cells.append(
                {
                    "row": row_names[r],
                    "col": col_names[c],
                    "value": val,
                    "row_idx": ri,
                    "col_idx": ci,
                }
            )
            all_values.append(val)

    color_range = (min(all_values), max(all_values)) if all_values else (0.0, 1.0)

    dendrograms = {
        "rows": {"order": row_order, "clustered": cluster_rows, "linkage": row_linkage},
        "cols": {"order": col_order, "clustered": cluster_cols, "linkage": col_linkage},
    }

    logger.info("Created interactive heatmap: %dx%d", n_rows, n_cols)

    return {
        "plot_data": cells,
        "row_order": row_order,
        "col_order": col_order,
        "dendrograms": dendrograms,
        "color_range": color_range,
    }


# ---------------------------------------------------------------------------
# Genome Browser Track
# ---------------------------------------------------------------------------


def create_genome_browser_track(
    features: list[dict],
    track_type: str = "gene",
) -> dict:
    """Generate genome browser track data.

    Args:
        features: List of feature dicts, each with ``start``, ``end``, and
            ``name`` keys. Optional: ``strand`` (default ``"+"``), ``score``
            (default 0), ``chromosome`` (default None: unknown, not a guessed
            chromosome), ``color`` (default by track type).
        track_type: Track type: ``"gene"``, ``"variant"``, ``"peak"``,
            ``"coverage"``.

    Returns:
        Dictionary with keys:
            - track_data: Processed feature list.
            - track_type: Track type.
            - coordinate_range: (min_start, max_end); (0, 0) for a track with no features.
            - n_features: Number of features.

    Raises:
        ValueError: If a feature lacks ``start``, ``end`` or ``name``. They used to be replaced by position 0,
            ``start + 1`` and "unnamed", which drew features that were never in the input.
    """
    processed = []
    min_start = float("inf")
    max_end = float("-inf")

    for index, feat in enumerate(features):
        missing = [key for key in ("start", "end", "name") if key not in feat]
        if missing:
            raise ValueError(f"feature {index} is missing the required key(s) {missing}")
        start = feat["start"]
        end = feat["end"]
        name = feat["name"]

        min_start = min(min_start, start)
        max_end = max(max_end, end)

        entry = {
            "start": start,
            "end": end,
            "name": name,
            "strand": feat.get("strand", "+"),
            "score": feat.get("score", 0),
            "chromosome": feat.get("chromosome"),
            "color": feat.get("color", _track_color(track_type)),
            "width": end - start,
        }
        processed.append(entry)

    if not processed:
        min_start = 0
        max_end = 0

    logger.info(
        "Created genome browser track (%s): %d features, range %d-%d",
        track_type,
        len(processed),
        int(min_start),
        int(max_end),
    )

    return {
        "track_data": processed,
        "track_type": track_type,
        "coordinate_range": (int(min_start), int(max_end)),
        "n_features": len(processed),
    }


def _track_color(track_type: str) -> str:
    """Default color for a track type."""
    colors = {
        "gene": "#2171b5",
        "variant": "#d94801",
        "peak": "#238b45",
        "coverage": "#6a51a3",
    }
    return colors.get(track_type, "#636363")


# ---------------------------------------------------------------------------
# Interactive Volcano
# ---------------------------------------------------------------------------


def create_interactive_volcano(
    de_results: list[dict],
    fc_threshold: float = 1.0,
    p_threshold: float = 0.05,
) -> dict:
    """Generate interactive volcano plot data.

    Args:
        de_results: List of differential expression result dicts, each with
            ``gene``, ``log2_fold_change`` (or ``log2fc``), and ``p_value``
            (or ``pvalue`` or ``adjusted_p``) keys.
        fc_threshold: Absolute log2 fold change threshold for significance.
        p_threshold: P-value threshold for significance.

    Returns:
        Dictionary with keys:
            - plot_data: Point data with x (log2FC), y (-log10 p), gene, status.
            - significance_counts: Dict with up, down, ns counts.
            - thresholds: Applied threshold values.

    Raises:
        ValueError: If a result has no fold change or no p-value. A missing p-value used to count as 1.0 and a
            missing fold change as 0.0, which reported genes as "not significant" that were never tested.
    """
    points = []
    counts = {"up": 0, "down": 0, "not_significant": 0}

    for index, result in enumerate(de_results):
        gene = result.get("gene", result.get("name", "unnamed"))
        fc_keys = [key for key in ("log2_fold_change", "log2fc") if key in result]
        p_keys = [key for key in ("p_value", "pvalue", "adjusted_p") if key in result]
        if not fc_keys:
            raise ValueError(f"result {index} ({gene!r}) has no fold change: expected 'log2_fold_change' or 'log2fc'")
        if not p_keys:
            raise ValueError(f"result {index} ({gene!r}) has no p-value: expected 'p_value', 'pvalue' or 'adjusted_p'")
        fc = result[fc_keys[0]]
        pval = result[p_keys[0]]

        if pval <= 0:
            neg_log_p = 300.0
        elif pval >= 1:
            neg_log_p = 0.0
        else:
            neg_log_p = -math.log10(pval)

        if pval < p_threshold and abs(fc) > fc_threshold:
            if fc > 0:
                status = "up"
                color = "#e41a1c"
                counts["up"] += 1
            else:
                status = "down"
                color = "#377eb8"
                counts["down"] += 1
        else:
            status = "not_significant"
            color = "#999999"
            counts["not_significant"] += 1

        points.append(
            {
                "gene": gene,
                "x": fc,
                "y": neg_log_p,
                "p_value": pval,
                "log2fc": fc,
                "status": status,
                "color": color,
            }
        )

    logger.info(
        "Created interactive volcano: %d genes (up=%d, down=%d, ns=%d)",
        len(points),
        counts["up"],
        counts["down"],
        counts["not_significant"],
    )

    return {
        "plot_data": points,
        "significance_counts": counts,
        "thresholds": {
            "fc_threshold": fc_threshold,
            "p_threshold": p_threshold,
            "neg_log_p_threshold": -math.log10(p_threshold) if p_threshold > 0 else 0,
        },
    }


# ---------------------------------------------------------------------------
# HTML Export
# ---------------------------------------------------------------------------


def export_to_html(
    plot_data: dict,
    output_path: str,
    title: str = "Interactive Plot",
) -> str:
    """Export interactive plot to standalone HTML file.

    Generates a self-contained UTF-8 HTML file. Scatter-like data (points with ``x`` and ``y``, such as
    :func:`create_interactive_scatter` and :func:`create_interactive_volcano`) is drawn with embedded SVG and vanilla
    JavaScript for hover tooltips. Heatmap and genome-browser data are drawn as static SVG with native hover
    titles. Any other structure is embedded as data but not drawn, and a warning says so.

    The title and every label are HTML-escaped, and the data is embedded so that no value can close the script
    block: untrusted gene names or titles cannot inject markup.

    Args:
        plot_data: Plot data structure from any ``create_interactive_*`` function.
        output_path: Output file path.
        title: HTML page title.

    Returns:
        Absolute path of the generated HTML file.
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    data_json = _json_for_script(plot_data)
    safe_title = html.escape(title)

    points = plot_data.get("plot_data")
    scatter_like = bool(
        isinstance(points, list) and points and isinstance(points[0], dict) and "x" in points[0] and "y" in points[0]
    )
    static_svg = "" if scatter_like else (_panel_svg(plot_data, 800, 600) or "")
    if not scatter_like and not static_svg:
        logger.warning("export_to_html cannot draw this plot structure; the file embeds the data only")

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{safe_title}</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
               margin: 0; padding: 20px; background: #fafafa; }}
        .container {{ max-width: 1200px; margin: 0 auto; }}
        h1 {{ color: #333; font-size: 1.5em; margin-bottom: 10px; }}
        .plot-area {{ background: white; border: 1px solid #ddd; border-radius: 4px;
                      padding: 20px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
        svg {{ width: 100%; height: auto; }}
        .tooltip {{ position: absolute; background: rgba(0,0,0,0.8); color: white;
                    padding: 6px 10px; border-radius: 4px; font-size: 12px;
                    pointer-events: none; display: none; z-index: 1000; }}
        .point:hover {{ opacity: 0.7; cursor: pointer; }}
        .stats {{ margin-top: 10px; color: #666; font-size: 0.9em; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{safe_title}</h1>
        <div class="plot-area" id="plot">{static_svg}</div>
        <div class="stats" id="stats"></div>
    </div>
    <div class="tooltip" id="tooltip"></div>
    <script>
        const plotData = {data_json};
        const container = document.getElementById('plot');
        const tooltip = document.getElementById('tooltip');
        const stats = document.getElementById('stats');

        // Escape text and attribute values: labels and colours come from the data
        const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
                                  .replace(/"/g, '&quot;').replace(/'/g, '&#39;');

        // Render based on data structure
        if (plotData.plot_data && Array.isArray(plotData.plot_data)) {{
            const points = plotData.plot_data;
            const n = points.length;
            stats.textContent = 'Data points: ' + n;

            if (points[0] && 'x' in points[0] && 'y' in points[0]) {{
                // Scatter-like plot
                const w = 800, h = 600, m = {{t:40,r:30,b:50,l:60}};
                const pw = w - m.l - m.r, ph = h - m.t - m.b;

                const xs = points.map(p => p.x);
                const ys = points.map(p => p.y);
                const xMin = Math.min(...xs), xMax = Math.max(...xs);
                const yMin = Math.min(...ys), yMax = Math.max(...ys);
                const xRange = xMax - xMin || 1;
                const yRange = yMax - yMin || 1;

                let svg = '<svg viewBox="0 0 '+w+' '+h+'" xmlns="http://www.w3.org/2000/svg">';
                svg += '<rect width="'+w+'" height="'+h+'" fill="white"/>';

                // Points
                points.forEach((p, i) => {{
                    const cx = m.l + ((p.x - xMin) / xRange) * pw;
                    const cy = m.t + ph - ((p.y - yMin) / yRange) * ph;
                    const color = esc(p.color || '#1f77b4');
                    const label = esc(p.gene || p.label || ('Point ' + i));
                    svg += '<circle class="point" cx="'+cx+'" cy="'+cy+'" r="4" fill="'+color+'"' +
                           ' data-label="'+label+'" data-x="'+p.x.toFixed(3)+'" data-y="'+p.y.toFixed(3)+'"/>';
                }});

                // Axes
                svg += (
                    '<line x1="'+m.l+'" y1="'+(m.t+ph)+'" x2="'+(m.l+pw)+'" '
                    'y2="'+(m.t+ph)+'" stroke="#333" stroke-width="1"/>'
                );
                svg += (
                    '<line x1="'+m.l+'" y1="'+m.t+'" x2="'+m.l+'" '
                    'y2="'+(m.t+ph)+'" stroke="#333" stroke-width="1"/>'
                );
                svg += '</svg>';

                container.innerHTML = svg;

                // Hover tooltips
                container.querySelectorAll('.point').forEach(el => {{
                    el.addEventListener('mouseover', e => {{
                        tooltip.style.display = 'block';
                        tooltip.textContent = el.dataset.label + ' (' + el.dataset.x + ', ' + el.dataset.y + ')';
                        tooltip.style.left = (e.pageX + 10) + 'px';
                        tooltip.style.top = (e.pageY - 25) + 'px';
                    }});
                    el.addEventListener('mouseout', () => {{ tooltip.style.display = 'none'; }});
                }});
            }}
        }}
    </script>
</body>
</html>"""

    out.write_text(html_content, encoding="utf-8")
    abs_path = str(out.resolve())
    logger.info("Exported interactive plot to %s", abs_path)
    return abs_path


# ---------------------------------------------------------------------------
# Dashboard Composition
# ---------------------------------------------------------------------------

_DASHBOARD_LAYOUTS = ("grid", "rows", "columns")


def create_dashboard(
    panels: list[dict],
    layout: str = "grid",
    title: str = "Dashboard",
) -> dict:
    """Compose multiple plots into a dashboard layout.

    Each panel is drawn into the returned HTML: scatter and volcano data as circles, heatmap data as coloured
    cells, genome-browser data as bars (see :func:`_panel_svg`), each with hover titles. A panel whose
    ``plot_data`` has any other structure is shown as its data, truncated with a visible marker and labelled "not
    drawn". The title, panel titles and all labels are HTML-escaped.

    Args:
        panels: List of panel dicts, each with ``title`` and ``plot_data``
            keys (from any ``create_interactive_*`` function). Optional
            ``width`` (pixels, default: an equal share of the row) and ``height`` (pixels) keys.
        layout: Layout strategy: ``"grid"`` (auto-arrange), ``"rows"``
            (one per row), or ``"columns"`` (side-by-side).
        title: Dashboard title.

    Returns:
        Dictionary with keys:
            - dashboard_data: Structured panel data with positions.
            - html: Standalone HTML string for the dashboard.
            - n_panels: Number of panels.
            - layout: Layout used.

    Raises:
        ValueError: If ``layout`` is not ``"grid"``, ``"rows"`` or ``"columns"`` (it used to fall back to
            ``"grid"`` without a word).
    """
    if layout not in _DASHBOARD_LAYOUTS:
        raise ValueError(f"Unknown dashboard layout {layout!r}; expected one of {list(_DASHBOARD_LAYOUTS)}")

    n = len(panels)
    if n == 0:
        cols, rows = 1, 1
    elif layout == "grid":
        cols = math.ceil(math.sqrt(n))
        rows = math.ceil(n / cols)
    elif layout == "rows":
        cols = 1
        rows = n
    else:  # columns
        cols = n
        rows = 1

    panel_width = 100.0 / cols
    panel_height = max(300, 600 // rows)

    positioned_panels = []
    for idx, panel in enumerate(panels):
        row = idx // cols
        col = idx % cols
        size = {
            "width_pct": panel_width,
            "height_px": panel.get("height", panel_height),
        }
        if "width" in panel:
            size["width_px"] = panel["width"]
        positioned_panels.append(
            {
                "title": panel.get("title", f"Panel {idx + 1}"),
                "plot_data": panel.get("plot_data", {}),
                "position": {"row": row, "col": col},
                "size": size,
            }
        )

    # Simple HTML dashboard: every panel is drawn, not dumped as JSON
    html_panels = []
    for pp in positioned_panels:
        svg = _panel_svg(pp["plot_data"], int(pp["size"].get("width_px", 480)), int(pp["size"]["height_px"]))
        if svg is None:
            dumped = json.dumps(pp["plot_data"], indent=1, default=str)
            if len(dumped) > 2000:
                dumped = dumped[:2000] + "\n... (truncated)"
            content = (
                '<p style="margin:0 0 5px 0;font-size:11px;color:#666;">not drawn: unsupported plot structure</p>'
                f'<pre style="font-size:11px;overflow:auto;max-height:400px;">{html.escape(dumped)}</pre>'
            )
        else:
            content = svg
        width_css = f"{pp['size']['width_px']}px" if "width_px" in pp["size"] else f"{pp['size']['width_pct']:.1f}%"
        panel_html = f"""<div style="display:inline-block;width:{width_css};
            vertical-align:top;padding:10px;box-sizing:border-box;">
            <h3 style="margin:0 0 5px 0;font-size:14px;color:#333;">{html.escape(str(pp['title']))}</h3>
            <div style="background:#f5f5f5;border:1px solid #ddd;border-radius:4px;
                padding:10px;min-height:{pp['size']['height_px']}px;">
                {content}
            </div>
        </div>"""
        html_panels.append(panel_html)

    safe_title = html.escape(title)
    html_doc = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>{safe_title}</title>
<style>body{{font-family:sans-serif;margin:20px;}}h1{{color:#333;}}svg{{width:100%;height:auto;}}</style></head>
<body><h1>{safe_title}</h1><div>{''.join(html_panels)}</div></body></html>"""

    logger.info("Created dashboard: %d panels, layout=%s", n, layout)

    return {
        "dashboard_data": positioned_panels,
        "html": html_doc,
        "n_panels": n,
        "layout": layout,
    }
