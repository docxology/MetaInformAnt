"""Interactive dashboard data structures for web-based visualization.

Generates JSON-serializable plot data structures for interactive scatter
plots, heatmaps, genome browser tracks, volcano plots, and composed
dashboards. Includes standalone HTML export using embedded SVG and
vanilla JavaScript (no external dependencies).

All functions produce pure Python dictionaries ready for JSON serialization
and rendering by any web framework or standalone HTML viewer.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


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


def _hierarchical_cluster_order(
    data: list[list[float]],
    axis: int = 0,
) -> list[int]:
    """Simple hierarchical clustering for row/column ordering.

    Uses a greedy nearest-neighbour approach for ordering.
    """
    if axis == 1:
        # Transpose for column clustering
        n_cols = len(data[0]) if data else 0
        vectors = [[data[i][j] for i in range(len(data))] for j in range(n_cols)]
    else:
        vectors = data

    n = len(vectors)
    if n <= 1:
        return list(range(n))

    remaining = list(range(n))
    order = [remaining.pop(0)]

    while remaining:
        last = order[-1]
        best_idx = 0
        best_dist = float("inf")
        for idx, r in enumerate(remaining):
            d = sum((vectors[last][k] - vectors[r][k]) ** 2 for k in range(len(vectors[0])))
            if d < best_dist:
                best_dist = d
                best_idx = idx
        order.append(remaining.pop(best_idx))

    return order


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
        cluster_rows: Whether to cluster rows.
        cluster_cols: Whether to cluster columns.

    Returns:
        Dictionary with keys:
            - plot_data: Heatmap cell data with row, col, value.
            - row_order: Ordered row indices.
            - col_order: Ordered column indices.
            - dendrograms: Row/col dendrogram data (structure for rendering).
            - color_range: (min_value, max_value).
    """
    n_rows = len(data)
    n_cols = len(data[0]) if data else 0

    if len(row_names) != n_rows:
        raise ValueError(f"row_names ({len(row_names)}) != data rows ({n_rows})")
    if len(col_names) != n_cols:
        raise ValueError(f"col_names ({len(col_names)}) != data cols ({n_cols})")

    row_order = _hierarchical_cluster_order(data, axis=0) if cluster_rows else list(range(n_rows))
    col_order = _hierarchical_cluster_order(data, axis=1) if cluster_cols else list(range(n_cols))

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
        "rows": {"order": row_order, "clustered": cluster_rows},
        "cols": {"order": col_order, "clustered": cluster_cols},
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
        features: List of feature dicts, each with at least ``start``,
            ``end``, and ``name`` keys. Optional: ``strand``, ``score``,
            ``chromosome``, ``color``.
        track_type: Track type: ``"gene"``, ``"variant"``, ``"peak"``,
            ``"coverage"``.

    Returns:
        Dictionary with keys:
            - track_data: Processed feature list.
            - track_type: Track type.
            - coordinate_range: (min_start, max_end).
            - n_features: Number of features.
    """
    processed = []
    min_start = float("inf")
    max_end = float("-inf")

    for feat in features:
        start = feat.get("start", 0)
        end = feat.get("end", start + 1)
        name = feat.get("name", "unnamed")

        min_start = min(min_start, start)
        max_end = max(max_end, end)

        entry = {
            "start": start,
            "end": end,
            "name": name,
            "strand": feat.get("strand", "+"),
            "score": feat.get("score", 0),
            "chromosome": feat.get("chromosome", "chr1"),
            "color": feat.get("color", _track_color(track_type)),
            "width": end - start,
        }
        processed.append(entry)

    if not processed:
        min_start = 0
        max_end = 1000

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
    """
    points = []
    counts = {"up": 0, "down": 0, "not_significant": 0}

    for result in de_results:
        gene = result.get("gene", result.get("name", "unnamed"))
        fc = result.get("log2_fold_change", result.get("log2fc", 0.0))
        pval = result.get("p_value", result.get("pvalue", result.get("adjusted_p", 1.0)))

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

    Generates a self-contained HTML file with embedded SVG graphics and
    vanilla JavaScript for interactivity (hover tooltips, zoom).

    Args:
        plot_data: Plot data structure from any ``create_interactive_*`` function.
        output_path: Output file path.
        title: HTML page title.

    Returns:
        Absolute path of the generated HTML file.
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    data_json = json.dumps(plot_data, default=str)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
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
        <h1>{title}</h1>
        <div class="plot-area" id="plot"></div>
        <div class="stats" id="stats"></div>
    </div>
    <div class="tooltip" id="tooltip"></div>
    <script>
        const plotData = {data_json};
        const container = document.getElementById('plot');
        const tooltip = document.getElementById('tooltip');
        const stats = document.getElementById('stats');

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
                    const color = p.color || '#1f77b4';
                    const label = p.gene || p.label || ('Point ' + i);
                    svg += '<circle class="point" cx="'+cx+'" cy="'+cy+'" r="4" fill="'+color+'"' +
                           ' data-label="'+label+'" data-x="'+p.x.toFixed(3)+'" data-y="'+p.y.toFixed(3)+'"/>';
                }});

                // Axes
                svg += '<line x1="'+m.l+'" y1="'+(m.t+ph)+'" x2="'+(m.l+pw)+'" y2="'+(m.t+ph)+'" stroke="#333" stroke-width="1"/>';
                svg += '<line x1="'+m.l+'" y1="'+m.t+'" x2="'+m.l+'" y2="'+(m.t+ph)+'" stroke="#333" stroke-width="1"/>';
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

    out.write_text(html_content)
    abs_path = str(out.resolve())
    logger.info("Exported interactive plot to %s", abs_path)
    return abs_path


# ---------------------------------------------------------------------------
# Dashboard Composition
# ---------------------------------------------------------------------------


def create_dashboard(
    panels: list[dict],
    layout: str = "grid",
    title: str = "Dashboard",
) -> dict:
    """Compose multiple plots into a dashboard layout.

    Args:
        panels: List of panel dicts, each with ``title`` and ``plot_data``
            keys (from any ``create_interactive_*`` function). Optional
            ``width`` and ``height`` keys.
        layout: Layout strategy: ``"grid"`` (auto-arrange), ``"rows"``
            (one per row), or ``"columns"`` (side-by-side).
        title: Dashboard title.

    Returns:
        Dictionary with keys:
            - dashboard_data: Structured panel data with positions.
            - html: Standalone HTML string for the dashboard.
            - n_panels: Number of panels.
            - layout: Layout used.
    """
    n = len(panels)

    # Compute positions
    if layout == "grid":
        cols = math.ceil(math.sqrt(n))
        rows = math.ceil(n / cols)
    elif layout == "rows":
        cols = 1
        rows = n
    elif layout == "columns":
        cols = n
        rows = 1
    else:
        cols = math.ceil(math.sqrt(n))
        rows = math.ceil(n / cols)

    panel_width = 100.0 / cols
    panel_height = max(300, 600 // rows)

    positioned_panels = []
    for idx, panel in enumerate(panels):
        row = idx // cols
        col = idx % cols
        positioned_panels.append(
            {
                "title": panel.get("title", f"Panel {idx + 1}"),
                "plot_data": panel.get("plot_data", {}),
                "position": {"row": row, "col": col},
                "size": {
                    "width_pct": panel_width,
                    "height_px": panel.get("height", panel_height),
                },
            }
        )

    # Generate simple HTML dashboard
    html_panels = []
    for pp in positioned_panels:
        panel_html = f"""<div style="display:inline-block;width:{pp['size']['width_pct']:.1f}%;
            vertical-align:top;padding:10px;box-sizing:border-box;">
            <h3 style="margin:0 0 5px 0;font-size:14px;color:#333;">{pp['title']}</h3>
            <div style="background:#f5f5f5;border:1px solid #ddd;border-radius:4px;
                padding:10px;min-height:{pp['size']['height_px']}px;">
                <pre style="font-size:11px;overflow:auto;max-height:400px;">{json.dumps(pp['plot_data'], indent=1, default=str)[:2000]}</pre>
            </div>
        </div>"""
        html_panels.append(panel_html)

    html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>{title}</title>
<style>body{{font-family:sans-serif;margin:20px;}}h1{{color:#333;}}</style></head>
<body><h1>{title}</h1><div>{''.join(html_panels)}</div></body></html>"""

    logger.info("Created dashboard: %d panels, layout=%s", n, layout)

    return {
        "dashboard_data": positioned_panels,
        "html": html,
        "n_panels": n,
        "layout": layout,
    }
