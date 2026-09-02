"""Visualization MCP tools: chart-convention rendering to explicit output dirs."""

from __future__ import annotations

from typing import Any

from metainformant.mcp.tools._spec import validate_output_dir


def _handle_chart(
    kind: str,
    x: list[float],
    y: list[float] | None,
    output_dir: str,
    title: str = "",
    x_label: str = "",
    y_label: str = "",
) -> dict:
    """Render a basic chart (line|scatter|bar) with repo chart conventions."""
    import numpy as np

    from metainformant.visualization.plots import basic

    out_dir = validate_output_dir(output_dir)
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float) if y is not None else None
    out_path = out_dir / f"chart_{kind}.png"
    if kind == "line":
        basic.lineplot(x_arr, y_arr, output_path=out_path)
    elif kind == "scatter":
        if y_arr is None:
            return {"error": "scatter requires y"}
        basic.scatter_plot(x_arr, y_arr, output_path=out_path)
    elif kind == "bar":
        if y_arr is None:
            return {"error": "bar requires y"}
        basic.bar_plot(x_arr, y_arr, output_path=out_path)
    else:
        return {"error": f"unknown kind: {kind} (use line|scatter|bar)"}
    return {"output": str(out_path), "kind": kind, "n_points": int(len(x_arr))}


TOOL_SPEC: dict[str, Any] = {
    "name": "viz_chart",
    "description": "Render a line/scatter/bar chart with repository chart conventions to an explicit output dir.",
    "input_schema": {
        "type": "object",
        "properties": {
            "kind": {"type": "string", "enum": ["line", "scatter", "bar"]},
            "x": {"type": "array", "items": {"type": "number"}},
            "y": {"type": "array", "items": {"type": "number"}},
            "output_dir": {"type": "string"},
            "title": {"type": "string"},
            "x_label": {"type": "string"},
            "y_label": {"type": "string"},
        },
        "required": ["kind", "x", "output_dir"],
    },
    "handler": _handle_chart,
    "writes": "output-dir-only",
}

ALL_SPECS: list[dict[str, Any]] = [TOOL_SPEC]
