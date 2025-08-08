from __future__ import annotations

import matplotlib


def test_lineplot_basic():
    matplotlib.use("Agg")  # non-interactive backend
    from metainformant.visualization import lineplot

    ax = lineplot(None, [1, 2, 3], label="y")
    assert ax.get_lines(), "Expected at least one line"


def test_heatmap_basic():
    matplotlib.use("Agg")
    from metainformant.visualization import heatmap

    ax = heatmap([[1, 0], [0, 1]])
    # Heatmap draws a QuadMesh
    assert any(hasattr(coll, "get_array") for coll in ax.collections)


def test_animation_builds():
    matplotlib.use("Agg")
    from metainformant.visualization import animate_time_series

    fig, anim = animate_time_series([0, 1, 0, 1], interval_ms=10)
    assert fig is not None and anim is not None


