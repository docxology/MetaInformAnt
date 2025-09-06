import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402

from metainformant.visualization import animate_time_series, heatmap, lineplot


def test_lineplot_and_heatmap_render():
    ax1 = lineplot(None, [0.1, 0.2, 0.3], label="a")
    assert ax1 is not None
    plt.close(ax1.figure)

    ax2 = heatmap([[1, 0], [0, 1]])
    assert ax2 is not None
    plt.close(ax2.figure)


def test_animation_constructs():
    fig, anim = animate_time_series([0, 1, 0], interval_ms=50)
    assert fig is not None
    assert anim is not None
    plt.close(fig)


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
