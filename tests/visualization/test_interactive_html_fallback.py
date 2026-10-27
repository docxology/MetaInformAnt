"""Without Plotly, an HTML ``output_path`` must fail clearly, before any work is done.

The interactive plot wrappers fall back to matplotlib when Plotly is missing, and the fallback saved the figure with
``fig.savefig(output_path)``. matplotlib has no HTML writer, so asking for ``fig.html`` raised
``ValueError: Format 'html' is not supported`` from deep inside ``savefig``, after a figure and the output directory
had already been created. Saving a PNG under a different name would silently hand the caller something other than the
interactive file they asked for, so the wrappers now raise an ``ImportError`` that names Plotly and leaves no figure,
file or directory behind; static image paths keep working.

The fallback is only reachable without Plotly. The in-process tests skip when Plotly is installed, and a child
interpreter that cannot import Plotly runs the same checks in every environment. All data are real frames drawn
onto real matplotlib Agg figures; no test doubles (real-implementation policy).
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable, Iterator

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes

from metainformant.visualization.dashboards import interactive

matplotlib.use("Agg")

without_plotly = pytest.mark.skipif(
    interactive.HAS_PLOTLY, reason="the matplotlib fallback is only used without Plotly"
)

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
MARKER = "HTML-FALLBACK-REPORT:"


def _scatter(path: Path | None) -> Any:
    frame = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [4.0, 5.0, 6.0]})
    return interactive.interactive_scatter(frame, "x", "y", output_path=path)


def _heatmap(path: Path | None) -> Any:
    return interactive.interactive_heatmap(np.arange(9.0).reshape(3, 3), output_path=path)


def _volcano(path: Path | None) -> Any:
    rng = np.random.default_rng(1)
    frame = pd.DataFrame(
        {
            "log2FoldChange": rng.uniform(-3, 3, 40),
            "pvalue": rng.uniform(1e-4, 1, 40),
            "gene": [f"G{i}" for i in range(40)],
        }
    )
    return interactive.interactive_volcano(frame, output_path=path)


def _manhattan(path: Path | None) -> Any:
    frame = pd.DataFrame(
        {
            "chromosome": np.repeat(["1", "2", "3"], 10),
            "position": np.tile(np.arange(1000, 11000, 1000), 3),
            "pvalue": np.linspace(1e-6, 0.9, 30),
        }
    )
    return interactive.interactive_manhattan(frame, output_path=path)


PLOTS: dict[str, Callable[[Path | None], Any]] = {
    "scatter": _scatter,
    "heatmap": _heatmap,
    "volcano": _volcano,
    "manhattan": _manhattan,
}


@pytest.fixture(autouse=True)
def _close_figures() -> Iterator[None]:
    yield
    plt.close("all")


@without_plotly
@pytest.mark.parametrize("name", sorted(PLOTS))
@pytest.mark.parametrize("filename", ["fig.html", "fig.htm", "FIG.HTML"])
def test_html_output_raises_an_import_error_that_names_plotly(tmp_path: Path, name: str, filename: str) -> None:
    target = tmp_path / filename

    with pytest.raises(ImportError) as excinfo:
        PLOTS[name](target)

    message = str(excinfo.value)
    assert "Plotly" in message
    assert str(target) in message
    assert "uv pip install plotly" in message
    assert ".png" in message  # says what to pass instead


@without_plotly
@pytest.mark.parametrize("name", sorted(PLOTS))
def test_a_refused_html_export_leaves_no_figure_file_or_directory(tmp_path: Path, name: str) -> None:
    figures_before = plt.get_fignums()
    target = tmp_path / "exports" / "fig.html"

    with pytest.raises(ImportError):
        PLOTS[name](target)

    assert plt.get_fignums() == figures_before  # no figure was created and abandoned
    assert not (tmp_path / "exports").exists()  # no output directory was created
    assert list(tmp_path.iterdir()) == []


@without_plotly
@pytest.mark.parametrize("name", sorted(PLOTS))
@pytest.mark.parametrize("suffix", [".png", ".svg"])
def test_static_image_export_still_works(tmp_path: Path, name: str, suffix: str) -> None:
    target = tmp_path / "nested" / f"fig{suffix}"

    result = PLOTS[name](target)

    assert isinstance(result, Axes)
    data = target.read_bytes()
    if suffix == ".png":
        assert data.startswith(PNG_SIGNATURE)
    else:
        assert b"<svg" in data


@without_plotly
@pytest.mark.parametrize("name", sorted(PLOTS))
def test_no_output_path_returns_the_axes_and_writes_nothing(tmp_path: Path, name: str) -> None:
    assert isinstance(PLOTS[name](None), Axes)
    assert list(tmp_path.iterdir()) == []


CHILD_SCRIPT = """
import json
import sys
from pathlib import Path

sys.modules["plotly"] = None  # a None entry makes ``import plotly`` raise ImportError, as if it were not installed

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from metainformant.visualization.dashboards import interactive

assert interactive.HAS_PLOTLY is False

workdir = Path(sys.argv[1])
frame = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [4.0, 5.0, 6.0]})
volcano = pd.DataFrame({"log2FoldChange": [-2.0, 0.1, 2.5, 1.2], "pvalue": [0.001, 0.5, 0.0004, 0.2]})
manhattan = pd.DataFrame(
    {"chromosome": ["1", "1", "2", "2"], "position": [10, 20, 10, 20], "pvalue": [0.5, 0.001, 0.2, 1e-6]}
)
plots = {
    "scatter": lambda path: interactive.interactive_scatter(frame, "x", "y", output_path=path),
    "heatmap": lambda path: interactive.interactive_heatmap(np.eye(3), output_path=path),
    "volcano": lambda path: interactive.interactive_volcano(volcano, gene_col=None, output_path=path),
    "manhattan": lambda path: interactive.interactive_manhattan(manhattan, output_path=path),
}

report = {}
for name, plot in plots.items():
    figures_before = len(plt.get_fignums())
    html_path = workdir / name / "fig.html"
    try:
        plot(html_path)
        outcome = ["no error", ""]
    except Exception as exc:
        outcome = [type(exc).__name__, str(exc)]
    png_path = workdir / name / "fig.png"
    plot(png_path)
    report[name] = {
        "html": outcome,
        "figures_leaked_by_html": len(plt.get_fignums()) - figures_before - 1,  # the PNG export keeps one open
        "png_is_png": png_path.read_bytes()[:8] == b"\\x89PNG\\r\\n\\x1a\\n",
        "html_file_exists": html_path.exists(),
    }
    plt.close("all")

print(sys.argv[2] + json.dumps(report))
"""


@pytest.fixture(scope="module")
def child_report(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """Run the four wrappers in an interpreter that cannot import Plotly, whether or not it is installed."""
    workdir = tmp_path_factory.mktemp("no_plotly")
    proc = subprocess.run(
        [sys.executable, "-c", CHILD_SCRIPT, str(workdir), MARKER],
        capture_output=True,
        text=True,
        timeout=180,
        cwd=workdir,
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    line = next(line for line in proc.stdout.splitlines() if line.startswith(MARKER))
    report: dict[str, Any] = json.loads(line[len(MARKER) :])
    return report


@pytest.mark.parametrize("name", sorted(PLOTS))
def test_without_plotly_html_export_is_refused_with_import_error(child_report: dict[str, Any], name: str) -> None:
    error_type, message = child_report[name]["html"]

    assert error_type == "ImportError", f"{name}: {error_type}: {message}"
    assert "Plotly" in message
    assert child_report[name]["html_file_exists"] is False
    assert child_report[name]["figures_leaked_by_html"] == 0


@pytest.mark.parametrize("name", sorted(PLOTS))
def test_without_plotly_static_images_are_still_written(child_report: dict[str, Any], name: str) -> None:
    assert child_report[name]["png_is_png"] is True
