"""Regression tests: GWAS plot functions must create the output directory.

Defect: ``manhattan_plot``, ``qq_plot`` and ``regional_plot`` saved through a
helper that never created the parent directory, so
``output_path="nodir/m.png"`` raised ``FileNotFoundError``. Their siblings in the
same module (PCA, kinship, effect-size, missingness and enrichment plots) called
``savefig`` directly and had the same gap; ``kinship_heatmap`` swallowed the
error and returned ``None`` instead of a figure.

Every case draws a real figure with matplotlib and writes a real image file.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from metainformant.gwas.visualization import general  # noqa: E402


def _results(n: int = 60) -> list[dict[str, Any]]:
    return [
        {
            "chrom": "1" if i < n // 2 else "2",
            "pos": 1000 * (i + 1),
            "p_value": 10 ** -(1 + i % 8),
            "beta": 0.05 * ((i % 7) - 3) + 0.01,
        }
        for i in range(n)
    ]


def _pca() -> tuple[Any, ...]:
    rng = np.random.default_rng(3)
    components = [rng.standard_normal(15) for _ in range(3)]
    variance = [0.5, 0.3, 0.2]
    return (components, variance, rng.standard_normal((3, 15))), variance


def _kinship() -> np.ndarray:
    rng = np.random.default_rng(4)
    matrix = rng.random((8, 8))
    matrix = (matrix + matrix.T) / 2
    np.fill_diagonal(matrix, 1.0)
    return matrix


def _genotypes() -> dict[str, Any]:
    rng = np.random.default_rng(5)
    calls = rng.integers(0, 3, size=(12, 40))
    calls[rng.random(calls.shape) < 0.1] = -1  # some missing calls
    return {"variants": [{"id": f"rs{i}"} for i in range(40)], "genotypes": calls}


def _plot_calls(tmp_path: Path) -> dict[str, Callable[[Path], Any]]:
    gff = tmp_path / "annotation.gff3"
    gff.write_text("##gff-version 3\nchr1\tsrc\tgene\t1\t5000\t.\t+\t.\tID=gene1\n", encoding="utf-8")
    pca_result, variance = _pca()
    significant = [{"chrom": "1", "pos": 1000 * (i + 1), "p_value": 1e-9 / (i + 1)} for i in range(5)]
    weak = [{"chrom": "1", "pos": 1000 * (i + 1), "p_value": 0.01} for i in range(5)]  # none genome-wide significant
    return {
        "manhattan_plot": lambda path: general.manhattan_plot(_results(), output_path=path),
        "qq_plot": lambda path: general.qq_plot([r["p_value"] for r in _results()], output_path=path),
        "regional_plot": lambda path: general.regional_plot(_results(), "1", 0, 100000, output_path=path),
        "pca_plot": lambda path: general.pca_plot(pca_result, output_path=path, explained_var=variance),
        "kinship_heatmap": lambda path: general.kinship_heatmap(_kinship(), output_path=path),
        "effect_size_plot": lambda path: general.effect_size_plot(_results(), output_path=path),
        "missingness_plot": lambda path: general.missingness_plot(_genotypes(), output_path=path),
        "functional_enrichment_plot": lambda path: general.functional_enrichment_plot(
            significant, gff, output_path=path
        ),
        "functional_enrichment_plot_no_hits": lambda path: general.functional_enrichment_plot(
            weak, gff, output_path=path
        ),
    }


PLOT_NAMES = [
    "manhattan_plot",
    "qq_plot",
    "regional_plot",
    "pca_plot",
    "kinship_heatmap",
    "effect_size_plot",
    "missingness_plot",
    "functional_enrichment_plot",
    "functional_enrichment_plot_no_hits",
]


@pytest.mark.parametrize("plot_name", PLOT_NAMES)
def test_plot_creates_missing_parent_directories(tmp_path: Path, plot_name: str) -> None:
    """``output_path="nodir/deeper/plot.png"`` must work without the caller creating ``nodir``."""
    target = tmp_path / "nodir" / "deeper" / "plot.png"
    assert not target.parent.exists()

    figure = _plot_calls(tmp_path)[plot_name](target)

    try:
        assert figure is not None, f"{plot_name} returned None instead of a figure"
        assert target.is_file()
        assert target.stat().st_size > 0
        assert target.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    finally:
        if figure is not None:
            plt.close(figure)


@pytest.mark.parametrize("plot_name", PLOT_NAMES)
def test_plot_still_saves_into_an_existing_directory_and_accepts_str_paths(tmp_path: Path, plot_name: str) -> None:
    existing = tmp_path / "already_there"
    existing.mkdir()
    target = existing / "plot.png"

    figure = _plot_calls(tmp_path)[plot_name](str(target))

    try:
        assert figure is not None
        assert target.is_file()
    finally:
        if figure is not None:
            plt.close(figure)


def test_bare_filename_needs_no_directory(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A path with no directory part has the working directory as its parent, which already exists."""
    monkeypatch.chdir(tmp_path)

    figure = general.qq_plot([r["p_value"] for r in _results()], output_path="qq.png")

    try:
        assert (tmp_path / "qq.png").is_file()
    finally:
        plt.close(figure)


def test_unwritable_parent_is_reported_not_swallowed(tmp_path: Path) -> None:
    """A parent that is a regular file cannot become a directory; that must stay an error."""
    blocker = tmp_path / "blocker"
    blocker.write_text("not a directory", encoding="utf-8")
    open_before = set(plt.get_fignums())

    try:
        with pytest.raises(OSError):
            general.qq_plot([r["p_value"] for r in _results()], output_path=blocker / "qq.png")
    finally:
        for number in set(plt.get_fignums()) - open_before:  # the figure drawn before the failed save
            plt.close(number)
