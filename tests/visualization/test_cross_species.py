"""Regression tests for cross-species divergence figures."""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from metainformant.visualization.plots import cross_species


def _divergence_matrix() -> pd.DataFrame:
    """Return a small matrix containing the full 1 - correlation range."""
    return pd.DataFrame(
        [[0.0, 0.4, 1.8], [0.4, 0.0, 1.2], [1.8, 1.2, 0.0]],
        index=["sp_a", "sp_b", "sp_c"],
        columns=["sp_a", "sp_b", "sp_c"],
    )


def test_divergence_heatmap_uses_full_correlation_divergence_scale(tmp_path: Path, monkeypatch) -> None:
    """Heatmaps must not clip valid divergences above 0.75."""
    captured: dict[str, object] = {}
    original_heatmap = cross_species.sns.heatmap

    def spy_heatmap(*args, **kwargs):
        captured.update(kwargs)
        return original_heatmap(*args, **kwargs)

    monkeypatch.setattr(cross_species.sns, "heatmap", spy_heatmap)
    output_path = tmp_path / "divergence_heatmap.png"

    cross_species.plot_divergence_heatmap(_divergence_matrix(), output_path)

    assert output_path.exists()
    assert captured["vmin"] == 0.0
    assert captured["vmax"] == 2.0
    assert np.array_equal(captured["cbar_kws"]["ticks"], np.linspace(0.0, 2.0, 5))
    mask = np.asarray(captured["mask"])
    assert mask.shape == (3, 3)
    assert np.array_equal(np.diag(mask), np.ones(3, dtype=bool))
    plt.close("all")


def test_cividis_annotations_select_readable_contrast() -> None:
    """Numeric cells must remain legible at both ends of the color scale."""

    assert cross_species._cividis_annotation_color(0.0) == "#FFFFFF"
    assert cross_species._cividis_annotation_color(2.0) == "#111111"

    for value in np.linspace(0.0, 2.0, 41):
        position = value / 2.0
        background = cross_species.matplotlib.colormaps["cividis"](position)[:3]
        foreground_hex = cross_species._cividis_annotation_color(float(value))
        foreground = tuple(int(foreground_hex[index : index + 2], 16) / 255 for index in (1, 3, 5))
        background_luminance = cross_species._relative_luminance(background)
        foreground_luminance = cross_species._relative_luminance(foreground)
        contrast = (max(background_luminance, foreground_luminance) + 0.05) / (
            min(background_luminance, foreground_luminance) + 0.05
        )
        assert contrast >= 4.5


def test_combined_summary_uses_full_correlation_divergence_scale(tmp_path: Path, monkeypatch) -> None:
    """The combined manuscript figure must use the same fixed scale."""
    captured: dict[str, object] = {}
    original_heatmap = cross_species.sns.heatmap

    def spy_heatmap(*args, **kwargs):
        captured.update(kwargs)
        return original_heatmap(*args, **kwargs)

    monkeypatch.setattr(cross_species.sns, "heatmap", spy_heatmap)
    output_path = tmp_path / "combined_summary.png"

    cross_species.plot_combined_summary(_divergence_matrix(), output_path)

    assert output_path.exists()
    assert captured["vmin"] == 0.0
    assert captured["vmax"] == 2.0
    plt.close("all")


def test_cross_species_clustering_uses_average_linkage(tmp_path: Path, monkeypatch) -> None:
    """Correlation-derived dissimilarities must not be passed to Ward linkage."""

    methods: list[str] = []
    original_linkage = cross_species.linkage

    def spy_linkage(values, method, *args, **kwargs):
        methods.append(method)
        return original_linkage(values, method=method, *args, **kwargs)

    monkeypatch.setattr(cross_species, "linkage", spy_linkage)
    cross_species.plot_divergence_heatmap(_divergence_matrix(), tmp_path / "heatmap.png")
    cross_species.plot_dendrogram(_divergence_matrix(), tmp_path / "dendrogram.png")

    assert methods == ["average", "average"]
    plt.close("all")


def test_cross_species_plot_rejects_incomplete_divergence_matrix(
    tmp_path: Path,
) -> None:
    """Missing pairwise distances must be resolved before figure generation."""

    matrix = _divergence_matrix().copy()
    matrix.loc["sp_a", "sp_b"] = np.nan
    matrix.loc["sp_b", "sp_a"] = np.nan

    with pytest.raises(ValueError, match="non-finite"):
        cross_species.plot_dendrogram(matrix, tmp_path / "dendrogram.png")


def test_profile_quality_plot_uses_source_table_columns(tmp_path: Path) -> None:
    """Profile quality figures render from explicit validity counts."""

    quality = pd.DataFrame(
        {
            "species": ["sp_a", "sp_b"],
            "positive_features": [90, 80],
            "zero_features": [10, 20],
            "nonfinite_features": [0, 1],
        }
    )
    output = tmp_path / "profile_quality.png"
    cross_species.plot_profile_quality(quality, output)
    assert output.is_file()
    plt.close("all")


def test_divergence_stability_plot_has_fixed_descriptive_scale(tmp_path: Path) -> None:
    """Sensitivity intervals use the same bounded descriptive distance scale."""

    stability = pd.DataFrame(
        {
            "species_a": ["sp_a", "sp_a"],
            "species_b": ["sp_b", "sp_c"],
            "point_estimate": [0.4, 1.1],
            "sensitivity_lower": [0.2, 0.8],
            "sensitivity_upper": [0.7, 1.4],
            "sensitivity_iqr": [0.5, 0.6],
        }
    )
    output = tmp_path / "divergence_stability.png"
    cross_species.plot_divergence_stability(stability, output)
    assert output.is_file()


def _capture_saved_figure(monkeypatch) -> dict[str, object]:
    """Capture the figure passed to the deterministic saver for text checks."""
    captured: dict[str, object] = {}
    original_save = cross_species.save_figure_deterministic

    def spy_save(fig, path, **kwargs):
        captured["fig"] = fig
        return original_save(fig, path, **kwargs)

    monkeypatch.setattr(cross_species, "save_figure_deterministic", spy_save)
    return captured


def _figure_texts(fig) -> list[str]:
    """Collect every rendered string: titles, annotations, ticks, legends."""
    texts = [text.get_text() for text in fig.texts]
    for ax in fig.axes:
        for loc in ("left", "center", "right"):
            texts.append(ax.get_title(loc))
        texts.extend(text.get_text() for text in ax.texts)
        texts.extend(label.get_text() for label in ax.get_xticklabels())
        texts.extend(label.get_text() for label in ax.get_yticklabels())
        legend = ax.get_legend()
        if legend:
            texts.extend(text.get_text() for text in legend.get_texts())
    return texts


def _assert_no_inferential_language(texts: list[str]) -> None:
    """Native pairwise figures must not carry stars or p-value claims."""
    blob = "\n".join(texts)
    assert not re.search(r"\bp\s*[<=]\s*\d", blob, flags=re.IGNORECASE), blob
    assert "***" not in blob, blob
    for text in texts:
        for sentence in re.split(r"[.;\n]", text):
            if "confidence interval" in sentence.lower():
                lowered = sentence.lower()
                assert "not confidence interval" in lowered or "no " in lowered, sentence


def test_divergence_heatmap_annotates_descriptive_summary_with_denominator(
    tmp_path: Path,
    monkeypatch,
) -> None:
    """The heatmap must label itself descriptive and carry a data-derived denominator."""
    captured = _capture_saved_figure(monkeypatch)
    cross_species.plot_divergence_heatmap(_divergence_matrix(), tmp_path / "heatmap.png")

    fig = captured["fig"]
    texts = _figure_texts(fig)
    blob = "\n".join(texts)
    assert "n=3 species in plotted matrix" in blob
    assert "Descriptive summary of native pairwise expression divergence" in blob
    assert "no significance tests, p-values, or confidence intervals" in blob
    _assert_no_inferential_language(texts)
    # Cell annotations are numeric only: no significance stars are possible.
    for ax in fig.axes:
        for annotation in ax.texts:
            assert "*" not in annotation.get_text()
    plt.close("all")


def test_divergence_heatmap_denominator_follows_plotted_data(tmp_path: Path, monkeypatch) -> None:
    """The denominator must be computed from the matrix, not hard-coded."""
    captured = _capture_saved_figure(monkeypatch)
    matrix = _divergence_matrix().iloc[:2, :2]

    cross_species.plot_divergence_heatmap(matrix, tmp_path / "heatmap_two_species.png")

    blob = "\n".join(_figure_texts(captured["fig"]))
    assert "n=2 species in plotted matrix" in blob
    assert "n=3" not in blob
    plt.close("all")


def test_dendrogram_declares_expression_profile_clustering(tmp_path: Path, monkeypatch) -> None:
    """Dendrogram titles/labels must disclaim a species-tree reading."""
    captured = _capture_saved_figure(monkeypatch)
    cross_species.plot_dendrogram(_divergence_matrix(), tmp_path / "dendrogram.png")

    fig = captured["fig"]
    texts = _figure_texts(fig)
    blob = "\n".join(texts)
    assert "not a species tree or phylogeny" in blob
    assert "expression profiles" in blob
    assert "n=3 species in plotted matrix" in blob
    _assert_no_inferential_language(texts)
    plt.close("all")


def test_divergence_stability_reports_descriptive_sensitivity_and_denominator(
    tmp_path: Path,
    monkeypatch,
) -> None:
    """Stability intervals stay sensitivity diagnostics with a species denominator."""
    captured = _capture_saved_figure(monkeypatch)
    stability = pd.DataFrame(
        {
            "species_a": ["sp_a", "sp_a"],
            "species_b": ["sp_b", "sp_c"],
            "point_estimate": [0.4, 1.1],
            "sensitivity_lower": [0.2, 0.8],
            "sensitivity_upper": [0.7, 1.4],
            "sensitivity_iqr": [0.5, 0.6],
        }
    )

    cross_species.plot_divergence_stability(stability, tmp_path / "stability.png")

    texts = _figure_texts(captured["fig"])
    blob = "\n".join(texts)
    assert "not confidence intervals" in blob
    assert "no p-values or significance stars" in blob
    assert "n=3 species in plotted pairs" in blob
    _assert_no_inferential_language(texts)
    plt.close("all")


def test_combined_summary_states_denominator_and_not_a_species_tree(
    tmp_path: Path,
    monkeypatch,
) -> None:
    """The manuscript panel figure carries the denominator and the disclaimer."""
    captured = _capture_saved_figure(monkeypatch)
    cross_species.plot_combined_summary(_divergence_matrix(), tmp_path / "combined.png")

    fig = captured["fig"]
    texts = _figure_texts(fig)
    blob = "\n".join(texts)
    assert "n=3 species in plotted matrix" in blob
    assert "not a species tree" in blob
    _assert_no_inferential_language(texts)
    plt.close("all")


def test_species_level_figures_report_data_derived_denominators(
    tmp_path: Path,
    monkeypatch,
) -> None:
    """Coverage and pair-extreme figures derive denominators from plotted rows."""
    captured = _capture_saved_figure(monkeypatch)
    coverage = pd.Series([10, 20, 30], index=["sp_a", "sp_b", "sp_c"])

    cross_species.plot_coverage(coverage, total_groups=40, output_path=tmp_path / "coverage.png")
    blob = "\n".join(_figure_texts(captured["fig"]))
    assert "n=3 species shown" in blob

    cross_species.plot_top_pairs(_divergence_matrix(), tmp_path / "top_pairs.png")
    blob = "\n".join(_figure_texts(captured["fig"]))
    assert "n=3 species in plotted matrix" in blob
    _assert_no_inferential_language(_figure_texts(captured["fig"]))
    plt.close("all")
