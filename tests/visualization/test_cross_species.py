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


def _family_map() -> dict[str, str]:
    return {"sp_a": "Fam1", "sp_b": "Fam1", "sp_c": "Fam2"}


def _family_colors() -> dict[str, str]:
    return {"Fam1": "#0072B2", "Fam2": "#D55E00"}


def _stability_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "species_a": ["sp_a", "sp_a"],
            "species_b": ["sp_b", "sp_c"],
            "point_estimate": [0.4, 1.1],
            "sensitivity_lower": [0.2, 0.8],
            "sensitivity_upper": [0.7, 1.4],
            "sensitivity_iqr": [0.5, 0.6],
        }
    )


def _species_summary_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "species": ["sp_a", "sp_b", "sp_c"],
            "total_features": [100, 100, 100],
            "expressed_features": [90, 80, 70],
            "mean_expression": [1.5, 2.5, 0.5],
        }
    )


def _profile_quality_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "species": ["sp_a", "sp_b", "sp_c"],
            "positive_features": [90, 80, 70],
            "zero_features": [10, 20, 30],
            "nonfinite_features": [0, 0, 0],
        }
    )


def test_all_cross_species_figures_close_after_saving(tmp_path: Path) -> None:
    """Every figure-producing helper must close its figure; none may leak."""
    matrix = _divergence_matrix()
    coverage = pd.Series([10, 20, 30], index=["sp_a", "sp_b", "sp_c"])
    calls = [
        lambda p: cross_species.plot_divergence_heatmap(matrix, p),
        lambda p: cross_species.plot_dendrogram(matrix, p),
        lambda p: cross_species.plot_coverage(coverage, total_groups=40, output_path=p),
        lambda p: cross_species.plot_top_pairs(matrix, p),
        lambda p: cross_species.plot_family_violin(matrix, p, _family_map()),
        lambda p: cross_species.plot_method_comparison(matrix, matrix, p),
        lambda p: cross_species.plot_mean_divergence_rank(matrix, p, _family_map(), _family_colors()),
        lambda p: cross_species.plot_species_summary(_species_summary_table(), p),
        lambda p: cross_species.plot_profile_quality(_profile_quality_table(), p),
        lambda p: cross_species.plot_divergence_stability(_stability_table(), p),
        lambda p: cross_species.plot_combined_summary(matrix, p),
    ]
    for index, plot_call in enumerate(calls):
        plot_call(tmp_path / f"figure_{index}.png")
        assert plt.get_fignums() == [], f"figure leaked by call index {index}"


def test_species_summary_requires_all_source_columns(tmp_path: Path) -> None:
    """Missing mean_expression/species must raise a clear ValueError, not KeyError."""
    incomplete = _species_summary_table().drop(columns=["mean_expression"])

    with pytest.raises(ValueError, match="requires columns"):
        cross_species.plot_species_summary(incomplete, tmp_path / "summary.png")


def test_species_summary_and_mean_rank_render(tmp_path: Path) -> None:
    """Untested public figure helpers must render their output files."""
    summary_path = tmp_path / "summary.png"
    cross_species.plot_species_summary(_species_summary_table(), summary_path)
    assert summary_path.is_file()

    rank_path = tmp_path / "mean_rank.png"
    cross_species.plot_mean_divergence_rank(_divergence_matrix(), rank_path, _family_map(), _family_colors())
    assert rank_path.is_file()

    comparison_path = tmp_path / "method_comparison.png"
    cross_species.plot_method_comparison(_divergence_matrix(), _divergence_matrix(), comparison_path)
    assert comparison_path.is_file()
    plt.close("all")


def test_family_violin_renders_and_skips_empty(tmp_path: Path) -> None:
    """Family violin renders with enough pairs and skips silently without any."""
    output = tmp_path / "family_violin.png"
    cross_species.plot_family_violin(_divergence_matrix(), output, _family_map())
    assert output.is_file()

    empty_output = tmp_path / "family_violin_empty.png"
    cross_species.plot_family_violin(_divergence_matrix().iloc[:1, :1], empty_output, _family_map())
    assert not empty_output.exists()
    plt.close("all")


def test_validated_condensed_rejects_invalid_matrices() -> None:
    """Structural and numerical violations must each be rejected explicitly."""
    with pytest.raises(ValueError, match="row and column labels must match"):
        cross_species._validated_condensed(pd.DataFrame([[0.0, 0.5], [0.5, 0.0]], index=["a", "b"], columns=["a", "c"]))
    with pytest.raises(ValueError, match="symmetric"):
        cross_species._validated_condensed(pd.DataFrame([[0.0, 0.5], [0.7, 0.0]], index=["a", "b"], columns=["a", "b"]))
    with pytest.raises(ValueError, match="diagonal must be zero"):
        cross_species._validated_condensed(pd.DataFrame([[0.1, 0.5], [0.5, 0.0]], index=["a", "b"], columns=["a", "b"]))
    with pytest.raises(ValueError, match="0--2 range"):
        cross_species._validated_condensed(pd.DataFrame([[0.0, 2.5], [2.5, 0.0]], index=["a", "b"], columns=["a", "b"]))
    with pytest.raises(ValueError, match="at least two species"):
        cross_species._validated_condensed(pd.DataFrame([[0.0]], index=["a"], columns=["a"]))


def test_get_family_color_defaults() -> None:
    """Unmapped species/families fall back to sensible defaults."""
    assert cross_species._get_family_color("sp_a") == "#34495e"
    assert cross_species._get_family_color("sp_unknown", _family_map(), _family_colors()) == "#95a5a6"
    assert cross_species._get_family_color("sp_a", _family_map(), _family_colors()) == "#0072B2"
