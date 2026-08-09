"""Regression tests for cross-species divergence figures."""

from __future__ import annotations

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
    plt.close("all")
