"""Tests for atlas-style visualization module (descriptive-only)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from PIL import Image

from metainformant.rna.analysis.atlas_plots import (
    OKABE_ITO,
    compute_tau,
    plot_orthogroup_small_multiples,
    plot_tau_heatmap,
    plot_tau_orthology_strips,
)

RNG_SEED = 42


def _synthetic_tau_frame() -> pd.DataFrame:
    rng = np.random.default_rng(RNG_SEED)
    species = ["apis_mellifera", "bombus_terrestris", "ceratina_calcarata"]
    tissues = ["brain", "ovary", "midgut", "antenna"]
    data = rng.uniform(0.0, 1.0, size=(len(species), len(tissues)))
    return pd.DataFrame(data, index=species, columns=tissues)


def _synthetic_long_expression() -> pd.DataFrame:
    rng = np.random.default_rng(RNG_SEED)
    rows = []
    tissues = ["brain", "ovary", "midgut"]
    for og in ["OG001", "OG002", "OG003", "OG004"]:
        base = rng.uniform(10.0, 1000.0)
        for sp in ["apis_mellifera", "bombus_terrestris"]:
            for t in tissues:
                rows.append(
                    {
                        "orthogroup": og,
                        "species": sp,
                        "tissue": t,
                        "expression": float(base * rng.uniform(0.2, 5.0)),
                    }
                )
    return pd.DataFrame(rows)


def _synthetic_tau_classes(n_per_class: int = 25) -> pd.DataFrame:
    rng = np.random.default_rng(RNG_SEED)
    frames = []
    class_means = {
        "many_to_many": 0.2,
        "many_to_one": 0.4,
        "one_to_many": 0.6,
        "one_to_one": 0.8,
    }
    for cls, mean in class_means.items():
        values = np.clip(rng.normal(mean, 0.15, size=n_per_class), 0.0, 1.0)
        frames.append(pd.DataFrame({"tau": values, "orthology_class": cls}))
    return pd.concat(frames, ignore_index=True)


def _assert_png_sane(path: Path, min_dimension: int = 50) -> None:
    assert path.exists()
    assert path.stat().st_size > 0
    with Image.open(path) as img:
        width, height = img.size
        assert width >= min_dimension and height >= min_dimension
        assert img.format == "PNG"


class TestComputeTau:
    def test_tau_bounds_and_known_value(self) -> None:
        # Gene expressed in exactly one tissue: tau == 1. Ubiquitous: tau == 0.
        frame = pd.DataFrame(
            {
                "t1": [10.0, 5.0],
                "t2": [0.0, 5.0],
                "t3": [0.0, 5.0],
            },
            index=["specific", "ubiquitous"],
        )
        tau = compute_tau(frame)
        assert tau["specific"] == pytest.approx(1.0)
        assert tau["ubiquitous"] == pytest.approx(0.0)
        assert ((tau >= 0) & (tau <= 1)).all()

    def test_rejects_nan_and_negative(self) -> None:
        bad_nan = pd.DataFrame({"t1": [1.0, 1.0], "t2": [np.nan, 1.0]})
        with pytest.raises(ValueError):
            compute_tau(bad_nan)
        bad_neg = pd.DataFrame({"t1": [-1.0, 1.0], "t2": [1.0, 1.0]})
        with pytest.raises(ValueError):
            compute_tau(bad_neg)

    def test_all_zero_rows_dropped(self) -> None:
        frame = pd.DataFrame({"t1": [0.0, 1.0], "t2": [0.0, 1.0]}, index=["zero", "ok"])
        tau = compute_tau(frame)
        assert "zero" not in tau.index
        assert "ok" in tau.index


class TestTauHeatmap:
    def test_writes_sane_png(self, tmp_path: Path) -> None:
        out = tmp_path / "tau_heatmap.png"
        result = plot_tau_heatmap(_synthetic_tau_frame(), out)
        assert result == out
        _assert_png_sane(out)

    def test_deterministic_ordering_independent_of_input_order(
        self, tmp_path: Path
    ) -> None:
        frame = _synthetic_tau_frame()
        shuffled = frame.loc[reversed(list(frame.index)), reversed(list(frame.columns))]
        out_a = tmp_path / "a.png"
        out_b = tmp_path / "b.png"
        plot_tau_heatmap(frame, out_a)
        plot_tau_heatmap(shuffled, out_b)
        assert out_a.read_bytes() == out_b.read_bytes()

    def test_rejects_out_of_range(self, tmp_path: Path) -> None:
        bad = _synthetic_tau_frame()
        bad.iloc[0, 0] = 1.5
        with pytest.raises(ValueError):
            plot_tau_heatmap(bad, tmp_path / "x.png")


class TestOrthogroupSmallMultiples:
    def test_writes_sane_png(self, tmp_path: Path) -> None:
        out = tmp_path / "small_multiples.png"
        result = plot_orthogroup_small_multiples(_synthetic_long_expression(), out)
        _assert_png_sane(out)
        assert result == out

    def test_max_panels_limits_panels(self, tmp_path: Path) -> None:
        out = tmp_path / "two_panels.png"
        plot_orthogroup_small_multiples(_synthetic_long_expression(), out, max_panels=2)
        _assert_png_sane(out)

    def test_missing_column_raises(self, tmp_path: Path) -> None:
        frame = _synthetic_long_expression().drop(columns=["species"])
        with pytest.raises(ValueError, match="species"):
            plot_orthogroup_small_multiples(frame, tmp_path / "x.png")

    def test_negative_expression_raises(self, tmp_path: Path) -> None:
        frame = _synthetic_long_expression()
        frame.loc[frame.index[0], "expression"] = -1.0
        with pytest.raises(ValueError):
            plot_orthogroup_small_multiples(frame, tmp_path / "x.png")


class TestTauOrthologyStrips:
    def test_writes_sane_png(self, tmp_path: Path) -> None:
        out = tmp_path / "strips.png"
        result = plot_tau_orthology_strips(_synthetic_tau_classes(), out)
        _assert_png_sane(out)
        assert result == out

    def test_rejects_out_of_range_tau(self, tmp_path: Path) -> None:
        frame = _synthetic_tau_classes()
        frame.loc[frame.index[0], "tau"] = 2.0
        with pytest.raises(ValueError):
            plot_tau_orthology_strips(frame, tmp_path / "x.png")

    def test_missing_column_raises(self, tmp_path: Path) -> None:
        frame = _synthetic_tau_classes().drop(columns=["orthology_class"])
        with pytest.raises(ValueError, match="orthology_class"):
            plot_tau_orthology_strips(frame, tmp_path / "x.png")


class TestDescriptiveOnlyGuarantee:
    def test_no_inferential_statistics_in_module_source(self) -> None:
        import inspect

        from metainformant.rna.analysis import atlas_plots

        source = inspect.getsource(atlas_plots)
        for forbidden in (
            "pvalue",
            "p_value",
            "wilcoxon",
            "chi2",
            "ttest",
            "t_test",
            "scipy.stats",
        ):
            assert forbidden not in source.lower(), forbidden

    def test_palette_is_okabe_ito(self) -> None:
        assert OKABE_ITO[0] == "#E69F00"
        assert OKABE_ITO[4] == "#0072B2"
        assert len(OKABE_ITO) == 8
