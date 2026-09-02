"""Zero-mocks tests for visualization chart-convention utilities.

All tests use real matplotlib Agg figures, real files in tmp_path, and real
byte comparisons of saved output. Real-implementation policy: no test doubles.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from metainformant.visualization.config.conventions import (
    CATEGORY_HATCHES,
    OKABE_ITO,
    add_log_zero_mark,
    apply_deterministic_rcparams,
    category_style,
    okabe_ito,
    save_figure_deterministic,
)


class TestOkabeIto:
    def test_palette_is_fixed_order_and_complete(self) -> None:
        assert OKABE_ITO == (
            "#E69F00",
            "#56B4E9",
            "#009E73",
            "#F0E442",
            "#0072B2",
            "#D55E00",
            "#CC79A7",
            "#000000",
        )
        assert len(set(OKABE_ITO)) == 8

    def test_okabe_ito_returns_requested_count(self) -> None:
        colors = okabe_ito(5)
        assert colors == list(OKABE_ITO[:5])
        assert len(okabe_ito(8)) == 8

    def test_okabe_ito_cycles_for_large_n(self) -> None:
        colors = okabe_ito(10)
        assert len(colors) == 10
        assert colors[8] == OKABE_ITO[0]
        assert colors[9] == OKABE_ITO[1]

    def test_okabe_ito_rejects_nonpositive_n(self) -> None:
        with pytest.raises(ValueError):
            okabe_ito(0)
        with pytest.raises(ValueError):
            okabe_ito(-3)


class TestCategoryStyle:
    def test_mapping_is_deterministic_and_sorted(self) -> None:
        labels = ["many_to_many", "one_to_one", "one_to_many"]
        first = category_style(labels)
        second = category_style(list(reversed(labels)))
        assert first == second
        # Sorted unique order: many_to_many, one_to_many, one_to_one
        assert list(first) == ["many_to_many", "one_to_many", "one_to_one"]

    def test_every_style_pairs_color_with_hatch(self) -> None:
        styles = category_style(["a", "b", "c", "d", "e"])
        for color, hatch in styles.values():
            assert color in OKABE_ITO
            assert hatch in CATEGORY_HATCHES

    def test_duplicate_labels_collapse(self) -> None:
        styles = category_style(["x", "x", "y"])
        assert set(styles) == {"x", "y"}


class TestAddLogZeroMark:
    def test_draws_reference_line_on_log_y_axis(self) -> None:
        fig, ax = plt.subplots()
        ax.set_yscale("log")
        ax.set_ylim(1e-6, 10)
        ax.set_ylabel("tau")
        add_log_zero_mark(ax, "y")
        # One new axhline at the axis floor
        floors = [line.get_ydata()[0] for line in ax.get_lines() if len(line.get_ydata()) == 2]
        assert 1e-6 in floors
        assert "zero at axis floor" in ax.get_ylabel()
        plt.close(fig)

    def test_draws_reference_line_on_log_x_axis(self) -> None:
        fig, ax = plt.subplots()
        ax.set_xscale("log")
        ax.set_xlim(1e-3, 10)
        add_log_zero_mark(ax, "x", label="zero count")
        assert "zero count at axis floor" in ax.get_xlabel()
        plt.close(fig)

    def test_appends_to_existing_label_without_duplicating(self) -> None:
        fig, ax = plt.subplots()
        ax.set_yscale("log")
        ax.set_ylabel("divergence")
        add_log_zero_mark(ax, "y")
        assert ax.get_ylabel() == "divergence (zero at axis floor)"
        plt.close(fig)

    def test_refuses_linear_axis(self) -> None:
        fig, ax = plt.subplots()
        ax.set_yscale("linear")
        with pytest.raises(ValueError, match="not log-scaled"):
            add_log_zero_mark(ax, "y")
        plt.close(fig)

    def test_rejects_unknown_axis(self) -> None:
        fig, ax = plt.subplots()
        ax.set_yscale("log")
        with pytest.raises(ValueError, match="axis must be"):
            add_log_zero_mark(ax, "z")
        plt.close(fig)


class TestDeterministicExport:
    def setup_method(self) -> None:
        plt.rcdefaults()

    def teardown_method(self) -> None:
        plt.rcdefaults()

    def test_saved_pngs_are_byte_identical_across_runs(self, tmp_path: Path) -> None:
        paths = []
        for run in range(2):
            apply_deterministic_rcparams()
            fig, ax = plt.subplots()
            ax.plot([1, 2, 3], [4, 5, 6], color=OKABE_ITO[4])
            ax.set_title("deterministic")
            path = tmp_path / f"run{run}.png"
            save_figure_deterministic(fig, path)
            paths.append(path)
            plt.close(fig)
        assert paths[0].read_bytes() == paths[1].read_bytes()

    def test_apply_deterministic_rcparams_set_real_params(self) -> None:
        import matplotlib as mpl

        apply_deterministic_rcparams()
        assert mpl.rcParams["ps.useafm"] is True
        assert mpl.rcParams["svg.fonttype"] == "none"

    def test_save_without_metadata_still_embeds_defaults(self, tmp_path: Path) -> None:
        # Negative control: default matplotlib savefig embeds producer
        # metadata, which is exactly what deterministic mode removes.
        fig, ax = plt.subplots()
        ax.plot([1], [1])
        path = tmp_path / "default.png"
        fig.savefig(path)
        plt.close(fig)
        # Default PNG must differ from a deterministic save of the same figure
        fig, ax = plt.subplots()
        ax.plot([1], [1])
        det_path = tmp_path / "det.png"
        save_figure_deterministic(fig, det_path)
        plt.close(fig)
        assert path.read_bytes() != det_path.read_bytes()
