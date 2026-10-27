"""Spatial plots put the right coordinate column on each axis.

``SpatialDataset.coordinates`` (Visium) is ``[row, col]`` but the Xenium and MERFISH loaders return
``[x, y]``; every plot used to treat column 0 as the vertical axis, so Xenium and MERFISH tissue
came out transposed. The tests read the positions matplotlib actually drew. Real matplotlib (Agg) only.
"""

from __future__ import annotations

from dataclasses import fields
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

from metainformant.spatial.io.merfish import MERFISHDataset
from metainformant.spatial.io.visium import SpatialDataset
from metainformant.spatial.io.xenium import XeniumDataset
from metainformant.spatial.visualization.plots import (
    plot_cell_type_map,
    plot_deconvolution_pie,
    plot_gene_expression_map,
    plot_neighborhood_graph,
    plot_spatial_autocorrelation,
    plot_spatial_scatter,
    plot_tissue_overlay,
)

matplotlib.use("Agg")  # Non-interactive backend for tests; selected before any figure is created

# Three spots whose two columns differ everywhere, so a swap cannot go unnoticed
COORDS = np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]])
VALUES = np.array([0.1, 0.5, 0.9])
GENES = ["G0", "G1"]
EXPRESSION = np.column_stack([VALUES, 1.0 - VALUES])
NAMES = ["a", "b", "c"]


@pytest.fixture(autouse=True)
def _close_plots():
    yield
    plt.close("all")


def _drawn_points(fig: Any) -> np.ndarray:
    """Every scatter position on the main axes, sorted by horizontal then vertical coordinate."""
    points = np.vstack([c.get_offsets() for c in fig.axes[0].collections])
    return points[np.lexsort((points[:, 1], points[:, 0]))]


def _sorted_rows(points: np.ndarray) -> np.ndarray:
    return points[np.lexsort((points[:, 1], points[:, 0]))]


class TestPlotSpatialScatter:
    def test_xy_order_puts_column_zero_on_the_horizontal_axis(self, tmp_path: Path) -> None:
        fig = plot_spatial_scatter(COORDS, VALUES, tmp_path / "xy.png", coordinate_order="xy")

        np.testing.assert_allclose(_drawn_points(fig), COORDS)

    def test_row_col_order_puts_column_one_on_the_horizontal_axis(self, tmp_path: Path) -> None:
        fig = plot_spatial_scatter(COORDS, VALUES, tmp_path / "rc.png", coordinate_order="row_col")

        np.testing.assert_allclose(_drawn_points(fig), _sorted_rows(COORDS[:, ::-1]))

    def test_the_default_stays_row_col_so_visium_arrays_plot_as_before(self, tmp_path: Path) -> None:
        default = plot_spatial_scatter(COORDS, VALUES, tmp_path / "default.png")
        explicit = plot_spatial_scatter(COORDS, VALUES, tmp_path / "explicit.png", coordinate_order="row_col")

        np.testing.assert_allclose(_drawn_points(default), _drawn_points(explicit))

    def test_categorical_values_follow_the_same_convention(self, tmp_path: Path) -> None:
        fig = plot_spatial_scatter(COORDS, np.array(NAMES), tmp_path / "cat.png", coordinate_order="xy")

        np.testing.assert_allclose(_drawn_points(fig), COORDS)

    def test_an_unknown_order_is_rejected(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="coordinate_order"):
            plot_spatial_scatter(COORDS, VALUES, tmp_path / "bad.png", coordinate_order="yx")  # type: ignore[arg-type]


class TestDatasetsDeclareTheirConvention:
    def test_the_declared_orders(self) -> None:
        assert SpatialDataset.coordinate_order == "row_col"
        assert XeniumDataset.coordinate_order == "xy"
        assert MERFISHDataset.coordinate_order == "xy"

    def test_the_declaration_is_not_a_constructor_field(self) -> None:
        dataset = XeniumDataset(expression=EXPRESSION, coordinates=COORDS, cell_ids=NAMES, gene_names=GENES)

        assert dataset.coordinate_order == "xy"
        assert "coordinate_order" not in {field.name for field in fields(XeniumDataset)}


class TestDatasetDrivenPlots:
    @staticmethod
    def _datasets() -> dict[str, Any]:
        return {
            "xenium": XeniumDataset(expression=EXPRESSION, coordinates=COORDS, cell_ids=NAMES, gene_names=GENES),
            "merfish": MERFISHDataset(expression=EXPRESSION, coordinates=COORDS, cell_ids=NAMES, gene_names=GENES),
            "visium": SpatialDataset(expression=EXPRESSION, coordinates=COORDS, barcodes=NAMES, gene_names=GENES),
        }

    @pytest.mark.parametrize("platform, expected_is_xy", [("xenium", True), ("merfish", True), ("visium", False)])
    def test_the_gene_map_reads_the_order_from_the_dataset(
        self, tmp_path: Path, platform: str, expected_is_xy: bool
    ) -> None:
        dataset = self._datasets()[platform]

        fig = plot_gene_expression_map(dataset, "G0", tmp_path / f"{platform}.png")

        expected = COORDS if expected_is_xy else COORDS[:, ::-1]
        np.testing.assert_allclose(_drawn_points(fig), _sorted_rows(expected))

    @pytest.mark.parametrize("platform, expected_is_xy", [("xenium", True), ("merfish", True), ("visium", False)])
    def test_the_cell_type_map_reads_the_order_from_the_dataset(
        self, tmp_path: Path, platform: str, expected_is_xy: bool
    ) -> None:
        dataset = self._datasets()[platform]

        fig = plot_cell_type_map(dataset, np.array(NAMES), tmp_path / f"{platform}_types.png")

        expected = COORDS if expected_is_xy else COORDS[:, ::-1]
        np.testing.assert_allclose(_drawn_points(fig), _sorted_rows(expected))

    def test_an_explicit_order_overrides_the_dataset(self, tmp_path: Path) -> None:
        dataset = self._datasets()["xenium"]

        fig = plot_gene_expression_map(dataset, "G0", tmp_path / "override.png", coordinate_order="row_col")

        np.testing.assert_allclose(_drawn_points(fig), _sorted_rows(COORDS[:, ::-1]))

    def test_a_bare_array_keeps_the_legacy_row_col_reading(self, tmp_path: Path) -> None:
        fig = plot_cell_type_map(COORDS, np.array(NAMES), tmp_path / "bare.png")

        np.testing.assert_allclose(_drawn_points(fig), _sorted_rows(COORDS[:, ::-1]))


class TestTheOtherPlotsShareTheConvention:
    def test_neighborhood_graph_nodes_and_edges(self, tmp_path: Path) -> None:
        graph = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=float)

        fig = plot_neighborhood_graph(COORDS, graph, tmp_path / "graph.png", coordinate_order="xy")

        ax = fig.axes[0]
        np.testing.assert_allclose(_drawn_points(fig), COORDS)
        np.testing.assert_allclose(ax.lines[0].get_xdata(), [1.0, 2.0])  # x of spots 0 and 1
        np.testing.assert_allclose(ax.lines[0].get_ydata(), [10.0, 20.0])  # y of spots 0 and 1

    def test_deconvolution_pies_are_centred_on_the_spots(self, tmp_path: Path) -> None:
        fractions = np.array([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])

        fig = plot_deconvolution_pie(COORDS, fractions, tmp_path / "pie.png", pie_radius=0.2, coordinate_order="xy")

        centres = np.array([patch.center for patch in fig.axes[0].patches])
        np.testing.assert_allclose(_sorted_rows(centres), COORDS)

    def test_lisa_map_with_labels_and_with_scores(self, tmp_path: Path) -> None:
        labelled = plot_spatial_autocorrelation(
            COORDS, VALUES, tmp_path / "lisa_labels.png", cluster_labels=["HH", "LL", "NS"], coordinate_order="xy"
        )
        continuous = plot_spatial_autocorrelation(COORDS, VALUES, tmp_path / "lisa.png", coordinate_order="xy")

        np.testing.assert_allclose(_drawn_points(labelled), COORDS)
        np.testing.assert_allclose(_drawn_points(continuous), COORDS)

    def test_tissue_overlay_scales_then_orders(self, tmp_path: Path) -> None:
        image = np.zeros((70, 70, 3), dtype=np.uint8)

        fig = plot_tissue_overlay(
            COORDS, VALUES, image, tmp_path / "overlay.png", scale_factor=2.0, coordinate_order="xy"
        )

        np.testing.assert_allclose(_drawn_points(fig), 2.0 * COORDS)
