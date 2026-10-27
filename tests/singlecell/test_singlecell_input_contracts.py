"""Input-mutation contract of the single-cell analysis functions.

Like the rest of the module, ``select_hvgs``, ``compute_neighbors`` and ``compute_diffusion_map`` return
a modified copy and leave their input untouched; ``inplace=True`` opts into writing into the input.
Real data and real library calls only.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

from metainformant.singlecell.analysis.nonlinear_methods import compute_diffusion_map, compute_neighbors
from metainformant.singlecell.analysis.pca_methods import select_hvgs
from metainformant.singlecell.data.preprocessing import SingleCellData


def _expression(n_cells: int = 60, n_genes: int = 12) -> SingleCellData:
    rng = np.random.default_rng(11)
    X = rng.gamma(2.0, 1.0, size=(n_cells, n_genes))
    var = pd.DataFrame({"symbol": [f"G{i}" for i in range(n_genes)]}, index=[f"g{i}" for i in range(n_genes)])
    return SingleCellData(X, var=var)


def _state(data: SingleCellData) -> dict[str, Any]:
    return {
        "X": data.X.copy(),
        "obs": data.obs.copy(),
        "var": data.var.copy(),
        "uns": sorted(data.uns),
        "obsm": sorted(data.obsm),
    }


def _assert_unchanged(data: SingleCellData, before: dict[str, Any]) -> None:
    np.testing.assert_array_equal(data.X, before["X"])
    pd.testing.assert_frame_equal(data.obs, before["obs"])
    pd.testing.assert_frame_equal(data.var, before["var"])
    assert sorted(data.uns) == before["uns"]
    assert sorted(data.obsm) == before["obsm"]


class TestSelectHvgs:
    def test_default_leaves_the_input_untouched(self) -> None:
        data = _expression()
        before = _state(data)

        result = select_hvgs(data, n_top_genes=5, flavor="variance")

        _assert_unchanged(data, before)
        assert result is not data
        assert int(result.var["highly_variable"].sum()) == 5
        assert result.uns["hvg"]["n_selected"] == 5

    def test_inplace_writes_into_the_input_and_returns_it(self) -> None:
        data = _expression()

        result = select_hvgs(data, n_top_genes=5, flavor="variance", inplace=True)

        assert result is data
        assert int(data.var["highly_variable"].sum()) == 5
        assert "hvg" in data.uns

    def test_copy_and_inplace_agree(self) -> None:
        copied = select_hvgs(_expression(), n_top_genes=4, flavor="seurat")
        in_place = select_hvgs(_expression(), n_top_genes=4, flavor="seurat", inplace=True)

        pd.testing.assert_frame_equal(copied.var, in_place.var)


class TestComputeNeighbors:
    def test_default_leaves_the_input_untouched(self) -> None:
        data = _expression()
        before = _state(data)

        result = compute_neighbors(data, n_neighbors=8)

        _assert_unchanged(data, before)
        assert "neighbors" not in data.uns
        assert result is not data
        assert result.uns["neighbors"]["connectivities"].shape == (60, 60)

    def test_inplace_writes_into_the_input_and_returns_it(self) -> None:
        data = _expression()

        result = compute_neighbors(data, n_neighbors=8, inplace=True)

        assert result is data
        assert data.uns["neighbors"]["params"]["n_neighbors"] == 8


class TestComputeDiffusionMap:
    def test_default_leaves_the_input_untouched_even_when_it_builds_the_graph(self) -> None:
        data = _expression()
        before = _state(data)

        result = compute_diffusion_map(data, n_components=4, n_neighbors=8)

        _assert_unchanged(data, before)
        assert {"neighbors", "diffmap"} <= set(result.uns)
        assert result.obsm["X_diffmap"].shape == (60, 4)

    def test_default_reuses_an_existing_graph_without_replacing_it(self) -> None:
        data = compute_neighbors(_expression(), n_neighbors=8)
        graph = data.uns["neighbors"]

        result = compute_diffusion_map(data, n_components=4)

        assert data.uns["neighbors"] is graph
        assert "diffmap" not in data.uns and "X_diffmap" not in data.obsm
        assert result.obsm["X_diffmap"].shape == (60, 4)

    def test_inplace_writes_into_the_input_and_returns_it(self) -> None:
        data = _expression()

        result = compute_diffusion_map(data, n_components=4, n_neighbors=8, inplace=True)

        assert result is data
        assert {"neighbors", "diffmap"} <= set(data.uns)
        assert data.obsm["X_diffmap"].shape == (60, 4)

    def test_copy_and_inplace_agree(self) -> None:
        copied = compute_diffusion_map(_expression(), n_components=4, n_neighbors=8)
        in_place = compute_diffusion_map(_expression(), n_components=4, n_neighbors=8, inplace=True)

        np.testing.assert_allclose(copied.obsm["X_diffmap"], in_place.obsm["X_diffmap"])

    def test_repeated_calls_give_identical_components(self) -> None:
        """ARPACK draws its start vector from state that advances per call; the components must not flip."""
        runs = [compute_diffusion_map(_expression(), n_components=4, n_neighbors=8) for _ in range(3)]

        for later in runs[1:]:
            np.testing.assert_allclose(later.obsm["X_diffmap"], runs[0].obsm["X_diffmap"], rtol=1e-10, atol=1e-12)


class TestCopyIndependence:
    def test_copy_does_not_share_nested_uns_dictionaries(self) -> None:
        data = _expression()
        data.uns["qc_summary"] = {"total_cells": 60}

        duplicate = data.copy()
        duplicate.uns["qc_summary"]["added_on_the_copy"] = True
        duplicate.uns["qc_summary"]["total_cells"] = -1

        assert data.uns["qc_summary"] == {"total_cells": 60}

    @pytest.mark.parametrize("function", [select_hvgs, compute_neighbors])
    def test_nested_uns_entries_of_the_input_survive_the_default_copy(self, function: Any) -> None:
        data = _expression()
        data.uns["provenance"] = {"steps": ["load"]}

        result = function(data)
        result.uns["provenance"]["steps"].append("edited on the result")

        assert data.uns["provenance"] == {"steps": ["load"]}
