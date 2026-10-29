"""``filter_cells`` never leaves a neighbor graph that describes cells which are no longer there.

``uns["neighbors"]`` (indices, distances and connectivities, one row per cell) used to be copied unchanged into
the filtered object, so its rows no longer matched its cells: the diffusion map then crashed or, with luck,
quietly used the wrong cells' neighbors. A k-nearest-neighbor graph of the retained cells is not a subset of the
old graph, so when any cell is removed the graph is dropped and the reason recorded; consumers rebuild it. Real
data and real calls only.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from metainformant.core.utils import errors
from metainformant.singlecell.analysis.nonlinear_methods import compute_diffusion_map, compute_neighbors
from metainformant.singlecell.analysis.trajectory import compute_diffusion_pseudotime, paga_trajectory
from metainformant.singlecell.data.preprocessing import SingleCellData, filter_cells

N_CELLS = 60
N_NEIGHBORS = 8


def _counts_with_empty_cells(n_empty: int = 6, seed: int = 0) -> SingleCellData:
    """Poisson counts in which the last ``n_empty`` cells have almost nothing, so ``min_counts`` removes them."""
    rng = np.random.default_rng(seed)
    X = rng.poisson(5.0, (N_CELLS, 30)).astype(float)
    X[N_CELLS - n_empty :] = 0.0
    X[N_CELLS - n_empty :, 0] = 1.0
    obs = pd.DataFrame({"cluster": ["A", "B"] * (N_CELLS // 2)}, index=[f"cell{i}" for i in range(N_CELLS)])
    return SingleCellData(X, obs=obs)


def _with_graph(n_empty: int = 6) -> SingleCellData:
    return compute_neighbors(_counts_with_empty_cells(n_empty), n_neighbors=N_NEIGHBORS)


class TestTheGraphIsDroppedWhenCellsAreRemoved:
    def test_no_stale_graph_survives(self) -> None:
        data = _with_graph()

        filtered = filter_cells(data, min_counts=10)

        assert filtered.n_obs == N_CELLS - 6
        assert "neighbors" not in filtered.uns

    def test_the_reason_the_counts_and_the_parameters_are_recorded(self) -> None:
        data = _with_graph()

        filtered = filter_cells(data, min_counts=10)

        record = filtered.uns["neighbors_dropped"]
        assert record["n_cells_before"] == N_CELLS and record["n_cells_after"] == N_CELLS - 6
        assert "removed 6 of 60 cells" in record["reason"] and "compute_neighbors" in record["reason"]
        assert record["params"]["n_neighbors"] == N_NEIGHBORS  # enough to recompute the same graph
        assert record["params"] == data.uns["neighbors"]["params"]

    def test_the_input_object_keeps_its_graph(self) -> None:
        data = _with_graph()
        before = data.uns["neighbors"]["indices"].copy()

        filter_cells(data, min_counts=10)

        np.testing.assert_array_equal(data.uns["neighbors"]["indices"], before)
        assert "neighbors_dropped" not in data.uns

    def test_the_drop_is_logged(self, caplog: pytest.LogCaptureFixture) -> None:
        data = _with_graph()

        with caplog.at_level("WARNING"):
            filter_cells(data, min_counts=10)

        assert any("dropped uns['neighbors']" in record.getMessage() for record in caplog.records)

    @pytest.mark.parametrize("keep_fraction", [0.5, 0.9])
    def test_the_graph_is_dropped_however_many_cells_go(self, keep_fraction: float) -> None:
        data = _with_graph()
        threshold = np.quantile(data.X.sum(axis=1), 1.0 - keep_fraction)

        filtered = filter_cells(data, min_counts=threshold)

        assert filtered.n_obs < N_CELLS
        assert "neighbors" not in filtered.uns


class TestEveryConsumerWorksOnTheFilteredData:
    def test_the_diffusion_map_rebuilds_its_graph_for_the_retained_cells(self) -> None:
        filtered = filter_cells(_with_graph(), min_counts=10)

        result = compute_diffusion_map(filtered, n_components=4, n_neighbors=N_NEIGHBORS)

        assert result.obsm["X_diffmap"].shape == (N_CELLS - 6, 4)
        assert np.all(np.isfinite(result.obsm["X_diffmap"]))
        assert result.uns["neighbors"]["indices"].shape[0] == N_CELLS - 6  # a fresh graph of the retained cells
        assert result.uns["neighbors"]["indices"].max() < N_CELLS - 6

    def test_diffusion_pseudotime_covers_the_retained_cells(self) -> None:
        filtered = filter_cells(_with_graph(), min_counts=10)

        result = compute_diffusion_pseudotime(filtered, root_cell=0, n_components=3)

        assert len(result.obs["dpt_pseudotime"]) == N_CELLS - 6
        assert result.uns["dpt"]["neighbors_source"] == "X"

    def test_paga_covers_the_retained_cells(self) -> None:
        filtered = filter_cells(_with_graph(), min_counts=10)

        result = paga_trajectory(filtered, groups="cluster", n_neighbors=N_NEIGHBORS)

        assert result.uns["paga"]["neighbors_source"] == "X"
        assert len(result.obs["paga_pseudotime"]) == N_CELLS - 6

    def test_a_graph_that_slipped_through_is_refused_by_the_diffusion_map(self) -> None:
        """Defence in depth: an object assembled by hand with a stale graph raises instead of using wrong cells."""
        data = _with_graph()
        stale = SingleCellData(data.X[: N_CELLS - 6].copy(), uns=dict(data.uns))

        with pytest.raises(errors.ValidationError, match=r"60 cells but data has 54"):
            compute_diffusion_map(stale, n_components=4)


class TestWhatStaysValid:
    def test_filtering_nothing_keeps_the_graph_untouched(self) -> None:
        data = _with_graph(n_empty=0)
        before = {key: data.uns["neighbors"][key] for key in ("indices", "distances")}

        filtered = filter_cells(data, min_counts=0)

        assert filtered.n_obs == N_CELLS
        np.testing.assert_array_equal(filtered.uns["neighbors"]["indices"], before["indices"])
        np.testing.assert_array_equal(filtered.uns["neighbors"]["distances"], before["distances"])
        assert (filtered.uns["neighbors"]["connectivities"] != data.uns["neighbors"]["connectivities"]).nnz == 0
        assert "neighbors_dropped" not in filtered.uns

    def test_a_parameter_only_neighbors_entry_is_kept_and_obsp_is_subset(self) -> None:
        """scanpy keeps only parameters in ``uns["neighbors"]``; its matrices live in ``obsp`` and are subset."""
        data = _counts_with_empty_cells()
        connectivities = sparse.random(N_CELLS, N_CELLS, density=0.2, format="csr", random_state=1)
        data.obsp["connectivities"] = connectivities
        data.uns["neighbors"] = {"connectivities_key": "connectivities", "params": {"n_neighbors": 15}}

        filtered = filter_cells(data, min_counts=10)

        assert filtered.uns["neighbors"] == {"connectivities_key": "connectivities", "params": {"n_neighbors": 15}}
        assert "neighbors_dropped" not in filtered.uns
        assert filtered.obsp["connectivities"].shape == (N_CELLS - 6, N_CELLS - 6)

    def test_data_without_a_graph_gets_no_record(self) -> None:
        filtered = filter_cells(_counts_with_empty_cells(), min_counts=10)

        assert "neighbors" not in filtered.uns and "neighbors_dropped" not in filtered.uns
