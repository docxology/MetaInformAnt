"""PAGA connectivity is computed from the kNN graph, not from which row follows which.

``paga_trajectory`` used to count transitions between adjacent rows of the data, so its "connectivity" was a
property of how the rows happened to be sorted. The connectivity of Wolf et al. (2019) compares the number of
kNN edges between two groups with the number expected if the same edges ended at random cells. These tests
check the formula on a graph small enough to count by hand, then planted groups on a trajectory whose cells
are shuffled. Real data and real kNN searches only.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.core.utils import errors
from metainformant.singlecell.analysis import trajectory
from metainformant.singlecell.analysis.nonlinear_methods import compute_neighbors
from metainformant.singlecell.analysis.trajectory import paga_trajectory
from metainformant.singlecell.data.preprocessing import SingleCellData

GROUP_NAMES = ("A", "B", "C")


def _chain(n_per_group: int = 30, seed: int = 0, shuffle: bool = True, groups: tuple[str, ...] = GROUP_NAMES):
    """Cells along a 1-D trajectory cut into consecutive groups, so only neighboring groups touch."""
    rng = np.random.default_rng(seed)
    n_cells = len(groups) * n_per_group
    t = np.linspace(0.0, float(len(groups)), n_cells, endpoint=False)
    X = np.column_stack([t, 0.002 * rng.standard_normal((n_cells, 4))])
    labels = np.repeat(np.array(groups), n_per_group)
    if shuffle:
        order = rng.permutation(n_cells)
        X, labels = X[order], labels[order]
    obs = pd.DataFrame({"group": labels}, index=[f"cell{i}" for i in range(n_cells)])
    return SingleCellData(X, obs=obs)


def _matrix(result: SingleCellData) -> np.ndarray:
    return np.array(result.uns["paga"]["connectivity_matrix"])


class TestTheFormulaOnAGraphCountedByHand:
    # Cells 0-2 are group 0 and cells 3-5 group 1. Row v lists the two neighbors of cell v.
    GROUPS = np.array([0, 0, 0, 1, 1, 1])
    NEIGHBORS = np.array([[1, 3], [0, 2], [1, 3], [4, 2], [3, 5], [4, 3]])

    def test_observed_expected_and_connectivity(self) -> None:
        connectivity, observed, expected = trajectory._compute_cluster_connectivity(self.GROUPS, 2, self.NEIGHBORS)

        # Directed edges 0 -> 3 and 2 -> 3 leave group 0 for group 1, edge 3 -> 2 comes back: 3 in all
        assert observed[0, 1] == observed[1, 0] == 3.0
        # Each group starts 3 cells x 2 = 6 edges and has 3 cells: (6 * 3 + 6 * 3) / (6 - 1)
        assert expected[0, 1] == expected[1, 0] == pytest.approx(36.0 / 5.0)
        assert connectivity[0, 1] == connectivity[1, 0] == pytest.approx(3.0 / 7.2)
        assert connectivity[0, 0] == connectivity[1, 1] == 0.0

    def test_the_connectivity_is_capped_at_one(self) -> None:
        # Every cell lists only cells of the other group: 6 observed edges from each side against 36/5 expected...
        groups = np.array([0, 0, 1, 1])
        neighbors = np.array([[2, 3], [2, 3], [0, 1], [0, 1]])

        connectivity, observed, expected = trajectory._compute_cluster_connectivity(groups, 2, neighbors)

        assert observed[0, 1] == 8.0
        assert expected[0, 1] == pytest.approx((4 * 2 + 4 * 2) / 3.0)  # 16/3 < 8
        assert connectivity[0, 1] == 1.0

    def test_a_cells_own_entry_in_its_neighbor_row_is_not_an_edge(self) -> None:
        with_self = np.array([[0, 1, 3], [1, 0, 2], [2, 1, 3], [3, 4, 2], [4, 3, 5], [5, 4, 3]])
        without_self = np.array([[1, 3], [0, 2], [1, 3], [4, 2], [3, 5], [4, 3]])

        a = trajectory._compute_cluster_connectivity(self.GROUPS, 2, with_self)
        b = trajectory._compute_cluster_connectivity(self.GROUPS, 2, without_self)

        for left, right in zip(a, b):
            np.testing.assert_allclose(left, right)

    def test_it_is_symmetric_with_a_zero_diagonal_and_in_the_unit_interval(self) -> None:
        rng = np.random.default_rng(1)
        groups = rng.integers(0, 4, size=60)
        neighbors = np.array([rng.choice(np.delete(np.arange(60), v), size=5, replace=False) for v in range(60)])

        connectivity, observed, expected = trajectory._compute_cluster_connectivity(groups, 4, neighbors)

        for matrix in (connectivity, observed, expected):
            np.testing.assert_array_equal(matrix, matrix.T)
            assert np.all(np.diag(matrix) == 0.0)
        assert connectivity.min() >= 0.0 and connectivity.max() <= 1.0

    def test_relabelling_the_cells_does_not_change_the_connectivity(self) -> None:
        """The old count read adjacent rows; a graph quantity cannot depend on the order of the cells."""
        rng = np.random.default_rng(2)
        groups = rng.integers(0, 3, size=40)
        neighbors = np.array([rng.choice(np.delete(np.arange(40), v), size=4, replace=False) for v in range(40)])
        order = rng.permutation(40)
        position = np.argsort(order)  # old cell -> new row

        original = trajectory._compute_cluster_connectivity(groups, 3, neighbors)[0]
        shuffled = trajectory._compute_cluster_connectivity(groups[order], 3, position[neighbors[order]])[0]

        np.testing.assert_allclose(original, shuffled)


class TestPlantedGroupsOnAShuffledTrajectory:
    def test_only_neighboring_groups_are_connected(self) -> None:
        result = paga_trajectory(_chain(), groups="group", n_neighbors=10)

        connectivity = _matrix(result)
        assert result.uns["paga"]["unique_groups"] == list(GROUP_NAMES)
        assert connectivity[0, 1] > 0.1 and connectivity[1, 2] > 0.1  # A-B and B-C touch
        assert connectivity[0, 2] == 0.0  # A and C share no edge
        assert connectivity.max() < 1.0  # a boundary is far from a random mixture
        np.testing.assert_array_equal(connectivity, connectivity.T)

    def test_the_order_of_the_rows_makes_no_difference(self) -> None:
        sorted_rows = _matrix(paga_trajectory(_chain(shuffle=False), groups="group", n_neighbors=10))
        shuffled_rows = _matrix(paga_trajectory(_chain(shuffle=True), groups="group", n_neighbors=10))

        # A shuffled copy of the same cells has the same kNN graph, hence the same connectivity
        np.testing.assert_allclose(sorted_rows, shuffled_rows, atol=1e-12)

    def test_groups_unrelated_to_the_graph_are_connected_as_often_as_chance(self) -> None:
        """Random labels are the null model, so observed ~ expected and the connectivity reaches about 1."""
        data = _chain()
        data.obs["group"] = np.random.default_rng(3).choice(["A", "B", "C"], size=data.n_obs)

        connectivity = _matrix(paga_trajectory(data, groups="group", n_neighbors=10))

        off_diagonal = connectivity[~np.eye(3, dtype=bool)]
        assert off_diagonal.min() > 0.8

    def test_the_observed_and_expected_edges_are_stored(self) -> None:
        result = paga_trajectory(_chain(), groups="group", n_neighbors=10)

        info = result.uns["paga"]
        observed, expected = np.array(info["observed_edges"]), np.array(info["expected_edges"])
        np.testing.assert_allclose(np.minimum(observed / np.where(expected > 0, expected, 1.0), 1.0), _matrix(result))
        assert info["n_groups"] == 3 and info["n_neighbors"] == 10 and info["neighbors_source"] == "X"

    def test_the_stored_neighbor_graph_is_used_when_there_is_one(self) -> None:
        data = compute_neighbors(_chain(), n_neighbors=7)

        result = paga_trajectory(data, groups="group")

        assert result.uns["paga"]["neighbors_source"] == "uns"
        assert result.uns["paga"]["n_neighbors"] == 7
        assert result.uns["neighbors"]["params"]["n_neighbors"] == 7  # not replaced
        assert "paga" not in data.uns  # the input is left alone

    def test_a_stale_neighbor_graph_is_refused(self) -> None:
        data = compute_neighbors(_chain(), n_neighbors=7)
        fewer = SingleCellData(data.X[:60].copy(), obs=data.obs.iloc[:60].copy(), uns=dict(data.uns))

        with pytest.raises(errors.ValidationError, match="90 cells but data has 60"):
            paga_trajectory(fewer, groups="group")


class TestGroupPseudotime:
    def test_hops_from_the_chosen_root_along_the_chain(self) -> None:
        result = paga_trajectory(_chain(), groups="group", n_neighbors=10, root_group="A")

        by_group = result.obs.groupby("group")["paga_pseudotime"].first()
        assert by_group["A"] == 0.0 and by_group["B"] == 0.5 and by_group["C"] == 1.0
        assert result.uns["paga"]["root_group"] == "A" and result.uns["paga"]["root_selection"] == "user"

        reversed_root = paga_trajectory(_chain(), groups="group", n_neighbors=10, root_group="C")
        by_group = reversed_root.obs.groupby("group")["paga_pseudotime"].first()
        assert by_group["C"] == 0.0 and by_group["B"] == 0.5 and by_group["A"] == 1.0

    def test_without_a_root_the_heuristic_picks_an_end_of_the_chain_and_says_so(self) -> None:
        result = paga_trajectory(_chain(), groups="group", n_neighbors=10)

        info = result.uns["paga"]
        assert info["root_selection"] == "heuristic"
        assert info["root_group"] in ("A", "C")
        by_group = result.obs.groupby("group")["paga_pseudotime"].first()
        assert by_group["B"] == 0.5

    def test_a_group_the_graph_does_not_link_to_the_root_gets_infinity(self) -> None:
        data = _chain(groups=("A", "B", "C", "D"))
        # Move group D far from everything: no edge reaches it
        far = (data.obs["group"] == "D").to_numpy()
        data.X[far, 0] += 1000.0

        result = paga_trajectory(data, groups="group", n_neighbors=10, root_group="A")

        by_group = result.obs.groupby("group")["paga_pseudotime"].first()
        assert np.isinf(by_group["D"])
        assert by_group["A"] == 0.0 and by_group["C"] == 1.0

    def test_a_raised_threshold_cuts_weak_links(self) -> None:
        connectivity = _matrix(paga_trajectory(_chain(), groups="group", n_neighbors=10))
        strongest_edge = float(connectivity.max())

        result = paga_trajectory(_chain(), groups="group", n_neighbors=10, threshold=1.0, root_group="A")

        assert strongest_edge < 1.0  # so no edge reaches threshold 1.0
        by_group = result.obs.groupby("group")["paga_pseudotime"].first()
        assert by_group["A"] == 0.0 and np.isinf(by_group["B"]) and np.isinf(by_group["C"])


class TestValidation:
    def test_unknown_groups_column(self) -> None:
        with pytest.raises(errors.ValidationError, match="not found"):
            paga_trajectory(_chain(), groups="nope")

    def test_missing_group_labels(self) -> None:
        data = _chain()
        data.obs["group"] = data.obs["group"].astype(object)
        data.obs.iloc[3, 0] = None

        with pytest.raises(errors.ValidationError, match="missing values"):
            paga_trajectory(data, groups="group")

    def test_a_root_that_is_not_a_group(self) -> None:
        with pytest.raises(errors.ValidationError, match="root_group"):
            paga_trajectory(_chain(), groups="group", root_group="Z")

    @pytest.mark.parametrize("threshold", [-0.1, 1.5])
    def test_threshold_out_of_range(self, threshold: float) -> None:
        with pytest.raises(errors.ValidationError, match="threshold"):
            paga_trajectory(_chain(), groups="group", threshold=threshold)

    def test_a_single_cell_has_no_graph(self) -> None:
        data = SingleCellData(np.zeros((1, 3)), obs=pd.DataFrame({"group": ["A"]}))

        with pytest.raises(errors.ValidationError, match="at least 2 cells"):
            paga_trajectory(data, groups="group")
