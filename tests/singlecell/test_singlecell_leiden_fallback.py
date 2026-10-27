"""leiden_clustering never substitutes k-means silently.

Real data and real library calls only; nothing is replaced by a test double. The two classes are
gated on whether the Leiden stack (networkx, python-igraph, leidenalg) is installed, so each runs in
the environments where its contract can actually be observed.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from metainformant.singlecell.analysis.clustering import (
    HAS_IGRAPH,
    HAS_LEIDEN,
    HAS_NETWORKX,
    leiden_clustering,
)
from metainformant.singlecell.data.preprocessing import SingleCellData

LEIDEN_STACK_AVAILABLE = bool(HAS_NETWORKX and HAS_IGRAPH and HAS_LEIDEN)


def _three_blobs(per_blob: int = 20, n_genes: int = 6) -> np.ndarray:
    rng = np.random.default_rng(1)
    return np.vstack([rng.normal(center, 0.2, (per_blob, n_genes)) for center in (0.0, 5.0, 10.0)])


@pytest.mark.skipif(LEIDEN_STACK_AVAILABLE, reason="needs an environment without the Leiden stack")
class TestWithoutTheLeidenStack:
    def test_feature_matrix_raises_instead_of_running_kmeans(self) -> None:
        with pytest.raises(ImportError, match="allow_kmeans_fallback"):
            leiden_clustering(feature_matrix=_three_blobs())

    def test_positional_matrix_raises_too(self) -> None:
        with pytest.raises(ImportError, match="leidenalg"):
            leiden_clustering(_three_blobs())

    def test_single_cell_data_is_pointed_at_kmeans_clustering(self) -> None:
        with pytest.raises(ImportError, match="kmeans_clustering"):
            leiden_clustering(SingleCellData(_three_blobs()))

    def test_explicit_fallback_runs_kmeans_and_says_so(self) -> None:
        with pytest.warns(UserWarning, match="k-means"):
            labels = leiden_clustering(_three_blobs(), resolution=0.75, random_state=3, allow_kmeans_fallback=True)

        # resolution 0.75 asks for round(0.75 * 4) = 3 clusters; the three blobs are recovered exactly
        assert len(labels) == 60
        assert len(set(labels)) == 3
        for start in (0, 20, 40):
            assert len(set(labels[start : start + 20])) == 1

    def test_the_fallback_flag_does_not_rescue_single_cell_data(self) -> None:
        with pytest.raises(ImportError, match="kmeans_clustering"):
            leiden_clustering(SingleCellData(_three_blobs()), allow_kmeans_fallback=True)


@pytest.mark.skipif(not LEIDEN_STACK_AVAILABLE, reason="needs networkx, python-igraph and leidenalg")
class TestWithTheLeidenStack:
    def test_matrix_input_is_partitioned_by_leiden(self) -> None:
        X = _three_blobs()

        labels = leiden_clustering(feature_matrix=X, resolution=1.0, n_neighbors=8, random_state=0)
        data = leiden_clustering(SingleCellData(X), resolution=1.0, n_neighbors=8, random_state=0)

        np.testing.assert_array_equal(labels, data.obs["leiden_cluster"].to_numpy())
        assert data.uns["leiden_clustering"]["method"] == "leiden"

    def test_the_fallback_flag_is_inert_when_leiden_is_available(self) -> None:
        X = _three_blobs()

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with_flag = leiden_clustering(X, n_neighbors=8, random_state=0, allow_kmeans_fallback=True)

        assert not [w for w in caught if "k-means" in str(w.message)]
        np.testing.assert_array_equal(with_flag, leiden_clustering(X, n_neighbors=8, random_state=0))


class TestMatrixValidation:
    """Holds with or without the Leiden stack: the shape check happens before any dependency check."""

    def test_one_dimensional_input_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="2D"):
            leiden_clustering(feature_matrix=np.zeros(5))
