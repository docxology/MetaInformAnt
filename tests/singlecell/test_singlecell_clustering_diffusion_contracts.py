"""Input-validation and fallback contracts of hierarchical clustering and the diffusion map.

Real data and real scipy/numpy calls only; nothing is replaced by a test double.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.core.utils import errors
from metainformant.singlecell.analysis.clustering import hierarchical_clustering
from metainformant.singlecell.analysis.nonlinear_methods import compute_diffusion_map, compute_neighbors
from metainformant.singlecell.data.preprocessing import SingleCellData

_LINKAGE_METHODS = ["single", "complete", "average", "weighted", "centroid", "median", "ward"]


def _two_blob_data(cells_per_blob: int = 15, n_genes: int = 6) -> SingleCellData:
    rng = np.random.default_rng(3)
    blob_a = rng.normal(0.0, 0.2, (cells_per_blob, n_genes))
    blob_b = rng.normal(5.0, 0.2, (cells_per_blob, n_genes))
    return SingleCellData(X=np.vstack([blob_a, blob_b]))


class TestHierarchicalClusteringValidation:
    @pytest.mark.parametrize("linkage_method", _LINKAGE_METHODS)
    def test_every_supported_linkage_separates_the_blobs(self, linkage_method: str) -> None:
        result = hierarchical_clustering(_two_blob_data(), n_clusters=2, linkage_method=linkage_method)

        labels = result.obs["hierarchical_cluster"].to_numpy()
        assert len(set(labels[:15])) == 1
        assert len(set(labels[15:])) == 1
        assert labels[0] != labels[-1]
        assert result.uns["hierarchical_clustering"]["linkage_method"] == linkage_method

    def test_unknown_linkage_method_is_rejected_before_clustering(self) -> None:
        with pytest.raises(errors.ValidationError, match="Invalid linkage method: bogus"):
            hierarchical_clustering(_two_blob_data(), n_clusters=2, linkage_method="bogus")

    def test_unknown_metric_is_rejected_by_scipy(self) -> None:
        with pytest.raises(ValueError):
            hierarchical_clustering(_two_blob_data(), n_clusters=2, linkage_method="average", metric="bogus")

    def test_ward_linkage_always_uses_the_euclidean_metric(self) -> None:
        result = hierarchical_clustering(_two_blob_data(), n_clusters=2, linkage_method="ward", metric="cosine")

        assert result.uns["hierarchical_clustering"]["metric"] == "euclidean"


class TestDiffusionMapOperators:
    @staticmethod
    def _neighbors(n_cells: int, n_neighbors: int) -> SingleCellData:
        rng = np.random.default_rng(5)
        data = SingleCellData(X=rng.random((n_cells, 8)) + 0.1)
        return compute_neighbors(data, n_neighbors=n_neighbors)

    @pytest.mark.parametrize("alpha", [0.0, 0.5, 1.0])
    def test_each_normalisation_yields_finite_components(self, alpha: float) -> None:
        data = compute_diffusion_map(self._neighbors(40, 8), n_components=4, n_neighbors=8, alpha=alpha)

        components = data.obsm["X_diffmap"]
        assert components.shape == (40, 4)
        assert np.all(np.isfinite(components))
        eigenvalues = np.abs(data.uns["diffmap"]["eigenvalues"])
        assert np.all(np.diff(eigenvalues) <= 1e-10)  # largest magnitude first

    def test_dense_eigendecomposition_handles_a_matrix_too_small_for_the_sparse_solver(self) -> None:
        # With one cell the sparse solver rejects k = 0, so the dense symmetric solver
        # must receive a densified operator (not a scipy sparse matrix) and succeed.
        data = compute_diffusion_map(self._neighbors(1, 1), n_components=3, n_neighbors=1)

        assert data.obsm["X_diffmap"].shape[0] == 1
        assert np.all(np.isfinite(data.obsm["X_diffmap"]))
        assert "diffmap" in data.uns
