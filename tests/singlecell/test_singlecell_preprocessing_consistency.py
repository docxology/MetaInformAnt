"""Consistency of filter_cells / filter_genes / normalize_counts bookkeeping.

- filter_cells recomputes the per-gene QC columns on the cells that remain;
- filter_cells / filter_genes (and SingleCellData.copy) never share nested ``uns`` dictionaries with their input;
- normalize_counts stores size factors that point the same way for every method.

Real data and real library calls only.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from metainformant.singlecell.data.preprocessing import (
    SingleCellData,
    calculate_qc_metrics,
    filter_cells,
    filter_genes,
    normalize_counts,
)

# Cells c2 and c3 have 2 and 1 counts in total; the other four carry 17, 11, 13 and 10.
COUNTS = np.array(
    [
        [10, 0, 2, 0],
        [8, 0, 3, 0],
        [0, 1, 0, 1],
        [0, 1, 0, 0],
        [9, 0, 4, 0],
        [7, 1, 2, 0],
    ],
    dtype=float,
)


def _data(matrix: np.ndarray | sparse.spmatrix = COUNTS) -> SingleCellData:
    var = pd.DataFrame(index=["g0", "g1", "g2", "g3"])
    obs = pd.DataFrame(index=[f"c{i}" for i in range(COUNTS.shape[0])])
    return calculate_qc_metrics(SingleCellData(matrix, obs=obs, var=var))


@pytest.fixture(params=["dense", "sparse"])
def counts_data(request: pytest.FixtureRequest) -> SingleCellData:
    return _data(COUNTS if request.param == "dense" else sparse.csr_matrix(COUNTS))


class TestFilterCellsRefreshesGeneMetrics:
    def test_the_metrics_before_filtering_describe_all_six_cells(self, counts_data: SingleCellData) -> None:
        assert counts_data.var["n_cells"].tolist() == [4, 3, 4, 1]

    def test_per_gene_metrics_describe_the_cells_that_remain(self, counts_data: SingleCellData) -> None:
        kept = filter_cells(counts_data, min_counts=3)

        assert kept.obs.index.tolist() == ["c0", "c1", "c4", "c5"]
        assert kept.var["n_cells"].tolist() == [4, 1, 4, 0]
        assert kept.var["total_counts"].tolist() == [34, 1, 11, 0]
        assert kept.var["mean_expression"].tolist() == pytest.approx([8.5, 0.25, 2.75, 0.0])
        assert kept.var["pct_dropout"].tolist() == pytest.approx([0.0, 75.0, 0.0, 100.0])

    def test_a_later_gene_filter_sees_the_refreshed_counts(self, counts_data: SingleCellData) -> None:
        kept = filter_cells(counts_data, min_counts=3)

        genes = filter_genes(kept, min_cells=1)

        assert genes.var.index.tolist() == ["g0", "g1", "g2"]  # g3 is not expressed in any remaining cell

    def test_columns_that_were_never_computed_are_not_invented(self) -> None:
        bare = SingleCellData(COUNTS.copy(), obs=pd.DataFrame({"total_counts": COUNTS.sum(axis=1)}))

        kept = filter_cells(bare, min_counts=3)

        assert "n_cells" not in kept.var.columns  # filter_cells never computed gene metrics for this input

    def test_removing_every_cell_leaves_nan_ratios_without_warnings(self, counts_data: SingleCellData) -> None:
        kept = filter_cells(counts_data, min_counts=1000)  # warnings are errors under this suite's config

        assert kept.n_obs == 0
        assert kept.var["n_cells"].tolist() == [0, 0, 0, 0]
        assert kept.var["mean_expression"].isna().all()


class TestUnsIsNotShared:
    def test_filter_cells_leaves_the_input_qc_summary_alone(self, counts_data: SingleCellData) -> None:
        summary_before = dict(counts_data.uns["qc_summary"])

        kept = filter_cells(counts_data, min_counts=3)

        assert kept.uns["qc_summary"]["cells_after_filtering"] == 4
        assert counts_data.uns["qc_summary"] == summary_before
        assert "cells_after_filtering" not in counts_data.uns["qc_summary"]

    def test_filter_genes_leaves_the_input_qc_summary_alone(self, counts_data: SingleCellData) -> None:
        kept = filter_genes(counts_data, min_cells=2)

        assert kept.uns["qc_summary"]["genes_after_filtering"] == 3
        assert "genes_after_filtering" not in counts_data.uns["qc_summary"]

    def test_copy_does_not_share_nested_uns(self, counts_data: SingleCellData) -> None:
        duplicate = counts_data.copy()
        duplicate.uns["qc_summary"]["total_cells"] = -1

        assert counts_data.uns["qc_summary"]["total_cells"] == 6


class TestSizeFactorDirection:
    @pytest.fixture()
    def counts(self) -> np.ndarray:
        rng = np.random.default_rng(2)
        return rng.poisson(lam=rng.uniform(1.0, 20.0, size=(30, 1)) * rng.uniform(0.5, 2.0, size=(1, 15))) + 1.0

    @pytest.mark.parametrize("method", ["total_count", "total", "median", "size_factors"])
    def test_every_method_normalises_as_x_over_size_factor_times_target(self, counts: np.ndarray, method: str) -> None:
        result = normalize_counts(SingleCellData(counts.copy()), method=method)

        info = result.uns["normalization"]
        size_factors = np.asarray(info["size_factors"])
        expected = counts / size_factors[:, None] * info["target_sum"]
        np.testing.assert_allclose(result.X, expected, rtol=1e-10)

    @pytest.mark.parametrize("method", ["total_count", "median"])
    def test_library_size_methods_store_the_library_size(self, counts: np.ndarray, method: str) -> None:
        result = normalize_counts(SingleCellData(counts.copy()), method=method)

        np.testing.assert_allclose(result.uns["normalization"]["size_factors"], counts.sum(axis=1))

    @pytest.mark.parametrize("method", ["total_count", "median"])
    def test_a_larger_library_has_a_larger_size_factor_with_every_method(self, method: str) -> None:
        counts = np.array([[1.0, 1.0], [4.0, 4.0], [9.0, 9.0]])

        result = normalize_counts(SingleCellData(counts), method=method)

        size_factors = result.uns["normalization"]["size_factors"]
        assert size_factors == sorted(size_factors)
        assert size_factors[0] < size_factors[-1]

    def test_all_methods_agree_on_the_direction(self, counts: np.ndarray) -> None:
        depth = counts.sum(axis=1)
        for method in ("total_count", "median", "size_factors"):
            factors = normalize_counts(SingleCellData(counts.copy()), method=method).uns["normalization"][
                "size_factors"
            ]
            assert np.corrcoef(factors, depth)[0, 1] > 0.9, method

    @pytest.mark.parametrize("method", ["total_count", "median"])
    def test_an_empty_cell_stays_empty_with_a_zero_size_factor(self, method: str) -> None:
        counts = np.array([[3.0, 1.0], [0.0, 0.0], [2.0, 6.0]])

        result = normalize_counts(SingleCellData(counts), method=method)

        assert result.uns["normalization"]["size_factors"][1] == 0.0
        assert np.all(result.X[1] == 0.0)
        assert np.all(np.isfinite(result.X))
