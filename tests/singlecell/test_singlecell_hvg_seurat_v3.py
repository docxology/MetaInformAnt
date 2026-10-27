"""``select_hvgs`` implements the ``seurat_v3`` flavor its docstring advertised.

The docstring listed ``'seurat_v3'`` among the flavors, but the function raised ``ValueError: Unknown HVG
selection method`` for it. The flavor is the variance-stabilising selection as scanpy defines it: raw counts, a
degree-2 local regression of log10 variance on log10 mean (span 0.3), clipping and standardisation with the
fitted standard deviation, and a ranking by the variance of the standardised counts.

scanpy cannot be the oracle here (its ``seurat_v3`` needs scikit-misc, which is not installed), so the pieces
are checked independently: the smoother against statsmodels' lowess (degree 1) and against a brute-force loop
(degrees 1 and 2) and its exact reproduction of quadratics; the statistic against an explicit per-gene loop of
the documented formula; and the selection against genes planted with real extra-Poisson variance. Real
computation only.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from metainformant.singlecell.analysis.pca_methods import select_hvgs
from metainformant.singlecell.data.preprocessing import SingleCellData

N_CELLS, N_GENES = 400, 600
SPAN = 0.3


def _brute_force_fit(x: np.ndarray, y: np.ndarray, span: float, degree: int) -> np.ndarray:
    """Local polynomial regression by an explicit loop: sort distances, tricube weights, weighted lstsq."""
    n = len(x)
    q = int(n * span + 1e-5)
    fitted = np.empty(n)
    for i in range(n):
        distance = np.abs(x - x[i])
        nearest = np.argsort(distance, kind="stable")[:q]
        width = distance[nearest].max()
        u = (x[nearest] - x[i]) / width
        weight = np.clip(1.0 - np.abs(u) ** 3, 0.0, None) ** 3
        design = np.vander(u, degree + 1, increasing=True)
        root = np.sqrt(weight)
        coefficients, *_ = np.linalg.lstsq(design * root[:, None], y[nearest] * root, rcond=None)
        fitted[i] = coefficients[0]
    return fitted


def _reference_statistic(counts: np.ndarray) -> np.ndarray:
    """The seurat_v3 statistic gene by gene.

    Clip each count at ``sd * sqrt(N) + mean``, then take ``sum((clipped - mean)^2) / ((N - 1) * sd^2)``.
    """
    n_cells, n_genes = counts.shape
    mean = counts.mean(axis=0)
    variance = np.array([counts[:, g].var(ddof=1) for g in range(n_genes)])
    varying = np.flatnonzero(variance > 0)
    expected_std = np.ones(n_genes)
    expected_std[varying] = np.sqrt(
        10.0 ** _brute_force_fit(np.log10(mean[varying]), np.log10(variance[varying]), SPAN, 2)
    )
    statistic = np.zeros(n_genes)
    for g in range(n_genes):
        clipped = np.minimum(counts[:, g], expected_std[g] * np.sqrt(n_cells) + mean[g])
        statistic[g] = np.sum((clipped - mean[g]) ** 2) / ((n_cells - 1) * expected_std[g] ** 2)
    return statistic


def _poisson_counts_with_planted_genes(seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Raw counts of Poisson genes (variance equals mean) over a wide range of means; 30 genes are planted.

    A planted gene switches between a high and a low rate in a quarter of the cells: a real extra-Poisson
    variance. Returns the counts and the indices of the planted genes.
    """
    rng = np.random.default_rng(seed)
    means = np.exp(rng.normal(0.0, 1.2, N_GENES))
    counts = rng.poisson(means, size=(N_CELLS, N_GENES)).astype(float)
    planted = np.arange(0, N_GENES, 20)
    state = rng.random(N_CELLS) < 0.25
    for gene in planted:
        counts[:, gene] = rng.poisson(np.where(state, means[gene] * 6.0, means[gene] * 0.4))
    return counts, planted


def _flagged(data: SingleCellData) -> np.ndarray:
    return data.var["highly_variable"].to_numpy()


class TestLocalPolynomialFit:
    @staticmethod
    def _points(n: int = 317, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
        rng = np.random.default_rng(seed)
        x = rng.normal(0.0, 0.7, n)
        y = 1.3 * x + 0.4 * x**2 + rng.normal(0.0, 0.25, n)
        return x, y

    def test_degree_one_is_the_statsmodels_lowess_fit(self) -> None:
        from statsmodels.nonparametric.smoothers_lowess import lowess

        from metainformant.singlecell.analysis.pca_methods import _local_polynomial_fit

        x, y = self._points()

        ours = _local_polynomial_fit(x, y, SPAN, degree=1)
        reference = lowess(y, x, frac=SPAN, it=0, delta=0.0, return_sorted=False)

        np.testing.assert_allclose(ours, reference, rtol=0, atol=1e-9)

    @pytest.mark.parametrize("degree", [1, 2])
    def test_it_matches_a_brute_force_loop(self, degree: int) -> None:
        from metainformant.singlecell.analysis.pca_methods import _local_polynomial_fit

        x, y = self._points()

        np.testing.assert_allclose(
            _local_polynomial_fit(x, y, SPAN, degree), _brute_force_fit(x, y, SPAN, degree), rtol=0, atol=1e-9
        )

    def test_a_quadratic_is_reproduced_exactly(self) -> None:
        from metainformant.singlecell.analysis.pca_methods import _local_polynomial_fit

        x, _ = self._points()
        exact = 0.5 - 2.0 * x + 1.5 * x**2

        np.testing.assert_allclose(_local_polynomial_fit(x, exact, SPAN, degree=2), exact, rtol=0, atol=1e-8)

    def test_the_order_of_the_points_does_not_matter(self) -> None:
        from metainformant.singlecell.analysis.pca_methods import _local_polynomial_fit

        x, y = self._points()
        order = np.random.default_rng(1).permutation(len(x))

        straight = _local_polynomial_fit(x, y, SPAN, degree=2)
        shuffled = _local_polynomial_fit(x[order], y[order], SPAN, degree=2)

        np.testing.assert_allclose(shuffled, straight[order], rtol=0, atol=1e-9)

    def test_more_repeated_values_than_a_neighbourhood_hold_are_refused(self) -> None:
        from metainformant.singlecell.analysis.pca_methods import _local_polynomial_fit

        x = np.concatenate([np.zeros(40), np.linspace(1.0, 2.0, 20)])  # 40 of 60 points share x = 0; q = 18

        with pytest.raises(ValueError, match="zero width"):
            _local_polynomial_fit(x, np.arange(60.0), SPAN, degree=2)

    def test_a_span_that_leaves_too_few_points_is_refused(self) -> None:
        from metainformant.singlecell.analysis.pca_methods import _local_polynomial_fit

        x16, y16 = self._points(n=16)
        x17, y17 = self._points(n=17)

        with pytest.raises(ValueError, match="at least 17 points are required"):
            _local_polynomial_fit(x16, y16, SPAN, degree=2)
        assert len(_local_polynomial_fit(x17, y17, SPAN, degree=2)) == 17  # int(17 * 0.3) = 5 points per fit


class TestStatistic:
    def test_it_is_the_documented_formula_gene_by_gene(self) -> None:
        from metainformant.singlecell.analysis.pca_methods import _seurat_v3_normalized_variance

        counts, _ = _poisson_counts_with_planted_genes()
        counts[:, 5] = 0.0  # an all-zero gene
        counts[:, 6] = 3.0  # a constant gene

        np.testing.assert_allclose(
            _seurat_v3_normalized_variance(counts), _reference_statistic(counts), rtol=1e-7, atol=1e-9
        )

    def test_genes_without_variance_score_zero(self) -> None:
        from metainformant.singlecell.analysis.pca_methods import _seurat_v3_normalized_variance

        counts, _ = _poisson_counts_with_planted_genes()
        counts[:, 5] = 0.0
        counts[:, 6] = 3.0

        statistic = _seurat_v3_normalized_variance(counts)

        assert statistic[5] == pytest.approx(0.0, abs=1e-12)
        assert statistic[6] == pytest.approx(0.0, abs=1e-12)

    def test_a_gene_that_follows_the_trend_scores_about_one(self) -> None:
        from metainformant.singlecell.analysis.pca_methods import _seurat_v3_normalized_variance

        counts, planted = _poisson_counts_with_planted_genes()

        statistic = _seurat_v3_normalized_variance(counts)

        ordinary = np.setdiff1d(np.arange(N_GENES), planted)
        assert np.median(statistic[ordinary]) == pytest.approx(1.0, abs=0.1)
        assert statistic[planted].min() > statistic[ordinary].max()  # the planted genes stand clear of the rest

    def test_a_single_extreme_count_is_clipped(self) -> None:
        from metainformant.singlecell.analysis.pca_methods import _seurat_v3_normalized_variance

        counts, _ = _poisson_counts_with_planted_genes()
        gene = 7
        counts[0, gene] = 500.0
        mean = counts[:, gene].mean()
        unclipped = np.sum((counts[:, gene] - mean) ** 2) / (N_CELLS - 1)

        statistic = _seurat_v3_normalized_variance(counts)

        assert statistic[gene] == pytest.approx(_reference_statistic(counts)[gene], rel=1e-7)
        assert unclipped > 600.0  # a variance of several hundred, from one cell
        assert statistic[gene] < 0.05 * unclipped  # clipped at about sd * sqrt(N) + mean = 20 or so

    def test_negative_values_are_refused(self) -> None:
        from metainformant.singlecell.analysis.pca_methods import _seurat_v3_normalized_variance

        counts, _ = _poisson_counts_with_planted_genes()
        counts[0, 0] = -1.0

        with pytest.raises(ValueError, match="raw counts"):
            _seurat_v3_normalized_variance(counts)

    def test_non_finite_values_are_refused(self) -> None:
        from metainformant.singlecell.analysis.pca_methods import _seurat_v3_normalized_variance

        counts, _ = _poisson_counts_with_planted_genes()
        counts[0, 0] = np.nan

        with pytest.raises(ValueError, match="raw counts"):
            _seurat_v3_normalized_variance(counts)

    def test_a_single_cell_is_refused(self) -> None:
        from metainformant.singlecell.analysis.pca_methods import _seurat_v3_normalized_variance

        with pytest.raises(ValueError, match="at least two cells"):
            _seurat_v3_normalized_variance(np.ones((1, 30)))

    def test_non_integer_values_warn_that_raw_counts_are_expected(self) -> None:
        from metainformant.singlecell.analysis.pca_methods import _seurat_v3_normalized_variance

        counts, _ = _poisson_counts_with_planted_genes()

        with pytest.warns(UserWarning, match="expects raw count data"):
            _seurat_v3_normalized_variance(np.log1p(counts))

    def test_many_genes_with_the_same_mean_are_refused(self) -> None:
        from metainformant.singlecell.analysis.pca_methods import _seurat_v3_normalized_variance

        rng = np.random.default_rng(2)
        counts = np.zeros((200, 100))
        counts[:, :40] = rng.poisson(np.exp(rng.normal(1.0, 1.0, 40)), size=(200, 40))
        for gene in range(40, 100):  # 60 genes with a single count in one cell: identical mean and variance
            counts[gene, gene] = 1.0

        with pytest.raises(ValueError, match="zero width"):
            _seurat_v3_normalized_variance(counts)


class TestSelectHvgsSeuratV3:
    def test_the_flavor_is_accepted_and_flags_the_requested_number_of_genes(self) -> None:
        counts, _ = _poisson_counts_with_planted_genes()

        result = select_hvgs(SingleCellData(counts), n_top_genes=30, flavor="seurat_v3")

        assert int(_flagged(result).sum()) == 30
        assert result.uns["hvg"] == {"flavor": "seurat_v3", "n_top_genes": 30, "n_selected": 30}

    def test_planted_extra_poisson_genes_are_found_where_plain_variance_misses_them(self) -> None:
        counts, planted = _poisson_counts_with_planted_genes()

        v3 = _flagged(select_hvgs(SingleCellData(counts), n_top_genes=30, flavor="seurat_v3"))
        plain = _flagged(select_hvgs(SingleCellData(counts), n_top_genes=30, flavor="variance"))

        # Poisson genes' variance grows with their mean, so plain variance prefers the highest-mean genes
        # (15 of the 30 planted genes here); the variance-stabilised statistic finds all of them.
        assert int(v3[planted].sum()) == len(planted)
        assert int(plain[planted].sum()) < int(v3[planted].sum())

    def test_variances_norm_holds_the_variance_stabilised_statistic(self) -> None:
        counts, _ = _poisson_counts_with_planted_genes()

        result = select_hvgs(SingleCellData(counts), n_top_genes=30, flavor="seurat_v3")

        np.testing.assert_allclose(
            result.var["variances_norm"].to_numpy(), _reference_statistic(counts), rtol=1e-7, atol=1e-9
        )
        assert {"highly_variable", "means", "variances", "variances_norm", "dispersions"} <= set(result.var.columns)

    def test_the_method_alias_gives_the_same_selection(self) -> None:
        counts, _ = _poisson_counts_with_planted_genes()

        by_flavor = select_hvgs(SingleCellData(counts), n_top_genes=30, flavor="seurat_v3")
        by_method = select_hvgs(SingleCellData(counts), n_top_genes=30, method="seurat_v3")

        np.testing.assert_array_equal(_flagged(by_flavor), _flagged(by_method))

    def test_sparse_counts_give_the_dense_selection(self) -> None:
        counts, _ = _poisson_counts_with_planted_genes()

        dense = select_hvgs(SingleCellData(counts), n_top_genes=30, flavor="seurat_v3")
        from_sparse = select_hvgs(SingleCellData(sparse.csr_matrix(counts)), n_top_genes=30, flavor="seurat_v3")

        np.testing.assert_array_equal(_flagged(from_sparse), _flagged(dense))

    def test_the_selection_does_not_depend_on_the_gene_order(self) -> None:
        counts, _ = _poisson_counts_with_planted_genes()
        order = np.random.default_rng(4).permutation(N_GENES)

        original = _flagged(select_hvgs(SingleCellData(counts), n_top_genes=30, flavor="seurat_v3"))
        shuffled = _flagged(select_hvgs(SingleCellData(counts[:, order]), n_top_genes=30, flavor="seurat_v3"))

        np.testing.assert_array_equal(shuffled, original[order])

    def test_the_input_is_left_alone_unless_inplace(self) -> None:
        counts, _ = _poisson_counts_with_planted_genes()
        data = SingleCellData(counts, var=pd.DataFrame({"gene_name": [f"g{i}" for i in range(N_GENES)]}))

        copied = select_hvgs(data, n_top_genes=30, flavor="seurat_v3")
        assert "highly_variable" not in data.var.columns
        assert copied is not data

        in_place = select_hvgs(data, n_top_genes=30, flavor="seurat_v3", inplace=True)
        assert in_place is data and "highly_variable" in data.var.columns

    def test_asking_for_more_genes_than_exist_flags_them_all(self) -> None:
        counts, _ = _poisson_counts_with_planted_genes()

        result = select_hvgs(SingleCellData(counts), n_top_genes=N_GENES + 100, flavor="seurat_v3")

        assert bool(_flagged(result).all())

    def test_genes_without_variance_are_flagged_last(self) -> None:
        counts, _ = _poisson_counts_with_planted_genes()
        counts[:, 5] = 0.0
        counts[:, 6] = 3.0
        varying = np.setdiff1d(np.arange(N_GENES), [5, 6])

        result = _flagged(select_hvgs(SingleCellData(counts), n_top_genes=len(varying), flavor="seurat_v3"))

        np.testing.assert_array_equal(np.flatnonzero(result), varying)

    def test_normalised_data_are_flagged_with_a_warning_not_silently(self) -> None:
        counts, _ = _poisson_counts_with_planted_genes()

        with pytest.warns(UserWarning, match="expects raw count data"):
            select_hvgs(SingleCellData(np.log1p(counts)), n_top_genes=30, flavor="seurat_v3")

    def test_an_unknown_flavor_still_raises_and_now_names_seurat_v3(self) -> None:
        counts, _ = _poisson_counts_with_planted_genes()

        with pytest.raises(ValueError, match="Unknown HVG selection method: seurat_v4.*'seurat_v3'"):
            select_hvgs(SingleCellData(counts), flavor="seurat_v4")
