"""``validate_deconvolution`` reports the correlation between estimated proportions and marker expression.

Its ``marker_correlations`` used to be ``min(1, 10 x the variance of the type's proportions)``: no expression
data reached the function, so the number was not a correlation and said nothing about whether the markers
agree with the estimate (it was never negative, and was the same for any markers). It is now the Pearson
(default) or Spearman correlation across spots between the estimated proportion of a cell type and the mean
expression of that type's markers, and the expected values below come from scipy. Real computation only: the
spots are real mixtures of reference profiles, deconvolved by the real NNLS routine.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse, stats

from metainformant.spatial.deconvolution.spatial_deconvolution import deconvolve_spots, validate_deconvolution
from metainformant.spatial.io.visium import SpatialDataset

TYPES = ["astrocyte", "microglia", "neuron"]
GENES = ["GFAP", "AQP4", "S100B", "CX3CR1", "P2RY12", "TMEM119", "RBFOX3", "SNAP25", "SYT1", "ACTB", "GAPDH", "B2M"]
MARKERS = {
    "astrocyte": ["GFAP", "AQP4", "S100B"],
    "microglia": ["CX3CR1", "P2RY12", "TMEM119"],
    "neuron": ["RBFOX3", "SNAP25", "SYT1"],
}


def _scenario(n_spots: int = 80, seed: int = 0) -> tuple[np.ndarray, dict[str, object]]:
    """Spots that mix three cell types; returns the spots x genes counts and the real deconvolution of them."""
    rng = np.random.default_rng(seed)
    reference = np.full((len(TYPES), len(GENES)), 1.0)  # a shared low level in every gene
    for k, cell_type in enumerate(TYPES):
        for gene in MARKERS[cell_type]:
            reference[k, GENES.index(gene)] = 20.0  # markers are high in their own type only
    mixture = rng.dirichlet([0.6, 0.6, 0.6], size=n_spots)
    counts = mixture @ reference + rng.normal(0.0, 0.8, size=(n_spots, len(GENES)))
    counts = np.clip(counts, 0.0, None)
    estimated = deconvolve_spots(counts, {t: reference[k].tolist() for k, t in enumerate(TYPES)})
    return counts, estimated


def _expected(estimated: dict[str, object], counts: np.ndarray, cell_type: str, method: str) -> float:
    """The correlation between a type's proportion and the mean of its markers, computed with scipy."""
    proportion = np.asarray(estimated["proportions_matrix"])[:, list(estimated["cell_types"]).index(cell_type)]
    signal = counts[:, [GENES.index(g) for g in MARKERS[cell_type]]].mean(axis=1)
    if method == "pearson":
        return float(stats.pearsonr(proportion, signal).statistic)
    return float(stats.spearmanr(proportion, signal).statistic)


class TestMarkerCorrelation:
    @pytest.mark.parametrize("method", ["pearson", "spearman"])
    def test_it_is_the_scipy_correlation_between_proportion_and_marker_mean(self, method: str) -> None:
        counts, estimated = _scenario()

        report = validate_deconvolution(
            estimated, MARKERS, expression=counts, gene_names=GENES, correlation_method=method  # type: ignore[arg-type]
        )

        assert set(report["marker_correlations"]) == set(TYPES)
        for cell_type in TYPES:
            assert report["marker_correlations"][cell_type] == pytest.approx(
                _expected(estimated, counts, cell_type, method), abs=1e-9
            )
        assert report["marker_correlation_method"] == method

    def test_pearson_is_the_default_and_the_two_methods_differ(self) -> None:
        counts, estimated = _scenario()

        default = validate_deconvolution(estimated, MARKERS, expression=counts, gene_names=GENES)
        pearson = validate_deconvolution(
            estimated, MARKERS, expression=counts, gene_names=GENES, correlation_method="pearson"
        )
        spearman = validate_deconvolution(
            estimated, MARKERS, expression=counts, gene_names=GENES, correlation_method="spearman"
        )

        assert default["marker_correlation_method"] == "pearson"
        assert default["marker_correlations"] == pearson["marker_correlations"]
        assert default["marker_correlations"] != spearman["marker_correlations"]

    def test_markers_that_agree_with_the_estimate_correlate_strongly(self) -> None:
        counts, estimated = _scenario()

        report = validate_deconvolution(estimated, MARKERS, expression=counts, gene_names=GENES)

        assert all(value > 0.9 for value in report["marker_correlations"].values())
        assert all(-1.0 <= value <= 1.0 for value in report["marker_correlations"].values())

    def test_the_markers_of_another_type_do_not_agree(self) -> None:
        counts, estimated = _scenario()
        swapped = {"astrocyte": MARKERS["neuron"], "microglia": MARKERS["astrocyte"], "neuron": MARKERS["microglia"]}

        right = validate_deconvolution(estimated, MARKERS, expression=counts, gene_names=GENES)
        wrong = validate_deconvolution(estimated, swapped, expression=counts, gene_names=GENES)

        for cell_type in TYPES:
            assert wrong["marker_correlations"][cell_type] < right["marker_correlations"][cell_type] - 0.8
            assert wrong["marker_correlations"][cell_type] < 0.0  # the old stand-in could never be negative

    def test_the_value_depends_on_the_expression_and_not_only_on_the_proportions(self) -> None:
        _, estimated = _scenario()
        proportion = np.asarray(estimated["proportions_matrix"])[:, TYPES.index("microglia")]
        markers = {"microglia": ["MARKER"]}
        aligned = (3.0 * proportion + 1.0).reshape(-1, 1)
        anti_aligned = (-3.0 * proportion + 9.0).reshape(-1, 1)

        up = validate_deconvolution(estimated, markers, expression=aligned, gene_names=["MARKER"])
        down = validate_deconvolution(estimated, markers, expression=anti_aligned, gene_names=["MARKER"])

        assert up["marker_correlations"]["microglia"] == pytest.approx(1.0)
        assert down["marker_correlations"]["microglia"] == pytest.approx(-1.0)

    def test_several_markers_are_averaged_before_correlating(self) -> None:
        counts, estimated = _scenario()
        only_one = {"astrocyte": ["GFAP"]}
        three = {"astrocyte": MARKERS["astrocyte"]}

        one_gene = validate_deconvolution(estimated, only_one, expression=counts, gene_names=GENES)
        three_genes = validate_deconvolution(estimated, three, expression=counts, gene_names=GENES)

        proportion = np.asarray(estimated["proportions_matrix"])[:, TYPES.index("astrocyte")]
        assert one_gene["marker_correlations"]["astrocyte"] == pytest.approx(
            stats.pearsonr(proportion, counts[:, GENES.index("GFAP")]).statistic
        )
        assert three_genes["marker_correlations"]["astrocyte"] != one_gene["marker_correlations"]["astrocyte"]


class TestExpressionInputs:
    def test_a_dataset_supplies_expression_and_gene_names(self) -> None:
        counts, estimated = _scenario()
        dataset = SpatialDataset(
            expression=counts,
            coordinates=np.zeros((len(counts), 2)),
            barcodes=[f"spot{i}" for i in range(len(counts))],
            gene_names=list(GENES),
        )

        from_dataset = validate_deconvolution(estimated, MARKERS, expression=dataset)
        from_matrix = validate_deconvolution(estimated, MARKERS, expression=counts, gene_names=GENES)

        assert from_dataset["marker_correlations"] == from_matrix["marker_correlations"]

    def test_sparse_and_list_expression_give_the_dense_result(self) -> None:
        counts, estimated = _scenario()
        dense = validate_deconvolution(estimated, MARKERS, expression=counts, gene_names=GENES)

        from_sparse = validate_deconvolution(estimated, MARKERS, expression=sparse.csr_matrix(counts), gene_names=GENES)
        from_lists = validate_deconvolution(estimated, MARKERS, expression=counts.tolist(), gene_names=GENES)

        assert from_sparse["marker_correlations"] == dense["marker_correlations"]
        assert from_lists["marker_correlations"] == dense["marker_correlations"]

    def test_a_dataframe_supplies_its_columns_as_gene_names(self) -> None:
        pandas = pytest.importorskip("pandas")
        counts, estimated = _scenario()

        from_frame = validate_deconvolution(estimated, MARKERS, expression=pandas.DataFrame(counts, columns=GENES))
        from_matrix = validate_deconvolution(estimated, MARKERS, expression=counts, gene_names=GENES)

        assert from_frame["marker_correlations"] == from_matrix["marker_correlations"]

    def test_column_order_does_not_matter(self) -> None:
        counts, estimated = _scenario()
        order = np.random.default_rng(5).permutation(len(GENES))

        shuffled = validate_deconvolution(
            estimated, MARKERS, expression=counts[:, order], gene_names=[GENES[i] for i in order]
        )
        original = validate_deconvolution(estimated, MARKERS, expression=counts, gene_names=GENES)

        for cell_type in TYPES:
            assert shuffled["marker_correlations"][cell_type] == pytest.approx(
                original["marker_correlations"][cell_type]
            )


class TestWhatIsLeftOut:
    def test_markers_missing_from_the_data_are_ignored_and_reported(self) -> None:
        counts, estimated = _scenario()
        markers = {"astrocyte": ["GFAP", "NOT_A_GENE", "AQP4"]}

        report = validate_deconvolution(estimated, markers, expression=counts, gene_names=GENES)

        assert report["marker_genes_used"] == {"astrocyte": ["GFAP", "AQP4"]}
        proportion = np.asarray(estimated["proportions_matrix"])[:, TYPES.index("astrocyte")]
        signal = counts[:, [GENES.index("GFAP"), GENES.index("AQP4")]].mean(axis=1)
        assert report["marker_correlations"]["astrocyte"] == pytest.approx(stats.pearsonr(proportion, signal).statistic)

    def test_a_type_without_any_marker_in_the_data_has_no_number_and_a_reason(self) -> None:
        counts, estimated = _scenario()

        report = validate_deconvolution(
            estimated, {"astrocyte": ["NOT_A_GENE"], "neuron": MARKERS["neuron"]}, expression=counts, gene_names=GENES
        )

        assert "astrocyte" not in report["marker_correlations"]
        assert "neuron" in report["marker_correlations"]
        assert "none of its markers" in report["marker_correlations_skipped"]["astrocyte"]

    def test_a_type_that_is_not_in_the_result_has_no_number_and_a_reason(self) -> None:
        counts, estimated = _scenario()

        report = validate_deconvolution(estimated, {"oligodendrocyte": ["GFAP"]}, expression=counts, gene_names=GENES)

        assert report["marker_correlations"] == {}
        assert "not in the deconvolution result" in report["marker_correlations_skipped"]["oligodendrocyte"]

    def test_constant_marker_expression_has_no_number_and_a_reason(self) -> None:
        counts, estimated = _scenario()
        flat = counts.copy()
        flat[:, GENES.index("GFAP")] = 5.0

        report = validate_deconvolution(estimated, {"astrocyte": ["GFAP"]}, expression=flat, gene_names=GENES)

        assert "astrocyte" not in report["marker_correlations"]
        assert "marker expression is constant" in report["marker_correlations_skipped"]["astrocyte"]

    def test_a_constant_estimated_proportion_has_no_number_and_a_reason(self) -> None:
        counts, estimated = _scenario()
        constant = dict(estimated)
        proportions = np.asarray(estimated["proportions_matrix"]).copy()
        proportions[:, TYPES.index("astrocyte")] = 0.25
        constant["proportions_matrix"] = proportions

        report = validate_deconvolution(constant, {"astrocyte": ["GFAP"]}, expression=counts, gene_names=GENES)

        assert "astrocyte" not in report["marker_correlations"]
        assert "proportion is constant" in report["marker_correlations_skipped"]["astrocyte"]

    def test_nothing_is_reported_without_markers(self) -> None:
        _, estimated = _scenario()

        report = validate_deconvolution(estimated)

        assert report["marker_correlations"] == {}
        assert report["marker_correlation_method"] is None
        assert report["marker_genes_used"] == {}
        assert report["marker_correlations_skipped"] == {}
        assert report["normalization_check"] is True
        assert report["type_diversity"] >= 2


class TestWhatIsRefused:
    def test_markers_without_expression_raise(self) -> None:
        _, estimated = _scenario()

        with pytest.raises(ValueError, match="needs the spatial expression matrix"):
            validate_deconvolution(estimated, MARKERS)

    def test_a_plain_matrix_without_gene_names_raises(self) -> None:
        counts, estimated = _scenario()

        with pytest.raises(ValueError, match="gene_names is required"):
            validate_deconvolution(estimated, MARKERS, expression=counts)

    def test_an_unknown_correlation_method_raises(self) -> None:
        counts, estimated = _scenario()

        with pytest.raises(ValueError, match="Unknown correlation_method"):
            validate_deconvolution(
                estimated,
                MARKERS,
                expression=counts,
                gene_names=GENES,
                correlation_method="kendall",  # type: ignore[arg-type]
            )

    def test_expression_must_have_one_row_per_spot(self) -> None:
        counts, estimated = _scenario()

        with pytest.raises(ValueError, match="expression has 79 rows but the deconvolution has 80 spots"):
            validate_deconvolution(estimated, MARKERS, expression=counts[:-1], gene_names=GENES)

    def test_gene_names_must_name_every_column(self) -> None:
        counts, estimated = _scenario()

        with pytest.raises(ValueError, match="gene_names has 11 entries but the expression matrix has 12 columns"):
            validate_deconvolution(estimated, MARKERS, expression=counts, gene_names=GENES[:-1])

    def test_a_used_marker_that_labels_two_columns_is_ambiguous(self) -> None:
        counts, estimated = _scenario()
        names = list(GENES)
        names[GENES.index("ACTB")] = "GFAP"

        with pytest.raises(ValueError, match="more than one expression column"):
            validate_deconvolution(estimated, MARKERS, expression=counts, gene_names=names)
