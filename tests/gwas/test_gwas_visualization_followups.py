"""Regression tests for the follow-up GWAS plotting fixes.

Each class covers one defect found while fixing the first round:

* ``functional_enrichment_plot`` invented categories (``POS % 8``) and a forced "intergenic" count;
* ``calculate_population_enrichment`` carried an invented ``p = 0.5`` fallback and crashed on degenerate tables;
* ``analyze_genetic_architecture`` reported lambda / (lambda + 1) as "heritability";
* ``interactive_pca`` rejected the NumPy array its docstring accepts;
* the regional plots and ``variant_density_plot`` could not read the pipeline's own ``chrom`` / ``pos`` records;
* ``generate_all_plots`` left every figure it returned open;
* ``effect_size_manhattan`` used the row rank as x; ``plot_gene_track`` ignored the chromosome;
* ``hwe_deviation_plot`` drew random numbers as the expected values of its Q-Q panel;
* ``recombination_rate_plot`` raised ``ValueError`` from its smoothing line for 11 to 19 positions.

As elsewhere, assertions read the artists (bar widths, offsets, labels, text) or the returned metadata.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PathCollection  # noqa: E402

from metainformant.gwas.visualization import utils as viz_utils  # noqa: E402
from metainformant.gwas.visualization.genomic import regional, variants  # noqa: E402
from metainformant.gwas.visualization.interactive import interactive, suite  # noqa: E402
from metainformant.gwas.visualization.statistical import comparison, effects  # noqa: E402


def _offsets(collection) -> np.ndarray:
    return np.asarray(collection.get_offsets(), dtype=float)


def _scatters(ax) -> list[PathCollection]:
    return [c for c in ax.collections if isinstance(c, PathCollection)]


# ---------------------------------------------------------------------------
# 1. functional_enrichment_plot: real annotations only
# ---------------------------------------------------------------------------


def _enrichment_variants() -> list[dict]:
    return [
        {"CHROM": "1", "POS": 1000, "ID": "rs1", "p_value": 1e-10},
        {"CHROM": "1", "POS": 2000, "ID": "rs2", "p_value": 1e-12},
        {"CHROM": "2", "POS": 5000, "ID": "rs3", "p_value": 1e-8},
        {"CHROM": "2", "POS": 6000, "ID": "rs4", "p_value": 0.01},  # not significant
    ]


def _with_field(categories: dict, key: str = "annotation") -> list[dict]:
    return [
        {**variant, key: categories[variant["ID"]]} for variant in _enrichment_variants() if variant["ID"] in categories
    ]


class TestFunctionalEnrichmentAnnotations:
    CATEGORIES = {"rs1": "missense", "rs2": "intronic", "rs3": "missense", "rs4": "exonic"}

    def test_field_annotations_are_counted_as_given(self, tmp_path):
        result = effects.functional_enrichment_plot(_with_field(self.CATEGORIES), tmp_path / "f.png")

        assert result["status"] == "success"
        assert result["category_counts"] == {"missense": 2, "intronic": 1}  # rs4 is not significant
        assert result["total_significant_variants"] == 3
        assert result["categories_found"] == 2
        assert result["n_unannotated"] == 0
        assert (tmp_path / "f.png").stat().st_size > 0

    def test_no_category_is_invented(self):
        rows = _with_field({"rs1": "class_x", "rs2": "class_y", "rs3": "class_x"})

        result = effects.functional_enrichment_plot(rows)

        assert set(result["category_counts"]) == {"class_x", "class_y"}  # no "intergenic", "exonic", ...

    def test_bars_show_the_counts_largest_first(self):
        effects.functional_enrichment_plot(_with_field(self.CATEGORIES))

        ax = plt.gcf().axes[0]
        assert [t.get_text() for t in ax.get_yticklabels()] == ["missense", "intronic"]
        assert [p.get_width() for p in ax.patches] == [2, 1]

    def test_mapping_by_variant_id(self):
        result = effects.functional_enrichment_plot(_enrichment_variants(), annotations=self.CATEGORIES)

        assert result["category_counts"] == {"missense": 2, "intronic": 1}

    def test_mapping_by_chromosome_and_position(self):
        by_position = {("1", 1000): "a", (1, 2000): "b", ("2", 5000): "a", ("2", 6000): "c"}

        result = effects.functional_enrichment_plot(_enrichment_variants(), annotations=by_position)

        assert result["category_counts"] == {"a": 2, "b": 1}

    def test_annotations_may_be_a_pandas_series(self):
        series = pd.Series(self.CATEGORIES)  # variant id -> category; a Series has no unambiguous truth value

        result = effects.functional_enrichment_plot(_enrichment_variants(), annotations=series)

        assert result["category_counts"] == {"missense": 2, "intronic": 1}

    def test_a_nan_position_leaves_that_variant_unannotated_instead_of_aborting(self):
        rows = _enrichment_variants()[:3]
        rows[0] = {"CHROM": "1", "POS": float("nan"), "p_value": 1e-10}  # no id, no field, no usable position

        result = effects.functional_enrichment_plot(rows, annotations={("1", 2000): "intronic", ("2", 5000): "exonic"})

        assert result["status"] == "success"
        assert result["category_counts"] == {"intronic": 1, "exonic": 1}
        assert result["n_unannotated"] == 1

    def test_custom_field_name(self):
        rows = _with_field(self.CATEGORIES, key="consequence")

        result = effects.functional_enrichment_plot(rows, annotation_key="consequence")

        assert result["category_counts"] == {"missense": 2, "intronic": 1}

    def test_a_field_on_the_variant_wins_over_the_mapping(self):
        rows = _with_field({"rs1": "from_field"})  # only rs1 carries a field

        result = effects.functional_enrichment_plot(rows + _enrichment_variants()[1:], annotations=self.CATEGORIES)

        assert result["category_counts"] == {"from_field": 1, "intronic": 1, "missense": 1}

    def test_unannotated_significant_variants_are_reported_not_assigned(self):
        rows = _with_field({"rs2": "intronic", "rs3": "missense"}) + [_enrichment_variants()[0]]  # rs1 unannotated

        result = effects.functional_enrichment_plot(rows)

        assert result["category_counts"] == {"intronic": 1, "missense": 1}
        assert result["total_significant_variants"] == 3
        assert result["n_unannotated"] == 1
        notes = [t.get_text() for t in plt.gcf().axes[0].texts]
        assert any("1 significant variant" in note and "annotation" in note for note in notes)

    def test_blank_and_nan_annotations_count_as_missing(self):
        rows = [
            {"CHROM": "1", "POS": 1, "p_value": 1e-10, "annotation": "missense"},
            {"CHROM": "1", "POS": 2, "p_value": 1e-10, "annotation": float("nan")},
            {"CHROM": "1", "POS": 3, "p_value": 1e-10, "annotation": "  "},
            {"CHROM": "1", "POS": 4, "p_value": 1e-10, "annotation": None},
        ]

        result = effects.functional_enrichment_plot(rows)

        assert result["category_counts"] == {"missense": 1}
        assert result["n_unannotated"] == 3

    def test_without_any_annotation_the_result_is_unavailable_not_invented(self, tmp_path):
        out = tmp_path / "f.png"

        result = effects.functional_enrichment_plot(_enrichment_variants(), out)

        assert result["status"] == "skipped"
        assert "annotation" in result["message"]
        assert result["total_significant_variants"] == 3
        assert result["n_unannotated"] == 3
        assert "category_counts" not in result
        assert not out.exists()

    def test_equal_counts_keep_first_seen_order(self):
        rows = _with_field({"rs3": "beta_class", "rs1": "alpha_class"})
        rows = sorted(rows, key=lambda v: v["ID"], reverse=True)  # rs3 first

        result = effects.functional_enrichment_plot(rows)

        assert list(result["category_counts"]) == ["beta_class", "alpha_class"]
        assert [t.get_text() for t in plt.gcf().axes[0].get_yticklabels()] == ["beta_class", "alpha_class"]

    def test_top_n_limits_the_bars(self):
        rows = [
            {"CHROM": "1", "POS": i, "p_value": 1e-10, "annotation": f"cat{i}"} for i in range(1, 6) for _ in range(i)
        ]

        result = effects.functional_enrichment_plot(rows, top_n=2)

        assert result["categories_found"] == 5
        assert [p.get_width() for p in plt.gcf().axes[0].patches] == [5, 4]

    def test_dictionary_of_enrichment_statistics_is_unchanged(self, tmp_path):
        stats = {"missense": {"enrichment_ratio": 3.1, "p_value": 1e-5}, "synonymous": {"enrichment_ratio": 0.5}}

        result = effects.functional_enrichment_plot(stats, tmp_path / "f.png")

        assert result["status"] == "success"
        assert result["top_functions"] == ["missense"]


# ---------------------------------------------------------------------------
# 2. calculate_population_enrichment: no invented p-value
# ---------------------------------------------------------------------------


def _p_frame(n_significant: int, n_total: int) -> pd.DataFrame:
    return pd.DataFrame({"P": [1e-9] * n_significant + list(np.linspace(0.05, 0.95, n_total - n_significant))})


class TestPopulationEnrichmentPValue:
    def test_p_value_equals_scipys_chi_square_on_random_tables(self):
        scipy_stats = pytest.importorskip("scipy.stats")
        rng = np.random.default_rng(0)
        checked = 0
        for _ in range(300):
            table = rng.integers(0, 600, size=(2, 2)).tolist()
            if min(map(sum, table)) == 0 or min(map(sum, zip(*table))) == 0:
                continue
            expected = scipy_stats.chi2_contingency(table)[1]
            assert comparison._chi_square_p_value_2x2(table) == pytest.approx(expected, rel=1e-9, abs=1e-300)
            checked += 1
        assert checked > 250

    @pytest.mark.parametrize("table", [[[1, 2], [2, 1]], [[0, 5], [5, 0]], [[3, 3], [3, 3]], [[40, 60], [55, 45]]])
    def test_small_and_symmetric_tables_match_scipy(self, table):
        scipy_stats = pytest.importorskip("scipy.stats")

        assert comparison._chi_square_p_value_2x2(table) == pytest.approx(scipy_stats.chi2_contingency(table)[1])

    @pytest.mark.parametrize("table", [[[0, 5], [0, 7]], [[0, 0], [3, 4]], [[0, 0], [0, 0]], [[4, 0], [6, 0]]])
    def test_a_degenerate_table_has_no_p_value(self, table):
        assert math.isnan(comparison._chi_square_p_value_2x2(table))

    def test_population_results_match_scipy(self):
        scipy_stats = pytest.importorskip("scipy.stats")
        pops = {"EUR": {"results": _p_frame(10, 200)}, "AFR": {"results": _p_frame(0, 100)}}

        out = comparison.calculate_population_enrichment(pops)

        for name, (n_sig, n_total) in {"EUR": (10, 200), "AFR": (0, 100)}.items():
            table = [[n_sig, n_total - n_sig], [int(0.5 * n_total), int(0.5 * n_total)]]
            assert out[name]["p_value"] == pytest.approx(scipy_stats.chi2_contingency(table)[1], rel=1e-9)
            assert out[name]["enrichment"] == pytest.approx((n_sig / n_total) / 0.5)
            assert out[name]["observed_sig"] == n_sig
            assert out[name]["expected_sig"] == int(0.5 * n_total)

    def test_a_lone_population_without_hits_gets_nan_instead_of_a_crash(self):
        out = comparison.calculate_population_enrichment({"only": {"results": _p_frame(0, 50)}})

        assert out["only"]["enrichment"] == 1.0
        assert math.isnan(out["only"]["p_value"])  # never an invented 0.5

    def test_the_comparison_figure_marks_an_undefined_p_value(self):
        pytest.importorskip("seaborn")

        fig = comparison.compare_populations({"only": {"results": _p_frame(0, 50)}})

        labels = [t.get_text() for t in fig.axes[2].texts]
        assert labels == ["p=n/a"]

    def test_the_comparison_figure_prints_defined_p_values(self):
        pytest.importorskip("seaborn")
        pops = {"EUR": {"results": _p_frame(10, 200)}, "AFR": {"results": _p_frame(0, 100)}}

        fig = comparison.compare_populations(pops)

        labels = [t.get_text() for t in fig.axes[2].texts]
        assert len(labels) == 2 and all(label.startswith("p=") and label != "p=n/a" for label in labels)


# ---------------------------------------------------------------------------
# 3. analyze_genetic_architecture: no "heritability"
# ---------------------------------------------------------------------------


class TestGeneticArchitecture:
    FRAME = pd.DataFrame({"P": [0.5] * 100 + [1e-9]})

    def test_no_heritability_is_reported(self):
        architecture = comparison.analyze_genetic_architecture({"PopA": {"results": self.FRAME}})

        assert set(architecture["PopA"]) == {"polygenicity", "lambda_gc"}
        assert "heritability" not in architecture["PopA"]
        assert architecture["PopA"]["lambda_gc"] == pytest.approx(1.0, abs=0.01)
        assert architecture["PopA"]["polygenicity"] == pytest.approx(1 / 101)

    def test_the_comparison_figure_has_one_honestly_labelled_series(self):
        pytest.importorskip("seaborn")
        pops = {"A": {"results": _p_frame(2, 100)}, "B": {"results": _p_frame(5, 100)}}

        fig = comparison.compare_populations(pops)

        panel = fig.axes[3]
        assert [p.get_height() for p in panel.patches] == pytest.approx([0.02, 0.05])
        assert panel.get_legend() is None
        assert "heritab" not in (panel.get_title() + panel.get_ylabel()).lower()
        assert "significant" in panel.get_ylabel().lower()
        assert [t.get_text() for t in panel.get_xticklabels()] == ["A", "B"]


# ---------------------------------------------------------------------------
# 5. interactive_pca accepts arrays
# ---------------------------------------------------------------------------


def _association_tsv(path: Path, n: int = 30, seed: int = 1) -> Path:
    rng = np.random.default_rng(seed)
    lines = ["CHROM\tPOS\tp_value\tBETA\tSE"]
    for i in range(n):
        lines.append(f"{1 + i % 2}\t{1000 * (i + 1)}\t{rng.uniform(1e-6, 1):.6g}\t{rng.normal():.4f}\t0.1")
    path.write_text("\n".join(lines) + "\n")
    return path


class TestInteractivePcaArrays:
    PCS = np.random.default_rng(0).normal(size=(8, 3))

    def test_a_numpy_array_is_accepted(self, tmp_path):
        out = tmp_path / "pca.html"

        result = interactive.interactive_pca(
            {"pcs": self.PCS, "explained_variance_ratio": [0.5, 0.3, 0.1]}, output_file=out
        )

        assert result["status"] == "success"
        assert result["n_samples"] == 8
        assert out.stat().st_size > 0

    def test_lists_still_work(self, tmp_path):
        result = interactive.interactive_pca({"pcs": self.PCS.tolist()}, output_file=tmp_path / "pca.html")

        assert (result["status"], result["n_samples"]) == ("success", 8)

    @pytest.mark.parametrize("empty", [[], np.empty((0, 3)), None])
    def test_empty_input_is_skipped(self, empty):
        result = interactive.interactive_pca({"pcs": empty})

        assert result["status"] == "skipped"
        assert result["n_samples"] == 0

    def test_a_one_dimensional_array_is_a_failure_not_an_exception(self):
        assert interactive.interactive_pca({"pcs": np.arange(5.0)})["status"] == "failed"

    def test_fewer_than_three_components_fail(self):
        result = interactive.interactive_pca({"pcs": self.PCS[:, :2]})

        assert (result["status"], result["n_samples"]) == ("failed", 8)

    def test_the_suite_draws_the_interactive_pca_from_an_in_memory_array(self, tmp_path):
        result = suite.generate_all_plots(
            _association_tsv(tmp_path / "assoc.tsv"), tmp_path / "plots", pca_data={"pcs": self.PCS}
        )

        assert result["status"] == "completed"
        assert Path(result["plots"]["pca_interactive"]).stat().st_size > 0

    def test_the_file_loader_hands_over_an_array(self, tmp_path):
        path = tmp_path / "pca.tsv"
        np.savetxt(path, self.PCS, delimiter="\t")

        loaded = suite._load_pca_file(path)

        assert isinstance(loaded["pcs"], np.ndarray)
        np.testing.assert_allclose(loaded["pcs"], self.PCS)


# ---------------------------------------------------------------------------
# 6. regional plots and variant density read the pipeline's own record keys
# ---------------------------------------------------------------------------

PIPELINE_RECORDS = [
    {"chrom": "1", "pos": 1000, "p_value": 1e-8},
    {"chrom": "1", "pos": 1500, "p_value": 1e-6},
    {"chrom": "1", "pos": 2000, "p_value": 0.01},
    {"chrom": "2", "pos": 1000, "p_value": 0.5},
]


class TestPipelineRecordKeys:
    def test_regional_plot_reads_chrom_pos_p_value(self, tmp_path):
        result = regional.regional_plot(PIPELINE_RECORDS, tmp_path / "r.png", chrom="1", start=500, end=2500)

        assert result["status"] == "success"
        assert result["n_variants"] == 3

    def test_gene_annotation_plot_reads_chrom_pos_p_value(self, tmp_path):
        result = regional.gene_annotation_plot(PIPELINE_RECORDS, tmp_path / "g.png", chrom="1", start=500, end=2500)

        assert result["status"] == "success"
        assert result["n_variants"] == 3

    @pytest.mark.parametrize("plot", ["regional_plot", "gene_annotation_plot"])
    def test_dataframes_with_the_pipelines_columns(self, plot):
        frame = pd.DataFrame(PIPELINE_RECORDS)

        result = getattr(regional, plot)(frame, chrom="1", start=500, end=2500)

        assert result["n_variants"] == 3

    @pytest.mark.parametrize("key", ["chromosome", "chrom", "chr", "CHROM", "CHR"])
    @pytest.mark.parametrize("plot", ["regional_plot", "gene_annotation_plot"])
    def test_every_chromosome_spelling(self, plot, key):
        rows = [{key: "1", "POS": 1000, "p_value": 1e-8}, {key: "2", "POS": 1000, "p_value": 1e-3}]

        result = getattr(regional, plot)(rows, chrom="1", start=500, end=2500)

        assert result["n_variants"] == 1

    @pytest.mark.parametrize("key", ["position", "pos", "POS", "BP"])
    def test_every_position_spelling(self, key):
        rows = [{"CHROM": "1", key: 1000, "p_value": 1e-8}, {"CHROM": "1", key: 9000, "p_value": 1e-3}]

        result = regional.regional_plot(rows, chrom="1", start=500, end=2500)

        assert result["n_variants"] == 1

    @pytest.mark.parametrize("key", ["p_value", "pvalue", "pval", "P"])
    def test_every_p_value_spelling(self, key):
        rows = [{"CHROM": "1", "POS": 1000, key: 1e-8}, {"CHROM": "1", "POS": 1200, key: 1e-3}]

        result = regional.regional_plot(rows, chrom="1", start=500, end=2500)

        assert result["n_variants"] == 2

    def test_standardize_renames_the_first_alias_of_each_family(self):
        frame = pd.DataFrame({"chrom": ["1"], "pos": [5], "pvalue": [0.1], "other": [1]})

        out = viz_utils.standardize_result_columns(frame)

        assert list(out.columns) == ["CHR", "BP", "P", "other"]
        assert list(frame.columns) == ["chrom", "pos", "pvalue", "other"]  # the input is not modified

    def test_standardize_never_overwrites_a_canonical_column(self):
        frame = pd.DataFrame({"CHR": ["7"], "chrom": ["1"], "BP": [5], "pos": [9], "P": [0.1], "p_value": [0.9]})

        out = viz_utils.standardize_result_columns(frame)

        assert out["CHR"].tolist() == ["7"] and out["BP"].tolist() == [5] and out["P"].tolist() == [0.1]
        assert {"chrom", "pos", "p_value"} <= set(out.columns)  # the unused aliases are left alone

    def test_alias_tuples_cover_the_interactive_manhattans(self):
        assert set(interactive._POSITION_KEYS) <= set(viz_utils.POSITION_KEYS)
        assert set(interactive._P_VALUE_KEYS) <= set(viz_utils.P_VALUE_KEYS)

    def test_get_position(self):
        for key in ("position", "pos", "POS", "BP"):
            assert viz_utils.get_position({key: 42}) == 42
        assert viz_utils.get_position({"pos": None, "POS": 7}) == 7
        assert viz_utils.get_position({}) is None
        assert viz_utils.get_position({}, 0) == 0


class TestVariantDensityAliases:
    @pytest.mark.parametrize(
        ("chrom_key", "pos_key"),
        [("CHROM", "POS"), ("chrom", "pos"), ("chromosome", "position"), ("chr", "BP"), ("CHR", "POS")],
    )
    def test_every_spelling_is_read(self, chrom_key, pos_key):
        rows = [{chrom_key: "1", pos_key: 100}, {chrom_key: "1", pos_key: 200}, {chrom_key: "2", pos_key: 50}]

        result = variants.variant_density_plot(rows)

        assert result["status"] == "success"
        assert result["n_chromosomes"] == 2
        assert result["n_variants"] == 3

    def test_numeric_chromosomes_that_are_floats_group_with_their_integers(self):
        rows = [{"CHROM": 1, "POS": 100}, {"CHROM": 1.0, "POS": 200}, {"CHROM": "1", "POS": 300}]

        assert variants.variant_density_plot(rows)["n_chromosomes"] == 1


# ---------------------------------------------------------------------------
# 7. generate_all_plots closes the figures it opens
# ---------------------------------------------------------------------------


class TestSuiteClosesItsFigures:
    def test_no_pyplot_figure_is_left_open(self, tmp_path):
        result = suite.generate_all_plots(_association_tsv(tmp_path / "assoc.tsv"), tmp_path / "plots")

        assert result["status"] == "completed"
        assert plt.get_fignums() == []

    def test_files_are_saved_before_the_figures_are_closed(self, tmp_path):
        rng = np.random.default_rng(2)
        pcs = rng.normal(size=(10, 3))
        np.savetxt(tmp_path / "pca.tsv", pcs, delimiter="\t")
        np.savetxt(tmp_path / "kinship.tsv", np.corrcoef(rng.normal(size=(10, 30))), delimiter="\t")

        result = suite.generate_all_plots(
            _association_tsv(tmp_path / "assoc.tsv"),
            tmp_path / "plots",
            pca_file=tmp_path / "pca.tsv",
            kinship_file=tmp_path / "kinship.tsv",
        )

        for name in ("manhattan", "qq", "regional_1", "regional_5", "effect_size", "pca", "kinship"):
            assert Path(result["plots"][name]).stat().st_size > 0, name
        assert plt.get_fignums() == []

    def test_repeated_calls_neither_accumulate_figures_nor_lose_plots(self, tmp_path):
        # without closing, the 21st open figure turns into an error under this repository's warning filter and
        # later plots (the summary panel among them) silently disappear
        keys = []
        for run in range(4):
            result = suite.generate_all_plots(_association_tsv(tmp_path / f"assoc{run}.tsv"), tmp_path / f"plots{run}")
            assert result["status"] == "completed"
            keys.append(sorted(result["plots"]))

        assert all(k == keys[0] for k in keys)
        assert "summary_panel" in keys[0]
        assert plt.get_fignums() == []

    def test_a_call_with_every_optional_input_leaves_nothing_open(self, tmp_path):
        rng = np.random.default_rng(1)
        n_samples, n_variants = 12, 40
        pcs = rng.normal(size=(n_samples, 3))
        kinship = np.corrcoef(rng.normal(size=(n_samples, 60)))
        np.savetxt(tmp_path / "pca.tsv", pcs, delimiter="\t")
        np.savetxt(tmp_path / "kinship.tsv", kinship, delimiter="\t")
        (tmp_path / "calls.vcf").write_text("#CHROM\tPOS\tID\tREF\tALT\n1\t100\t.\tA\tG\n2\t100\t.\tG\tT\n")
        phenotype = rng.normal(size=n_samples).tolist()

        result = suite.generate_all_plots(
            _association_tsv(tmp_path / "assoc.tsv", n=n_variants),
            tmp_path / "plots",
            pca_file=tmp_path / "pca.tsv",
            kinship_file=tmp_path / "kinship.tsv",
            vcf_file=tmp_path / "calls.vcf",
            metadata={
                f"S{i}": {"latitude": 30 + i, "longitude": -100 + i, "population": "AB"[i % 2]}
                for i in range(n_samples)
            },
            genotype_matrix=[[int(g) for g in rng.integers(0, 3, n_samples)] for _ in range(n_variants)],
            kinship_matrix=kinship,
            pca_data={"pcs": pcs},
            phenotypes_dict={"t1": phenotype, "t2": rng.normal(size=n_samples).tolist()},
            phenotype_values=phenotype,
        )

        assert result["status"] == "completed"
        assert result["num_plots_generated"] > 25
        assert plt.get_fignums() == []


# ---------------------------------------------------------------------------
# 8a. effect_size_manhattan: genomic positions
# ---------------------------------------------------------------------------


class TestEffectSizeManhattanPositions:
    @staticmethod
    def _frame() -> pd.DataFrame:
        return pd.DataFrame(
            {
                "CHR": [1, 1, 1, 2, 2],
                "BP": [100, 200, 300, 50, 150],
                "P": [1e-3, 1e-4, 1e-5, 1e-2, 1e-6],
                "BETA": [0.1, -0.2, 0.3, -0.1, 0.2],
            }
        )

    def test_x_is_the_cumulative_genomic_position(self):
        fig = effects.effect_size_manhattan(self._frame())

        points = _offsets(fig.axes[0].collections[0])
        assert points[:, 0].tolist() == [100, 200, 300, 351, 451]  # chromosome 2 starts one base after chromosome 1
        assert points[:, 1] == pytest.approx([3, 4, 5, 2, 6])

    def test_x_follows_the_position_not_the_row_order(self):
        frame = self._frame().iloc[[2, 0, 1, 4, 3]].reset_index(drop=True)

        fig = effects.effect_size_manhattan(frame)

        assert _offsets(fig.axes[0].collections[0])[:, 0].tolist() == [300, 100, 200, 451, 351]

    def test_gaps_between_variants_are_kept(self):
        frame = pd.DataFrame({"CHR": [1, 1], "BP": [1_000, 5_000_000], "P": [1e-3, 1e-4], "BETA": [0.1, 0.2]})

        fig = effects.effect_size_manhattan(frame)

        assert _offsets(fig.axes[0].collections[0])[:, 0].tolist() == [1_000, 5_000_000]

    def test_chromosomes_are_laid_out_in_natural_order(self):
        frame = pd.DataFrame(
            {
                "CHR": ["10", "2", "1"],
                "BP": [100, 100, 100],
                "P": [1e-3, 1e-3, 1e-3],
                "BETA": [0.1, 0.1, 0.1],
            }
        )

        fig = effects.effect_size_manhattan(frame)

        # offsets: chromosome 1 -> 0, chromosome 2 -> 101, chromosome 10 -> 202
        assert _offsets(fig.axes[0].collections[0])[:, 0].tolist() == [302, 201, 100]

    def test_labels_sit_at_the_centre_of_each_chromosome(self):
        fig = effects.effect_size_manhattan(self._frame())

        labels = {t.get_text(): t.get_position()[0] for t in fig.axes[0].texts}
        assert labels == {"1": pytest.approx(200.0), "2": pytest.approx(401.0)}

    def test_colour_values_are_the_absolute_effect_sizes(self):
        fig = effects.effect_size_manhattan(self._frame())

        assert list(fig.axes[0].collections[0].get_array()) == pytest.approx([0.1, 0.2, 0.3, 0.1, 0.2])

    def test_dictionary_input_is_accepted(self):
        data = self._frame().to_dict(orient="list")

        fig = effects.effect_size_manhattan(data)

        assert _offsets(fig.axes[0].collections[0])[:, 0].tolist() == [100, 200, 300, 351, 451]

    def test_missing_columns_return_none(self):
        assert effects.effect_size_manhattan(pd.DataFrame({"CHR": [1], "BP": [1], "P": [0.1]})) is None


# ---------------------------------------------------------------------------
# 8c. plot_gene_track honours the chromosome
# ---------------------------------------------------------------------------


class TestGeneTrackChromosome:
    @staticmethod
    def _track(chrom, genes):
        _, ax = plt.subplots()
        regional.plot_gene_track(ax, chrom, 0, 100, {"genes": genes})
        return ax

    GENES = [
        {"name": "on_1", "chrom": "1", "start": 10, "end": 50},
        {"name": "on_2", "chrom": "2", "start": 10, "end": 50},
        {"name": "unlabelled", "start": 60, "end": 90},
    ]

    def test_genes_on_other_chromosomes_are_not_drawn(self):
        ax = self._track("1", self.GENES)

        assert sorted((p.get_x(), p.get_width()) for p in ax.patches) == [(10, 40), (60, 30)]
        assert sorted(t.get_text() for t in ax.texts) == ["on_1", "unlabelled"]

    def test_chromosome_types_are_compared_as_text(self):
        for chrom in (1, 1.0, "1"):
            ax = self._track(chrom, self.GENES)
            assert sorted(t.get_text() for t in ax.texts) == ["on_1", "unlabelled"]
            plt.close("all")

    def test_every_chromosome_spelling_on_the_gene_is_read(self):
        for key in ("chromosome", "chrom", "chr", "CHROM", "CHR"):
            ax = self._track("1", [{"name": "g", key: "2", "start": 10, "end": 50}])
            assert len(ax.patches) == 0
            plt.close("all")

    def test_only_foreign_genes_leaves_an_honest_note(self):
        ax = self._track("1", [{"name": "on_2", "chrom": "2", "start": 10, "end": 50}])

        assert len(ax.patches) == 0
        assert [t.get_text() for t in ax.texts] == ["No annotated genes overlap this region"]

    def test_genes_without_a_chromosome_are_assumed_to_be_on_the_track(self):
        ax = self._track("7", [{"name": "plain", "start": 10, "end": 50}])

        assert [(p.get_x(), p.get_width()) for p in ax.patches] == [(10, 40)]


# ---------------------------------------------------------------------------
# hwe_deviation_plot: the Q-Q panel's expected values were random draws
# ---------------------------------------------------------------------------


class TestHweQQExpectedValues:
    P_VALUES = [0.9, 0.5, 0.2, 0.05, 1e-3, 1e-6]

    def test_expected_values_are_the_theoretical_null_quantiles(self):
        n = len(self.P_VALUES)

        expected = variants._null_neg_log10_quantiles(n)

        assert list(expected) == sorted(expected)
        assert expected.tolist() == pytest.approx(sorted(-math.log10((i - 0.5) / n) for i in range(1, n + 1)))

    def test_the_figure_is_identical_from_run_to_run(self, tmp_path):
        first, second = tmp_path / "first.png", tmp_path / "second.png"

        variants.hwe_deviation_plot(self.P_VALUES, first)
        variants.hwe_deviation_plot(self.P_VALUES, second)

        assert first.read_bytes() == second.read_bytes()

    def test_metadata_is_unchanged(self, tmp_path):
        result = variants.hwe_deviation_plot(self.P_VALUES, tmp_path / "hwe.png")

        assert result == {"status": "success", "n_variants": 6}


# ---------------------------------------------------------------------------
# recombination_rate_plot: the smoothing window must be longer than the polynomial order
# ---------------------------------------------------------------------------


class TestRecombinationSmoothing:
    """A chromosome of 11 to 19 positions got a window of 3, which a cubic ``savgol_filter`` rejects (ValueError)."""

    @staticmethod
    def _series(n_points: int) -> tuple[list[float], list[int]]:
        rates = np.random.default_rng(n_points).uniform(0.1, 3.0, n_points).tolist()
        return rates, list(range(1000, 1000 + 1000 * n_points, 1000))

    @pytest.mark.parametrize("n_points", [5, 10, 11, 12, 15, 19, 20, 25, 99, 100, 250])
    def test_every_series_length_is_plotted(self, n_points, tmp_path):
        pytest.importorskip("scipy")  # the smoothing line needs SciPy
        rates, positions = self._series(n_points)
        out = tmp_path / "recombination.png"

        result = regional.recombination_rate_plot([], out, recombination_rates={"1": rates}, positions={"1": positions})

        assert result["status"] == "success"
        assert out.stat().st_size > 0

    def test_several_chromosomes_with_different_lengths_are_plotted_together(self, tmp_path):
        pytest.importorskip("scipy")
        lengths = {"1": 8, "2": 15, "3": 40}
        series = {chrom: self._series(n) for chrom, n in lengths.items()}

        result = regional.recombination_rate_plot(
            [],
            tmp_path / "recombination.png",
            recombination_rates={chrom: rates for chrom, (rates, _) in series.items()},
            positions={chrom: positions for chrom, (_, positions) in series.items()},
        )

        assert (result["status"], result["n_chromosomes"]) == ("success", 3)

    @pytest.mark.parametrize("n_points", [11, 12, 15, 19])
    def test_short_series_are_smoothed_with_a_window_of_five(self, n_points):
        scipy_signal = pytest.importorskip("scipy.signal")
        rates, _ = self._series(n_points)

        smoothed = regional._smoothed_rates(rates)

        assert len(smoothed) == n_points
        assert smoothed == pytest.approx(scipy_signal.savgol_filter(rates, 5, 3))

    @pytest.mark.parametrize("n_points", [20, 25, 39, 40, 99, 100, 250])
    def test_longer_series_keep_the_window_they_always_had(self, n_points):
        scipy_signal = pytest.importorskip("scipy.signal")
        rates, _ = self._series(n_points)
        window = min(21, n_points // 10 * 2 + 1)

        smoothed = regional._smoothed_rates(rates)

        assert smoothed == pytest.approx(scipy_signal.savgol_filter(rates, window, 3))

    @pytest.mark.parametrize("n_points", [0, 1, 5, 10])
    def test_series_of_ten_or_fewer_get_no_smoothing_line(self, n_points):
        assert regional._smoothed_rates([1.0] * n_points) is None
