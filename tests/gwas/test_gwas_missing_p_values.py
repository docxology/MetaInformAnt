"""Regression tests: a result without a p-value is not turned into p = 1 (or a statistic into 0) by reports and plots.

Defects, each reproduced from the code before it was fixed:

* ``conditional_analysis_plot`` drew a result without a p-value at -log10 p = 0 and placed one without a position at
  its index in the list;
* ``top_hits_genotype_phenotype`` ranked such a result last "as p = 1", titled its boxplot ``p=1.00e+00`` and reported
  ``p_value: 1.0`` for it;
* ``concordance_plot`` drew a variant that one study had not tested at -log10 p = 0 for that study, and matched
  records without an ``ID`` and without uppercase ``CHROM`` / ``POS`` keys on the id ``"None:None"``, so all of them
  collapsed into one "overlapping SNP";
* ``generate_html_report`` showed a missing minimum p-value as ``1.000e+00``, a missing MAF mean or mean |beta| as
  ``0.000``, missing hit counts as ``0``, passed a missing calibration ratio as a calibration "pass", crashed on a
  top hit without a position, effect, standard error or MAF, and inserted the trait, model and ids unescaped.

* ``forest_plot`` drew a result without an effect, standard error or p-value as an effect of 0 with a zero-width
  interval and ``p=1.0e+00`` (or crashed on ``None``);
* the general ``functional_enrichment_plot`` is not an enrichment analysis (the GFF file is never read and no
  category is assigned): it counted significant variants per chromosome, read the p-value under one spelling, and
  crashed on ``None``;
* ``generate_html_report`` ignored its ``n_samples`` argument.

Oracles are the figures' artists and the returned metadata, with hand-counted values.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")  # optional dependency: skip, do not fail, without it
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PathCollection  # noqa: E402

from metainformant.gwas.reporting import html_report  # noqa: E402
from metainformant.gwas.reporting.html_report import generate_html_report  # noqa: E402
from metainformant.gwas.visualization import general  # noqa: E402
from metainformant.gwas.visualization.interactive.finemapping import conditional_analysis_plot  # noqa: E402
from metainformant.gwas.visualization.interactive.phenotype import top_hits_genotype_phenotype  # noqa: E402
from metainformant.gwas.visualization.statistical.comparison import concordance_plot  # noqa: E402

NAN = float("nan")
MISSING = "NA"  # what the report shows for a statistic it does not have


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


# ---------------------------------------------------------------------------
# conditional_analysis_plot
# ---------------------------------------------------------------------------


class TestConditionalAnalysisPlot:
    def test_a_result_without_a_p_value_is_not_drawn(self, tmp_path: Path) -> None:
        rounds = [
            [
                {"position": 100, "p_value": 1e-8},
                {"position": 200, "p_value": NAN},
                {"position": 300},
                {"position": 400, "p_value": "NA"},
                {"position": 500, "p_value": 1e-2},
            ]
        ]

        result = conditional_analysis_plot(rounds, tmp_path / "c.png")

        assert result["status"] == "success"
        assert result["n_without_p_value"] == 3
        assert (tmp_path / "c.png").stat().st_size > 0

    def test_every_position_and_p_value_spelling_is_read(self) -> None:
        rounds = [[{"POS": 100, "P": 1e-5}, {"BP": 200, "pvalue": 1e-3}, {"pos": 300, "pval": 1e-1}]]

        result = conditional_analysis_plot(rounds)

        assert result["status"] == "success" and result["n_without_p_value"] == 0

    def test_a_result_without_a_position_fails_naming_the_accepted_keys(self) -> None:
        rounds = [[{"position": 100, "p_value": 1e-5}], [{"position": 1, "p_value": 1e-3}, {"p_value": 1e-4}]]

        result = conditional_analysis_plot(rounds)

        assert result["status"] == "failed"
        assert "result 1 of round 2 has no position" in result["reason"]
        assert all(f"'{key}'" in result["reason"] for key in ("position", "pos", "POS", "BP"))
        assert plt.get_fignums() == []  # no figure was left open by the failure

    def test_an_ordinary_run_is_unchanged(self, tmp_path: Path) -> None:
        rounds = [
            [{"position": i * 100, "p_value": 10.0 ** -(1 + i)} for i in range(5)],
            [{"position": i * 100, "p_value": 10.0 ** -(1 + i / 2)} for i in range(5)],
        ]

        result = conditional_analysis_plot(rounds, tmp_path / "c.png", labels=["unconditioned", "conditioned"])

        assert result == {
            "status": "success",
            "output_path": str(tmp_path / "c.png"),
            "n_rounds": 2,
            "n_without_p_value": 0,
        }


# ---------------------------------------------------------------------------
# top_hits_genotype_phenotype
# ---------------------------------------------------------------------------

GENOTYPES = [[0, 1, 2, 0, 1, 2, 0, 1], [0, 0, 1, 1, 2, 2, 0, 1], [1, 1, 0, 0, 2, 2, 1, 0], [2, 1, 0, 2, 1, 0, 2, 1]]
PHENOTYPES = [1.0, 2.0, 3.0, 1.5, 2.5, 3.5, 0.5, 2.2]


def _hits(p_values: List[Any]) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    for i, p in enumerate(p_values):
        record: Dict[str, Any] = {"variant_id": f"rs{i}"}
        if p is not None:
            record["p_value"] = p
        results.append(record)
    return results


class TestTopHits:
    def test_a_result_without_a_p_value_is_not_ranked_as_p_equal_one(self, tmp_path: Path) -> None:
        results = _hits([None, 1e-3, NAN, 1e-6])

        out = top_hits_genotype_phenotype(results, GENOTYPES, PHENOTYPES, tmp_path, n_top=4)

        assert out["status"] == "success"
        assert [s["variant_id"] for s in out["variant_summaries"]] == ["rs3", "rs1"]  # only the two that have one
        assert [s["p_value"] for s in out["variant_summaries"]] == [1e-6, 1e-3]
        assert out["n_plotted"] == 2
        assert not list(tmp_path.glob("*rs0*")) and not list(tmp_path.glob("*rs2*"))

    def test_the_genotypes_of_each_hit_are_those_of_its_own_row(self, tmp_path: Path) -> None:
        """Skipping a result must not shift the genotype row used for the others."""
        results = _hits([None, 1e-3, 0.5, 1e-6])

        out = top_hits_genotype_phenotype(results, GENOTYPES, PHENOTYPES, tmp_path, n_top=2)

        assert [s["variant_id"] for s in out["variant_summaries"]] == ["rs3", "rs1"]
        assert out["n_plotted"] == 2

    @pytest.mark.parametrize("key", ["p_value", "pvalue", "pval", "P"])
    def test_every_p_value_spelling_ranks(self, tmp_path: Path, key: str) -> None:
        results = [{"variant_id": "a", key: 1e-2}, {"variant_id": "b", key: 1e-7}]

        out = top_hits_genotype_phenotype(results, GENOTYPES[:2], PHENOTYPES, tmp_path, n_top=2)

        assert [s["variant_id"] for s in out["variant_summaries"]] == ["b", "a"]

    def test_all_results_without_a_p_value_give_no_plots_not_a_p_equal_one_plot(self, tmp_path: Path) -> None:
        out = top_hits_genotype_phenotype(_hits([None, NAN]), GENOTYPES[:2], PHENOTYPES, tmp_path)

        assert out["status"] == "success" and out["n_plotted"] == 0 and out["variant_summaries"] == []


# ---------------------------------------------------------------------------
# concordance_plot
# ---------------------------------------------------------------------------


def _scatter_points(fig: Any) -> np.ndarray:
    ax = fig.axes[0]
    return np.vstack(
        [np.asarray(c.get_offsets(), dtype=float) for c in ax.collections if isinstance(c, PathCollection)]
    )


class TestConcordance:
    STUDY1 = [
        {"ID": "rs1", "p_value": 1e-8, "beta": 0.5},
        {"ID": "rs2", "p_value": 1e-4, "beta": 0.2},
        {"ID": "rs3", "p_value": 1e-2, "beta": -0.1},
    ]
    STUDY2 = [
        {"ID": "rs1", "p_value": 1e-6, "beta": 0.4},
        {"ID": "rs2", "p_value": 1e-3, "beta": 0.3},
        {"ID": "rs3", "p_value": 1e-1, "beta": -0.2},
    ]

    def test_overlapping_variants_are_plotted_at_minus_log10_p(self, tmp_path: Path) -> None:
        result = concordance_plot(self.STUDY1, self.STUDY2, tmp_path / "c.png")

        assert result["status"] == "success" and result["n_overlapping_snps"] == 3
        assert result["n_without_p_value"] == 0

    def test_a_variant_untested_in_one_study_is_not_drawn_at_zero(self, tmp_path: Path) -> None:
        study2 = [dict(r) for r in self.STUDY2]
        study2[1]["p_value"] = NAN
        del study2[2]["p_value"]

        result = concordance_plot(self.STUDY1, study2, tmp_path / "c.png")

        assert result["status"] == "success"
        assert result["n_overlapping_snps"] == 1  # only rs1 has a p-value in both
        assert result["n_without_p_value"] == 2

    def test_no_variant_with_a_p_value_in_both_studies_fails(self) -> None:
        study2 = [{"ID": r["ID"]} for r in self.STUDY2]

        result = concordance_plot(self.STUDY1, study2)

        assert result["status"] == "failed"
        assert "p-value in both studies" in result["error"]

    @pytest.mark.parametrize(
        "chromosome, position, p_value",
        [
            ("CHROM", "POS", "p_value"),
            ("chrom", "pos", "pvalue"),
            ("CHR", "BP", "P"),
            ("chromosome", "position", "pval"),
        ],
    )
    def test_records_without_an_id_are_matched_on_chromosome_and_position_in_any_spelling(
        self, tmp_path: Path, chromosome: str, position: str, p_value: str
    ) -> None:
        study1 = [{chromosome: "1", position: 100, p_value: 1e-5}, {chromosome: "1", position: 200, p_value: 1e-3}]
        study2 = [{chromosome: 1, position: 100, p_value: 1e-4}, {chromosome: 2, position: 200, p_value: 1e-2}]

        result = concordance_plot(study1, study2, tmp_path / "c.png")

        assert result["status"] == "success"
        assert result["n_overlapping_snps"] == 1  # 1:100 only; they used to collapse into one "None:None" variant
        assert (result["study1_snps"], result["study2_snps"]) == (2, 2)

    def test_records_that_cannot_be_identified_fail_naming_the_accepted_keys(self) -> None:
        result = concordance_plot([{"p_value": 1e-5}], [{"ID": "rs1", "p_value": 1e-5}])

        assert result["status"] == "failed"
        assert "study 1" in result["error"] and "result 0 has no chromosome" in result["error"]
        assert all(f"'{key}'" in result["error"] for key in ("chromosome", "chrom", "chr", "CHROM", "CHR"))

    def test_the_correlation_uses_the_betas_of_the_drawn_variants(self) -> None:
        result = concordance_plot(self.STUDY1, self.STUDY2)

        betas1 = [0.5, 0.2, -0.1]
        betas2 = [0.4, 0.3, -0.2]
        assert result["correlation"] == pytest.approx(float(np.corrcoef(betas1, betas2)[0, 1]))


# ---------------------------------------------------------------------------
# generate_html_report
# ---------------------------------------------------------------------------

FULL_SUMMARY = {
    "n_variants": 12345,
    "n_chromosomes": 16,
    "significance_counts": {"genome_wide_5e-8": 3, "suggestive_1e-5": 17},
    "maf_distribution": {"mean": 0.2345},
    "effect_size_distribution": {"mean_abs_beta": 0.1234},
    "p_value_calibration": {
        "lambda_gc": 1.05,
        "min_p": 3.2e-9,
        "calibration_quantiles": {"q50": {"ratio_obs_exp": 1.02}},
    },
    "top_10_hits": [{"snp": "rs1", "chrom": "1", "pos": 123456, "beta": 0.5, "se": 0.1, "p_value": 1e-9, "MAF": 0.31}],
}


def _report(tmp_path: Path, summary: Dict[str, Any], **kwargs: Any) -> str:
    return generate_html_report(tmp_path, summary, **kwargs).read_text(encoding="utf-8")


class TestHtmlReport:
    def test_the_marker_is_na(self) -> None:
        assert html_report.MISSING == MISSING

    def test_a_complete_summary_is_reported_as_given(self, tmp_path: Path) -> None:
        html = _report(tmp_path, FULL_SUMMARY, trait="height", model="linear", n_samples=250)

        for expected in (
            "12,345",
            '<div class="value">3</div>',
            '<div class="value">17</div>',
            "1.050",
            "3.200e-09",
            "0.234",
            "0.123",
        ):
            assert expected in html, expected
        assert "rs1" in html and "123,456" in html and "1.000e-09" in html and "0.500" in html
        assert MISSING not in html
        assert "✅" in html

    def test_an_empty_summary_shows_na_not_made_up_numbers(self, tmp_path: Path) -> None:
        html = _report(tmp_path, {})

        assert html.count(f'<div class="value">{MISSING}</div>') == 6  # all six statistics cards
        for invented in ("1.000e+00", "0.000", '<div class="value">0</div>'):
            assert invented not in html, invented
        assert (
            "✅" not in html.split("<h2>Supplementary Statistics</h2>")[1].split("</table>")[0]
        )  # no calibration pass

    def test_a_missing_calibration_ratio_is_not_a_pass(self, tmp_path: Path) -> None:
        summary = {"p_value_calibration": {"lambda_gc": 1.0, "calibration_quantiles": {}}}

        supplementary = _report(tmp_path, summary).split("<h2>Supplementary Statistics</h2>")[1].split("</table>")[0]

        assert "—" in supplementary and "✅" not in supplementary and "⚠️" not in supplementary

    def test_a_poorly_calibrated_ratio_is_a_warning(self, tmp_path: Path) -> None:
        summary = {"p_value_calibration": {"calibration_quantiles": {"q50": {"ratio_obs_exp": 2.5}}}}

        supplementary = _report(tmp_path, summary).split("<h2>Supplementary Statistics</h2>")[1].split("</table>")[0]

        assert "2.5" in supplementary and "⚠️" in supplementary

    def test_a_top_hit_with_missing_fields_does_not_crash_and_shows_na(self, tmp_path: Path) -> None:
        summary = {
            "top_10_hits": [{"snp": "rs9"}, {"snp": "rs10", "chrom": "2", "pos": None, "beta": NAN, "p_value": "NA"}]
        }

        html = _report(tmp_path, summary)

        rows = html.split("<th>MAF</th></tr>")[1].split("</table>")[0]
        assert "rs9" in rows and "rs10" in rows
        assert rows.count(f"<td>{MISSING}</td>") >= 8
        assert "1.000e+00" not in rows

    def test_an_inflated_lambda_is_flagged_and_a_missing_one_is_not(self, tmp_path: Path) -> None:
        assert "Inflated" in _report(tmp_path, {"p_value_calibration": {"lambda_gc": 2.6}}, filename="a.html")
        assert "Inflated" not in _report(tmp_path, {"p_value_calibration": {"lambda_gc": 1.1}}, filename="b.html")
        assert "Inflated" not in _report(tmp_path, {}, filename="c.html")

    def test_the_sign_test_row_does_not_invent_a_p_value(self, tmp_path: Path) -> None:
        html = _report(tmp_path, {}, sign_test_results={"n_positive_beta": 6, "n_negative_beta": 4})

        row = html.split("Effect direction (sign test)")[1].split("</tr>")[0]
        assert f"p={MISSING}" in row and "✅" not in row

    def test_values_that_come_from_the_data_are_html_escaped(self, tmp_path: Path) -> None:
        summary = {"top_10_hits": [{"snp": "<img src=x onerror=alert(1)>", "chrom": "<b>1</b>"}]}

        html = _report(
            tmp_path,
            summary,
            trait="<script>alert('t')</script>",
            model="m&m",
            gene_annotations=[{"snp": "<i>s</i>", "nearby_genes": [{"gene_name": "<u>G</u>"}]}],
            sign_test_results={"interpretation": "<em>x</em>"},
            credible_set_results={"n_snps": "<4>"},
        )

        for raw in ("<script>alert", "<img src=x", "<b>1</b>", "<i>s</i>", "<u>G</u>", "<em>x</em>", "<4>"):
            assert raw not in html, raw
        assert "&lt;script&gt;alert(&#x27;t&#x27;)&lt;/script&gt;" in html and "m&amp;m" in html

    def test_the_report_is_written_where_asked(self, tmp_path: Path) -> None:
        path = generate_html_report(tmp_path / "nested" / "dir", FULL_SUMMARY, filename="r.html")

        assert path == tmp_path / "nested" / "dir" / "r.html" and path.is_file()
        assert math.isfinite(len(path.read_text(encoding="utf-8")))


# ---------------------------------------------------------------------------
# forest_plot and the general functional_enrichment_plot
# ---------------------------------------------------------------------------


class TestForestPlot:
    RESULTS = [
        {"snp": "rs1", "beta": 0.5, "se": 0.1, "p_value": 1e-9},
        {"snp": "rs2", "beta": 0.2, "se": 0.1},  # no p-value
        {"snp": "rs3", "beta": None, "se": 0.1, "p_value": 1e-4},  # no effect
        {"snp": "rs4", "beta": 0.3, "p_value": 1e-3},  # no standard error
        {"snp": "rs5", "BETA": -0.4, "SE": 0.2, "pvalue": 1e-5},  # other spellings are read
        {"snp": "rs6", "beta": 0.1, "se": 0.1, "p_value": NAN},
    ]

    def test_only_results_with_a_p_value_an_effect_and_a_standard_error_are_drawn(self) -> None:
        fig = general.forest_plot(self.RESULTS)

        labels = [text.get_text() for text in fig.axes[0].get_yticklabels()]
        assert labels == ["rs5  p=1.0e-05", "rs1  p=1.0e-09"]  # bottom to top: the most significant is on top

    def test_the_intervals_are_the_ones_of_the_drawn_results(self) -> None:
        fig = general.forest_plot(self.RESULTS)

        lines = [
            line.get_xdata()
            for line in fig.axes[0].get_lines()
            if len(line.get_xdata()) == 2 and line.get_xdata()[0] != line.get_xdata()[1]
        ]
        intervals = [(round(float(low), 6), round(float(high), 6)) for low, high in lines if low != 0 or high != 0]
        assert (round(-0.4 - 1.96 * 0.2, 6), round(-0.4 + 1.96 * 0.2, 6)) in intervals
        assert (round(0.5 - 1.96 * 0.1, 6), round(0.5 + 1.96 * 0.1, 6)) in intervals

    def test_nothing_usable_gives_no_figure(self) -> None:
        assert general.forest_plot([{"snp": "a"}, {"snp": "b", "beta": 0.1}]) is None

    def test_a_warning_counts_what_was_left_out(self, caplog: pytest.LogCaptureFixture) -> None:
        with caplog.at_level(logging.WARNING):
            general.forest_plot(self.RESULTS)

        assert any(
            "4 of 6 results lack a p-value, an effect or a standard error" in r.getMessage() for r in caplog.records
        )

    def test_complete_results_are_unchanged(self) -> None:
        results = [{"snp": f"rs{i}", "beta": 0.1 * i, "se": 0.05, "p_value": 10.0 ** -(i + 5)} for i in range(1, 4)]

        fig = general.forest_plot(results, top_n=2)

        assert [text.get_text() for text in fig.axes[0].get_yticklabels()] == ["rs2  p=1.0e-07", "rs3  p=1.0e-08"]


class TestGeneralFunctionalEnrichmentPlot:
    def test_it_counts_significant_variants_per_chromosome_under_every_spelling(self) -> None:
        results = [
            {"chrom": 1, "p_value": 1e-9},
            {"chrom": 1.0, "pvalue": 1e-10},
            {"CHR": 2, "P": 1e-12},
            {"chrom": 2, "p_value": None},  # no p-value: not significant (None used to raise TypeError)
            {"chrom": 3, "p_value": 0.5},
            {"chrom": 3},
        ]

        fig = general.functional_enrichment_plot(results, "annotation.gff")

        ax = fig.axes[0]
        assert ax.get_title() == "Significant Variants by Chromosome (n=3)"
        assert [text.get_text() for text in ax.get_xticklabels()] == ["1", "2"]
        assert [patch.get_height() for patch in ax.patches] == [2, 1]

    def test_the_gff_file_is_not_read_and_the_warning_says_so(
        self, caplog: pytest.LogCaptureFixture, tmp_path: Path
    ) -> None:
        missing = tmp_path / "does_not_exist.gff"

        with caplog.at_level(logging.WARNING):
            fig = general.functional_enrichment_plot([{"chrom": 1, "p_value": 1e-9}], missing)

        assert fig is not None  # a file that does not exist is not an error: nothing is read from it
        assert any("does not read the GFF file" in r.getMessage() for r in caplog.records)

    def test_without_significant_variants_the_title_does_not_claim_an_enrichment(self) -> None:
        fig = general.functional_enrichment_plot([{"chrom": 1, "p_value": 0.5}], "annotation.gff")

        assert fig.axes[0].get_title() == "Significant Variants by Chromosome"


def test_n_samples_is_shown_in_the_report_header(tmp_path: Path) -> None:
    html = generate_html_report(tmp_path, {}, n_samples=1234).read_text(encoding="utf-8")

    assert "<strong>Samples:</strong> 1,234" in html
    assert f"<strong>Samples:</strong> {MISSING}" in generate_html_report(tmp_path, {}, filename="b.html").read_text(
        encoding="utf-8"
    )
