"""Regression tests: a failed association test is missing data, not an effect of 0 with p = 1.

Defect: ``association_test_linear`` (and the logistic and mixed-model tests) returned ``beta = 0``, ``se = 0`` and
``p_value = 1`` for a test that could not be fitted. Everything downstream then treated the failure as a real null
result: the multiple-testing corrections counted it in the family (inflating Bonferroni's ``m`` and every
Benjamini-Hochberg rank), lambda_GC took its chi-square of 0 into the median, the summary tables wrote
``0 / 0 / 1.000000e+00`` as if they were estimates, and ``_apply_full_correction_outputs`` even filtered rows without a
p-value out of the family and then wrote the corrected values back to the wrong rows.

A failed test now carries NaN for its statistics plus ``status``/``error``; the corrections leave it out of the
family, the summary writers print ``NA`` and the summaries, plots, top hit and fine-mapping skip it.

Oracles: ``scipy.stats.linregress`` for the fitted tests, ``scipy.stats.false_discovery_control`` (Benjamini-Hochberg)
and hand arithmetic for the corrections, and ``scipy.stats.chi2`` for lambda_GC.
"""

from __future__ import annotations

import json
import logging
import math
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pytest
from scipy import stats

from metainformant.gwas.analysis.association import (
    association_test_linear,
    association_test_logistic,
    run_linear_model_gwas,
    run_logistic_model_gwas,
)
from metainformant.gwas.analysis.correction import (
    adjust_p_values,
    bonferroni_correction,
    fdr_correction,
    genomic_control,
    lambda_gc_from_p_values,
    qvalue_estimation,
)
from metainformant.gwas.analysis.mixed_model import association_test_mixed, run_mixed_model_gwas
from metainformant.gwas.analysis.power import _compute_statistic, subsample_convergence
from metainformant.gwas.analysis.summary_stats import (
    compute_comprehensive_summary,
    create_results_summary,
    per_chromosome_summary,
    sign_test,
    write_significant_hits,
    write_summary_statistics,
)
from metainformant.gwas.workflow.workflow_execution import (
    _apply_full_correction_outputs,
    _apply_qvalue_and_bonferroni_annotations,
    _lambda_gc_from_associations,
    execute_gwas_workflow,
    run_gwas,
)

NAN = float("nan")
N_SAMPLES = 12

# Genotype dosages per variant over twelve samples. rsM is monomorphic and rsN heterozygous in every sample: neither
# can be tested (the genotype is constant), so these two tests fail, and they sit between and after valid variants.
GENOTYPES: Dict[str, List[int]] = {
    "rsA": [0, 1, 2, 0, 1, 2, 0, 1, 2, 1, 0, 2],
    "rsM": [0] * N_SAMPLES,
    "rsB": [2, 2, 1, 0, 0, 1, 2, 1, 0, 0, 1, 1],
    "rsC": [0, 0, 1, 1, 2, 2, 1, 0, 2, 1, 2, 0],
    "rsN": [1] * N_SAMPLES,
    "rsD": [1, 0, 2, 1, 0, 1, 2, 2, 0, 1, 1, 0],
    "rsE": [0, 2, 1, 1, 0, 2, 2, 0, 1, 0, 2, 1],
}
FAILING = ["rsM", "rsN"]
VALID = [v for v in GENOTYPES if v not in FAILING]
# rsA carries a real effect (about +0.9 per allele), rsC a weaker one; the rest of the phenotype is fixed noise.
NOISE = [0.12, -0.31, 0.25, 0.4, -0.22, 0.05, -0.37, 0.19, -0.08, 0.33, -0.15, 0.28]
PHENOTYPE = [
    1.5 + 0.9 * a + 0.3 * c + e
    for a, c, e in zip(GENOTYPES["rsA"], GENOTYPES["rsC"], NOISE)  # noqa: B905 - lengths are equal
]
VARIANT_INFO = [{"id": v, "chrom": "chr1" if i < 4 else "chr2", "pos": 100 * (i + 1)} for i, v in enumerate(GENOTYPES)]
GENOTYPE_ROWS = list(GENOTYPES.values())


def is_nan(value: Any) -> bool:
    return isinstance(value, float) and math.isnan(value)


def linear_rows() -> List[Dict[str, Any]]:
    return run_linear_model_gwas(GENOTYPE_ROWS, PHENOTYPE, variant_info=VARIANT_INFO)


# The p-values of the tests that can be fitted, from scipy, in the order of VALID (independent of the code under test)
VALID_P = [float(stats.linregress(GENOTYPES[v], PHENOTYPE).pvalue) for v in VALID]


def bh(p_values: List[float]) -> List[float]:
    return [float(q) for q in stats.false_discovery_control(p_values, method="bh")]


def lambda_gc_oracle(p_values: List[float]) -> float:
    chi2 = np.array([stats.chi2.isf(p, 1) for p in p_values])
    return float(np.median(chi2) / stats.chi2.ppf(0.5, 1))


def read_table(path: Path) -> List[Dict[str, str]]:
    lines = path.read_text().strip().split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


# ---------------------------------------------------------------------------
# The tests themselves
# ---------------------------------------------------------------------------


def assert_missing_statistics(result: Dict[str, Any], fields: List[str]) -> None:
    for field in fields:
        assert is_nan(result[field]), f"{field} should be NaN, got {result[field]!r}"
    assert result["converged"] is False
    assert isinstance(result["error"], str) and result["error"]


LINEAR_FIELDS = ["beta", "se", "t_stat", "p_value", "r_squared"]
LOGISTIC_FIELDS = ["beta", "se", "z_stat", "p_value", "odds_ratio"]


@pytest.mark.parametrize(
    ("genotypes", "phenotypes", "covariates", "message"),
    [
        ([1] * 8, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0], None, "constant"),
        ([-1, -1, -1, -1], [1.0, 2.0, 3.0, 4.0], None, "No complete cases"),
        ([0, 1], [1.0, 2.0], None, "at least 3 complete cases"),
        (
            [0, 1, 2, 0, 1],
            [1.0, 2.0, 3.0, 1.5, 2.5],
            [[1.0, 2.0, 3.0, 4.0, 5.0], [2.0, 1.0, 4.0, 3.0, 6.0], [0.0, 1.0, 0.0, 1.0, 1.0]],
            "degrees of freedom",
        ),
        (
            [0, 1, 2, 0, 1, 2, 1, 0],
            [1.0, 2.0, 3.0, 1.5, 2.5, 3.5, 2.0, 1.0],
            [[3.0 * g for g in [0, 1, 2, 0, 1, 2, 1, 0]]],  # a covariate collinear with the genotype
            "singular",
        ),
    ],
)
def test_a_linear_test_that_cannot_be_fitted_has_missing_statistics(
    genotypes: List[int], phenotypes: List[float], covariates: Any, message: str
) -> None:
    result = association_test_linear(genotypes, phenotypes, covariates=covariates)

    assert result["status"] == "error"
    assert message in result["error"]
    assert_missing_statistics(result, LINEAR_FIELDS)


def test_a_linear_test_that_is_fitted_is_unchanged_and_has_no_error_field() -> None:
    result = association_test_linear(GENOTYPES["rsA"], PHENOTYPE)

    oracle = stats.linregress(GENOTYPES["rsA"], PHENOTYPE)
    assert result["status"] == "success" and result["converged"] is True and "error" not in result
    assert result["beta"] == pytest.approx(oracle.slope, rel=1e-9)
    assert result["se"] == pytest.approx(oracle.stderr, rel=1e-9)
    assert result["p_value"] == pytest.approx(oracle.pvalue, rel=1e-6)


def test_a_logistic_test_without_enough_cases_has_missing_statistics() -> None:
    result = association_test_logistic([0, 1, 2, 0, 1, 2], [0, 0, 0, 0, 0, 1])  # one case

    assert result["status"] == "failed"
    assert result["n_cases"] == 1 and result["n_controls"] == 5
    assert_missing_statistics(result, LOGISTIC_FIELDS)


def test_a_logistic_fit_that_fails_has_missing_statistics() -> None:
    result = association_test_logistic([0] * 8, [0, 1, 0, 1, 0, 1, 0, 1])  # constant genotype: singular Hessian

    assert result["status"] == "error"
    assert_missing_statistics(result, LOGISTIC_FIELDS)


def test_a_mixed_model_test_of_a_monomorphic_snp_has_missing_statistics() -> None:
    kinship = np.eye(N_SAMPLES).tolist()

    result = association_test_mixed(GENOTYPES["rsM"], PHENOTYPE, kinship)

    assert result["status"] == "error"
    assert "singular" in result["error"]
    assert_missing_statistics(result, ["beta", "se", "t_stat", "p_value"])
    # Only the SNP failed: the null model's variance components are real numbers
    assert all(math.isfinite(result[key]) for key in ("sigma_g", "sigma_e", "heritability"))


def test_a_mixed_model_test_that_is_fitted_reports_success() -> None:
    result = association_test_mixed(GENOTYPES["rsA"], PHENOTYPE, np.eye(N_SAMPLES).tolist())

    assert result["status"] == "success" and result["converged"] is True and "error" not in result
    assert 0 < result["p_value"] <= 1


def test_the_mixed_model_batch_marks_failed_rows_and_keeps_the_others() -> None:
    kinship = np.eye(N_SAMPLES).tolist()

    rows = run_mixed_model_gwas(GENOTYPE_ROWS, PHENOTYPE, kinship, variant_info=VARIANT_INFO)

    assert [row["variant_id"] for row in rows] == list(GENOTYPES)  # one row per variant, in order
    for row in rows:
        if row["variant_id"] in FAILING:
            assert row["status"] == "error" and row["converged"] is False and "singular" in row["error"]
            assert all(is_nan(row[key]) for key in ("beta", "se", "t_stat", "p_value"))
        else:
            single = association_test_mixed(GENOTYPES[row["variant_id"]], PHENOTYPE, kinship)
            assert row["status"] == "success" and "error" not in row
            assert row["p_value"] == pytest.approx(single["p_value"], rel=1e-9)


def test_the_linear_batch_keeps_one_row_per_variant_with_failed_rows_missing() -> None:
    rows = linear_rows()

    assert [row["variant_id"] for row in rows] == list(GENOTYPES)
    for row in rows:
        if row["variant_id"] in FAILING:
            assert row["status"] == "error"
            assert all(is_nan(row[key]) for key in LINEAR_FIELDS)
        else:
            oracle = stats.linregress(GENOTYPES[row["variant_id"]], PHENOTYPE)
            assert row["status"] == "success"
            assert row["beta"] == pytest.approx(oracle.slope, rel=1e-9)
            assert row["p_value"] == pytest.approx(oracle.pvalue, rel=1e-6)


def test_the_logistic_batch_keeps_one_row_per_variant_with_failed_rows_missing() -> None:
    cases = [0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0]

    rows = run_logistic_model_gwas(GENOTYPE_ROWS, cases, variant_info=VARIANT_INFO)

    assert [row["variant_id"] for row in rows] == list(GENOTYPES)
    assert [row["status"] for row in rows if row["variant_id"] in FAILING] == ["error", "error"]
    assert all(is_nan(row["p_value"]) for row in rows if row["variant_id"] in FAILING)
    assert all(0 < row["p_value"] <= 1 for row in rows if row["variant_id"] not in FAILING)


# ---------------------------------------------------------------------------
# Corrections leave missing p-values out of the family
# ---------------------------------------------------------------------------

# Three real p-values and two failed tests, in the order a pipeline would hold them.
P_WITH_FAILURES = [0.001, NAN, 0.01, NAN, 0.04]
P_PRESENT = [0.001, 0.01, 0.04]


@pytest.mark.parametrize("missing", [NAN, None])
def test_fdr_ignores_missing_p_values_instead_of_counting_them_as_one(missing: Any) -> None:
    p_values = [0.001, missing, 0.01, missing, 0.04]

    result = fdr_correction(p_values)

    adjusted = result["adjusted_p_values"]
    assert [adjusted[0], adjusted[2], adjusted[4]] == pytest.approx(bh(P_PRESENT))  # 0.003, 0.015, 0.04
    assert is_nan(adjusted[1]) and is_nan(adjusted[3])
    assert result["significant"] == [True, False, True, False, True]
    assert (result["n_tests"], result["n_excluded"]) == (3, 2)
    # What the old failed-as-p=1 family gave: five tests, so every adjusted value was inflated
    old_family = bh(P_PRESENT + [1.0, 1.0])[:3]
    assert adjusted[0] != pytest.approx(old_family[0])


def test_fdr_by_method_ignores_missing_p_values() -> None:
    result = fdr_correction(P_WITH_FAILURES, method="by")

    oracle = [float(q) for q in stats.false_discovery_control(P_PRESENT, method="by")]
    adjusted = result["adjusted_p_values"]
    assert [adjusted[0], adjusted[2], adjusted[4]] == pytest.approx(oracle)


@pytest.mark.parametrize("missing", [NAN, None])
def test_bonferroni_divides_by_the_tests_that_have_a_p_value(missing: Any) -> None:
    result = bonferroni_correction([0.001, missing, 0.01, missing, 0.04], alpha=0.05)

    assert result["n_tests"] == 3 and result["n_excluded"] == 2
    assert result["corrected_alpha"] == pytest.approx(0.05 / 3)  # not 0.05 / 5
    assert result["significant"] == [True, False, True, False, False]  # 0.01 <= 0.01666, 0.04 > 0.01666
    flags, corrected_alpha = bonferroni_correction([0.001, missing, 0.01], alpha=0.05, return_dict=False)
    assert flags == [True, False, True] and corrected_alpha == pytest.approx(0.025)


def test_qvalues_ignore_missing_p_values() -> None:
    q_values, pi0 = qvalue_estimation(P_WITH_FAILURES, pi0=1.0)

    oracle = bh(P_PRESENT)  # with pi0 = 1 the q-value procedure is Benjamini-Hochberg
    assert [q_values[0], q_values[2], q_values[4]] == pytest.approx(oracle)
    assert is_nan(q_values[1]) and is_nan(q_values[3]) and pi0 == 1.0


def test_adjust_p_values_ignores_missing_p_values() -> None:
    bonferroni = adjust_p_values(P_WITH_FAILURES, "bonferroni")
    fdr = adjust_p_values(P_WITH_FAILURES, "fdr")

    assert [bonferroni[0], bonferroni[2], bonferroni[4]] == pytest.approx([0.003, 0.03, 0.12])  # m = 3
    assert [fdr[0], fdr[2], fdr[4]] == pytest.approx(bh(P_PRESENT))
    assert all(is_nan(x) for x in (bonferroni[1], bonferroni[3], fdr[1], fdr[3]))


def test_genomic_control_takes_lambda_from_the_p_values_that_are_present() -> None:
    p_values = [0.9, NAN, 0.5, 0.2, NAN, 0.05, 0.7, 0.35, 0.6]

    result = genomic_control(p_values=p_values)

    present = [p for p in p_values if not math.isnan(p)]
    assert result["lambda_gc"] == pytest.approx(lambda_gc_oracle(present), rel=1e-9)
    assert lambda_gc_from_p_values(p_values) == pytest.approx(lambda_gc_oracle(present), rel=1e-9)
    assert is_nan(result["corrected_p_values"][1]) and is_nan(result["corrected_p_values"][4])  # passed through


def test_every_correction_handles_a_family_without_any_p_value() -> None:
    assert fdr_correction([NAN, None])["n_tests"] == 0
    assert fdr_correction([NAN, None])["significant"] == [False, False]
    assert bonferroni_correction([NAN, None])["n_significant"] == 0
    assert qvalue_estimation([NAN, NAN])[0][0] != qvalue_estimation([NAN, NAN])[0][0]  # NaN


def test_corrections_without_missing_values_are_unchanged() -> None:
    p_values = [0.5, 0.001, 0.03, 0.2, 0.0004]

    assert fdr_correction(p_values)["adjusted_p_values"] == pytest.approx(bh(p_values))
    assert fdr_correction(p_values)["n_excluded"] == 0
    assert bonferroni_correction(p_values)["corrected_alpha"] == pytest.approx(0.01)


# ---------------------------------------------------------------------------
# Workflow annotations
# ---------------------------------------------------------------------------


def test_qvalue_and_bonferroni_annotations_leave_failed_rows_out_of_the_family() -> None:
    rows = linear_rows()
    present = VALID_P

    _apply_qvalue_and_bonferroni_annotations(rows)

    expected_q = iter(bh(present))
    for row in rows:
        if row["variant_id"] in FAILING:
            assert is_nan(row["q_value"]) and row["bonferroni_significant"] is False
        else:
            assert row["q_value"] == pytest.approx(next(expected_q))
            assert row["bonferroni_significant"] == (row["p_value"] <= 0.05 / len(present))


def test_full_correction_outputs_annotate_the_right_rows_and_exclude_failed_tests() -> None:
    rows = linear_rows()
    present = VALID_P

    outputs = _apply_full_correction_outputs(rows)

    assert outputs["fdr"]["n_tests"] == outputs["bonferroni"]["n_tests"] == len(present)
    assert outputs["fdr"]["n_excluded"] == len(FAILING)
    expected_q = iter(bh(present))
    for row in rows:
        if row["variant_id"] in FAILING:
            assert is_nan(row["fdr_p_value"]) and row["fdr_significant"] is False
            assert row["bonferroni_significant"] is False
        else:
            assert row["fdr_p_value"] == pytest.approx(next(expected_q))
    assert outputs["genomic_control"]["lambda_gc"] == pytest.approx(lambda_gc_oracle(present), rel=1e-9)


def test_full_correction_outputs_do_not_shift_results_onto_the_wrong_rows() -> None:
    """Rows without a p-value used to be filtered out of the list that the results were then indexed against."""
    rows = [{"p_value": 0.01}, {"p_value": None}, {"p_value": 0.04}, {"p_value": 0.03}]

    _apply_full_correction_outputs(rows, correction_method="fdr")

    oracle = bh([0.01, 0.04, 0.03])
    assert rows[0]["fdr_p_value"] == pytest.approx(oracle[0])
    assert is_nan(rows[1]["fdr_p_value"]) and rows[1]["fdr_significant"] is False
    assert rows[2]["fdr_p_value"] == pytest.approx(oracle[1])  # its own value, not row 1's
    assert rows[3]["fdr_p_value"] == pytest.approx(oracle[2])


def test_full_correction_outputs_report_nothing_when_every_test_failed() -> None:
    assert _apply_full_correction_outputs([{"p_value": NAN}, {"p_value": NAN}]) == {}


def test_lambda_gc_of_the_workflow_ignores_failed_tests() -> None:
    rows = linear_rows()

    assert _lambda_gc_from_associations(rows) == pytest.approx(lambda_gc_oracle(VALID_P), abs=1e-4)


# ---------------------------------------------------------------------------
# Summary tables and summaries
# ---------------------------------------------------------------------------


def test_summary_statistics_write_failed_tests_as_na_and_leave_them_out_of_q_values(tmp_path: Path) -> None:
    rows = linear_rows()

    path = write_summary_statistics(rows, VARIANT_INFO, tmp_path / "summary.tsv", genotypes=GENOTYPE_ROWS)

    table = {row["SNP"]: row for row in read_table(path)}
    assert list(table) == list(GENOTYPES)
    for snp in FAILING:
        assert [table[snp][col] for col in ("BETA", "SE", "P", "Q_FDR")] == ["NA"] * 4
        assert table[snp]["MAF"] != "NA"  # the allele frequency is a property of the genotypes, not of the test
    oracle_q = dict(zip(VALID, bh(VALID_P)))
    for snp in VALID:
        oracle = stats.linregress(GENOTYPES[snp], PHENOTYPE)
        assert float(table[snp]["BETA"]) == pytest.approx(oracle.slope, rel=1e-5)
        assert float(table[snp]["P"]) == pytest.approx(oracle.pvalue, rel=1e-5)
        assert float(table[snp]["Q_FDR"]) == pytest.approx(oracle_q[snp], rel=1e-5)  # m counts the valid tests only


def test_significant_hits_never_contain_a_failed_test(tmp_path: Path) -> None:
    rows = linear_rows()

    path = write_significant_hits(rows, VARIANT_INFO, tmp_path / "hits.tsv", threshold=1.0, genotypes=GENOTYPE_ROWS)

    hits = read_table(path)
    assert sorted(row["SNP"] for row in hits) == sorted(VALID)
    assert [float(row["P"]) for row in hits] == sorted(float(row["P"]) for row in hits)


def test_the_results_summary_counts_only_tests_that_ran(tmp_path: Path) -> None:
    rows = linear_rows()
    present = VALID_P

    summary = create_results_summary(rows, tmp_path / "summary.json", threshold=0.05)

    assert summary["n_variants_tested"] == len(present) and summary["n_tests_failed"] == len(FAILING)
    assert summary["n_significant_genome_wide"] == sum(1 for p in present if p < 0.05)
    assert summary["lambda_gc"] == pytest.approx(lambda_gc_oracle(present), rel=1e-6)
    top = summary["top_hits"]
    assert len(top) == len(present)
    assert [hit["p_value"] for hit in top] == pytest.approx(sorted(present), rel=1e-9)
    assert all(rows[hit["variant_index"]]["variant_id"] in VALID for hit in top)  # index into the full list
    text = (tmp_path / "summary.json").read_text()
    assert "NaN" not in text
    assert json.loads(text)["n_tests_failed"] == len(FAILING)


def tabular_rows() -> List[Dict[str, Any]]:
    """Results in the shape ``compute_comprehensive_summary`` and friends take (MAF, chrom, pos and snp keys)."""
    rows = linear_rows()
    for row, genotypes in zip(rows, GENOTYPE_ROWS):
        alt_frequency = sum(genotypes) / (2 * len(genotypes))
        row["MAF"] = min(alt_frequency, 1 - alt_frequency)
        row["snp"] = row["variant_id"]
    return rows


def test_the_comprehensive_summary_ignores_failed_tests() -> None:
    rows = tabular_rows()
    valid_rows = [row for row in rows if row["variant_id"] in VALID]

    summary = compute_comprehensive_summary(rows, significance_threshold=0.05)

    reference = compute_comprehensive_summary(valid_rows, significance_threshold=0.05)
    assert summary.pop("n_tests_failed") == len(FAILING)
    assert reference.pop("n_tests_failed") == 0
    assert summary == reference
    assert summary["n_variants"] == len(VALID)
    assert "nan" not in json.dumps(summary).lower()


def test_the_per_chromosome_summary_and_sign_test_ignore_failed_tests() -> None:
    rows = tabular_rows()
    valid_rows = [row for row in rows if row["variant_id"] in VALID]

    assert per_chromosome_summary(rows, significance_threshold=0.05) == per_chromosome_summary(
        valid_rows, significance_threshold=0.05
    )
    assert sign_test(rows, significance_threshold=0.05) == sign_test(valid_rows, significance_threshold=0.05)
    for chromosome_row in per_chromosome_summary(rows):
        assert not math.isnan(chromosome_row["MIN_P"]) and not math.isnan(chromosome_row["MEAN_ABS_BETA"])


# ---------------------------------------------------------------------------
# Power analysis
# ---------------------------------------------------------------------------


def test_the_subsample_convergence_metrics_leave_failed_variants_out() -> None:
    with_failures = subsample_convergence(GENOTYPE_ROWS, PHENOTYPE, fractions=[1.0], n_replicates=1, alpha=0.05)
    clean = subsample_convergence([GENOTYPES[v] for v in VALID], PHENOTYPE, fractions=[1.0], n_replicates=1, alpha=0.05)

    for metric in ("mean_abs_beta", "mean_neg_log10_p", "n_significant", "lambda_gc"):
        got = with_failures["metrics"][metric][0]["mean"]
        assert not math.isnan(got)
        assert got == pytest.approx(clean["metrics"][metric][0]["mean"], rel=1e-9), metric


def test_the_jackknife_statistic_leaves_failed_variants_out() -> None:
    indices = list(range(len(GENOTYPE_ROWS)))
    valid_indices = [i for i, v in enumerate(GENOTYPES) if v in VALID]

    for statistic in ("mean_abs_beta", "n_significant", "lambda_gc"):
        with_failures = _compute_statistic(GENOTYPE_ROWS, PHENOTYPE, indices, statistic, 0.05)
        clean = _compute_statistic(GENOTYPE_ROWS, PHENOTYPE, valid_indices, statistic, 0.05)
        assert with_failures == pytest.approx(clean, rel=1e-9), statistic


# ---------------------------------------------------------------------------
# The pipelines
# ---------------------------------------------------------------------------


def vcf_text() -> str:
    calls = {0: "0/0", 1: "0/1", 2: "1/1"}
    samples = [f"S{i + 1}" for i in range(N_SAMPLES)]
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    for info, genotypes in zip(VARIANT_INFO, GENOTYPE_ROWS):
        fields = [info["chrom"], str(info["pos"]), info["id"], "A", "G", "60", "PASS", ".", "GT"]
        lines.append("\t".join(fields + [calls[g] for g in genotypes]))
    return "\n".join(lines) + "\n"


@pytest.fixture
def dataset(tmp_path: Path) -> Dict[str, Path]:
    vcf = tmp_path / "failing.vcf"
    vcf.write_text(vcf_text())
    phenotypes = tmp_path / "phenotypes.tsv"
    rows = ["sample_id\tyield"] + [f"S{i + 1}\t{y}" for i, y in enumerate(PHENOTYPE)]
    phenotypes.write_text("\n".join(rows) + "\n")
    return {"vcf": vcf, "phenotypes": phenotypes}


KEEP_ALL = {"min_maf": 0.0, "max_missing": 1.0, "min_hwe_p": 0.0}


def test_run_gwas_writes_failed_tests_as_na_and_excludes_them_from_corrections(
    tmp_path: Path, dataset: Dict[str, Path], caplog: pytest.LogCaptureFixture
) -> None:
    out = tmp_path / "run_gwas"

    with caplog.at_level(logging.WARNING):
        result = run_gwas(
            dataset["vcf"],
            dataset["phenotypes"],
            {"model": "linear", "trait": "yield", "quality_control": KEEP_ALL, "significance_threshold": 1.0},
            output_dir=out,
        )

    assert result["status"] == "success"
    assert result["results"]["n_tests_failed"] == len(FAILING)
    table = {row["SNP"]: row for row in read_table(out / "summary_statistics.tsv")}
    assert list(table) == list(GENOTYPES)
    for snp in FAILING:
        assert [table[snp][col] for col in ("BETA", "SE", "P", "Q_FDR")] == ["NA"] * 4
    oracle_p = [stats.linregress(GENOTYPES[snp], PHENOTYPE).pvalue for snp in VALID]
    for snp, q in zip(VALID, bh(oracle_p)):
        assert float(table[snp]["Q_FDR"]) == pytest.approx(q, rel=1e-5)
    assert result["results"]["lambda_gc"] == pytest.approx(lambda_gc_oracle(oracle_p), abs=1e-4)
    # the q-values written to the rows of the returned results agree with the file
    annotated = {row["variant_id"]: row for row in result["results"]["association_results"]}
    assert all(is_nan(annotated[f"variant_{i}"]["q_value"]) for i, snp in enumerate(GENOTYPES) if snp in FAILING)
    assert {row["SNP"] for row in read_table(out / "significant_hits.tsv")} == set(VALID)
    assert json.loads((out / "results_summary.json").read_text())["n_variants_tested"] == len(VALID)
    assert not [r for r in caplog.records if "plot failed" in r.getMessage().lower()]


def test_execute_gwas_workflow_writes_failed_tests_as_na_and_excludes_them_from_corrections(
    tmp_path: Path, dataset: Dict[str, Path]
) -> None:
    out = tmp_path / "execute"
    config = {
        "vcf_path": str(dataset["vcf"]),
        "phenotype_path": str(dataset["phenotypes"]),
        "output_dir": str(out),
        "results_dir": str(out),
        "trait": "yield",
        "model": "linear",
        "quality_control": KEEP_ALL,
        "ld_pruning": {"r2_threshold": 1.0},  # keep every variant
        "structure": {"pca_components": 1},
        "significance_threshold": 1.0,
    }

    result = execute_gwas_workflow(config)

    assert result["status"] == "completed"
    assert result["outputs"]["n_tests_failed"] == len(FAILING)
    table = {row["SNP"]: row for row in read_table(out / "summary_statistics.tsv")}
    assert list(table) == list(GENOTYPES)
    for snp in FAILING:
        assert [table[snp][col] for col in ("BETA", "SE", "P", "Q_FDR")] == ["NA"] * 4
    present_p = [float(table[snp]["P"]) for snp in VALID]
    for snp, q in zip(VALID, bh(present_p)):  # the family is the tests that ran
        assert float(table[snp]["Q_FDR"]) == pytest.approx(q, rel=1e-4)
    outputs = result["outputs"]
    assert outputs["fdr"]["n_tests"] == len(VALID) and outputs["fdr"]["n_excluded"] == len(FAILING)
    assert outputs["bonferroni"]["n_tests"] == len(VALID)
    assert outputs["genomic_control"]["lambda_gc"] == pytest.approx(lambda_gc_oracle(present_p), rel=1e-4)
    assert {row["SNP"] for row in read_table(out / "significant_hits.tsv")} == set(VALID)
