"""Regression tests: an eQTL fit that cannot be made is missing data, not an effect of 0 with p = 1.

Defect: ``gwas.finemapping.eqtl._linear_regression`` returned ``p = 1.0`` for a fit it could not make (fewer than 3
complete observations, a constant genotype, a zero standard error), and ``cis_eqtl_scan`` recorded that as a test.
``eqtl_summary_stats`` then counted every failed fit in its Benjamini-Hochberg family as p = 1, which inflates every
adjusted p-value and can hide real eQTLs; ``eqtl_effect_sizes`` reported ``r_squared = 0.0`` for them.

A failed fit now carries NaN for beta, se and the p-value plus ``status == "error"`` and the reason; the summary
leaves it out of the family and reports it in ``n_tests_failed``; ``write_eqtl_results`` writes it as ``NA``.

Oracles: ``scipy.stats.linregress`` for every fit that can be made, ``scipy.stats.false_discovery_control``
(Benjamini-Hochberg) for the family, and plain numpy arithmetic for the residuals of the conditional analysis.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from metainformant.gwas.finemapping.eqtl import (
    cis_eqtl_scan,
    conditional_eqtl,
    eqtl_effect_sizes,
    eqtl_summary_stats,
)

N_SAMPLES = 12
SAMPLES = [f"s{i}" for i in range(N_SAMPLES)]
NAN = float("nan")


def make_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Two genes and four variants, all on one chromosome and inside one cis window.

    * ``gOK`` varies; ``gFLAT`` is constant (zero standard error for any variant that varies).
    * ``vOK1`` and ``vOK2`` vary; ``vCONST`` is heterozygous in every sample (MAF 0.5, but constant);
      ``vSPARSE`` has just two called samples.
    """
    rng = np.random.default_rng(5)
    expression = pd.DataFrame([rng.normal(10.0, 2.0, N_SAMPLES), np.full(N_SAMPLES, 5.0)], index=["gOK", "gFLAT"])
    expression.columns = SAMPLES

    sparse = np.full(N_SAMPLES, np.nan)
    sparse[[0, 1]] = 1.0
    genotypes = pd.DataFrame(
        [rng.integers(0, 3, N_SAMPLES), rng.integers(0, 3, N_SAMPLES), np.ones(N_SAMPLES), sparse],
        index=["vOK1", "vOK2", "vCONST", "vSPARSE"],
        columns=SAMPLES,
    ).astype(float)

    gene_positions = pd.DataFrame(
        {"gene_id": ["gOK", "gFLAT"], "chrom": ["1", "1"], "tss_position": [100_000, 150_000]}
    )
    variant_positions = pd.DataFrame(
        {"variant_id": list(genotypes.index), "chrom": ["1"] * 4, "position": [101_000, 102_000, 103_000, 104_000]}
    )
    return expression, genotypes, gene_positions, variant_positions


def scan() -> pd.DataFrame:
    expression, genotypes, gene_positions, variant_positions = make_inputs()
    return cis_eqtl_scan(expression, genotypes, gene_positions, variant_positions, maf_threshold=0.05)


def by_pair(frame: pd.DataFrame) -> Dict[tuple[str, str], pd.Series]:
    return {(row["gene_id"], row["variant_id"]): row for _, row in frame.iterrows()}


# ---------------------------------------------------------------------------
# cis_eqtl_scan
# ---------------------------------------------------------------------------

EXPECTED_FAILURES = {
    ("gOK", "vCONST"): "genotype is constant",
    ("gOK", "vSPARSE"): "fewer than 3 complete observations (2)",
    ("gFLAT", "vOK1"): "zero standard error",
    ("gFLAT", "vOK2"): "zero standard error",
    ("gFLAT", "vCONST"): "genotype is constant",
    ("gFLAT", "vSPARSE"): "fewer than 3 complete observations (2)",
}


def test_failed_fits_have_missing_statistics_a_status_and_a_reason() -> None:
    pairs = by_pair(scan())

    for pair, reason in EXPECTED_FAILURES.items():
        row = pairs[pair]
        assert row["status"] == "error", pair
        assert reason in row["error"], (pair, row["error"])
        assert all(math.isnan(row[column]) for column in ("beta", "se", "pvalue")), pair


def test_fits_that_can_be_made_match_linregress_and_carry_no_error() -> None:
    expression, genotypes, _, _ = make_inputs()
    pairs = by_pair(scan())

    for variant in ("vOK1", "vOK2"):
        row = pairs[("gOK", variant)]
        oracle = stats.linregress(genotypes.loc[variant].to_numpy(), expression.loc["gOK"].to_numpy())
        assert row["status"] == "success" and pd.isna(row["error"])
        assert row["beta"] == pytest.approx(oracle.slope, rel=1e-9)
        assert row["se"] == pytest.approx(oracle.stderr, rel=1e-9)
        assert row["pvalue"] == pytest.approx(oracle.pvalue, rel=1e-6)
    assert len(pairs) == 8  # every pair is still a row: failures are reported, not dropped


# ---------------------------------------------------------------------------
# eqtl_summary_stats
# ---------------------------------------------------------------------------


def test_the_summary_leaves_failed_fits_out_of_the_family_and_reports_them() -> None:
    results = scan()
    successes = results[results["status"] == "success"]
    oracle_q = stats.false_discovery_control(successes["pvalue"].to_numpy(), method="bh")

    summary = eqtl_summary_stats(results, fdr_threshold=0.05)

    assert summary["n_tests"] == len(successes) == 2
    assert summary["n_tests_failed"] == 6
    assert summary["n_eqtls"] == int((oracle_q < 0.05).sum())
    assert summary["mean_effect_size"] == pytest.approx(float(successes["beta"].abs().mean()))
    assert not math.isnan(summary["mean_effect_size"])


def test_failed_fits_do_not_inflate_the_benjamini_hochberg_family() -> None:
    present = [0.012, 0.02, 0.03, 0.04]
    results = pd.DataFrame(
        {
            "gene_id": [f"g{i}" for i in range(10)],
            "variant_id": [f"v{i}" for i in range(10)],
            "beta": [0.5] * 4 + [NAN] * 6,
            "pvalue": present + [NAN] * 6,
        }
    )

    summary = eqtl_summary_stats(results, fdr_threshold=0.05)

    # BH over the four tests that ran: every adjusted value is 0.048 or 0.04, so all four are eQTLs ...
    assert summary["n_eqtls"] == 4 and summary["n_egenes"] == 4
    # ... whereas counting the six failures as p = 1 (a family of ten) leaves none of them below 0.05
    assert not (stats.false_discovery_control(present + [1.0] * 6, method="bh")[:4] < 0.05).any()
    assert (summary["n_tests"], summary["n_tests_failed"]) == (4, 6)


def test_a_scan_without_any_successful_fit_has_an_empty_family() -> None:
    results = pd.DataFrame(
        {"gene_id": ["a", "b"], "variant_id": ["x", "y"], "beta": [NAN, NAN], "pvalue": [NAN, NAN], "status": "error"}
    )

    assert eqtl_summary_stats(results) == {"n_tests": 0, "n_tests_failed": 2, "n_egenes": 0, "n_eqtls": 0}


def test_summary_of_results_without_failures_is_unchanged() -> None:
    results = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3", "g4"],
            "variant_id": ["v1", "v2", "v3", "v4"],
            "beta": [1.0, -1.0, 1.0, -1.0],
            "pvalue": [0.9, 0.001, 0.8, 0.7],
        }
    )

    summary = eqtl_summary_stats(results, fdr_threshold=0.05)

    assert (summary["n_tests"], summary["n_tests_failed"], summary["n_eqtls"], summary["n_egenes"]) == (4, 0, 1, 1)
    assert summary["mean_effect_size"] == pytest.approx(1.0)
    assert eqtl_summary_stats(pd.DataFrame())["n_tests_failed"] == 0


# ---------------------------------------------------------------------------
# eqtl_effect_sizes and conditional_eqtl
# ---------------------------------------------------------------------------


def test_effect_sizes_of_a_failed_fit_are_missing_not_zero() -> None:
    expression, genotypes, _, _ = make_inputs()
    associations = pd.DataFrame({"gene_id": ["gOK", "gOK", "gFLAT"], "variant_id": ["vOK1", "vCONST", "vOK1"]})

    result = eqtl_effect_sizes(expression, genotypes, associations)

    rows = by_pair(result)
    good = rows[("gOK", "vOK1")]
    oracle = stats.linregress(genotypes.loc["vOK1"].to_numpy(), expression.loc["gOK"].to_numpy())
    assert good["status"] == "success" and pd.isna(good["error"])
    assert good["r_squared"] == pytest.approx(oracle.rvalue**2, rel=1e-9)
    for pair in (("gOK", "vCONST"), ("gFLAT", "vOK1")):
        row = rows[pair]
        assert row["status"] == "error" and row["error"]
        assert all(math.isnan(row[column]) for column in ("beta", "se", "pvalue", "r_squared")), pair


def test_conditional_analysis_marks_the_variants_it_could_not_test() -> None:
    expression, genotypes, _, _ = make_inputs()
    y = expression.loc["gOK"].to_numpy()
    lead = genotypes.loc["vOK1"].to_numpy()

    result = conditional_eqtl(y, genotypes, lead_variant_id="vOK1")

    rows = {row["variant_id"]: row for _, row in result.iterrows()}
    # oracle: residuals of the OLS fit of the expression on the lead variant, then a simple regression on vOK2
    slope = stats.linregress(lead, y).slope
    residuals = y - (y.mean() + slope * (lead - lead.mean()))
    oracle = stats.linregress(genotypes.loc["vOK2"].to_numpy(), residuals)
    assert rows["vOK2"]["status"] == "success"
    assert rows["vOK2"]["pvalue"] == pytest.approx(oracle.pvalue, rel=1e-6)
    for variant in ("vCONST", "vSPARSE"):
        assert rows[variant]["status"] == "error" and rows[variant]["error"]
        assert all(math.isnan(rows[variant][column]) for column in ("beta", "se", "pvalue"))


# ---------------------------------------------------------------------------
# The writer
# ---------------------------------------------------------------------------


def test_the_writer_prints_missing_values_as_na(tmp_path: Path) -> None:
    from metainformant.gwas.finemapping.eqtl import write_eqtl_results

    results = scan()
    path = write_eqtl_results(results, tmp_path / "nested" / "cis.tsv")

    lines = path.read_text().strip().split("\n")
    header = lines[0].split("\t")
    table = {(r["gene_id"], r["variant_id"]): r for r in (dict(zip(header, line.split("\t"))) for line in lines[1:])}
    assert len(table) == 8
    for pair in EXPECTED_FAILURES:
        assert [table[pair][column] for column in ("beta", "se", "pvalue")] == ["NA", "NA", "NA"], pair
        assert table[pair]["status"] == "error" and table[pair]["error"] != "NA"
    for variant in ("vOK1", "vOK2"):
        row = table[("gOK", variant)]
        assert all(math.isfinite(float(row[column])) for column in ("beta", "se", "pvalue"))
        assert row["status"] == "success"
    assert all(field != "" for line in lines for field in line.split("\t"))  # no empty field stands for "missing"


def test_a_table_written_with_na_reads_back_as_missing(tmp_path: Path) -> None:
    from metainformant.gwas.finemapping.eqtl import write_eqtl_results

    path = write_eqtl_results(scan(), tmp_path / "cis.tsv")

    frame = pd.read_csv(path, sep="\t")
    assert frame["pvalue"].dtype == float
    assert int(frame["pvalue"].isna().sum()) == 6 and int(frame["beta"].isna().sum()) == 6
    assert (frame.loc[frame["pvalue"].isna(), "status"] == "error").all()
    # a family built from the re-read table is the same as the one built from the in-memory scan
    assert eqtl_summary_stats(frame) == eqtl_summary_stats(scan())
