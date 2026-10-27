"""Regression tests: summary-statistics output must not fabricate values.

Defect: ``write_summary_statistics`` and ``write_significant_hits`` wrote
``Q_FDR = 1`` and ``MAF = 0`` (and ``N = 0``, ``BETA = 0``, ``SE = 0``, ``P = 1``)
whenever a result lacked those keys. Results from ``run_*_gwas`` never carry a
q-value or an allele frequency, so every output file reported a genome-wide
significant hit with Q_FDR = 1 and a MAF of 0.

Independent oracles: Benjamini-Hochberg q-values worked out by hand, a separate
reference implementation of the step-up procedure, and
``scipy.stats.false_discovery_control``; minor allele frequencies counted by hand.
"""

from __future__ import annotations

import math
import random
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from metainformant.gwas.analysis.summary_stats import write_significant_hits, write_summary_statistics

MISSING = "NA"


def _read_table(path: Path) -> list[dict[str, str]]:
    """Parse a written TSV into one dict per row."""
    lines = path.read_text().strip().split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def _gwas_results(p_values: list[Any], n_samples: int = 100) -> list[dict[str, Any]]:
    """Results shaped like ``run_*_gwas`` output: beta, se, p_value, n_samples; no q-value, no MAF."""
    return [
        {"beta": 0.1 * (i + 1), "se": 0.05 + 0.01 * i, "p_value": p, "n_samples": n_samples}
        for i, p in enumerate(p_values)
    ]


def _variants(n: int) -> list[dict[str, Any]]:
    return [{"chrom": "1", "pos": 1000 * (i + 1), "id": f"rs{i}", "ref": "A", "alt": ["G"]} for i in range(n)]


def _benjamini_hochberg_reference(p_values: list[float]) -> list[float]:
    """Step-up BH written out longhand: q_(i) = min over j >= i of p_(j) * m / j."""
    m = len(p_values)
    order = sorted(range(m), key=lambda i: p_values[i])
    q = [0.0] * m
    running = 1.0
    for rank in range(m, 0, -1):
        index = order[rank - 1]
        running = min(running, p_values[index] * m / rank)
        q[index] = running
    return q


def _floats(rows: list[dict[str, str]], column: str) -> list[float]:
    return [float(row[column]) for row in rows]


class TestQFdr:
    """Q_FDR is Benjamini-Hochberg over the p-values, never a constant 1."""

    def test_run_gwas_shaped_results_get_hand_computed_q_values(self, tmp_path: Path) -> None:
        # Sorted: 1e-10, 1e-6, 0.04, 0.3, 0.5 with m = 5 -> q = 5e-10, 2.5e-6, 0.0666667, 0.375, 0.5
        p_values = [1e-6, 0.5, 1e-10, 0.04, 0.3]
        out = write_summary_statistics(_gwas_results(p_values), _variants(5), tmp_path / "stats.tsv")

        q = _floats(_read_table(out), "Q_FDR")
        assert q == pytest.approx([2.5e-6, 0.5, 5e-10, 0.2 / 3.0, 0.375], rel=1e-5)

    def test_highly_significant_hit_is_not_given_q_one(self, tmp_path: Path) -> None:
        """Regression: p = 1e-10 used to be written with Q_FDR = 1."""
        out = write_summary_statistics(_gwas_results([1e-10, 0.4, 0.8]), _variants(3), tmp_path / "stats.tsv")
        rows = _read_table(out)

        assert float(rows[0]["Q_FDR"]) < 1e-9
        assert all(row["Q_FDR"] != "1.000000e+00" for row in rows[:2])

    def test_monotonicity_is_enforced_in_sorted_order(self, tmp_path: Path) -> None:
        # Raw p * m / rank = 0.04, 0.06, 0.0533, 0.9 for sorted p = 0.01, 0.03, 0.04, 0.9;
        # the step-up minimum caps the second at the third.
        out = write_summary_statistics(_gwas_results([0.01, 0.04, 0.03, 0.9]), _variants(4), tmp_path / "stats.tsv")

        q = _floats(_read_table(out), "Q_FDR")
        assert q == pytest.approx([0.04, 0.16 / 3.0, 0.16 / 3.0, 0.9], rel=1e-5)

    def test_matches_reference_and_scipy_on_many_p_values(self, tmp_path: Path) -> None:
        scipy_stats = pytest.importorskip("scipy.stats")
        rng = random.Random(7)
        p_values = [10 ** -rng.uniform(0.0, 12.0) for _ in range(150)]
        p_values += [0.25, 0.25, 0.25, 1.0, 1e-300]  # ties, the boundary value 1, an extreme p

        out = write_summary_statistics(_gwas_results(p_values), _variants(len(p_values)), tmp_path / "stats.tsv")
        written = _floats(_read_table(out), "Q_FDR")

        assert written == pytest.approx(_benjamini_hochberg_reference(p_values), rel=2e-6)
        assert written == pytest.approx(list(scipy_stats.false_discovery_control(p_values, method="bh")), rel=2e-6)

    def test_missing_or_invalid_p_values_get_na_and_stay_out_of_the_family(self, tmp_path: Path) -> None:
        # Only 1e-6, 0.5 and 1e-10 are valid: m = 3, so q = 1.5e-6, 0.5, 3e-10.
        p_values = [1e-6, None, 0.5, float("nan"), 2.0, 1e-10]
        out = write_summary_statistics(_gwas_results(p_values), _variants(6), tmp_path / "stats.tsv")
        rows = _read_table(out)

        assert [rows[i]["Q_FDR"] for i in (1, 3, 4)] == [MISSING, MISSING, MISSING]  # None, nan, 2.0
        assert float(rows[0]["Q_FDR"]) == pytest.approx(1.5e-6, rel=1e-5)
        assert float(rows[2]["Q_FDR"]) == pytest.approx(0.5, rel=1e-5)
        assert float(rows[5]["Q_FDR"]) == pytest.approx(3e-10, rel=1e-5)
        assert rows[1]["P"] == MISSING
        assert rows[3]["P"] == MISSING

    @pytest.mark.parametrize("key", ["q_value", "q_fdr", "Q_FDR"])
    def test_a_q_value_the_result_carries_is_kept(self, tmp_path: Path, key: str) -> None:
        results = _gwas_results([1e-6, 0.5])
        results[0][key] = 0.123
        out = write_summary_statistics(results, _variants(2), tmp_path / "stats.tsv")

        rows = _read_table(out)
        assert float(rows[0]["Q_FDR"]) == pytest.approx(0.123)
        # The other result has none, so it gets the BH value over both p-values: 0.5 * 2 / 2.
        assert float(rows[1]["Q_FDR"]) == pytest.approx(0.5)

    def test_unusable_q_value_does_not_crash_and_falls_back(self, tmp_path: Path) -> None:
        results = _gwas_results([1e-6, 0.5])
        results[0]["q_value"] = None
        out = write_summary_statistics(results, _variants(2), tmp_path / "stats.tsv")

        assert float(_read_table(out)[0]["Q_FDR"]) == pytest.approx(2e-6, rel=1e-5)


class TestMaf:
    """MAF comes from the result, from the genotypes, or is NA; it is never 0 by default."""

    def test_missing_maf_is_na_not_zero(self, tmp_path: Path) -> None:
        out = write_summary_statistics(_gwas_results([1e-6, 0.5, 0.9]), _variants(3), tmp_path / "stats.tsv")

        assert [row["MAF"] for row in _read_table(out)] == [MISSING] * 3

    def test_maf_is_counted_from_genotypes_when_given(self, tmp_path: Path) -> None:
        genotypes = [
            [0, 0, 0, 0, 1, 2],  # 3 alt alleles of 12 -> 0.25
            [2, 2, 2, 1, -1, 2],  # 5 called: 9 alt alleles of 10 -> alt freq 0.9 -> MAF 0.1
            [-1, -1, -1, -1, -1, -1],  # nothing called -> no MAF
            [1, 1, 0, 0, 0, 0],  # supplied by the result; the genotypes must not override it
        ]
        results = _gwas_results([1e-6, 0.5, 0.9, 0.3])
        results[3]["maf"] = 0.33
        out = write_summary_statistics(results, _variants(4), tmp_path / "stats.tsv", genotypes=genotypes)

        assert [row["MAF"] for row in _read_table(out)] == ["0.2500", "0.1000", MISSING, "0.3300"]

    def test_numpy_genotype_matrix_is_accepted(self, tmp_path: Path) -> None:
        genotypes = np.array([[0, 1, 2, 1], [2, 2, 2, 2]])
        out = write_summary_statistics(
            _gwas_results([1e-3, 0.2]), _variants(2), tmp_path / "s.tsv", genotypes=genotypes
        )

        assert [row["MAF"] for row in _read_table(out)] == ["0.5000", "0.0000"]  # fixed allele: MAF is truly 0

    def test_uppercase_maf_key_is_honoured(self, tmp_path: Path) -> None:
        results = _gwas_results([1e-3])
        results[0]["MAF"] = 0.2
        out = write_summary_statistics(results, _variants(1), tmp_path / "stats.tsv", genotypes=[[2, 2, 2]])

        assert _read_table(out)[0]["MAF"] == "0.2000"

    def test_genotype_rows_must_match_results(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="same length"):
            write_summary_statistics(_gwas_results([0.1, 0.2]), _variants(2), tmp_path / "s.tsv", genotypes=[[0, 1]])
        with pytest.raises(ValueError, match="same length"):
            write_significant_hits(_gwas_results([0.1, 0.2]), _variants(2), tmp_path / "h.tsv", genotypes=[[0, 1]])


class TestOtherColumns:
    """No other column may invent a number either."""

    def test_missing_sample_size_is_na_not_zero(self, tmp_path: Path) -> None:
        results = [{"beta": 0.1, "se": 0.05, "p_value": 0.3}]
        out = write_summary_statistics(results, _variants(1), tmp_path / "stats.tsv")

        assert _read_table(out)[0]["N"] == MISSING

    def test_sample_size_is_read_from_either_key_and_written_as_an_integer(self, tmp_path: Path) -> None:
        results = [
            {"beta": 0.1, "se": 0.05, "p_value": 0.3, "n_samples": 100},
            {"beta": 0.1, "se": 0.05, "p_value": 0.3, "N": 250.0},
        ]
        out = write_summary_statistics(results, _variants(2), tmp_path / "stats.tsv")

        assert [row["N"] for row in _read_table(out)] == ["100", "250"]

    def test_missing_beta_se_and_p_are_na_not_zero_zero_one(self, tmp_path: Path) -> None:
        out = write_summary_statistics([{"n_samples": 10}], _variants(1), tmp_path / "stats.tsv")
        row = _read_table(out)[0]

        assert (row["BETA"], row["SE"], row["P"], row["Q_FDR"], row["MAF"]) == (MISSING,) * 5

    def test_none_and_non_finite_values_are_na(self, tmp_path: Path) -> None:
        results = [
            {"beta": None, "se": float("nan"), "p_value": 0.1, "n_samples": float("inf"), "maf": float("nan")},
            {"beta": float("-inf"), "se": 0.1, "p_value": 0.2, "n_samples": 5, "maf": 0.1},
        ]
        out = write_summary_statistics(results, _variants(2), tmp_path / "stats.tsv")
        rows = _read_table(out)

        assert (rows[0]["BETA"], rows[0]["SE"], rows[0]["N"], rows[0]["MAF"]) == (MISSING,) * 4
        assert rows[1]["BETA"] == MISSING
        assert rows[1]["SE"] == "0.1"

    def test_missing_position_is_na_and_missing_text_uses_the_vcf_dot(self, tmp_path: Path) -> None:
        out = write_summary_statistics(_gwas_results([0.2]), [{"alt": []}], tmp_path / "stats.tsv")
        row = _read_table(out)[0]

        assert row["POS"] == MISSING
        assert (row["CHR"], row["SNP"], row["REF"], row["ALT"]) == (".", ".", ".", ".")

    def test_real_values_are_still_written_exactly(self, tmp_path: Path) -> None:
        results = [{"beta": 0.5, "se": 0.1, "p_value": 1e-6, "q_value": 2e-6, "n_samples": 100, "maf": 0.3}]
        out = write_summary_statistics(results, _variants(1), tmp_path / "stats.tsv")

        assert _read_table(out)[0] == {
            "CHR": "1",
            "POS": "1000",
            "SNP": "rs0",
            "REF": "A",
            "ALT": "G",
            "BETA": "0.5",
            "SE": "0.1",
            "P": "1.000000e-06",
            "Q_FDR": "2.000000e-06",
            "N": "100",
            "MAF": "0.3000",
        }

    def test_na_is_read_back_as_missing_by_pandas(self, tmp_path: Path) -> None:
        pd = pytest.importorskip("pandas")
        out = write_summary_statistics(_gwas_results([1e-6, None, 0.5]), _variants(3), tmp_path / "stats.tsv")
        table = pd.read_csv(out, sep="\t")

        assert table["MAF"].isna().all()
        assert table["P"].isna().tolist() == [False, True, False]
        assert table["Q_FDR"].isna().tolist() == [False, True, False]
        assert table["N"].tolist() == [100, 100, 100]


class TestSignificantHits:
    """write_significant_hits follows the same conventions, computed over every test."""

    def test_q_values_reflect_all_tests_not_just_the_hits(self, tmp_path: Path) -> None:
        # BH over all three p-values: q = 3e-9, 1.5e-8 for the two hits (over only the hits it would be 2e-9, 1e-8).
        results = _gwas_results([1e-9, 1e-8, 0.5])
        out = write_significant_hits(results, _variants(3), tmp_path / "hits.tsv", threshold=1e-5)
        rows = _read_table(out)

        assert [row["SNP"] for row in rows] == ["rs0", "rs1"]
        assert _floats(rows, "Q_FDR") == pytest.approx([3e-9, 1.5e-8], rel=1e-5)
        assert [row["MAF"] for row in rows] == [MISSING, MISSING]

    def test_maf_comes_from_genotypes_and_hits_are_sorted_by_p(self, tmp_path: Path) -> None:
        genotypes = [[0, 1, 2, 1], [0, 0, 0, 1], [2, 2, 2, 2]]
        results = _gwas_results([1e-7, 1e-12, 0.4])
        out = write_significant_hits(results, _variants(3), tmp_path / "hits.tsv", threshold=1e-5, genotypes=genotypes)
        rows = _read_table(out)

        assert [row["SNP"] for row in rows] == ["rs1", "rs0"]
        assert [row["MAF"] for row in rows] == ["0.1250", "0.5000"]

    def test_results_without_a_usable_p_value_are_never_hits(self, tmp_path: Path) -> None:
        results = _gwas_results([None, float("nan"), 1e-9])
        results.append({"beta": 1.0, "se": 0.1})
        out = write_significant_hits(results, _variants(4), tmp_path / "hits.tsv", threshold=1e-5)

        assert [row["SNP"] for row in _read_table(out)] == ["rs2"]
        assert math.isclose(float(_read_table(out)[0]["Q_FDR"]), 1e-9, rel_tol=1e-5)  # one valid test: q = p
