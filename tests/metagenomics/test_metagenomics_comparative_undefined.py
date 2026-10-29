"""Undefined tests are NaN (never p = 1.0), ANCOM has no p-value, and missing cross-validation is None.

Covers the Kruskal-Wallis and Welch tests, Benjamini-Hochberg with failed tests, the ANCOM-like W
statistic, the cross-validation accuracy of ``biomarker_discovery`` and the Monte Carlo of
``aldex2_like`` when NumPy is missing. Real implementations; SciPy and scikit-learn are oracles
guarded with ``pytest.importorskip``.
"""

from __future__ import annotations

import json
import math
import os
import random
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from metainformant.metagenomics.comparative import differential_abundance as da
from metainformant.metagenomics.comparative.differential_abundance import (
    biomarker_discovery,
    differential_abundance,
    effect_size_analysis,
)


def _package_root() -> str:
    import metainformant

    return os.path.dirname(os.path.dirname(os.path.abspath(metainformant.__file__)))


def _run_child(code: str, cwd: Path) -> str:
    """Run ``code`` in a child interpreter that imports exactly the package under test."""
    guarded = (
        "import os, metainformant\n"
        f"assert os.path.dirname(os.path.dirname(os.path.abspath(metainformant.__file__))) == {_package_root()!r}\n"
    ) + code
    run = subprocess.run(
        [sys.executable, "-c", guarded],
        env={**os.environ, "PYTHONPATH": _package_root()},
        cwd=cwd,
        capture_output=True,
        text=True,
        check=True,
        timeout=240,
    )
    return run.stdout.strip().splitlines()[-1]


def _bh_oracle(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg over the defined p-values only, written independently of the module."""
    defined = [(p, i) for i, p in enumerate(p_values) if not math.isnan(p)]
    defined.sort()
    m = len(defined)
    out = [math.nan] * len(p_values)
    running = 1.0
    for rank in range(m, 0, -1):
        p, i = defined[rank - 1]
        running = min(running, p * m / rank)
        out[i] = running
    return out


def _two_group_counts(seed: int = 0, n_per_group: int = 5, n_taxa: int = 7) -> tuple[list[list[int]], list[int]]:
    rng = random.Random(seed)
    counts, groups = [], []
    for group in (0, 1):
        for _ in range(n_per_group):
            row = [rng.randint(8, 60) for _ in range(n_taxa)]
            if group == 1:
                row[0] += 90
                row[1] = max(1, row[1] - 6)
            counts.append(row)
            groups.append(group)
    return counts, groups


class TestUndefinedKruskalWallis:
    """The Kruskal-Wallis statistic does not exist for identical values or fewer than two groups."""

    UNDEFINED = [
        [[1.0, 1.0], [1.0, 1.0]],  # every value identical
        [[1.0, 2.0]],  # a single group
        [[1.0, 2.0], []],  # an empty group
    ]

    @pytest.mark.parametrize("groups", UNDEFINED)
    def test_pure_implementation_is_nan(self, groups: list[list[float]]) -> None:
        assert math.isnan(da._kruskal_wallis_pure(groups))

    @pytest.mark.parametrize("groups", UNDEFINED)
    def test_dispatching_implementations_are_nan(self, groups: list[list[float]]) -> None:
        assert math.isnan(da._kruskal_wallis(groups))
        if len(groups) == 2:
            assert math.isnan(da._kruskal_wallis_two_group(groups[0], groups[1]))

    def test_effect_size_analysis_reports_nan_for_identical_samples(self) -> None:
        rows = effect_size_analysis([[5, 10, 15]] * 4, [0, 0, 1, 1], ["a", "b", "c"])
        assert len(rows) == 3
        assert all(math.isnan(r["p_value"]) for r in rows)
        assert all(r["lda_score"] == 0.0 for r in rows)


class TestBenjaminiHochbergWithFailedTests:
    """NaN p-values stay NaN and are left out of the family."""

    def test_nan_is_excluded_from_ranks_and_family_size(self) -> None:
        adjusted = da._benjamini_hochberg([0.01, math.nan, 0.04, 0.03, 0.5])
        assert math.isnan(adjusted[1])
        assert [adjusted[i] for i in (0, 2, 3, 4)] == pytest.approx([0.04, 0.16 / 3, 0.16 / 3, 0.5], rel=1e-12)

    def test_matches_independent_oracle(self) -> None:
        rng = random.Random(3)
        for _ in range(50):
            p_values = [math.nan if rng.random() < 0.25 else rng.random() ** 2 for _ in range(rng.randint(1, 30))]
            got = da._benjamini_hochberg(p_values)
            expected = _bh_oracle(p_values)
            assert len(got) == len(expected)
            for a, b in zip(got, expected):
                assert (math.isnan(a) and math.isnan(b)) or a == pytest.approx(b, rel=1e-12)

    def test_all_nan_and_empty(self) -> None:
        assert da._benjamini_hochberg([]) == []
        assert all(math.isnan(v) for v in da._benjamini_hochberg([math.nan, math.nan]))


class TestDifferentialAbundanceUndefinedTests:
    """A taxon that cannot be tested is NaN and last; it never becomes p = 1."""

    def test_absent_taxon_is_nan_and_sorted_last(self) -> None:
        counts, groups = _two_group_counts()
        zero_column = len(counts[0])
        counts = [row + [0] for row in counts]
        taxa = [f"t{j}" for j in range(zero_column)] + ["absent"]
        rows = differential_abundance(counts, groups, taxa, method="simple_deseq")
        by_taxon = {r["taxon"]: r for r in rows}
        assert rows[-1]["taxon"] == "absent"
        assert math.isnan(by_taxon["absent"]["p_value"])
        assert math.isnan(by_taxon["absent"]["adjusted_p"])
        assert math.isnan(by_taxon["absent"]["effect_size"])
        # the family is the tests that ran: 7 of 8 taxa
        defined = [r for r in rows if not math.isnan(r["p_value"])]
        assert len(defined) == zero_column
        expected = _bh_oracle([r["p_value"] for r in defined])
        assert [r["adjusted_p"] for r in defined] == pytest.approx(expected, rel=1e-12)
        assert [r["adjusted_p"] for r in defined] == sorted(r["adjusted_p"] for r in defined)

    def test_group_with_a_single_sample_is_not_reported_as_p_one(self) -> None:
        counts, _ = _two_group_counts(n_per_group=3)
        groups = [0, 1, 1, 1, 1, 1]
        for method in ("aldex2_like", "simple_deseq"):
            rows = differential_abundance(counts, groups, [f"t{j}" for j in range(7)], method=method, seed=1)
            assert all(math.isnan(r["p_value"]) for r in rows), method
            assert all(math.isnan(r["adjusted_p"]) for r in rows), method
            assert all(r["p_value"] != 1.0 for r in rows), method


class TestAncomLikeReportsAWStatistic:
    """ANCOM has no p-value; the old ``1 - W`` was passed through Benjamini-Hochberg as if it were one."""

    def test_rows_carry_w_and_no_p_value(self) -> None:
        counts, groups = _two_group_counts(seed=2)
        taxa = [f"t{j}" for j in range(7)]
        rows = differential_abundance(counts, groups, taxa, method="ancom_like")
        assert len(rows) == 7
        for row in rows:
            assert math.isnan(row["p_value"])
            assert math.isnan(row["adjusted_p"])
            assert isinstance(row["w_statistic"], int)
            assert row["n_ratio_tests"] == 6
            assert row["effect_size"] == pytest.approx(row["w_statistic"] / row["n_ratio_tests"])
            assert row["detected"] == (row["w_statistic"] >= 0.7 * row["n_ratio_tests"])

    def test_w_statistic_matches_an_independent_count(self) -> None:
        stats = pytest.importorskip("scipy.stats")
        counts, groups = _two_group_counts(seed=2)
        taxa = [f"t{j}" for j in range(7)]
        rows = {r["taxon"]: r for r in differential_abundance(counts, groups, taxa, method="ancom_like")}
        g0 = [i for i, g in enumerate(groups) if g == 0]
        g1 = [i for i, g in enumerate(groups) if g == 1]
        for i, taxon in enumerate(taxa):
            expected = 0
            for j in range(7):
                if i == j:
                    continue
                lr = [math.log(row[i] + 0.5) - math.log(row[j] + 0.5) for row in counts]
                p = float(stats.ttest_ind([lr[k] for k in g0], [lr[k] for k in g1], equal_var=False).pvalue)
                expected += p < 0.05
            assert rows[taxon]["w_statistic"] == expected

    def test_the_shifted_taxon_wins_and_results_are_ranked_by_w(self) -> None:
        counts, groups = _two_group_counts(seed=2)
        rows = differential_abundance(counts, groups, [f"t{j}" for j in range(7)], method="ancom_like")
        assert rows[0]["taxon"] == "t0"
        assert rows[0]["detected"] is True
        w = [r["w_statistic"] for r in rows]
        assert w == sorted(w, reverse=True)
        assert sum(r["detected"] for r in rows) < 7  # not everything is "detected"


class TestBiomarkerCrossValidation:
    """``cv_accuracy`` is None, with a reason, when no cross-validation ran (it used to be 0.0)."""

    @pytest.fixture(autouse=True)
    def _require_sklearn(self) -> None:
        pytest.importorskip("sklearn")

    TAXA = ["t0", "t1", "t2"]

    def test_one_sample_per_group_has_no_cross_validation(self) -> None:
        result = biomarker_discovery([[30, 2, 5], [3, 28, 6]], [0, 1], self.TAXA, seed=1)
        assert result["cv_accuracy"] is None
        assert "at least 2 samples" in result["cv_unavailable_reason"]
        assert result["model_summary"]["cv_folds"] is None
        assert result["model_summary"]["cv_scores"] == []
        assert set(result["importances"]) == set(self.TAXA)  # the forest itself is still fitted

    def test_a_singleton_group_has_no_cross_validation(self) -> None:
        counts = [[30, 2, 5], [28, 3, 6], [31, 1, 4], [29, 2, 5], [27, 4, 7], [30, 3, 6], [3, 28, 6]]
        result = biomarker_discovery(counts, [0, 0, 0, 0, 0, 0, 1], self.TAXA, seed=1)
        assert result["cv_accuracy"] is None
        assert "smallest group has 1" in result["cv_unavailable_reason"]
        assert result["model_summary"]["cv_scores"] == []

    def test_cross_validation_that_ran_has_a_value_and_no_reason(self) -> None:
        counts, groups = _two_group_counts(seed=4, n_taxa=3)
        result = biomarker_discovery(counts, groups, self.TAXA, cv_folds=4, seed=1)
        assert result["cv_unavailable_reason"] is None
        assert 0.0 <= result["cv_accuracy"] <= 1.0
        assert result["cv_accuracy"] == pytest.approx(sum(result["model_summary"]["cv_scores"]) / 4)
        assert result["model_summary"]["cv_folds"] == 4

    def test_fewer_than_two_folds_is_rejected(self) -> None:
        counts, groups = _two_group_counts(n_taxa=3)
        for folds in (1, 0, -3):
            with pytest.raises(ValueError, match="cv_folds"):
                biomarker_discovery(counts, groups, self.TAXA, cv_folds=folds)

    def test_a_single_group_has_nothing_to_discriminate(self) -> None:
        with pytest.raises(ValueError, match="at least two distinct"):
            biomarker_discovery([[1, 2, 3], [2, 3, 4], [3, 4, 5]], [1, 1, 1], self.TAXA)

    def test_without_scikit_learn_the_fallback_reports_no_accuracy(self, tmp_path: Path) -> None:
        counts, groups = _two_group_counts(seed=4, n_taxa=3)
        code = textwrap.dedent(f"""
            import json, sys
            sys.modules["sklearn"] = None  # makes `import sklearn` raise ImportError
            from metainformant.metagenomics.comparative import differential_abundance as da
            assert da.HAS_SKLEARN is False
            result = da.biomarker_discovery({counts!r}, {groups!r}, {self.TAXA!r})
            print(json.dumps(result))
            """)
        result = json.loads(_run_child(code, tmp_path))
        assert result["cv_accuracy"] is None
        assert "scikit-learn" in result["cv_unavailable_reason"]
        assert result["model_summary"]["cv_folds"] is None
        assert result["model_summary"]["cv_scores"] == []
        assert result["model_summary"]["method"] == "effect_size_ranking"
        assert max(result["importances"], key=result["importances"].get) == result["selected_taxa"][0]


class TestAldexWithoutNumpy:
    """Without NumPy the Dirichlet Monte Carlo is still run (it used to become the deterministic CLR)."""

    COUNTS = [
        [500, 10, 30],
        [600, 12, 28],
        [550, 9, 31],
        [520, 11, 29],
        [100, 40, 25],
        [90, 45, 27],
        [110, 38, 26],
        [95, 42, 24],
    ]
    GROUPS = [0, 0, 0, 0, 1, 1, 1, 1]

    def test_seeded_monte_carlo_without_numpy(self, tmp_path: Path) -> None:
        code = textwrap.dedent(f"""
            import json, sys
            sys.modules["numpy"] = None  # makes `import numpy` raise ImportError
            from metainformant.metagenomics.comparative import differential_abundance as da
            assert da.HAS_NUMPY is False
            counts, groups, taxa = {self.COUNTS!r}, {self.GROUPS!r}, ["a", "b", "c"]
            def run(seed, draws):
                rows = da.differential_abundance(counts, groups, taxa, n_monte_carlo=draws, seed=seed)
                return {{r["taxon"]: [r["p_value"], r["adjusted_p"]] for r in rows}}
            strong = run(7, 128)
            print(json.dumps({{
                "same_a": run(11, 2),
                "same_b": run(11, 2),
                "other": run(12, 2),
                "strong": {{taxon: values[1] for taxon, values in strong.items()}},
            }}))
            """)
        out = json.loads(_run_child(code, tmp_path))
        assert out["same_a"] == out["same_b"]  # seeded: bit-identical
        assert out["same_a"] != out["other"]  # a Monte Carlo: another seed perturbs the p-values
        assert out["strong"]["a"] < 0.05  # the strongly shifted taxon is still found
