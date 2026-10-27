"""Zero-mocks tests for gwas.analysis.ldsr (LD score computation + regression)."""

from __future__ import annotations

import json
import random
import subprocess
import sys
from typing import Any, Dict, List

import pytest

from metainformant.gwas.analysis import ldsr


class TestComputeLDScores:
    def test_identical_variants_full_ld_score(self) -> None:
        # Two identical dosage rows: r²=1 with each other → each LD score = 2 (self + other)
        geno = [[0, 1, 2, 1, 0], [0, 1, 2, 1, 0]]
        scores = ldsr.compute_ld_scores(geno)
        assert len(scores) == 2
        assert scores[0] == pytest.approx(2.0)
        assert scores[1] == pytest.approx(2.0)

    def test_uncorrelated_variants_ld_score_one(self) -> None:
        rng = random.Random(42)
        geno = [rng.choice([0, 1, 2]) for _ in range(500)]
        rows = [list(geno), [rng.choice([0, 1, 2]) for _ in range(500)]]
        scores = ldsr.compute_ld_scores(rows)
        # r²≈0 between rows → score ≈ 1.0 (self term only)
        assert scores[0] < 1.1

    def test_window_excludes_distant_variants(self) -> None:
        geno = [[0, 1, 2, 1, 0], [0, 1, 2, 1, 0]]
        positions = [1000, 5_000_000]  # 5 Mb apart
        scores = ldsr.compute_ld_scores(geno, positions=positions, window_kb=1000)
        # Outside window: only the self term contributes
        assert scores[0] == pytest.approx(1.0)
        assert scores[1] == pytest.approx(1.0)

    def test_single_variant_score_one(self) -> None:
        assert ldsr.compute_ld_scores([[0, 1, 2]]) == [1.0]

    def test_scores_nonnegative(self) -> None:
        rng = random.Random(7)
        geno = [[rng.choice([0, 1, 2]) for _ in range(30)] for _ in range(10)]
        scores = ldsr.compute_ld_scores(geno)
        assert all(s >= 1.0 - 1e-9 for s in scores)


class TestLDSRRegression:
    def test_recovers_h2_from_synthetic_chi2(self) -> None:
        """Generate chi2 = 1 + (N h² / M) ℓ_j + noise; regression must recover h²."""
        rng = random.Random(13)
        n_samples, n_variants = 10_000, 2000
        h2_true = 0.3
        ld_scores = [rng.uniform(1.0, 30.0) for _ in range(n_variants)]
        slope = n_samples * h2_true / n_variants
        chi2 = [1.0 + slope * s + rng.gauss(0, 0.5) for s in ld_scores]
        result = ldsr.ldsr_regression(chi2, ld_scores, n_samples=n_samples, n_variants=n_variants)
        assert result["method"] != "insufficient_data"
        assert abs(result["h2_snp"] - h2_true) < 0.05
        assert abs(result["intercept"] - 1.0) < 0.2
        assert result["h2_snp_se"] > 0

    def test_null_statistics_intercept_one_zero_h2(self) -> None:
        rng = random.Random(17)
        n_variants = 500
        ld_scores = [rng.uniform(1.0, 20.0) for _ in range(n_variants)]
        chi2 = [1.0 + rng.gauss(0, 0.3) for _ in range(n_variants)]  # pure confound-free null
        result = ldsr.ldsr_regression(chi2, ld_scores, n_samples=5000, n_variants=n_variants)
        assert abs(result["h2_snp"]) < 0.05
        assert abs(result["intercept"] - 1.0) < 0.15

    def test_intercept_constraint(self) -> None:
        rng = random.Random(23)
        n_variants = 300
        ld_scores = [rng.uniform(1.0, 15.0) for _ in range(n_variants)]
        chi2 = [1.2 + 0.01 * s for s in ld_scores]
        result = ldsr.ldsr_regression(chi2, ld_scores, n_samples=8000, n_variants=n_variants, intercept_constraint=1.2)
        assert result["intercept"] == pytest.approx(1.2)

    def test_genomic_control_lambda_uses_the_true_median(self) -> None:
        # 40 statistics 1.0, 1.1, ..., 4.9: the middle pair is 2.9 and 3.0, so the median is 2.95.
        chi2 = [1.0 + 0.1 * i for i in range(40)]
        ld_scores = [1.0 + i for i in range(40)]
        result = ldsr.ldsr_regression(chi2, ld_scores, n_samples=1000, n_variants=40)
        assert result["lambda_gc_median"] == pytest.approx(round(2.95 / 0.4549, 4))

    def test_insufficient_data(self) -> None:
        result = ldsr.ldsr_regression([2.0, 1.0], [5.0, 3.0], n_samples=1000, n_variants=100)
        assert result["method"] == "insufficient_data"
        assert result["h2_snp"] == 0.0

    def test_output_keys(self) -> None:
        rng = random.Random(29)
        ld = [rng.uniform(1.0, 10.0) for _ in range(50)]
        chi2 = [1.0 + 0.02 * s + rng.gauss(0, 0.2) for s in ld]
        result = ldsr.ldsr_regression(chi2, ld, n_samples=2000, n_variants=50)
        for key in ("h2_snp", "intercept", "mean_chi2", "method"):
            assert key in result
        assert result["mean_chi2"] == pytest.approx(sum(chi2) / len(chi2), rel=1e-3)


# Runs ``ldsr_regression`` in a fresh interpreter in which ``statsmodels.api`` cannot be imported,
# i.e. an environment without the optional dependency. ``ldsr_regression`` then takes its
# ``scipy.stats.linregress`` branch, which the in-process tests never reach because statsmodels is
# installed alongside the test suite.
_WITHOUT_STATSMODELS = """
import json
import sys

sys.modules["statsmodels"] = None
sys.modules["statsmodels.api"] = None

from metainformant.gwas.analysis import ldsr

assert not ldsr.HAS_STATSMODELS and ldsr.HAS_SCIPY
payload = json.load(sys.stdin)
result = ldsr.ldsr_regression(
    payload["chi2"],
    payload["ld"],
    n_samples=payload["n_samples"],
    n_variants=payload["n_variants"],
)
print(json.dumps(result))
"""


def _scipy_branch_result(chi2: List[float], ld_scores: List[float], n_samples: int, n_variants: int) -> Dict[str, Any]:
    payload = {"chi2": chi2, "ld": ld_scores, "n_samples": n_samples, "n_variants": n_variants}
    proc = subprocess.run(
        [sys.executable, "-c", _WITHOUT_STATSMODELS],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert proc.returncode == 0, proc.stderr
    parsed: Dict[str, Any] = json.loads(proc.stdout.strip().splitlines()[-1])
    return parsed


class TestLDSRScipyFallback:
    """The scipy branch must report the same OLS standard errors as an independent fit."""

    def test_scipy_path_intercept_se_matches_independent_ols(self) -> None:
        np = pytest.importorskip("numpy")
        sm = pytest.importorskip("statsmodels.api")

        rng = random.Random(13)
        n_samples, n_variants = 10_000, 400
        ld_scores = [rng.uniform(1.0, 30.0) for _ in range(n_variants)]
        slope = n_samples * 0.3 / n_variants
        chi2 = [1.0 + slope * s + rng.gauss(0, 0.5) for s in ld_scores]

        result = _scipy_branch_result(chi2, ld_scores, n_samples, n_variants)
        assert result["method"] == "scipy_linregress"

        x = np.asarray(ld_scores)
        y = np.maximum(np.asarray(chi2), 0.0)
        reference = sm.OLS(y, sm.add_constant(x)).fit()  # column 0 = intercept, column 1 = slope

        assert result["intercept"] == pytest.approx(float(reference.params[0]), abs=1e-4)
        assert result["slope"] == pytest.approx(float(reference.params[1]), abs=1e-6)
        assert result["slope_se"] == pytest.approx(float(reference.bse[1]), abs=1e-6)
        # SE(intercept) = SE(slope) * sqrt(mean(x^2)). Scaling by std(x) instead would drop the mean
        # term and roughly halve the value for strictly positive LD scores (0.025 vs 0.0509 here).
        assert result["intercept_se"] == pytest.approx(float(reference.bse[0]), abs=1e-4)
        assert result["intercept_se"] == pytest.approx(result["slope_se"] * float(np.sqrt(np.mean(x**2))), abs=1e-4)

    def test_scipy_and_statsmodels_branches_agree(self) -> None:
        pytest.importorskip("statsmodels.api")
        assert ldsr.HAS_STATSMODELS

        rng = random.Random(29)
        n_samples, n_variants = 5_000, 150
        ld_scores = [rng.uniform(1.0, 12.0) for _ in range(n_variants)]
        chi2 = [1.1 + 0.04 * s + rng.gauss(0, 0.3) for s in ld_scores]

        scipy_result = _scipy_branch_result(chi2, ld_scores, n_samples, n_variants)
        statsmodels_result = ldsr.ldsr_regression(chi2, ld_scores, n_samples=n_samples, n_variants=n_variants)

        assert statsmodels_result["method"] == "statsmodels_OLS"
        assert scipy_result["method"] == "scipy_linregress"
        for key in ("h2_snp", "h2_snp_se", "intercept", "intercept_se", "slope", "slope_se"):
            assert scipy_result[key] == pytest.approx(statsmodels_result[key], abs=1e-4), key
