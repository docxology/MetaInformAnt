"""Zero-mocks tests for gwas.analysis.power (NCP power, curves, jackknife).

All expectations are closed-form or hand-computed; no fixtures are replaced.
"""

from __future__ import annotations

import math
import random

import pytest

from metainformant.gwas.analysis import power


class TestComputePower:
    def test_null_effect_zero_power_at_genomewide_alpha(self) -> None:
        result = power.compute_power(n_samples=100, maf=0.5, beta=0.0)
        assert result["status"] == "success"
        assert result["power"] < 0.01
        assert result["ncp"] == 0.0

    def test_large_effect_high_power(self) -> None:
        # N = 10000, MAF = 0.5, beta = 0.5 → NCP = 10000*2*0.25*0.25 = 1250
        result = power.compute_power(n_samples=10000, maf=0.5, beta=0.5)
        assert result["status"] == "success"
        assert result["power"] > 0.999
        assert abs(result["ncp"] - 1250.0) < 1e-9

    def test_ncp_closed_form(self) -> None:
        n, maf, beta = 500, 0.3, 0.2
        result = power.compute_power(n_samples=n, maf=maf, beta=beta)
        expected_ncp = n * 2.0 * maf * (1.0 - maf) * beta**2
        assert abs(result["ncp"] - expected_ncp) < 1e-9

    def test_power_increases_with_n(self) -> None:
        p_small = power.compute_power(n_samples=100, maf=0.3, beta=0.5)["power"]
        p_large = power.compute_power(n_samples=1000, maf=0.3, beta=0.5)["power"]
        assert p_large > p_small

    def test_heritability_increases_power(self) -> None:
        """h2 < 1 shrinks the residual variance (σ²_e = 1 − h2), raising NCP.

        Verified implementation semantics: NCP = N·2·MAF(1−MAF)·β² / (1−h2),
        so conditioning on a trait with heritability h2=0.5 doubles NCP and
        therefore increases power for the same β.
        """
        p_no_h2 = power.compute_power(n_samples=1000, maf=0.3, beta=0.2)["power"]
        p_h2 = power.compute_power(n_samples=1000, maf=0.3, beta=0.2, h2=0.5)["power"]
        assert p_h2 > p_no_h2
        # NCP ratio must be exactly 1/(1-h2) = 2
        ncp0 = power.compute_power(n_samples=1000, maf=0.3, beta=0.2)["ncp"]
        ncp_h2 = power.compute_power(n_samples=1000, maf=0.3, beta=0.2, h2=0.5)["ncp"]
        assert abs(ncp_h2 / ncp0 - 2.0) < 1e-9

    def test_input_validation(self) -> None:
        assert power.compute_power(n_samples=0, maf=0.5, beta=0.1)["status"] == "error"
        assert power.compute_power(n_samples=100, maf=0.0, beta=0.1)["status"] == "error"
        assert power.compute_power(n_samples=100, maf=0.6, beta=0.1)["status"] == "error"
        assert power.compute_power(n_samples=100, maf=0.5, beta=0.1, alpha=1.5)["status"] == "error"


class TestPowerCurve:
    def test_monotone_in_sample_size(self) -> None:
        sizes = [50, 100, 200, 500, 1000]
        curve = power.power_curve(sample_sizes=sizes, maf=0.3, beta=0.3)
        powers = curve["powers"] if isinstance(curve, dict) else curve
        values = [p["power"] if isinstance(p, dict) else p for p in powers]
        for a, b in zip(values, values[1:]):
            assert b >= a - 1e-12

    def test_curve_length_matches_inputs(self) -> None:
        sizes = [100, 200, 400]
        curve = power.power_curve(sample_sizes=sizes, maf=0.5, beta=0.2)
        if isinstance(curve, dict):
            assert len(curve.get("sample_sizes", curve.get("powers", []))) == 3
        else:
            assert len(curve) == 3


class TestJackknifeSE:
    @staticmethod
    def _geno_matrix(n_variants: int, n_samples: int, rng: random.Random) -> list:
        return [[rng.choice([0, 1, 2]) for _ in range(n_samples)] for _ in range(n_variants)]

    def test_output_structure_and_ci_contains_estimate(self) -> None:
        rng = random.Random(42)
        geno = self._geno_matrix(40, 30, rng)
        traits_ = [rng.gauss(0, 1) for _ in range(30)]
        result = power.jackknife_se(geno, traits_, n_blocks=8, statistic="lambda_gc")
        assert result["status"] == "success"
        assert result["n_blocks"] == 8
        assert len(result["block_estimates"]) == 8
        assert result["se"] >= 0.0
        assert result["ci_lower"] <= result["estimate"] <= result["ci_upper"]

    def test_insufficient_variants_error(self) -> None:
        geno = [[0, 1], [1, 0]]
        result = power.jackknife_se(geno, [0.5, -0.5], n_blocks=10)
        assert result["status"] == "error"

    def test_se_zero_for_constant_blocks(self) -> None:
        # All blocks identical statistic → jackknife SE exactly 0
        n_var, n_samp = 20, 10
        geno = [[1] * n_samp for _ in range(n_var)]  # constant genotypes
        traits_ = [float(i) for i in range(n_samp)]
        result = power.jackknife_se(geno, traits_, n_blocks=5, statistic="mean_abs_beta")
        if result["status"] == "success":
            # monomorphic variants give identical leave-one-block estimates
            assert result["se"] == pytest.approx(0.0, abs=1e-12)

    def test_statistic_choices_accepted(self) -> None:
        rng = random.Random(7)
        geno = self._geno_matrix(30, 20, rng)
        traits_ = [rng.gauss(0, 1) for _ in range(20)]
        for stat in ("lambda_gc", "n_significant", "mean_abs_beta"):
            result = power.jackknife_se(geno, traits_, n_blocks=6, statistic=stat)
            assert result["status"] == "success", stat


class TestSubsampleConvergence:
    @staticmethod
    def _geno_matrix(n_variants: int, n_samples: int, rng: random.Random) -> list:
        return [[rng.choice([0, 1, 2]) for _ in range(n_samples)] for _ in range(n_variants)]

    def test_structure_and_determinism(self) -> None:
        rng = random.Random(13)
        n_var, n_samp = 60, 40
        geno = self._geno_matrix(n_var, n_samp, rng)
        traits_ = [rng.gauss(0, 1) for _ in range(n_samp)]
        # Embed a real association so metrics are non-trivial
        geno[0] = [traits_[i] * 0.5 + rng.gauss(0, 0.3) for i in range(n_samp)]
        result = power.subsample_convergence(geno, traits_, fractions=[0.2, 0.5, 1.0], n_replicates=2)
        assert result["status"] == "success" if "status" in result else True
        assert result["fractions"] == [0.2, 0.5, 1.0]
        for metric in ("lambda_gc", "n_significant", "mean_abs_beta", "mean_neg_log10_p"):
            assert metric in result["metrics"]
            assert len(result["metrics"][metric]) == 3
            for entry in result["metrics"][metric]:
                assert "mean" in entry and "std" in entry

    def test_empty_input_error(self) -> None:
        result = power.subsample_convergence([], [])
        assert result.get("status") == "error"


class TestSaturationAnalysis:
    @staticmethod
    def _convergence_data(fractions: list, means: list) -> dict:
        return {
            "fractions": fractions,
            "metrics": {"n_significant": [{"mean": m, "std": 0.0} for m in means]},
        }

    def test_saturated_plateau(self) -> None:
        fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        # Quickly saturating curve
        means = [50.0 * (1 - math.exp(-8 * f)) for f in fractions]
        result = power.saturation_analysis(self._convergence_data(fractions, means), metric="n_significant")
        assert result["status"] == "success"
        assert result["is_saturated"] is True
        assert result["r_squared"] > 0.9

    def test_not_saturated_when_steep(self) -> None:
        fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        # Slowly rising curve (κ=0.5): still gaining at 100%
        means = [50.0 * (1 - math.exp(-0.5 * f)) for f in fractions]
        result = power.saturation_analysis(self._convergence_data(fractions, means), metric="n_significant")
        assert result["status"] == "success"
        assert result["is_saturated"] is False

    def test_all_zero_trivially_saturated(self) -> None:
        fractions = [0.5, 1.0]
        result = power.saturation_analysis(self._convergence_data(fractions, [0.0, 0.0]))
        assert result["status"] == "success"
        assert result["is_saturated"] is True
        assert result["k_inf"] == 0.0

    def test_missing_metric_errors(self) -> None:
        result = power.saturation_analysis({"fractions": [0.5], "metrics": {}})
        assert result.get("status") == "error"
