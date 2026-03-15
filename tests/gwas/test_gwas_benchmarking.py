"""Tests for GWAS compute-time benchmarking module.

Tests scaling models, extrapolation, formatting, and data-driven estimate_runtime().
All tests use real functional methods — zero mocks.
"""

from __future__ import annotations

import pytest

from metainformant.gwas.analysis.benchmarking import (
    ComputeTimeEstimate,
    StepTiming,
    _estimate_step_proportions,
    _format_duration,
    extrapolate_full_genome_time,
    scaling_model,
)
from metainformant.gwas.data.config import estimate_runtime, _compute_scaling_factor


# ---------------------------------------------------------------------------
# scaling_model() unit tests
# ---------------------------------------------------------------------------


class TestScalingModel:
    """Test the scaling_model function for each complexity class."""

    def test_linear_nm_scaling(self) -> None:
        """O(n·m): doubling both dimensions → 4× time."""
        factor = scaling_model("n_m", pilot_n=100, pilot_m=1000, target_n=200, target_m=2000)
        assert factor == pytest.approx(4.0)

    def test_linear_nm_same_dimensions(self) -> None:
        """Same dimensions → 1× scaling."""
        factor = scaling_model("n_m", pilot_n=100, pilot_m=1000, target_n=100, target_m=1000)
        assert factor == pytest.approx(1.0)

    def test_quadratic_n2m_scaling(self) -> None:
        """O(n²·m): doubling samples → 4× from n², doubling variants → 2× from m → 8× total."""
        factor = scaling_model("n2_m", pilot_n=100, pilot_m=1000, target_n=200, target_m=2000)
        assert factor == pytest.approx(8.0)

    def test_linear_m_scaling(self) -> None:
        """O(m): only variant count matters."""
        factor = scaling_model("m", pilot_n=100, pilot_m=1000, target_n=5000, target_m=3000)
        assert factor == pytest.approx(3.0)

    def test_mk2_scaling(self) -> None:
        """O(m·k²): scaling with PCA components."""
        # 2× variants, 2× components → 2 × 4 = 8× ... but then divided by pilot k²=100
        factor = scaling_model("m_k2", pilot_n=100, pilot_m=1000, target_n=100, target_m=2000, k_pilot=10, k_target=20)
        # = (2000 * 400) / (1000 * 100) = 800000 / 100000 = 8.0
        assert factor == pytest.approx(8.0)

    def test_k3_scaling_subcubic(self) -> None:
        """O(k³) approximated as m^1.5 ratio."""
        factor = scaling_model("k3", pilot_n=100, pilot_m=100, target_n=100, target_m=400)
        # (400/100)^1.5 = 4^1.5 = 8.0
        assert factor == pytest.approx(8.0)

    def test_unknown_model_fallback(self) -> None:
        """Unknown model falls back to O(m)."""
        factor = scaling_model("unknown_model", pilot_n=100, pilot_m=1000, target_n=100, target_m=2000)
        assert factor == pytest.approx(2.0)

    def test_zero_pilot_raises(self) -> None:
        """Zero pilot dimensions should raise ValueError."""
        with pytest.raises(ValueError, match="positive"):
            scaling_model("n_m", pilot_n=0, pilot_m=1000, target_n=100, target_m=2000)


# ---------------------------------------------------------------------------
# extrapolate_full_genome_time() tests
# ---------------------------------------------------------------------------


class TestExtrapolateFullGenomeTime:
    """Test extrapolation from pilot timings to full genome."""

    def test_basic_extrapolation(self) -> None:
        """Extrapolation should scale times by the appropriate factor."""
        pilot_timings = [
            StepTiming(step_name="association_testing", elapsed_seconds=10.0, n_samples=100, n_variants=1000),
            StepTiming(step_name="qc_filters", elapsed_seconds=2.0, n_samples=100, n_variants=1000),
        ]

        est = extrapolate_full_genome_time(pilot_timings, target_n_samples=200, target_n_variants=2000)

        assert isinstance(est, ComputeTimeEstimate)
        assert est.total_seconds > 0
        assert "association_testing" in est.per_step
        assert "qc_filters" in est.per_step
        # Association: O(n·m) → 4× from doubling both
        assert est.per_step["association_testing"] == pytest.approx(40.0)
        # QC: same O(n·m) → 4×
        assert est.per_step["qc_filters"] == pytest.approx(8.0)

    def test_empty_timings_returns_zero(self) -> None:
        """Empty pilot timings → zero estimate."""
        est = extrapolate_full_genome_time([], target_n_samples=1000, target_n_variants=1000000)
        assert est.total_seconds == 0.0

    def test_same_dimensions_identity(self) -> None:
        """Same pilot and target dimensions → same times."""
        pilot_timings = [
            StepTiming(step_name="parse_vcf", elapsed_seconds=5.0, n_samples=100, n_variants=500),
        ]
        est = extrapolate_full_genome_time(pilot_timings, target_n_samples=100, target_n_variants=500)
        assert est.per_step["parse_vcf"] == pytest.approx(5.0)

    def test_custom_models_override(self) -> None:
        """Custom models should override defaults."""
        pilot_timings = [
            StepTiming(step_name="association_testing", elapsed_seconds=10.0, n_samples=100, n_variants=1000),
        ]
        # Force O(m) instead of default O(n·m)
        est = extrapolate_full_genome_time(
            pilot_timings, target_n_samples=200, target_n_variants=2000,
            custom_models={"association_testing": "m"},
        )
        # O(m): 2000/1000 = 2× → 20.0s
        assert est.per_step["association_testing"] == pytest.approx(20.0)

    def test_summary_output(self) -> None:
        """Summary should be a readable string."""
        pilot_timings = [
            StepTiming(step_name="association_testing", elapsed_seconds=10.0, n_samples=100, n_variants=1000),
        ]
        est = extrapolate_full_genome_time(pilot_timings, target_n_samples=200, target_n_variants=2000)
        summary = est.summary()
        assert "Estimated total runtime" in summary
        assert "association_testing" in summary


# ---------------------------------------------------------------------------
# Helper tests
# ---------------------------------------------------------------------------


class TestFormatDuration:
    """Test duration formatting."""

    def test_seconds(self) -> None:
        assert _format_duration(45.2) == "45.2s"

    def test_minutes(self) -> None:
        assert _format_duration(90.0) == "1m 30s"

    def test_hours(self) -> None:
        assert _format_duration(3661.0) == "1h 1m 1s"

    def test_zero(self) -> None:
        assert _format_duration(0.0) == "0.0s"


class TestStepProportions:
    """Test step proportion estimation fallback."""

    def test_proportions_sum_to_one(self) -> None:
        steps = ["parse_vcf", "qc_filters", "association_testing"]
        proportions = _estimate_step_proportions(steps)
        assert sum(proportions.values()) == pytest.approx(1.0)

    def test_unknown_step_handled(self) -> None:
        proportions = _estimate_step_proportions(["parse_vcf", "unknown_step"])
        assert "unknown_step" in proportions
        assert sum(proportions.values()) == pytest.approx(1.0)


class TestStepTiming:
    """Test StepTiming dataclass."""

    def test_throughput(self) -> None:
        t = StepTiming(step_name="test", elapsed_seconds=10.0, n_samples=100, n_variants=500)
        assert t.throughput == pytest.approx(5000.0)

    def test_zero_elapsed_throughput(self) -> None:
        t = StepTiming(step_name="test", elapsed_seconds=0.0, n_samples=100, n_variants=500)
        assert t.throughput == 0.0


# ---------------------------------------------------------------------------
# estimate_runtime() data-driven tests (config.py)
# ---------------------------------------------------------------------------


class TestEstimateRuntime:
    """Test data-driven estimate_runtime in config.py."""

    def test_default_config_returns_estimate(self) -> None:
        """Default config (no dimensions) should still return a valid estimate."""
        config = {"threads": 4, "pca_components": 10}
        est = estimate_runtime(config)
        assert est["estimated_hours"] > 0
        assert est["data_driven"] is False
        assert "per_step_seconds" in est

    def test_data_driven_with_dimensions(self) -> None:
        """Providing dimensions should enable data-driven mode."""
        config = {"threads": 4, "pca_components": 10}
        est = estimate_runtime(config, n_samples=500, n_variants=100_000)
        assert est["data_driven"] is True
        assert est["estimated_seconds"] > 0
        assert "association_testing" in est["per_step_seconds"]

    def test_larger_data_takes_longer(self) -> None:
        """More samples/variants should yield higher runtime estimates."""
        config = {"threads": 4}
        est_small = estimate_runtime(config, n_samples=100, n_variants=1000)
        est_large = estimate_runtime(config, n_samples=1000, n_variants=100_000)
        assert est_large["estimated_seconds"] > est_small["estimated_seconds"]

    def test_more_threads_reduces_time(self) -> None:
        """More threads should produce shorter estimates via Amdahl's law."""
        config_1t = {"threads": 1, "pca_components": 10}
        config_8t = {"threads": 8, "pca_components": 10}
        est_1 = estimate_runtime(config_1t, n_samples=500, n_variants=10000)
        est_8 = estimate_runtime(config_8t, n_samples=500, n_variants=10000)
        assert est_8["estimated_seconds"] < est_1["estimated_seconds"]

    def test_limiting_factors_reported(self) -> None:
        """Large data should report limiting factors."""
        config = {"threads": 4}
        est = estimate_runtime(config, n_samples=2000, n_variants=5_000_000)
        assert len(est["limiting_factors"]) > 0

    def test_plots_disabled_reduces_time(self) -> None:
        """Disabling plots should reduce total time."""
        config_plots = {"threads": 1, "create_plots": True}
        config_noplots = {"threads": 1, "create_plots": False}
        est_p = estimate_runtime(config_plots, n_samples=100, n_variants=1000)
        est_np = estimate_runtime(config_noplots, n_samples=100, n_variants=1000)
        assert est_np["estimated_seconds"] < est_p["estimated_seconds"]


class TestComputeScalingFactor:
    """Test _compute_scaling_factor from config.py."""

    def test_nm_model(self) -> None:
        factor = _compute_scaling_factor("n_m", 100, 1000, 200, 2000)
        assert factor == pytest.approx(4.0)

    def test_n2m_model(self) -> None:
        factor = _compute_scaling_factor("n2_m", 100, 1000, 200, 2000)
        assert factor == pytest.approx(8.0)

    def test_m_model(self) -> None:
        factor = _compute_scaling_factor("m", 100, 1000, 200, 3000)
        assert factor == pytest.approx(3.0)

    def test_zero_pilot_returns_one(self) -> None:
        factor = _compute_scaling_factor("n_m", 0, 1000, 200, 2000)
        assert factor == 1.0
