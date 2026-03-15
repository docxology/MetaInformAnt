"""Tests for phenotype morphological: Measurement, MorphometricProfile, allometry, comparison.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import math

import pytest

from metainformant.core.utils.errors import ValidationError
from metainformant.phenotype.morphological.measurement import Measurement
from metainformant.phenotype.morphological.profile import (
    MorphometricProfile,
    allometric_regression,
    coefficient_of_variation,
    compare_profiles,
    summary_statistics,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_measurement(name="head_length", value=2.5, unit="mm"):
    return Measurement(name=name, value=value, unit=unit)


def _make_profile(specimen_id="ANT001", measurements=None):
    if measurements is None:
        measurements = [
            Measurement(name="head_length", value=2.5, unit="mm"),
            Measurement(name="head_width", value=2.0, unit="mm"),
            Measurement(name="thorax_length", value=3.0, unit="mm"),
            Measurement(name="femur_length", value=1.8, unit="mm"),
        ]
    return MorphometricProfile(specimen_id=specimen_id, measurements=measurements)


def _make_profiles(n=10, seed=42):
    """Create n profiles with slightly varying measurements."""
    import random

    random.seed(seed)
    profiles = []
    for i in range(n):
        scale = 1.0 + random.gauss(0, 0.1)
        measurements = [
            Measurement(name="head_length", value=2.5 * scale, unit="mm"),
            Measurement(name="head_width", value=2.0 * scale, unit="mm"),
            Measurement(name="thorax_length", value=3.0 * scale, unit="mm"),
            Measurement(name="femur_length", value=1.8 * scale, unit="mm"),
        ]
        profiles.append(MorphometricProfile(specimen_id=f"ANT{i:03d}", measurements=measurements))
    return profiles


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------


class TestMeasurement:
    def test_basic_creation(self):
        m = Measurement(name="head_length", value=2.5)
        assert m.name == "head_length"
        assert m.value == 2.5
        assert m.unit == "mm"
        assert m.uncertainty is None
        assert m.description is None

    def test_with_uncertainty(self):
        m = Measurement(name="head_length", value=2.5, uncertainty=0.1, description="HL")
        assert m.uncertainty == 0.1
        assert m.description == "HL"

    def test_convert_mm_to_cm(self):
        m = Measurement(name="length", value=25.0, unit="mm")
        converted = m.convert("cm")
        assert converted.unit == "cm"
        assert converted.value == pytest.approx(2.5)

    def test_convert_cm_to_mm(self):
        m = Measurement(name="length", value=2.5, unit="cm")
        converted = m.convert("mm")
        assert converted.unit == "mm"
        assert converted.value == pytest.approx(25.0)

    def test_convert_mm_to_m(self):
        m = Measurement(name="length", value=1000.0, unit="mm")
        converted = m.convert("m")
        assert converted.value == pytest.approx(1.0)

    def test_convert_mm_to_um(self):
        m = Measurement(name="length", value=1.0, unit="mm")
        converted = m.convert("um")
        assert converted.value == pytest.approx(1000.0)

    def test_convert_um_to_mm(self):
        m = Measurement(name="length", value=500.0, unit="um")
        converted = m.convert("mm")
        assert converted.value == pytest.approx(0.5)

    def test_convert_same_unit(self):
        m = Measurement(name="length", value=5.0, unit="mm")
        converted = m.convert("mm")
        assert converted is m  # Same object returned

    def test_convert_preserves_uncertainty(self):
        m = Measurement(name="length", value=10.0, unit="mm", uncertainty=0.5)
        converted = m.convert("cm")
        assert converted.uncertainty == pytest.approx(0.05)

    def test_convert_unsupported_unit(self):
        m = Measurement(name="length", value=5.0, unit="mm")
        with pytest.raises(ValueError, match="Unsupported unit"):
            m.convert("inches")

    def test_convert_preserves_name(self):
        m = Measurement(name="head_length", value=2.5, unit="mm")
        converted = m.convert("cm")
        assert converted.name == "head_length"


# ---------------------------------------------------------------------------
# MorphometricProfile - basics
# ---------------------------------------------------------------------------


class TestMorphometricProfile:
    def test_create_empty(self):
        p = MorphometricProfile(specimen_id="ANT001")
        assert p.specimen_id == "ANT001"
        assert p.n_measurements == 0
        assert p.measurement_names == []

    def test_create_with_measurements(self):
        p = _make_profile()
        assert p.specimen_id == "ANT001"
        assert p.n_measurements == 4

    def test_add_measurement(self):
        p = MorphometricProfile(specimen_id="ANT001")
        p.add(Measurement(name="head_length", value=2.5))
        assert p.n_measurements == 1
        assert "head_length" in p.measurement_names

    def test_add_replaces_existing(self):
        p = MorphometricProfile(specimen_id="ANT001")
        p.add(Measurement(name="head_length", value=2.5))
        p.add(Measurement(name="head_length", value=3.0))
        assert p.n_measurements == 1
        assert p.get("head_length").value == 3.0

    def test_get_existing(self):
        p = _make_profile()
        m = p.get("head_length")
        assert m is not None
        assert m.value == 2.5

    def test_get_missing(self):
        p = _make_profile()
        assert p.get("nonexistent") is None

    def test_values_dict(self):
        p = _make_profile()
        vals = p.values_dict()
        assert isinstance(vals, dict)
        assert vals["head_length"] == 2.5
        assert vals["head_width"] == 2.0
        assert len(vals) == 4

    def test_metadata(self):
        p = MorphometricProfile(
            specimen_id="ANT001",
            metadata={"species": "Formica rufa", "collector": "EO Wilson"},
        )
        assert p.metadata["species"] == "Formica rufa"

    def test_default_metadata(self):
        p = MorphometricProfile(specimen_id="ANT001")
        assert p.metadata == {}


# ---------------------------------------------------------------------------
# MorphometricProfile - calculate_index
# ---------------------------------------------------------------------------


class TestCalculateIndex:
    def test_cephalic_index(self):
        p = _make_profile()
        ci = p.calculate_index("CI", "head_width", "head_length")
        assert ci == pytest.approx(80.0)  # (2.0 / 2.5) * 100

    def test_index_missing_measurement(self):
        p = _make_profile()
        with pytest.raises(ValidationError, match="Missing measurements"):
            p.calculate_index("test", "head_width", "nonexistent")

    def test_index_zero_denominator(self):
        p = MorphometricProfile(
            specimen_id="test",
            measurements=[
                Measurement(name="a", value=5.0),
                Measurement(name="b", value=0.0),
            ],
        )
        with pytest.raises(ValueError, match="Denominator is zero"):
            p.calculate_index("test", "a", "b")

    def test_index_different_units_auto_converts(self):
        p = MorphometricProfile(
            specimen_id="test",
            measurements=[
                Measurement(name="a", value=10.0, unit="mm"),
                Measurement(name="b", value=1.0, unit="cm"),  # = 10mm
            ],
        )
        idx = p.calculate_index("test", "a", "b")
        # 10mm / 10mm * 100 = 100
        assert idx == pytest.approx(100.0)


# ---------------------------------------------------------------------------
# MorphometricProfile - geometric_mean_size
# ---------------------------------------------------------------------------


class TestGeometricMeanSize:
    def test_positive_values(self):
        p = _make_profile()
        gm = p.geometric_mean_size()
        assert gm > 0
        # GM of 2.5, 2.0, 3.0, 1.8
        expected = math.exp((math.log(2.5) + math.log(2.0) + math.log(3.0) + math.log(1.8)) / 4)
        assert gm == pytest.approx(expected)

    def test_empty_profile(self):
        p = MorphometricProfile(specimen_id="empty")
        assert p.geometric_mean_size() == 0.0

    def test_profile_with_zero_value(self):
        p = MorphometricProfile(
            specimen_id="test",
            measurements=[
                Measurement(name="a", value=0.0),
                Measurement(name="b", value=5.0),
            ],
        )
        gm = p.geometric_mean_size()
        # Only positive values used: just 5.0
        assert gm == pytest.approx(5.0)


# ---------------------------------------------------------------------------
# MorphometricProfile - size_corrected_values
# ---------------------------------------------------------------------------


class TestSizeCorrectedValues:
    def test_basic_correction(self):
        p = _make_profile()
        corrected = p.size_corrected_values()
        assert len(corrected) == 4
        # All values should be around 1.0 (for similar-sized measurements)
        for v in corrected.values():
            assert v > 0

    def test_empty_profile_returns_zeros(self):
        p = MorphometricProfile(specimen_id="empty")
        corrected = p.size_corrected_values()
        assert corrected == {}

    def test_uniform_measurements_give_one(self):
        p = MorphometricProfile(
            specimen_id="uniform",
            measurements=[
                Measurement(name="a", value=5.0),
                Measurement(name="b", value=5.0),
                Measurement(name="c", value=5.0),
            ],
        )
        corrected = p.size_corrected_values()
        for v in corrected.values():
            assert v == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# coefficient_of_variation
# ---------------------------------------------------------------------------


class TestCoefficientOfVariation:
    def test_basic_cv(self):
        profiles = _make_profiles(n=10)
        cv = coefficient_of_variation(profiles, "head_length")
        assert cv is not None
        assert cv > 0

    def test_identical_values_zero_cv(self):
        profiles = [
            MorphometricProfile(
                specimen_id=f"s{i}",
                measurements=[Measurement(name="x", value=5.0)],
            )
            for i in range(5)
        ]
        cv = coefficient_of_variation(profiles, "x")
        assert cv == pytest.approx(0.0)

    def test_missing_measurement_returns_none(self):
        profiles = _make_profiles(n=3)
        cv = coefficient_of_variation(profiles, "nonexistent")
        assert cv is None

    def test_single_profile_returns_none(self):
        profiles = [_make_profile()]
        cv = coefficient_of_variation(profiles, "head_length")
        assert cv is None  # Need at least 2

    def test_cv_is_percentage(self):
        profiles = _make_profiles(n=20)
        cv = coefficient_of_variation(profiles, "head_length")
        # CV should be expressed as percentage (typically 5-15% for biological data)
        assert 0 < cv < 100


# ---------------------------------------------------------------------------
# allometric_regression
# ---------------------------------------------------------------------------


class TestAllometricRegression:
    def test_basic_regression(self):
        profiles = _make_profiles(n=20)
        result = allometric_regression(profiles, "head_length", "head_width")
        assert "slope" in result
        assert "intercept" in result
        assert "r_squared" in result
        assert "n" in result
        assert "allometry_type" in result
        assert result["n"] == 20

    def test_isometric_scaling(self):
        # When y scales 1:1 with x, slope ~= 1
        profiles = _make_profiles(n=20)
        result = allometric_regression(profiles, "head_length", "head_width")
        # Since both scale with the same factor, slope should be near 1
        assert result["allometry_type"] == "isometric"
        assert abs(result["slope"] - 1.0) < 0.1

    def test_r_squared_high_for_correlated(self):
        profiles = _make_profiles(n=20)
        result = allometric_regression(profiles, "head_length", "head_width")
        assert result["r_squared"] > 0.8

    def test_insufficient_data(self):
        profiles = _make_profiles(n=2)
        with pytest.raises(ValidationError, match="at least 3"):
            allometric_regression(profiles, "head_length", "head_width")

    def test_missing_measurement_reduces_n(self):
        profiles = _make_profiles(n=5)
        # Add a profile without the x measurement
        p = MorphometricProfile(
            specimen_id="missing",
            measurements=[Measurement(name="head_width", value=2.0)],
        )
        profiles.append(p)
        result = allometric_regression(profiles, "head_length", "head_width")
        assert result["n"] == 5  # The incomplete profile is excluded

    def test_positive_allometry(self):
        # Create data where y scales faster than x
        import random

        random.seed(42)
        profiles = []
        for i in range(20):
            x_val = 1.0 + i * 0.2
            y_val = x_val**2.0 * (1 + random.gauss(0, 0.02))
            profiles.append(
                MorphometricProfile(
                    specimen_id=f"s{i}",
                    measurements=[
                        Measurement(name="body_size", value=x_val),
                        Measurement(name="head_size", value=y_val),
                    ],
                )
            )
        result = allometric_regression(profiles, "body_size", "head_size")
        assert result["slope"] > 1.0
        assert result["allometry_type"] == "positive"

    def test_negative_allometry(self):
        # Create data where y scales slower than x
        import random

        random.seed(42)
        profiles = []
        for i in range(20):
            x_val = 1.0 + i * 0.2
            y_val = x_val**0.5 * (1 + random.gauss(0, 0.02))
            profiles.append(
                MorphometricProfile(
                    specimen_id=f"s{i}",
                    measurements=[
                        Measurement(name="body_size", value=x_val),
                        Measurement(name="appendage", value=y_val),
                    ],
                )
            )
        result = allometric_regression(profiles, "body_size", "appendage")
        assert result["slope"] < 1.0
        assert result["allometry_type"] == "negative"


# ---------------------------------------------------------------------------
# compare_profiles
# ---------------------------------------------------------------------------


class TestCompareProfiles:
    def test_basic_comparison(self):
        p_a = _make_profile(specimen_id="A")
        p_b = _make_profile(specimen_id="B")
        result = compare_profiles(p_a, p_b)
        assert result["specimen_a"] == "A"
        assert result["specimen_b"] == "B"
        assert result["shared_measurements"] == 4

    def test_identical_profiles_zero_diff(self):
        p_a = _make_profile(specimen_id="A")
        p_b = _make_profile(specimen_id="B")
        result = compare_profiles(p_a, p_b)
        for name, comp in result["comparisons"].items():
            assert comp["difference"] == pytest.approx(0.0)
            assert comp["percent_difference"] == pytest.approx(0.0)

    def test_different_profiles(self):
        p_a = MorphometricProfile(
            specimen_id="A",
            measurements=[Measurement(name="head_length", value=2.5, unit="mm")],
        )
        p_b = MorphometricProfile(
            specimen_id="B",
            measurements=[Measurement(name="head_length", value=3.0, unit="mm")],
        )
        result = compare_profiles(p_a, p_b)
        comp = result["comparisons"]["head_length"]
        assert comp["value_a"] == 2.5
        assert comp["value_b"] == 3.0
        assert comp["difference"] == pytest.approx(-0.5)

    def test_unique_measurements(self):
        p_a = MorphometricProfile(
            specimen_id="A",
            measurements=[
                Measurement(name="head_length", value=2.5),
                Measurement(name="antenna_length", value=1.0),
            ],
        )
        p_b = MorphometricProfile(
            specimen_id="B",
            measurements=[
                Measurement(name="head_length", value=3.0),
                Measurement(name="leg_length", value=4.0),
            ],
        )
        result = compare_profiles(p_a, p_b)
        assert "antenna_length" in result["unique_to_a"]
        assert "leg_length" in result["unique_to_b"]
        assert result["shared_measurements"] == 1

    def test_empty_profiles(self):
        p_a = MorphometricProfile(specimen_id="A")
        p_b = MorphometricProfile(specimen_id="B")
        result = compare_profiles(p_a, p_b)
        assert result["shared_measurements"] == 0
        assert result["comparisons"] == {}


# ---------------------------------------------------------------------------
# summary_statistics
# ---------------------------------------------------------------------------


class TestSummaryStatistics:
    def test_basic_summary(self):
        profiles = _make_profiles(n=10)
        stats = summary_statistics(profiles)
        assert "head_length" in stats
        assert "head_width" in stats
        assert "thorax_length" in stats
        assert "femur_length" in stats

    def test_stat_keys(self):
        profiles = _make_profiles(n=10)
        stats = summary_statistics(profiles)
        for name, s in stats.items():
            assert "mean" in s
            assert "median" in s
            assert "std" in s
            assert "min" in s
            assert "max" in s
            assert "n" in s
            assert "cv" in s

    def test_n_equals_profiles(self):
        profiles = _make_profiles(n=8)
        stats = summary_statistics(profiles)
        for name, s in stats.items():
            assert s["n"] == 8

    def test_min_less_than_max(self):
        profiles = _make_profiles(n=10)
        stats = summary_statistics(profiles)
        for name, s in stats.items():
            assert s["min"] <= s["max"]

    def test_mean_between_min_max(self):
        profiles = _make_profiles(n=10)
        stats = summary_statistics(profiles)
        for name, s in stats.items():
            assert s["min"] <= s["mean"] <= s["max"]

    def test_cv_non_negative(self):
        profiles = _make_profiles(n=10)
        stats = summary_statistics(profiles)
        for name, s in stats.items():
            assert s["cv"] >= 0

    def test_single_profile(self):
        profiles = [_make_profile()]
        stats = summary_statistics(profiles)
        assert stats["head_length"]["n"] == 1
        assert stats["head_length"]["std"] == 0.0

    def test_empty_profiles(self):
        profiles = [MorphometricProfile(specimen_id="empty")]
        stats = summary_statistics(profiles)
        assert stats == {}
