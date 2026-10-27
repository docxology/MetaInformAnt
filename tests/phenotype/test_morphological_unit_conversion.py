"""Regression tests: ``Measurement.convert`` is exact for every supported unit pair and rejects other units.

Defect: ``convert`` carried the note "Simple conversion (placeholder for comprehensive unit lib)" and multiplied
through inexact float constants (``0.001``, ``0.01``, ``0.000001``) into metres and back out again. Most
conversions therefore came out one or two units in the last place off: 1 mm was 1000.0000000000001 um and
1 um was 9.999999999999999e-05 cm. The docstring promised only mm <-> cm <-> m although micrometres were also
converted, and the micrometre's SI symbol (micro sign) was rejected as an unknown unit.

The oracle is ``decimal`` arithmetic: a float is an exact binary fraction, a power of ten is exact in decimal, and
``float(Decimal)`` rounds correctly, so the expected value of every conversion is the exactly computed product
rounded once.
"""

from __future__ import annotations

import itertools
import math
from decimal import Decimal

import pytest

from metainformant.phenotype.morphological.measurement import Measurement

# Decimal exponent of one unit in metres: 1 um = 1e-6 m, 1 mm = 1e-3 m, 1 cm = 1e-2 m, 1 m = 1e0 m.
EXPONENT_IN_METRES = {"um": -6, "mm": -3, "cm": -2, "m": 0}
UNITS = list(EXPONENT_IN_METRES)
VALUES = [1.0, 0.1, 0.3, 2.5, 7.0, 123.456, 1e-3, 3.3, 1e15, 1e-9, 0.7, 25.0, 1000.0, -4.2, 0.0]
MICRO_SIGN = "µm"
GREEK_MU = "μm"


def expected(value: float, source: str, target: str) -> float:
    """value * 10**(exp_source - exp_target), computed exactly in decimal and rounded once to a float."""
    exponent = EXPONENT_IN_METRES[source] - EXPONENT_IN_METRES[target]
    return float(Decimal(value) * Decimal(10) ** exponent)


@pytest.mark.parametrize(("source", "target"), list(itertools.permutations(UNITS, 2)))
def test_every_supported_pair_is_converted_exactly(source: str, target: str) -> None:
    for value in VALUES:
        converted = Measurement("length", value, source).convert(target)

        assert converted.value == expected(value, source, target), (value, source, target)
        assert converted.unit == target


def test_the_conversions_that_used_to_leave_rounding_error_are_exact() -> None:
    assert Measurement("l", 1.0, "mm").convert("um").value == 1000.0  # was 1000.0000000000001
    assert Measurement("l", 1.0, "um").convert("cm").value == 1e-4  # was 9.999999999999999e-05
    assert Measurement("l", 0.1, "mm").convert("um").value == 100.0  # was 100.00000000000001
    assert Measurement("l", 7.0, "mm").convert("um").value == 7000.0  # was 7000.000000000001
    assert Measurement("l", 2.5, "um").convert("mm").value == 0.0025  # was 0.0024999999999999996
    assert Measurement("l", 25.0, "mm").convert("cm").value == 2.5
    assert Measurement("l", 1000.0, "mm").convert("m").value == 1.0


def test_the_uncertainty_is_converted_exactly_too() -> None:
    converted = Measurement("l", 10.0, "mm", uncertainty=0.1).convert("um")

    assert converted.value == 10_000.0
    assert converted.uncertainty == 100.0  # 0.1 * 1000, was 100.00000000000001
    for source, target in itertools.permutations(UNITS, 2):
        assert Measurement("l", 3.0, source, uncertainty=0.3).convert(target).uncertainty == expected(
            0.3, source, target
        )


def test_a_missing_uncertainty_stays_missing() -> None:
    assert Measurement("l", 10.0, "mm").convert("cm").uncertainty is None


@pytest.mark.parametrize("spelling", [MICRO_SIGN, GREEK_MU])
def test_the_micrometre_symbol_is_accepted_as_a_source_and_as_a_target(spelling: str) -> None:
    from_micro = Measurement("l", 500.0, spelling).convert("mm")
    to_micro = Measurement("l", 0.5, "mm").convert(spelling)

    assert from_micro.value == 0.5
    assert to_micro.value == 500.0
    assert to_micro.unit == spelling  # the unit label is the one that was asked for


@pytest.mark.parametrize(("source", "target"), [("um", MICRO_SIGN), (MICRO_SIGN, GREEK_MU), (GREEK_MU, "um")])
def test_spellings_of_the_same_unit_leave_the_value_alone(source: str, target: str) -> None:
    converted = Measurement("l", 12.5, source, uncertainty=0.25).convert(target)

    assert (converted.value, converted.uncertainty, converted.unit) == (12.5, 0.25, target)


@pytest.mark.parametrize("target", ["inches", "in", "MM", "Mm", "CM", "km", "nm", "", "mm ", "meter", "feet"])
def test_an_unsupported_target_unit_raises(target: str) -> None:
    with pytest.raises(ValueError, match="Unsupported unit conversion: mm to") as excinfo:
        Measurement("l", 5.0, "mm").convert(target)

    assert "supported units: um, mm, cm, m" in str(excinfo.value)


@pytest.mark.parametrize("source", ["furlong", "MM", "µ", ""])
def test_an_unsupported_source_unit_raises(source: str) -> None:
    with pytest.raises(ValueError, match="Unsupported unit conversion"):
        Measurement("l", 5.0, source).convert("mm")


def test_two_different_unsupported_units_raise() -> None:
    with pytest.raises(ValueError, match="Unsupported unit conversion: furlong to chain"):
        Measurement("l", 5.0, "furlong").convert("chain")


def test_converting_to_the_own_unit_returns_the_same_object() -> None:
    measurement = Measurement("l", 5.0, "mm")
    odd = Measurement("l", 5.0, "furlong")

    assert measurement.convert("mm") is measurement
    assert odd.convert("furlong") is odd  # nothing to convert, so nothing to reject


def test_name_and_description_are_preserved() -> None:
    converted = Measurement("head_length", 2.5, "mm", uncertainty=0.1, description="HL").convert("cm")

    assert (converted.name, converted.description) == ("head_length", "HL")


def test_the_original_measurement_is_not_modified() -> None:
    original = Measurement("l", 2.5, "mm", uncertainty=0.1)

    original.convert("um")

    assert (original.value, original.unit, original.uncertainty) == (2.5, "mm", 0.1)


def test_a_round_trip_returns_the_original_within_one_rounding() -> None:
    for source, target in itertools.permutations(UNITS, 2):
        for value in VALUES:
            there = Measurement("l", value, source).convert(target)
            back = there.convert(source)
            assert math.isclose(back.value, value, rel_tol=4e-16, abs_tol=0.0), (value, source, target)


def test_non_finite_values_pass_through() -> None:
    assert Measurement("l", math.inf, "mm").convert("um").value == math.inf
    assert Measurement("l", -math.inf, "m").convert("mm").value == -math.inf
    assert math.isnan(Measurement("l", math.nan, "mm").convert("cm").value)
    assert Measurement("l", 1.0, "mm", uncertainty=math.inf).convert("cm").uncertainty == math.inf


def test_a_result_beyond_the_float_range_overflows_to_infinity_of_the_same_sign() -> None:
    assert Measurement("l", 1e308, "m").convert("um").value == math.inf
    assert Measurement("l", -1e308, "m").convert("um").value == -math.inf


def test_a_result_below_the_float_range_underflows_to_zero() -> None:
    assert Measurement("l", 5e-324, "um").convert("m").value == 0.0


def test_integer_values_are_accepted() -> None:
    converted = Measurement("l", 3, "mm").convert("um")  # type: ignore[arg-type]

    assert converted.value == 3000.0


def test_the_docstring_lists_the_supported_units_and_is_no_longer_a_placeholder_note() -> None:
    doc = Measurement.convert.__doc__ or ""

    assert "placeholder" not in doc.lower()
    for unit in ("um", "mm", "cm", "``m``"):
        assert unit in doc
    assert "ValueError" in doc
