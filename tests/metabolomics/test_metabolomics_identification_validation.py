"""``ion_mode`` and ``ppm_tolerance`` of the identification functions are validated.

Regression tests for ``metainformant.metabolomics.analysis.identification``:

* ``identify_with_adducts(..., ion_mode="POSITIVE")`` fell through the ``== "positive"`` test and silently used the
  *negative* adduct set (as did ``"bogus"``), so a protonated ion was reported as unidentified;
* ``ppm_tolerance`` was never checked: ``0`` scored an exact match as 0.0 (0/0), a negative or NaN tolerance matched
  nothing without saying why.

Real computation on small databases; no test doubles.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from metainformant.metabolomics.analysis.identification import (
    COMMON_ADDUCTS,
    identify_metabolites,
    identify_with_adducts,
)

GLUCOSE_MASS = 180.0634
DATABASE = {"glucose": GLUCOSE_MASS}
PROTONATED = np.array([GLUCOSE_MASS + COMMON_ADDUCTS["[M+H]+"]])
DEPROTONATED = np.array([GLUCOSE_MASS + COMMON_ADDUCTS["[M-H]-"]])


class TestIonMode:
    @pytest.mark.parametrize("mode", ["positive", "POSITIVE", "Positive", "  positive  "])
    def test_positive_mode_in_any_case_selects_the_positive_adducts(self, mode: str) -> None:
        reference = identify_with_adducts(PROTONATED, DATABASE, ion_mode="positive")
        matches = identify_with_adducts(PROTONATED, DATABASE, ion_mode=mode)

        assert matches == reference
        assert any(m.matched_name == "glucose" and m.adduct_type == "[M+H]+" for m in matches[0])
        assert all(m.adduct_type.endswith("]+") for m in matches[0])

    @pytest.mark.parametrize("mode", ["negative", "NEGATIVE", "Negative", "\tnegative\n"])
    def test_negative_mode_in_any_case_selects_the_negative_adducts(self, mode: str) -> None:
        matches = identify_with_adducts(DEPROTONATED, DATABASE, ion_mode=mode)

        assert any(m.matched_name == "glucose" and m.adduct_type == "[M-H]-" for m in matches[0])
        assert all(m.adduct_type.endswith("]-") for m in matches[0])

    @pytest.mark.parametrize("bad", ["", "both", "pos", "neg", "neutral", "positve", None, 1])
    def test_an_unknown_ion_mode_raises_instead_of_falling_back_to_the_negative_set(self, bad: Any) -> None:
        with pytest.raises(ValueError, match=r"ion_mode must be 'positive' or 'negative'"):
            identify_with_adducts(PROTONATED, DATABASE, ion_mode=bad)

    def test_the_ion_mode_is_checked_even_when_explicit_adducts_make_it_unused(self) -> None:
        with pytest.raises(ValueError, match=r"ion_mode"):
            identify_with_adducts(PROTONATED, DATABASE, adducts={"[M+H]+": 1.007276}, ion_mode="typo")

    def test_explicit_adducts_are_used_as_given(self) -> None:
        matches = identify_with_adducts(PROTONATED, DATABASE, adducts={"[M+H]+": COMMON_ADDUCTS["[M+H]+"]})

        assert [m.adduct_type for m in matches[0]] == ["[M+H]+"]


BAD_TOLERANCES = [0, 0.0, -1.0, -10, float("nan"), float("inf"), float("-inf"), "ten", None]


class TestPpmTolerance:
    @pytest.mark.parametrize("bad", BAD_TOLERANCES)
    def test_identify_metabolites_rejects_a_meaningless_tolerance(self, bad: Any) -> None:
        with pytest.raises(ValueError, match=r"ppm_tolerance must be a positive finite number"):
            identify_metabolites(np.array([GLUCOSE_MASS]), DATABASE, ppm_tolerance=bad)

    @pytest.mark.parametrize("bad", BAD_TOLERANCES)
    def test_identify_with_adducts_rejects_a_meaningless_tolerance(self, bad: Any) -> None:
        with pytest.raises(ValueError, match=r"ppm_tolerance must be a positive finite number"):
            identify_with_adducts(PROTONATED, DATABASE, ppm_tolerance=bad)

    @pytest.mark.parametrize("tolerance", [5, 5.0, np.float64(5.0), np.int64(5)])
    def test_any_positive_real_number_type_is_accepted(self, tolerance: Any) -> None:
        by_mass = identify_metabolites(np.array([GLUCOSE_MASS]), DATABASE, ppm_tolerance=tolerance)
        by_adduct = identify_with_adducts(PROTONATED, DATABASE, ppm_tolerance=tolerance)

        assert by_mass[0][0].matched_name == "glucose"
        assert by_adduct[0][0].matched_name == "glucose"

    # ppm error is measured relative to the observed m/z, so this observed value is exactly 5 ppm from the database
    FIVE_PPM_OFF = np.array([GLUCOSE_MASS / (1 - 5e-6)])

    def test_a_valid_tolerance_keeps_the_documented_linear_score(self) -> None:
        # 5 ppm error against a 10 ppm window scores 1 - 5/10
        (matches,) = identify_metabolites(self.FIVE_PPM_OFF, DATABASE, ppm_tolerance=10.0)

        assert matches[0].delta_ppm == pytest.approx(5.0, rel=1e-9)
        assert matches[0].score == pytest.approx(0.5, rel=1e-9)

    def test_a_tolerance_just_outside_the_window_still_excludes_the_match(self) -> None:
        (matches,) = identify_metabolites(self.FIVE_PPM_OFF, DATABASE, ppm_tolerance=4.0)

        assert matches == []
