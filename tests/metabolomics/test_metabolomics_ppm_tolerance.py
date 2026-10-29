"""``extract_chromatogram`` and ``prepare_spectrum_plot`` validate ``ppm_tolerance`` (and the arguments beside it).

Regression tests for ``metainformant.metabolomics.io.formats.extract_chromatogram`` and
``metainformant.metabolomics.visualization.plots.prepare_spectrum_plot``. Both accepted any ``ppm_tolerance``:

* a negative or NaN tolerance matched nothing, so ``extract_chromatogram`` returned an all-zero chromatogram that
  looks like an ion that is simply absent, and ``prepare_spectrum_plot`` returned an unannotated spectrum;
* an infinite tolerance matched every peak;
* ``target_mz`` was not checked either (NaN or a negative value also gave an all-zero chromatogram).

A tolerance of exactly zero is meaningful (only a peak at exactly the target m/z is inside the window) and must
keep working. Real computation on small spectra; no test doubles (real-implementation policy).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from metainformant.metabolomics.io.formats import MassSpectrum, extract_chromatogram
from metainformant.metabolomics.visualization.plots import prepare_spectrum_plot

BAD_TOLERANCES: list[Any] = [-1.0, -1e-12, float("nan"), float("inf"), float("-inf"), "ten", "5", b"5", None, [5.0]]
MESSAGE = r"ppm_tolerance must be a finite number >= 0"

TARGET = 500.0
# one ppm of 500 is 0.0005, so these peaks are exactly on, 0.2 ppm off, 4 ppm off and 40 ppm off the target
PEAK_MZ = np.array([500.0, 500.0001, 500.002, 500.02])
PEAK_INTENSITY = np.array([10.0, 20.0, 30.0, 40.0])


def _scans() -> list[MassSpectrum]:
    # four MS1 scans, each with the peak at the target m/z shifted by a different amount
    return [
        MassSpectrum(i, 10.0 * i, 1, np.array([400.0, mz]), np.array([5.0, 100.0 + i])) for i, mz in enumerate(PEAK_MZ)
    ]


class TestExtractChromatogramTolerance:
    @pytest.mark.parametrize("bad", BAD_TOLERANCES)
    def test_a_meaningless_tolerance_raises_instead_of_returning_zeros(self, bad: Any) -> None:
        with pytest.raises(ValueError, match=MESSAGE):
            extract_chromatogram(_scans(), TARGET, ppm_tolerance=bad)

    def test_a_negative_tolerance_no_longer_looks_like_an_absent_ion(self) -> None:
        _, intensities = extract_chromatogram(_scans(), TARGET, ppm_tolerance=1.0)
        assert intensities.max() > 0  # the ion is there ...

        with pytest.raises(ValueError, match=MESSAGE):  # ... and a negative window does not claim otherwise
            extract_chromatogram(_scans(), TARGET, ppm_tolerance=-1.0)

    def test_zero_keeps_only_peaks_at_exactly_the_target(self) -> None:
        rts, intensities = extract_chromatogram(_scans(), TARGET, ppm_tolerance=0)

        assert rts.tolist() == [0.0, 10.0, 20.0, 30.0]
        assert intensities.tolist() == [100.0, 0.0, 0.0, 0.0]  # scan 0 holds the exact peak, 1-3 are off by >= 0.2 ppm

    def test_zero_with_no_exact_peak_is_an_honest_zero_chromatogram(self) -> None:
        scans = [MassSpectrum(0, 1.0, 1, np.array([500.0000001]), np.array([50.0]))]

        _, intensities = extract_chromatogram(scans, TARGET, ppm_tolerance=0.0)

        assert intensities.tolist() == [0.0]

    @pytest.mark.parametrize(
        ("tolerance", "expected"),
        [
            (0.0, [100.0, 0.0, 0.0, 0.0]),
            (1.0, [100.0, 101.0, 0.0, 0.0]),
            (5.0, [100.0, 101.0, 102.0, 0.0]),
            (50.0, None),
        ],
    )
    def test_the_window_grows_with_the_tolerance(self, tolerance: float, expected: list[float] | None) -> None:
        _, intensities = extract_chromatogram(_scans(), TARGET, ppm_tolerance=tolerance)

        assert intensities.tolist() == (expected if expected is not None else [100.0, 101.0, 102.0, 103.0])

    @pytest.mark.parametrize("tolerance", [5, 5.0, np.float64(5.0), np.int64(5), np.float32(5.0)])
    def test_any_real_number_type_is_accepted(self, tolerance: Any) -> None:
        _, intensities = extract_chromatogram(_scans(), TARGET, ppm_tolerance=tolerance)

        assert intensities.tolist() == [100.0, 101.0, 102.0, 0.0]

    def test_the_default_tolerance_is_ten_ppm(self) -> None:
        by_default = extract_chromatogram(_scans(), TARGET)[1]
        explicit = extract_chromatogram(_scans(), TARGET, ppm_tolerance=10.0)[1]

        np.testing.assert_array_equal(by_default, explicit)

    @pytest.mark.parametrize("bad", [0, -500.0, float("nan"), float("inf"), "500", None])
    def test_the_target_must_be_a_positive_finite_number(self, bad: Any) -> None:
        with pytest.raises(ValueError, match=r"target_mz must be a positive finite number"):
            extract_chromatogram(_scans(), bad)

    def test_the_tolerance_is_checked_even_when_there_is_no_ms1_scan(self) -> None:
        ms2_only = [MassSpectrum(0, 1.0, 2, np.array([500.0]), np.array([1.0]))]

        with pytest.raises(ValueError, match=MESSAGE):
            extract_chromatogram(ms2_only, TARGET, ppm_tolerance=float("nan"))


class TestPrepareSpectrumPlotTolerance:
    @pytest.mark.parametrize("bad", BAD_TOLERANCES)
    def test_a_meaningless_tolerance_raises(self, bad: Any) -> None:
        with pytest.raises(ValueError, match=MESSAGE):
            prepare_spectrum_plot(PEAK_MZ, PEAK_INTENSITY, annotation_db={500.0: "x"}, ppm_tolerance=bad)

    @pytest.mark.parametrize("bad", BAD_TOLERANCES)
    def test_it_is_checked_even_when_no_annotation_database_uses_it(self, bad: Any) -> None:
        with pytest.raises(ValueError, match=MESSAGE):
            prepare_spectrum_plot(PEAK_MZ, PEAK_INTENSITY, ppm_tolerance=bad)

    @pytest.mark.parametrize("bad", BAD_TOLERANCES)
    def test_it_is_checked_even_for_an_empty_spectrum(self, bad: Any) -> None:
        with pytest.raises(ValueError, match=MESSAGE):
            prepare_spectrum_plot(np.array([]), np.array([]), ppm_tolerance=bad)

    def test_zero_annotates_only_peaks_at_exactly_a_database_mass(self) -> None:
        result = prepare_spectrum_plot(
            PEAK_MZ, PEAK_INTENSITY, annotation_db={500.0: "exact", 500.002: "four ppm"}, ppm_tolerance=0.0
        )

        assert result.annotations == {500.0: "exact", 500.002: "four ppm"}  # each equals a database mass exactly

    def test_zero_does_not_annotate_a_peak_that_is_off_by_a_fraction_of_a_ppm(self) -> None:
        result = prepare_spectrum_plot(PEAK_MZ, PEAK_INTENSITY, annotation_db={500.0: "exact"}, ppm_tolerance=0.0)

        assert result.annotations == {500.0: "exact"}  # 500.0001 is 0.2 ppm away

    @pytest.mark.parametrize(
        ("tolerance", "annotated"),
        [
            (0.0, [500.0]),
            (0.5, [500.0, 500.0001]),
            (5.0, [500.0, 500.0001, 500.002]),
            (100.0, [500.0, 500.0001, 500.002, 500.02]),
        ],
    )
    def test_the_annotation_window_grows_with_the_tolerance(self, tolerance: float, annotated: list[float]) -> None:
        result = prepare_spectrum_plot(PEAK_MZ, PEAK_INTENSITY, annotation_db={500.0: "x"}, ppm_tolerance=tolerance)

        assert sorted(result.annotations) == annotated
        assert set(result.annotations.values()) == {"x"}

    def test_the_default_tolerance_is_ten_ppm(self) -> None:
        by_default = prepare_spectrum_plot(PEAK_MZ, PEAK_INTENSITY, annotation_db={500.0: "x"})
        explicit = prepare_spectrum_plot(PEAK_MZ, PEAK_INTENSITY, annotation_db={500.0: "x"}, ppm_tolerance=10.0)

        assert by_default.annotations == explicit.annotations
        assert sorted(by_default.annotations) == [500.0, 500.0001, 500.002]

    def test_a_peak_at_zero_mass_is_not_a_division_by_zero(self) -> None:
        # pytest turns every RuntimeWarning into an error, so a 0/0 here would fail the test
        result = prepare_spectrum_plot(
            np.array([0.0, 100.0]), np.array([5.0, 10.0]), annotation_db={0.0: "origin", 100.0: "x"}
        )

        assert result.annotations == {0.0: "origin", 100.0: "x"}
