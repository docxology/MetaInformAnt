"""``read_mgf`` must not lose a spectrum because its block is not closed, and must not invent a retention time.

Regression tests for ``metainformant.metabolomics.io.formats``:

* ``read_mgf`` split the text on ``BEGIN IONS`` and skipped, without a word, every piece that had no ``END IONS``.
  A file cut off inside its last spectrum therefore came back one spectrum short and looked complete; a block
  that was never closed in the middle of a file was dropped the same way. Both are now a ``ValueError`` that names
  the spectrum, and so is an ``END IONS`` that closes nothing.
* The markers were searched for anywhere in the text, so ``BEGIN IONS`` or ``END IONS`` inside a ``TITLE=`` header
  cut a spectrum in two. They are recognised as whole lines now.
* A spectrum without ``RTINSECONDS`` got ``retention_time == 0.0``, a measured-looking value that
  ``filter_spectra(rt_range=...)`` then used; it is NaN now.
* ``MassSpectrum.total_ion_current`` defaulted to ``0.0`` although it is documented as the sum of the intensities,
  so ``filter_spectra(min_tic=...)`` dropped every hand-built spectrum.
* ``write_mgf`` wrote intensities with one decimal, so ``0.0123`` came back as ``0.0``.

Real files in ``tmp_path``; no test doubles (real-implementation policy).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from metainformant.metabolomics.io.formats import MassSpectrum, filter_spectra, read_mgf, write_mgf


def _block(title: str, peaks: str = "100.0 10\n200.0 20", extra: str = "") -> str:
    return f"BEGIN IONS\nTITLE={title}\nPEPMASS=300.0\nRTINSECONDS=12.5\n{extra}{peaks}\nEND IONS\n"


# every default block above is 7 lines long, so whatever follows one starts on line 8
LINE_AFTER_ONE_BLOCK = _block("x").count("\n") + 1


class TestTruncatedFiles:
    @pytest.mark.parametrize(
        "tail",
        [
            "BEGIN IONS\n",
            "BEGIN IONS\nTITLE=cut\n",
            "BEGIN IONS\nTITLE=cut\nPEPMASS=300.0\nRTINSECONDS=3\n",
            "BEGIN IONS\nTITLE=cut\nPEPMASS=300.0\n100.0 10\n200.0 2",
            "BEGIN IONS\nTITLE=cut\nPEPMASS=300.0\n100.0 10\n\n\n",
            "BEGIN IONS",
        ],
    )
    def test_a_final_block_without_end_ions_raises_and_names_the_spectrum(self, tmp_path: Path, tail: str) -> None:
        path = tmp_path / "truncated.mgf"
        path.write_text(_block("one") + _block("two") + tail)

        with pytest.raises(ValueError, match=r"spectrum 3: BEGIN IONS \(line \d+\).*not closed by END IONS") as error:
            read_mgf(path)

        assert str(path) in str(error.value)
        assert "truncated" in str(error.value)

    def test_the_message_carries_the_title_of_the_cut_spectrum_and_the_line_it_started_on(self, tmp_path: Path) -> None:
        path = tmp_path / "truncated.mgf"
        path.write_text(_block("first") + "BEGIN IONS\nTITLE=lactate_3\nPEPMASS=89.02\n90.0 5\n")

        with pytest.raises(ValueError) as error:
            read_mgf(path)

        assert f"spectrum 2: BEGIN IONS (line {LINE_AFTER_ONE_BLOCK}) with TITLE='lactate_3'" in str(error.value)
        assert LINE_AFTER_ONE_BLOCK == 8

    def test_a_truncated_file_returns_nothing_rather_than_the_spectra_before_the_cut(self, tmp_path: Path) -> None:
        path = tmp_path / "truncated.mgf"
        path.write_text(_block("one") + "BEGIN IONS\n100.0 1\n")

        with pytest.raises(ValueError):
            read_mgf(path)

    def test_a_file_that_is_cut_inside_the_only_spectrum_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "only.mgf"
        path.write_text("BEGIN IONS\nPEPMASS=300.0\n100.0 10\n")

        with pytest.raises(ValueError, match=r"spectrum 1: BEGIN IONS \(line 1\)"):
            read_mgf(path)

    def test_a_block_that_is_never_closed_in_the_middle_of_the_file_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "middle.mgf"
        path.write_text(_block("one") + "BEGIN IONS\nTITLE=lost\n100.0 1\n" + _block("three"))

        with pytest.raises(ValueError, match=r"spectrum 2: BEGIN IONS \(line 8\) with TITLE='lost' is not closed") as e:
            read_mgf(path)

        assert "before the next BEGIN IONS (line 11)" in str(e.value)

    def test_an_end_ions_without_a_begin_ions_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "orphan.mgf"
        # the BEGIN IONS line of the second spectrum is damaged, so its peaks and its END IONS stand outside any block
        path.write_text(_block("one") + "BEGIN ION\nTITLE=two\n100.0 1\nEND IONS\n")

        with pytest.raises(ValueError, match=r"line 11: END IONS without a matching BEGIN IONS"):
            read_mgf(path)

    def test_a_doubled_end_ions_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "doubled.mgf"
        path.write_text(_block("one") + "END IONS\n")

        with pytest.raises(ValueError, match=r"line 8: END IONS without a matching BEGIN IONS"):
            read_mgf(path)


class TestWellFormedFilesAreStillRead:
    def test_a_preamble_and_text_between_blocks_are_ignored(self, tmp_path: Path) -> None:
        path = tmp_path / "global.mgf"
        path.write_text(
            "COM=run 7\nCHARGE=2+\nMASS=Monoisotopic\n"
            + _block("one")
            + "this line sits between two spectra\nTITLE=not a header of any spectrum\n"
            + _block("two")
        )

        first, second = read_mgf(path)

        assert first.mz_array.tolist() == [100.0, 200.0] and second.mz_array.tolist() == [100.0, 200.0]
        assert (first.scan_number, second.scan_number) == (0, 1)

    def test_crlf_line_endings_blank_lines_and_a_missing_final_newline_are_fine(self, tmp_path: Path) -> None:
        path = tmp_path / "crlf.mgf"
        body = (_block("one") + "\n\n" + _block("two")).replace("\n", "\r\n")
        path.write_bytes(body.rstrip().encode())  # no newline after the last END IONS

        assert [s.mz_array.size for s in read_mgf(path)] == [2, 2]

    def test_indented_and_padded_marker_lines_are_markers(self, tmp_path: Path) -> None:
        path = tmp_path / "padded.mgf"
        path.write_text("  BEGIN IONS  \n100.0 1\n\tEND IONS \n")

        (spectrum,) = read_mgf(path)

        assert spectrum.mz_array.tolist() == [100.0]

    def test_marker_text_inside_a_title_does_not_cut_a_spectrum_in_two(self, tmp_path: Path) -> None:
        path = tmp_path / "title.mgf"
        path.write_text(
            "BEGIN IONS\nTITLE=lysate BEGIN IONS replicate\nPEPMASS=300.0\n100.0 1\n"
            "# END IONS is only a comment here\n200.0 2\nEND IONS\n"
        )

        (spectrum,) = read_mgf(path)

        assert spectrum.mz_array.tolist() == [100.0, 200.0]

    def test_an_empty_block_is_an_empty_spectrum(self, tmp_path: Path) -> None:
        path = tmp_path / "empty_block.mgf"
        path.write_text("BEGIN IONS\nPEPMASS=300.0\nEND IONS\n")

        (spectrum,) = read_mgf(path)

        assert spectrum.mz_array.size == 0 and spectrum.total_ion_current == 0.0


class TestRetentionTimeIsNeverInvented:
    def test_a_spectrum_without_rtinseconds_has_a_nan_retention_time(self, tmp_path: Path) -> None:
        path = tmp_path / "no_rt.mgf"
        path.write_text("BEGIN IONS\nPEPMASS=300.0\n100.0 10\nEND IONS\n" + _block("with"))

        without, with_rt = read_mgf(path)

        assert math.isnan(without.retention_time)
        assert with_rt.retention_time == 12.5

    def test_filter_spectra_does_not_treat_an_unknown_retention_time_as_zero(self, tmp_path: Path) -> None:
        path = tmp_path / "no_rt.mgf"
        path.write_text("BEGIN IONS\nPEPMASS=300.0\n100.0 10\nEND IONS\n" + _block("with"))
        spectra = read_mgf(path)

        inside_zero_window = filter_spectra(spectra, min_peaks=1, rt_range=(0.0, 1.0))
        around_known_rt = filter_spectra(spectra, min_peaks=1, rt_range=(10.0, 15.0))
        unfiltered = filter_spectra(spectra, min_peaks=1)

        assert inside_zero_window == []  # the unknown RT is not 0.0 s
        assert [s.retention_time for s in around_known_rt] == [12.5]
        assert len(unfiltered) == 2  # without a range nothing depends on the RT

    def test_write_mgf_leaves_out_rtinseconds_for_a_nan_retention_time_and_read_mgf_gives_nan_back(
        self, tmp_path: Path
    ) -> None:
        spectrum = MassSpectrum(3, math.nan, 2, np.array([100.0]), np.array([10.0]), precursor_mz=250.0)
        path = tmp_path / "nan_rt.mgf"

        write_mgf([spectrum], path)
        (reloaded,) = read_mgf(path)

        assert "RTINSECONDS" not in path.read_text()
        assert math.isnan(reloaded.retention_time)
        assert reloaded.scan_number == 3 and reloaded.precursor_mz == 250.0


class TestTotalIonCurrent:
    def test_it_is_computed_from_the_intensities_when_left_out(self) -> None:
        spectrum = MassSpectrum(0, 1.0, 1, np.array([100.0, 200.0, 300.0]), np.array([1.5, 2.5, 6.0]))

        assert spectrum.total_ion_current == 10.0

    def test_a_value_given_explicitly_is_kept_even_when_it_is_zero(self) -> None:
        intensities = np.array([1.0, 2.0])

        assert (
            MassSpectrum(0, 1.0, 1, np.array([1.0, 2.0]), intensities, total_ion_current=0.0).total_ion_current == 0.0
        )
        assert (
            MassSpectrum(0, 1.0, 1, np.array([1.0, 2.0]), intensities, total_ion_current=99.0).total_ion_current == 99.0
        )

    def test_an_empty_spectrum_has_a_total_ion_current_of_zero(self) -> None:
        assert MassSpectrum(0, 1.0, 1, np.array([]), np.array([])).total_ion_current == 0.0

    def test_filter_spectra_min_tic_sees_the_computed_value(self) -> None:
        strong = MassSpectrum(0, 1.0, 1, np.arange(1.0, 7.0), np.full(6, 100.0))
        weak = MassSpectrum(1, 2.0, 1, np.arange(1.0, 7.0), np.full(6, 0.5))

        kept = filter_spectra([strong, weak], min_peaks=1, min_tic=50.0)

        assert kept == [strong]  # a TIC left at its old default of 0.0 would have dropped both

    def test_read_mgf_still_reports_the_sum_of_the_peaks_it_read(self, tmp_path: Path) -> None:
        path = tmp_path / "tic.mgf"
        path.write_text(_block("one", peaks="100.0 10\n200.0 25.5"))

        (spectrum,) = read_mgf(path)

        assert spectrum.total_ion_current == 35.5


class TestWriteMgfKeepsSmallIntensities:
    def test_a_normalised_spectrum_survives_a_write_and_a_read(self, tmp_path: Path) -> None:
        spectrum = MassSpectrum(
            1, 5.0, 2, np.array([100.123456, 200.5, 301.25]), np.array([0.0123, 0.5, 1.0]), precursor_mz=400.0
        )
        path = tmp_path / "normalised.mgf"

        write_mgf([spectrum], path)
        (reloaded,) = read_mgf(path)

        assert reloaded.intensity_array.tolist() == [0.0123, 0.5, 1.0]  # was [0.0, 0.5, 1.0]
        np.testing.assert_allclose(reloaded.mz_array, spectrum.mz_array, rtol=0, atol=5e-7)
        assert reloaded.total_ion_current == pytest.approx(1.5123)

    def test_large_counts_round_trip_exactly(self, tmp_path: Path) -> None:
        spectrum = MassSpectrum(0, 1.0, 2, np.array([150.0]), np.array([123456789.0]))
        path = tmp_path / "counts.mgf"

        write_mgf([spectrum], path)

        assert read_mgf(path)[0].intensity_array.tolist() == [123456789.0]
