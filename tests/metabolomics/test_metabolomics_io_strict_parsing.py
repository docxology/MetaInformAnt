"""``read_csv`` and ``read_mgf`` must neither invent nor silently drop data.

Regression tests for defects in ``metainformant.metabolomics.io.formats``:

* ``read_csv`` turned every cell it could not parse (``NA``, a blank, a decimal comma) into ``0.0``, so a missing
  measurement became a measured zero;
* ``read_csv`` never filled ``MetabolomicsDataset.metadata``, so annotation columns (m/z, retention time, ...) were
  thrown away;
* ``read_mgf`` skipped every peak line that does not start with a digit (``.5 99``, ``+200.0 20``) and an empty
  ``PEPMASS=`` line raised a bare ``IndexError``.

Real files in ``tmp_path``; no test doubles (real-implementation policy).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from metainformant.metabolomics.analysis.identification import missing_value_imputation
from metainformant.metabolomics.io.formats import MetabolomicsDataset, read_csv, read_mgf, write_csv


class TestCsvMissingValues:
    """A missing cell stays missing (NaN); a cell that is not a number and not a missing marker is an error."""

    @pytest.mark.parametrize("token", ["NA", "na", "N/A", "#N/A", "NaN", "nan", "null", "NULL", "None", ""])
    def test_missing_markers_are_nan_not_zero(self, tmp_path: Path, token: str) -> None:
        path = tmp_path / "missing.csv"
        path.write_text(f"met,S1,S2\nA,1.5,{token}\nB,{token},4\n")

        data = read_csv(path).intensities

        assert np.isnan(data[0, 1])
        assert np.isnan(data[1, 0])
        assert data[0, 0] == 1.5
        assert data[1, 1] == 4.0

    def test_a_measured_zero_stays_zero_and_differs_from_a_missing_value(self, tmp_path: Path) -> None:
        path = tmp_path / "zero.csv"
        path.write_text("met,S1,S2,S3\nA,0,0.0,NA\n")

        row = read_csv(path).intensities[0]

        assert row[0] == 0.0 and row[1] == 0.0
        assert not np.isnan(row[0]) and not np.isnan(row[1])
        assert np.isnan(row[2])

    def test_text_that_is_not_a_number_raises_and_says_where(self, tmp_path: Path) -> None:
        path = tmp_path / "text.csv"
        path.write_text("met,S1,S2\nglucose,1,2\nfructose,3,oops\n")

        with pytest.raises(ValueError, match=r"line 3, metabolite 'fructose', sample 'S2'.*'oops'"):
            read_csv(path)

    def test_a_decimal_comma_is_reported_not_read_as_zero(self, tmp_path: Path) -> None:
        path = tmp_path / "european.csv"
        path.write_text("met;S1;S2\nA;1,5;2\n")

        with pytest.raises(ValueError, match=r"'1,5'.*decimal comma"):
            read_csv(path, delimiter=";")

    def test_a_decimal_comma_in_a_comma_delimited_file_is_a_width_error_not_a_misaligned_row(
        self, tmp_path: Path
    ) -> None:
        # "1,5" splits into two cells, so the row is one cell too wide: it must not be aligned to the header
        path = tmp_path / "split.csv"
        path.write_text("met,S1,S2\nA,1,5,2\n")

        with pytest.raises(ValueError, match=r"line 2 has 4 cells but the header has 3"):
            read_csv(path)

    def test_a_short_row_raises_instead_of_producing_a_ragged_matrix(self, tmp_path: Path) -> None:
        path = tmp_path / "short.csv"
        path.write_text("met,S1,S2\nA,1\nB,3,4\n")

        with pytest.raises(ValueError, match=r"line 2 has 2 cells but the header has 3"):
            read_csv(path)

    def test_nan_survives_a_write_read_round_trip(self, tmp_path: Path) -> None:
        dataset = MetabolomicsDataset(
            metabolites=["A", "B"],
            samples=["S1", "S2"],
            intensities=np.array([[1.0, np.nan], [np.nan, 0.0]]),
        )
        path = tmp_path / "round_trip.csv"

        write_csv(dataset, path)
        loaded = read_csv(path)

        np.testing.assert_array_equal(loaded.intensities, dataset.intensities)  # NaN positions compare equal

    def test_missing_values_flow_into_the_documented_imputation_step(self, tmp_path: Path) -> None:
        path = tmp_path / "pipeline.csv"
        path.write_text("met,S1,S2,S3\nA,10,NA,30\nB,4,6,8\n")

        imputed = missing_value_imputation(read_csv(path).intensities, method="min_half")

        # the gap in row A is filled with half the smallest observed value of that row; the other cells are untouched
        np.testing.assert_allclose(imputed, [[10.0, 5.0, 30.0], [4.0, 6.0, 8.0]])


class TestCsvMetadata:
    """Annotation columns before the intensity columns are returned as per-metabolite metadata."""

    def test_annotation_columns_fill_metadata(self, tmp_path: Path) -> None:
        path = tmp_path / "annotated.csv"
        path.write_text(
            "metabolite,mz,rt,formula,S1,S2\n"
            "glucose,180.0634,1.5,C6H12O6,10,20\n"
            "lactate,89.0244,NA,C3H6O3,30,40\n"
        )

        dataset = read_csv(path, data_start_col=4)

        assert dataset.metabolites == ["glucose", "lactate"]
        assert dataset.samples == ["S1", "S2"]
        np.testing.assert_array_equal(dataset.intensities, [[10.0, 20.0], [30.0, 40.0]])
        assert dataset.metadata["glucose"] == {"mz": 180.0634, "rt": 1.5, "formula": "C6H12O6"}
        lactate: dict[str, Any] = dataset.metadata["lactate"]
        assert lactate["mz"] == 89.0244
        assert math.isnan(lactate["rt"])  # a missing annotation is NaN, not 0.0
        assert lactate["formula"] == "C3H6O3"

    def test_the_default_layout_has_no_annotation_columns_and_no_metadata(self, tmp_path: Path) -> None:
        path = tmp_path / "plain.csv"
        path.write_text("metabolite,S1,S2\nA,1,2\n")

        assert read_csv(path).metadata == {}

    def test_the_metabolite_column_need_not_be_the_first(self, tmp_path: Path) -> None:
        path = tmp_path / "middle.csv"
        path.write_text("mz,metabolite,S1,S2\n180.0634,glucose,10,20\n")

        dataset = read_csv(path, metabolite_col=1, data_start_col=2)

        assert dataset.metabolites == ["glucose"]
        assert dataset.metadata == {"glucose": {"mz": 180.0634}}

    def test_an_unnamed_annotation_column_gets_a_positional_name(self, tmp_path: Path) -> None:
        path = tmp_path / "unnamed.csv"
        path.write_text("metabolite,,S1\nA,7.5,1\n")

        assert read_csv(path, data_start_col=2).metadata == {"A": {"column_1": 7.5}}

    def test_a_repeated_metabolite_with_annotations_is_rejected_not_overwritten(self, tmp_path: Path) -> None:
        path = tmp_path / "duplicate.csv"
        path.write_text("metabolite,mz,S1\nA,100.0,1\nA,200.0,2\n")

        with pytest.raises(ValueError, match=r"metabolite 'A' occurs more than once"):
            read_csv(path, data_start_col=2)

    def test_repeated_annotation_column_names_are_rejected(self, tmp_path: Path) -> None:
        path = tmp_path / "duplicate_columns.csv"
        path.write_text("metabolite,mz,mz,S1\nA,100.0,101.0,1\n")

        with pytest.raises(ValueError, match=r"duplicate annotation column names"):
            read_csv(path, data_start_col=3)


class TestCsvLayout:
    """Header and line handling that used to misalign samples silently."""

    def test_a_header_with_an_empty_first_cell_keeps_every_sample(self, tmp_path: Path) -> None:
        # R's write.table emits such a header; stripping the whole text used to eat the leading tab and drop S1
        path = tmp_path / "r_style.tsv"
        path.write_text("\tS1\tS2\nA\t1\t2\n")

        dataset = read_csv(path, delimiter="\t")

        assert dataset.samples == ["S1", "S2"]
        np.testing.assert_array_equal(dataset.intensities, [[1.0, 2.0]])

    def test_blank_lines_and_windows_line_endings_are_tolerated(self, tmp_path: Path) -> None:
        path = tmp_path / "crlf.csv"
        path.write_bytes(b"met,S1,S2\r\n\r\nA,1,2\r\n\r\nB,3,4\r\n\r\n")

        dataset = read_csv(path)

        assert dataset.metabolites == ["A", "B"]
        assert dataset.samples == ["S1", "S2"]
        np.testing.assert_array_equal(dataset.intensities, [[1.0, 2.0], [3.0, 4.0]])

    def test_a_utf8_byte_order_mark_does_not_leak_into_an_annotation_name(self, tmp_path: Path) -> None:
        # Excel's "CSV UTF-8" export starts with a BOM; with the metabolite in column 1 the first header cell is an
        # annotation name and must read "mz", not "\ufeffmz"
        path = tmp_path / "bom.csv"
        path.write_bytes(b"\xef\xbb\xbfmz,metabolite,S1\n180.0634,glucose,1\n")

        dataset = read_csv(path, metabolite_col=1, data_start_col=2)

        assert dataset.metadata == {"glucose": {"mz": 180.0634}}

    def test_an_empty_file_gives_an_empty_dataset(self, tmp_path: Path) -> None:
        path = tmp_path / "empty.csv"
        path.write_text("")

        dataset = read_csv(path)

        assert dataset.metabolites == [] and dataset.samples == []
        assert dataset.intensities.size == 0

    @pytest.mark.parametrize("kwargs", [{"metabolite_col": 5}, {"data_start_col": 9}, {"metabolite_col": -1}])
    def test_column_indices_outside_the_header_are_rejected(self, tmp_path: Path, kwargs: dict[str, int]) -> None:
        path = tmp_path / "narrow.csv"
        path.write_text("met,S1\nA,1\n")

        with pytest.raises(ValueError, match=r"do not fit the 2 header columns"):
            read_csv(path, **kwargs)


def _mgf(*blocks: str) -> str:
    return "".join(f"BEGIN IONS\n{block}\nEND IONS\n" for block in blocks)


class TestMgfPeakLines:
    """Every peak line is read or reported; none is skipped silently."""

    def test_peaks_that_do_not_start_with_a_digit_are_read(self, tmp_path: Path) -> None:
        path = tmp_path / "peaks.mgf"
        path.write_text(_mgf("PEPMASS=300.0\nRTINSECONDS=10\n100.0 10\n.5 99\n+200.0 20\n3.0E+02\t30"))

        (spectrum,) = read_mgf(path)

        assert spectrum.mz_array.tolist() == [100.0, 0.5, 200.0, 300.0]
        assert spectrum.intensity_array.tolist() == [10.0, 99.0, 20.0, 30.0]
        assert spectrum.total_ion_current == 159.0

    @pytest.mark.parametrize("comment", ["# note", "; note", "! note", "/ note"])
    def test_comment_lines_are_skipped(self, tmp_path: Path, comment: str) -> None:
        path = tmp_path / "comment.mgf"
        path.write_text(_mgf(f"{comment}\n100.0 1\n{comment}\n200.0 2"))

        (spectrum,) = read_mgf(path)

        assert spectrum.mz_array.tolist() == [100.0, 200.0]

    def test_headers_other_than_the_three_that_are_read_are_ignored(self, tmp_path: Path) -> None:
        path = tmp_path / "headers.mgf"
        path.write_text(_mgf("TITLE=sample=1\nCHARGE=2+\nMSLEVEL=2\nPEPMASS=250.5 99.0\nSCANS=7\n100.0 1"))

        (spectrum,) = read_mgf(path)

        assert spectrum.precursor_mz == 250.5  # the first number of "mass intensity"
        assert spectrum.scan_number == 7
        assert spectrum.mz_array.tolist() == [100.0]

    def test_a_block_that_is_never_closed_is_an_error_not_a_dropped_spectrum(self, tmp_path: Path) -> None:
        # it used to be dropped without a word; see test_metabolomics_mgf_truncation.py for the full contract
        path = tmp_path / "truncated.mgf"
        path.write_text(_mgf("100.0 1") + "BEGIN IONS\n200.0 2\n")

        with pytest.raises(ValueError, match=r"spectrum 2: .*not closed by END IONS"):
            read_mgf(path)

    def test_a_third_peak_column_such_as_a_charge_is_ignored(self, tmp_path: Path) -> None:
        path = tmp_path / "charge.mgf"
        path.write_text(_mgf("100.0 10 2+\n200.0 20 1-"))

        (spectrum,) = read_mgf(path)

        assert spectrum.mz_array.tolist() == [100.0, 200.0]
        assert spectrum.intensity_array.tolist() == [10.0, 20.0]

    @pytest.mark.parametrize("line", ["abc 5", "100.0", "100.0 abc", "nan 5", "100.0 inf", "5,5 3"])
    def test_a_malformed_peak_line_is_an_error_that_names_the_spectrum_and_the_line(
        self, tmp_path: Path, line: str
    ) -> None:
        path = tmp_path / "bad_peak.mgf"
        path.write_text(_mgf(f"PEPMASS=300.0\n100.0 10\n{line}"))

        with pytest.raises(ValueError, match=r"spectrum 1:.*" + line.split()[0].replace(".", r"\.")):
            read_mgf(path)

    def test_the_spectrum_number_counts_blocks_in_file_order(self, tmp_path: Path) -> None:
        path = tmp_path / "second.mgf"
        path.write_text(_mgf("100.0 1", "100.0 1", "nope 1"))

        with pytest.raises(ValueError, match=r"spectrum 3:"):
            read_mgf(path)


class TestMgfHeaderLines:
    """An unreadable header value is a clear ValueError, never an IndexError."""

    @pytest.mark.parametrize("line", ["PEPMASS=", "PEPMASS=   "])
    def test_an_empty_pepmass_is_a_clear_error(self, tmp_path: Path, line: str) -> None:
        path = tmp_path / "empty_pepmass.mgf"
        path.write_text(_mgf(f"{line}\n100.0 10"))

        with pytest.raises(ValueError, match=r"spectrum 1: PEPMASS has no value"):
            read_mgf(path)

    @pytest.mark.parametrize(
        ("line", "message"),
        [
            ("PEPMASS=abc", r"PEPMASS must be a number, got 'abc'"),
            ("RTINSECONDS=", r"RTINSECONDS must be a number, got ''"),
            ("RTINSECONDS=10-20", r"RTINSECONDS must be a number, got '10-20'"),
            ("SCANS=12-15", r"SCANS must be an integer, got '12-15'"),
            ("SCANS=", r"SCANS must be an integer, got ''"),
        ],
    )
    def test_a_non_numeric_header_value_names_the_key_and_the_value(
        self, tmp_path: Path, line: str, message: str
    ) -> None:
        path = tmp_path / "bad_header.mgf"
        path.write_text(_mgf(f"{line}\n100.0 10"))

        with pytest.raises(ValueError, match=message):
            read_mgf(path)

    def test_a_spectrum_without_a_pepmass_line_still_has_no_precursor(self, tmp_path: Path) -> None:
        path = tmp_path / "no_pepmass.mgf"
        path.write_text(_mgf("RTINSECONDS=5\n100.0 10"))

        (spectrum,) = read_mgf(path)

        assert spectrum.precursor_mz is None
        assert spectrum.retention_time == 5.0
