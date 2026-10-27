"""Round trips and column layout of the peak-file readers and writers (ChIP-seq and ATAC-seq).

``save_*_peaks`` used to write the summit in the seventh column and the signal value in the eighth, while
``load_*_peaks`` read the signal value from the seventh and the summit offset from the tenth, so a saved narrowPeak
file loaded back as an empty list (and a broadPeak file loaded back with the summit as the signal value and the
signal as the p-value). Both sides now use the ENCODE layouts:

    narrowPeak: chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak (offset from start)
    broadPeak:  the first nine of those
    BED:        the first six

The pValue and qValue columns hold ``-log10`` values, as the specification and MACS2 write them, while the peaks
carry probabilities: a file column of 2 is a q-value of 0.01. The conversion is not bit-exact in either direction, so
the round trips compare probabilities with ``pytest.approx``.

Every test builds real peaks, writes real files in ``tmp_path`` and reads them back. REAL IMPLEMENTATION.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Callable, List

import pytest

from metainformant.epigenome.assays.atacseq import ATACPeak, load_atac_peaks, save_atac_peaks
from metainformant.epigenome.assays.chipseq import ChIPPeak, load_chip_peaks, save_chip_peaks


def _chip_peaks() -> List[ChIPPeak]:
    return [
        ChIPPeak("chr1", 100, 500, summit=300, score=50.0, strand="+", signal_value=12.5, p_value=1e-5, q_value=1e-3),
        # summit at the very start (offset 0)
        ChIPPeak(
            "chr2", 1000, 1400, summit=1000, score=80.0, strand="-", signal_value=7.25, p_value=0.0123, q_value=0.05
        ),
        # summit at the very end, fractional score, no p-value or q-value
        ChIPPeak("chrX", 5, 25, summit=25, score=1234.5, strand=".", signal_value=0.0, p_value=None, q_value=None),
        # everything at its default
        ChIPPeak("chr3", 10, 20),
    ]


def _atac_peaks() -> List[ATACPeak]:
    return [
        ATACPeak("chr1", 100, 500, score=50.0, strand="+", signal_value=12.5, p_value=1e-5, q_value=1e-3, summit=300),
        ATACPeak(
            "chr2", 1000, 1400, score=80.0, strand="-", signal_value=7.25, p_value=0.0123, q_value=0.05, summit=1000
        ),
        ATACPeak("chrX", 5, 25, score=1234.5, strand=".", signal_value=0.0, p_value=None, q_value=None, summit=25),
        ATACPeak("chr3", 10, 20),
    ]


def _data_rows(path: Path) -> List[List[str]]:
    """Tab-split data lines of a peak file (comment lines dropped)."""
    return [line.split("\t") for line in path.read_text().splitlines() if line and not line.startswith("#")]


def _assert_same_peak(back: Any, original: Any, *, keep_summit: bool = True) -> None:
    """Every field equal; the p-value and q-value only up to the rounding of the -log10 round trip."""
    got, want = back.to_dict(), original.to_dict()
    if not keep_summit:
        want["summit"] = None  # a broadPeak line has no summit column
    for key in ("p_value", "q_value"):
        got_value, want_value = got.pop(key), want.pop(key)
        if want_value is None:
            assert got_value is None, key
        else:
            assert got_value == pytest.approx(want_value, rel=1e-12), key
    assert got == want


# Each assay module exposes the same three functions under different names
CHIP = (_chip_peaks, save_chip_peaks, load_chip_peaks)
ATAC = (_atac_peaks, save_atac_peaks, load_atac_peaks)

# A MACS2-style narrowPeak line: -log10 p-value 12.3 and -log10 q-value 9.1; the summit is 200 bases from the start
MACS2_LINE = "chr1\t100\t500\tpeak1\t50\t.\t7.25\t12.3\t9.1\t200\n"


class TestNarrowPeakRoundTrip:
    @pytest.mark.parametrize("assay", [CHIP, ATAC], ids=["chipseq", "atacseq"])
    def test_every_field_survives(self, tmp_path: Path, assay: Any) -> None:
        make_peaks, save, load = assay
        peaks = make_peaks()
        path = tmp_path / "peaks.narrowPeak"

        save(peaks, path, format="narrowpeak")
        loaded = load(path, format="narrowpeak")

        assert len(loaded) == len(peaks) == 4
        for original, back in zip(peaks, loaded):
            _assert_same_peak(back, original)

    @pytest.mark.parametrize("assay", [CHIP, ATAC], ids=["chipseq", "atacseq"])
    def test_default_format_is_narrowpeak(self, tmp_path: Path, assay: Any) -> None:
        make_peaks, save, load = assay
        path = tmp_path / "peaks.narrowPeak"

        save(make_peaks(), path)

        loaded = load(path)
        assert len(loaded) == 4
        for original, back in zip(make_peaks(), loaded):
            _assert_same_peak(back, original)

    def test_chipseq_columns_are_the_standard_layout(self, tmp_path: Path) -> None:
        path = tmp_path / "chip.narrowPeak"
        save_chip_peaks(_chip_peaks(), path)

        rows = _data_rows(path)

        assert all(len(row) == 10 for row in rows)
        # chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak
        # (the -log10 columns: p = 1e-05 is written 5 and q = 0.001 is written 3)
        assert rows[0] == ["chr1", "100", "500", "peak_100_500", "50", "+", "12.5", "5", "3", "200"]
        # the last column is the summit as an offset from the start, 0 for a summit at the start
        assert rows[1][9] == "0"
        assert rows[1][:3] == ["chr2", "1000", "1400"]
        assert float(rows[1][7]) == pytest.approx(-math.log10(0.0123), rel=1e-12)
        assert float(rows[1][8]) == pytest.approx(-math.log10(0.05), rel=1e-12)
        # unavailable p-value, q-value and summit are -1; a fractional score is not rounded
        assert rows[2] == ["chrX", "5", "25", "peak_5_25", "1234.5", ".", "0", "-1", "-1", "20"]
        assert rows[3][6:] == ["0", "-1", "-1", "-1"]

    def test_atacseq_columns_are_the_standard_layout(self, tmp_path: Path) -> None:
        path = tmp_path / "atac.narrowPeak"
        save_atac_peaks(_atac_peaks(), path)

        rows = _data_rows(path)

        assert all(len(row) == 10 for row in rows)
        assert rows[0] == ["chr1", "100", "500", "atac_peak_100_500", "50", "+", "12.5", "5", "3", "200"]
        assert rows[2][6:] == ["0", "-1", "-1", "20"]

    def test_to_bed_format_is_one_narrowpeak_line(self) -> None:
        peak = ChIPPeak(
            "chr1", 100, 500, summit=300, score=50.0, strand="+", signal_value=12.5, p_value=1e-5, q_value=1e-3
        )

        assert peak.to_bed_format().split("\t") == [
            "chr1",
            "100",
            "500",
            "peak_100_500",
            "50",
            "+",
            "12.5",
            "5",
            "3",
            "200",
        ]
        assert peak.to_bed_format("broadpeak").split("\t")[6:] == ["12.5", "5", "3"]
        assert peak.to_bed_format("bed").split("\t") == ["chr1", "100", "500", "peak_100_500", "50", "+"]

    @pytest.mark.parametrize("loader", [load_chip_peaks, load_atac_peaks], ids=["chipseq", "atacseq"])
    def test_hand_written_standard_file_is_read_by_column(self, tmp_path: Path, loader: Callable[..., Any]) -> None:
        # A MACS2-style narrowPeak: the summit is column 10, an offset from the start; -1 and . mean "none"
        path = tmp_path / "macs2.narrowPeak"
        path.write_text(
            "track name=peaks\n"
            "chr1\t1000\t1400\tpeak_1\t350\t.\t8.5\t12.3\t9.1\t200\n"
            "chr1\t2000\t2300\tpeak_2\t10\t+\t1.5\t-1\t.\t-1\n"
        )

        first, second = loader(path, format="narrowpeak")

        assert (first.chromosome, first.start, first.end) == ("chr1", 1000, 1400)
        assert first.score == 350.0
        assert first.signal_value == 8.5
        # the columns hold -log10 values; the peak carries probabilities
        assert first.p_value == pytest.approx(10**-12.3, rel=1e-12)
        assert first.q_value == pytest.approx(10**-9.1, rel=1e-12)
        assert first.summit == 1200
        assert second.strand == "+"
        assert (second.p_value, second.q_value, second.summit) == (None, None, None)


class TestSignificanceColumnsAreMinusLog10:
    """Peaks carry probabilities; the files store -log10 values (ENCODE specification, MACS2)."""

    @pytest.mark.parametrize("loader", [load_chip_peaks, load_atac_peaks], ids=["chipseq", "atacseq"])
    def test_a_macs2_line_loads_with_probabilities(self, tmp_path: Path, loader: Callable[..., Any]) -> None:
        path = tmp_path / "macs2.narrowPeak"
        path.write_text(MACS2_LINE)

        (peak,) = loader(path, format="narrowpeak")

        assert (peak.chromosome, peak.start, peak.end, peak.score, peak.signal_value) == ("chr1", 100, 500, 50.0, 7.25)
        assert peak.summit == 300
        assert peak.p_value == pytest.approx(10**-12.3, rel=1e-12)
        assert peak.q_value == pytest.approx(10**-9.1, rel=1e-12)
        # a strongly significant MACS2 peak passes an FDR limit of 0.05 (read verbatim, 9.1 would have failed it)
        assert peak.q_value is not None and peak.q_value <= 0.05

    @pytest.mark.parametrize(
        ("probability", "column"),
        [(1e-5, "5"), (1e-3, "3"), (0.01, "2"), (1.0, "0"), (0.5, repr(-math.log10(0.5)))],
    )
    @pytest.mark.parametrize("peak_class", [ChIPPeak, ATACPeak], ids=["chipseq", "atacseq"])
    def test_a_probability_is_written_as_minus_log10(self, peak_class: Any, probability: float, column: str) -> None:
        peak = peak_class("chr1", 0, 100, p_value=probability, q_value=probability)

        fields = peak.to_bed_format().split("\t")

        assert (fields[7], fields[8]) == (column, column)

    @pytest.mark.parametrize(
        ("column", "probability"),
        [("5", 1e-5), ("3", 1e-3), ("2", 0.01), ("0", 1.0), ("9.1", 10**-9.1), ("0.30103", 10**-0.30103)],
    )
    @pytest.mark.parametrize("loader", [load_chip_peaks, load_atac_peaks], ids=["chipseq", "atacseq"])
    def test_a_column_is_read_as_a_probability(
        self, tmp_path: Path, loader: Callable[..., Any], column: str, probability: float
    ) -> None:
        path = tmp_path / "peaks.narrowPeak"
        path.write_text(f"chr1\t0\t100\tp\t5\t.\t1.0\t{column}\t{column}\t50\n")

        (peak,) = loader(path)

        assert peak.p_value == pytest.approx(probability, rel=1e-12)
        assert peak.q_value == pytest.approx(probability, rel=1e-12)

    @pytest.mark.parametrize(("save", "load"), [(save_chip_peaks, load_chip_peaks), (save_atac_peaks, load_atac_peaks)])
    def test_a_macs2_file_keeps_its_columns_through_a_load_and_a_save(
        self, tmp_path: Path, save: Callable[..., Any], load: Callable[..., Any]
    ) -> None:
        source = tmp_path / "macs2.narrowPeak"
        source.write_text(MACS2_LINE)
        copy = tmp_path / "copy.narrowPeak"

        save(load(source), copy)

        (row,) = _data_rows(copy)
        assert float(row[7]) == pytest.approx(12.3, rel=1e-12)
        assert float(row[8]) == pytest.approx(9.1, rel=1e-12)
        assert row[6] == "7.25"
        assert row[9] == "200"

    @pytest.mark.parametrize("peak_class", [ChIPPeak, ATACPeak], ids=["chipseq", "atacseq"])
    def test_probability_zero_is_written_as_the_largest_finite_column(self, peak_class: Any) -> None:
        peak = peak_class("chr1", 0, 100, p_value=0.0, q_value=0.0)

        fields = peak.to_bed_format().split("\t")

        assert float(fields[7]) == pytest.approx(-math.log10(5e-324), rel=1e-15)  # about 323.3
        assert fields[7] == fields[8]

    @pytest.mark.parametrize(
        ("save", "load", "peak_class"),
        [(save_chip_peaks, load_chip_peaks, ChIPPeak), (save_atac_peaks, load_atac_peaks, ATACPeak)],
        ids=["chipseq", "atacseq"],
    )
    def test_probability_zero_reads_back_as_the_smallest_float(
        self, tmp_path: Path, save: Callable[..., Any], load: Callable[..., Any], peak_class: Any
    ) -> None:
        path = tmp_path / "zero.narrowPeak"
        save([peak_class("chr1", 0, 100, p_value=0.0, q_value=0.0)], path)

        (back,) = load(path)

        assert back.p_value == 5e-324 and back.q_value == 5e-324

    @pytest.mark.parametrize("loader", [load_chip_peaks, load_atac_peaks], ids=["chipseq", "atacseq"])
    def test_a_column_beyond_the_float_range_reads_as_zero(self, tmp_path: Path, loader: Callable[..., Any]) -> None:
        # 10 ** -400 underflows: the strongest peaks of a deep experiment are probability 0.0
        path = tmp_path / "strong.narrowPeak"
        path.write_text("chr1\t0\t100\tp\t5\t.\t1.0\t400\t400\t50\n")

        (peak,) = loader(path)

        assert (peak.p_value, peak.q_value) == (0.0, 0.0)

    @pytest.mark.parametrize("loader", [load_chip_peaks, load_atac_peaks], ids=["chipseq", "atacseq"])
    def test_negative_columns_other_than_minus_one_are_skipped(
        self, tmp_path: Path, loader: Callable[..., Any]
    ) -> None:
        # -2 and nan are not -log10 of a probability; only -1 (and .) mean "not available"
        path = tmp_path / "bad.narrowPeak"
        path.write_text(
            "chr1\t0\t100\tgood1\t5\t.\t1.0\t2\t2\t50\n"
            "chr1\t100\t200\tneg\t5\t.\t1.0\t-2\t2\t50\n"
            "chr1\t200\t300\tnan\t5\t.\t1.0\t2\tnan\t50\n"
            "chr1\t300\t400\tgood2\t5\t.\t1.0\t-1\t.\t50\n"
        )

        peaks = loader(path)

        assert [(p.start, p.p_value, p.q_value) for p in peaks] == [
            (0, pytest.approx(0.01), pytest.approx(0.01)),
            (300, None, None),
        ]

    @pytest.mark.parametrize(("save", "peak_class"), [(save_chip_peaks, ChIPPeak), (save_atac_peaks, ATACPeak)])
    @pytest.mark.parametrize("bad", [5.2, -0.1, 1.5, float("nan")])
    @pytest.mark.parametrize("field", ["p_value", "q_value"])
    def test_a_value_outside_zero_one_is_rejected(
        self, tmp_path: Path, save: Callable[..., Any], peak_class: Any, bad: float, field: str
    ) -> None:
        # A -log10 value (5.2) given where a probability belongs would otherwise be written as a nonsense column
        path = tmp_path / "bad.narrowPeak"
        peak = peak_class("chr1", 0, 100, **{field: bad})

        with pytest.raises(ValueError, match=field):
            save([peak], path, format="narrowpeak")
        with pytest.raises(ValueError, match=field):
            peak.to_bed_format("broadpeak")

    @pytest.mark.parametrize(("save", "peak_class"), [(save_chip_peaks, ChIPPeak), (save_atac_peaks, ATACPeak)])
    def test_a_failed_save_leaves_no_partial_file_and_keeps_an_existing_one(
        self, tmp_path: Path, save: Callable[..., Any], peak_class: Any
    ) -> None:
        peaks = [peak_class("chr1", 0, 100, q_value=0.01), peak_class("chr1", 200, 300, q_value=5.2)]

        fresh = tmp_path / "fresh.narrowPeak"
        with pytest.raises(ValueError):
            save(peaks, fresh)
        assert not fresh.exists()

        existing = tmp_path / "existing.narrowPeak"
        existing.write_text("keep me\n")
        with pytest.raises(ValueError):
            save(peaks, existing)
        assert existing.read_text() == "keep me\n"

    @pytest.mark.parametrize(("save", "peak_class"), [(save_chip_peaks, ChIPPeak), (save_atac_peaks, ATACPeak)])
    def test_bed_does_not_write_the_significance_columns(
        self, tmp_path: Path, save: Callable[..., Any], peak_class: Any
    ) -> None:
        path = tmp_path / "peaks.bed"

        save([peak_class("chr1", 0, 100, q_value=5.2)], path, format="bed")

        assert [len(row) for row in _data_rows(path)] == [6]


class TestBroadPeakRoundTrip:
    @pytest.mark.parametrize("assay", [CHIP, ATAC], ids=["chipseq", "atacseq"])
    def test_everything_but_the_summit_survives(self, tmp_path: Path, assay: Any) -> None:
        make_peaks, save, load = assay
        peaks = make_peaks()
        path = tmp_path / "peaks.broadPeak"

        save(peaks, path, format="broadpeak")
        loaded = load(path, format="broadpeak")

        assert len(loaded) == 4
        for original, back in zip(peaks, loaded):
            _assert_same_peak(back, original, keep_summit=False)

    @pytest.mark.parametrize("assay", [CHIP, ATAC], ids=["chipseq", "atacseq"])
    def test_lines_have_nine_columns(self, tmp_path: Path, assay: Any) -> None:
        make_peaks, save, _ = assay
        path = tmp_path / "peaks.broadPeak"

        save(make_peaks(), path, format="broadpeak")

        rows = _data_rows(path)
        assert all(len(row) == 9 for row in rows)
        assert rows[0][6:] == ["12.5", "5", "3"]
        assert rows[2][6:] == ["0", "-1", "-1"]

    def test_signal_is_not_read_from_the_summit_column(self, tmp_path: Path) -> None:
        # The old writer put the summit position in column 7: it came back as the signal value
        path = tmp_path / "peaks.broadPeak"
        save_chip_peaks([ChIPPeak("chr1", 100, 500, summit=300, signal_value=12.5)], path, format="broadpeak")

        (back,) = load_chip_peaks(path, format="broadpeak")

        assert back.signal_value == 12.5

    @pytest.mark.parametrize("loader", [load_chip_peaks, load_atac_peaks], ids=["chipseq", "atacseq"])
    def test_a_hand_written_broadpeak_line_is_read_as_probabilities(
        self, tmp_path: Path, loader: Callable[..., Any]
    ) -> None:
        path = tmp_path / "hand.broadPeak"
        path.write_text("chr1\t100\t900\tb1\t40\t-\t3.5\t6\t2\n")

        (peak,) = loader(path, format="broadpeak")

        assert (peak.signal_value, peak.strand, peak.summit) == (3.5, "-", None)
        assert peak.p_value == pytest.approx(1e-6, rel=1e-12)
        assert peak.q_value == pytest.approx(0.01, rel=1e-12)


class TestBedRoundTrip:
    def test_chipseq_keeps_position_score_and_strand(self, tmp_path: Path) -> None:
        path = tmp_path / "peaks.bed"
        save_chip_peaks(_chip_peaks(), path, format="bed")

        assert all(len(row) == 6 for row in _data_rows(path))
        loaded = load_chip_peaks(path, format="bed")

        assert [(p.chromosome, p.start, p.end, p.score, p.strand) for p in loaded] == [
            ("chr1", 100, 500, 50.0, "+"),
            ("chr2", 1000, 1400, 80.0, "-"),
            ("chrX", 5, 25, 1234.5, "."),
            ("chr3", 10, 20, 0.0, "."),
        ]

    def test_atacseq_uses_the_score_as_the_signal(self, tmp_path: Path) -> None:
        path = tmp_path / "peaks.bed"
        save_atac_peaks(_atac_peaks(), path, format="bed")

        assert all(len(row) == 6 for row in _data_rows(path))
        loaded = load_atac_peaks(path, format="bed")

        assert [(p.start, p.end, p.score, p.signal_value) for p in loaded] == [
            (100, 500, 50.0, 50.0),
            (1000, 1400, 80.0, 80.0),
            (5, 25, 1234.5, 1234.5),
            (10, 20, 0.0, 0.0),
        ]

    @pytest.mark.parametrize("loader", [load_chip_peaks, load_atac_peaks], ids=["chipseq", "atacseq"])
    def test_a_bed3_file_loads(self, tmp_path: Path, loader: Callable[..., Any]) -> None:
        path = tmp_path / "regions.bed"
        path.write_text("chr1\t100\t200\nchr2\t300\t450\n")

        loaded = loader(path, format="bed")

        assert [(p.chromosome, p.start, p.end, p.score, p.strand) for p in loaded] == [
            ("chr1", 100, 200, 0.0, "."),
            ("chr2", 300, 450, 0.0, "."),
        ]

    def test_atacseq_reads_a_short_line_as_bed(self, tmp_path: Path) -> None:
        # A BED6 file read with the default layout still loads (the score stands in as the signal value)
        path = tmp_path / "bed6.bed"
        path.write_text("chr1\t100\t200\tp1\t7\t+\n")

        (peak,) = load_atac_peaks(path)

        assert (peak.start, peak.end, peak.score, peak.signal_value, peak.strand) == (100, 200, 7.0, 7.0, "+")


class TestFormatNames:
    @pytest.mark.parametrize("assay", [CHIP, ATAC], ids=["chipseq", "atacseq"])
    def test_names_are_case_insensitive(self, tmp_path: Path, assay: Any) -> None:
        # "narrowPeak" is the standard spelling; it used to match no parser and give an empty list
        make_peaks, save, load = assay
        path = tmp_path / "peaks.narrowPeak"

        save(make_peaks(), path, format="NarrowPeak")

        assert len(load(path, format="narrowPeak")) == 4

    @pytest.mark.parametrize("assay", [CHIP, ATAC], ids=["chipseq", "atacseq"])
    def test_unknown_format_is_an_error_not_an_empty_file(self, tmp_path: Path, assay: Any) -> None:
        make_peaks, save, load = assay
        path = tmp_path / "peaks.txt"

        with pytest.raises(ValueError, match="Unsupported peak format"):
            save(make_peaks(), path, format="gff")
        assert not path.exists()

        path.write_text("chr1\t1\t2\n")
        with pytest.raises(ValueError, match="Unsupported peak format"):
            load(path, format="gff")
