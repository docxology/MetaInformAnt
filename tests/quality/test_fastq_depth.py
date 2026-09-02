"""Depth tests for FASTQ IO: record parsing, per-read metrics, and filtering.

All tests operate on real FASTQ files written to temporary directories and
parsed back through the public API.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.quality.io.fastq import (
    FastqRecord,
    filter_reads,
    read_fastq_records,
)


def _write_fastq(path: Path, records: list[tuple[str, str, str]]) -> Path:
    lines: list[str] = []
    for header, seq, qual in records:
        lines.extend(["@" + header, seq, "+", qual])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _record(name: str, seq: str, qual: str) -> tuple[str, str, str]:
    assert len(seq) == len(qual)
    return name, seq, qual


class TestFastqRecord:
    def test_name_length_quality_scores(self) -> None:
        rec = FastqRecord("@read1 description", "ACGT", "+", "IIII")
        assert rec.name == "read1"
        assert rec.length == 4
        assert rec.quality_scores() == [40, 40, 40, 40]
        assert rec.mean_quality() == pytest.approx(40.0)

    def test_gc_content_percent_scale(self) -> None:
        rec = FastqRecord("@r", "GGCC", "+", "IIII")
        assert rec.gc_content() == pytest.approx(100.0)
        at_only = FastqRecord("@r", "AATT", "+", "IIII")
        assert at_only.gc_content() == pytest.approx(0.0)
        mixed = FastqRecord("@r", "AGCT", "+", "IIII")
        assert mixed.gc_content() == pytest.approx(50.0)
        empty = FastqRecord("@r", "", "+", "")
        assert empty.gc_content() == 0.0

    def test_mean_quality_mixed_scores(self) -> None:
        # ASCII: I=40, #=2, 5=20
        rec = FastqRecord("@r", "ACGT", "+", "I#5I")
        assert rec.mean_quality() == pytest.approx((40 + 2 + 20 + 40) / 4)


class TestReadFastqRecords:
    def test_roundtrip_multi_record(self, tmp_path: Path) -> None:
        records = [
            _record("r1", "ACGTACGTAC", "IIIIIIIIII"),
            _record("r2", "GGGGCCCCAA", "##########"),
            _record("r3", "TTTTAAAACC", "5" * 10),
        ]
        path = _write_fastq(tmp_path / "in.fastq", records)
        parsed = list(read_fastq_records(path))
        assert [r.name for r in parsed] == ["r1", "r2", "r3"]
        assert parsed[1].quality_scores() == [2] * 10
        assert parsed[2].mean_quality() == pytest.approx(20.0)

    def test_max_records_limit(self, tmp_path: Path) -> None:
        records = [_record(f"r{i}", "ACGT", "IIII") for i in range(10)]
        path = _write_fastq(tmp_path / "in.fastq", records)
        parsed = list(read_fastq_records(path, max_records=3))
        assert len(parsed) == 3

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(Exception):
            list(read_fastq_records(tmp_path / "nope.fastq"))


class TestFilterReads:
    def _make_file(self, tmp_path: Path) -> tuple[Path, list[tuple[str, str, str]]]:
        # r1: high quality, long, no N -> passes everything
        # r2: low mean quality
        # r3: too short
        # r4: too many N bases
        records = [
            _record("r1", "ACGTACGTACGT", "IIIIIIIIIIII"),
            _record("r2", "ACGTACGTACGT", "##########ii"),
            _record("r3", "ACGT", "IIII"),
            _record("r4", "ACGTNACGTNAC", "IIIIIIIIIIII"),
        ]
        path = _write_fastq(tmp_path / "in.fastq", records)
        return path, records

    def test_filters_by_quality_length_and_n(self, tmp_path: Path) -> None:
        path, _ = self._make_file(tmp_path)
        out = tmp_path / "out.fastq"
        stats = filter_reads(path, out, min_quality=25.0, min_length=10, max_n_bases=1)
        assert stats["total_reads"] == 4
        assert stats["passed_reads"] == 1
        assert stats["filtered_reads"] == 3
        assert stats["pass_rate"] == pytest.approx(25.0)
        assert stats["filters_applied"] == {"min_quality": 25.0, "min_length": 10, "max_n_bases": 1}

        survivors = list(read_fastq_records(out))
        assert [r.name for r in survivors] == ["r1"]
        # Roundtrip fidelity
        assert survivors[0].sequence == "ACGTACGTACGT"

    def test_no_length_filter(self, tmp_path: Path) -> None:
        path, _ = self._make_file(tmp_path)
        out = tmp_path / "out.fastq"
        stats = filter_reads(path, out, min_quality=25.0, min_length=None, max_n_bases=1)
        # No length floor: short r3 passes; r4 still fails the N cap.
        assert stats["filters_applied"]["min_length"] is None
        assert [r.name for r in read_fastq_records(out)] == ["r1", "r3"]
        assert stats["passed_reads"] == 2

    def test_negative_max_n_disables_n_filter(self, tmp_path: Path) -> None:
        path, _ = self._make_file(tmp_path)
        out = tmp_path / "out.fastq"
        stats = filter_reads(path, out, min_quality=25.0, min_length=10, max_n_bases=-1)
        assert stats["passed_reads"] == 2
        assert [r.name for r in read_fastq_records(out)] == ["r1", "r4"]

    def test_empty_input(self, tmp_path: Path) -> None:
        path = tmp_path / "empty.fastq"
        path.write_text("", encoding="utf-8")
        out = tmp_path / "out.fastq"
        stats = filter_reads(path, out, min_quality=20.0, min_length=10, max_n_bases=0)
        assert stats == {
            "total_reads": 0,
            "passed_reads": 0,
            "filtered_reads": 0,
            "pass_rate": 0,
            "filters_applied": {"min_quality": 20.0, "min_length": 10, "max_n_bases": 0},
        }
        assert out.exists()
