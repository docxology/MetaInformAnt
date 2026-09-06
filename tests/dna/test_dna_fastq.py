"""Tests for DNA FASTQ file handling."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from metainformant.dna.io import fastq


def test_average_phred_by_position(tmp_path: Path) -> None:
    """Test average Phred quality calculation by position."""
    content = """@r1
ACGTN
+
IIIII
@r2
ACGTN
+
IIIII
"""
    p = tmp_path / "reads.fq"
    p.write_text(content)
    avgs = fastq.average_phred_by_position(p)
    assert set(avgs.keys()) == {0, 1, 2, 3, 4}
    assert all(abs(x - 40.0) < 1e-9 for x in avgs.values())


def test_iter_fastq_and_head(tmp_path: Path) -> None:
    """Test FASTQ iteration and head functions."""
    content = """@r1 first read
ACGTN
+
IIIII
@r2 second read
GGCCTA
+
HHHHHH
"""
    p = tmp_path / "reads.fq"
    p.write_text(content)
    it = list(fastq.iter_fastq(p))
    assert len(it) == 2
    assert it[0][0] == "r1"
    assert it[0][1] == "ACGTN"
    assert it[0][2] == "IIIII"


def test_iter_fastq_gz(tmp_path: Path) -> None:
    """Test FASTQ iteration from gzipped file."""
    content = """@a
ACGT
+
IIII
@b
TTTT
+
HHHH
"""
    p = tmp_path / "reads.fastq.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(content)
    records = list(fastq.iter_fastq(p))
    assert len(records) == 2


def test_summarize_fastq(tmp_path: Path) -> None:
    """Test FASTQ file summarization."""
    content = """@x
ACGT
+
IIII
@y
AGGT
+
HHHH
"""
    p = tmp_path / "reads.fq"
    p.write_text(content)
    summary = fastq.summarize_fastq(p)
    assert summary["total_reads"] == 2
    assert summary["total_bases"] == 8
    assert summary["min_length"] == 4
    assert summary["max_length"] == 4
    assert summary["mean_length"] == 4.0


def test_read_write_fastq_roundtrip(tmp_path: Path) -> None:
    reads = {
        "r1": ("ACGTACGT", "IIIIIIII"),
        "r2": ("TTTTGGGG", "HHHHHHHH"),
    }
    out = tmp_path / "sub" / "out.fq"
    fastq.write_fastq(reads, out)
    parsed = fastq.read_fastq(out)
    assert parsed == reads


def test_read_fastq_gzipped(tmp_path: Path) -> None:
    p = tmp_path / "reads.fastq.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("@r1\nACGT\n+\nIIII\n")
    assert fastq.read_fastq(p) == {"r1": ("ACGT", "IIII")}


def test_read_fastq_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        fastq.read_fastq(tmp_path / "absent.fq")


def test_read_fastq_rejects_bad_header(tmp_path: Path) -> None:
    p = tmp_path / "bad.fq"
    p.write_text("NOT_A_HEADER\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError):
        fastq.read_fastq(p)


def test_read_fastq_rejects_length_mismatch(tmp_path: Path) -> None:
    p = tmp_path / "mismatch.fq"
    p.write_text("@r1\nACGT\n+\nII\n")
    with pytest.raises(ValueError):
        fastq.read_fastq(p)


def test_write_fastq_rejects_length_mismatch(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        fastq.write_fastq({"r1": ("ACGT", "II")}, tmp_path / "out.fq")


def test_assess_quality_statistics(tmp_path: Path) -> None:
    # Q40 ('I') for ACGT, Q10 ('+') for the second read.
    p = tmp_path / "reads.fq"
    p.write_text("@r1\nACGT\n+\nIIII\n@r2\nGGCC\n+\n++++\n")
    stats = fastq.assess_quality(p)
    assert stats["total_reads"] == 2
    assert stats["mean_quality"] == pytest.approx(25.0)
    assert stats["gc_content"] == pytest.approx(75.0)  # 6 GC of 8 bases
    assert stats["read_lengths"] == [4, 4]


def test_assess_quality_empty_file(tmp_path: Path) -> None:
    p = tmp_path / "empty.fq"
    p.write_text("")
    stats = fastq.assess_quality(p)
    assert stats["total_reads"] == 0
    assert stats["mean_quality"] == 0.0


def test_filter_reads_by_min_quality(tmp_path: Path) -> None:
    p = tmp_path / "reads.fq"
    p.write_text("@good\nACGT\n+\nIIII\n@bad\nACGT\n+\n####\n")
    kept = list(fastq.filter_reads(p, min_quality=20))
    assert len(kept) == 1
    assert kept[0].split("\n")[0] == "@good"


def test_convert_fastq_to_fasta(tmp_path: Path) -> None:
    p = tmp_path / "reads.fq"
    p.write_text("@r1\nACGT\n+\nIIII\n")
    out = tmp_path / "out.fasta"
    fastq.convert_fastq_to_fasta(p, out)
    assert out.read_text() == ">r1\nACGT\n"


def test_trim_reads(tmp_path: Path) -> None:
    p = tmp_path / "reads.fq"
    p.write_text("@r1\nAAAAAACCCCCC\n+\nIIIIIIIIIIII\n")
    out = tmp_path / "trimmed.fq"
    fastq.trim_reads(p, out, min_length=4, trim_5p=2, trim_3p=2)
    trimmed = fastq.read_fastq(out)
    assert trimmed["r1"][0] == "AAAACCCC"
    assert trimmed["r1"][1] == "IIIIIIII"


def test_trim_reads_drops_too_short(tmp_path: Path) -> None:
    p = tmp_path / "reads.fq"
    p.write_text("@short\nAAA\n+\nIII\n")
    out = tmp_path / "trimmed.fq"
    fastq.trim_reads(p, out, min_length=50)
    assert fastq.read_fastq(out) == {}


def test_fastq_record_properties() -> None:
    rec = fastq.FastqRecord("r1", "ACGT", "IIII")
    assert rec.mean_quality == pytest.approx(40.0)
    assert rec.min_quality == 40
    assert rec.max_quality == 40
    assert rec.median_quality == 40
    assert rec.length() == 4
    assert len(rec) == 4
    assert rec.to_string() == "@r1\nACGT\n+\nIIII\n"


def test_fastq_record_rejects_mismatched_lengths() -> None:
    with pytest.raises(ValueError):
        fastq.FastqRecord("r1", "ACGT", "II")


def test_fastq_record_even_length_median() -> None:
    # Phred scores for "I!": [40, 0] -> median is the average (20).
    rec = fastq.FastqRecord("r1", "AC", "I!")
    assert rec.median_quality == 20


def test_fastq_record_gc_content() -> None:
    rec = fastq.FastqRecord("r1", "GGCC", "IIII")
    assert rec.gc_content() == pytest.approx(1.0)
    rec2 = fastq.FastqRecord("r2", "GAAT", "IIII")
    assert rec2.gc_content() == pytest.approx(0.25)


def test_fastq_record_trim_low_quality() -> None:
    # Phred scores for "II5!": [40, 40, 20, 0]. Trailing bases below the
    # threshold are trimmed from the 3' end.
    rec = fastq.FastqRecord("r1", "ACGT", "II5!")
    trimmed = rec.trim_low_quality(min_quality=20)
    assert trimmed.sequence == "ACG"
    assert trimmed.quality_string == "II5"
    # Nothing to trim when the last base passes the threshold.
    rec2 = fastq.FastqRecord("r2", "ACGT", "IIII")
    assert rec2.trim_low_quality(min_quality=20) is rec2


def test_fastq_record_from_string_roundtrip() -> None:
    text = "@r1\nACGT\n+\nIIII\n"
    rec = fastq.FastqRecord.from_string(text)
    assert rec.header == "r1"
    assert rec.to_string() == text


def test_fastq_record_from_string_rejects_malformed() -> None:
    with pytest.raises(ValueError):
        fastq.FastqRecord.from_string("@r1\nACGT\n+\n")
    with pytest.raises(ValueError):
        fastq.FastqRecord.from_string("@r1\nACGT\n-\nIIII\n")
