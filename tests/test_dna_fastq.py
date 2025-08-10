from __future__ import annotations

from pathlib import Path

from metainformant.dna import fastq


def test_average_phred_by_position(tmp_path: Path) -> None:
    content = """@r1
ACGT
+
IIII
@r2
ACGA
+
IIII
"""
    p = tmp_path / "reads.fq"
    p.write_text(content)
    avgs = fastq.average_phred_by_position(p)
    # 'I' is Phred+33 -> 40
    assert all(abs(x - 40.0) < 1e-9 for x in avgs)


def test_iter_fastq_and_head(tmp_path: Path) -> None:
    content = """@r1 first read
ACGTN
+
IIIII
@r2 second read
ACGAN
+
IIIII
"""
    p = tmp_path / "reads.fq"
    p.write_text(content)
    it = list(fastq.iter_fastq(p))
    assert len(it) == 2
    rid1, s1, q1 = it[0]
    assert rid1 == "r1"
    assert s1 == "ACGTN"
    assert q1 == "IIIII"
    head = fastq.head(p, n=1)
    assert len(head) == 1
    assert head[0].read_id == "r1"


def test_iter_fastq_gz(tmp_path: Path) -> None:
    import gzip

    content = """@a
ACGT
+
IIII
@b
ACGT
+
IIII
"""
    p = tmp_path / "reads.fastq.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(content)
    records = list(fastq.iter_fastq(p))
    assert len(records) == 2
    assert records[0][0] == "a"
    assert records[1][0] == "b"


def test_summarize_fastq(tmp_path: Path) -> None:
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
    assert summary["num_reads"] == 2
    assert summary["length_min"] == 4
    assert summary["length_max"] == 4
    assert abs(summary["gc_mean"] - ((2/4)+(3/4))/2) < 1e-9
    assert len(summary["avg_phred_by_pos"]) == 4

