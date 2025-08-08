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


