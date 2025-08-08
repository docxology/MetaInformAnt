from __future__ import annotations

from pathlib import Path


def _phred33(qchar: str) -> int:
    return max(0, ord(qchar) - 33)


def average_phred_by_position(path: str | Path) -> list[float]:
    """Compute average Phred+33 quality per position from a FASTQ file.

    Reads are assumed to be uniform length; if they differ, positions beyond
    shorter reads are ignored.
    """
    p = Path(path)
    lines = p.read_text().splitlines()
    # FASTQ blocks of 4 lines: header, seq, plus, qual
    qual_strings: list[str] = []
    for i in range(0, len(lines), 4):
        if i + 3 >= len(lines):
            break
        qual = lines[i + 3].rstrip("\n")
        if qual:
            qual_strings.append(qual)
    if not qual_strings:
        return []
    L = min(len(q) for q in qual_strings)
    sums = [0] * L
    counts = [0] * L
    for q in qual_strings:
        for j in range(L):
            sums[j] += _phred33(q[j])
            counts[j] += 1
    return [s / c if c else 0.0 for s, c in zip(sums, counts)]



from __future__ import annotations

from pathlib import Path
from typing import List


def average_phred_by_position(path: str | Path) -> List[float]:
    """Compute average Phred quality per position (Sanger, +33).

    Reads FASTQ naive four-line chunks; for tests only.
    """
    qualities: list[list[int]] = []
    with open(path, "r") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            plus = fh.readline()
            qual = fh.readline().rstrip("\n")
            if not qual:
                break
            for i, ch in enumerate(qual):
                q = ord(ch) - 33
                if i >= len(qualities):
                    qualities.append([])
                qualities[i].append(q)
    return [sum(col) / len(col) if col else 0.0 for col in qualities]


