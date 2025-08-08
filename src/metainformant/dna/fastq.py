from __future__ import annotations

from pathlib import Path


def _phred33_to_q(ch: str) -> int:
    return ord(ch) - 33


def average_phred_by_position(path: Path | str) -> list[float]:
    """Compute average Phred+33 quality score at each position across reads.

    Assumes all reads are equal length; truncates to the shortest among reads
    for safety.
    """
    p = Path(path)
    text = p.read_text().splitlines()
    # FASTQ comes in 4-line records
    quals: list[str] = []
    for i in range(0, len(text), 4):
        if i + 3 >= len(text):
            break
        qual = text[i + 3].rstrip("\n")
        quals.append(qual)
    if not quals:
        return []
    L = min(len(q) for q in quals)
    sums = [0.0] * L
    for qline in quals:
        for j in range(L):
            sums[j] += float(_phred33_to_q(qline[j]))
    return [s / len(quals) for s in sums]


