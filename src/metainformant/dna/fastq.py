from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

from metainformant.core.io import open_text_auto


def _phred33_to_q(ch: str) -> int:
    return ord(ch) - 33


def average_phred_by_position(path: Path | str) -> list[float]:
    """Compute average Phred+33 quality score at each position across reads.

    - Streams the file (supports ``.gz``) via ``open_text_auto``
    - Truncates to the shortest read length encountered for safety
    - Ignores incomplete trailing records
    """
    sums: List[float] = []
    count = 0
    min_len: Optional[int] = None
    for _id, seq, qual in iter_fastq(path):
        qlen = len(qual)
        if min_len is None:
            min_len = qlen
            sums = [0.0] * min_len
        else:
            if qlen < min_len:
                # Truncate sums to the new shortest
                sums = sums[:qlen]
                min_len = qlen
            # If qlen > min_len, we only accumulate up to min_len
        # Accumulate
        upto = min_len if min_len is not None else 0
        for j in range(upto):
            sums[j] += float(_phred33_to_q(qual[j]))
        count += 1
    if count == 0:
        return []
    return [s / float(count) for s in sums]


@dataclass(frozen=True)
class FastqRecord:
    read_id: str
    sequence: str
    quality: str


def iter_fastq(path: Path | str) -> Iterator[Tuple[str, str, str]]:
    """Yield FASTQ records as (read_id, sequence, quality).

    - Supports plain text and ``.gz`` files
    - Minimal validation (``@`` on header, ``+`` separator); skips incomplete trailing records
    - Quality string is not validated beyond length consistency with sequence
    """
    with open_text_auto(path, mode="rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                return  # incomplete trailing record
            header = header.rstrip("\n")
            seq = seq.rstrip("\n")
            plus = plus.rstrip("\n")
            qual = qual.rstrip("\n")
            if not header.startswith("@"):
                # Skip malformed block
                continue
            if not plus.startswith("+"):
                # Skip malformed block
                continue
            # Trim spaces in ID, keep up to first whitespace
            read_id = header[1:].split()[0]
            # Ensure lengths match; if not, take the min for safety
            if len(seq) != len(qual):
                L = min(len(seq), len(qual))
                yield (read_id, seq[:L], qual[:L])
            else:
                yield (read_id, seq, qual)


def head(path: Path | str, n: int = 5) -> List[FastqRecord]:
    """Return the first ``n`` records as ``FastqRecord`` objects."""
    out: List[FastqRecord] = []
    for read_id, seq, qual in iter_fastq(path):
        out.append(FastqRecord(read_id, seq, qual))
        if len(out) >= n:
            break
    return out


def gc_content(seq: str) -> float:
    """Compute GC fraction of a nucleotide sequence (A/C/G/T/N-insensitive)."""
    if not seq:
        return 0.0
    g = seq.count("G") + seq.count("g")
    c = seq.count("C") + seq.count("c")
    atgc = sum(seq.upper().count(x) for x in ("A", "C", "G", "T"))
    if atgc == 0:
        return 0.0
    return (g + c) / float(atgc)


def summarize_fastq(path: Path | str, max_reads: Optional[int] = None) -> Dict[str, object]:
    """Compute lightweight summary statistics over a FASTQ file.

    Returns a dictionary with keys:
    - ``num_reads``: total reads processed (bounded by ``max_reads`` if provided)
    - ``length_min``, ``length_max``, ``length_mean``
    - ``gc_mean``: average GC fraction across reads
    - ``n_content_mean``: average fraction of 'N' bases across reads
    - ``avg_phred_by_pos``: list[float] average Phred+33 per position (truncated to shortest)

    Streaming implementation; processes ``.gz`` transparently.
    """
    num_reads = 0
    min_len: Optional[int] = None
    max_len: int = 0
    sum_len: int = 0
    sum_gc: float = 0.0
    sum_n_frac: float = 0.0
    phred_sums: List[float] = []
    phred_min_len: Optional[int] = None

    for _id, seq, qual in iter_fastq(path):
        L = len(seq)
        num_reads += 1
        sum_len += L
        max_len = max(max_len, L)
        min_len = L if min_len is None else min(min_len, L)
        # GC and N content
        sum_gc += gc_content(seq)
        n_frac = 0.0 if L == 0 else (seq.upper().count("N") / float(L))
        sum_n_frac += n_frac
        # Phred accumulation (truncate to shortest)
        if phred_min_len is None:
            phred_min_len = len(qual)
            phred_sums = [0.0] * phred_min_len
        else:
            if len(qual) < phred_min_len:
                phred_sums = phred_sums[: len(qual)]
                phred_min_len = len(qual)
        upto = phred_min_len if phred_min_len is not None else 0
        for j in range(upto):
            phred_sums[j] += float(_phred33_to_q(qual[j]))
        if max_reads is not None and num_reads >= max_reads:
            break

    if num_reads == 0:
        return {
            "num_reads": 0,
            "length_min": 0,
            "length_max": 0,
            "length_mean": 0.0,
            "gc_mean": 0.0,
            "n_content_mean": 0.0,
            "avg_phred_by_pos": [],
        }

    length_mean = sum_len / float(num_reads)
    gc_mean = sum_gc / float(num_reads)
    n_mean = sum_n_frac / float(num_reads)
    avg_phred = [s / float(num_reads) for s in phred_sums]
    return {
        "num_reads": num_reads,
        "length_min": int(min_len if min_len is not None else 0),
        "length_max": int(max_len),
        "length_mean": float(length_mean),
        "gc_mean": float(gc_mean),
        "n_content_mean": float(n_mean),
        "avg_phred_by_pos": avg_phred,
    }
