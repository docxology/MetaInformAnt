from __future__ import annotations

from typing import Dict
from collections import Counter

from Bio import SeqIO


def read_fasta(path: str) -> Dict[str, str]:
    """Read a FASTA file into a dictionary of id -> sequence string."""
    records = SeqIO.parse(path, "fasta")
    seqs: Dict[str, str] = {}
    for rec in records:
        seqs[rec.id] = str(rec.seq)
    return seqs


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    upper = seq.upper()
    gc = upper.count("G") + upper.count("C")
    return gc / len(upper)


def kmer_counts(seq: str, k: int) -> Dict[str, int]:
    if k <= 0 or len(seq) < k:
        return {}
    return dict(Counter(seq[i : i + k] for i in range(0, len(seq) - k + 1)))


def kmer_frequencies(seq: str, k: int) -> Dict[str, float]:
    counts = kmer_counts(seq, k)
    total = sum(counts.values())
    if total == 0:
        return {}
    return {kmer: cnt / total for kmer, cnt in counts.items()}


