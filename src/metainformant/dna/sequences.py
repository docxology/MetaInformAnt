from __future__ import annotations

from collections import Counter
from typing import Dict

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
    """Calculate GC content of a DNA sequence."""
    if not seq:
        return 0.0

    # Convert to uppercase for consistent counting
    upper = seq.upper()

    # Count GC bases efficiently
    gc_count = upper.count("G") + upper.count("C")

    # Handle non-standard characters by counting only ACGT
    total_valid = sum(upper.count(base) for base in "ACGT")

    return gc_count / total_valid if total_valid > 0 else 0.0


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
