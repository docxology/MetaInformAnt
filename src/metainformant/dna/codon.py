from __future__ import annotations

from collections import Counter


def codon_counts(seq: str) -> dict[str, int]:
    """Count complete codons in the sequence (ignore trailing partial triplets).

    Returns a mapping of codon -> count. Input is uppercased; only full triplets are counted.
    """
    s = (seq or "").upper()
    counts: Counter[str] = Counter()
    for i in range(0, len(s) - len(s) % 3, 3):
        codon = s[i : i + 3]
        if len(codon) == 3:
            counts[codon] += 1
    return dict(counts)


def codon_frequencies(seq: str) -> dict[str, float]:
    """Compute relative frequencies of codons. Empty dict if there are no full codons."""
    counts = codon_counts(seq)
    total = sum(counts.values())
    if total == 0:
        return {}
    return {codon: count / total for codon, count in counts.items()}
