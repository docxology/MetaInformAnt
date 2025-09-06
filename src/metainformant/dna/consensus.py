from __future__ import annotations

from collections import Counter


def consensus_from_alignment(id_to_seq: dict[str, str]) -> str:
    """Return a simple majority-rule consensus sequence.

    Gaps ('-') are ignored when determining the majority at each position.
    If no non-gap characters exist at a position, a '-' is emitted.
    Assumes all sequences are of equal length; if not, truncates to the
    minimum length.
    """
    if not id_to_seq:
        return ""
    lengths = [len(s) for s in id_to_seq.values()]
    if not lengths:
        return ""
    L = min(lengths)
    consensus_chars: list[str] = []
    for pos in range(L):
        column_chars = [seq[pos] for seq in id_to_seq.values() if seq[pos] != "-"]
        if not column_chars:
            consensus_chars.append("-")
            continue
        counts = Counter(column_chars)
        # deterministic tie-break by lexicographic order
        best_char = max(sorted(counts.keys()), key=lambda c: counts[c])
        consensus_chars.append(best_char)
    return "".join(consensus_chars)
