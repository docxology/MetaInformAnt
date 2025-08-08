from __future__ import annotations

from collections import Counter
from typing import Dict


def consensus_from_alignment(id_to_aligned_seq: Dict[str, str]) -> str:
    if not id_to_aligned_seq:
        return ""
    seqs = list(id_to_aligned_seq.values())
    L = min(len(s) for s in seqs)
    out: list[str] = []
    for pos in range(L):
        col = [s[pos] for s in seqs]
        col = [c for c in col if c != "-"]
        if not col:
            out.append("-")
            continue
        c, _ = Counter(col).most_common(1)[0]
        out.append(c)
    return "".join(out)

from __future__ import annotations

from collections import Counter
from typing import Dict


def consensus_from_alignment(id_to_seq: Dict[str, str]) -> str:
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




