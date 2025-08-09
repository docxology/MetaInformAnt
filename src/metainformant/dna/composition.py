from __future__ import annotations

from typing import List


def gc_skew(seq: str) -> float:
    """(G - C) / (G + C); returns 0.0 when there is no G/C content or empty input."""
    if not seq:
        return 0.0
    s = seq.upper()
    g = s.count("G")
    c = s.count("C")
    denom = g + c
    if denom == 0:
        return 0.0
    return (g - c) / denom


def cumulative_gc_skew(seq: str) -> List[float]:
    """Cumulative GC skew walk: +1 for G, -1 for C, 0 otherwise; returns floats."""
    total = 0
    out: List[float] = []
    for ch in (seq or "").upper():
        if ch == "G":
            total += 1
        elif ch == "C":
            total -= 1
        out.append(float(total))
    return out


def melting_temperature(seq: str) -> float:
    """Very rough Tm estimate.

    - Short sequences (<= 14): Wallace rule 2*(A+T) + 4*(G+C)
    - Otherwise: 64.9 + 41*(G+C-16.4)/N
    """
    s = (seq or "").upper()
    n = len(s)
    if n == 0:
        return 0.0
    a = s.count("A")
    t = s.count("T")
    g = s.count("G")
    c = s.count("C")
    if n <= 14:
        return 2.0 * (a + t) + 4.0 * (g + c)
    return 64.9 + 41.0 * (g + c - 16.4) / n


