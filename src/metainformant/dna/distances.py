from __future__ import annotations

import math


def p_distance(s1: str, s2: str) -> float:
    """Proportion of differing sites between two sequences (Hamming/L)."""
    if not s1 or not s2:
        return 0.0
    L = min(len(s1), len(s2))
    if L == 0:
        return 0.0
    diff = sum(1 for a, b in zip(s1[:L], s2[:L]) if a != b)
    return diff / L


def jc69_distance(s1: str, s2: str) -> float:
    """Jukes-Cantor 69 distance from p-distance.

    d = -3/4 * ln(1 - 4/3 * p)
    """
    p = p_distance(s1, s2)
    x = 1.0 - (4.0 / 3.0) * p
    if x <= 0:
        return float("inf")
    return -0.75 * math.log(x)



from __future__ import annotations

import math


def p_distance(s1: str, s2: str) -> float:
    L = min(len(s1), len(s2))
    if L == 0:
        return 0.0
    diffs = sum(1 for a, b in zip(s1[:L], s2[:L]) if a != b)
    return diffs / L


def jc69_distance(s1: str, s2: str) -> float:
    p = p_distance(s1, s2)
    if p >= 0.75:
        # cap to avoid math domain error; returns a large distance
        return float("inf")
    return -3.0 / 4.0 * math.log(1 - 4 * p / 3)


