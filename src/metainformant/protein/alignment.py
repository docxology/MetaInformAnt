from __future__ import annotations

from typing import Tuple


def pairwise_identity(a: str, b: str) -> float:
    """Compute simple pairwise identity for equal-length strings.

    If lengths differ, compare up to min length.
    """
    if not a or not b:
        return 0.0
    L = min(len(a), len(b))
    if L == 0:
        return 0.0
    same = sum(1 for x, y in zip(a[:L], b[:L]) if x == y)
    return same / L


def needleman_wunsch(a: str, b: str, *, match: int = 1, mismatch: int = -1, gap: int = -1) -> Tuple[int, str, str]:
    """Very small Needleman–Wunsch global alignment implementation returning (score, align_a, align_b)."""
    n = len(a)
    m = len(b)
    # DP matrices
    score = [[0] * (m + 1) for _ in range(n + 1)]
    pointer = [[0] * (m + 1) for _ in range(n + 1)]  # 0 diag, 1 up, 2 left
    for i in range(1, n + 1):
        score[i][0] = i * gap
        pointer[i][0] = 1
    for j in range(1, m + 1):
        score[0][j] = j * gap
        pointer[0][j] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap
            best = max(diag, up, left)
            score[i][j] = best
            pointer[i][j] = 0 if best == diag else (1 if best == up else 2)
    # traceback
    i, j = n, m
    align_a = []
    align_b = []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and pointer[i][j] == 0:
            align_a.append(a[i - 1])
            align_b.append(b[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or pointer[i][j] == 1):
            align_a.append(a[i - 1])
            align_b.append("-")
            i -= 1
        else:
            align_a.append("-")
            align_b.append(b[j - 1])
            j -= 1
    align_a.reverse()
    align_b.reverse()
    return score[n][m], "".join(align_a), "".join(align_b)


