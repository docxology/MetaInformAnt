from __future__ import annotations

import math
from collections import Counter
from typing import Dict, List
from collections import Counter
from typing import Dict, List


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


def _kmer_counts(seq: str, k: int) -> Dict[str, int]:
    if k <= 0:
        raise ValueError("k must be positive")
    seq_u = seq.upper()
    counts: Dict[str, int] = Counter()
    for i in range(0, max(0, len(seq_u) - k + 1)):
        kmer = seq_u[i : i + k]
        if any(ch not in "ACGT" for ch in kmer):
            continue
        counts[kmer] += 1
    return counts


def kmer_distance(a: str, b: str, *, k: int = 3, metric: str = "cosine") -> float:
    """Distance between two sequences using k-mer vectors.

    Supports metric = 'cosine' or 'euclidean'.
    """
    ca = _kmer_counts(a, k)
    cb = _kmer_counts(b, k)
    keys = set(ca) | set(cb)
    if not keys:
        return 0.0
    va = [ca.get(km, 0) for km in keys]
    vb = [cb.get(km, 0) for km in keys]
    if metric == "euclidean":
        return math.sqrt(sum((x - y) ** 2 for x, y in zip(va, vb)))
    # cosine
    dot = sum(x * y for x, y in zip(va, vb))
    na = math.sqrt(sum(x * x for x in va))
    nb = math.sqrt(sum(y * y for y in vb))
    if na == 0 or nb == 0:
        return 1.0
    cos_sim = dot / (na * nb)
    return max(0.0, 1.0 - cos_sim)


def kmer_distance_matrix(id_to_seq: Dict[str, str], *, k: int = 3, metric: str = "cosine") -> List[List[float]]:
    ids = list(id_to_seq.keys())
    n = len(ids)
    dm = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i):
            d = kmer_distance(id_to_seq[ids[i]], id_to_seq[ids[j]], k=k, metric=metric)
            dm[i][j] = dm[j][i] = d
    return dm

def _kmer_counts(seq: str, k: int) -> Dict[str, int]:
    if k <= 0 or k > len(seq):
        return {}
    counts: Dict[str, int] = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def kmer_distance(a: str, b: str, *, k: int = 3, metric: str = "cosine") -> float:
    """Compute distance between two sequences using k-mer frequency vectors.

    Supported metric: cosine distance = 1 - cosine similarity.
    """
    ca = _kmer_counts(a, k)
    cb = _kmer_counts(b, k)
    if not ca and not cb:
        return 0.0
    # cosine similarity
    vocab = set(ca) | set(cb)
    dot = sum(ca.get(w, 0) * cb.get(w, 0) for w in vocab)
    na = math.sqrt(sum(v * v for v in ca.values()))
    nb = math.sqrt(sum(v * v for v in cb.values()))
    if na == 0 or nb == 0:
        return 1.0
    cos = dot / (na * nb)
    return 1.0 - cos


def kmer_distance_matrix(id_to_seq: Dict[str, str], *, k: int = 3, metric: str = "cosine") -> List[List[float]]:
    ids = list(id_to_seq.keys())
    n = len(ids)
    matrix: List[List[float]] = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i):
            d = kmer_distance(id_to_seq[ids[i]], id_to_seq[ids[j]], k=k, metric=metric)
            matrix[i][j] = matrix[j][i] = d
    return matrix


