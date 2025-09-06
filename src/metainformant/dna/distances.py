from __future__ import annotations

import math
from collections import Counter
from typing import Dict, List, Sequence


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


def kimura_2p_distance(s1: str, s2: str) -> float:
    """Kimura 2-parameter distance model accounting for transitions/transversions.

    K = -0.5 * ln((1 - 2P - Q) * sqrt(1 - 2Q))
    where P = transition frequency, Q = transversion frequency
    """
    if not s1 or not s2:
        return 0.0
    L = min(len(s1), len(s2))
    if L == 0:
        return 0.0

    transitions = 0  # A<->G, C<->T
    transversions = 0  # A<->C, A<->T, G<->C, G<->T

    for a, b in zip(s1[:L].upper(), s2[:L].upper()):
        if a != b:
            if (a in "AG" and b in "AG") or (a in "CT" and b in "CT"):
                transitions += 1
            else:
                transversions += 1

    P = transitions / L
    Q = transversions / L

    if 1 - 2 * P - Q <= 0 or 1 - 2 * Q <= 0:
        return float("inf")

    return -0.5 * math.log((1 - 2 * P - Q) * math.sqrt(1 - 2 * Q))


def _kmer_counts(seq: str, k: int, *, dna_only: bool = True) -> Dict[str, int]:
    """Count k-mers in a sequence with optional DNA filtering."""
    if k <= 0:
        raise ValueError("k must be positive")
    seq_u = seq.upper()
    counts: Dict[str, int] = Counter()

    for i in range(max(0, len(seq_u) - k + 1)):
        kmer = seq_u[i : i + k]
        if dna_only and any(ch not in "ACGT" for ch in kmer):
            continue
        counts[kmer] += 1
    return dict(counts)


def kmer_distance(a: str, b: str, *, k: int = 3, metric: str = "cosine", dna_only: bool = True) -> float:
    """Compute distance between two sequences using k-mer frequency vectors.

    Args:
        a, b: Input sequences
        k: K-mer length
        metric: Distance metric ('cosine', 'euclidean', 'manhattan', 'jaccard')
        dna_only: If True, filter out non-ACGT k-mers

    Returns:
        Distance value (0 = identical, higher = more different)
    """
    ca = _kmer_counts(a, k, dna_only=dna_only)
    cb = _kmer_counts(b, k, dna_only=dna_only)

    if not ca and not cb:
        return 0.0

    vocab = set(ca.keys()) | set(cb.keys())
    va = [ca.get(km, 0) for km in vocab]
    vb = [cb.get(km, 0) for km in vocab]

    if metric == "euclidean":
        return math.sqrt(sum((x - y) ** 2 for x, y in zip(va, vb)))
    elif metric == "manhattan":
        return sum(abs(x - y) for x, y in zip(va, vb))
    elif metric == "jaccard":
        intersection = sum(min(x, y) for x, y in zip(va, vb))
        union = sum(max(x, y) for x, y in zip(va, vb))
        return 1.0 - (intersection / union if union > 0 else 0.0)
    else:  # cosine
        dot = sum(x * y for x, y in zip(va, vb))
        na = math.sqrt(sum(x * x for x in va))
        nb = math.sqrt(sum(y * y for y in vb))
        if na == 0 or nb == 0:
            return 1.0
        cos_sim = dot / (na * nb)
        return max(0.0, 1.0 - cos_sim)


def kmer_distance_matrix(id_to_seq: Dict[str, str], *, k: int = 3, metric: str = "cosine") -> List[List[float]]:
    """Compute pairwise k-mer distance matrix for a set of sequences."""
    ids = list(id_to_seq.keys())
    n = len(ids)
    matrix: List[List[float]] = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i):
            d = kmer_distance(id_to_seq[ids[i]], id_to_seq[ids[j]], k=k, metric=metric)
            matrix[i][j] = matrix[j][i] = d

    return matrix


def tajima_nei_distance(s1: str, s2: str) -> float:
    """Tajima-Nei distance accounting for unequal base frequencies."""
    if not s1 or not s2:
        return 0.0

    L = min(len(s1), len(s2))
    if L == 0:
        return 0.0

    # Count base frequencies in both sequences
    counts1 = Counter(s1[:L].upper())
    counts2 = Counter(s2[:L].upper())

    # Calculate average base frequencies
    bases = "ACGT"
    freq = {}
    for base in bases:
        freq[base] = (counts1.get(base, 0) + counts2.get(base, 0)) / (2 * L)

    # Calculate p-distance
    p = p_distance(s1, s2)

    # Calculate b (correction factor)
    b = sum(freq[i] * freq[j] for i in bases for j in bases if i != j)

    if b == 0 or p >= b:
        return float("inf")

    return -b * math.log(1 - p / b)


def sequence_identity_matrix(sequences: Sequence[str]) -> List[List[float]]:
    """Compute sequence identity matrix (1 - p_distance)."""
    n = len(sequences)
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        matrix[i][i] = 1.0  # Identity with self
        for j in range(i + 1, n):
            identity = 1.0 - p_distance(sequences[i], sequences[j])
            matrix[i][j] = matrix[j][i] = identity

    return matrix
