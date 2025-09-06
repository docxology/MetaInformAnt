from __future__ import annotations

import random
from typing import Dict


def apply_point_mutations(seq: str, changes: Dict[int, str]) -> str:
    """Apply point mutations specified by a mapping index->new_base.

    Indices are 0-based. Out-of-range indices are ignored. The new base is
    inserted as given (case preserved).
    """
    if not seq or not changes:
        return seq
    arr = list(seq)
    n = len(arr)
    for idx, base in changes.items():
        if 0 <= idx < n and base:
            arr[idx] = str(base)[0]
    return "".join(arr)


def hamming_distance(a: str, b: str) -> int:
    """Hamming distance up to the shorter length."""
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] != b[i])


def random_point_mutations(seq: str, *, num_mutations: int, seed: int | None = None) -> str:
    """Introduce a given number of random point mutations into the sequence.

    Only A/C/G/T are used as replacement bases when possible. Positions are
    unique; if num_mutations exceeds length, it is clamped.
    """
    if not seq or num_mutations <= 0:
        return seq
    rng = random.Random(seed)
    n = len(seq)
    k = min(num_mutations, n)
    positions = rng.sample(range(n), k)
    dna_bases = "ACGT"
    arr = list(seq)
    for pos in positions:
        current_upper = arr[pos].upper()
        options = [b for b in dna_bases if b != current_upper] or list(dna_bases)
        new_base = rng.choice(options)
        arr[pos] = new_base if arr[pos].isupper() else new_base.lower()
    return "".join(arr)
