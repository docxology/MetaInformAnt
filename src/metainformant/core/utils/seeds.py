"""Deterministic seed derivation shared across domains.

Pure stdlib utilities with no dependencies beyond the standard library, so any
domain module (simulation, math, rna, ...) can import them from
``metainformant.core`` without crossing domain boundaries.
"""

from __future__ import annotations

import hashlib
from typing import List


def deterministic_replicate_seeds(base_seed: int, n_replicates: int) -> List[int]:
    """Stable, order-independent replicate seeds derived from one base seed.

    Campaign simulation sweeps need each replicate to have a reproducible but
    distinct seed, and the mapping must not change when replicates are run in
    a different order or distribution. Seeds are derived with SHA-256 over
    ``f"{base_seed}:{replicate_index}"`` so the assignment is deterministic
    across Python processes and platforms (unlike ``random.Random`` streams).

    Args:
        base_seed: Root seed for the sweep.
        n_replicates: Number of replicate seeds to derive (>= 1).

    Returns:
        List of ``n_replicates`` non-negative 63-bit integer seeds.

    Raises:
        ValueError: If ``n_replicates`` < 1.
    """
    if n_replicates < 1:
        raise ValueError("n_replicates must be >= 1")
    seeds: List[int] = []
    for index in range(n_replicates):
        digest = hashlib.sha256(f"{base_seed}:{index}".encode("utf-8")).digest()
        seeds.append(int.from_bytes(digest[:8], "big") & 0x7FFFFFFFFFFFFFFF)
    return seeds
