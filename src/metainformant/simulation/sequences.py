from __future__ import annotations

import random
from typing import Iterable

_DNA = "ACGT"
_AA = "ACDEFGHIKLMNPQRSTVWY"


def generate_random_dna(length: int, *, gc_content: float = 0.5, rng: random.Random | None = None) -> str:
    """Generate random DNA sequence with specified GC content.

    Args:
        length: Length of sequence to generate
        gc_content: Target GC content (0.0-1.0)
        rng: Random number generator

    Returns:
        Random DNA sequence string

    Raises:
        ValueError: If gc_content is not in [0, 1] or length is negative
    """
    if not (0.0 <= gc_content <= 1.0):
        raise ValueError(f"GC content must be between 0.0 and 1.0, got {gc_content}")

    if length < 0:
        raise ValueError(f"Length must be non-negative, got {length}")

    if length == 0:
        return ""

    r = rng or random

    # Adjust base probabilities for GC content
    gc_prob = gc_content / 2  # Equal probability for G and C
    at_prob = (1 - gc_content) / 2  # Equal probability for A and T

    bases = ["A", "T", "G", "C"]
    cumulative_weights = [at_prob, at_prob + at_prob, at_prob + at_prob + gc_prob, 1.0]

    result = []
    for _ in range(length):
        rand_val = r.random()  # Get random float [0, 1)
        for i, weight in enumerate(cumulative_weights):
            if rand_val < weight:
                result.append(bases[i])
                break

    return "".join(result)


def mutate_sequence(seq: str, n_mut: int, *, rng: random.Random | None = None) -> str:
    r = rng or random
    if not seq or n_mut <= 0:
        return seq
    seq_list = list(seq)

    # Handle numpy random state vs standard random for sampling positions
    if hasattr(r, "choice") and hasattr(r, "seed") and not hasattr(r, "sample"):  # numpy RandomState
        # numpy doesn't have sample, use choice without replacement
        indices = list(range(len(seq_list)))
        positions = []
        for _ in range(min(n_mut, len(seq_list))):
            if not indices:
                break
            idx = r.choice(len(indices))
            positions.append(indices.pop(idx))
    else:
        # standard random module
        positions = r.sample(range(len(seq_list)), k=min(n_mut, len(seq_list)))

    for pos in positions:
        alphabet = _DNA if set(seq_list).issubset(set(_DNA)) else _AA
        choices = [c for c in alphabet if c != seq_list[pos]]
        if choices:  # Make sure there are alternatives
            seq_list[pos] = r.choice(choices)
    return "".join(seq_list)


def generate_random_protein(length: int, *, rng: random.Random | None = None) -> str:
    r = rng or random
    # Handle numpy random state vs standard random
    if hasattr(r, "choice") and hasattr(r, "seed"):  # numpy RandomState
        # numpy choice needs a sequence/array, not string
        aa_bases = list(_AA)
        return "".join(r.choice(aa_bases) for _ in range(max(0, length)))
    else:
        # standard random module
        return "".join(r.choice(_AA) for _ in range(max(0, length)))
