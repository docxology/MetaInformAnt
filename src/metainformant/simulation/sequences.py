from __future__ import annotations

import random
from typing import Iterable


_DNA = "ACGT"
_AA = "ACDEFGHIKLMNPQRSTVWY"


def generate_random_dna(length: int, *, rng: random.Random | None = None) -> str:
    r = rng or random
    return "".join(r.choice(_DNA) for _ in range(max(0, length)))


def mutate_sequence(seq: str, n_mut: int, *, rng: random.Random | None = None) -> str:
    r = rng or random
    if not seq or n_mut <= 0:
        return seq
    seq_list = list(seq)
    positions = r.sample(range(len(seq_list)), k=min(n_mut, len(seq_list)))
    for pos in positions:
        alphabet = _DNA if set(seq_list).issubset(set(_DNA)) else _AA
        choices = [c for c in alphabet if c != seq_list[pos]]
        seq_list[pos] = r.choice(choices)
    return "".join(seq_list)


def generate_random_protein(length: int, *, rng: random.Random | None = None) -> str:
    r = rng or random
    return "".join(r.choice(_AA) for _ in range(max(0, length)))


