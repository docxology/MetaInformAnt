from __future__ import annotations

import random
from typing import Iterable

from ..core import logging, validation

logger = logging.get_logger(__name__)

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
        ValidationError: If gc_content is not in [0, 1] or length is negative
    """
    validation.validate_type(length, int, "length")
    validation.validate_range(length, min_val=0, name="length")
    validation.validate_range(gc_content, min_val=0.0, max_val=1.0, name="gc_content")
    
    logger.debug(f"Generating DNA sequence: length={length}, gc_content={gc_content}")

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
    """Introduce random mutations into a sequence.
    
    Args:
        seq: DNA or protein sequence string
        n_mut: Number of mutations to introduce
        rng: Random number generator (default: random module)
        
    Returns:
        Mutated sequence with random substitutions
        
    Raises:
        ValidationError: If parameters are invalid
    """
    validation.validate_type(seq, str, "seq")
    validation.validate_type(n_mut, int, "n_mut")
    validation.validate_range(n_mut, min_val=0, name="n_mut")
    
    r = rng or random
    if not seq or n_mut <= 0:
        return seq
    
    logger.debug(f"Mutating sequence: length={len(seq)}, n_mutations={n_mut}")
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
    """Generate random protein sequence of specified length.
    
    Args:
        length: Desired sequence length
        rng: Random number generator (default: random module)
        
    Returns:
        Random protein sequence using standard 20 amino acids
        
    Raises:
        ValidationError: If length is negative
    """
    validation.validate_type(length, int, "length")
    validation.validate_range(length, min_val=0, name="length")
    
    logger.debug(f"Generating protein sequence: length={length}")
    r = rng or random
    # Handle numpy random state vs standard random
    if hasattr(r, "choice") and hasattr(r, "seed"):  # numpy RandomState
        # numpy choice needs a sequence/array, not string
        aa_bases = list(_AA)
        return "".join(r.choice(aa_bases) for _ in range(max(0, length)))
    else:
        # standard random module
        return "".join(r.choice(_AA) for _ in range(max(0, length)))
