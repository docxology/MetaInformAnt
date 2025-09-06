from __future__ import annotations

from typing import List

# Very simple helix propensity scale (normalized roughly 0..1)
HELIX_PROP = {
    "A": 1.00,
    "C": 0.40,
    "D": 0.35,
    "E": 0.85,
    "F": 0.55,
    "G": 0.10,
    "H": 0.50,
    "I": 0.50,
    "K": 0.80,
    "L": 0.90,
    "M": 0.75,
    "N": 0.35,
    "P": 0.05,
    "Q": 0.60,
    "R": 0.55,
    "S": 0.35,
    "T": 0.35,
    "V": 0.45,
    "W": 0.45,
    "Y": 0.40,
}


def simple_helix_coil_propensity(seq: str) -> List[float]:
    """Return per-residue helix propensity in [0,1] using a simple scale.

    This is a placeholder deterministic method that avoids external models,
    providing a quick baseline suitable for unit tests and examples.
    """
    out: List[float] = []
    for c in seq.upper():
        out.append(HELIX_PROP.get(c, 0.0))
    return out
