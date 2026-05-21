"""Compatibility exports for population genetics math helpers."""

from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

from metainformant.math.population_genetics.core import hardy_weinberg_genotype_freqs
from metainformant.math.population_genetics.fst import fst_from_allele_freqs


def hardy_weinberg_allele_freqs(p: float, q: float | None = None) -> Tuple[float, float, float]:
    """Calculate Hardy-Weinberg genotype frequencies from allele frequencies."""
    if q is not None and abs((p + q) - 1.0) > 1e-6:
        raise ValueError("Allele frequencies must sum to 1")
    return hardy_weinberg_genotype_freqs(p)


def fst_from_freqs(freq1: Sequence[float], freq2: Sequence[float]) -> float:
    """Compatibility wrapper for FST from allele-frequency vectors."""
    return fst_from_allele_freqs(np.asarray(freq1, dtype=float), np.asarray(freq2, dtype=float))


__all__ = [
    "hardy_weinberg_genotype_freqs",
    "hardy_weinberg_allele_freqs",
    "fst_from_freqs",
]
