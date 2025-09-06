from __future__ import annotations

from typing import Iterable


def fst_from_heterozygosity(Hs: float, Ht: float) -> float:
    """Wright's Fst = (Ht - Hs) / Ht, clamped to [0, 1] when valid."""
    if Ht <= 0:
        return 0.0
    fst = (Ht - max(0.0, Hs)) / Ht
    if fst < 0:
        return 0.0
    if fst > 1:
        return 1.0
    return fst


def fst_from_allele_freqs(subpop_allele_freqs: Iterable[float]) -> float:
    """Compute Fst for a bi-allelic locus from subpopulation allele-A frequencies.

    Uses Hs = mean(2 p (1-p)) and Ht = 2 p_bar (1-p_bar).
    """
    ps = [p for p in subpop_allele_freqs if 0.0 <= p <= 1.0]
    if not ps:
        return 0.0
    pbar = sum(ps) / len(ps)
    Hs = sum(2.0 * p * (1.0 - p) for p in ps) / len(ps)
    Ht = 2.0 * pbar * (1.0 - pbar)
    return fst_from_heterozygosity(Hs, Ht)
