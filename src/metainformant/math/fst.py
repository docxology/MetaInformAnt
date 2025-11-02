from __future__ import annotations

from typing import Iterable


def fst_from_heterozygosity(Hs: float, Ht: float) -> float:
    """Calculate Wright's Fst from heterozygosity values.
    
    Fst measures population differentiation: the proportion of genetic
    diversity that is between subpopulations rather than within them.
    
    Args:
        Hs: Average heterozygosity within subpopulations
        Ht: Total heterozygosity in the entire population
        
    Returns:
        Fst value in [0, 1]. Returns 0.0 if Ht <= 0.
        Formula: Fst = (Ht - Hs) / Ht
        
    Examples:
        >>> fst_from_heterozygosity(Hs=0.3, Ht=0.4)
        0.25
        >>> fst_from_heterozygosity(Hs=0.4, Ht=0.4)
        0.0  # No differentiation
        
    References:
        Wright, S. (1949). The genetical structure of populations.
        Annals of Eugenics, 15(1), 323-354.
    """
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
