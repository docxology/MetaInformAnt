"""Functional trait diversity and community-weighted mean analysis.

Implements trait-based ecology metrics for understanding community assembly
rules, environmental filtering, and niche partitioning.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class TraitData:
    """Container for species trait data.

    Attributes:
        species: List of species names.
        trait_names: List of trait names.
        values: 2D array (species × traits) of trait values.
    """

    species: list[str]
    trait_names: list[str]
    values: np.ndarray


@dataclass
class TraitDiversityResult:
    """Result of functional trait diversity calculation.

    Attributes:
        fric: Functional richness (convex hull volume in trait space).
        feve: Functional evenness (regularity of spacing in trait space).
        fdiv: Functional divergence (deviation from centroid).
        cwm: Community-weighted mean for each trait.
        rao_q: Rao's quadratic entropy (abundance-weighted trait dissimilarity).
    """

    fric: float
    feve: float
    fdiv: float
    cwm: dict[str, float]
    rao_q: float


def community_weighted_mean(
    abundances: np.ndarray,
    trait_values: np.ndarray,
    trait_names: list[str] | None = None,
) -> dict[str, float]:
    """Compute community-weighted mean (CWM) for each trait.

    CWM_j = Σ_i (p_i × t_ij) where p_i is relative abundance of species i
    and t_ij is the trait value j of species i.

    Args:
        abundances: 1D array of species abundances (will be normalized).
        trait_values: 2D array (species × traits).
        trait_names: Optional names for traits; defaults to Trait_0, Trait_1, ...

    Returns:
        Dictionary mapping trait name to community-weighted mean value.
    """
    if abundances.ndim != 1:
        raise ValueError("abundances must be a 1D array")
    if trait_values.ndim != 2:
        raise ValueError("trait_values must be a 2D array (species × traits)")
    if len(abundances) != trait_values.shape[0]:
        raise ValueError("Number of species in abundances and trait_values must match")

    rel_abund = abundances / abundances.sum()
    n_traits = trait_values.shape[1]
    names = trait_names or [f"Trait_{i}" for i in range(n_traits)]

    cwm_values = rel_abund @ trait_values
    return dict(zip(names, cwm_values.tolist()))


def functional_richness(trait_values: np.ndarray) -> float:
    """Compute functional richness as convex hull volume in trait space.

    Uses the quick-hull algorithm via scipy. Falls back to range for 1D traits.

    Args:
        trait_values: 2D array (species × traits), at least 3 species required.

    Returns:
        Convex hull volume (or area in 2D, length in 1D).
    """
    n_species, n_traits = trait_values.shape
    if n_species < n_traits + 1:
        return 0.0

    if n_traits == 1:
        return float(trait_values.max() - trait_values.min())

    try:
        from scipy.spatial import ConvexHull

        hull = ConvexHull(trait_values)
        return float(hull.volume)
    except Exception:
        return 0.0


def rao_quadratic_entropy(
    abundances: np.ndarray,
    trait_values: np.ndarray,
) -> float:
    """Compute Rao's quadratic entropy (Q).

    Q = Σ_i Σ_j d_ij × p_i × p_j

    where d_ij is the Euclidean distance between species i and j in trait space,
    and p_i is the relative abundance of species i.

    Args:
        abundances: 1D array of species abundances.
        trait_values: 2D array (species × traits).

    Returns:
        Rao's Q value.
    """
    rel_abund = abundances / abundances.sum()
    n = len(rel_abund)

    q = 0.0
    for i in range(n):
        for j in range(n):
            d_ij = float(np.sqrt(np.sum((trait_values[i] - trait_values[j]) ** 2)))
            q += d_ij * rel_abund[i] * rel_abund[j]
    return q


def trait_diversity(
    abundances: np.ndarray,
    trait_data: TraitData,
) -> TraitDiversityResult:
    """Compute comprehensive functional trait diversity metrics.

    Calculates FRic, FEve, FDiv, CWM, and Rao's Q in one call.

    Args:
        abundances: 1D array of species abundances.
        trait_data: TraitData container with species, trait_names, and values.

    Returns:
        TraitDiversityResult with all computed metrics.
    """
    cwm = community_weighted_mean(abundances, trait_data.values, trait_data.trait_names)
    fric = functional_richness(trait_data.values)
    rao_q = rao_quadratic_entropy(abundances, trait_data.values)

    # Functional evenness: regularity of spacing along minimum spanning tree
    rel_abund = abundances / abundances.sum()
    n = len(rel_abund)
    if n < 3:
        feve = 0.0
    else:
        dists = np.sqrt(
            np.sum((trait_data.values[:, None, :] - trait_data.values[None, :, :]) ** 2, axis=2)
        )
        sorted_idx = np.argsort(rel_abund)[::-1]
        ew = np.array([dists[sorted_idx[i], sorted_idx[i + 1]] for i in range(n - 1)])
        if ew.sum() > 0:
            pew = ew / ew.sum()
            feve = float((1.0 / (n - 1)) * np.sum(np.minimum(pew, 1.0 / (n - 1))))
        else:
            feve = 0.0

    # Functional divergence: deviation from centroid weighted by abundance
    centroid = rel_abund @ trait_data.values
    deviations = np.sqrt(np.sum((trait_data.values - centroid) ** 2, axis=1))
    mean_dev = float(np.mean(deviations))
    if mean_dev > 0:
        weighted_dev = float(rel_abund @ deviations)
        fdiv = weighted_dev / mean_dev
    else:
        fdiv = 0.0

    return TraitDiversityResult(
        fric=fric,
        feve=feve,
        fdiv=fdiv,
        cwm=cwm,
        rao_q=rao_q,
    )
