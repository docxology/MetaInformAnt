"""Community ecology analysis tools."""

from __future__ import annotations

import math
from collections import Counter
from typing import Dict, List, Sequence, Tuple

import numpy as np


def shannon_diversity(abundances: Sequence[float]) -> float:
    """Calculate Shannon diversity index H' = -Σ(pi * ln(pi)).

    Args:
        abundances: Species abundance values (counts or proportions)

    Returns:
        Shannon diversity index
    """
    if len(abundances) == 0:
        return 0.0

    # Convert to proportions if needed
    total = sum(abundances)
    if total == 0:
        return 0.0

    props = [x / total for x in abundances if x > 0]
    return -sum(p * math.log(p) for p in props)


def simpson_diversity(abundances: Sequence[float]) -> float:
    """Calculate Simpson's diversity index D = 1 - Σ(pi^2).

    Args:
        abundances: Species abundance values

    Returns:
        Simpson's diversity index
    """
    if len(abundances) == 0:
        return 0.0

    total = sum(abundances)
    if total == 0:
        return 0.0

    props = [x / total for x in abundances]
    simpson_index = sum(p * p for p in props)
    return 1.0 - simpson_index


def pielou_evenness(abundances: Sequence[float]) -> float:
    """Calculate Pielou's evenness index J' = H'/ln(S).

    Args:
        abundances: Species abundance values

    Returns:
        Pielou's evenness index
    """
    if len(abundances) == 0:
        return 0.0

    S = sum(1 for x in abundances if x > 0)  # Species richness
    if S <= 1:
        return 1.0 if S == 1 else 0.0

    H = shannon_diversity(abundances)
    return H / math.log(S)


def chao1_estimator(abundances: Sequence[int]) -> float:
    """Estimate total species richness using Chao1 estimator.

    Sest = Sobs + (f1^2)/(2*f2)
    where f1 = singletons, f2 = doubletons

    Args:
        abundances: Species abundance counts (integers)

    Returns:
        Estimated total species richness
    """
    if len(abundances) == 0:
        return 0.0

    counts = [int(x) for x in abundances if x > 0]
    S_obs = len(counts)

    f1 = sum(1 for x in counts if x == 1)  # Singletons
    f2 = sum(1 for x in counts if x == 2)  # Doubletons

    if f2 == 0:
        return S_obs + (f1 * (f1 - 1)) / 2.0

    return S_obs + (f1 * f1) / (2.0 * f2)


def bray_curtis_dissimilarity(site1: Sequence[float], site2: Sequence[float]) -> float:
    """Calculate Bray-Curtis dissimilarity between two sites.

    BC = Σ|x1i - x2i| / Σ(x1i + x2i)

    Args:
        site1, site2: Species abundance vectors for two sites

    Returns:
        Bray-Curtis dissimilarity (0 = identical, 1 = completely different)
    """
    if not site1 or not site2:
        return 1.0

    # Pad shorter vector with zeros
    max_len = max(len(site1), len(site2))
    s1 = list(site1) + [0] * (max_len - len(site1))
    s2 = list(site2) + [0] * (max_len - len(site2))

    numerator = sum(abs(a - b) for a, b in zip(s1, s2))
    denominator = sum(a + b for a, b in zip(s1, s2))

    return numerator / denominator if denominator > 0 else 0.0


def jaccard_similarity(site1: Sequence[float], site2: Sequence[float], presence_threshold: float = 0.0) -> float:
    """Calculate Jaccard similarity coefficient between two sites.

    J = |A ∩ B| / |A ∪ B|

    Args:
        site1, site2: Species abundance vectors
        presence_threshold: Minimum abundance to consider species present

    Returns:
        Jaccard similarity coefficient (0-1)
    """
    if not site1 or not site2:
        return 0.0

    max_len = max(len(site1), len(site2))
    s1 = list(site1) + [0] * (max_len - len(site1))
    s2 = list(site2) + [0] * (max_len - len(site2))

    present1 = set(i for i, x in enumerate(s1) if x > presence_threshold)
    present2 = set(i for i, x in enumerate(s2) if x > presence_threshold)

    intersection = len(present1 & present2)
    union = len(present1 | present2)

    return intersection / union if union > 0 else 0.0


def sorensen_similarity(site1: Sequence[float], site2: Sequence[float], presence_threshold: float = 0.0) -> float:
    """Calculate Sørensen similarity coefficient between two sites.

    S = 2|A ∩ B| / (|A| + |B|)

    Args:
        site1, site2: Species abundance vectors
        presence_threshold: Minimum abundance to consider species present

    Returns:
        Sørensen similarity coefficient (0-1)
    """
    if not site1 or not site2:
        return 0.0

    max_len = max(len(site1), len(site2))
    s1 = list(site1) + [0] * (max_len - len(site1))
    s2 = list(site2) + [0] * (max_len - len(site2))

    present1 = set(i for i, x in enumerate(s1) if x > presence_threshold)
    present2 = set(i for i, x in enumerate(s2) if x > presence_threshold)

    intersection = len(present1 & present2)

    return (2 * intersection) / (len(present1) + len(present2)) if (len(present1) + len(present2)) > 0 else 0.0


def rarefaction_curve(abundances: Sequence[int], max_samples: int | None = None) -> List[Tuple[int, float]]:
    """Generate rarefaction curve for species richness.

    Args:
        abundances: Species abundance counts
        max_samples: Maximum number of samples to rarify to (default: total abundance)

    Returns:
        List of (sample_size, expected_species_richness) tuples
    """
    if len(abundances) == 0:
        return []

    counts = [int(x) for x in abundances if x > 0]
    if not counts:
        return []

    total_individuals = sum(counts)
    if max_samples is None:
        max_samples = total_individuals

    curve = []
    for n in range(1, min(max_samples + 1, total_individuals + 1)):
        # Expected number of species in sample of size n
        expected_species = 0
        for ni in counts:
            # Probability that species i is NOT in sample of size n
            prob_absent = (
                math.comb(total_individuals - ni, n) / math.comb(total_individuals, n)
                if n <= total_individuals - ni
                else 0
            )
            expected_species += 1 - prob_absent

        curve.append((n, expected_species))

    return curve


def species_accumulation_curve(site_abundances: List[Sequence[float]]) -> List[Tuple[int, float]]:
    """Generate species accumulation curve across multiple sites.

    Args:
        site_abundances: List of abundance vectors, one per site

    Returns:
        List of (num_sites, cumulative_species_richness) tuples
    """
    if not site_abundances:
        return []

    curve = []
    cumulative_species = set()

    for i, site in enumerate(site_abundances, 1):
        # Add species present in this site
        for j, abundance in enumerate(site):
            if abundance > 0:
                cumulative_species.add(j)

        curve.append((i, len(cumulative_species)))

    return curve


def rank_abundance_distribution(abundances: Sequence[float]) -> List[Tuple[int, float]]:
    """Generate rank-abundance distribution (Whittaker plot data).

    Args:
        abundances: Species abundance values

    Returns:
        List of (rank, abundance) tuples, sorted by decreasing abundance
    """
    if len(abundances) == 0:
        return []

    sorted_abundances = sorted([x for x in abundances if x > 0], reverse=True)
    return [(i + 1, abundance) for i, abundance in enumerate(sorted_abundances)]


def beta_diversity_partitioning(site_abundances: List[Sequence[float]]) -> Dict[str, float]:
    """Partition beta diversity into turnover and nestedness components.

    Following Baselga (2010) framework.

    Args:
        site_abundances: List of abundance vectors, one per site

    Returns:
        Dictionary with beta diversity components
    """
    if len(site_abundances) < 2:
        return {"beta_total": 0.0, "beta_turnover": 0.0, "beta_nestedness": 0.0}

    # Convert to presence/absence
    presence_data = []
    for site in site_abundances:
        presence_data.append([1 if x > 0 else 0 for x in site])

    # Calculate pairwise components
    n_sites = len(presence_data)
    beta_components = []

    for i in range(n_sites):
        for j in range(i + 1, n_sites):
            site_i = presence_data[i]
            site_j = presence_data[j]

            # Ensure same length
            max_len = max(len(site_i), len(site_j))
            si = site_i + [0] * (max_len - len(site_i))
            sj = site_j + [0] * (max_len - len(site_j))

            a = sum(1 for x, y in zip(si, sj) if x == 1 and y == 1)  # Shared species
            b = sum(1 for x, y in zip(si, sj) if x == 1 and y == 0)  # Unique to i
            c = sum(1 for x, y in zip(si, sj) if x == 0 and y == 1)  # Unique to j

            if a + b + c > 0:
                beta_sim = min(b, c) / (a + min(b, c))  # Turnover component
                beta_sne = abs(b - c) / (2 * a + b + c) * a / (a + min(b, c)) if a + min(b, c) > 0 else 0  # Nestedness
                beta_sor = (b + c) / (2 * a + b + c)  # Total beta diversity

                beta_components.append({"turnover": beta_sim, "nestedness": beta_sne, "total": beta_sor})

    # Average across all pairwise comparisons
    if not beta_components:
        return {"beta_total": 0.0, "beta_turnover": 0.0, "beta_nestedness": 0.0}

    return {
        "beta_total": sum(x["total"] for x in beta_components) / len(beta_components),
        "beta_turnover": sum(x["turnover"] for x in beta_components) / len(beta_components),
        "beta_nestedness": sum(x["nestedness"] for x in beta_components) / len(beta_components),
    }


def functional_diversity_metrics(
    abundances: Sequence[float], trait_matrix: Sequence[Sequence[float]]
) -> Dict[str, float]:
    """Calculate functional diversity metrics.

    Args:
        abundances: Species abundance values
        trait_matrix: Matrix of trait values (species x traits)

    Returns:
        Dictionary with functional diversity metrics
    """
    if not abundances or not trait_matrix:
        return {"functional_richness": 0.0, "functional_evenness": 0.0, "functional_dispersion": 0.0}

    # Filter to present species
    present_species = [(i, ab) for i, ab in enumerate(abundances) if ab > 0 and i < len(trait_matrix)]
    if not present_species:
        return {"functional_richness": 0.0, "functional_evenness": 0.0, "functional_dispersion": 0.0}

    traits = [trait_matrix[i] for i, _ in present_species]
    weights = [ab / sum(x[1] for x in present_species) for _, ab in present_species]

    if not traits or not all(len(t) == len(traits[0]) for t in traits):
        return {"functional_richness": 0.0, "functional_evenness": 0.0, "functional_dispersion": 0.0}

    # Functional richness (convex hull volume - approximated by trait range)
    n_traits = len(traits[0])
    trait_ranges = []
    for j in range(n_traits):
        trait_values = [traits[i][j] for i in range(len(traits))]
        trait_ranges.append(max(trait_values) - min(trait_values) if len(set(trait_values)) > 1 else 0)

    functional_richness = math.prod(trait_ranges) if all(r > 0 for r in trait_ranges) else 0.0

    # Functional evenness (evenness of species in trait space)
    if len(traits) <= 1:
        functional_evenness = 1.0
    else:
        # Calculate pairwise distances in trait space
        distances = []
        for i in range(len(traits)):
            for j in range(i + 1, len(traits)):
                dist = math.sqrt(sum((traits[i][k] - traits[j][k]) ** 2 for k in range(n_traits)))
                distances.append(dist)

        mean_dist = sum(distances) / len(distances) if distances else 0
        functional_evenness = 1.0 - (2 / math.pi) * math.atan(np.var(distances) / mean_dist) if mean_dist > 0 else 1.0

    # Functional dispersion (abundance-weighted mean distance to centroid)
    if n_traits == 0:
        functional_dispersion = 0.0
    else:
        # Calculate abundance-weighted centroid
        centroid = [sum(weights[i] * traits[i][j] for i in range(len(traits))) for j in range(n_traits)]

        # Calculate distances to centroid
        distances_to_centroid = []
        for i in range(len(traits)):
            dist = math.sqrt(sum((traits[i][j] - centroid[j]) ** 2 for j in range(n_traits)))
            distances_to_centroid.append(dist * weights[i])

        functional_dispersion = sum(distances_to_centroid)

    return {
        "functional_richness": functional_richness,
        "functional_evenness": functional_evenness,
        "functional_dispersion": functional_dispersion,
    }
