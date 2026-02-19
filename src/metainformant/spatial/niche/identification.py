"""Tissue niche identification and characterization.

Identifies spatially coherent tissue niches (microenvironments) from
cell type composition data. Uses clustering on local cell type proportions
to define niches, then characterizes each niche by its composition,
spatial extent, and boundary properties.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class NicheResult:
    """Result of spatial niche identification.

    Attributes:
        niche_labels: Per-spot niche assignment.
        n_niches: Number of identified niches.
        compositions: Per-niche mean cell type composition (niche × cell_type).
        cell_type_names: Names of cell types.
        niche_sizes: Number of spots per niche.
        niche_diversity: Shannon diversity of cell types within each niche.
    """

    niche_labels: np.ndarray
    n_niches: int
    compositions: np.ndarray
    cell_type_names: list[str]
    niche_sizes: list[int]
    niche_diversity: list[float]


def identify_niches(
    cell_type_proportions: np.ndarray,
    coordinates: np.ndarray,
    cell_type_names: list[str] | None = None,
    n_niches: int = 5,
    n_neighbors: int = 15,
    spatial_weight: float = 0.3,
    random_state: int = 42,
) -> NicheResult:
    """Identify tissue niches from local cell type composition.

    1. Build a neighborhood graph from spatial coordinates.
    2. Smooth cell type proportions over the spatial neighborhood.
    3. Cluster smoothed proportions using K-Means to define niches.
    4. Characterize each niche by its mean composition and diversity.

    Args:
        cell_type_proportions: 2D array (spots × cell_types), rows sum to 1.
        coordinates: 2D array (spots × 2) of spatial coordinates.
        cell_type_names: Optional names for cell types.
        n_niches: Number of niches to identify.
        n_neighbors: Neighbors for spatial smoothing.
        spatial_weight: Weight for spatial smoothing (0 = no smoothing, 1 = full).
        random_state: Random seed.

    Returns:
        NicheResult with niche assignments and characterizations.
    """
    rng = np.random.RandomState(random_state)
    n_spots, n_types = cell_type_proportions.shape
    names = cell_type_names or [f"CellType_{i}" for i in range(n_types)]

    # Spatial smoothing
    smoothed = _spatial_smooth(
        cell_type_proportions, coordinates, n_neighbors, spatial_weight
    )

    # K-Means clustering on smoothed proportions
    labels = _kmeans(smoothed, n_niches, rng)

    # Characterize niches
    compositions = np.zeros((n_niches, n_types))
    niche_sizes = []
    niche_diversity = []

    for k in range(n_niches):
        mask = labels == k
        n_k = int(mask.sum())
        niche_sizes.append(n_k)

        if n_k > 0:
            mean_comp = cell_type_proportions[mask].mean(axis=0)
            compositions[k] = mean_comp
            # Shannon diversity
            p = mean_comp[mean_comp > 0]
            h = float(-np.sum(p * np.log(p)))
            niche_diversity.append(h)
        else:
            niche_diversity.append(0.0)

    return NicheResult(
        niche_labels=labels,
        n_niches=n_niches,
        compositions=compositions,
        cell_type_names=names,
        niche_sizes=niche_sizes,
        niche_diversity=niche_diversity,
    )


def _spatial_smooth(
    values: np.ndarray,
    coords: np.ndarray,
    k: int,
    weight: float,
) -> np.ndarray:
    """Smooth values over spatial neighbors.

    Args:
        values: 2D array (spots × features).
        coords: 2D array (spots × 2).
        k: Number of nearest neighbors.
        weight: Smoothing weight (0-1).

    Returns:
        Smoothed values.
    """
    n = len(values)
    k = min(k, n - 1)
    smoothed = values.copy()

    for i in range(n):
        dists = np.sqrt(np.sum((coords - coords[i]) ** 2, axis=1))
        dists[i] = np.inf
        nn_idx = np.argpartition(dists, k)[:k]
        neighbor_mean = values[nn_idx].mean(axis=0)
        smoothed[i] = (1 - weight) * values[i] + weight * neighbor_mean

    return smoothed


def _kmeans(
    data: np.ndarray,
    k: int,
    rng: np.random.RandomState,
    max_iter: int = 100,
) -> np.ndarray:
    """Simple K-Means clustering.

    Args:
        data: 2D array (samples × features).
        k: Number of clusters.
        rng: Random state.
        max_iter: Maximum iterations.

    Returns:
        1D array of cluster labels.
    """
    n = len(data)
    centers = data[rng.choice(n, size=k, replace=False)]
    labels = np.zeros(n, dtype=int)

    for _ in range(max_iter):
        # Assign
        dists = np.array([
            np.sqrt(np.sum((data - c) ** 2, axis=1)) for c in centers
        ])
        new_labels = np.argmin(dists, axis=0)

        if np.array_equal(new_labels, labels):
            break
        labels = new_labels

        # Update
        for j in range(k):
            mask = labels == j
            if mask.sum() > 0:
                centers[j] = data[mask].mean(axis=0)

    return labels
