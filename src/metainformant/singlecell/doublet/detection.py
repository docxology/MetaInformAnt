"""Doublet detection for single-cell RNA-seq data.

Implements a simulation-based doublet detection approach inspired by
Scrublet and DoubletFinder. Generates synthetic doublets from observed
expression profiles, then scores each observed cell by its proximity
to synthetic doublets in a low-dimensional embedding.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class DoubletResult:
    """Result of doublet detection.

    Attributes:
        scores: Per-cell doublet score (0 to 1, higher = more likely doublet).
        predicted_doublets: Boolean mask of predicted doublets.
        threshold: Score threshold used for classification.
        n_doublets: Number of predicted doublets.
        doublet_rate: Fraction of cells predicted as doublets.
        synthetic_scores: Scores for synthetic doublets (for calibration).
    """

    scores: np.ndarray
    predicted_doublets: np.ndarray
    threshold: float
    n_doublets: int
    doublet_rate: float
    synthetic_scores: np.ndarray


def detect_doublets(
    counts: np.ndarray,
    expected_doublet_rate: float = 0.06,
    n_synthetic: int | None = None,
    n_neighbors: int = 30,
    n_pcs: int = 30,
    random_state: int = 42,
) -> DoubletResult:
    """Detect doublets using simulation-based scoring.

    1. Normalize and select highly variable genes.
    2. Generate synthetic doublets by averaging random cell pairs.
    3. Project observed + synthetic cells into PCA space.
    4. Score each observed cell by the fraction of its KNN neighbors
       that are synthetic doublets.

    Args:
        counts: 2D array (cells × genes) of raw counts.
        expected_doublet_rate: Expected fraction of doublets (for threshold calibration).
        n_synthetic: Number of synthetic doublets to generate. Defaults to n_cells.
        n_neighbors: Number of nearest neighbors for KNN scoring.
        n_pcs: Number of principal components for embedding.
        random_state: Random seed for reproducibility.

    Returns:
        DoubletResult with scores and predicted doublets.
    """
    rng = np.random.RandomState(random_state)
    n_cells, n_genes = counts.shape

    if n_synthetic is None:
        n_synthetic = n_cells

    # Normalize (library size + log1p)
    lib_sizes = counts.sum(axis=1, keepdims=True)
    lib_sizes = np.where(lib_sizes > 0, lib_sizes, 1.0)
    normed = np.log1p(counts / lib_sizes * 1e4)

    # Select top variable genes
    gene_var = normed.var(axis=0)
    gene_mean = normed.mean(axis=0)
    # Coefficient of variation-based selection
    cv = np.where(gene_mean > 0, np.sqrt(gene_var) / gene_mean, 0.0)
    n_hvg = min(2000, n_genes)
    hvg_idx = np.argsort(cv)[-n_hvg:]
    normed_hvg = normed[:, hvg_idx]

    # Generate synthetic doublets
    idx1 = rng.randint(0, n_cells, size=n_synthetic)
    idx2 = rng.randint(0, n_cells, size=n_synthetic)
    synthetic = (normed_hvg[idx1] + normed_hvg[idx2]) / 2.0

    # Combine observed and synthetic
    combined = np.vstack([normed_hvg, synthetic])
    is_synthetic = np.zeros(len(combined), dtype=bool)
    is_synthetic[n_cells:] = True

    # PCA
    centered = combined - combined.mean(axis=0)
    n_components = min(n_pcs, combined.shape[0] - 1, combined.shape[1])
    try:
        u, s, vt = np.linalg.svd(centered, full_matrices=False)
        pca_coords = u[:, :n_components] * s[:n_components]
    except np.linalg.LinAlgError:
        pca_coords = centered[:, :n_components]

    # KNN doublet scoring
    k = min(n_neighbors, len(combined) - 1)
    scores = _knn_doublet_scores(pca_coords, is_synthetic, k, n_cells)

    # Synthetic scores for calibration
    synthetic_scores = _knn_doublet_scores(
        pca_coords, is_synthetic, k, n_cells, score_synthetic=True
    )

    # Threshold: find score that yields expected doublet rate
    sorted_scores = np.sort(scores)
    threshold_idx = max(0, int((1.0 - expected_doublet_rate) * n_cells) - 1)
    threshold = float(sorted_scores[threshold_idx])

    predicted = scores > threshold
    n_doublets = int(predicted.sum())

    return DoubletResult(
        scores=scores,
        predicted_doublets=predicted,
        threshold=threshold,
        n_doublets=n_doublets,
        doublet_rate=n_doublets / n_cells if n_cells > 0 else 0.0,
        synthetic_scores=synthetic_scores,
    )


def _knn_doublet_scores(
    coords: np.ndarray,
    is_synthetic: np.ndarray,
    k: int,
    n_observed: int,
    score_synthetic: bool = False,
) -> np.ndarray:
    """Compute doublet scores using brute-force KNN.

    For each target cell, computes fraction of K nearest neighbors
    that are synthetic doublets.

    Args:
        coords: 2D PCA coordinates for all cells (observed + synthetic).
        is_synthetic: Boolean mask indicating synthetic cells.
        k: Number of nearest neighbors.
        n_observed: Number of observed (non-synthetic) cells.
        score_synthetic: If True, score synthetic cells; otherwise score observed.

    Returns:
        1D array of doublet scores for observed (or synthetic) cells.
    """
    if score_synthetic:
        target_range = range(n_observed, len(coords))
    else:
        target_range = range(n_observed)

    scores = np.zeros(len(target_range))
    for idx, i in enumerate(target_range):
        dists = np.sqrt(np.sum((coords - coords[i]) ** 2, axis=1))
        dists[i] = np.inf  # exclude self
        nn_idx = np.argpartition(dists, k)[:k]
        scores[idx] = is_synthetic[nn_idx].sum() / k

    return scores
