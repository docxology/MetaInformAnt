"""Batch effect detection and correction for sequencing data.

Implements statistical methods to identify batch effects, assess their
magnitude, and apply corrections using empirical Bayes (ComBat-like)
and linear model approaches.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass
class BatchEffectReport:
    """Summary of detected batch effects.

    Attributes:
        n_samples: Total number of samples.
        n_batches: Number of distinct batches.
        batch_labels: List of batch labels per sample.
        pvca_variance: Dict mapping source (batch, residual) to variance proportion.
        silhouette_score: Score measuring batch separation in PCA space (-1 to 1).
        n_significant_features: Number of features significantly associated with batch.
        severity: Overall severity (low, moderate, high, critical).
    """

    n_samples: int
    n_batches: int
    batch_labels: list[str]
    pvca_variance: dict[str, float]
    silhouette_score: float
    n_significant_features: int
    severity: str


def detect_batch_effects(
    data: np.ndarray,
    batch_labels: list[str],
    alpha: float = 0.05,
) -> BatchEffectReport:
    """Detect batch effects via PVCA and silhouette analysis.

    Performs Principal Variance Component Analysis (PVCA) to estimate the
    proportion of variance attributable to batch, and computes a silhouette
    score in PCA space to quantify batch separation.

    Args:
        data: 2D array (samples × features).
        batch_labels: List of batch identifiers, one per sample.
        alpha: Significance threshold for per-feature batch association test.

    Returns:
        BatchEffectReport with severity assessment.
    """
    n_samples, n_features = data.shape
    unique_batches = sorted(set(batch_labels))
    n_batches = len(unique_batches)
    batch_map = {b: i for i, b in enumerate(unique_batches)}
    batch_idx = np.array([batch_map[b] for b in batch_labels])

    # Per-feature one-way ANOVA-like F-test for batch association
    grand_mean = data.mean(axis=0)
    ss_between = np.zeros(n_features)
    ss_within = np.zeros(n_features)
    for b in range(n_batches):
        mask = batch_idx == b
        n_b = mask.sum()
        if n_b == 0:
            continue
        batch_mean = data[mask].mean(axis=0)
        ss_between += n_b * (batch_mean - grand_mean) ** 2
        ss_within += ((data[mask] - batch_mean) ** 2).sum(axis=0)

    df_between = max(n_batches - 1, 1)
    df_within = max(n_samples - n_batches, 1)
    ms_between = ss_between / df_between
    ms_within = ss_within / df_within
    f_stats = np.where(ms_within > 0, ms_between / ms_within, 0.0)

    # Approximate p-values using F-distribution
    try:
        from scipy.stats import f as f_dist

        p_values = 1.0 - f_dist.cdf(f_stats, df_between, df_within)
    except ImportError:
        # Rough threshold: F > 3.0 considered significant
        p_values = np.where(f_stats > 3.0, 0.01, 0.5)

    n_significant = int((p_values < alpha).sum())

    # PVCA-like variance decomposition
    total_ss = ss_between.sum() + ss_within.sum()
    batch_var_prop = float(ss_between.sum() / total_ss) if total_ss > 0 else 0.0
    residual_var_prop = 1.0 - batch_var_prop

    # PCA-based silhouette score
    centered = data - data.mean(axis=0)
    try:
        u, s, vt = np.linalg.svd(centered, full_matrices=False)
        n_pcs = min(10, n_features, n_samples)
        pca_coords = u[:, :n_pcs] * s[:n_pcs]
    except np.linalg.LinAlgError:
        pca_coords = centered[:, :min(10, n_features)]

    # Simplified silhouette
    silhouette = _compute_silhouette(pca_coords, batch_idx, n_batches)

    # Severity classification
    if batch_var_prop > 0.3 or silhouette > 0.5:
        severity = "critical"
    elif batch_var_prop > 0.15 or silhouette > 0.3:
        severity = "high"
    elif batch_var_prop > 0.05 or silhouette > 0.1:
        severity = "moderate"
    else:
        severity = "low"

    return BatchEffectReport(
        n_samples=n_samples,
        n_batches=n_batches,
        batch_labels=batch_labels,
        pvca_variance={"batch": batch_var_prop, "residual": residual_var_prop},
        silhouette_score=silhouette,
        n_significant_features=n_significant,
        severity=severity,
    )


def correct_batch_combat(
    data: np.ndarray,
    batch_labels: list[str],
) -> np.ndarray:
    """Apply ComBat-like empirical Bayes batch correction.

    Centers each batch to the grand mean and shrinks batch-specific
    variance estimates using an empirical Bayes prior.

    Args:
        data: 2D array (samples × features).
        batch_labels: List of batch identifiers, one per sample.

    Returns:
        Corrected data array with batch effects removed.
    """
    corrected = data.copy().astype(float)
    unique_batches = sorted(set(batch_labels))
    batch_map = {b: i for i, b in enumerate(unique_batches)}
    batch_idx = np.array([batch_map[b] for b in batch_labels])

    grand_mean = data.mean(axis=0)

    for b in range(len(unique_batches)):
        mask = batch_idx == b
        if mask.sum() == 0:
            continue
        batch_mean = corrected[mask].mean(axis=0)
        batch_std = corrected[mask].std(axis=0)
        grand_std = corrected.std(axis=0)

        # Shift to grand mean
        corrected[mask] -= (batch_mean - grand_mean)

        # Scale variance (empirical Bayes shrinkage toward grand variance)
        safe_batch_std = np.where(batch_std > 0, batch_std, 1.0)
        safe_grand_std = np.where(grand_std > 0, grand_std, 1.0)
        scale = safe_grand_std / safe_batch_std

        # Shrink scale toward 1.0 (prior)
        shrinkage = 0.5
        scale = shrinkage * scale + (1.0 - shrinkage) * 1.0

        corrected[mask] = (corrected[mask] - grand_mean) * scale + grand_mean

    return corrected


def _compute_silhouette(
    coords: np.ndarray,
    labels: np.ndarray,
    n_groups: int,
) -> float:
    """Compute mean silhouette score for cluster separation.

    Args:
        coords: 2D array of coordinates.
        labels: Integer cluster labels.
        n_groups: Number of distinct groups.

    Returns:
        Mean silhouette score (-1 to 1).
    """
    n = len(labels)
    if n < 2 or n_groups < 2:
        return 0.0

    silhouettes = []
    for i in range(n):
        own_mask = labels == labels[i]
        own_count = own_mask.sum()
        if own_count <= 1:
            silhouettes.append(0.0)
            continue

        # Mean intra-cluster distance
        a_i = float(np.sqrt(np.sum((coords[own_mask] - coords[i]) ** 2, axis=1)).sum()) / (
            own_count - 1
        )

        # Mean nearest-cluster distance
        b_i = float("inf")
        for g in range(n_groups):
            if g == labels[i]:
                continue
            g_mask = labels == g
            if g_mask.sum() == 0:
                continue
            dist_g = float(np.sqrt(np.sum((coords[g_mask] - coords[i]) ** 2, axis=1)).mean())
            b_i = min(b_i, dist_g)

        if b_i == float("inf"):
            silhouettes.append(0.0)
        else:
            s_i = (b_i - a_i) / max(a_i, b_i) if max(a_i, b_i) > 0 else 0.0
            silhouettes.append(s_i)

    return float(np.mean(silhouettes))
