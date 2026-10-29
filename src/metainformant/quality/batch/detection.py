"""Batch effect detection and correction for sequencing data.

Implements statistical methods to identify batch effects, assess their
magnitude, and correct them with parametric empirical Bayes ComBat
(Johnson, Li & Rabinovic 2007, Biostatistics 8(1):118-127).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import f as f_dist

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Convergence threshold and safety cap of the empirical Bayes posterior iteration in ComBat.
_COMBAT_TOLERANCE = 1e-10
_COMBAT_MAX_ITERATIONS = 1000


@dataclass
class BatchEffectReport:
    """Summary of detected batch effects.

    Attributes:
        n_samples: Total number of samples.
        n_batches: Number of distinct batches.
        batch_labels: List of batch labels per sample.
        pvca_variance: Dict mapping source (batch, residual) to variance proportion: the between-batch and
            within-batch sums of squares of all features pooled, each as a share of their total.
        silhouette_score: Score measuring batch separation in PCA space (-1 to 1).
        n_significant_features: Number of features whose one-way F-test against the batch labels has a raw
            (uncorrected for multiple testing) p-value below ``alpha``.
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
    """Detect batch effects via a PVCA-style variance split, per-feature F-tests and a silhouette score.

    Three quantities are computed:

    * the variance split ``pvca_variance``: the between-batch sum of squares of all features, pooled,
      as a share of the total sum of squares. This is a PVCA-style summary, not Li et al.'s PVCA: it
      works on the features themselves, not on principal components, and fits no mixed model, so
      features with a large variance dominate it;
    * one one-way F-test per feature against the batch labels, whose exact p-values come from the
      survival function of the F distribution, so a very small p-value is not rounded to zero. A feature
      that does not vary within a batch but differs between batches is perfectly explained by the batch
      (F is infinite, p = 0); a feature that does not vary at all has no association (p = 1);
    * the silhouette score of the batches in the space of the first (up to 10) principal components.

    ``severity`` is "critical" when the batch share of the variance exceeds 0.3 or the silhouette exceeds 0.5,
    otherwise "high" above 0.15 / 0.3, "moderate" above 0.05 / 0.1, and "low" below.

    Args:
        data: 2D array (samples × features), finite.
        batch_labels: List of batch identifiers, one per sample.
        alpha: Significance threshold for per-feature batch association test, in (0, 1). The p-values are
            raw: they are not corrected for the number of features.

    Returns:
        BatchEffectReport with severity assessment. With a single batch there is nothing to detect: the
        batch share of the variance is 0 and no feature is significant.

    Raises:
        ValueError: If ``data`` is not 2D, has no samples or contains NaN/inf (a missing value would
            otherwise make its feature look unassociated with the batch), ``batch_labels`` does not have
            one entry per sample, ``alpha`` is not in (0, 1), or every sample is a batch of its own (no
            feature has within-batch variation to test against; the labels may be sample identifiers).
    """
    values = np.asarray(data, dtype=float)
    if values.ndim != 2:
        raise ValueError(f"data must be a 2D samples x features array, got {values.ndim}D")
    n_samples, n_features = values.shape
    if n_samples == 0:
        raise ValueError("data has no samples")
    if len(batch_labels) != n_samples:
        raise ValueError(f"batch_labels has {len(batch_labels)} entries for {n_samples} samples")
    if not np.all(np.isfinite(values)):
        raise ValueError("data contains NaN or infinite values; impute or remove them before batch detection")
    if not 0.0 < alpha < 1.0:  # False for NaN as well
        raise ValueError(f"alpha must be in (0, 1), got {alpha!r}")
    data = values

    unique_batches = sorted(set(batch_labels))
    n_batches = len(unique_batches)
    if n_batches > 1 and n_batches == n_samples:
        raise ValueError(
            f"all {n_samples} samples are in batches of their own, so there is no within-batch variation to "
            "test against; check that batch_labels are batch identifiers and not sample identifiers"
        )
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

    if n_batches > 1:
        ms_between = ss_between / (n_batches - 1)
        ms_within = ss_within / (n_samples - n_batches)
        varies_within = ms_within > 0
        f_stats = np.zeros(n_features)
        np.divide(ms_between, ms_within, out=f_stats, where=varies_within)
        # constant within every batch yet different between batches: the batch explains the feature exactly
        f_stats[~varies_within & (ms_between > 0)] = np.inf
        # Exact p-values from the F-distribution. sf, not 1 - cdf: the cdf rounds to 1.0 for every p-value
        # below about 1e-16, so 1 - cdf reports 0.0 for p-values of 1e-17 and of 1e-300 alike.
        p_values = f_dist.sf(f_stats, n_batches - 1, n_samples - n_batches)
    else:
        p_values = np.ones(n_features)  # one batch: nothing to be associated with

    n_significant = int((p_values < alpha).sum())

    # PVCA-style variance decomposition
    total_ss = ss_between.sum() + ss_within.sum()
    batch_var_prop = float(ss_between.sum() / total_ss) if total_ss > 0 else 0.0
    residual_var_prop = 1.0 - batch_var_prop

    # PCA-based silhouette score (a failure of the SVD propagates: a silhouette in some other space would
    # not be the PCA-space score the report promises)
    centered = data - data.mean(axis=0)
    u, s, _ = np.linalg.svd(centered, full_matrices=False)
    n_pcs = min(10, n_features, n_samples)
    pca_coords = u[:, :n_pcs] * s[:n_pcs]

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
    """Remove batch effects with parametric empirical Bayes ComBat.

    Implements Johnson, Li & Rabinovic (2007), without covariates. For sample ``j`` of batch ``i`` and
    feature ``g`` the model is ``Y_ijg = alpha_g + gamma_ig + delta_ig * sigma_g * e_ijg``: a
    batch-specific location ``gamma_ig`` and scale ``delta_ig`` per feature. The steps are:

    1. Standardise once, with statistics of the whole data set: ``Z_ijg = (Y_ijg - alpha_g) / sigma_g``,
       where ``alpha_g`` is the grand mean and ``sigma_g**2`` the pooled within-batch variance (divided by
       the number of samples, as in the paper and in ``sva``).
    2. Per batch, estimate the batch effects by the mean ``gamma_hat_ig`` and the sample variance
       ``delta_hat_ig**2`` of ``Z``.
    3. Per batch, borrow strength across features: the method of moments fits a normal prior
       ``N(gamma_bar_i, tau_i**2)`` to ``gamma_hat_ig`` and an inverse-gamma prior to ``delta_hat_ig**2``
       (hyperparameters ``lambda_i`` and ``theta_i``).
    4. Iterate the two posterior means to convergence:
       ``gamma*_ig = (n_i tau_i**2 gamma_hat_ig + delta*_ig**2 gamma_bar_i) / (n_i tau_i**2 + delta*_ig**2)``
       and ``delta*_ig**2 = (theta_i + sum_j (Z_ijg - gamma*_ig)**2 / 2) / (n_i / 2 + lambda_i - 1)``.
    5. Adjust: ``Y*_ijg = sigma_g / delta*_ig * (Z_ijg - gamma*_ig) + alpha_g``.

    Every batch is estimated independently from the same once-standardised data, so the result does not
    depend on the order or the labels of the batches, nor on the order of the samples. The amount of
    shrinkage is not a constant: it follows from the spread of ``gamma_hat`` and ``delta_hat**2`` across
    features and from the batch size.

    Features that are constant within some batch have no estimable batch variance and are returned
    unchanged (the behaviour of ``sva``). If there is a single batch there is nothing to correct and a
    copy of ``data`` is returned.

    Args:
        data: 2D array (samples × features), finite.
        batch_labels: List of batch identifiers, one per sample.

    Returns:
        Corrected data array of the same shape; ``data`` itself is not modified.

    Raises:
        ValueError: If ``data`` is not 2D or contains NaN/inf, ``batch_labels`` does not match the
            number of samples, a batch has a single sample (its variance cannot be estimated, so drop
            or merge it), or exactly one feature is estimable (the empirical Bayes priors need at
            least two features).
    """
    values = np.asarray(data, dtype=float)
    if values.ndim != 2:
        raise ValueError(f"data must be a 2D samples x features array, got {values.ndim}D")
    n_samples = values.shape[0]
    if len(batch_labels) != n_samples:
        raise ValueError(f"batch_labels has {len(batch_labels)} entries for {n_samples} samples")
    if not np.all(np.isfinite(values)):
        raise ValueError("data contains NaN or infinite values; impute or remove them before batch correction")

    members: dict[str, list[int]] = {}
    for index, label in enumerate(batch_labels):
        members.setdefault(label, []).append(index)
    batches = [np.array(indices) for indices in members.values()]

    corrected = values.copy()
    if len(batches) < 2:
        return corrected
    if min(len(batch) for batch in batches) < 2:
        raise ValueError(
            "every batch needs at least 2 samples to estimate its variance; drop or merge the batches "
            f"with one sample ({[label for label, idx in members.items() if len(idx) < 2]})"
        )

    # Without within-batch variation in some batch the scale of a feature cannot be estimated
    estimable = ~np.any([np.ptp(values[batch], axis=0) == 0 for batch in batches], axis=0)
    n_estimable = int(estimable.sum())
    if n_estimable == 0:
        logger.warning("ComBat: every feature is constant within some batch; returning the data unchanged")
        return corrected
    if n_estimable == 1:
        raise ValueError(
            "ComBat needs at least 2 features that vary within every batch to fit its empirical Bayes priors"
        )

    # Step 1: standardise once
    x = values[:, estimable]
    grand_mean = x.mean(axis=0)
    batch_mean = np.empty_like(x)
    for batch in batches:
        batch_mean[batch] = x[batch].mean(axis=0)
    pooled_sd = np.sqrt(((x - batch_mean) ** 2).sum(axis=0) / n_samples)
    z = (x - grand_mean) / pooled_sd

    # Steps 2-5, independently per batch
    adjusted = np.empty_like(z)
    for batch in batches:
        gamma_star, delta_star_sq = _combat_posterior(z[batch])
        adjusted[batch] = (z[batch] - gamma_star) / np.sqrt(delta_star_sq)

    corrected[:, estimable] = adjusted * pooled_sd + grand_mean
    return corrected


def _combat_posterior(z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Empirical Bayes posterior means of one batch's location and scale effects (ComBat steps 2-4).

    Args:
        z: Standardised data of one batch, shape (batch samples, features), at least 2 samples and
            2 features, no feature constant.

    Returns:
        ``(gamma_star, delta_star_sq)``: per-feature posterior mean of the batch location effect and of
        the batch variance effect.
    """
    n = z.shape[0]
    gamma_hat = z.mean(axis=0)
    delta_hat_sq = z.var(axis=0, ddof=1)

    # Hyperparameters by the method of moments, over features
    gamma_bar = gamma_hat.mean()
    tau_sq = gamma_hat.var(ddof=1)
    delta_mean = delta_hat_sq.mean()
    delta_var = delta_hat_sq.var(ddof=1)
    lam = theta = 0.0
    if delta_var > 0:
        lam = (2.0 * delta_var + delta_mean**2) / delta_var
        theta = (delta_mean * delta_var + delta_mean**3) / delta_var

    gamma_star = gamma_hat
    delta_star_sq = delta_hat_sq
    for _ in range(_COMBAT_MAX_ITERATIONS):
        new_gamma = (n * tau_sq * gamma_hat + delta_star_sq * gamma_bar) / (n * tau_sq + delta_star_sq)
        if delta_var > 0:
            new_delta = (theta + 0.5 * ((z - new_gamma) ** 2).sum(axis=0)) / (n / 2.0 + lam - 1.0)
        else:
            # all features share one variance: the inverse-gamma prior degenerates to a point mass
            new_delta = np.full_like(delta_star_sq, delta_mean)
        change = max(
            float(np.max(np.abs(new_gamma - gamma_star))),
            float(np.max(np.abs(new_delta - delta_star_sq) / delta_star_sq)),
        )
        gamma_star, delta_star_sq = new_gamma, new_delta
        if change < _COMBAT_TOLERANCE:
            break
    else:
        logger.warning("ComBat: posterior iteration did not converge in %d steps", _COMBAT_MAX_ITERATIONS)

    return gamma_star, delta_star_sq


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
        a_i = float(np.sqrt(np.sum((coords[own_mask] - coords[i]) ** 2, axis=1)).sum()) / (own_count - 1)

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
