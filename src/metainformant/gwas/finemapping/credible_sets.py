"""Statistical fine-mapping: credible sets, SuSiE, Bayes factors, and colocalization.

Implements core fine-mapping algorithms for identifying causal variants from
GWAS summary statistics. Includes approximate Bayes factor credible sets,
a pure-Python SuSiE (Sum of Single Effects) regression, Wakefield approximate
Bayes factors, coloc-style colocalization, stepwise conditional analysis,
and functional annotation enrichment of credible sets.

References:
    - Wakefield (2009) Genetic Epidemiology 33:79-86 (ABF)
    - Wang et al. (2020) JRSS-B 82:1273-1300 (SuSiE)
    - Giambartolomei et al. (2014) PLoS Genetics 10:e1004383 (coloc)
    - Yang et al. (2012) Nature Genetics 44:369-375 (conditional analysis)
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_credible_set(
    z_scores: list[float],
    ld_matrix: list[list[float]] | None = None,
    coverage: float = 0.95,
    prior: str = "uniform",
) -> dict:
    """Compute an approximate Bayes factor credible set from Z-scores.

    Calculates posterior inclusion probabilities (PIPs) from GWAS Z-scores
    using Wakefield approximate Bayes factors, optionally adjusting for
    linkage disequilibrium. Returns the smallest set of variants whose
    cumulative PIP exceeds the requested coverage.

    Args:
        z_scores: Z-scores for each variant in the region.
        ld_matrix: Optional LD (r-squared) matrix. If provided, used to
            weight Bayes factors by LD structure. Square matrix of size
            len(z_scores) x len(z_scores).
        coverage: Target cumulative posterior probability for the credible
            set (default 0.95).
        prior: Prior type on causal status. One of "uniform" (equal prior
            per variant) or "distance" (decaying prior from lead variant).

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - snps_in_set: List of indices in the credible set
            - pips: List of posterior inclusion probabilities for all variants
            - coverage_achieved: Actual cumulative PIP of the credible set
            - n_snps: Number of variants in the credible set
            - bayes_factors: List of approximate Bayes factors
    """
    if not z_scores:
        return {"status": "error", "message": "No Z-scores provided"}

    n_variants = len(z_scores)

    if coverage <= 0.0 or coverage > 1.0:
        return {"status": "error", "message": f"Coverage must be in (0, 1], got {coverage}"}

    if prior not in ("uniform", "distance"):
        return {"status": "error", "message": f"Unknown prior type: {prior}"}

    # Validate LD matrix dimensions if provided
    if ld_matrix is not None:
        if len(ld_matrix) != n_variants:
            return {
                "status": "error",
                "message": (f"LD matrix rows ({len(ld_matrix)}) does not match " f"number of Z-scores ({n_variants})"),
            }
        for row_idx, row in enumerate(ld_matrix):
            if len(row) != n_variants:
                return {
                    "status": "error",
                    "message": (f"LD matrix row {row_idx} has {len(row)} columns, " f"expected {n_variants}"),
                }

    logger.debug(f"Computing credible set for {n_variants} variants, " f"coverage={coverage}, prior={prior}")

    # Step 1: Compute approximate Bayes factors
    bayes_factors = compute_bayes_factors(z_scores)

    # Step 2: Apply LD adjustment if matrix provided
    if ld_matrix is not None and HAS_NUMPY:
        bayes_factors = _apply_ld_adjustment(bayes_factors, ld_matrix)

    # Step 3: Set prior probabilities
    prior_probs = _compute_prior_probabilities(z_scores, prior, n_variants)

    # Step 4: Compute posterior inclusion probabilities
    pips = _compute_pips(bayes_factors, prior_probs)

    # Step 5: Build credible set by ranking PIPs
    ranked_indices = sorted(range(n_variants), key=lambda i: pips[i], reverse=True)
    cumulative = 0.0
    snps_in_set: list[int] = []

    for idx in ranked_indices:
        snps_in_set.append(idx)
        cumulative += pips[idx]
        if cumulative >= coverage:
            break

    return {
        "status": "success",
        "snps_in_set": sorted(snps_in_set),
        "pips": pips,
        "coverage_achieved": float(cumulative),
        "n_snps": len(snps_in_set),
        "bayes_factors": bayes_factors,
    }


def susie_regression(
    X: Any,
    y: Any,
    L: int = 10,
    max_iter: int = 100,
    tol: float = 1e-3,
) -> dict:
    """Simplified SuSiE (Sum of Single Effects) iterative Bayesian regression.

    Pure-Python implementation of the core SuSiE algorithm. Fits a sparse
    regression model as a sum of L single-effect regression components,
    each contributing at most one non-zero coefficient. Iteratively updates
    each component via Bayesian single-effect regression and tracks the
    evidence lower bound (ELBO) for convergence.

    Args:
        X: Genotype matrix (n_samples x n_variants). Accepts numpy array
            or list of lists.
        y: Phenotype vector (n_samples,). Accepts numpy array or list.
        L: Maximum number of causal effects to fit (default 10).
        max_iter: Maximum number of EM iterations (default 100).
        tol: Convergence tolerance on ELBO change (default 1e-3).

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - alpha: L x p matrix of posterior inclusion probabilities per effect
            - mu: L x p matrix of posterior means per effect
            - sigma2: Residual variance estimate
            - pip: Per-variant PIP (1 - product of (1 - alpha) across effects)
            - credible_sets: List of credible set dicts per effect
            - converged: Whether the algorithm converged
            - elbo: Final evidence lower bound
    """
    if not HAS_NUMPY:
        return {"status": "error", "message": "numpy is required for SuSiE regression"}

    try:
        X_arr = np.asarray(X, dtype=float)
        y_arr = np.asarray(y, dtype=float).ravel()
    except (ValueError, TypeError) as exc:
        return {"status": "error", "message": f"Failed to convert inputs: {exc}"}

    n, p = X_arr.shape
    if n != len(y_arr):
        return {
            "status": "error",
            "message": f"X has {n} rows but y has {len(y_arr)} elements",
        }

    if n < 3:
        return {"status": "error", "message": f"Need at least 3 samples, got {n}"}

    if L < 1:
        return {"status": "error", "message": f"L must be >= 1, got {L}"}

    logger.debug(f"Running SuSiE with n={n}, p={p}, L={L}")

    # Center y
    y_mean = float(np.mean(y_arr))
    y_centered = y_arr - y_mean

    # Pre-compute X'X diagonal and X'y
    XtX_diag = np.sum(X_arr**2, axis=0)  # (p,)
    Xty = X_arr.T @ y_centered  # (p,)

    # Initialize
    alpha = np.full((L, p), 1.0 / p)  # Uniform initial inclusion probs
    mu = np.zeros((L, p))  # Posterior means
    mu2 = np.zeros((L, p))  # Posterior second moments
    sigma2 = float(np.var(y_centered))  # Residual variance
    if sigma2 < 1e-12:
        sigma2 = 1.0

    prior_variance = 0.2 * sigma2  # Prior on effect size variance

    converged = False
    prev_elbo = -np.inf

    for iteration in range(max_iter):
        # Compute current fitted values from all effects
        # fitted = sum over l of X @ (alpha[l] * mu[l])
        fitted = np.zeros(n)
        for l_idx in range(L):
            b_l = alpha[l_idx] * mu[l_idx]  # Expected coefficient vector
            fitted += X_arr @ b_l

        for l_idx in range(L):
            # Remove contribution of effect l
            b_l = alpha[l_idx] * mu[l_idx]
            residual = y_centered - fitted + X_arr @ b_l

            # Single effect regression for effect l
            Xtr = X_arr.T @ residual  # (p,)

            # Posterior variance for each variant
            s2 = sigma2 / (XtX_diag + sigma2 / max(prior_variance, 1e-12))
            s2 = np.maximum(s2, 1e-20)

            # Posterior mean for each variant
            mu_l = s2 * Xtr / sigma2

            # Log Bayes factor for each variant
            log_bf = 0.5 * np.log(s2 / max(prior_variance, 1e-12))
            log_bf += 0.5 * mu_l**2 / s2

            # Numerical stability for softmax
            log_bf_max = np.max(log_bf)
            log_bf_shifted = log_bf - log_bf_max

            # Posterior inclusion probabilities (softmax of log BF)
            alpha_l = np.exp(log_bf_shifted)
            alpha_sum = np.sum(alpha_l)
            if alpha_sum > 0:
                alpha_l = alpha_l / alpha_sum
            else:
                alpha_l = np.full(p, 1.0 / p)

            # Update
            alpha[l_idx] = alpha_l
            mu[l_idx] = mu_l
            mu2[l_idx] = mu_l**2 + s2

            # Update fitted values with new effect l
            b_l_new = alpha_l * mu_l
            fitted = fitted - X_arr @ b_l + X_arr @ b_l_new

        # Update residual variance (sigma2)
        residual_final = y_centered - fitted
        sigma2 = float(np.sum(residual_final**2) / n)
        if sigma2 < 1e-12:
            sigma2 = 1e-12

        # Update prior variance (empirical Bayes)
        expected_b2 = 0.0
        for l_idx in range(L):
            expected_b2 += float(np.sum(alpha[l_idx] * mu2[l_idx]))
        if L > 0:
            prior_variance = max(expected_b2 / L, 1e-12)

        # Compute ELBO (simplified)
        elbo = _compute_susie_elbo(y_centered, fitted, sigma2, alpha, mu, mu2, prior_variance)

        # Check convergence
        if abs(elbo - prev_elbo) < tol:
            converged = True
            logger.debug(f"SuSiE converged at iteration {iteration + 1}")
            break

        prev_elbo = elbo

    # Compute per-variant PIP: 1 - prod(1 - alpha[l])
    pip = np.ones(p)
    for l_idx in range(L):
        pip *= 1.0 - alpha[l_idx]
    pip = 1.0 - pip

    # Build credible sets for each effect
    credible_sets = []
    for l_idx in range(L):
        cs = _build_single_effect_cs(alpha[l_idx].tolist(), coverage=0.95)
        credible_sets.append(cs)

    return {
        "status": "success",
        "alpha": alpha.tolist(),
        "mu": mu.tolist(),
        "sigma2": float(sigma2),
        "pip": pip.tolist(),
        "credible_sets": credible_sets,
        "converged": converged,
        "elbo": float(elbo),
    }


def compute_bayes_factors(
    z_scores: list[float],
    prior_variance: float = 0.04,
) -> list[float]:
    """Compute Wakefield approximate Bayes factors from Z-scores.

    The approximate Bayes factor (ABF) for variant j is:

        ABF_j = sqrt(1 + W / V) * exp(-z_j^2 * W / (2 * V * (V + W)))

    where W is the prior variance on the effect size and V = 1 (assuming
    standardized Z-scores with unit variance).

    Args:
        z_scores: Z-scores for each variant.
        prior_variance: Prior variance on the true effect size (W).
            Default 0.04 corresponds to a prior SD of 0.2 on the log-OR
            scale, typical for common-variant GWAS.

    Returns:
        List of approximate Bayes factors, one per variant.

    Raises:
        ValueError: If prior_variance is not positive.
    """
    if prior_variance <= 0:
        raise ValueError(f"prior_variance must be positive, got {prior_variance}")

    if not z_scores:
        return []

    W = prior_variance
    V = 1.0  # Variance of Z-scores under the null

    # Precompute constant terms
    shrinkage = W / (V + W)
    log_scale = 0.5 * math.log(V / (V + W))

    bayes_factors: list[float] = []
    for z in z_scores:
        log_abf = log_scale + 0.5 * z * z * shrinkage
        # Clamp to avoid overflow
        if log_abf > 500:
            abf = math.exp(500)
        elif log_abf < -500:
            abf = 0.0
        else:
            abf = math.exp(log_abf)
        bayes_factors.append(abf)

    return bayes_factors


def colocalization(
    z_scores_1: list[float],
    z_scores_2: list[float],
    ld_matrix: list[list[float]] | None = None,
    prior_p1: float = 1e-4,
    prior_p2: float = 1e-4,
    prior_p12: float = 1e-5,
) -> dict:
    """Coloc-style colocalization analysis between two traits.

    Tests five hypotheses about shared causal variants between two traits
    at a genomic locus:
        H0: No association with either trait
        H1: Association with trait 1 only
        H2: Association with trait 2 only
        H3: Association with both, different causal variants
        H4: Association with both, shared causal variant

    Uses approximate Bayes factors from Z-scores and combines evidence
    across all variants in the region.

    Args:
        z_scores_1: Z-scores for trait 1 at each variant.
        z_scores_2: Z-scores for trait 2 at each variant.
        ld_matrix: Optional LD matrix. Currently used for future LD-aware
            extensions but does not modify the core coloc computation.
        prior_p1: Prior probability a variant is causal for trait 1 only.
        prior_p2: Prior probability a variant is causal for trait 2 only.
        prior_p12: Prior probability a variant is causal for both traits.

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - PP_H0: Posterior probability of hypothesis 0
            - PP_H1: Posterior probability of hypothesis 1
            - PP_H2: Posterior probability of hypothesis 2
            - PP_H3: Posterior probability of hypothesis 3
            - PP_H4: Posterior probability of hypothesis 4
            - summary: String describing the most supported hypothesis
    """
    if not z_scores_1 or not z_scores_2:
        return {"status": "error", "message": "Both trait Z-scores are required"}

    if len(z_scores_1) != len(z_scores_2):
        return {
            "status": "error",
            "message": (f"Z-score lengths differ: trait 1 has {len(z_scores_1)}, " f"trait 2 has {len(z_scores_2)}"),
        }

    n_variants = len(z_scores_1)

    # Validate priors
    if any(p <= 0 for p in (prior_p1, prior_p2, prior_p12)):
        return {"status": "error", "message": "Prior probabilities must be positive"}

    prior_p0 = 1.0 - prior_p1 - prior_p2 - prior_p12
    if prior_p0 < 0:
        return {"status": "error", "message": "Prior probabilities sum to more than 1"}

    logger.debug(f"Running colocalization for {n_variants} variants")

    # Compute per-variant Bayes factors for each trait
    bf1 = compute_bayes_factors(z_scores_1)
    bf2 = compute_bayes_factors(z_scores_2)

    # Sum of Bayes factors across variants
    sum_bf1 = sum(bf1)
    sum_bf2 = sum(bf2)

    # Sum of product of Bayes factors (shared causal variant)
    sum_bf12 = sum(b1 * b2 for b1, b2 in zip(bf1, bf2))

    # Sum of cross-products (different causal variants)
    # H3: sum_i sum_{j != i} BF1_i * BF2_j = sum_BF1 * sum_BF2 - sum_BF12
    cross_bf = sum_bf1 * sum_bf2 - sum_bf12

    # Log-space computation for numerical stability
    # Unnormalized log-posteriors for each hypothesis
    # H0: no association
    log_h0 = math.log(max(prior_p0, 1e-300)) + math.log(max(n_variants, 1))

    # H1: causal for trait 1 only
    log_h1 = math.log(max(prior_p1, 1e-300)) + math.log(max(sum_bf1, 1e-300))

    # H2: causal for trait 2 only
    log_h2 = math.log(max(prior_p2, 1e-300)) + math.log(max(sum_bf2, 1e-300))

    # H3: causal for both, different variants
    log_h3_scale = prior_p1 * prior_p2
    log_h3 = math.log(max(log_h3_scale, 1e-300)) + math.log(max(cross_bf, 1e-300))

    # H4: shared causal variant
    log_h4 = math.log(max(prior_p12, 1e-300)) + math.log(max(sum_bf12, 1e-300))

    # Normalize to posterior probabilities
    log_values = [log_h0, log_h1, log_h2, log_h3, log_h4]
    log_max = max(log_values)

    posteriors = [math.exp(lv - log_max) for lv in log_values]
    total = sum(posteriors)
    if total > 0:
        posteriors = [p / total for p in posteriors]
    else:
        posteriors = [0.2] * 5

    pp_h0, pp_h1, pp_h2, pp_h3, pp_h4 = posteriors

    # Determine most supported hypothesis
    hypothesis_names = {
        0: "H0: No association with either trait",
        1: "H1: Association with trait 1 only",
        2: "H2: Association with trait 2 only",
        3: "H3: Both traits associated, different causal variants",
        4: "H4: Both traits associated, shared causal variant",
    }
    best_idx = posteriors.index(max(posteriors))
    summary = f"Most supported: {hypothesis_names[best_idx]} (PP={posteriors[best_idx]:.4f})"

    return {
        "status": "success",
        "PP_H0": float(pp_h0),
        "PP_H1": float(pp_h1),
        "PP_H2": float(pp_h2),
        "PP_H3": float(pp_h3),
        "PP_H4": float(pp_h4),
        "summary": summary,
    }


def conditional_analysis(
    z_scores: list[float],
    ld_matrix: list[list[float]],
    max_signals: int = 10,
    p_threshold: float = 5e-8,
) -> list[dict]:
    """Stepwise conditional analysis to identify independent association signals.

    Iteratively identifies the most significant variant, conditions on it
    by adjusting Z-scores using the LD matrix, and repeats until no variant
    exceeds the significance threshold or the maximum number of signals
    is reached.

    Args:
        z_scores: Z-scores for each variant.
        ld_matrix: LD (correlation) matrix, size n_variants x n_variants.
        max_signals: Maximum number of independent signals to identify
            (default 10).
        p_threshold: Genome-wide significance threshold (default 5e-8).

    Returns:
        List of dicts, one per identified signal, each containing:
            - index: Variant index of the lead SNP
            - z_score: Original Z-score
            - z_conditional: Conditional Z-score
            - p_value: P-value from the conditional Z-score
            - signal_number: 1-indexed signal number
    """
    if not z_scores:
        return []

    n_variants = len(z_scores)

    # Validate LD matrix
    if len(ld_matrix) != n_variants:
        logger.error(f"LD matrix rows ({len(ld_matrix)}) does not match " f"Z-score count ({n_variants})")
        return []

    if not HAS_NUMPY:
        logger.warning("numpy not available; using pure Python conditional analysis")
        return _conditional_analysis_pure(z_scores, ld_matrix, max_signals, p_threshold)

    logger.debug(f"Running conditional analysis: {n_variants} variants, " f"max_signals={max_signals}")

    z_arr = np.array(z_scores, dtype=float)
    ld_arr = np.array(ld_matrix, dtype=float)

    # Current (conditional) Z-scores
    z_cond = z_arr.copy()
    conditioned_indices: list[int] = []
    signals: list[dict] = []

    # Z-score threshold from p-value (two-tailed)
    z_threshold = _p_to_z(p_threshold)

    for signal_num in range(1, max_signals + 1):
        # Find the variant with the largest absolute conditional Z-score
        abs_z = np.abs(z_cond)
        lead_idx = int(np.argmax(abs_z))
        lead_z = float(z_cond[lead_idx])

        # Check significance
        if abs(lead_z) < z_threshold:
            logger.debug(
                f"No more significant signals at step {signal_num} "
                f"(max |Z|={abs(lead_z):.2f}, threshold={z_threshold:.2f})"
            )
            break

        p_val = _z_to_p(lead_z)

        signals.append(
            {
                "index": lead_idx,
                "z_score": float(z_arr[lead_idx]),
                "z_conditional": float(lead_z),
                "p_value": float(p_val),
                "signal_number": signal_num,
            }
        )

        conditioned_indices.append(lead_idx)

        # Condition out this signal: adjust Z-scores
        # z_cond = z_cond - ld[., lead] * z_cond[lead] / ld[lead, lead]
        ld_col = ld_arr[:, lead_idx]
        ld_diag = ld_arr[lead_idx, lead_idx]
        if abs(ld_diag) > 1e-12:
            adjustment = ld_col * (lead_z / ld_diag)
            z_cond = z_cond - adjustment
            # Zero out the lead variant itself
            z_cond[lead_idx] = 0.0
        else:
            z_cond[lead_idx] = 0.0

    logger.info(f"Conditional analysis identified {len(signals)} independent signals")
    return signals


def annotate_credible_set(
    credible_set: dict,
    annotations: dict | None = None,
) -> dict:
    """Enrich credible set SNPs with functional annotations.

    Takes a credible set result (from compute_credible_set) and overlays
    functional annotation categories. Computes enrichment of causal
    posterior probability in each annotation category relative to baseline
    expectation.

    Args:
        credible_set: Output from compute_credible_set, must contain
            "snps_in_set" and "pips" keys.
        annotations: Optional dict mapping annotation category names to
            lists of variant indices belonging to that category. Example:
            {"coding": [0, 5, 12], "enhancer": [1, 3, 7, 15]}.
            If None, returns the credible set with empty annotation fields.

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - credible_set: The original credible set data
            - annotated_snps: List of dicts with index, pip, and categories
            - enrichment: Dict mapping category to enrichment fold and p_pip
            - n_annotated: Number of annotated variants in the set
    """
    if credible_set.get("status") != "success":
        return {
            "status": "error",
            "message": "Input credible set does not have success status",
        }

    snps_in_set = credible_set.get("snps_in_set", [])
    pips = credible_set.get("pips", [])

    if not snps_in_set:
        return {
            "status": "error",
            "message": "Credible set contains no variants",
        }

    n_total = len(pips)
    set_indices = set(snps_in_set)

    # Build per-SNP annotation lookup
    snp_categories: dict[int, list[str]] = {idx: [] for idx in snps_in_set}

    if annotations:
        for category, indices in annotations.items():
            for idx in indices:
                if idx in snp_categories:
                    snp_categories[idx].append(category)

    # Build annotated SNP list
    annotated_snps: list[dict] = []
    for idx in sorted(snps_in_set):
        pip_val = pips[idx] if idx < len(pips) else 0.0
        annotated_snps.append(
            {
                "index": idx,
                "pip": float(pip_val),
                "categories": snp_categories.get(idx, []),
            }
        )

    # Compute enrichment per annotation category
    enrichment: dict[str, dict] = {}
    if annotations:
        total_pip_in_set = sum(pips[i] for i in snps_in_set if i < len(pips))

        for category, cat_indices in annotations.items():
            cat_set = set(cat_indices)

            # PIP in this category within the credible set
            pip_in_cat = sum(pips[i] for i in snps_in_set if i < len(pips) and i in cat_set)

            # Expected fraction: proportion of all variants in this category
            n_in_cat_total = len(cat_set)
            expected_fraction = n_in_cat_total / max(n_total, 1)

            # Observed fraction of PIP
            observed_fraction = pip_in_cat / max(total_pip_in_set, 1e-12)

            # Enrichment fold
            fold = observed_fraction / max(expected_fraction, 1e-12)

            enrichment[category] = {
                "pip_in_category": float(pip_in_cat),
                "n_variants_in_category": len(cat_set & set_indices),
                "expected_fraction": float(expected_fraction),
                "observed_fraction": float(observed_fraction),
                "enrichment_fold": float(fold),
            }

    n_annotated = sum(1 for snp in annotated_snps if snp["categories"])

    return {
        "status": "success",
        "credible_set": credible_set,
        "annotated_snps": annotated_snps,
        "enrichment": enrichment,
        "n_annotated": n_annotated,
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _apply_ld_adjustment(
    bayes_factors: list[float],
    ld_matrix: list[list[float]],
) -> list[float]:
    """Adjust Bayes factors using LD structure.

    Applies a penalty to variants in high LD with stronger signals, reducing
    the effective Bayes factor of tagging variants.

    Args:
        bayes_factors: Raw Bayes factors.
        ld_matrix: LD correlation matrix.

    Returns:
        Adjusted Bayes factors.
    """
    n = len(bayes_factors)
    if not HAS_NUMPY:
        return bayes_factors

    bf_arr = np.array(bayes_factors, dtype=float)
    ld_arr = np.array(ld_matrix, dtype=float)

    # For each variant, reduce BF proportionally to LD with higher-BF variants
    adjusted = bf_arr.copy()
    ranked = np.argsort(bf_arr)[::-1]  # Descending BF order

    for rank, idx in enumerate(ranked):
        if rank == 0:
            continue
        # Find maximum LD with any higher-ranked variant
        max_ld_with_stronger = 0.0
        for prev_rank in range(rank):
            prev_idx = ranked[prev_rank]
            ld_val = abs(ld_arr[idx, prev_idx])
            if ld_val > max_ld_with_stronger:
                max_ld_with_stronger = ld_val

        # Penalize: BF * (1 - r^2) to down-weight tagging variants
        penalty = 1.0 - max_ld_with_stronger**2
        adjusted[idx] = bf_arr[idx] * max(penalty, 0.01)

    return adjusted.tolist()


def _compute_prior_probabilities(
    z_scores: list[float],
    prior: str,
    n_variants: int,
) -> list[float]:
    """Compute prior probabilities for each variant.

    Args:
        z_scores: Z-scores (used for "distance" prior to find lead).
        prior: Prior type ("uniform" or "distance").
        n_variants: Number of variants.

    Returns:
        List of prior probabilities summing to 1.
    """
    if prior == "uniform" or n_variants == 0:
        return [1.0 / n_variants] * n_variants

    # Distance-from-lead prior: higher prior near the lead variant
    lead_idx = max(range(n_variants), key=lambda i: abs(z_scores[i]))

    priors: list[float] = []
    for i in range(n_variants):
        distance = abs(i - lead_idx)
        # Exponential decay with distance
        weight = math.exp(-0.1 * distance)
        priors.append(weight)

    total = sum(priors)
    if total > 0:
        priors = [p / total for p in priors]
    else:
        priors = [1.0 / n_variants] * n_variants

    return priors


def _compute_pips(
    bayes_factors: list[float],
    prior_probs: list[float],
) -> list[float]:
    """Compute posterior inclusion probabilities from BFs and priors.

    PIP_i = BF_i * prior_i / sum_j(BF_j * prior_j)

    Args:
        bayes_factors: Per-variant Bayes factors.
        prior_probs: Per-variant prior probabilities.

    Returns:
        List of PIPs summing to 1.
    """
    n = len(bayes_factors)
    if n == 0:
        return []

    products = [bf * pr for bf, pr in zip(bayes_factors, prior_probs)]
    total = sum(products)

    if total > 0:
        return [p / total for p in products]
    else:
        return [1.0 / n] * n


def _build_single_effect_cs(
    alpha_l: list[float],
    coverage: float = 0.95,
) -> dict:
    """Build a credible set for a single SuSiE effect.

    Args:
        alpha_l: Posterior inclusion probabilities for one effect.
        coverage: Target coverage (default 0.95).

    Returns:
        Dict with indices, coverage_achieved, and size.
    """
    n = len(alpha_l)
    if n == 0:
        return {"indices": [], "coverage_achieved": 0.0, "size": 0}

    ranked = sorted(range(n), key=lambda i: alpha_l[i], reverse=True)
    cumulative = 0.0
    indices: list[int] = []

    for idx in ranked:
        indices.append(idx)
        cumulative += alpha_l[idx]
        if cumulative >= coverage:
            break

    return {
        "indices": sorted(indices),
        "coverage_achieved": float(cumulative),
        "size": len(indices),
    }


def _compute_susie_elbo(
    y: Any,
    fitted: Any,
    sigma2: float,
    alpha: Any,
    mu: Any,
    mu2: Any,
    prior_variance: float,
) -> float:
    """Compute the evidence lower bound (ELBO) for SuSiE.

    Simplified ELBO: data log-likelihood minus KL divergences.

    Args:
        y: Centered phenotype vector.
        fitted: Current fitted values.
        sigma2: Residual variance.
        alpha: L x p inclusion probability matrix.
        mu: L x p posterior mean matrix.
        mu2: L x p posterior second moment matrix.
        prior_variance: Prior variance on effects.

    Returns:
        ELBO value.
    """
    n = len(y)
    residual = y - fitted
    rss = float(np.sum(residual**2))

    # Data log-likelihood
    ll = -0.5 * n * math.log(2 * math.pi * max(sigma2, 1e-12))
    ll -= 0.5 * rss / max(sigma2, 1e-12)

    # KL divergence for each effect (simplified)
    L = alpha.shape[0]
    kl_total = 0.0
    for l_idx in range(L):
        # Entropy of the categorical distribution
        for j in range(alpha.shape[1]):
            a = alpha[l_idx, j]
            if a > 1e-30:
                kl_total -= a * math.log(a)
                # Prior contribution
                kl_total += a * 0.5 * math.log(max(prior_variance, 1e-12) / max(sigma2, 1e-12))

    return ll + kl_total


def _conditional_analysis_pure(
    z_scores: list[float],
    ld_matrix: list[list[float]],
    max_signals: int,
    p_threshold: float,
) -> list[dict]:
    """Pure Python fallback for conditional analysis without numpy.

    Args:
        z_scores: Z-scores for each variant.
        ld_matrix: LD matrix.
        max_signals: Maximum signals to find.
        p_threshold: Significance threshold.

    Returns:
        List of signal dicts.
    """
    n = len(z_scores)
    z_cond = list(z_scores)
    signals: list[dict] = []
    z_threshold = _p_to_z(p_threshold)

    for signal_num in range(1, max_signals + 1):
        # Find lead variant
        lead_idx = 0
        max_abs_z = abs(z_cond[0])
        for i in range(1, n):
            if abs(z_cond[i]) > max_abs_z:
                max_abs_z = abs(z_cond[i])
                lead_idx = i

        lead_z = z_cond[lead_idx]

        if abs(lead_z) < z_threshold:
            break

        p_val = _z_to_p(lead_z)

        signals.append(
            {
                "index": lead_idx,
                "z_score": z_scores[lead_idx],
                "z_conditional": float(lead_z),
                "p_value": float(p_val),
                "signal_number": signal_num,
            }
        )

        # Condition out this signal
        ld_diag = ld_matrix[lead_idx][lead_idx]
        if abs(ld_diag) > 1e-12:
            for i in range(n):
                adjustment = ld_matrix[i][lead_idx] * (lead_z / ld_diag)
                z_cond[i] -= adjustment
            z_cond[lead_idx] = 0.0
        else:
            z_cond[lead_idx] = 0.0

    return signals


def _p_to_z(p: float) -> float:
    """Convert a two-tailed p-value to a Z-score threshold.

    Args:
        p: Two-tailed p-value.

    Returns:
        Positive Z-score threshold.
    """
    if p <= 0 or p >= 1:
        return 0.0
    # Use inverse of normal CDF approximation (Beasley-Springer-Moro)
    return _normal_quantile(1.0 - p / 2.0)


def _z_to_p(z: float) -> float:
    """Convert a Z-score to a two-tailed p-value.

    Args:
        z: Z-score.

    Returns:
        Two-tailed p-value.
    """
    return 2.0 * (1.0 - _normal_cdf(abs(z)))


def _normal_cdf(x: float) -> float:
    """Approximate standard normal CDF using Abramowitz & Stegun.

    Args:
        x: Input value.

    Returns:
        CDF value.
    """
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    sign = 1 if x >= 0 else -1
    x = abs(x) / math.sqrt(2.0)

    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-x * x)

    return 0.5 * (1.0 + sign * y)


def _normal_quantile(p: float) -> float:
    """Approximate quantile function of the standard normal distribution.

    Uses the rational approximation from Abramowitz & Stegun 26.2.23.

    Args:
        p: Probability in (0, 1).

    Returns:
        Z-score such that P(Z <= z) = p.
    """
    if p <= 0.0:
        return -10.0
    if p >= 1.0:
        return 10.0
    if p == 0.5:
        return 0.0

    if p < 0.5:
        return -_normal_quantile(1.0 - p)

    # Rational approximation for p > 0.5
    t = math.sqrt(-2.0 * math.log(1.0 - p))
    c0 = 2.515517
    c1 = 0.802853
    c2 = 0.010328
    d1 = 1.432788
    d2 = 0.189269
    d3 = 0.001308

    z = t - (c0 + c1 * t + c2 * t * t) / (1.0 + d1 * t + d2 * t * t + d3 * t * t * t)
    return z
