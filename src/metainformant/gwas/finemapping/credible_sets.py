"""Statistical fine-mapping: credible sets, SuSiE, Bayes factors, and colocalization.

Implements core fine-mapping algorithms for identifying causal variants from
GWAS summary statistics. Includes approximate Bayes factor credible sets,
a pure-Python SuSiE (Sum of Single Effects) regression, Wakefield approximate
Bayes factors, coloc-style colocalization, stepwise conditional analysis,
and functional annotation enrichment of credible sets.

Approximate Bayes factors (ABFs) follow Wakefield (2009): the estimate
``beta_hat ~ N(beta, V)`` with ``V = SE**2`` is compared under the prior
``beta ~ N(0, W)``, giving ``ABF = sqrt(V / (V + W)) * exp(z**2 * W / (2 * (V + W)))``.
``V`` and ``W`` must be on the same scale, so z-scores alone do not fix an ABF:
callers supply standard errors, or a sample size (with an optional MAF and case
fraction) from which ``V`` is modelled. The prior defaults are coloc's
(``W = 0.15**2`` quantitative, ``W = 0.2**2`` case-control). All Bayes factors,
posterior inclusion probabilities and colocalization posteriors are accumulated
in log space, because correctly scaled ABFs reach ``exp(z**2 / 2)`` and overflow
double precision for ``|z|`` above about 38.

References:
    - Wakefield (2009) Genetic Epidemiology 33:79-86 (ABF)
    - Wang et al. (2020) JRSS-B 82:1273-1300 (SuSiE)
    - Giambartolomei et al. (2014) PLoS Genetics 10:e1004383 (coloc)
    - Yang et al. (2012) Nature Genetics 44:369-375 (conditional analysis)
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from numbers import Real
from typing import TYPE_CHECKING, Any, NamedTuple, cast

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    if not TYPE_CHECKING:  # runtime-only fallback; type checkers keep the real module type
        np = None


# ---------------------------------------------------------------------------
# Prior defaults for approximate Bayes factors
# ---------------------------------------------------------------------------

#: Prior SD of the true effect for a quantitative trait, in phenotype-SD units
#: (coloc's ``sd.prior`` default).
DEFAULT_PRIOR_SD_QUANTITATIVE = 0.15
#: Prior SD of the true log odds ratio for a case-control trait (coloc's default).
DEFAULT_PRIOR_SD_CASE_CONTROL = 0.2
#: Default prior variance ``W`` on the beta scale, quantitative trait.
DEFAULT_PRIOR_VARIANCE_QUANTITATIVE = DEFAULT_PRIOR_SD_QUANTITATIVE**2
#: Default prior variance ``W`` on the beta scale, case-control trait.
DEFAULT_PRIOR_VARIANCE_CASE_CONTROL = DEFAULT_PRIOR_SD_CASE_CONTROL**2
#: Prior SD of the true z-score noncentrality. Used only on the z-score scale
#: (no standard errors and no sample size), where it replaces the beta-scale
#: prior. 5.2 is close to the z-score at genome-wide significance (5.45 for 5e-8).
DEFAULT_PRIOR_SD_Z_SCALE = 5.2
#: Default prior variance of the noncentrality on the z-score scale.
DEFAULT_PRIOR_VARIANCE_Z_SCALE = DEFAULT_PRIOR_SD_Z_SCALE**2

#: Largest log Bayes factor reported on the linear scale (``exp(500)`` is finite).
_MAX_LOG_BAYES_FACTOR = 500.0

#: How the effect variance ``V`` of an ABF was obtained.
VARIANCE_MODEL_STANDARD_ERRORS = "standard_errors"
VARIANCE_MODEL_SAMPLE_SIZE = "sample_size"
VARIANCE_MODEL_Z_SCALE = "z_scale"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_credible_set(
    z_scores: list[float],
    ld_matrix: list[list[float]] | None = None,
    coverage: float = 0.95,
    prior: str = "uniform",
    *,
    prior_variance: float | None = None,
    standard_errors: float | Sequence[float] | None = None,
    sample_size: float | Sequence[float] | None = None,
    maf: float | Sequence[float] | None = None,
    case_fraction: float | None = None,
) -> dict:
    """Compute an approximate Bayes factor credible set from Z-scores.

    Calculates posterior inclusion probabilities (PIPs) from GWAS Z-scores
    using Wakefield (2009) approximate Bayes factors, optionally adjusting for
    linkage disequilibrium. Returns the smallest set of variants whose
    cumulative PIP reaches the requested coverage, assuming one causal variant
    in the region.

    The Bayes factors depend on the scale of the effect estimates, which
    z-scores alone do not carry. Pass ``standard_errors`` (V = SE**2), or
    ``sample_size`` (with ``maf`` and ``case_fraction`` as available) so that V
    is modelled as in coloc. See :func:`compute_log_bayes_factors` for the exact
    model and the default priors (W = 0.15**2 quantitative, 0.2**2 case-control).
    With neither, the Bayes factors are computed on the z-score scale with a
    documented noncentrality prior (SD 5.2), and the result reports
    ``variance_model == "z_scale"``.

    Args:
        z_scores: Z-scores for each variant in the region.
        ld_matrix: Optional LD correlation matrix r (signed entries in
            [-1, 1], ones on the diagonal), square with len(z_scores) rows and
            columns. It is not an r-squared matrix: the adjustment squares it.
            Each variant's Bayes factor is multiplied by ``1 - r**2``, where r
            is its largest absolute correlation with a variant that has a
            larger Bayes factor (the factor never drops below 0.01), so a
            variant that tags a stronger signal loses posterior mass to it.
        coverage: Target cumulative posterior probability for the credible
            set (default 0.95).
        prior: Prior type on causal status. One of "uniform" (equal prior
            per variant) or "distance" (decaying prior from lead variant).
        prior_variance: Prior variance W of the true effect (beta scale when
            standard errors or a sample size are given, z-score scale otherwise).
            Defaults to the coloc priors described above.
        standard_errors: Standard error of each effect estimate (scalar or one
            per variant). Mutually exclusive with ``sample_size``.
        sample_size: Total sample size N (scalar or one per variant) used to
            model V when standard errors are unavailable.
        maf: Allele frequency (scalar or one per variant) for the sample-size
            variance model. Without it, genotypes are standardised to unit
            variance.
        case_fraction: Fraction of cases in (0, 1) for a case-control trait.
            Leave as None for a quantitative trait.

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - snps_in_set: List of indices in the credible set
            - pips: List of posterior inclusion probabilities for all variants
            - coverage_achieved: Actual cumulative PIP of the credible set
            - n_snps: Number of variants in the credible set
            - bayes_factors: List of approximate Bayes factors (capped at
              exp(500) so they stay finite; PIPs use the uncapped logs)
            - log_bayes_factors: Natural-log approximate Bayes factors
            - variance_model: "standard_errors", "sample_size" or "z_scale"
            - prior_variance: The prior variance W that was used
    """
    if len(z_scores) == 0:
        return {"status": "error", "message": "No Z-scores provided"}

    n_variants = len(z_scores)

    if coverage <= 0.0 or coverage > 1.0:
        return {
            "status": "error",
            "message": f"Coverage must be in (0, 1], got {coverage}",
        }

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

    # Step 1: Compute Wakefield approximate Bayes factors (log scale)
    try:
        scaled = _wakefield_log_abf(
            z_scores,
            prior_variance,
            standard_errors=standard_errors,
            sample_size=sample_size,
            maf=maf,
            case_fraction=case_fraction,
        )
    except ValueError as exc:
        return {"status": "error", "message": str(exc)}
    log_bayes_factors = scaled.log_bayes_factors

    # Step 2: Apply LD adjustment if matrix provided
    if ld_matrix is not None and HAS_NUMPY:
        log_bayes_factors = _apply_ld_adjustment(log_bayes_factors, ld_matrix)

    # Step 3: Set prior probabilities
    prior_probs = _compute_prior_probabilities(z_scores, prior, n_variants)

    # Step 4: Compute posterior inclusion probabilities (log-space, overflow-safe)
    pips = _compute_pips(log_bayes_factors, prior_probs)

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
        "bayes_factors": [_exp_capped(value) for value in log_bayes_factors],
        "log_bayes_factors": log_bayes_factors,
        "variance_model": scaled.variance_model,
        "prior_variance": scaled.prior_variance,
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
    X_arr.T @ y_centered  # (p,)

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


def compute_log_bayes_factors(
    z_scores: Sequence[float],
    prior_variance: float | None = None,
    *,
    standard_errors: float | Sequence[float] | None = None,
    sample_size: float | Sequence[float] | None = None,
    maf: float | Sequence[float] | None = None,
    case_fraction: float | None = None,
) -> list[float]:
    """Compute natural-log Wakefield (2009) approximate Bayes factors, one per variant.

    Wakefield models the effect estimate as ``beta_hat ~ N(beta, V)`` with
    ``V = SE**2`` and the true effect as ``beta ~ N(0, W)``. The Bayes factor
    for association (``beta != 0``) against the null is

        ABF = sqrt(V / (V + W)) * exp(z**2 * W / (2 * (V + W))),   z = beta_hat / sqrt(V)

    or, with the shrinkage ``r = W / (V + W)``, ``log ABF = (log(1 - r) + r * z**2) / 2``.
    The ABF depends on ``z`` and on ``V / W`` only, so a z-score is enough once the
    scale of ``V`` is known. ``V`` and ``W`` must be on the same scale: pairing a
    beta-scale prior with ``V = 1`` makes the Bayes factors nearly flat, so that a
    z = 7.5 hit among 300 variants cannot be told from the nulls.

    Exactly one source of ``V`` is used:

    * ``standard_errors``: ``V = SE**2``, exact.
    * ``sample_size``: ``V`` is modelled as in coloc (``Var.data`` and
      ``Var.data.cc``) for a phenotype with unit variance,
      ``V = 1 / (2 N f (1 - f))`` for a quantitative trait and
      ``V = 1 / (2 N f (1 - f) s (1 - s))`` for a case-control trait, with allele
      frequency ``f = maf`` and case fraction ``s = case_fraction``. Without
      ``maf`` the genotype is standardised to unit variance (the ``2 f (1 - f)``
      factor becomes 1), and ``W`` then applies per SD of genotype.
    * neither: the z-score scale. A z-score has unit sampling variance, so
      ``V = 1`` is exact on that scale, and ``W`` is the prior variance of the
      noncentrality ``beta / SE``. This is a coarse fallback because it fixes the
      effect size relative to the noise instead of measuring it; prefer one of
      the two sources above.

    The prior variance ``W`` defaults to coloc's: ``0.15**2`` for a quantitative
    trait in phenotype-SD units (pass ``prior_variance=(0.15 * sd_y) ** 2`` for
    raw units), ``0.2**2`` on the log-odds scale when ``case_fraction`` marks a
    case-control trait, and ``5.2**2`` on the z-score scale.

    Args:
        z_scores: Z-scores for each variant (``beta_hat / SE``).
        prior_variance: Prior variance ``W`` of the true effect. It lives on the
            beta scale when ``standard_errors`` or ``sample_size`` is given and on
            the z-score scale otherwise. ``None`` selects the default above.
        standard_errors: Standard error of each effect estimate, a scalar or one
            value per variant. Mutually exclusive with ``sample_size``.
        sample_size: Total sample size ``N``, a scalar or one value per variant.
        maf: Allele frequency in (0, 1), a scalar or one value per variant. Used
            only by the sample-size variance model.
        case_fraction: Fraction of cases in (0, 1) for a case-control trait;
            ``None`` means a quantitative trait. It selects the default prior and
            enters the sample-size variance model.

    Returns:
        Natural-log approximate Bayes factors (association against null), one per
        variant. Empty for empty input.

    Raises:
        ValueError: If ``prior_variance`` is not positive, if both
            ``standard_errors`` and ``sample_size`` are given, or if an input is
            non-finite, out of range or has the wrong length.
    """
    return _wakefield_log_abf(
        z_scores,
        prior_variance,
        standard_errors=standard_errors,
        sample_size=sample_size,
        maf=maf,
        case_fraction=case_fraction,
    ).log_bayes_factors


def compute_bayes_factors(
    z_scores: Sequence[float],
    prior_variance: float | None = None,
    *,
    standard_errors: float | Sequence[float] | None = None,
    sample_size: float | Sequence[float] | None = None,
    maf: float | Sequence[float] | None = None,
    case_fraction: float | None = None,
) -> list[float]:
    """Compute Wakefield (2009) approximate Bayes factors on the linear scale.

    ``ABF = sqrt(V / (V + W)) * exp(z**2 * W / (2 * (V + W)))``, where ``V`` is
    the squared standard error and ``W`` the prior variance of the true effect.
    Takes the same arguments as :func:`compute_log_bayes_factors`, which documents
    how ``V`` and ``W`` are obtained and their defaults.

    Correctly scaled ABFs grow like ``exp(z**2 / 2)``, so values are capped at
    ``exp(500)`` to stay finite. Use :func:`compute_log_bayes_factors` when
    ranking or combining strong signals.

    Returns:
        List of approximate Bayes factors, one per variant.

    Raises:
        ValueError: As for :func:`compute_log_bayes_factors`.
    """
    log_bayes_factors = compute_log_bayes_factors(
        z_scores,
        prior_variance,
        standard_errors=standard_errors,
        sample_size=sample_size,
        maf=maf,
        case_fraction=case_fraction,
    )
    return [_exp_capped(value) for value in log_bayes_factors]


def log_sum_exp(values: Iterable[float]) -> float:
    """Return ``log(sum(exp(v)))`` without overflow or underflow.

    Args:
        values: Natural-log quantities; ``-inf`` entries contribute zero.

    Returns:
        The log of the sum, or ``-inf`` for empty input or all ``-inf`` values.
    """
    items = list(values)
    if not items:
        return -math.inf
    peak = max(items)
    if math.isinf(peak):
        return peak
    return peak + math.log(math.fsum(math.exp(value - peak) for value in items))


def colocalization(
    z_scores_1: list[float],
    z_scores_2: list[float],
    ld_matrix: list[list[float]] | None = None,
    prior_p1: float = 1e-4,
    prior_p2: float = 1e-4,
    prior_p12: float = 1e-5,
    *,
    prior_variance_1: float | None = None,
    prior_variance_2: float | None = None,
    standard_errors_1: float | Sequence[float] | None = None,
    standard_errors_2: float | Sequence[float] | None = None,
    sample_size_1: float | Sequence[float] | None = None,
    sample_size_2: float | Sequence[float] | None = None,
    maf: float | Sequence[float] | None = None,
    case_fraction_1: float | None = None,
    case_fraction_2: float | None = None,
) -> dict:
    """Coloc-style colocalization analysis between two traits.

    Tests five hypotheses about shared causal variants between two traits
    at a genomic locus:
        H0: No association with either trait
        H1: Association with trait 1 only
        H2: Association with trait 2 only
        H3: Association with both, different causal variants
        H4: Association with both, shared causal variant

    Uses Wakefield approximate Bayes factors from Z-scores and combines
    evidence across all variants in the region (Giambartolomei et al. 2014),
    entirely in log space so strong signals cannot overflow.

    The Bayes factors need the scale of each trait's effect estimates. Pass
    ``standard_errors_*`` or ``sample_size_*`` (with ``maf`` and ``case_fraction_*``
    as available) for each trait; see :func:`compute_log_bayes_factors` for the
    model and the default priors. A trait without either is scored on the
    z-score scale (noncentrality prior SD 5.2) and reported as such.

    Args:
        z_scores_1: Z-scores for trait 1 at each variant.
        z_scores_2: Z-scores for trait 2 at each variant.
        ld_matrix: Optional LD matrix. Currently used for future LD-aware
            extensions but does not modify the core coloc computation.
        prior_p1: Prior probability a variant is causal for trait 1 only.
        prior_p2: Prior probability a variant is causal for trait 2 only.
        prior_p12: Prior probability a variant is causal for both traits.
        prior_variance_1: Prior variance W of the true effect for trait 1.
        prior_variance_2: Prior variance W of the true effect for trait 2.
        standard_errors_1: Standard errors of trait 1 (scalar or per variant).
        standard_errors_2: Standard errors of trait 2 (scalar or per variant).
        sample_size_1: Sample size of trait 1 (scalar or per variant).
        sample_size_2: Sample size of trait 2 (scalar or per variant).
        maf: Allele frequency shared by both traits, for the sample-size model.
        case_fraction_1: Case fraction of trait 1; None for a quantitative trait.
        case_fraction_2: Case fraction of trait 2; None for a quantitative trait.

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - PP_H0: Posterior probability of hypothesis 0
            - PP_H1: Posterior probability of hypothesis 1
            - PP_H2: Posterior probability of hypothesis 2
            - PP_H3: Posterior probability of hypothesis 3
            - PP_H4: Posterior probability of hypothesis 4
            - summary: String describing the most supported hypothesis
            - snp_pp_h4: Per-variant probability of being the shared causal
              variant given H4 (coloc's SNP.PP.H4); sums to 1
            - bayes_factor_model: Per trait ("trait_1", "trait_2"), the
              ``variance_model`` and ``prior_variance`` behind its Bayes factors
    """
    if len(z_scores_1) == 0 or len(z_scores_2) == 0:
        return {"status": "error", "message": "Both trait Z-scores are required"}

    if len(z_scores_1) != len(z_scores_2):
        return {
            "status": "error",
            "message": (f"Z-score lengths differ: trait 1 has {len(z_scores_1)}, " f"trait 2 has {len(z_scores_2)}"),
        }

    n_variants = len(z_scores_1)
    logger.debug(f"Running colocalization for {n_variants} variants")

    try:
        trait_1 = _wakefield_log_abf(
            z_scores_1,
            prior_variance_1,
            standard_errors=standard_errors_1,
            sample_size=sample_size_1,
            maf=maf,
            case_fraction=case_fraction_1,
        )
        trait_2 = _wakefield_log_abf(
            z_scores_2,
            prior_variance_2,
            standard_errors=standard_errors_2,
            sample_size=sample_size_2,
            maf=maf,
            case_fraction=case_fraction_2,
        )
    except ValueError as exc:
        return {"status": "error", "message": str(exc)}

    result = colocalization_from_log_bayes_factors(
        trait_1.log_bayes_factors,
        trait_2.log_bayes_factors,
        prior_p1=prior_p1,
        prior_p2=prior_p2,
        prior_p12=prior_p12,
    )
    if result.get("status") == "success":
        result["bayes_factor_model"] = {
            "trait_1": {"variance_model": trait_1.variance_model, "prior_variance": trait_1.prior_variance},
            "trait_2": {"variance_model": trait_2.variance_model, "prior_variance": trait_2.prior_variance},
        }
    return result


def colocalization_from_log_bayes_factors(
    log_bayes_factors_1: Sequence[float],
    log_bayes_factors_2: Sequence[float],
    prior_p1: float = 1e-4,
    prior_p2: float = 1e-4,
    prior_p12: float = 1e-5,
) -> dict:
    """Coloc posteriors from per-variant log Bayes factors (``coloc.bf_bf``).

    Combines the evidence of two traits over the variants of one region. With
    per-variant Bayes factors ``BF1_i`` and ``BF2_i`` against each trait's null,
    the unnormalised weights of the five hypotheses are

        H0: 1 - p1 - p2 - p12
        H1: p1 * sum_i BF1_i
        H2: p2 * sum_i BF2_i
        H3: p1 * p2 * sum_{i != j} BF1_i * BF2_j
        H4: p12 * sum_i BF1_i * BF2_i

    Everything is accumulated in log space, including the ``i != j`` sum, which
    is built from leave-one-out sums instead of the difference
    ``sum(BF1) * sum(BF2) - sum(BF1 * BF2)`` that cancels catastrophically when one
    variant dominates both traits.

    Use this directly when the Bayes factors come from another source (SuSiE,
    a different variance model per trait, ...).

    Args:
        log_bayes_factors_1: Natural-log Bayes factors of trait 1, per variant.
        log_bayes_factors_2: Natural-log Bayes factors of trait 2, per variant.
        prior_p1: Prior probability a variant is causal for trait 1 only.
        prior_p2: Prior probability a variant is causal for trait 2 only.
        prior_p12: Prior probability a variant is causal for both traits.

    Returns:
        Dictionary with keys status, PP_H0 .. PP_H4, summary and snp_pp_h4, as
        for :func:`colocalization`.
    """
    if len(log_bayes_factors_1) == 0 or len(log_bayes_factors_2) == 0:
        return {"status": "error", "message": "Both trait Bayes-factor vectors are required"}

    if len(log_bayes_factors_1) != len(log_bayes_factors_2):
        return {
            "status": "error",
            "message": (
                f"Bayes-factor lengths differ: trait 1 has {len(log_bayes_factors_1)}, "
                f"trait 2 has {len(log_bayes_factors_2)}"
            ),
        }

    log_bf1 = [float(value) for value in log_bayes_factors_1]
    log_bf2 = [float(value) for value in log_bayes_factors_2]
    if not all(math.isfinite(value) for value in log_bf1 + log_bf2):
        return {"status": "error", "message": "Log Bayes factors must be finite"}

    # Validate priors
    if any(p <= 0 for p in (prior_p1, prior_p2, prior_p12)):
        return {"status": "error", "message": "Prior probabilities must be positive"}

    prior_p0 = 1.0 - prior_p1 - prior_p2 - prior_p12
    if prior_p0 < 0:
        return {"status": "error", "message": "Prior probabilities sum to more than 1"}

    log_bf12 = [a + b for a, b in zip(log_bf1, log_bf2)]

    # Log of each evidence sum; every term stays in log space.
    log_sum_1 = log_sum_exp(log_bf1)
    log_sum_2 = log_sum_exp(log_bf2)
    log_sum_12 = log_sum_exp(log_bf12)
    log_sum_cross = _log_sum_exp_distinct_pairs(log_bf1, log_bf2)

    # Unnormalized log-weights for H0..H4
    log_weights = [
        _log_or_neg_inf(prior_p0),
        _log_or_neg_inf(prior_p1) + log_sum_1,
        _log_or_neg_inf(prior_p2) + log_sum_2,
        _log_or_neg_inf(prior_p1) + _log_or_neg_inf(prior_p2) + log_sum_cross,
        _log_or_neg_inf(prior_p12) + log_sum_12,
    ]

    log_total = log_sum_exp(log_weights)
    if math.isinf(log_total):
        posteriors = [0.2] * 5
    else:
        posteriors = [math.exp(w - log_total) for w in log_weights]

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

    # Which variant is the shared causal one, given H4 (coloc's SNP.PP.H4)
    snp_pp_h4 = [math.exp(value - log_sum_12) for value in log_bf12]

    return {
        "status": "success",
        "PP_H0": float(pp_h0),
        "PP_H1": float(pp_h1),
        "PP_H2": float(pp_h2),
        "PP_H3": float(pp_h3),
        "PP_H4": float(pp_h4),
        "summary": summary,
        "snp_pp_h4": snp_pp_h4,
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


class _WakefieldABF(NamedTuple):
    """Log approximate Bayes factors plus how their scale was fixed."""

    log_bayes_factors: list[float]
    variance_model: str
    prior_variance: float


def _exp_capped(log_value: float) -> float:
    """Exponentiate a log Bayes factor, capped at ``exp(500)`` so it stays finite."""
    return math.exp(min(log_value, _MAX_LOG_BAYES_FACTOR))


def _log_or_neg_inf(value: float) -> float:
    """Natural log of a non-negative number, with ``log(0) = -inf``."""
    return math.log(value) if value > 0 else -math.inf


def _log_add_exp(a: float, b: float) -> float:
    """Return ``log(exp(a) + exp(b))`` without overflow."""
    if a == -math.inf:
        return b
    if b == -math.inf:
        return a
    peak = max(a, b)
    return peak + math.log1p(math.exp(-abs(a - b)))


def _log_sum_exp_distinct_pairs(log_a: Sequence[float], log_b: Sequence[float]) -> float:
    """Return ``log(sum_{i != j} exp(log_a[i] + log_b[j]))`` stably.

    Built from prefix/suffix log-sums of ``log_b`` so that each leave-one-out sum
    is accumulated directly. The algebraic shortcut ``sum(a) * sum(b) - sum(a * b)``
    loses every significant digit when one variant dominates both vectors.

    Returns ``-inf`` when fewer than two variants exist (no distinct pair).
    """
    n = len(log_a)
    if n < 2:
        return -math.inf

    prefix = [-math.inf] * n
    running = -math.inf
    for i in range(n):
        prefix[i] = running
        running = _log_add_exp(running, log_b[i])

    suffix = [-math.inf] * n
    running = -math.inf
    for i in range(n - 1, -1, -1):
        suffix[i] = running
        running = _log_add_exp(running, log_b[i])

    return log_sum_exp(log_a[i] + _log_add_exp(prefix[i], suffix[i]) for i in range(n))


def _finite_values(name: str, values: float | Sequence[float], n_variants: int) -> list[float]:
    """Return ``values`` as ``n_variants`` finite floats (a scalar is broadcast)."""
    if isinstance(values, Real):
        items = [float(cast(float, values))] * n_variants
    else:
        items = [float(v) for v in cast("Sequence[float]", values)]
        if len(items) != n_variants:
            raise ValueError(f"{name} has {len(items)} values but {n_variants} z-scores were given")
    if not all(math.isfinite(v) for v in items):
        raise ValueError(f"{name} must be finite")
    return items


def _effect_variances(
    n_variants: int,
    standard_errors: float | Sequence[float] | None,
    sample_size: float | Sequence[float] | None,
    maf: float | Sequence[float] | None,
    case_fraction: float | None,
) -> tuple[list[float] | None, str]:
    """Resolve the sampling variance ``V`` of each effect estimate.

    Returns ``(variances, variance_model)``; ``variances`` is ``None`` on the
    z-score scale, where ``V = 1`` for every variant.
    """
    if standard_errors is not None and sample_size is not None:
        raise ValueError("Provide standard_errors or sample_size, not both")

    if case_fraction is not None and not 0.0 < case_fraction < 1.0:
        raise ValueError(f"case_fraction must be in (0, 1), got {case_fraction}")

    if standard_errors is not None:
        se = _finite_values("standard_errors", standard_errors, n_variants)
        if any(value <= 0.0 for value in se):
            raise ValueError("standard_errors must be positive")
        return [value * value for value in se], VARIANCE_MODEL_STANDARD_ERRORS

    if sample_size is not None:
        n_obs = _finite_values("sample_size", sample_size, n_variants)
        if any(value <= 0.0 for value in n_obs):
            raise ValueError("sample_size must be positive")
        freqs = _finite_values("maf", maf, n_variants) if maf is not None else None
        if freqs is not None and any(not 0.0 < f < 1.0 for f in freqs):
            raise ValueError("maf must be in (0, 1)")
        outcome_variance = case_fraction * (1.0 - case_fraction) if case_fraction is not None else 1.0
        variances = []
        for j in range(n_variants):
            # Var(genotype) = 2 f (1 - f) under HWE, or 1 for a standardised genotype
            genotype_variance = 2.0 * freqs[j] * (1.0 - freqs[j]) if freqs is not None else 1.0
            variances.append(1.0 / (n_obs[j] * genotype_variance * outcome_variance))
        return variances, VARIANCE_MODEL_SAMPLE_SIZE

    return None, VARIANCE_MODEL_Z_SCALE


def _wakefield_log_abf(
    z_scores: Sequence[float],
    prior_variance: float | None,
    *,
    standard_errors: float | Sequence[float] | None = None,
    sample_size: float | Sequence[float] | None = None,
    maf: float | Sequence[float] | None = None,
    case_fraction: float | None = None,
) -> _WakefieldABF:
    """Evaluate ``log ABF = (log(V / (V + W)) + z**2 * W / (V + W)) / 2`` per variant.

    See :func:`compute_log_bayes_factors` for the meaning and sources of ``V`` and
    ``W``. ``log(V / (V + W))`` is evaluated as ``-log1p(W / V)`` so it is accurate
    whether ``V`` is far smaller or far larger than ``W``.
    """
    if prior_variance is not None and not (math.isfinite(prior_variance) and prior_variance > 0):
        raise ValueError(f"prior_variance must be positive, got {prior_variance}")

    z_values = [float(z) for z in z_scores]
    if not all(math.isfinite(z) for z in z_values):
        raise ValueError("z_scores must be finite")

    variances, variance_model = _effect_variances(len(z_values), standard_errors, sample_size, maf, case_fraction)

    if prior_variance is not None:
        prior = float(prior_variance)
    elif variance_model == VARIANCE_MODEL_Z_SCALE:
        prior = DEFAULT_PRIOR_VARIANCE_Z_SCALE
        logger.debug(
            "No standard_errors or sample_size given: Bayes factors use the z-score scale "
            f"(V = 1, noncentrality prior SD {DEFAULT_PRIOR_SD_Z_SCALE})"
        )
    elif case_fraction is not None:
        prior = DEFAULT_PRIOR_VARIANCE_CASE_CONTROL
    else:
        prior = DEFAULT_PRIOR_VARIANCE_QUANTITATIVE

    log_abf: list[float] = []
    for j, z in enumerate(z_values):
        variance = 1.0 if variances is None else variances[j]
        value = 0.5 * (-math.log1p(prior / variance) + z * z * prior / (variance + prior))
        if not math.isfinite(value):
            raise ValueError(f"Bayes factor for variant {j} is not finite (z={z}, V={variance}, W={prior})")
        log_abf.append(value)

    return _WakefieldABF(log_abf, variance_model, prior)


def _apply_ld_adjustment(
    log_bayes_factors: list[float],
    ld_matrix: list[list[float]],
) -> list[float]:
    """Adjust log Bayes factors using LD structure.

    Applies a penalty to variants in high LD with stronger signals, reducing
    the effective Bayes factor of tagging variants. Works on log Bayes factors,
    so the penalty is added in log space and strong signals cannot overflow.

    Args:
        log_bayes_factors: Raw natural-log Bayes factors.
        ld_matrix: LD correlation matrix.

    Returns:
        Adjusted natural-log Bayes factors.
    """
    if not HAS_NUMPY:
        return log_bayes_factors

    bf_arr = np.array(log_bayes_factors, dtype=float)
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
        adjusted[idx] = bf_arr[idx] + math.log(max(penalty, 0.01))

    return cast(list[float], adjusted.tolist())


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
    log_bayes_factors: list[float],
    prior_probs: list[float],
) -> list[float]:
    """Compute posterior inclusion probabilities from log BFs and priors.

    PIP_i = BF_i * prior_i / sum_j(BF_j * prior_j), evaluated as a softmax of
    ``log BF_i + log prior_i`` so it is exact for any |z|.

    Args:
        log_bayes_factors: Per-variant natural-log Bayes factors.
        prior_probs: Per-variant prior probabilities.

    Returns:
        List of PIPs summing to 1.
    """
    n = len(log_bayes_factors)
    if n == 0:
        return []

    log_posts = [lbf + _log_or_neg_inf(pr) for lbf, pr in zip(log_bayes_factors, prior_probs)]
    log_total = log_sum_exp(log_posts)

    if math.isinf(log_total):
        return [1.0 / n] * n
    return [math.exp(lp - log_total) for lp in log_posts]


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
    # erfc avoids the cancellation in 1 - cdf(|z|), which reaches 0 beyond |z| = 8.3.
    return math.erfc(abs(z) / math.sqrt(2.0))


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
