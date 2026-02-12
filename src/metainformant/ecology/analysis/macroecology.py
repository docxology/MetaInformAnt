"""Macroecological analysis methods for species abundance distributions and scaling laws.

This module implements classical macroecological models including species abundance
distributions (SAD), species-area relationships (SAR), distance-decay relationships,
metabolic scaling, and related macroecological patterns. All implementations are
pure Python with no external numerical dependencies.

References:
    - Fisher, R.A., Corbet, A.S. & Williams, C.B. (1943). The relation between
      the number of species and the number of individuals.
    - Preston, F.W. (1948). The commonness, and rarity, of species.
    - MacArthur, R.H. (1957). On the relative abundance of bird species.
    - Arrhenius, O. (1921). Species and area. Journal of Ecology.
    - Nekola, J.C. & White, P.S. (1999). The distance decay of similarity.
    - Kleiber, M. (1932). Body size and metabolism.
    - Taylor, L.R. (1961). Aggregation, variance and the mean.
"""

from __future__ import annotations

import math
import random
import statistics
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.data import validation
from metainformant.core.utils import errors
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


# ---------------------------------------------------------------------------
# Helper / utility functions
# ---------------------------------------------------------------------------


def simple_linear_regression(x: List[float], y: List[float]) -> Dict[str, float]:
    """Ordinary least-squares linear regression (y = slope*x + intercept).

    Args:
        x: Independent variable values.
        y: Dependent variable values (same length as *x*).

    Returns:
        Dictionary with keys ``slope``, ``intercept``, and ``r_squared``.

    Raises:
        ValueError: If inputs have different lengths or fewer than 2 points.

    Example:
        >>> result = simple_linear_regression([1, 2, 3], [2, 4, 6])
        >>> abs(result["slope"] - 2.0) < 1e-9
        True
    """
    if len(x) != len(y):
        raise ValueError(f"x and y must have same length, got {len(x)} and {len(y)}")
    n = len(x)
    if n < 2:
        raise ValueError("Need at least 2 data points for linear regression")

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    ss_xx = sum((xi - mean_x) ** 2 for xi in x)
    ss_xy = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    ss_yy = sum((yi - mean_y) ** 2 for yi in y)

    if ss_xx == 0.0:
        return {"slope": 0.0, "intercept": mean_y, "r_squared": 0.0}

    slope = ss_xy / ss_xx
    intercept = mean_y - slope * mean_x

    ss_res = sum((yi - (slope * xi + intercept)) ** 2 for xi, yi in zip(x, y))
    r_squared = 1.0 - (ss_res / ss_yy) if ss_yy > 0.0 else 0.0

    return {"slope": slope, "intercept": intercept, "r_squared": r_squared}


def chi_squared_gof(observed: List[float], expected: List[float]) -> Dict[str, float]:
    """Chi-squared goodness-of-fit statistic.

    Args:
        observed: Observed frequency counts.
        expected: Expected frequency counts (same length).

    Returns:
        Dictionary with ``chi_squared`` value and ``degrees_of_freedom``.

    Raises:
        ValueError: If inputs differ in length or contain negative values.

    Example:
        >>> result = chi_squared_gof([10, 20, 30], [15, 15, 30])
        >>> result["degrees_of_freedom"]
        2
    """
    if len(observed) != len(expected):
        raise ValueError("observed and expected must have the same length")

    chi2 = 0.0
    k = 0
    for obs, exp in zip(observed, expected):
        if exp > 0.0:
            chi2 += (obs - exp) ** 2 / exp
            k += 1

    df = max(k - 1, 1)
    return {"chi_squared": chi2, "degrees_of_freedom": df}


def aic_from_gof(chi_squared: float, n_params: int, n_obs: int) -> float:
    """Approximate AIC from chi-squared goodness-of-fit.

    Uses the relationship AIC = chi2 + 2*k where k is the number of estimated
    parameters. This is an approximation suitable for model ranking.

    Args:
        chi_squared: Chi-squared statistic.
        n_params: Number of estimated parameters in the model.
        n_obs: Number of observations (used for AICc correction).

    Returns:
        AIC value (lower is better).
    """
    aic = chi_squared + 2.0 * n_params
    # Small-sample correction (AICc)
    if n_obs > n_params + 1:
        aic += (2.0 * n_params * (n_params + 1)) / (n_obs - n_params - 1)
    return aic


def bootstrap_ci(
    x: List[float],
    y: List[float],
    statistic_fn: Any,
    n_bootstrap: int = 1000,
    confidence: float = 0.95,
    seed: int = 42,
) -> Tuple[float, float]:
    """Bootstrap confidence interval for a regression statistic.

    Args:
        x: Independent variable values.
        y: Dependent variable values.
        statistic_fn: Callable(x_sample, y_sample) -> float returning the statistic.
        n_bootstrap: Number of bootstrap resamples.
        confidence: Confidence level (e.g. 0.95 for 95 %).
        seed: Random seed for reproducibility.

    Returns:
        Tuple of (lower_bound, upper_bound).
    """
    rng = random.Random(seed)
    n = len(x)
    if n < 2:
        val = statistic_fn(x, y)
        return (val, val)

    stats: List[float] = []
    for _ in range(n_bootstrap):
        indices = [rng.randint(0, n - 1) for _ in range(n)]
        x_sample = [x[i] for i in indices]
        y_sample = [y[i] for i in indices]
        try:
            stats.append(statistic_fn(x_sample, y_sample))
        except (ValueError, ZeroDivisionError):
            continue

    if not stats:
        return (0.0, 0.0)

    stats.sort()
    alpha = 1.0 - confidence
    lower_idx = max(0, int(math.floor(alpha / 2.0 * len(stats))))
    upper_idx = min(len(stats) - 1, int(math.ceil((1.0 - alpha / 2.0) * len(stats))) - 1)
    return (stats[lower_idx], stats[upper_idx])


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _sorted_abundances(abundances: List[float]) -> List[float]:
    """Return positive abundances sorted in descending order."""
    return sorted([a for a in abundances if a > 0], reverse=True)


def _validate_abundances(abundances: List[float], name: str = "abundances") -> List[float]:
    """Validate and clean an abundance vector, returning positive values sorted descending."""
    validation.validate_not_empty(abundances, name)
    pos = _sorted_abundances(abundances)
    if not pos:
        raise ValueError(f"All values in {name} are zero or negative")
    return pos


# ---------------------------------------------------------------------------
# 1. Fisher's Log-Series Distribution
# ---------------------------------------------------------------------------


def fit_logseries(abundances: List[float]) -> Dict[str, Any]:
    """Fit Fisher's log-series distribution to abundance data.

    Uses Newton's method to estimate Fisher's alpha from the relationship
    ``S = alpha * ln(1 + N/alpha)`` where S is species richness and N is
    the total number of individuals.

    Args:
        abundances: Species abundance counts (positive integers expected).

    Returns:
        Dictionary with:
            - ``alpha``: Fisher's alpha diversity parameter.
            - ``x``: The log-series x parameter (N / (N + alpha)).
            - ``expected_frequencies``: Expected number of species with 1, 2, ... individuals.
            - ``goodness_of_fit``: Chi-squared GOF result dict.
            - ``fisher_diversity``: Fisher's alpha (identical to alpha, for clarity).

    Raises:
        ValueError: If abundances are empty or all zero.

    Example:
        >>> result = fit_logseries([20, 15, 12, 8, 5, 3, 2, 1, 1, 1])
        >>> result["alpha"] > 0
        True
    """
    pos = _validate_abundances(abundances)
    s = len(pos)
    n = sum(pos)

    logger.info(f"Fitting Fisher's log-series: S={s}, N={n}")

    if s == 1:
        return {
            "alpha": 0.0,
            "x": 0.0,
            "expected_frequencies": [float(s)],
            "goodness_of_fit": {"chi_squared": 0.0, "degrees_of_freedom": 0},
            "fisher_diversity": 0.0,
        }

    # Newton's method to solve S = alpha * ln(1 + N/alpha) for alpha
    alpha = s / math.log(1.0 + n / s)  # initial guess
    for _ in range(200):
        ratio = n / alpha
        f_val = alpha * math.log(1.0 + ratio) - s
        # derivative: ln(1 + N/alpha) - N/(alpha + N)
        f_prime = math.log(1.0 + ratio) - ratio / (1.0 + ratio)
        if abs(f_prime) < 1e-15:
            break
        alpha_new = alpha - f_val / f_prime
        if alpha_new <= 0:
            alpha_new = alpha / 2.0
        if abs(alpha_new - alpha) < 1e-12 * alpha:
            alpha = alpha_new
            break
        alpha = alpha_new

    x = n / (n + alpha)

    # Expected frequencies: E(n_r) = alpha * x^r / r  for r = 1, 2, ...
    max_abund = int(max(pos))
    expected: List[float] = []
    for r in range(1, max_abund + 1):
        expected.append(alpha * (x**r) / r)

    # Build observed frequency histogram
    from collections import Counter

    counts = Counter(int(a) for a in pos)
    observed: List[float] = []
    exp_trimmed: List[float] = []
    for r in range(1, max_abund + 1):
        observed.append(float(counts.get(r, 0)))
        if r <= len(expected):
            exp_trimmed.append(expected[r - 1])
        else:
            exp_trimmed.append(0.0)

    gof = chi_squared_gof(observed, exp_trimmed)

    return {
        "alpha": alpha,
        "x": x,
        "expected_frequencies": expected,
        "goodness_of_fit": gof,
        "fisher_diversity": alpha,
    }


# ---------------------------------------------------------------------------
# 2. Preston's Lognormal Distribution
# ---------------------------------------------------------------------------


def fit_lognormal(abundances: List[float]) -> Dict[str, Any]:
    """Fit Preston's lognormal distribution to abundance data.

    Uses method of moments on log-transformed abundances to estimate the
    parameters of the underlying normal distribution on the log scale.

    Args:
        abundances: Species abundance counts (positive values).

    Returns:
        Dictionary with:
            - ``mu``: Mean of log-abundances.
            - ``sigma``: Standard deviation of log-abundances.
            - ``expected_frequencies``: Expected species per Preston octave.
            - ``goodness_of_fit``: Chi-squared GOF result dict.
            - ``s_star``: Estimated total species richness (including unseen).

    Raises:
        ValueError: If abundances are empty or all zero.

    Example:
        >>> result = fit_lognormal([100, 50, 25, 12, 6, 3, 1, 1])
        >>> result["mu"] > 0
        True
    """
    pos = _validate_abundances(abundances)
    s = len(pos)

    log_abund = [math.log(a) for a in pos]
    mu = statistics.mean(log_abund)
    sigma = statistics.stdev(log_abund) if s > 1 else 0.0

    logger.info(f"Fitting lognormal: mu={mu:.3f}, sigma={sigma:.3f}, S={s}")

    # Preston octaves: bin species into log2 abundance classes
    max_log2 = int(math.log2(max(pos))) + 1 if max(pos) > 0 else 1
    n_octaves = max_log2 + 1

    observed_octaves: List[float] = [0.0] * n_octaves
    for a in pos:
        octave = int(math.log2(a)) if a >= 1 else 0
        octave = min(octave, n_octaves - 1)
        observed_octaves[octave] += 1.0

    # Expected frequencies per octave from lognormal CDF
    # Each octave r spans [2^r, 2^(r+1)) on the abundance axis
    expected_octaves: List[float] = []
    for r in range(n_octaves):
        lower = math.log(2**r) if r > 0 else math.log(0.5)
        upper = math.log(2 ** (r + 1))
        if sigma > 0:
            p_lower = _normal_cdf((lower - mu) / sigma)
            p_upper = _normal_cdf((upper - mu) / sigma)
            expected_octaves.append(s * (p_upper - p_lower))
        else:
            # All species at the same abundance
            if lower <= mu < upper:
                expected_octaves.append(float(s))
            else:
                expected_octaves.append(0.0)

    gof = (
        chi_squared_gof(observed_octaves, expected_octaves)
        if n_octaves > 1
        else {
            "chi_squared": 0.0,
            "degrees_of_freedom": 0,
        }
    )

    # S* estimate (Veil line correction): total species including those behind the veil
    if sigma > 0:
        # Fraction of the normal curve visible (above abundance = 0.5)
        veil_z = (math.log(0.5) - mu) / sigma
        fraction_visible = 1.0 - _normal_cdf(veil_z)
        s_star = s / fraction_visible if fraction_visible > 0 else float(s)
    else:
        s_star = float(s)

    return {
        "mu": mu,
        "sigma": sigma,
        "expected_frequencies": expected_octaves,
        "goodness_of_fit": gof,
        "s_star": s_star,
    }


def _normal_cdf(z: float) -> float:
    """Approximate the standard normal CDF using the Abramowitz & Stegun formula.

    Accurate to about 1e-7 for all z.

    Args:
        z: Standard normal deviate.

    Returns:
        Probability P(Z <= z).
    """
    if z < -8.0:
        return 0.0
    if z > 8.0:
        return 1.0
    # Horner form of the rational approximation
    a1, a2, a3, a4, a5 = 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429
    p = 0.3275911
    sign = 1.0 if z >= 0 else -1.0
    t = 1.0 / (1.0 + p * abs(z))
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-z * z / 2.0)
    return 0.5 * (1.0 + sign * y)


# ---------------------------------------------------------------------------
# 3. MacArthur's Broken Stick Model
# ---------------------------------------------------------------------------


def fit_broken_stick(abundances: List[float]) -> Dict[str, Any]:
    """Fit MacArthur's broken stick model to abundance data.

    The expected proportional abundance for species of rank *k* out of *S*
    species is ``(N/S) * sum(1/i for i in range(k, S+1))``.

    Args:
        abundances: Species abundance counts (positive values).

    Returns:
        Dictionary with:
            - ``expected_abundances``: Predicted abundance at each rank.
            - ``observed_abundances``: Sorted observed abundances (descending).
            - ``goodness_of_fit``: Chi-squared GOF result dict.

    Raises:
        ValueError: If abundances are empty or all zero.

    Example:
        >>> result = fit_broken_stick([50, 30, 15, 5])
        >>> len(result["expected_abundances"]) == 4
        True
    """
    pos = _validate_abundances(abundances)
    s = len(pos)
    n = sum(pos)

    logger.info(f"Fitting broken stick model: S={s}, N={n}")

    expected: List[float] = []
    for k in range(1, s + 1):
        e_k = (n / s) * sum(1.0 / i for i in range(k, s + 1))
        expected.append(e_k)

    gof = chi_squared_gof(pos, expected)

    return {
        "expected_abundances": expected,
        "observed_abundances": list(pos),
        "goodness_of_fit": gof,
    }


# ---------------------------------------------------------------------------
# 4. Motomura's Geometric Series
# ---------------------------------------------------------------------------


def fit_geometric_series(abundances: List[float]) -> Dict[str, Any]:
    """Fit Motomura's geometric series (niche preemption) to abundance data.

    The model assumes each species pre-empts a constant proportion *k* of
    remaining resources. Expected abundance of the species at rank *r* is
    ``N * k * (1-k)^(r-1)`` (with a normalising constant).

    The parameter *k* is estimated by minimising the sum of squared residuals
    between observed and predicted log-abundances via a grid + refinement search.

    Args:
        abundances: Species abundance counts (positive values).

    Returns:
        Dictionary with:
            - ``k``: Estimated dominance (preemption) proportion.
            - ``expected_abundances``: Predicted abundances at each rank.
            - ``goodness_of_fit``: Chi-squared GOF result dict.

    Raises:
        ValueError: If abundances are empty or all zero.

    Example:
        >>> result = fit_geometric_series([100, 50, 25, 12, 6])
        >>> 0.0 < result["k"] < 1.0
        True
    """
    pos = _validate_abundances(abundances)
    s = len(pos)
    n = sum(pos)

    logger.info(f"Fitting geometric series: S={s}, N={n}")

    if s == 1:
        return {
            "k": 1.0,
            "expected_abundances": [float(n)],
            "goodness_of_fit": {"chi_squared": 0.0, "degrees_of_freedom": 0},
        }

    def _expected_geometric(k_val: float) -> List[float]:
        """Compute expected abundances for a given k."""
        if k_val <= 0.0 or k_val >= 1.0:
            return [n / s] * s
        # Normalising constant so expected total = N
        raw = [k_val * ((1.0 - k_val) ** (r - 1)) for r in range(1, s + 1)]
        raw_sum = sum(raw)
        if raw_sum == 0:
            return [n / s] * s
        return [n * r / raw_sum for r in raw]

    def _ss_residuals(k_val: float) -> float:
        exp = _expected_geometric(k_val)
        return sum((math.log(o + 1) - math.log(e + 1)) ** 2 for o, e in zip(pos, exp))

    # Coarse grid search
    best_k = 0.5
    best_ss = _ss_residuals(0.5)
    for k_int in range(1, 100):
        k_try = k_int / 100.0
        ss = _ss_residuals(k_try)
        if ss < best_ss:
            best_ss = ss
            best_k = k_try

    # Golden-section refinement
    lo = max(0.01, best_k - 0.02)
    hi = min(0.99, best_k + 0.02)
    gr = (math.sqrt(5) + 1) / 2
    for _ in range(60):
        c = hi - (hi - lo) / gr
        d = lo + (hi - lo) / gr
        if _ss_residuals(c) < _ss_residuals(d):
            hi = d
        else:
            lo = c
    best_k = (lo + hi) / 2.0

    expected = _expected_geometric(best_k)
    gof = chi_squared_gof(pos, expected)

    return {
        "k": best_k,
        "expected_abundances": expected,
        "goodness_of_fit": gof,
    }


# ---------------------------------------------------------------------------
# 5. Compare SAD models
# ---------------------------------------------------------------------------


def compare_sad_models(abundances: List[float]) -> Dict[str, Dict[str, Any]]:
    """Fit and compare all four species abundance distribution models.

    Models compared: log-series, lognormal, broken stick, geometric series.
    Ranking is by AIC (lower is better).

    Args:
        abundances: Species abundance counts (positive values).

    Returns:
        Dictionary keyed by model name, each containing:
            - ``parameters``: Fitted model parameters.
            - ``gof``: Goodness-of-fit chi-squared dict.
            - ``aic``: AIC value.
            - ``rank``: Integer rank (1 = best).

    Raises:
        ValueError: If abundances are empty or all zero.

    Example:
        >>> result = compare_sad_models([50, 30, 20, 10, 5, 3, 1, 1])
        >>> "logseries" in result
        True
    """
    pos = _validate_abundances(abundances)
    n_obs = len(pos)

    logger.info(f"Comparing SAD models for {n_obs} species")

    # Fit each model
    ls = fit_logseries(pos)
    ln = fit_lognormal(pos)
    bs = fit_broken_stick(pos)
    gs = fit_geometric_series(pos)

    models: Dict[str, Dict[str, Any]] = {
        "logseries": {
            "parameters": {"alpha": ls["alpha"], "x": ls["x"]},
            "gof": ls["goodness_of_fit"],
            "aic": aic_from_gof(ls["goodness_of_fit"]["chi_squared"], 1, n_obs),
        },
        "lognormal": {
            "parameters": {"mu": ln["mu"], "sigma": ln["sigma"]},
            "gof": ln["goodness_of_fit"],
            "aic": aic_from_gof(ln["goodness_of_fit"]["chi_squared"], 2, n_obs),
        },
        "broken_stick": {
            "parameters": {},
            "gof": bs["goodness_of_fit"],
            "aic": aic_from_gof(bs["goodness_of_fit"]["chi_squared"], 0, n_obs),
        },
        "geometric_series": {
            "parameters": {"k": gs["k"]},
            "gof": gs["goodness_of_fit"],
            "aic": aic_from_gof(gs["goodness_of_fit"]["chi_squared"], 1, n_obs),
        },
    }

    # Rank by AIC
    ranked = sorted(models.items(), key=lambda item: item[1]["aic"])
    for rank, (name, _) in enumerate(ranked, 1):
        models[name]["rank"] = rank

    return models


# ---------------------------------------------------------------------------
# 6. Species-Area Relationship (Power Law / Arrhenius)
# ---------------------------------------------------------------------------


def species_area_power(
    areas: List[float],
    species_counts: List[float],
    n_bootstrap: int = 1000,
) -> Dict[str, Any]:
    """Fit the Arrhenius power-law species-area relationship with bootstrap CI.

    Model: ``S = c * A^z`` (linearised as ``log(S) = log(c) + z * log(A)``).

    Args:
        areas: Area sizes (must be positive).
        species_counts: Number of species observed at each area.
        n_bootstrap: Number of bootstrap resamples for confidence intervals.

    Returns:
        Dictionary with:
            - ``c``: Intercept on the natural scale.
            - ``z``: Power-law exponent.
            - ``r_squared``: Coefficient of determination.
            - ``ci_z``: 95 % bootstrap confidence interval for z.
            - ``predicted_species``: Predicted S at each observed area.

    Raises:
        ValueError: If inputs differ in length, have < 2 points, or contain non-positive values.

    Example:
        >>> result = species_area_power([1, 10, 100, 1000], [10, 30, 80, 200])
        >>> result["z"] > 0
        True
    """
    if len(areas) != len(species_counts):
        raise ValueError("areas and species_counts must have the same length")
    if len(areas) < 2:
        raise ValueError("Need at least 2 data points")
    if any(a <= 0 for a in areas) or any(s <= 0 for s in species_counts):
        raise ValueError("All areas and species counts must be positive")

    log_a = [math.log(a) for a in areas]
    log_s = [math.log(s) for s in species_counts]

    reg = simple_linear_regression(log_a, log_s)
    z = reg["slope"]
    c = math.exp(reg["intercept"])
    r_sq = reg["r_squared"]

    # Bootstrap CI for z
    def _slope_fn(xb: List[float], yb: List[float]) -> float:
        r = simple_linear_regression(xb, yb)
        return r["slope"]

    ci_z = bootstrap_ci(log_a, log_s, _slope_fn, n_bootstrap=n_bootstrap)

    predicted = [c * (a**z) for a in areas]

    logger.info(f"SAR power law: c={c:.4f}, z={z:.4f}, R2={r_sq:.4f}")

    return {
        "c": c,
        "z": z,
        "r_squared": r_sq,
        "ci_z": ci_z,
        "predicted_species": predicted,
    }


# ---------------------------------------------------------------------------
# 7. Species-Area Relationship (Logarithmic / Gleason)
# ---------------------------------------------------------------------------


def species_area_logarithmic(
    areas: List[float],
    species_counts: List[float],
) -> Dict[str, Any]:
    """Fit Gleason's logarithmic species-area model.

    Model: ``S = a + b * ln(A)``.

    Args:
        areas: Area sizes (positive values).
        species_counts: Observed species counts.

    Returns:
        Dictionary with:
            - ``a``: Intercept.
            - ``b``: Slope coefficient for ln(area).
            - ``r_squared``: Coefficient of determination.
            - ``predicted_species``: Predicted S at each observed area.

    Raises:
        ValueError: If inputs differ in length or have < 2 points.

    Example:
        >>> result = species_area_logarithmic([1, 10, 100], [5, 15, 25])
        >>> result["b"] > 0
        True
    """
    if len(areas) != len(species_counts):
        raise ValueError("areas and species_counts must have the same length")
    if len(areas) < 2:
        raise ValueError("Need at least 2 data points")
    if any(a <= 0 for a in areas):
        raise ValueError("All areas must be positive")

    ln_a = [math.log(a) for a in areas]
    s_vals = [float(s) for s in species_counts]

    reg = simple_linear_regression(ln_a, s_vals)
    a = reg["intercept"]
    b = reg["slope"]
    r_sq = reg["r_squared"]

    predicted = [a + b * math.log(area) for area in areas]

    logger.info(f"SAR logarithmic: a={a:.4f}, b={b:.4f}, R2={r_sq:.4f}")

    return {
        "a": a,
        "b": b,
        "r_squared": r_sq,
        "predicted_species": predicted,
    }


# ---------------------------------------------------------------------------
# 8. Distance-Decay of Similarity
# ---------------------------------------------------------------------------


def distance_decay(
    distances: List[float],
    similarities: List[float],
) -> Dict[str, Any]:
    """Fit exponential and power-law distance-decay models.

    Exponential: ``similarity = a * exp(-b * distance)``
    Power-law:   ``similarity = a * distance^(-b)``

    Follows the Nekola & White (1999) approach.

    Args:
        distances: Pairwise geographic distances (positive).
        similarities: Corresponding community similarity values (0-1 range typical).

    Returns:
        Dictionary with:
            - ``exponential``: {a, b, r_squared} for the exponential model.
            - ``power_law``: {a, b, r_squared} for the power-law model.
            - ``best_model``: Name of the model with higher R-squared.

    Raises:
        ValueError: If inputs differ in length or have < 2 points.

    Example:
        >>> result = distance_decay([1, 5, 10, 20], [0.9, 0.6, 0.3, 0.1])
        >>> result["best_model"] in ("exponential", "power_law")
        True
    """
    if len(distances) != len(similarities):
        raise ValueError("distances and similarities must have the same length")
    if len(distances) < 2:
        raise ValueError("Need at least 2 data points")

    # Filter out non-positive distances/similarities for log transforms
    valid = [(d, s) for d, s in zip(distances, similarities) if d > 0 and s > 0]
    if len(valid) < 2:
        raise ValueError("Need at least 2 points with positive distance and similarity")

    dist_v = [v[0] for v in valid]
    sim_v = [v[1] for v in valid]

    # --- Exponential model: ln(sim) = ln(a) - b * distance ---
    log_sim = [math.log(s) for s in sim_v]
    exp_reg = simple_linear_regression(dist_v, log_sim)
    exp_a = math.exp(exp_reg["intercept"])
    exp_b = -exp_reg["slope"]  # negate because model has -b*d
    exp_r2 = exp_reg["r_squared"]

    # --- Power-law model: ln(sim) = ln(a) - b * ln(distance) ---
    log_dist = [math.log(d) for d in dist_v]
    pow_reg = simple_linear_regression(log_dist, log_sim)
    pow_a = math.exp(pow_reg["intercept"])
    pow_b = -pow_reg["slope"]
    pow_r2 = pow_reg["r_squared"]

    best = "exponential" if exp_r2 >= pow_r2 else "power_law"

    logger.info(f"Distance decay: best model = {best}")

    return {
        "exponential": {"a": exp_a, "b": exp_b, "r_squared": exp_r2},
        "power_law": {"a": pow_a, "b": pow_b, "r_squared": pow_r2},
        "best_model": best,
    }


# ---------------------------------------------------------------------------
# 9. Occupancy-Frequency Distribution
# ---------------------------------------------------------------------------


def occupancy_frequency(
    presence_absence_matrix: List[List[int]],
) -> Dict[str, Any]:
    """Analyse the occupancy-frequency distribution of a species-by-site matrix.

    Species are classified as:
        - **core**: present at > 66 % of sites
        - **common**: present at 33--66 % of sites
        - **satellite**: present at < 33 % of sites

    A bimodality index is computed as the ratio of species in the two
    extreme thirds (core + satellite) to total species.

    Args:
        presence_absence_matrix: Binary matrix where rows are sites and
            columns are species (1 = present, 0 = absent).

    Returns:
        Dictionary with:
            - ``core_species``: Number of core species.
            - ``common_species``: Number of common species.
            - ``satellite_species``: Number of satellite species.
            - ``occupancy_distribution``: List of occupancy fractions per species.
            - ``bimodality_index``: Ratio of extreme-class species to total.

    Raises:
        ValueError: If the matrix is empty.

    Example:
        >>> mat = [[1, 0, 1], [1, 1, 0], [1, 1, 0]]
        >>> result = occupancy_frequency(mat)
        >>> result["core_species"] >= 0
        True
    """
    validation.validate_not_empty(presence_absence_matrix, "presence_absence_matrix")

    n_sites = len(presence_absence_matrix)
    if n_sites == 0:
        raise ValueError("presence_absence_matrix must not be empty")

    n_species = len(presence_absence_matrix[0]) if presence_absence_matrix[0] else 0
    if n_species == 0:
        return {
            "core_species": 0,
            "common_species": 0,
            "satellite_species": 0,
            "occupancy_distribution": [],
            "bimodality_index": 0.0,
        }

    # Calculate occupancy fraction for each species
    occupancy: List[float] = []
    for sp in range(n_species):
        n_present = sum(1 for site in range(n_sites) if presence_absence_matrix[site][sp])
        occupancy.append(n_present / n_sites)

    core = sum(1 for o in occupancy if o > 2 / 3)
    satellite = sum(1 for o in occupancy if o < 1 / 3)
    common = n_species - core - satellite

    bimodality = (core + satellite) / n_species if n_species > 0 else 0.0

    logger.info(f"Occupancy: {core} core, {common} common, {satellite} satellite species")

    return {
        "core_species": core,
        "common_species": common,
        "satellite_species": satellite,
        "occupancy_distribution": occupancy,
        "bimodality_index": bimodality,
    }


# ---------------------------------------------------------------------------
# 10. Metabolic Theory of Ecology (MTE) scaling
# ---------------------------------------------------------------------------


def metabolic_scaling(
    body_masses: List[float],
    metabolic_rates: List[float],
) -> Dict[str, Any]:
    """Fit allometric metabolic scaling and test against Kleiber's law.

    Model: ``B = B0 * M^b``, linearised as ``log(B) = log(B0) + b * log(M)``.

    Kleiber's law predicts b ~ 0.75 (three-quarter power scaling).

    Args:
        body_masses: Body mass values (positive).
        metabolic_rates: Corresponding metabolic rate values (positive).

    Returns:
        Dictionary with:
            - ``b0``: Scaling constant (natural scale).
            - ``b_exponent``: Allometric exponent.
            - ``r_squared``: Coefficient of determination.
            - ``kleiber_deviation``: Absolute difference of b from 0.75.

    Raises:
        ValueError: If inputs differ in length, have < 2 points, or are non-positive.

    Example:
        >>> result = metabolic_scaling([1, 10, 100, 1000], [0.5, 3.0, 17.0, 100.0])
        >>> result["b_exponent"] > 0
        True
    """
    if len(body_masses) != len(metabolic_rates):
        raise ValueError("body_masses and metabolic_rates must have the same length")
    if len(body_masses) < 2:
        raise ValueError("Need at least 2 data points")
    if any(m <= 0 for m in body_masses) or any(r <= 0 for r in metabolic_rates):
        raise ValueError("All body masses and metabolic rates must be positive")

    log_m = [math.log(m) for m in body_masses]
    log_b = [math.log(r) for r in metabolic_rates]

    reg = simple_linear_regression(log_m, log_b)
    b_exp = reg["slope"]
    b0 = math.exp(reg["intercept"])
    r_sq = reg["r_squared"]

    kleiber_dev = abs(b_exp - 0.75)

    logger.info(f"Metabolic scaling: B0={b0:.4f}, b={b_exp:.4f}, deviation from 0.75 = {kleiber_dev:.4f}")

    return {
        "b0": b0,
        "b_exponent": b_exp,
        "r_squared": r_sq,
        "kleiber_deviation": kleiber_dev,
    }


# ---------------------------------------------------------------------------
# 11. Weighted Endemism Index
# ---------------------------------------------------------------------------


def endemism_index(
    species_ranges: List[float],
    area_of_interest: float,
) -> Dict[str, Any]:
    """Compute weighted and corrected weighted endemism indices.

    Weighted endemism (WE) sums the inverse range sizes of all species
    present. Corrected WE divides by species count to remove the effect
    of richness.

    A species is considered *endemic* if its range size is smaller than
    the area of interest.

    Args:
        species_ranges: Range sizes for each species (positive values, same
            units as *area_of_interest*).
        area_of_interest: Size of the focal area.

    Returns:
        Dictionary with:
            - ``weighted_endemism``: Sum of 1/range_size across species.
            - ``corrected_weighted_endemism``: WE / species count.
            - ``n_endemic_species``: Count of species with range < area_of_interest.

    Raises:
        ValueError: If species_ranges is empty or area_of_interest is non-positive.

    Example:
        >>> result = endemism_index([10, 50, 200, 500], 100)
        >>> result["n_endemic_species"]
        2
    """
    validation.validate_not_empty(species_ranges, "species_ranges")
    if area_of_interest <= 0:
        raise ValueError("area_of_interest must be positive")

    positive_ranges = [r for r in species_ranges if r > 0]
    if not positive_ranges:
        raise ValueError("No positive range sizes provided")

    we = sum(1.0 / r for r in positive_ranges)
    n_species = len(positive_ranges)
    cwe = we / n_species if n_species > 0 else 0.0
    n_endemic = sum(1 for r in positive_ranges if r < area_of_interest)

    logger.info(f"Endemism: WE={we:.4f}, CWE={cwe:.4f}, {n_endemic}/{n_species} endemic")

    return {
        "weighted_endemism": we,
        "corrected_weighted_endemism": cwe,
        "n_endemic_species": n_endemic,
    }


# ---------------------------------------------------------------------------
# 12. Taylor's Power Law
# ---------------------------------------------------------------------------


def taylors_power_law(
    means: List[float],
    variances: List[float],
) -> Dict[str, Any]:
    """Fit Taylor's power law relating variance to the mean.

    Model: ``V = a * M^b``, linearised as ``log(V) = log(a) + b * log(M)``.

    The exponent *b* indicates spatial aggregation:
        - b ~ 1: Poisson (random) distribution
        - b ~ 2: negative binomial (aggregated) distribution
        - b > 2: highly aggregated

    Args:
        means: Sample means (positive values).
        variances: Corresponding sample variances (positive values).

    Returns:
        Dictionary with:
            - ``a``: Scaling constant (natural scale).
            - ``b``: Power-law exponent.
            - ``r_squared``: Coefficient of determination.

    Raises:
        ValueError: If inputs differ in length, have < 2 points, or are non-positive.

    Example:
        >>> result = taylors_power_law([1, 5, 10, 50], [1.2, 28, 110, 2600])
        >>> result["b"] > 1
        True
    """
    if len(means) != len(variances):
        raise ValueError("means and variances must have the same length")
    if len(means) < 2:
        raise ValueError("Need at least 2 data points")

    # Filter pairs where both are positive
    valid = [(m, v) for m, v in zip(means, variances) if m > 0 and v > 0]
    if len(valid) < 2:
        raise ValueError("Need at least 2 pairs with positive mean and variance")

    log_m = [math.log(v[0]) for v in valid]
    log_v = [math.log(v[1]) for v in valid]

    reg = simple_linear_regression(log_m, log_v)
    a = math.exp(reg["intercept"])
    b = reg["slope"]
    r_sq = reg["r_squared"]

    logger.info(f"Taylor's power law: a={a:.4f}, b={b:.4f}, R2={r_sq:.4f}")

    return {
        "a": a,
        "b": b,
        "r_squared": r_sq,
    }
