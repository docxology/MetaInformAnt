"""Multi-omic survival analysis.

Implements Cox proportional hazards regression, Kaplan-Meier survival
curve estimation, log-rank tests, regularised multi-omic survival models,
and risk stratification.  All methods operate on plain Python lists with
optional numpy acceleration.

Core algorithms:
    - Cox PH: Newton-Raphson on the partial likelihood with Breslow ties.
    - Kaplan-Meier: product-limit estimator with Greenwood variance.
    - Log-rank: Mantel-Haenszel chi-squared comparing survival curves.
    - Lasso-Cox: coordinate descent on the L1-penalised partial likelihood.
    - Concordance index: Harrell's C for discrimination evaluation.
"""

from __future__ import annotations

import math
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional scientific dependencies
np: Any  # None when numpy is unavailable (optional dependency)
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None

sp_stats: Any  # None when scipy is unavailable (scipy is a declared dependency, so this is a safety net)
try:
    from scipy import stats as sp_stats
except ImportError:
    sp_stats = None


# ---------------------------------------------------------------------------
# Cox Proportional Hazards Regression
# ---------------------------------------------------------------------------


_COX_MAX_ITERATIONS = 50
_COX_MAX_HALVINGS = 30


def _breslow_partial_likelihood(time: Any, event: Any, covariates: Any, beta: Any) -> tuple[float, Any, Any]:
    """Log partial likelihood, its gradient and the information matrix at ``beta`` (Breslow ties).

    The three inputs are sorted by ascending time. For every distinct event time ``tau`` with ``d`` events
    and risk set ``R`` (everyone observed for at least ``tau``, tied subjects included), with
    ``S0 = sum_R exp(eta)``, ``S1 = sum_R x exp(eta)`` and ``S2 = sum_R x x' exp(eta)``::

        loglik   += sum_events eta - d log S0
        gradient += sum_events x - d S1 / S0
        info     += d (S2 / S0 - (S1 / S0)(S1 / S0)')

    The risk-set sums are accumulated while walking the distinct times from the latest to the earliest, so
    only ``p x p`` matrices are formed.

    Returns:
        ``(loglik, gradient, information)``; the information matrix is minus the Hessian.
    """
    n, p = covariates.shape
    eta = covariates @ beta
    shift = float(eta.max())  # exp(eta - shift) cannot overflow; the shift is added back into the likelihood
    weight = np.exp(eta - shift)

    loglik = 0.0
    gradient = np.zeros(p)
    information = np.zeros((p, p))
    s0 = 0.0
    s1 = np.zeros(p)
    s2 = np.zeros((p, p))

    stop = n
    while stop > 0:
        start = stop - 1
        while start > 0 and time[start - 1] == time[stop - 1]:
            start -= 1
        group_weight = weight[start:stop]
        group_x = covariates[start:stop]
        s0 += float(group_weight.sum())
        s1 += group_weight @ group_x
        s2 += (group_x * group_weight[:, None]).T @ group_x

        is_event = event[start:stop] == 1
        n_events = int(is_event.sum())
        if n_events:
            mean_x = s1 / s0
            loglik += float(eta[start:stop][is_event].sum()) - n_events * (math.log(s0) + shift)
            gradient += group_x[is_event].sum(axis=0) - n_events * mean_x
            information += n_events * (s2 / s0 - np.outer(mean_x, mean_x))
        stop = start

    return loglik, gradient, information


def _checked_information(information: Any, *, at_start: bool) -> None:
    """Raise unless the information matrix is numerically positive definite (a unique maximum exists).

    At the starting point (all coefficients zero) a singular matrix means the design itself is degenerate; later
    in the iteration it means the information is vanishing because a coefficient is running off to infinity.
    """
    eigenvalues = np.linalg.eigvalsh(information)
    tolerance = information.shape[0] * np.finfo(float).eps * max(float(eigenvalues[-1]), 1e-300)
    if float(eigenvalues[0]) > tolerance:
        return
    if at_start:
        raise ValueError(
            "The Cox information matrix is singular, so the coefficients are not identified: a covariate is "
            "constant or the covariates are collinear (or the sample has too few events for this many "
            "covariates). Drop or combine covariates, or use multi_omic_survival_model for penalised fits."
        )
    raise ValueError(
        "Cox regression did not converge: the information matrix vanished as the coefficients grew, so a "
        "coefficient is infinite because a covariate (or a combination of them) separates the events perfectly"
    )


def cox_regression(
    time: list[float],
    event: list[int],
    covariates: list[list[float]] | Any,
    covariate_names: list[str] | None = None,
) -> dict[str, Any]:
    """Cox proportional hazards regression.

    Fits a Cox PH model by maximising the partial likelihood using
    Newton-Raphson optimisation with exact gradient and Hessian (step halving
    whenever a step lowers the likelihood).  Handles tied event times with the
    Breslow approximation: every subject observed for at least the event time,
    tied subjects included, is in the risk set.

    The hazard for subject i is:
        h_i(t) = h_0(t) * exp(X_i @ beta)

    where h_0 is the baseline hazard and beta are the regression
    coefficients.

    Args:
        time: Observed survival times (n subjects), finite.
        event: Event indicators (1 = event/death, 0 = censored).
        covariates: Covariate matrix (n x p) as numpy array or
            list-of-lists.
        covariate_names: Optional names for the p covariates.

    Returns:
        Dictionary with keys:
            coefficients: Regression coefficients (list of float).
            hazard_ratios: exp(coefficients) per covariate.
            se: Standard errors of coefficients (square roots of the diagonal of the inverse information
                matrix at the maximum).
            p_values: Wald test p-values for each coefficient.
            concordance_index: Harrell's C-index for the fitted model.
            log_likelihood: Log partial likelihood (Breslow) at the fitted coefficients.
            covariate_names: Names of covariates.

    Raises:
        ValueError: If inputs have incompatible lengths, a time is not finite, an event indicator is not 0
            or 1, there is no event, the covariates are collinear or constant (the information matrix is
            singular), or the iteration does not converge (typically complete separation, where a
            coefficient is infinite).
        ImportError: If numpy is not available.
    """
    n = len(time)
    if n != len(event):
        raise ValueError(f"time ({n}) and event ({len(event)}) must have the same length")
    if n == 0:
        raise ValueError("Cannot fit Cox model on empty data")
    if any(e not in (0, 1) for e in event):
        raise ValueError("event must contain only 0 (censored) and 1 (event)")
    if sum(event) == 0:
        raise ValueError("No events in the data; Cox model cannot be fitted")

    if not HAS_NUMPY:
        raise ImportError("numpy is required for Cox regression. " "Install with: uv pip install numpy")

    logger.info("Fitting Cox PH model: n=%d, events=%d", n, sum(event))

    X = np.asarray(covariates, dtype=np.float64)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    if X.shape[0] != n:
        raise ValueError(f"covariates has {X.shape[0]} rows but time has {n} entries")

    p = X.shape[1]
    t = np.asarray(time, dtype=np.float64)
    d = np.asarray(event, dtype=np.float64)
    if not np.all(np.isfinite(t)):
        raise ValueError("time must be finite")
    if not np.all(np.isfinite(X)):
        raise ValueError("covariates must be finite")

    if covariate_names is None:
        covariate_names = [f"X{i}" for i in range(p)]
    if len(covariate_names) != p:
        raise ValueError(f"covariate_names has {len(covariate_names)} entries but need {p}")

    order = np.argsort(t, kind="stable")
    t_sorted, d_sorted, x_sorted = t[order], d[order], X[order]

    beta = np.zeros(p, dtype=np.float64)
    log_lik, gradient, information = _breslow_partial_likelihood(t_sorted, d_sorted, x_sorted, beta)
    converged = False
    for iteration in range(1, _COX_MAX_ITERATIONS + 1):
        _checked_information(information, at_start=iteration == 1)
        step = np.linalg.solve(information, gradient)

        # Newton step, halved while it lowers the likelihood
        for _ in range(_COX_MAX_HALVINGS):
            candidate = beta + step
            new_log_lik, new_gradient, new_information = _breslow_partial_likelihood(
                t_sorted, d_sorted, x_sorted, candidate
            )
            if new_log_lik >= log_lik - 1e-12 * abs(log_lik):
                break
            step = step / 2.0
        else:
            raise ValueError("Cox regression could not find an increasing step; the likelihood is not improvable")

        converged = (
            abs(new_log_lik - log_lik) <= 1e-10 * (abs(new_log_lik) + 0.1) and float(np.max(np.abs(step))) < 1e-7
        )
        beta, log_lik, gradient, information = candidate, new_log_lik, new_gradient, new_information
        if converged:
            logger.info("Cox regression converged at iteration %d", iteration)
            break
    if not converged:
        raise ValueError(
            f"Cox regression did not converge in {_COX_MAX_ITERATIONS} iterations; a coefficient is probably "
            "infinite because a covariate (or a combination of them) separates the events perfectly"
        )

    # Standard errors from the inverse information matrix at the maximum
    _checked_information(information, at_start=False)
    se = np.sqrt(np.diag(np.linalg.inv(information)))

    # Wald test p-values: two-sided normal tail, exact through erfc
    p_values = [math.erfc(abs(float(b / s)) / math.sqrt(2.0)) for b, s in zip(beta, se)]

    c_result = compute_concordance_index((X @ beta).tolist(), t.tolist(), d.astype(int).tolist())

    logger.info(
        "Cox regression complete: C-index=%.4f, log_lik=%.4f",
        c_result["c_index"],
        log_lik,
    )

    return {
        "coefficients": beta.tolist(),
        "hazard_ratios": np.exp(beta).tolist(),
        "se": se.tolist(),
        "p_values": p_values,
        "concordance_index": c_result["c_index"],
        "log_likelihood": log_lik,
        "covariate_names": covariate_names,
    }


# ---------------------------------------------------------------------------
# Kaplan-Meier Estimator
# ---------------------------------------------------------------------------


def kaplan_meier(
    time: list[float],
    event: list[int],
    groups: list[int] | None = None,
) -> dict[str, Any]:
    """Kaplan-Meier survival curve estimation.

    Computes the product-limit estimator for the survival function.  If
    *groups* is provided, computes separate curves for each group.

    Confidence intervals use the Greenwood formula:
        Var(S(t)) = S(t)^2 * sum(d_i / (n_i * (n_i - d_i)))

    Args:
        time: Observed survival times.
        event: Event indicators (1 = event, 0 = censored).
        groups: Optional group labels.  If provided, separate curves are
            computed per group.

    Returns:
        If groups is None, dictionary with keys:
            times: Sorted unique event times.
            survival_prob: Survival probability at each time.
            confidence_lower: Lower 95% confidence bound.
            confidence_upper: Upper 95% confidence bound.
            n_at_risk: Number at risk at each time.
            median_survival: Median survival time (or None if not reached).

        If groups is provided, dictionary with keys:
            groups: dict mapping group label to the per-group result dict.

    Raises:
        ValueError: If inputs have incompatible lengths.
    """
    n = len(time)
    if n != len(event):
        raise ValueError(f"time ({n}) and event ({len(event)}) must have same length")
    if groups is not None and len(groups) != n:
        raise ValueError(f"groups ({len(groups)}) must have same length as time ({n})")
    if n == 0:
        raise ValueError("Cannot compute KM on empty data")

    logger.info("Computing Kaplan-Meier: n=%d, events=%d", n, sum(event))

    if groups is not None:
        unique_groups = sorted(set(groups))
        result_groups: dict[int, dict[str, Any]] = {}
        for g in unique_groups:
            idx = [i for i in range(n) if groups[i] == g]
            sub_time = [time[i] for i in idx]
            sub_event = [event[i] for i in idx]
            result_groups[g] = _km_single(sub_time, sub_event)
        return {"groups": result_groups}
    else:
        return _km_single(time, event)


def _km_single(time: list[float], event: list[int]) -> dict[str, Any]:
    """Compute KM curve for a single group.

    Args:
        time: Observed times.
        event: Event indicators.

    Returns:
        Dictionary with KM curve data.
    """
    # Pair and sort by time
    paired = sorted(zip(time, event), key=lambda x: x[0])
    n = len(paired)

    times_out: list[float] = []
    surv_prob: list[float] = []
    conf_lower: list[float] = []
    conf_upper: list[float] = []
    n_at_risk_out: list[int] = []

    s = 1.0
    var_sum = 0.0  # Greenwood variance accumulator
    at_risk = n
    median_survival: float | None = None

    i = 0
    while i < n:
        t_i = paired[i][0]
        # Count events and censored at this time
        d_i = 0  # events
        c_i = 0  # censored
        while i < n and paired[i][0] == t_i:
            if paired[i][1] == 1:
                d_i += 1
            else:
                c_i += 1
            i += 1

        if d_i > 0:
            s *= (at_risk - d_i) / at_risk if at_risk > 0 else 0.0

            # Greenwood variance
            if at_risk > 0 and (at_risk - d_i) > 0:
                var_sum += d_i / (at_risk * (at_risk - d_i))

            var_s = s**2 * var_sum
            se = math.sqrt(max(var_s, 0.0))

            times_out.append(t_i)
            surv_prob.append(s)
            conf_lower.append(max(0.0, s - 1.96 * se))
            conf_upper.append(min(1.0, s + 1.96 * se))
            n_at_risk_out.append(at_risk)

            # Check median
            if median_survival is None and s <= 0.5:
                median_survival = t_i

        at_risk -= d_i + c_i

    return {
        "times": times_out,
        "survival_prob": surv_prob,
        "confidence_lower": conf_lower,
        "confidence_upper": conf_upper,
        "n_at_risk": n_at_risk_out,
        "median_survival": median_survival,
    }


# ---------------------------------------------------------------------------
# Log-Rank Test
# ---------------------------------------------------------------------------


def _solve_symmetric(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """Solve ``matrix @ x = rhs`` by Gaussian elimination with partial pivoting (pure Python).

    Raises:
        ValueError: If the matrix is numerically singular.
    """
    size = len(rhs)
    aug = [list(matrix[i]) + [rhs[i]] for i in range(size)]
    scale = max((abs(matrix[i][i]) for i in range(size)), default=0.0)
    for col in range(size):
        pivot = max(range(col, size), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot][col]) <= 1e-12 * scale:
            raise ValueError(
                "The log-rank covariance matrix is singular, so the groups cannot be compared: some groups "
                "are never separable at the observed event times (check for groups with identical at-risk sets)"
            )
        aug[col], aug[pivot] = aug[pivot], aug[col]
        for row in range(col + 1, size):
            factor = aug[row][col] / aug[col][col]
            for k in range(col, size + 1):
                aug[row][k] -= factor * aug[col][k]
    solution = [0.0] * size
    for row in range(size - 1, -1, -1):
        solution[row] = (aug[row][size] - sum(aug[row][k] * solution[k] for k in range(row + 1, size))) / aug[row][row]
    return solution


def log_rank_test(
    time: list[float],
    event: list[int],
    groups: list[int],
) -> dict[str, Any]:
    """Log-rank test comparing survival between groups.

    Implements the Mantel-Haenszel test statistic comparing observed vs
    expected events in each group at every distinct event time.  With ``d``
    events among ``N`` subjects at risk, ``n_g`` of them in group ``g``, the
    expected events of group ``g`` are ``d n_g / N`` and the hypergeometric
    covariance of the group scores is
    ``d (N - d) n_g (N delta_gh - n_h) / (N**2 (N - 1))``.  The statistic is
    ``(O - E)' V^-1 (O - E)`` over all groups but one, which is chi-squared
    with (number of groups - 1) degrees of freedom; for two groups it is
    ``(O - E)**2 / V``.  A group that is never at risk when an event occurs
    carries no information and is left out (and the degrees of freedom with it).

    Args:
        time: Observed survival times.
        event: Event indicators (1 = event, 0 = censored).
        groups: Group labels for each subject.

    Returns:
        Dictionary with keys:
            chi2: Chi-squared test statistic.
            p_value: P-value from chi-squared distribution with ``df`` degrees of freedom.
            df: Degrees of freedom (the number of compared groups minus one).
            observed_expected_per_group: {group: {observed, expected}}.

    Raises:
        ValueError: If inputs have incompatible lengths, fewer than 2 groups are given, no event occurs
            while two groups are at risk (the groups cannot be compared, and a p-value of 1.0 would
            claim they were), or the covariance matrix is singular.
    """
    n = len(time)
    if n != len(event) or n != len(groups):
        raise ValueError("time, event, and groups must have the same length")
    if n == 0:
        raise ValueError("Cannot perform log-rank test on empty data")

    unique_groups = sorted(set(groups))
    if len(unique_groups) < 2:
        raise ValueError("Need at least 2 groups for log-rank test")
    n_groups = len(unique_groups)
    slot = {g: i for i, g in enumerate(unique_groups)}

    logger.info("Running log-rank test: n=%d, groups=%d", n, n_groups)

    # Sort by time
    order = sorted(range(n), key=lambda i: time[i])
    t_sorted = [time[i] for i in order]
    e_sorted = [event[i] for i in order]
    g_sorted = [slot[groups[i]] for i in order]

    # Unique event times
    event_times = sorted(set(t_sorted[i] for i in range(n) if e_sorted[i] == 1))

    # Track at-risk counts per group
    at_risk = [0] * n_groups
    for g in g_sorted:
        at_risk[g] += 1
    observed = [0.0] * n_groups
    expected = [0.0] * n_groups
    covariance = [[0.0] * n_groups for _ in range(n_groups)]

    ptr = 0  # pointer into sorted arrays
    for t_event in event_times:
        # Remove subjects censored/evented before this time
        while ptr < n and t_sorted[ptr] < t_event:
            at_risk[g_sorted[ptr]] -= 1
            ptr += 1

        # Count events at this time
        total_at_risk = sum(at_risk)
        if total_at_risk == 0:
            continue

        events_at_t = [0] * n_groups
        temp_ptr = ptr
        while temp_ptr < n and t_sorted[temp_ptr] == t_event:
            if e_sorted[temp_ptr] == 1:
                events_at_t[g_sorted[temp_ptr]] += 1
            temp_ptr += 1

        d_j = sum(events_at_t)

        for g in range(n_groups):
            observed[g] += events_at_t[g]
            expected[g] += d_j * at_risk[g] / total_at_risk

        # Hypergeometric covariance of the group scores at this event time
        if total_at_risk > 1:
            factor = d_j * (total_at_risk - d_j) / (total_at_risk**2 * (total_at_risk - 1))
            for g in range(n_groups):
                for h in range(n_groups):
                    own = total_at_risk if g == h else 0
                    covariance[g][h] += factor * at_risk[g] * (own - at_risk[h])

        # Advance pointer past events at this time
        while ptr < n and t_sorted[ptr] == t_event:
            at_risk[g_sorted[ptr]] -= 1
            ptr += 1

    # Chi-squared statistic over the groups that were at risk at an event time, all but the last of them
    active = [g for g in range(n_groups) if expected[g] > 0]
    df = len(active) - 1
    if df < 1:
        raise ValueError(
            "No event occurs while two or more groups are at risk, so the groups cannot be compared "
            "(a chi-squared of 0 with p = 1.0 would claim there is no difference)"
        )
    kept = active[:-1]
    deviation = [observed[g] - expected[g] for g in kept]
    weighted = _solve_symmetric([[covariance[g][h] for h in kept] for g in kept], deviation)
    chi2 = max(0.0, sum(a * b for a, b in zip(deviation, weighted)))

    # P-value
    if sp_stats is not None:
        p_value = float(sp_stats.chi2.sf(chi2, df))
    else:
        logger.warning("scipy is not available: the log-rank p-value is a Wilson-Hilferty approximation")
        z = ((chi2 / df) ** (1.0 / 3) - (1 - 2.0 / (9 * df))) / math.sqrt(2.0 / (9 * df))
        p_value = 0.5 * math.erfc(z / math.sqrt(2.0))

    obs_exp: dict[int, dict[str, float]] = {}
    for g in unique_groups:
        obs_exp[g] = {"observed": observed[slot[g]], "expected": expected[slot[g]]}

    logger.info("Log-rank test: chi2=%.4f, df=%d, p=%.4e", chi2, df, p_value)

    return {
        "chi2": chi2,
        "p_value": p_value,
        "df": df,
        "observed_expected_per_group": obs_exp,
    }


# ---------------------------------------------------------------------------
# Multi-Omic Survival Model (Lasso-Cox)
# ---------------------------------------------------------------------------


def multi_omic_survival_model(
    omic_features: dict[str, Any],
    time: list[float],
    event: list[int],
    method: str = "lasso_cox",
) -> dict[str, Any]:
    """Build survival model from multi-omic features with regularisation.

    Concatenates features from multiple omic layers and fits a Cox model
    with L1 (lasso) penalty for automatic feature selection.

    The L1-penalised partial log-likelihood is optimised using coordinate
    descent.  The penalty strength is not tuned: it is fixed at one tenth of the
    smallest penalty that would remove every feature (lambda_max / 10), with the
    features standardized, so the number of selected features is a property of
    that choice and of the data, not the outcome of a cross-validated search.

    Args:
        omic_features: Mapping from omic name to feature matrix (n x p_m)
            as numpy array or list-of-lists.
        time: Survival times for n subjects.
        event: Event indicators for n subjects.
        method: Regularisation method.  Currently only ``"lasso_cox"`` is
            supported.

    Returns:
        Dictionary with keys:
            selected_features: Names of features with non-zero coefficients.
            coefficients: Coefficient values for selected features.
            c_index: Concordance index of the model.
            risk_scores: Predicted risk scores for each subject.

    Raises:
        ValueError: If inputs are incompatible or method is unsupported.
        ImportError: If numpy is not available.
    """
    if method != "lasso_cox":
        raise ValueError(f"Unsupported method: '{method}'. Use 'lasso_cox'.")
    if not omic_features:
        raise ValueError("omic_features must contain at least one entry")

    if not HAS_NUMPY:
        raise ImportError("numpy is required for multi-omic survival model. " "Install with: uv pip install numpy")

    n = len(time)
    if n != len(event):
        raise ValueError(f"time ({n}) and event ({len(event)}) must have same length")

    logger.info("Fitting multi-omic survival model: n=%d, method=%s", n, method)

    # Concatenate features with names
    feature_names: list[str] = []
    parts: list[Any] = []
    for omic_name, mat in omic_features.items():
        arr = np.asarray(mat, dtype=np.float64)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        if arr.shape[0] != n:
            raise ValueError(f"omic '{omic_name}' has {arr.shape[0]} samples but time has {n}")
        for j in range(arr.shape[1]):
            feature_names.append(f"{omic_name}_f{j}")
        parts.append(arr)

    X = np.hstack(parts)
    p = X.shape[1]

    # Standardise
    mean = X.mean(axis=0)
    std = X.std(axis=0) + 1e-10
    X_std = (X - mean) / std

    t_arr = np.asarray(time, dtype=np.float64)
    d_arr = np.asarray(event, dtype=np.float64)

    # Sort ascending
    order = np.argsort(t_arr)
    t_arr = t_arr[order]
    d_arr = d_arr[order]
    X_std = X_std[order]

    # Coordinate descent for Lasso-Cox
    beta = np.zeros(p, dtype=np.float64)
    # Choose lambda via heuristic: lambda_max / 10
    # lambda_max is the max absolute gradient at beta=0
    exp_eta = np.ones(n, dtype=np.float64)
    S0 = np.cumsum(exp_eta[::-1])[::-1]
    grad0 = np.zeros(p, dtype=np.float64)
    S1 = np.cumsum((X_std * exp_eta[:, None])[::-1], axis=0)[::-1]
    for i in range(n):
        if d_arr[i] == 0:
            continue
        s0_i = S0[i] + 1e-10
        grad0 += X_std[i] - S1[i] / s0_i

    lam_max = float(np.max(np.abs(grad0))) / n
    lam = lam_max * 0.1

    max_iter_cd = 100
    tol_cd = 1e-5

    for iteration in range(max_iter_cd):
        eta = X_std @ beta
        exp_eta = np.exp(eta - eta.max())

        S0 = np.cumsum(exp_eta[::-1])[::-1]
        S1 = np.cumsum((X_std * exp_eta[:, None])[::-1], axis=0)[::-1]

        beta_old = beta.copy()

        for j in range(p):
            # Partial residual gradient for coordinate j
            grad_j = 0.0
            hess_j = 0.0
            for i in range(n):
                if d_arr[i] == 0:
                    continue
                s0_i = S0[i] + 1e-10
                mean_xj = S1[i, j] / s0_i
                grad_j += X_std[i, j] - mean_xj
                # Approximate Hessian
                hess_j += (S1[i, j] ** 2) / (s0_i**2)

            grad_j /= n
            hess_j = max(hess_j / n, 1e-6)

            # Soft-thresholding
            z = beta[j] + grad_j / hess_j
            beta[j] = math.copysign(max(abs(z) - lam / hess_j, 0.0), z)

        if np.max(np.abs(beta - beta_old)) < tol_cd:
            logger.info("Lasso-Cox converged at iteration %d", iteration + 1)
            break

    # Selected features (non-zero coefficients)
    selected_mask = np.abs(beta) > 1e-8
    selected_indices = np.where(selected_mask)[0]
    selected_names = [feature_names[i] for i in selected_indices]
    selected_coeffs = beta[selected_indices].tolist()

    # Risk scores (in original order)
    risk_scores_ordered = X_std @ beta
    # Revert to original order
    inv_order = np.argsort(order)
    risk_scores = risk_scores_ordered[inv_order].tolist()

    # C-index
    c_result = compute_concordance_index(
        risk_scores,
        [float(time[i]) for i in range(n)],
        [int(event[i]) for i in range(n)],
    )

    logger.info(
        "Multi-omic survival model: %d/%d features selected, C-index=%.4f",
        len(selected_names),
        p,
        c_result["c_index"],
    )

    return {
        "selected_features": selected_names,
        "coefficients": selected_coeffs,
        "c_index": c_result["c_index"],
        "risk_scores": risk_scores,
    }


# ---------------------------------------------------------------------------
# Risk Stratification
# ---------------------------------------------------------------------------


def risk_stratification(
    risk_scores: list[float],
    time: list[float],
    event: list[int],
    n_groups: int = 2,
) -> dict[str, Any]:
    """Stratify patients into risk groups.

    Uses the median risk score (for n_groups=2) or quantiles to split
    patients into risk groups, then compares survival between groups.

    Args:
        risk_scores: Predicted risk scores per patient.
        time: Survival times.
        event: Event indicators.
        n_groups: Number of risk groups (default 2: high/low).

    Returns:
        Dictionary with keys:
            group_labels: Assigned risk group per patient (0 = lowest risk).
            survival_curves: Per-group Kaplan-Meier curves.
            log_rank_p: Log-rank test p-value comparing groups.
            hazard_ratio: Hazard ratio of highest vs lowest risk group
                (approximate, based on log-rank O/E ratio).

    Raises:
        ValueError: If inputs are incompatible or n_groups < 2, or the risk groups cannot be compared
            (from :func:`log_rank_test`: all scores fall in one group, or no event occurs while two groups
            are at risk).
    """
    n = len(risk_scores)
    if n != len(time) or n != len(event):
        raise ValueError("risk_scores, time, and event must have the same length")
    if n_groups < 2:
        raise ValueError(f"n_groups must be >= 2, got {n_groups}")
    if n < n_groups:
        raise ValueError(f"Need at least {n_groups} subjects, got {n}")

    logger.info("Stratifying patients into %d risk groups", n_groups)

    # Compute quantile thresholds
    sorted_scores = sorted(risk_scores)
    thresholds: list[float] = []
    for q in range(1, n_groups):
        idx = int(len(sorted_scores) * q / n_groups)
        thresholds.append(sorted_scores[min(idx, len(sorted_scores) - 1)])

    # Assign groups
    group_labels: list[int] = []
    for score in risk_scores:
        group = 0
        for t_val in thresholds:
            if score > t_val:
                group += 1
        group_labels.append(min(group, n_groups - 1))

    # Kaplan-Meier per group
    km_result = kaplan_meier(time, event, groups=group_labels)

    # Log-rank test
    lr_result = log_rank_test(time, event, group_labels)

    # Approximate hazard ratio (highest vs lowest group)
    obs_exp = lr_result["observed_expected_per_group"]
    groups_sorted = sorted(obs_exp.keys())
    low_group = groups_sorted[0]
    high_group = groups_sorted[-1]

    o_high = obs_exp[high_group]["observed"]
    e_high = obs_exp[high_group]["expected"]
    o_low = obs_exp[low_group]["observed"]
    e_low = obs_exp[low_group]["expected"]

    if e_high > 0 and e_low > 0:
        hazard_ratio = (o_high / e_high) / (o_low / e_low)
    elif o_high > 0:
        hazard_ratio = float("inf")
    else:
        hazard_ratio = 0.0

    logger.info(
        "Risk stratification: HR=%.4f, log-rank p=%.4e",
        hazard_ratio,
        lr_result["p_value"],
    )

    return {
        "group_labels": group_labels,
        "survival_curves": km_result.get("groups", km_result),
        "log_rank_p": lr_result["p_value"],
        "hazard_ratio": hazard_ratio,
    }


# ---------------------------------------------------------------------------
# Concordance Index (Harrell's C)
# ---------------------------------------------------------------------------


def compute_concordance_index(
    risk_scores: list[float],
    time: list[float],
    event: list[int],
) -> dict[str, Any]:
    """Harrell's concordance index for survival model discrimination.

    The C-index measures the probability that, for a randomly chosen pair
    of subjects where one experienced the event earlier, the subject with
    the higher risk score is the one who experienced the event first.

    A C-index of 0.5 indicates random prediction; 1.0 indicates perfect
    discrimination.

    Args:
        risk_scores: Predicted risk scores per subject.
        time: Observed survival times.
        event: Event indicators (1 = event, 0 = censored).

    Returns:
        Dictionary with keys:
            c_index: Concordance index value.  When no pair of subjects is comparable (for
                example no event at all) it is reported as 0.5, the value of a model that does
                not discriminate, and ``se`` as 0.0: a convention, not an estimate.
            se: Standard error estimate.
            n_concordant: Number of concordant pairs.
            n_discordant: Number of discordant pairs.
            n_tied: Number of tied pairs.

    Raises:
        ValueError: If inputs have incompatible lengths.
    """
    n = len(risk_scores)
    if n != len(time) or n != len(event):
        raise ValueError("risk_scores, time, and event must have the same length")

    concordant = 0
    discordant = 0
    tied = 0
    total_pairs = 0

    for i in range(n):
        if event[i] == 0:
            continue
        for j in range(n):
            if i == j:
                continue
            if time[j] <= time[i] and event[j] == 1:
                # j had event no later than i -- skip (i is the one with event)
                # Only count pairs where i's event time < j's time or j is censored after
                continue
            if time[i] < time[j] or (time[i] <= time[j] and event[j] == 0):
                # Valid pair: i had event, j survived longer or censored after
                total_pairs += 1
                if risk_scores[i] > risk_scores[j]:
                    concordant += 1
                elif risk_scores[i] < risk_scores[j]:
                    discordant += 1
                else:
                    tied += 1

    if total_pairs == 0:
        c_index = 0.5
    else:
        c_index = (concordant + 0.5 * tied) / total_pairs

    # SE estimate (Noether's formula)
    if total_pairs > 1:
        q1 = c_index / (2.0 - c_index)
        q2 = 2 * c_index**2 / (1.0 + c_index)
        n_events = sum(event)
        n_censored = n - n_events
        if n_events > 0 and n_censored > 0:
            se = math.sqrt(
                (c_index * (1 - c_index) + (n_events - 1) * (q1 - c_index**2) + (n_censored - 1) * (q2 - c_index**2))
                / (n_events * n_censored)
            )
        else:
            se = 0.0
    else:
        se = 0.0

    return {
        "c_index": c_index,
        "se": se,
        "n_concordant": concordant,
        "n_discordant": discordant,
        "n_tied": tied,
    }
