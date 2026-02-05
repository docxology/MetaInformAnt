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
import random
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional scientific dependencies
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Cox Proportional Hazards Regression
# ---------------------------------------------------------------------------


def cox_regression(
    time: list[float],
    event: list[int],
    covariates: list[list[float]] | Any,
    covariate_names: list[str] | None = None,
) -> dict[str, Any]:
    """Cox proportional hazards regression.

    Fits a Cox PH model by maximising the partial likelihood using
    Newton-Raphson optimisation.  Handles tied event times with the
    Breslow approximation.

    The hazard for subject i is:
        h_i(t) = h_0(t) * exp(X_i @ beta)

    where h_0 is the baseline hazard and beta are the regression
    coefficients.

    Args:
        time: Observed survival times (n subjects).
        event: Event indicators (1 = event/death, 0 = censored).
        covariates: Covariate matrix (n x p) as numpy array or
            list-of-lists.
        covariate_names: Optional names for the p covariates.

    Returns:
        Dictionary with keys:
            coefficients: Regression coefficients (list of float).
            hazard_ratios: exp(coefficients) per covariate.
            se: Standard errors of coefficients.
            p_values: Wald test p-values for each coefficient.
            concordance_index: Harrell's C-index for the fitted model.
            log_likelihood: Final log partial likelihood.
            covariate_names: Names of covariates.

    Raises:
        ValueError: If inputs have incompatible lengths or no events.
        ImportError: If numpy is not available.
    """
    n = len(time)
    if n != len(event):
        raise ValueError(f"time ({n}) and event ({len(event)}) must have the same length")
    if n == 0:
        raise ValueError("Cannot fit Cox model on empty data")
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

    if covariate_names is None:
        covariate_names = [f"X{i}" for i in range(p)]
    if len(covariate_names) != p:
        raise ValueError(f"covariate_names has {len(covariate_names)} entries but need {p}")

    # Sort by time (descending for risk set construction)
    order = np.argsort(-t)
    t = t[order]
    d = d[order]
    X = X[order]

    # Newton-Raphson
    beta = np.zeros(p, dtype=np.float64)
    max_iter = 50
    tol = 1e-6

    for iteration in range(max_iter):
        # Compute risk scores
        eta = X @ beta
        exp_eta = np.exp(eta - eta.max())  # numerical stability

        # Cumulative sums for risk set (already sorted descending by time,
        # so forward cumsum gives the risk set at each time)
        # We need cumsum from the end (latest time first after descending sort
        # means first entry is latest time; risk set at time t includes all
        # subjects with observed time >= t)
        # After descending sort: t[0] >= t[1] >= ... >= t[n-1]
        # Risk set at t[i] includes subjects 0..i (those with time >= t[i])
        # But with descending sort that reverses intuition -- let's use
        # ascending sort and reverse cumsum instead.

        # Re-sort ascending
        asc = np.argsort(t)
        t_asc = t[asc]
        d_asc = d[asc]
        X_asc = X[asc]
        exp_eta_asc = exp_eta[asc]

        # Risk set at time t[i]: indices j where t_asc[j] >= t_asc[i]
        # Reverse cumulative sums
        S0 = np.cumsum(exp_eta_asc[::-1])[::-1]  # sum exp_eta for risk set
        S1 = np.cumsum((X_asc * exp_eta_asc[:, None])[::-1], axis=0)[::-1]
        S2_diag = np.cumsum((X_asc**2 * exp_eta_asc[:, None])[::-1], axis=0)[::-1]

        # Gradient and Hessian
        gradient = np.zeros(p, dtype=np.float64)
        hessian = np.zeros((p, p), dtype=np.float64)
        log_lik = 0.0

        for i in range(n):
            if d_asc[i] == 0:
                continue
            s0_i = S0[i] + 1e-10
            s1_i = S1[i]
            mean_x = s1_i / s0_i
            gradient += X_asc[i] - mean_x
            log_lik += float(np.dot(X_asc[i], beta)) - math.log(s0_i)

            # Hessian contribution
            for j1 in range(p):
                for j2 in range(p):
                    hessian[j1, j2] -= (
                        (S2_diag[i, j1] / s0_i * (j1 == j2) - s1_i[j1] * s1_i[j2] / (s0_i**2))
                        if j1 == j2
                        else (-s1_i[j1] * s1_i[j2] / (s0_i**2))
                    )

        # Simplified Hessian (diagonal + off-diagonal from S1 outer product)
        hessian_full = np.zeros((p, p), dtype=np.float64)
        for i in range(n):
            if d_asc[i] == 0:
                continue
            s0_i = S0[i] + 1e-10
            s1_i = S1[i]
            mean_x = s1_i / s0_i
            var_x = S2_diag[i] / s0_i - mean_x**2
            for j in range(p):
                hessian_full[j, j] -= var_x[j]
                for j2 in range(j + 1, p):
                    cov_jj2 = s1_i[j] * s1_i[j2] / (s0_i**2)
                    # Cross-covariance approximation
                    hessian_full[j, j2] += cov_jj2
                    hessian_full[j2, j] += cov_jj2

        # Newton step
        try:
            delta = np.linalg.solve(hessian_full - 1e-6 * np.eye(p), gradient)
        except np.linalg.LinAlgError:
            logger.warning("Hessian singular at iteration %d; using gradient step", iteration)
            delta = 0.01 * gradient

        beta_new = beta - delta
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            logger.info("Cox regression converged at iteration %d", iteration + 1)
            break
        beta = beta_new

    # Standard errors from inverse Hessian
    try:
        info = -hessian_full
        cov_matrix = np.linalg.inv(info + 1e-8 * np.eye(p))
        se = np.sqrt(np.abs(np.diag(cov_matrix)))
    except np.linalg.LinAlgError:
        se = np.full(p, float("nan"))

    # Wald test p-values
    z_scores = beta / (se + 1e-10)
    p_values: list[float] = []
    for z in z_scores:
        pval = 2.0 * (0.5 * math.erfc(abs(float(z)) / math.sqrt(2.0)))
        p_values.append(pval)

    # Concordance index
    risk_scores_final = (X @ beta).tolist()
    # Need original order for concordance
    original_risk = [0.0] * n
    original_time = [0.0] * n
    original_event = [0] * n
    for i, idx in enumerate(asc):
        original_risk[idx] = risk_scores_final[i]
        original_time[idx] = float(t_asc[i])
        original_event[idx] = int(d_asc[i])

    c_result = compute_concordance_index(original_risk, original_time, original_event)

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
            s_prev = s
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


def log_rank_test(
    time: list[float],
    event: list[int],
    groups: list[int],
) -> dict[str, Any]:
    """Log-rank test comparing survival between groups.

    Implements the Mantel-Haenszel test statistic comparing observed vs
    expected events in each group at every distinct event time.

    Args:
        time: Observed survival times.
        event: Event indicators (1 = event, 0 = censored).
        groups: Group labels for each subject.

    Returns:
        Dictionary with keys:
            chi2: Chi-squared test statistic.
            p_value: P-value from chi-squared distribution with (g-1) df.
            df: Degrees of freedom.
            observed_expected_per_group: {group: {observed, expected}}.

    Raises:
        ValueError: If inputs have incompatible lengths or fewer than 2 groups.
    """
    n = len(time)
    if n != len(event) or n != len(groups):
        raise ValueError("time, event, and groups must have the same length")
    if n == 0:
        raise ValueError("Cannot perform log-rank test on empty data")

    unique_groups = sorted(set(groups))
    if len(unique_groups) < 2:
        raise ValueError("Need at least 2 groups for log-rank test")

    logger.info("Running log-rank test: n=%d, groups=%d", n, len(unique_groups))

    # Sort by time
    order = sorted(range(n), key=lambda i: time[i])
    t_sorted = [time[i] for i in order]
    e_sorted = [event[i] for i in order]
    g_sorted = [groups[i] for i in order]

    # Unique event times
    event_times = sorted(set(t_sorted[i] for i in range(n) if e_sorted[i] == 1))

    # Track at-risk counts per group
    at_risk: dict[int, int] = {g: sum(1 for gi in groups if gi == g) for g in unique_groups}
    observed: dict[int, float] = {g: 0.0 for g in unique_groups}
    expected: dict[int, float] = {g: 0.0 for g in unique_groups}
    variance: dict[int, float] = {g: 0.0 for g in unique_groups}

    ptr = 0  # pointer into sorted arrays
    for t_event in event_times:
        # Remove subjects censored/evented before this time
        while ptr < n and t_sorted[ptr] < t_event:
            at_risk[g_sorted[ptr]] -= 1
            ptr += 1

        # Count events at this time
        total_at_risk = sum(at_risk.values())
        if total_at_risk == 0:
            continue

        events_at_t: dict[int, int] = {g: 0 for g in unique_groups}
        temp_ptr = ptr
        while temp_ptr < n and t_sorted[temp_ptr] == t_event:
            if e_sorted[temp_ptr] == 1:
                events_at_t[g_sorted[temp_ptr]] += 1
            temp_ptr += 1

        d_j = sum(events_at_t.values())

        for g in unique_groups:
            n_g = at_risk[g]
            e_g = d_j * n_g / total_at_risk if total_at_risk > 0 else 0.0
            observed[g] += events_at_t[g]
            expected[g] += e_g

            # Variance contribution
            if total_at_risk > 1:
                v_g = (d_j * n_g * (total_at_risk - n_g) * (total_at_risk - d_j)) / (
                    total_at_risk**2 * (total_at_risk - 1)
                )
                variance[g] += v_g

        # Advance pointer past events at this time
        while ptr < n and t_sorted[ptr] == t_event:
            at_risk[g_sorted[ptr]] -= 1
            ptr += 1

    # Chi-squared statistic (use first g-1 groups)
    df = len(unique_groups) - 1
    chi2 = 0.0
    for g in unique_groups[:-1]:
        v = variance[g]
        if v > 0:
            chi2 += (observed[g] - expected[g]) ** 2 / v

    # P-value
    if HAS_NUMPY:
        from scipy import stats as _stats

        p_value = float(_stats.chi2.sf(chi2, df))
    else:
        # Approximate
        if df > 0:
            z = ((chi2 / df) ** (1.0 / 3) - (1 - 2.0 / (9 * df))) / math.sqrt(2.0 / (9 * df))
            p_value = 0.5 * math.erfc(z / math.sqrt(2.0))
        else:
            p_value = 1.0

    obs_exp: dict[int, dict[str, float]] = {}
    for g in unique_groups:
        obs_exp[g] = {"observed": observed[g], "expected": expected[g]}

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
    descent.

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
    for i in range(n):
        if d_arr[i] == 0:
            continue
        s0_i = S0[i] + 1e-10
        cum_weighted = np.cumsum((X_std * exp_eta[:, None])[::-1], axis=0)[::-1]
        s1_i = cum_weighted[i]
        grad0 += X_std[i] - s1_i / s0_i

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
        ValueError: If inputs are incompatible or n_groups < 2.
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
            c_index: Concordance index value.
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
            if time[j] <= time[i] and event[j] == 1 and j != i:
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
