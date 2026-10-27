"""Survival and time-to-event analysis methods.

Provides Kaplan-Meier survival curve estimation, Cox proportional hazards
modelling, competing risks analysis (cumulative incidence), recurrent event
analysis, and time-varying covariate data preparation.

The Kaplan-Meier, competing-risks, recurrent-event and time-varying-covariate routines are pure Python
step-function estimators. The Cox proportional hazards fit (``cox_ph_model``) runs on NumPy: Newton-Raphson on
the Breslow partial likelihood with the exact gradient and information matrix.
"""

from __future__ import annotations

import math
import warnings

import numpy as np

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _mean(values: list[float]) -> float:
    """Arithmetic mean."""
    return sum(values) / len(values) if values else 0.0


def _variance(values: list[float]) -> float:
    """Sample variance."""
    n = len(values)
    if n < 2:
        return 0.0
    mu = _mean(values)
    return sum((x - mu) ** 2 for x in values) / (n - 1)


def _quantile(sorted_data: list[float], q: float) -> float:
    """Quantile from sorted data."""
    n = len(sorted_data)
    if n == 0:
        return 0.0
    idx = q * (n - 1)
    low = int(math.floor(idx))
    high = min(int(math.ceil(idx)), n - 1)
    if low == high:
        return sorted_data[low]
    frac = idx - low
    return sorted_data[low] * (1 - frac) + sorted_data[high] * frac


# ---------------------------------------------------------------------------
# Kaplan-Meier
# ---------------------------------------------------------------------------


def kaplan_meier_estimator(
    times: list[float],
    events: list[int],
    groups: list[str] | None = None,
) -> dict:
    """Kaplan-Meier survival curve estimation.

    Computes the non-parametric survival function with confidence
    intervals using Greenwood's formula.

    Args:
        times: Observed survival times for each subject.
        events: Event indicators (1 = event occurred, 0 = censored).
        groups: Optional group labels for stratified analysis. If provided,
            separate curves are computed per group.

    Returns:
        Dictionary with keys:
            - time_points: Sorted unique event times.
            - survival_probability: Step-function survival probabilities.
            - confidence_interval: List of (lower, upper) 95% CI tuples.
            - n_at_risk: Number at risk at each time point.
            - median_survival: Median survival time (or None if not reached).
            - n_events: Total number of events.
            - n_censored: Total number of censored observations.
            - groups: If stratified, dict mapping group to its own KM result.

    Raises:
        ValueError: If times and events have different lengths.
    """
    if len(times) != len(events):
        raise ValueError(f"times ({len(times)}) and events ({len(events)}) must match")

    if groups is not None:
        return _stratified_km(times, events, groups)

    # Sort by time, events first at ties
    data = sorted(zip(times, events), key=lambda x: (x[0], -x[1]))
    n = len(data)

    unique_times: list[float] = []
    surv_probs: list[float] = []
    ci_list: list[tuple[float, float]] = []
    n_at_risk_list: list[int] = []

    current_surv = 1.0
    greenwood_sum = 0.0
    at_risk = n
    total_events = 0
    total_censored = 0

    i = 0
    while i < n:
        t = data[i][0]
        # Count events and censored at this time
        d_events = 0
        d_censored = 0
        while i < n and data[i][0] == t:
            if data[i][1] == 1:
                d_events += 1
            else:
                d_censored += 1
            i += 1

        if d_events > 0:
            surv_factor = 1.0 - d_events / at_risk
            current_surv *= surv_factor

            # Greenwood variance
            if at_risk > d_events:
                greenwood_sum += d_events / (at_risk * (at_risk - d_events))

            variance = current_surv**2 * greenwood_sum
            se = math.sqrt(max(0.0, variance))

            # 95% CI using log transform
            ci_lower = max(0.0, current_surv - 1.96 * se)
            ci_upper = min(1.0, current_surv + 1.96 * se)

            unique_times.append(t)
            surv_probs.append(current_surv)
            ci_list.append((ci_lower, ci_upper))
            n_at_risk_list.append(at_risk)

            total_events += d_events

        total_censored += d_censored
        at_risk -= d_events + d_censored

    # Median survival
    median_survival = None
    for tp, sp in zip(unique_times, surv_probs):
        if sp <= 0.5:
            median_survival = tp
            break

    logger.info(
        "KM estimator: %d subjects, %d events, %d censored, median=%.2f",
        n,
        total_events,
        total_censored,
        median_survival if median_survival is not None else float("nan"),
    )

    return {
        "time_points": unique_times,
        "survival_probability": surv_probs,
        "confidence_interval": ci_list,
        "n_at_risk": n_at_risk_list,
        "median_survival": median_survival,
        "n_events": total_events,
        "n_censored": total_censored,
    }


def _stratified_km(
    times: list[float],
    events: list[int],
    groups: list[str],
) -> dict:
    """Stratified Kaplan-Meier analysis."""
    group_data: dict[str, tuple[list[float], list[int]]] = {}
    for t, e, g in zip(times, events, groups):
        if g not in group_data:
            group_data[g] = ([], [])
        group_data[g][0].append(t)
        group_data[g][1].append(e)

    group_results = {}
    for g, (g_times, g_events) in group_data.items():
        group_results[g] = kaplan_meier_estimator(g_times, g_events)

    # Aggregate
    all_times: set[float] = set()
    for gr in group_results.values():
        all_times.update(gr["time_points"])

    return {
        "time_points": sorted(all_times),
        "survival_probability": [],
        "confidence_interval": [],
        "n_at_risk": [],
        "median_survival": None,
        "n_events": sum(gr["n_events"] for gr in group_results.values()),
        "n_censored": sum(gr["n_censored"] for gr in group_results.values()),
        "groups": group_results,
    }


# ---------------------------------------------------------------------------
# Cox Proportional Hazards
# ---------------------------------------------------------------------------


_COX_MAX_ITERATIONS = 50
_COX_TOLERANCE = 1e-9  # relative log-likelihood gain below which a full Newton step ends the fit


def cox_ph_model(
    times: list[float],
    events: list[int],
    covariates: list[list[float]],
    covariate_names: list[str] | None = None,
) -> dict:
    """Fit a Cox proportional hazards model.

    Maximises the Breslow partial log-likelihood by Newton-Raphson with the exact gradient and information
    matrix (negative Hessian) and step halving; this is the fit behind
    ``SurvivalPredictor(method="cox", penalizer=0)``. Event times that are tied share one risk set (Breslow's
    approximation) and right-censoring is supported. Standard errors are the square roots of the diagonal of
    the inverse information matrix, so correlated covariates are accounted for, and the p-values are Wald tests.

    A covariate that is not identified (constant, or a linear combination of other covariates) has an infinite
    standard error and a p-value of 1. A covariate that perfectly orders the event times makes the likelihood
    monotone: its coefficient grows until the likelihood stops improving and its standard error is huge.

    Args:
        times: Observed survival times.
        events: Event indicators (1 = event, 0 = censored); booleans are accepted.
        covariates: List of covariate vectors per subject, shape
            ``(n_subjects, n_covariates)``.
        covariate_names: Optional names for covariates.

    Returns:
        Dictionary with keys:
            - coefficients: List of estimated beta coefficients.
            - hazard_ratios: List of exp(beta) values.
            - se: List of standard errors (``inf`` for a covariate that is not identified).
            - p_values: List of Wald test p-values.
            - concordance: Harrell's C-index (pairs with equal risk scores are not counted).
            - covariate_names: Names of covariates.
            - n_subjects: Number of subjects.
            - n_events: Number of events.
            - log_likelihood: Maximised partial log-likelihood.
            - converged: Whether Newton-Raphson converged.
            - n_iterations: Number of Newton-Raphson iterations.

    Raises:
        ValueError: If input dimensions are inconsistent, there is no covariate, values are not finite
            numbers, or ``events`` holds anything but 0/1.

    Warns:
        RuntimeWarning: If Newton-Raphson does not converge within 50 iterations.
    """
    n = len(times)
    if len(events) != n or len(covariates) != n:
        raise ValueError("times, events, and covariates must have same length")

    try:
        p = len(covariates[0]) if n else 0
    except TypeError as exc:  # a flat list of numbers instead of one covariate vector per subject
        raise ValueError("covariates must hold one covariate vector per subject") from exc
    if p == 0:
        raise ValueError("At least one covariate is required")

    if covariate_names is None:
        covariate_names = [f"X{i}" for i in range(p)]

    try:
        time_array = np.asarray(times, dtype=float)
        covariate_array = np.asarray(covariates, dtype=float)
        event_array = np.asarray(events)
    except (TypeError, ValueError) as exc:
        raise ValueError("times and covariates must be numeric, and every subject needs the same covariates") from exc
    if covariate_array.shape != (n, p) or time_array.shape != (n,) or event_array.shape != (n,):
        raise ValueError("times and events must be one value per subject and covariates a subjects x covariates table")
    if not (np.all(np.isfinite(time_array)) and np.all(np.isfinite(covariate_array))):
        raise ValueError("times and covariates must be finite")
    if not np.all(np.isin(event_array, (0, 1))):
        raise ValueError("events must contain only 0 (censored) and 1 (event), or booleans")
    occurred = event_array.astype(bool)

    # The Cox fit is shared with SurvivalPredictor; imported here because models builds on this package's layer.
    from metainformant.life_events.models.statistical_models import _fit_cox

    fit = _fit_cox(covariate_array, time_array, occurred, 0.0, _COX_MAX_ITERATIONS, _COX_TOLERANCE)
    beta = fit["coefficients"]
    se = fit["std_errors"]
    if not fit["converged"]:
        warnings.warn(
            f"The Cox model did not converge in {_COX_MAX_ITERATIONS} iterations; the estimates are unreliable",
            RuntimeWarning,
            stacklevel=2,
        )

    p_values = [_wald_p_value(float(b), float(s)) for b, s in zip(beta, se)]
    with np.errstate(over="ignore"):
        hazard_ratios = np.exp(beta).tolist()  # inf rather than an OverflowError for a huge coefficient

    concordance = _compute_concordance(time_array, occurred, covariate_array, beta)

    logger.info(
        "Cox PH: %d subjects, %d events, %d covariates, C=%.3f",
        n,
        int(occurred.sum()),
        p,
        concordance,
    )

    return {
        "coefficients": beta.tolist(),
        "hazard_ratios": hazard_ratios,
        "se": se.tolist(),
        "p_values": p_values,
        "concordance": concordance,
        "covariate_names": covariate_names,
        "n_subjects": n,
        "n_events": int(occurred.sum()),
        "log_likelihood": float(fit["log_likelihood"]),
        "converged": bool(fit["converged"]),
        "n_iterations": int(fit["n_iter"]),
    }


def _wald_p_value(coefficient: float, standard_error: float) -> float:
    """Two-sided Wald test p-value; 1 when the standard error is zero or infinite (nothing can be tested)."""
    if 0.0 < standard_error < float("inf"):
        return math.erfc(abs(coefficient / standard_error) / math.sqrt(2.0))
    return 1.0


def _compute_concordance(
    times: np.ndarray,
    events: np.ndarray,
    covariates: np.ndarray,
    beta: np.ndarray,
) -> float:
    """Compute Harrell's concordance index.

    A pair is comparable when the subject with the event has the earlier time. It is concordant when that
    subject also has the higher risk score ``x . beta`` and discordant when the risk score is lower; pairs
    with equal risk scores are counted in neither group. Event rows are processed in blocks, so the
    comparisons are vectorised without holding all n^2 pairs in memory.
    """
    times = np.asarray(times, dtype=float)
    risk = np.asarray(covariates, dtype=float) @ np.asarray(beta, dtype=float)
    event_rows = np.flatnonzero(np.asarray(events) == 1)

    concordant = 0
    discordant = 0
    block = max(1, 2_000_000 // max(1, times.size))
    for start in range(0, event_rows.size, block):
        rows = event_rows[start : start + block]
        later = times[None, :] > times[rows][:, None]
        concordant += int(np.count_nonzero(later & (risk[rows][:, None] > risk[None, :])))
        discordant += int(np.count_nonzero(later & (risk[rows][:, None] < risk[None, :])))

    total = concordant + discordant
    return concordant / total if total > 0 else 0.5


# ---------------------------------------------------------------------------
# Competing Risks
# ---------------------------------------------------------------------------


def competing_risks(
    times: list[float],
    events: list[int],
    event_types: list[int],
) -> dict:
    """Competing risks analysis using cumulative incidence functions.

    Estimates the cumulative incidence of each event type accounting for
    other event types as competing risks using the Aalen-Johansen estimator.

    Args:
        times: Observed times for each subject.
        events: Event indicator (1 = event, 0 = censored).
        event_types: Type of event (integer) for subjects with events.
            Ignored when event == 0.

    Returns:
        Dictionary with keys:
            - cumulative_incidence_per_type: Dict mapping event type to
              list of (time, cumulative_incidence) tuples.
            - cause_specific_hazards: Dict mapping event type to
              list of (time, hazard) tuples.
            - overall_survival: List of (time, survival) tuples.
    """
    n = len(times)
    if len(events) != n or len(event_types) != n:
        raise ValueError("All inputs must have same length")

    # Sort by time
    order = sorted(range(n), key=lambda i: (times[i], -events[i]))
    t_sorted = [times[i] for i in order]
    e_sorted = [events[i] for i in order]
    et_sorted = [event_types[i] for i in order]

    unique_types = sorted(set(et_sorted[i] for i in range(n) if e_sorted[i] == 1))

    # Compute overall KM survival
    at_risk = n
    surv = 1.0

    ci_per_type: dict[int, list[tuple[float, float]]] = {t: [] for t in unique_types}
    hazard_per_type: dict[int, list[tuple[float, float]]] = {t: [] for t in unique_types}
    overall_surv: list[tuple[float, float]] = []
    cumulative_inc: dict[int, float] = {t: 0.0 for t in unique_types}

    i = 0
    while i < n and at_risk > 0:
        current_time = t_sorted[i]

        # Count events by type and censored at this time
        type_counts: dict[int, int] = {}
        n_censored = 0
        j = i
        while j < n and t_sorted[j] == current_time:
            if e_sorted[j] == 1:
                etype = et_sorted[j]
                type_counts[etype] = type_counts.get(etype, 0) + 1
            else:
                n_censored += 1
            j += 1

        total_events = sum(type_counts.values())

        if total_events > 0:
            # Cause-specific hazards
            for etype in unique_types:
                d_k = type_counts.get(etype, 0)
                if d_k > 0:
                    h_k = d_k / at_risk
                    hazard_per_type[etype].append((current_time, h_k))
                    # Cumulative incidence increment
                    cumulative_inc[etype] += surv * h_k

                ci_per_type[etype].append((current_time, cumulative_inc[etype]))

            # Update overall survival
            surv *= 1.0 - total_events / at_risk

        overall_surv.append((current_time, surv))
        at_risk -= total_events + n_censored
        i = j

    logger.info(
        "Competing risks: %d subjects, %d event types, %d events",
        n,
        len(unique_types),
        sum(e_sorted),
    )

    return {
        "cumulative_incidence_per_type": ci_per_type,
        "cause_specific_hazards": hazard_per_type,
        "overall_survival": overall_surv,
    }


# ---------------------------------------------------------------------------
# Recurrent Events
# ---------------------------------------------------------------------------


def recurrent_events(
    times: list[list[float]],
    events: list[list[int]],
    subject_ids: list[str],
) -> dict:
    """Analyse recurrent event patterns.

    Computes the mean cumulative function (MCF), per-subject event rates,
    and gap time distribution.

    Args:
        times: List of time vectors per subject (event and censoring times).
        events: List of event indicator vectors per subject.
        subject_ids: List of subject identifiers.

    Returns:
        Dictionary with keys:
            - mean_cumulative_function: List of (time, mcf) tuples.
            - rate_per_subject: Dict mapping subject to event rate.
            - gap_time_distribution: Dict with mean, median, std of gap times.
            - n_subjects: Number of subjects.
            - total_events: Total events across all subjects.
    """
    n = len(subject_ids)
    if len(times) != n or len(events) != n:
        raise ValueError("All inputs must have same length")

    # Collect all event times and per-subject info
    all_event_times: list[float] = []
    gap_times: list[float] = []
    rate_per_subject: dict[str, float] = {}
    total_events = 0

    for i in range(n):
        subj_times = times[i]
        subj_events = events[i]
        subj_event_times = sorted(subj_times[j] for j in range(len(subj_times)) if subj_events[j] == 1)
        n_events_subj = len(subj_event_times)
        total_events += n_events_subj

        max_time = max(subj_times) if subj_times else 0
        rate = n_events_subj / max_time if max_time > 0 else 0.0
        rate_per_subject[subject_ids[i]] = rate

        all_event_times.extend(subj_event_times)

        # Gap times
        for j in range(1, len(subj_event_times)):
            gap = subj_event_times[j] - subj_event_times[j - 1]
            gap_times.append(gap)

    # Mean Cumulative Function (Nelson-Aalen style for recurrent events)
    all_event_times.sort()
    unique_times = sorted(set(all_event_times))

    # At each time, count events across subjects / number at risk
    mcf_entries: list[tuple[float, float]] = []
    cumulative = 0.0

    for t in unique_times:
        # Count events at this time
        n_events_at_t = sum(1 for et in all_event_times if et == t)
        # Number at risk = subjects with max observation time >= t
        n_at_risk = sum(1 for i in range(n) if max(times[i]) >= t)
        if n_at_risk > 0:
            cumulative += n_events_at_t / n_at_risk
        mcf_entries.append((t, cumulative))

    # Gap time distribution
    gap_dist: dict[str, float] = {}
    if gap_times:
        sorted_gaps = sorted(gap_times)
        gap_dist = {
            "mean": _mean(gap_times),
            "median": _quantile(sorted_gaps, 0.5),
            "std": math.sqrt(_variance(gap_times)) if len(gap_times) > 1 else 0.0,
            "min": min(gap_times),
            "max": max(gap_times),
            "n_gaps": len(gap_times),
        }

    logger.info(
        "Recurrent events: %d subjects, %d total events",
        n,
        total_events,
    )

    return {
        "mean_cumulative_function": mcf_entries,
        "rate_per_subject": rate_per_subject,
        "gap_time_distribution": gap_dist,
        "n_subjects": n,
        "total_events": total_events,
    }


# ---------------------------------------------------------------------------
# Time-Varying Covariates
# ---------------------------------------------------------------------------


def time_varying_covariates(
    intervals: list[dict],
) -> dict:
    """Prepare data for time-varying covariate survival analysis.

    Expands subject-level data into counting-process format where each
    row represents an interval (start, stop] with the covariate values
    active during that interval.

    Args:
        intervals: List of interval dicts, each with keys:
            - subject_id: Subject identifier.
            - start: Interval start time.
            - stop: Interval end time.
            - event: Event indicator for this interval (0 or 1).
            - covariates: Dict of covariate values during this interval.

    Returns:
        Dictionary with keys:
            - expanded_data: List of processed interval dicts ready for
              Cox model, each with start, stop, event, subject_id, and
              flattened covariate values.
            - interval_counts: Dict mapping subject to number of intervals.
            - n_subjects: Number of unique subjects.
            - n_intervals: Total number of intervals.
            - covariate_names: List of covariate names encountered.
    """
    expanded: list[dict] = []
    interval_counts: dict[str, int] = {}
    all_covariates: set[str] = set()
    subjects: set[str] = set()

    for interval in intervals:
        subject = interval["subject_id"]
        subjects.add(subject)
        interval_counts[subject] = interval_counts.get(subject, 0) + 1

        covs = interval.get("covariates", {})
        all_covariates.update(covs.keys())

        entry = {
            "subject_id": subject,
            "start": interval["start"],
            "stop": interval["stop"],
            "event": interval.get("event", 0),
        }
        entry.update(covs)
        expanded.append(entry)

    # Sort by subject then start time
    expanded.sort(key=lambda x: (x["subject_id"], x["start"]))

    logger.info(
        "Time-varying covariates: %d subjects, %d intervals, %d covariates",
        len(subjects),
        len(expanded),
        len(all_covariates),
    )

    return {
        "expanded_data": expanded,
        "interval_counts": interval_counts,
        "n_subjects": len(subjects),
        "n_intervals": len(expanded),
        "covariate_names": sorted(all_covariates),
    }
