"""Comprehensive QC reporting, threshold checking, and trend analysis.

Provides functions for generating quality control reports from collected
metrics, aggregating multi-sample QC data, checking metrics against
configurable thresholds, and analysing QC metric trends over time or
across batches.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Default QC thresholds
# ---------------------------------------------------------------------------


def default_qc_thresholds() -> dict:
    """Return default QC thresholds for common sequencing metrics.

    Each metric maps to a dict with ``warn`` and ``fail`` values.
    For metrics where higher is better (e.g. mapping_rate), thresholds
    are lower bounds. For metrics where lower is better (e.g.
    duplication_rate), thresholds are upper bounds.

    Returns:
        Dictionary mapping metric name to threshold config with keys:
            - warn: Warning threshold value.
            - fail: Failure threshold value.
            - direction: ``"min"`` (value must be >= threshold) or
              ``"max"`` (value must be <= threshold).
            - description: Human-readable description.
    """
    return {
        "total_reads": {
            "warn": 1_000_000,
            "fail": 500_000,
            "direction": "min",
            "description": "Total number of reads",
        },
        "mapping_rate": {
            "warn": 0.85,
            "fail": 0.70,
            "direction": "min",
            "description": "Fraction of reads mapped to reference",
        },
        "duplication_rate": {
            "warn": 0.30,
            "fail": 0.50,
            "direction": "max",
            "description": "Fraction of duplicate reads",
        },
        "mean_quality": {
            "warn": 28.0,
            "fail": 20.0,
            "direction": "min",
            "description": "Mean base quality score (Phred)",
        },
        "gc_content": {
            "warn": 0.55,
            "fail": 0.65,
            "direction": "max",
            "description": "GC content fraction",
        },
        "contamination_rate": {
            "warn": 0.02,
            "fail": 0.05,
            "direction": "max",
            "description": "Estimated contamination fraction",
        },
        "adapter_content": {
            "warn": 0.05,
            "fail": 0.15,
            "direction": "max",
            "description": "Fraction of reads with adapter contamination",
        },
        "mean_insert_size": {
            "warn": 150.0,
            "fail": 100.0,
            "direction": "min",
            "description": "Mean insert size for paired-end libraries",
        },
        "on_target_rate": {
            "warn": 0.70,
            "fail": 0.50,
            "direction": "min",
            "description": "Fraction of reads on target (capture/amplicon)",
        },
        "mean_coverage": {
            "warn": 20.0,
            "fail": 10.0,
            "direction": "min",
            "description": "Mean sequencing depth across target regions",
        },
        "percent_bases_q30": {
            "warn": 0.80,
            "fail": 0.60,
            "direction": "min",
            "description": "Fraction of bases with quality >= Q30",
        },
        "percent_aligned": {
            "warn": 0.90,
            "fail": 0.75,
            "direction": "min",
            "description": "Percent of reads aligned",
        },
    }


# ---------------------------------------------------------------------------
# Threshold checking
# ---------------------------------------------------------------------------


def check_qc_thresholds(
    metrics: dict,
    thresholds: dict | None = None,
) -> dict:
    """Check QC metrics against configurable thresholds.

    Args:
        metrics: Dictionary mapping metric name to numeric value.
        thresholds: Optional threshold configuration. If ``None``, uses
            :func:`default_qc_thresholds`. Each metric entry should have
            ``warn``, ``fail``, and ``direction`` keys.

    Returns:
        Dictionary with keys:
            - status: Overall status ``"pass"`` | ``"warn"`` | ``"fail"``.
            - failed_metrics: List of dicts for metrics that failed.
            - warnings: List of dicts for metrics with warnings.
            - passed_metrics: List of metric names that passed.
            - n_checked: Number of metrics checked.
    """
    if thresholds is None:
        thresholds = default_qc_thresholds()

    failed: list[dict] = []
    warnings: list[dict] = []
    passed: list[str] = []

    for metric_name, value in metrics.items():
        if metric_name not in thresholds:
            passed.append(metric_name)
            continue

        thresh = thresholds[metric_name]
        direction = thresh.get("direction", "min")
        warn_val = thresh["warn"]
        fail_val = thresh["fail"]

        entry = {
            "metric": metric_name,
            "value": value,
            "warn_threshold": warn_val,
            "fail_threshold": fail_val,
            "direction": direction,
            "description": thresh.get("description", ""),
        }

        if direction == "min":
            if value < fail_val:
                entry["status"] = "fail"
                failed.append(entry)
            elif value < warn_val:
                entry["status"] = "warn"
                warnings.append(entry)
            else:
                passed.append(metric_name)
        else:  # direction == "max"
            if value > fail_val:
                entry["status"] = "fail"
                failed.append(entry)
            elif value > warn_val:
                entry["status"] = "warn"
                warnings.append(entry)
            else:
                passed.append(metric_name)

    if failed:
        overall = "fail"
    elif warnings:
        overall = "warn"
    else:
        overall = "pass"

    logger.info(
        "QC threshold check: %s (pass=%d, warn=%d, fail=%d)",
        overall,
        len(passed),
        len(warnings),
        len(failed),
    )

    return {
        "status": overall,
        "failed_metrics": failed,
        "warnings": warnings,
        "passed_metrics": passed,
        "n_checked": len(metrics),
    }


# ---------------------------------------------------------------------------
# Sample aggregation
# ---------------------------------------------------------------------------


def _mean(values: list[float]) -> float:
    """Arithmetic mean."""
    return sum(values) / len(values) if values else 0.0


def _std(values: list[float]) -> float:
    """Sample standard deviation."""
    n = len(values)
    if n < 2:
        return 0.0
    mu = _mean(values)
    return math.sqrt(sum((x - mu) ** 2 for x in values) / (n - 1))


def _median(values: list[float]) -> float:
    """Median of sorted values."""
    s = sorted(values)
    n = len(s)
    if n == 0:
        return 0.0
    if n % 2 == 1:
        return s[n // 2]
    return (s[n // 2 - 1] + s[n // 2]) / 2.0


def _iqr(values: list[float]) -> tuple[float, float]:
    """Interquartile range."""
    s = sorted(values)
    n = len(s)
    q1_idx = n // 4
    q3_idx = 3 * n // 4
    return s[q1_idx], s[min(q3_idx, n - 1)]


def aggregate_sample_qc(sample_metrics: list[dict]) -> dict:
    """Aggregate QC metrics across multiple samples.

    Computes summary statistics (mean, median, std) for each metric across
    samples and identifies outlier samples using the IQR method.

    Args:
        sample_metrics: List of dicts, each mapping metric name to value
            for one sample. Must also include ``"sample_id"`` key.

    Returns:
        Dictionary with keys:
            - summary_table: Dict mapping metric to {mean, median, std, min, max}.
            - outlier_samples: List of dicts describing outlier observations.
            - overall_pass_rate: Fraction of samples passing all default thresholds.
            - n_samples: Number of samples.
    """
    if not sample_metrics:
        return {
            "summary_table": {},
            "outlier_samples": [],
            "overall_pass_rate": 0.0,
            "n_samples": 0,
        }

    # Collect all metric keys (excluding sample_id)
    all_keys: set[str] = set()
    for sm in sample_metrics:
        all_keys.update(k for k in sm if k != "sample_id")

    summary_table: dict[str, dict[str, float]] = {}
    outliers: list[dict] = []

    for key in sorted(all_keys):
        values = []
        sample_ids = []
        for sm in sample_metrics:
            if key in sm:
                try:
                    v = float(sm[key])
                    if not math.isnan(v):
                        values.append(v)
                        sample_ids.append(sm.get("sample_id", "unknown"))
                except (TypeError, ValueError):
                    continue

        if not values:
            continue

        mu = _mean(values)
        sd = _std(values)
        med = _median(values)

        summary_table[key] = {
            "mean": mu,
            "median": med,
            "std": sd,
            "min": min(values),
            "max": max(values),
            "n_values": len(values),
        }

        # Detect outliers using IQR
        if len(values) >= 4:
            q1, q3 = _iqr(values)
            iqr_val = q3 - q1
            lower = q1 - 1.5 * iqr_val
            upper = q3 + 1.5 * iqr_val

            for val, sid in zip(values, sample_ids):
                if val < lower or val > upper:
                    outliers.append(
                        {
                            "sample_id": sid,
                            "metric": key,
                            "value": val,
                            "expected_range": (lower, upper),
                            "direction": "low" if val < lower else "high",
                        }
                    )

    # Pass rate
    thresholds = default_qc_thresholds()
    n_pass = 0
    for sm in sample_metrics:
        result = check_qc_thresholds({k: v for k, v in sm.items() if k != "sample_id"}, thresholds)
        if result["status"] == "pass":
            n_pass += 1

    pass_rate = n_pass / len(sample_metrics)

    logger.info(
        "Aggregated QC for %d samples: %d metrics, %d outliers, pass_rate=%.3f",
        len(sample_metrics),
        len(summary_table),
        len(outliers),
        pass_rate,
    )

    return {
        "summary_table": summary_table,
        "outlier_samples": outliers,
        "overall_pass_rate": pass_rate,
        "n_samples": len(sample_metrics),
    }


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def generate_qc_report(
    metrics: dict,
    output_path: str | None = None,
) -> dict:
    """Generate a comprehensive QC report from collected metrics.

    Combines threshold checking, summary statistics, and a narrative
    summary into a structured report.

    Args:
        metrics: Dictionary mapping metric name to value.
        output_path: Optional file path to write JSON report.

    Returns:
        Dictionary with keys:
            - summary: Narrative summary string.
            - pass_fail_status: Overall ``"pass"`` | ``"warn"`` | ``"fail"``.
            - warnings: List of warning descriptions.
            - threshold_results: Detailed threshold check results.
            - metrics: Input metrics echoed back.
            - report_path: Path where report was written (or None).
    """
    threshold_results = check_qc_thresholds(metrics)

    # Build narrative summary
    status = threshold_results["status"]
    n_pass = len(threshold_results["passed_metrics"])
    n_warn = len(threshold_results["warnings"])
    n_fail = len(threshold_results["failed_metrics"])

    lines = [f"QC Report: Overall status = {status.upper()}"]
    lines.append(f"Metrics checked: {threshold_results['n_checked']}")
    lines.append(f"Passed: {n_pass}, Warnings: {n_warn}, Failed: {n_fail}")

    if threshold_results["failed_metrics"]:
        lines.append("\nFailed metrics:")
        for fm in threshold_results["failed_metrics"]:
            lines.append(
                f"  - {fm['metric']}: {fm['value']} " f"(threshold: {fm['fail_threshold']}, {fm['direction']})"
            )

    if threshold_results["warnings"]:
        lines.append("\nWarning metrics:")
        for wm in threshold_results["warnings"]:
            lines.append(
                f"  - {wm['metric']}: {wm['value']} " f"(threshold: {wm['warn_threshold']}, {wm['direction']})"
            )

    summary = "\n".join(lines)
    warning_texts = [f"{w['metric']}: {w['value']} ({w.get('description', '')})" for w in threshold_results["warnings"]]

    report = {
        "summary": summary,
        "pass_fail_status": status,
        "warnings": warning_texts,
        "threshold_results": threshold_results,
        "metrics": metrics,
        "report_path": None,
    }

    if output_path is not None:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2, default=str))
        report["report_path"] = str(out)
        logger.info("QC report written to %s", out)

    logger.info("Generated QC report: status=%s", status)
    return report


# ---------------------------------------------------------------------------
# Trend analysis
# ---------------------------------------------------------------------------


def qc_trend_analysis(
    metrics_over_time: list[dict],
    metric_name: str,
) -> dict:
    """Analyse QC metric trends over time or across batches.

    Performs simple linear regression to detect drift and identifies
    change points where the metric shifts significantly.

    Args:
        metrics_over_time: List of dicts, each with at least the key
            ``metric_name`` and optionally ``"batch"`` or ``"date"``.
            Ordered chronologically.
        metric_name: Name of the metric to analyse.

    Returns:
        Dictionary with keys:
            - trend: ``"stable"`` | ``"increasing"`` | ``"decreasing"``.
            - is_drifting: ``True`` if trend is statistically significant.
            - slope: Linear regression slope.
            - r_squared: Coefficient of determination.
            - change_points: List of indices where shifts detected.
            - summary: Narrative description.
            - values: Extracted metric values.
    """
    values = []
    for entry in metrics_over_time:
        if metric_name in entry:
            try:
                values.append(float(entry[metric_name]))
            except (TypeError, ValueError):
                continue

    n = len(values)
    if n < 3:
        return {
            "trend": "insufficient_data",
            "is_drifting": False,
            "slope": 0.0,
            "r_squared": 0.0,
            "change_points": [],
            "summary": f"Insufficient data points ({n}) for trend analysis",
            "values": values,
        }

    # Simple linear regression: y = a + b*x
    x = list(range(n))
    x_mean = _mean(x)
    y_mean = _mean(values)

    ss_xx = sum((xi - x_mean) ** 2 for xi in x)
    ss_xy = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(x, values))

    if ss_xx == 0:
        slope = 0.0
    else:
        slope = ss_xy / ss_xx

    intercept = y_mean - slope * x_mean

    # R-squared
    ss_res = sum((yi - (intercept + slope * xi)) ** 2 for xi, yi in zip(x, values))
    ss_tot = sum((yi - y_mean) ** 2 for yi in values)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    # Significance test for slope
    if n > 2 and ss_xx > 0:
        mse = ss_res / (n - 2)
        se_slope = math.sqrt(mse / ss_xx) if mse > 0 else 0.0
        t_stat = abs(slope / se_slope) if se_slope > 0 else 0.0
        # Rough p-value approximation
        is_drifting = t_stat > 2.0  # Approximately p < 0.05 for moderate n
    else:
        is_drifting = False

    # Determine trend direction
    if abs(slope) < 1e-10:
        trend = "stable"
    elif slope > 0:
        trend = "increasing"
    else:
        trend = "decreasing"

    # Simple change point detection: find points where deviation from
    # running mean exceeds 2 SDs
    change_points: list[int] = []
    if n >= 5:
        window = max(3, n // 5)
        overall_std = _std(values)
        if overall_std > 0:
            for i in range(window, n):
                local_mean = _mean(values[i - window : i])
                if abs(values[i] - local_mean) > 2.0 * overall_std:
                    change_points.append(i)

    summary = f"Metric '{metric_name}' shows a {trend} trend " f"(slope={slope:.4f}, R^2={r_squared:.3f}). "
    if is_drifting:
        summary += "The drift is statistically significant. "
    else:
        summary += "No significant drift detected. "
    if change_points:
        summary += f"Change points detected at indices: {change_points}."
    else:
        summary += "No abrupt change points detected."

    logger.info(
        "Trend analysis for '%s': %s, slope=%.4f, R2=%.3f, %d change points",
        metric_name,
        trend,
        slope,
        r_squared,
        len(change_points),
    )

    return {
        "trend": trend,
        "is_drifting": is_drifting,
        "slope": slope,
        "r_squared": r_squared,
        "change_points": change_points,
        "summary": summary,
        "values": values,
    }
