"""Report generation for long-read analysis pipelines.

Generates structured QC reports, assembly quality reports, methylation
analysis reports, and cross-pipeline run summaries. Reports can be
exported as JSON, plain text, or HTML.

All report generators consume PipelineResult objects and produce
structured report dataclasses or dictionaries.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


@dataclass
class QCReport:
    """Structured QC report from a QC pipeline run.

    Attributes:
        sample_name: Name of the sequencing sample.
        total_reads: Total number of input reads.
        total_bases: Sum of all read lengths.
        n50: N50 read length.
        mean_length: Mean read length.
        mean_quality: Mean Phred quality score across reads.
        reads_passed_filter: Number of reads passing all filters.
        filter_pass_rate: Fraction of reads passing all filters.
        adapter_trimmed: Number of reads with adapters trimmed.
        length_distribution: Read length distribution statistics.
        quality_distribution: Quality score distribution statistics.
        timestamp: ISO-format timestamp when the report was generated.
    """

    sample_name: str = ""
    total_reads: int = 0
    total_bases: int = 0
    n50: int = 0
    mean_length: float = 0.0
    mean_quality: float = 0.0
    reads_passed_filter: int = 0
    filter_pass_rate: float = 0.0
    adapter_trimmed: int = 0
    length_distribution: dict[str, Any] = field(default_factory=dict)
    quality_distribution: dict[str, Any] = field(default_factory=dict)
    timestamp: str = ""


def generate_qc_report(pipeline_result: Any) -> QCReport:
    """Generate a structured QC report from a QC pipeline result.

    Extracts metrics from the pipeline summary and completed steps to
    build a comprehensive QC report.

    Args:
        pipeline_result: A PipelineResult from a QC pipeline run.

    Returns:
        QCReport with all available metrics populated.
    """
    summary = getattr(pipeline_result, "summary", {})
    if not isinstance(summary, dict):
        summary = {}

    steps = getattr(pipeline_result, "steps", [])
    output_dir = getattr(pipeline_result, "output_dir", Path("."))

    # Extract metrics from summary or step results
    metrics = summary.get("compute_metrics", summary)

    length_stats = metrics.get("length_stats", None)
    quality_dist = metrics.get("quality_distribution", None)

    total_reads = metrics.get("total_reads_input", 0)
    reads_passed = metrics.get("reads_after_processing", 0)

    # Build length distribution dict from dataclass if present
    length_distribution: dict[str, Any] = {}
    if length_stats is not None:
        if hasattr(length_stats, "count"):
            length_distribution = {
                "count": length_stats.count,
                "total_bases": length_stats.total_bases,
                "mean_length": length_stats.mean_length,
                "median_length": length_stats.median_length,
                "min_length": length_stats.min_length,
                "max_length": length_stats.max_length,
                "std_dev": length_stats.std_dev,
                "n50": length_stats.n50,
                "n90": length_stats.n90,
                "l50": length_stats.l50,
                "percentiles": length_stats.percentiles,
            }
        elif isinstance(length_stats, dict):
            length_distribution = length_stats

    # Build quality distribution dict
    quality_distribution: dict[str, Any] = {}
    if quality_dist is not None:
        if hasattr(quality_dist, "mean_quality"):
            quality_distribution = {
                "mean_quality": quality_dist.mean_quality,
                "median_quality": quality_dist.median_quality,
                "min_quality": quality_dist.min_quality,
                "max_quality": quality_dist.max_quality,
                "q7_fraction": quality_dist.q7_fraction,
                "q10_fraction": quality_dist.q10_fraction,
                "q15_fraction": quality_dist.q15_fraction,
                "q20_fraction": quality_dist.q20_fraction,
            }
        elif isinstance(quality_dist, dict):
            quality_distribution = quality_dist

    # Extract N50 and mean length from the right source
    n50 = 0
    mean_length = 0.0
    total_bases = 0
    mean_quality = 0.0

    if length_distribution:
        n50 = length_distribution.get("n50", 0)
        mean_length = length_distribution.get("mean_length", 0.0)
        total_bases = length_distribution.get("total_bases", 0)

    if quality_distribution:
        mean_quality = quality_distribution.get("mean_quality", 0.0)

    # Count adapter-trimmed reads
    adapter_trimmed = 0
    for step in steps:
        if step.name == "trim_adapters" and step.result is not None:
            if isinstance(step.result, list):
                adapter_trimmed = sum(
                    1
                    for r in step.result
                    if (isinstance(r, dict) and r.get("metadata", {}).get("trimmed_start", 0) > 0)
                    or (hasattr(r, "metadata") and r.metadata.get("trimmed_start", 0) > 0)
                )

    # Calculate filter pass rate
    filter_pass_rate = reads_passed / total_reads if total_reads > 0 else 0.0

    # Infer sample name from output directory
    sample_name = str(output_dir.name) if output_dir != Path(".") else "unknown_sample"

    timestamp = datetime.now(timezone.utc).isoformat()

    report = QCReport(
        sample_name=sample_name,
        total_reads=total_reads,
        total_bases=total_bases,
        n50=n50,
        mean_length=mean_length,
        mean_quality=mean_quality,
        reads_passed_filter=reads_passed,
        filter_pass_rate=filter_pass_rate,
        adapter_trimmed=adapter_trimmed,
        length_distribution=length_distribution,
        quality_distribution=quality_distribution,
        timestamp=timestamp,
    )

    logger.info(
        "Generated QC report: %d reads, N50=%d, mean_quality=%.1f",
        total_reads,
        n50,
        mean_quality,
    )
    return report


def generate_assembly_report(pipeline_result: Any) -> dict[str, Any]:
    """Generate an assembly quality report from an assembly pipeline result.

    Extracts consensus quality metrics, overlap graph statistics, and
    input read characteristics.

    Args:
        pipeline_result: A PipelineResult from an assembly pipeline run.

    Returns:
        Dictionary with assembly report metrics.
    """
    summary = getattr(pipeline_result, "summary", {})
    if not isinstance(summary, dict):
        summary = {}

    steps = getattr(pipeline_result, "steps", [])

    # Get stats from the calculate_stats step
    stats = summary.get("calculate_stats", {})
    if not stats:
        for step in steps:
            if step.name == "calculate_stats" and step.result is not None:
                stats = step.result if isinstance(step.result, dict) else {}
                break

    # Get polished consensus info
    polished = None
    for step in steps:
        if step.name == "polish" and step.result is not None:
            polished = step.result
            break

    consensus_info: dict[str, Any] = {}
    if polished is not None:
        if hasattr(polished, "sequence"):
            consensus_info = {
                "length": polished.length,
                "num_reads_used": polished.num_reads,
                "mean_quality": polished.mean_quality,
                "mean_coverage": polished.mean_coverage,
            }
        elif isinstance(polished, dict):
            consensus_info = polished

    report: dict[str, Any] = {
        "pipeline": "assembly",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "consensus": consensus_info,
        "statistics": stats,
        "success": getattr(pipeline_result, "success", False),
        "total_duration": getattr(pipeline_result, "total_duration", 0.0),
        "steps": [
            {
                "name": step.name,
                "status": step.status,
                "duration_seconds": step.duration_seconds,
                "error": step.error,
            }
            for step in steps
        ],
    }

    logger.info(
        "Generated assembly report: consensus_length=%d",
        consensus_info.get("length", 0),
    )
    return report


def generate_methylation_report(pipeline_result: Any) -> dict[str, Any]:
    """Generate a methylation analysis report from pipeline results.

    Summarizes modification calls, regional aggregation, and differential
    methylation findings.

    Args:
        pipeline_result: A PipelineResult from a methylation pipeline run.

    Returns:
        Dictionary with methylation analysis metrics.
    """
    summary = getattr(pipeline_result, "summary", {})
    if not isinstance(summary, dict):
        summary = {}

    steps = getattr(pipeline_result, "steps", [])

    # Extract call counts
    num_5mc = summary.get("5mc_calls", 0)
    num_6ma = summary.get("6ma_calls", 0)
    num_regions = summary.get("aggregated_regions", 0)

    # Extract differential analysis results
    diff_data = summary.get("differential", {})
    if not isinstance(diff_data, dict):
        diff_data = {}

    # Get actual call data from steps if available
    five_mc_calls: list[Any] = []
    six_ma_calls: list[Any] = []
    aggregated: list[Any] = []

    for step in steps:
        if step.name == "call_5mc" and step.result is not None:
            five_mc_calls = step.result if isinstance(step.result, list) else []
            num_5mc = len(five_mc_calls)
        elif step.name == "call_6ma" and step.result is not None:
            six_ma_calls = step.result if isinstance(step.result, list) else []
            num_6ma = len(six_ma_calls)
        elif step.name == "aggregate_regions" and step.result is not None:
            aggregated = step.result if isinstance(step.result, list) else []
            num_regions = len(aggregated)
        elif step.name == "differential_analysis" and step.result is not None:
            if isinstance(step.result, dict):
                diff_data = step.result

    # Summarize modification calls
    modification_summary: dict[str, Any] = {
        "5mC": {
            "num_calls": num_5mc,
        },
        "6mA": {
            "num_calls": num_6ma,
        },
    }

    # Compute average methylation from calls
    for calls, key in [(five_mc_calls, "5mC"), (six_ma_calls, "6mA")]:
        if calls:
            probs = []
            for c in calls:
                prob = c.probability if hasattr(c, "probability") else c.get("probability", 0.0)
                probs.append(prob)
            if probs:
                modification_summary[key]["mean_probability"] = sum(probs) / len(probs)

    # Summarize aggregated regions
    region_summary: list[dict[str, Any]] = []
    for reg in aggregated:
        if hasattr(reg, "chromosome"):
            region_summary.append(
                {
                    "chromosome": reg.chromosome,
                    "start": reg.start,
                    "end": reg.end,
                    "name": reg.name,
                    "mean_methylation": reg.mean_methylation,
                    "num_cpgs": reg.num_cpgs,
                    "mean_coverage": reg.mean_coverage,
                }
            )
        elif isinstance(reg, dict):
            region_summary.append(reg)

    report: dict[str, Any] = {
        "pipeline": "methylation",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "modifications": modification_summary,
        "total_calls": num_5mc + num_6ma,
        "aggregated_regions": {
            "count": num_regions,
            "regions": region_summary[:50],  # Limit to 50 for report size
        },
        "differential_methylation": {
            "num_tested": diff_data.get("num_tested", 0),
            "num_significant": diff_data.get("num_significant", 0),
        },
        "success": getattr(pipeline_result, "success", False),
        "total_duration": getattr(pipeline_result, "total_duration", 0.0),
        "steps": [
            {
                "name": step.name,
                "status": step.status,
                "duration_seconds": step.duration_seconds,
                "error": step.error,
            }
            for step in steps
        ],
    }

    logger.info(
        "Generated methylation report: %d 5mC calls, %d 6mA calls, %d regions",
        num_5mc,
        num_6ma,
        num_regions,
    )
    return report


def export_report(
    report: QCReport | dict[str, Any],
    output_path: Path | str,
    format: str = "json",
) -> Path:
    """Export a report to a file in the specified format.

    Supports JSON, plain text, and HTML output formats.

    Args:
        report: A QCReport dataclass or a dictionary report.
        output_path: Path for the output file.
        format: Output format - 'json', 'text', or 'html'.

    Returns:
        Path to the written report file.

    Raises:
        ValueError: If format is not supported.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert to dict if it is a dataclass
    if isinstance(report, QCReport):
        report_dict = asdict(report)
    elif isinstance(report, dict):
        report_dict = report
    else:
        raise TypeError(f"Report must be QCReport or dict, got {type(report).__name__}")

    if format == "json":
        _write_json_report(report_dict, output_path)
    elif format == "text":
        _write_text_report(report_dict, output_path)
    elif format == "html":
        _write_html_report(report_dict, output_path)
    else:
        raise ValueError(f"Unsupported report format '{format}'. Use 'json', 'text', or 'html'.")

    logger.info("Exported report to %s (format=%s)", output_path, format)
    return output_path


def generate_run_summary(pipeline_results: list[Any]) -> dict[str, Any]:
    """Generate a summary across multiple pipeline runs.

    Aggregates results from multiple pipeline executions into a single
    overview, useful for batch processing or multi-sample analyses.

    Args:
        pipeline_results: List of PipelineResult objects from different runs.

    Returns:
        Dictionary with cross-run summary metrics.
    """
    if not pipeline_results:
        return {
            "num_runs": 0,
            "runs": [],
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }

    runs: list[dict[str, Any]] = []
    total_duration = 0.0
    success_count = 0

    for result in pipeline_results:
        pipeline_name = getattr(result, "pipeline_name", "unknown")
        success = getattr(result, "success", False)
        duration = getattr(result, "total_duration", 0.0)
        output_dir = getattr(result, "output_dir", Path("."))
        steps = getattr(result, "steps", [])

        run_info: dict[str, Any] = {
            "pipeline_name": pipeline_name,
            "success": success,
            "total_duration": duration,
            "output_dir": str(output_dir),
            "steps_total": len(steps),
            "steps_completed": sum(1 for s in steps if s.status == "completed"),
            "steps_failed": sum(1 for s in steps if s.status == "failed"),
            "steps_skipped": sum(1 for s in steps if s.status == "skipped"),
        }

        runs.append(run_info)
        total_duration += duration
        if success:
            success_count += 1

    summary: dict[str, Any] = {
        "num_runs": len(pipeline_results),
        "success_count": success_count,
        "failure_count": len(pipeline_results) - success_count,
        "total_duration": total_duration,
        "mean_duration": total_duration / len(pipeline_results),
        "runs": runs,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    logger.info(
        "Generated run summary: %d runs, %d/%d successful, total %.1fs",
        len(pipeline_results),
        success_count,
        len(pipeline_results),
        total_duration,
    )
    return summary


# --- Internal report formatters ---


def _write_json_report(report_dict: dict[str, Any], output_path: Path) -> None:
    """Write report as formatted JSON."""

    def _default(obj: Any) -> Any:
        if isinstance(obj, Path):
            return str(obj)
        if hasattr(obj, "__dataclass_fields__"):
            return asdict(obj)
        return str(obj)

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(report_dict, f, indent=2, default=_default, ensure_ascii=False)


def _write_text_report(report_dict: dict[str, Any], output_path: Path) -> None:
    """Write report as human-readable plain text."""
    lines: list[str] = []

    # Header
    pipeline = report_dict.get("pipeline", report_dict.get("sample_name", "Report"))
    lines.append("=" * 60)
    lines.append(f"  Long-Read Analysis Report: {pipeline}")
    lines.append("=" * 60)
    lines.append(f"  Generated: {report_dict.get('timestamp', 'N/A')}")
    lines.append("")

    # Flatten top-level scalar values
    for key, value in report_dict.items():
        if isinstance(value, (str, int, float, bool)):
            label = key.replace("_", " ").title()
            if isinstance(value, float):
                lines.append(f"  {label}: {value:.4f}")
            else:
                lines.append(f"  {label}: {value}")

    lines.append("")

    # Sub-sections for nested dicts
    for key, value in report_dict.items():
        if isinstance(value, dict):
            lines.append(f"  --- {key.replace('_', ' ').title()} ---")
            for sub_key, sub_value in value.items():
                if isinstance(sub_value, (str, int, float, bool)):
                    label = sub_key.replace("_", " ").title()
                    if isinstance(sub_value, float):
                        lines.append(f"    {label}: {sub_value:.4f}")
                    else:
                        lines.append(f"    {label}: {sub_value}")
            lines.append("")

    # Step summaries if present
    steps_data = report_dict.get("steps", [])
    if steps_data:
        lines.append("  --- Pipeline Steps ---")
        for step in steps_data:
            if isinstance(step, dict):
                name = step.get("name", "?")
                status = step.get("status", "?")
                duration = step.get("duration_seconds", 0.0)
                error = step.get("error", "")
                status_symbol = "OK" if status == "completed" else status.upper()
                lines.append(f"    [{status_symbol}] {name} ({duration:.2f}s)")
                if error:
                    lines.append(f"           Error: {error}")
        lines.append("")

    lines.append("=" * 60)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def _write_html_report(report_dict: dict[str, Any], output_path: Path) -> None:
    """Write report as a self-contained HTML page."""
    pipeline = report_dict.get("pipeline", report_dict.get("sample_name", "Report"))
    timestamp = report_dict.get("timestamp", "N/A")

    # Build HTML content
    html_parts: list[str] = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        f"<title>Long-Read Report: {pipeline}</title>",
        "<style>",
        "body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; ",
        "  max-width: 900px; margin: 0 auto; padding: 20px; color: #333; }",
        "h1 { color: #1a237e; border-bottom: 2px solid #1a237e; padding-bottom: 10px; }",
        "h2 { color: #283593; margin-top: 30px; }",
        "table { border-collapse: collapse; width: 100%; margin: 10px 0; }",
        "th, td { padding: 8px 12px; text-align: left; border: 1px solid #ddd; }",
        "th { background-color: #e8eaf6; font-weight: 600; }",
        "tr:nth-child(even) { background-color: #f5f5f5; }",
        ".success { color: #2e7d32; font-weight: bold; }",
        ".failed { color: #c62828; font-weight: bold; }",
        ".skipped { color: #f57f17; }",
        ".metric-value { font-family: monospace; }",
        ".timestamp { color: #757575; font-size: 0.9em; }",
        "</style>",
        "</head>",
        "<body>",
        f"<h1>Long-Read Analysis Report: {pipeline}</h1>",
        f'<p class="timestamp">Generated: {timestamp}</p>',
    ]

    # Summary metrics table
    html_parts.append("<h2>Summary Metrics</h2>")
    html_parts.append("<table>")
    html_parts.append("<tr><th>Metric</th><th>Value</th></tr>")

    for key, value in report_dict.items():
        if isinstance(value, (str, int, float, bool)) and key not in ("pipeline", "timestamp"):
            label = key.replace("_", " ").title()
            if isinstance(value, float):
                formatted = f"{value:.4f}"
            elif isinstance(value, bool):
                formatted = "Yes" if value else "No"
            else:
                formatted = str(value)
            html_parts.append(f'<tr><td>{label}</td><td class="metric-value">{formatted}</td></tr>')

    html_parts.append("</table>")

    # Sub-sections
    for key, value in report_dict.items():
        if isinstance(value, dict) and key not in ("steps",):
            section_title = key.replace("_", " ").title()
            html_parts.append(f"<h2>{section_title}</h2>")
            html_parts.append("<table>")
            html_parts.append("<tr><th>Metric</th><th>Value</th></tr>")

            for sub_key, sub_value in value.items():
                if isinstance(sub_value, (str, int, float, bool)):
                    label = sub_key.replace("_", " ").title()
                    if isinstance(sub_value, float):
                        formatted = f"{sub_value:.4f}"
                    else:
                        formatted = str(sub_value)
                    html_parts.append(f'<tr><td>{label}</td><td class="metric-value">{formatted}</td></tr>')

            html_parts.append("</table>")

    # Step details
    steps_data = report_dict.get("steps", [])
    if steps_data:
        html_parts.append("<h2>Pipeline Steps</h2>")
        html_parts.append("<table>")
        html_parts.append("<tr><th>Step</th><th>Status</th><th>Duration (s)</th><th>Error</th></tr>")

        for step in steps_data:
            if isinstance(step, dict):
                name = step.get("name", "?")
                status = step.get("status", "?")
                duration = step.get("duration_seconds", 0.0)
                error = step.get("error", "")

                status_class = "success" if status == "completed" else ("failed" if status == "failed" else "skipped")
                html_parts.append(
                    f'<tr><td>{name}</td><td class="{status_class}">{status}</td>'
                    f'<td class="metric-value">{duration:.2f}</td><td>{error}</td></tr>'
                )

        html_parts.append("</table>")

    html_parts.extend(
        [
            "</body>",
            "</html>",
        ]
    )

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(html_parts))
