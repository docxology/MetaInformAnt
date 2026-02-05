"""Run summary generation for long-read sequencing analysis.

Aggregates results from multiple analysis steps into comprehensive
run summaries with statistics, comparisons, and exportable reports.
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "LR_"


@dataclass
class RunSummary:
    """Comprehensive summary of a long-read sequencing analysis run.

    Attributes:
        run_id: Unique identifier for this run.
        sample_name: Sample/experiment name.
        platform: Sequencing platform (ont, pacbio).
        timestamp: ISO format timestamp of when summary was generated.
        input_stats: Raw input data statistics.
        qc_stats: Quality control statistics.
        assembly_stats: Assembly statistics (if assembly was run).
        methylation_stats: Methylation analysis statistics (if run).
        sv_stats: Structural variant statistics (if run).
        phasing_stats: Phasing statistics (if run).
        pipeline_steps: List of pipeline steps executed.
        warnings: Any warnings generated during the run.
    """

    run_id: str = ""
    sample_name: str = ""
    platform: str = "ont"
    timestamp: str = ""
    input_stats: dict[str, Any] = field(default_factory=dict)
    qc_stats: dict[str, Any] = field(default_factory=dict)
    assembly_stats: dict[str, Any] = field(default_factory=dict)
    methylation_stats: dict[str, Any] = field(default_factory=dict)
    sv_stats: dict[str, Any] = field(default_factory=dict)
    phasing_stats: dict[str, Any] = field(default_factory=dict)
    pipeline_steps: list[dict[str, Any]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def generate_qc_summary(
    read_lengths: list[int],
    quality_scores: list[float] | None = None,
    filtered_count: int | None = None,
    adapter_trimmed: int = 0,
    chimeras_split: int = 0,
) -> dict[str, Any]:
    """Generate a QC summary from read data.

    Computes comprehensive quality statistics using the longread quality
    metrics functions.

    Args:
        read_lengths: List of read lengths (post-filtering if applicable).
        quality_scores: Optional list of per-read mean quality scores.
        filtered_count: Number of reads removed by filtering (None = not filtered).
        adapter_trimmed: Number of reads with adapters trimmed.
        chimeras_split: Number of chimeric reads split.

    Returns:
        Dict with QC summary statistics.
    """
    from metainformant.longread.quality.metrics import (
        calculate_n50,
        calculate_nx,
        read_length_stats,
    )

    if not read_lengths:
        return {
            "total_reads": 0,
            "total_bases": 0,
            "n50": 0,
            "mean_length": 0.0,
            "status": "no_data",
        }

    stats = read_length_stats(read_lengths)

    summary: dict[str, Any] = {
        "total_reads": stats.count,
        "total_bases": stats.total_bases,
        "mean_length": round(stats.mean_length, 1),
        "median_length": round(stats.median_length, 1),
        "min_length": stats.min_length,
        "max_length": stats.max_length,
        "n50": stats.n50,
        "n90": stats.n90,
        "l50": stats.l50,
        "std_dev": round(stats.std_dev, 1),
    }

    if quality_scores:
        mean_q = sum(quality_scores) / len(quality_scores)
        summary["mean_quality"] = round(mean_q, 2)
        summary["q7_fraction"] = round(sum(1 for q in quality_scores if q >= 7) / len(quality_scores), 4)
        summary["q10_fraction"] = round(sum(1 for q in quality_scores if q >= 10) / len(quality_scores), 4)
        summary["q20_fraction"] = round(sum(1 for q in quality_scores if q >= 20) / len(quality_scores), 4)

    if filtered_count is not None:
        original_count = stats.count + filtered_count
        summary["reads_before_filter"] = original_count
        summary["reads_removed"] = filtered_count
        summary["filter_pass_rate"] = round(stats.count / original_count, 4) if original_count > 0 else 0.0

    summary["adapter_trimmed"] = adapter_trimmed
    summary["chimeras_split"] = chimeras_split
    summary["status"] = "complete"

    return summary


def generate_assembly_summary(
    contig_lengths: list[int],
    num_reads_used: int = 0,
    polish_iterations: int = 0,
    coverage: float = 0.0,
) -> dict[str, Any]:
    """Generate assembly summary statistics.

    Args:
        contig_lengths: List of assembled contig lengths.
        num_reads_used: Number of input reads.
        polish_iterations: Number of polishing rounds performed.
        coverage: Mean coverage depth.

    Returns:
        Dict with assembly summary statistics.
    """
    from metainformant.longread.quality.metrics import calculate_n50, calculate_nx

    if not contig_lengths:
        return {"total_contigs": 0, "total_bases": 0, "status": "no_assembly"}

    sorted_lengths = sorted(contig_lengths, reverse=True)
    total_bases = sum(sorted_lengths)

    return {
        "total_contigs": len(sorted_lengths),
        "total_bases": total_bases,
        "largest_contig": sorted_lengths[0],
        "smallest_contig": sorted_lengths[-1],
        "mean_length": round(total_bases / len(sorted_lengths), 1),
        "n50": calculate_n50(sorted_lengths),
        "n90": calculate_nx(sorted_lengths, x=90),
        "gc_content": None,  # Computed separately if sequences available
        "reads_used": num_reads_used,
        "polish_iterations": polish_iterations,
        "mean_coverage": round(coverage, 1),
        "status": "complete",
    }


def generate_methylation_summary(
    total_sites: int,
    methylated_sites: int,
    modification_type: str = "5mC",
    mean_coverage: float = 0.0,
    regions_analyzed: int = 0,
    differential_sites: int = 0,
) -> dict[str, Any]:
    """Generate methylation analysis summary.

    Args:
        total_sites: Total CpG/modification sites analyzed.
        methylated_sites: Number of sites called as methylated.
        modification_type: Type of modification analyzed.
        mean_coverage: Mean coverage at analyzed sites.
        regions_analyzed: Number of genomic regions analyzed.
        differential_sites: Number of differentially methylated sites.

    Returns:
        Dict with methylation summary statistics.
    """
    methylation_rate = methylated_sites / total_sites if total_sites > 0 else 0.0

    return {
        "modification_type": modification_type,
        "total_sites": total_sites,
        "methylated_sites": methylated_sites,
        "unmethylated_sites": total_sites - methylated_sites,
        "global_methylation_rate": round(methylation_rate, 4),
        "mean_coverage": round(mean_coverage, 1),
        "regions_analyzed": regions_analyzed,
        "differential_sites": differential_sites,
        "status": "complete",
    }


def generate_sv_summary(
    variants: list[dict[str, Any]],
) -> dict[str, Any]:
    """Generate structural variant calling summary.

    Args:
        variants: List of SV dicts with at minimum 'sv_type' and 'size' keys.

    Returns:
        Dict with SV calling summary.
    """
    if not variants:
        return {"total_variants": 0, "status": "no_variants"}

    type_counts: dict[str, int] = {}
    sizes: list[int] = []

    for v in variants:
        sv_type = v.get("sv_type", "UNKNOWN")
        type_counts[sv_type] = type_counts.get(sv_type, 0) + 1
        size = v.get("size", 0)
        if size > 0:
            sizes.append(size)

    size_stats: dict[str, Any] = {}
    if sizes:
        sorted_sizes = sorted(sizes)
        size_stats = {
            "min_size": sorted_sizes[0],
            "max_size": sorted_sizes[-1],
            "median_size": sorted_sizes[len(sorted_sizes) // 2],
            "mean_size": round(sum(sorted_sizes) / len(sorted_sizes), 1),
        }

    return {
        "total_variants": len(variants),
        "by_type": type_counts,
        "size_distribution": size_stats,
        "phased_count": sum(1 for v in variants if v.get("haplotype", 0) != 0),
        "status": "complete",
    }


def build_run_summary(
    sample_name: str,
    platform: str = "ont",
    qc_stats: dict[str, Any] | None = None,
    assembly_stats: dict[str, Any] | None = None,
    methylation_stats: dict[str, Any] | None = None,
    sv_stats: dict[str, Any] | None = None,
    phasing_stats: dict[str, Any] | None = None,
) -> RunSummary:
    """Build a comprehensive run summary from component statistics.

    Args:
        sample_name: Name of the sample/experiment.
        platform: Sequencing platform ('ont' or 'pacbio').
        qc_stats: QC statistics dict (from generate_qc_summary).
        assembly_stats: Assembly statistics dict.
        methylation_stats: Methylation statistics dict.
        sv_stats: SV statistics dict.
        phasing_stats: Phasing statistics dict.

    Returns:
        RunSummary dataclass with all statistics.
    """
    import datetime

    summary = RunSummary(
        run_id=f"lr_{int(time.time())}_{sample_name}",
        sample_name=sample_name,
        platform=platform,
        timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    )

    if qc_stats:
        summary.qc_stats = qc_stats
    if assembly_stats:
        summary.assembly_stats = assembly_stats
    if methylation_stats:
        summary.methylation_stats = methylation_stats
    if sv_stats:
        summary.sv_stats = sv_stats
    if phasing_stats:
        summary.phasing_stats = phasing_stats

    logger.info(f"Built run summary for sample '{sample_name}' ({platform})")
    return summary


def export_run_summary(
    summary: RunSummary,
    output_path: Path | str,
    format: str = "json",
) -> Path:
    """Export run summary to file.

    Args:
        summary: RunSummary to export.
        output_path: Destination file path.
        format: Output format ('json', 'text').

    Returns:
        Path to the exported file.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if format == "json":
        data = asdict(summary)
        path.write_text(json.dumps(data, indent=2, default=str))

    elif format == "text":
        lines = [
            f"{'=' * 60}",
            f"LONG-READ SEQUENCING RUN SUMMARY",
            f"{'=' * 60}",
            f"Run ID:    {summary.run_id}",
            f"Sample:    {summary.sample_name}",
            f"Platform:  {summary.platform}",
            f"Timestamp: {summary.timestamp}",
            "",
        ]

        if summary.qc_stats:
            lines.extend(
                [
                    f"--- Quality Control ---",
                    f"  Total reads:   {summary.qc_stats.get('total_reads', 'N/A'):,}",
                    f"  Total bases:   {summary.qc_stats.get('total_bases', 'N/A'):,}",
                    f"  N50:           {summary.qc_stats.get('n50', 'N/A'):,}",
                    f"  Mean length:   {summary.qc_stats.get('mean_length', 'N/A'):,.1f}",
                    f"  Mean quality:  {summary.qc_stats.get('mean_quality', 'N/A')}",
                    "",
                ]
            )

        if summary.assembly_stats:
            lines.extend(
                [
                    f"--- Assembly ---",
                    f"  Contigs:       {summary.assembly_stats.get('total_contigs', 'N/A'):,}",
                    f"  Total bases:   {summary.assembly_stats.get('total_bases', 'N/A'):,}",
                    f"  N50:           {summary.assembly_stats.get('n50', 'N/A'):,}",
                    f"  Largest:       {summary.assembly_stats.get('largest_contig', 'N/A'):,}",
                    "",
                ]
            )

        if summary.methylation_stats:
            lines.extend(
                [
                    f"--- Methylation ---",
                    f"  Mod type:      {summary.methylation_stats.get('modification_type', 'N/A')}",
                    f"  Total sites:   {summary.methylation_stats.get('total_sites', 'N/A'):,}",
                    f"  Methylated:    {summary.methylation_stats.get('methylated_sites', 'N/A'):,}",
                    f"  Global rate:   {summary.methylation_stats.get('global_methylation_rate', 'N/A')}",
                    "",
                ]
            )

        if summary.sv_stats:
            lines.extend(
                [
                    f"--- Structural Variants ---",
                    f"  Total SVs:     {summary.sv_stats.get('total_variants', 'N/A')}",
                    f"  By type:       {summary.sv_stats.get('by_type', {})}",
                    "",
                ]
            )

        if summary.warnings:
            lines.extend(
                [
                    f"--- Warnings ---",
                    *[f"  - {w}" for w in summary.warnings],
                    "",
                ]
            )

        lines.append(f"{'=' * 60}")
        path.write_text("\n".join(lines))

    else:
        raise ValueError(f"Unsupported format: {format}. Use 'json' or 'text'.")

    logger.info(f"Exported run summary to {path} ({format})")
    return path


def compare_run_summaries(
    summaries: list[RunSummary],
) -> dict[str, Any]:
    """Compare statistics across multiple run summaries.

    Useful for comparing samples or tracking quality across runs.

    Args:
        summaries: List of RunSummary objects to compare.

    Returns:
        Dict with per-sample comparisons and aggregate statistics.
    """
    if not summaries:
        return {"samples": 0}

    comparison: dict[str, Any] = {
        "samples": len(summaries),
        "platforms": list({s.platform for s in summaries}),
        "per_sample": {},
    }

    n50_values: list[int] = []
    total_bases_values: list[int] = []
    mean_quality_values: list[float] = []

    for s in summaries:
        sample_data: dict[str, Any] = {
            "platform": s.platform,
        }

        if s.qc_stats:
            n50 = s.qc_stats.get("n50", 0)
            total_bases = s.qc_stats.get("total_bases", 0)
            mean_q = s.qc_stats.get("mean_quality")

            sample_data["n50"] = n50
            sample_data["total_bases"] = total_bases
            sample_data["total_reads"] = s.qc_stats.get("total_reads", 0)
            sample_data["mean_quality"] = mean_q

            if n50 > 0:
                n50_values.append(n50)
            if total_bases > 0:
                total_bases_values.append(total_bases)
            if mean_q is not None:
                mean_quality_values.append(mean_q)

        comparison["per_sample"][s.sample_name] = sample_data

    # Aggregate statistics
    if n50_values:
        comparison["aggregate_n50"] = {
            "min": min(n50_values),
            "max": max(n50_values),
            "mean": round(sum(n50_values) / len(n50_values)),
        }
    if total_bases_values:
        comparison["aggregate_bases"] = {
            "min": min(total_bases_values),
            "max": max(total_bases_values),
            "total": sum(total_bases_values),
        }
    if mean_quality_values:
        comparison["aggregate_quality"] = {
            "min": round(min(mean_quality_values), 2),
            "max": round(max(mean_quality_values), 2),
            "mean": round(sum(mean_quality_values) / len(mean_quality_values), 2),
        }

    return comparison


__all__ = [
    "RunSummary",
    "generate_qc_summary",
    "generate_assembly_summary",
    "generate_methylation_summary",
    "generate_sv_summary",
    "build_run_summary",
    "export_run_summary",
    "compare_run_summaries",
]
