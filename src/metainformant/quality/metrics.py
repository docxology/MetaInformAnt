"""Quality metrics and statistical analysis for biological data."""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np


def calculate_quality_metrics(quality_scores: List[List[int]]) -> Dict[str, float]:
    """
    Calculate comprehensive quality metrics from per-position quality scores.

    Args:
        quality_scores: List of lists containing quality scores per position

    Returns:
        Dict containing various quality metrics
    """
    if not quality_scores:
        return {}

    # Convert to numpy array for efficient computation
    scores_array = np.array(quality_scores)

    metrics = {}

    # Basic statistics
    metrics["mean_quality"] = float(np.mean(scores_array))
    metrics["median_quality"] = float(np.median(scores_array))
    metrics["std_quality"] = float(np.std(scores_array))
    metrics["min_quality"] = float(np.min(scores_array))
    metrics["max_quality"] = float(np.max(scores_array))

    # Quality distribution
    quality_range = np.arange(0, 43)  # Phred scores 0-42
    hist, _ = np.histogram(scores_array, bins=quality_range)

    # Q20, Q30 percentages
    metrics["pct_q20"] = float(np.sum(scores_array >= 20) / scores_array.size * 100)
    metrics["pct_q30"] = float(np.sum(scores_array >= 30) / scores_array.size * 100)

    # Quality score entropy
    hist = hist / hist.sum() if hist.sum() > 0 else hist
    entropy = -np.sum(hist * np.log2(hist + 1e-10))
    metrics["quality_entropy"] = float(entropy)

    return metrics


def calculate_gc_metrics(gc_content: List[float]) -> Dict[str, float]:
    """
    Calculate GC content related quality metrics.

    Args:
        gc_content: List of GC content percentages per sequence

    Returns:
        Dict containing GC-related quality metrics
    """
    if not gc_content:
        return {}

    gc_array = np.array(gc_content)

    metrics = {}

    # Basic statistics
    metrics["mean_gc"] = float(np.mean(gc_array))
    metrics["median_gc"] = float(np.median(gc_array))
    metrics["std_gc"] = float(np.std(gc_array))

    # GC content distribution
    metrics["gc_25th_percentile"] = float(np.percentile(gc_array, 25))
    metrics["gc_75th_percentile"] = float(np.percentile(gc_array, 75))

    # GC bias indicators
    expected_gc = 50.0  # Expected GC content for most organisms
    metrics["gc_bias"] = abs(metrics["mean_gc"] - expected_gc)

    # GC content uniformity
    cv = metrics["std_gc"] / metrics["mean_gc"] if metrics["mean_gc"] > 0 else 0
    metrics["gc_coefficient_variation"] = float(cv)

    return metrics


def calculate_length_metrics(sequence_lengths: List[int]) -> Dict[str, float]:
    """
    Calculate sequence length related quality metrics.

    Args:
        sequence_lengths: List of sequence lengths

    Returns:
        Dict containing length-related quality metrics
    """
    if not sequence_lengths:
        return {}

    lengths_array = np.array(sequence_lengths)

    metrics = {}

    # Basic statistics
    metrics["mean_length"] = float(np.mean(lengths_array))
    metrics["median_length"] = float(np.median(lengths_array))
    metrics["std_length"] = float(np.std(lengths_array))
    metrics["min_length"] = float(np.min(lengths_array))
    metrics["max_length"] = float(np.max(lengths_array))

    # Length distribution
    metrics["length_25th_percentile"] = float(np.percentile(lengths_array, 25))
    metrics["length_75th_percentile"] = float(np.percentile(lengths_array, 75))

    # Length uniformity
    cv = metrics["std_length"] / metrics["mean_length"] if metrics["mean_length"] > 0 else 0
    metrics["length_coefficient_variation"] = float(cv)

    # Length quality indicators
    metrics["pct_short_reads"] = float(np.sum(lengths_array < 50) / len(lengths_array) * 100)
    metrics["pct_long_reads"] = float(np.sum(lengths_array > 1000) / len(lengths_array) * 100)

    return metrics


def calculate_duplication_metrics(duplication_levels: Dict[int, int]) -> Dict[str, float]:
    """
    Calculate duplication-related quality metrics.

    Args:
        duplication_levels: Dict mapping duplication level to count

    Returns:
        Dict containing duplication quality metrics
    """
    if not duplication_levels:
        return {}

    total_reads = sum(duplication_levels.values())
    if total_reads == 0:
        return {}

    metrics = {}

    # Duplication statistics
    unique_reads = duplication_levels.get(1, 0)
    duplicated_reads = total_reads - unique_reads

    metrics["duplication_rate"] = float(duplicated_reads / total_reads * 100)
    metrics["unique_rate"] = float(unique_reads / total_reads * 100)

    # Duplication diversity (entropy)
    probs = np.array(list(duplication_levels.values())) / total_reads
    entropy = -np.sum(probs * np.log2(probs + 1e-10))
    metrics["duplication_entropy"] = float(entropy)

    # Average duplication level
    weighted_sum = sum(level * count for level, count in duplication_levels.items())
    metrics["mean_duplication_level"] = float(weighted_sum / total_reads)

    return metrics


def calculate_complexity_metrics(sequences: List[str]) -> Dict[str, float]:
    """
    Calculate sequence complexity metrics.

    Args:
        sequences: List of sequences to analyze

    Returns:
        Dict containing complexity metrics
    """
    if not sequences:
        return {}

    metrics = {}

    # Calculate complexity for each sequence
    complexities = []
    for seq in sequences:
        if len(seq) == 0:
            continue

        # Simple complexity: unique characters / total characters
        unique_chars = len(set(seq))
        complexity = unique_chars / len(seq)
        complexities.append(complexity)

    if not complexities:
        return {}

    complexity_array = np.array(complexities)

    metrics["mean_complexity"] = float(np.mean(complexity_array))
    metrics["median_complexity"] = float(np.median(complexity_array))
    metrics["std_complexity"] = float(np.std(complexity_array))

    # Complexity distribution
    metrics["low_complexity_rate"] = float(np.sum(complexity_array < 0.5) / len(complexities) * 100)

    return metrics


def calculate_coverage_metrics(
    coverage_values: List[float],
    target_coverage: float = 30.0,
) -> Dict[str, float]:
    """
    Calculate coverage-related quality metrics.

    Args:
        coverage_values: List of coverage values
        target_coverage: Expected target coverage

    Returns:
        Dict containing coverage quality metrics
    """
    if not coverage_values:
        return {}

    coverage_array = np.array(coverage_values)

    metrics = {}

    # Basic statistics
    metrics["mean_coverage"] = float(np.mean(coverage_array))
    metrics["median_coverage"] = float(np.median(coverage_array))
    metrics["std_coverage"] = float(np.std(coverage_array))

    # Coverage uniformity
    cv = metrics["std_coverage"] / metrics["mean_coverage"] if metrics["mean_coverage"] > 0 else 0
    metrics["coverage_coefficient_variation"] = float(cv)

    # Coverage bias
    metrics["coverage_bias"] = abs(metrics["mean_coverage"] - target_coverage)

    # Coverage distribution
    metrics["pct_low_coverage"] = float(np.sum(coverage_array < target_coverage * 0.5) / len(coverage_array) * 100)
    metrics["pct_high_coverage"] = float(np.sum(coverage_array > target_coverage * 1.5) / len(coverage_array) * 100)

    return metrics


def generate_quality_report(
    quality_data: Dict[str, Dict[str, float]],
    sample_name: str = "Unknown",
) -> str:
    """
    Generate a comprehensive quality report from quality metrics.

    Args:
        quality_data: Dict containing quality metrics from various analyses
        sample_name: Name of the sample being analyzed

    Returns:
        Formatted quality report
    """
    report_lines = []
    report_lines.append("METAINFORMANT Quality Control Report")
    report_lines.append("=" * 50)
    report_lines.append(f"Sample: {sample_name}")
    report_lines.append("")

    # Quality metrics section
    if "quality" in quality_data:
        qm = quality_data["quality"]
        report_lines.append("Quality Metrics:")
        report_lines.append(f"  Mean Quality: {qm.get('mean_quality', 0):.2f}")
        report_lines.append(f"  Median Quality: {qm.get('median_quality', 0):.2f}")
        report_lines.append(f"  Q20 Rate: {qm.get('pct_q20', 0):.1f}%")
        report_lines.append(f"  Q30 Rate: {qm.get('pct_q30', 0):.1f}%")
        report_lines.append("")

    # GC content metrics
    if "gc_content" in quality_data:
        gcm = quality_data["gc_content"]
        report_lines.append("GC Content Metrics:")
        report_lines.append(f"  Mean GC: {gcm.get('mean_gc', 0):.1f}%")
        report_lines.append(f"  GC Coefficient of Variation: {gcm.get('gc_coefficient_variation', 0):.3f}")
        report_lines.append("")

    # Length metrics
    if "length" in quality_data:
        lm = quality_data["length"]
        report_lines.append("Length Metrics:")
        report_lines.append(f"  Mean Length: {lm.get('mean_length', 0):.0f} bp")
        report_lines.append(f"  Length CV: {lm.get('length_coefficient_variation', 0):.3f}")
        report_lines.append("")

    # Duplication metrics
    if "duplication" in quality_data:
        dm = quality_data["duplication"]
        report_lines.append("Duplication Metrics:")
        report_lines.append(f"  Duplication Rate: {dm.get('duplication_rate', 0):.1f}%")
        report_lines.append(f"  Mean Duplication Level: {dm.get('mean_duplication_level', 0):.2f}")
        report_lines.append("")

    # Complexity metrics
    if "complexity" in quality_data:
        cm = quality_data["complexity"]
        report_lines.append("Complexity Metrics:")
        report_lines.append(f"  Mean Complexity: {cm.get('mean_complexity', 0):.3f}")
        report_lines.append(f"  Low Complexity Rate: {cm.get('low_complexity_rate', 0):.1f}%")
        report_lines.append("")

    # Coverage metrics
    if "coverage" in quality_data:
        covm = quality_data["coverage"]
        report_lines.append("Coverage Metrics:")
        report_lines.append(f"  Mean Coverage: {covm.get('mean_coverage', 0):.1f}x")
        report_lines.append(f"  Coverage CV: {covm.get('coverage_coefficient_variation', 0):.3f}")
        report_lines.append("")

    # Overall assessment
    report_lines.append("Quality Assessment:")
    overall_score = _calculate_overall_quality_score(quality_data)
    report_lines.append(f"  Overall Quality Score: {overall_score:.2f}/10.0")

    if overall_score >= 8.0:
        report_lines.append("  Assessment: EXCELLENT")
    elif overall_score >= 6.0:
        report_lines.append("  Assessment: GOOD")
    elif overall_score >= 4.0:
        report_lines.append("  Assessment: ACCEPTABLE")
    else:
        report_lines.append("  Assessment: POOR")

    return "\n".join(report_lines)


def _calculate_overall_quality_score(quality_data: Dict[str, Dict[str, float]]) -> float:
    """Calculate an overall quality score from 0-10."""
    score = 5.0  # Base score

    # Quality score contribution
    if "quality" in quality_data:
        qm = quality_data["quality"]
        mean_q = qm.get("mean_quality", 30)
        if mean_q >= 35:
            score += 2.0
        elif mean_q >= 30:
            score += 1.0
        elif mean_q < 25:
            score -= 1.0

    # GC content contribution
    if "gc_content" in quality_data:
        gcm = quality_data["gc_content"]
        gc_bias = gcm.get("gc_bias", 0)
        if gc_bias <= 5:
            score += 1.0
        elif gc_bias > 15:
            score -= 1.0

    # Duplication contribution
    if "duplication" in quality_data:
        dm = quality_data["duplication"]
        dup_rate = dm.get("duplication_rate", 50)
        if dup_rate <= 20:
            score += 1.0
        elif dup_rate > 60:
            score -= 1.0

    # Complexity contribution
    if "complexity" in quality_data:
        cm = quality_data["complexity"]
        low_complexity = cm.get("low_complexity_rate", 0)
        if low_complexity <= 5:
            score += 1.0
        elif low_complexity > 20:
            score -= 1.0

    return min(10.0, max(0.0, score))
