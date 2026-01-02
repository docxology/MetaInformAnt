"""Quality control metrics and scoring.

This module provides various quality control metrics for biological data,
including composite quality scores, data integrity checks, and quality
assessment algorithms.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
import numpy as np
import pandas as pd
from scipy import stats

from metainformant.core import logging, errors, validation

logger = logging.get_logger(__name__)


def calculate_quality_score(data: Dict[str, Any], data_type: str = "fastq") -> Dict[str, Any]:
    """Calculate an overall quality score for biological data.

    Args:
        data: Quality metrics data (e.g., from analyze_fastq_quality)
        data_type: Type of data ("fastq", "vcf", "bam", etc.)

    Returns:
        Dictionary with quality score and breakdown
    """
    if data_type == "fastq":
        return _calculate_fastq_quality_score(data)
    elif data_type == "vcf":
        return _calculate_vcf_quality_score(data)
    elif data_type == "bam":
        return _calculate_bam_quality_score(data)
    else:
        raise ValueError(f"Unsupported data type: {data_type}")


def _calculate_fastq_quality_score(data: Dict[str, Any]) -> Dict[str, Any]:
    """Calculate quality score for FASTQ data."""
    score_components = {}

    # Basic statistics component (40% weight)
    if "basic_statistics" in data:
        stats = data["basic_statistics"]
        mean_qual = stats.get("mean_quality", 0)
        # Quality score: 0-40 points based on mean quality
        qual_score = min(40, max(0, (mean_qual - 20) * 2))  # 20 = 0 points, 40 = 40 points
        score_components["basic_quality"] = {
            "score": qual_score,
            "weight": 0.4,
            "details": f"Mean quality: {mean_qual:.1f}",
        }

    # Per-base quality component (30% weight)
    if "per_base_quality" in data:
        pbq = data["per_base_quality"]
        positions = pbq.get("positions", [])
        if positions:
            # Check if quality drops significantly at end
            final_qualities = [p["mean"] for p in positions[-5:]]  # Last 5 positions
            initial_qualities = [p["mean"] for p in positions[:5]]  # First 5 positions

            if final_qualities and initial_qualities:
                final_mean = sum(final_qualities) / len(final_qualities)
                initial_mean = sum(initial_qualities) / len(initial_qualities)
                degradation = initial_mean - final_mean

                # Score: 0-30 points, penalty for >5 point degradation
                pbq_score = max(0, 30 - degradation * 2)
                score_components["per_base_quality"] = {
                    "score": pbq_score,
                    "weight": 0.3,
                    "details": f"Quality degradation: {degradation:.1f} points",
                }

    # GC content component (15% weight)
    if "gc_content_distribution" in data:
        gc_data = data["gc_content_distribution"]
        bins = gc_data.get("bins", [])
        if bins:
            # Check for abnormal GC distribution
            total_reads = sum(b["count"] for b in bins)
            gc_40_60 = sum(b["count"] for b in bins if 40 <= b["bin_start"] < 60)
            normal_gc_fraction = gc_40_60 / total_reads if total_reads > 0 else 0

            # Score: higher for more normal GC distribution
            gc_score = normal_gc_fraction * 15
            score_components["gc_content"] = {
                "score": gc_score,
                "weight": 0.15,
                "details": f"Normal GC fraction: {normal_gc_fraction:.2%}",
            }

    # Adapter content component (15% weight)
    if "adapter_content" in data:
        adapter_data = data["adapter_content"]
        adapters = adapter_data.get("adapters", {})
        max_adapter_percent = 0

        for adapter_info in adapters.values():
            positions = adapter_info.get("positions", [])
            for pos in positions:
                max_adapter_percent = max(max_adapter_percent, pos.get("percentage", 0))

        # Score: penalty for adapter contamination
        adapter_score = max(0, 15 - max_adapter_percent * 0.3)
        score_components["adapter_content"] = {
            "score": adapter_score,
            "weight": 0.15,
            "details": f"Max adapter content: {max_adapter_percent:.1f}%",
        }

    # Calculate total score
    total_score = 0
    max_score = 0

    for component in score_components.values():
        total_score += component["score"]
        max_score += component["score"] / component["weight"] * component["weight"]

    overall_score = (total_score / max_score * 100) if max_score > 0 else 0

    return {
        "overall_score": overall_score,
        "total_score": total_score,
        "max_possible_score": max_score,
        "components": score_components,
        "grade": _score_to_grade(overall_score),
    }


def _calculate_vcf_quality_score(data: Dict[str, Any]) -> Dict[str, Any]:
    """Calculate quality score for VCF data."""
    # Placeholder for VCF quality scoring
    return {
        "overall_score": 85.0,  # Placeholder
        "components": {},
        "grade": "B",
    }


def _calculate_bam_quality_score(data: Dict[str, Any]) -> Dict[str, Any]:
    """Calculate quality score for BAM data."""
    # Placeholder for BAM quality scoring
    return {
        "overall_score": 82.0,  # Placeholder
        "components": {},
        "grade": "B",
    }


def _score_to_grade(score: float) -> str:
    """Convert numeric score to letter grade."""
    if score >= 90:
        return "A"
    elif score >= 80:
        return "B"
    elif score >= 70:
        return "C"
    elif score >= 60:
        return "D"
    else:
        return "F"


def detect_outliers(data: List[float], method: str = "iqr", threshold: float = 1.5) -> Dict[str, Any]:
    """Detect outliers in quality metrics data.

    Args:
        data: List of numeric values
        method: Outlier detection method ("iqr", "zscore", "modified_zscore")
        threshold: Threshold for outlier detection

    Returns:
        Dictionary with outlier detection results
    """
    if not data:
        return {"outliers": [], "outlier_indices": []}

    data_array = np.array(data)

    if method == "iqr":
        q1 = np.percentile(data_array, 25)
        q3 = np.percentile(data_array, 75)
        iqr = q3 - q1
        lower_bound = q1 - threshold * iqr
        upper_bound = q3 + threshold * iqr

        outliers = []
        outlier_indices = []

        for i, value in enumerate(data):
            if value < lower_bound or value > upper_bound:
                outliers.append(value)
                outlier_indices.append(i)

    elif method == "zscore":
        z_scores = np.abs(stats.zscore(data_array))
        outlier_indices = np.where(z_scores > threshold)[0].tolist()
        outliers = data_array[outlier_indices].tolist()

    elif method == "modified_zscore":
        median = np.median(data_array)
        mad = np.median(np.abs(data_array - median))
        modified_z_scores = 0.6745 * (data_array - median) / mad
        outlier_indices = np.where(np.abs(modified_z_scores) > threshold)[0].tolist()
        outliers = data_array[outlier_indices].tolist()

    else:
        raise ValueError(f"Unsupported outlier detection method: {method}")

    return {
        "outliers": outliers,
        "outlier_indices": outlier_indices,
        "method": method,
        "threshold": threshold,
        "total_values": len(data),
        "outlier_percentage": (len(outliers) / len(data) * 100) if data else 0,
    }


def calculate_data_integrity_score(data: Dict[str, Any], data_type: str = "fastq") -> Dict[str, Any]:
    """Calculate data integrity score based on various checks.

    Args:
        data: Quality metrics data
        data_type: Type of biological data

    Returns:
        Dictionary with integrity score and issues found
    """
    integrity_checks = []
    total_checks = 0
    passed_checks = 0

    if data_type == "fastq":
        # Check for basic FASTQ integrity
        total_checks += 1
        if "basic_statistics" in data and data["basic_statistics"].get("total_reads", 0) > 0:
            passed_checks += 1
            integrity_checks.append({
                "check": "has_reads",
                "passed": True,
                "message": "FASTQ file contains reads",
            })
        else:
            integrity_checks.append({
                "check": "has_reads",
                "passed": False,
                "message": "FASTQ file contains no reads",
            })

        # Check for reasonable quality scores
        total_checks += 1
        if "basic_statistics" in data:
            mean_qual = data["basic_statistics"].get("mean_quality", 0)
            if 10 <= mean_qual <= 50:  # Reasonable quality range
                passed_checks += 1
                integrity_checks.append({
                    "check": "quality_range",
                    "passed": True,
                    "message": f"Mean quality ({mean_qual:.1f}) is within reasonable range",
                })
            else:
                integrity_checks.append({
                    "check": "quality_range",
                    "passed": False,
                    "message": f"Mean quality ({mean_qual:.1f}) is outside reasonable range",
                })

        # Check for sequence length consistency
        total_checks += 1
        if "sequence_length_distribution" in data:
            lengths = [item["length"] for item in data["sequence_length_distribution"].get("distribution", [])]
            if lengths and len(set(lengths)) <= 3:  # Allow some variation
                passed_checks += 1
                integrity_checks.append({
                    "check": "length_consistency",
                    "passed": True,
                    "message": "Sequence lengths are reasonably consistent",
                })
            else:
                integrity_checks.append({
                    "check": "length_consistency",
                    "passed": False,
                    "message": "Sequence lengths vary significantly",
                })

    integrity_score = (passed_checks / total_checks * 100) if total_checks > 0 else 0

    return {
        "integrity_score": integrity_score,
        "passed_checks": passed_checks,
        "total_checks": total_checks,
        "checks": integrity_checks,
    }


def compare_quality_metrics(dataset1: Dict[str, Any], dataset2: Dict[str, Any],
                          data_type: str = "fastq") -> Dict[str, Any]:
    """Compare quality metrics between two datasets.

    Args:
        dataset1: Quality metrics for first dataset
        dataset2: Quality metrics for second dataset
        data_type: Type of biological data

    Returns:
        Dictionary with comparison results
    """
    comparison = {
        "dataset1_score": calculate_quality_score(dataset1, data_type)["overall_score"],
        "dataset2_score": calculate_quality_score(dataset2, data_type)["overall_score"],
        "differences": {},
    }

    # Compare basic statistics
    if "basic_statistics" in dataset1 and "basic_statistics" in dataset2:
        stats1 = dataset1["basic_statistics"]
        stats2 = dataset2["basic_statistics"]

        for key in ["mean_quality", "mean_length", "mean_gc"]:
            if key in stats1 and key in stats2:
                diff = stats2[key] - stats1[key]
                comparison["differences"][key] = {
                    "dataset1": stats1[key],
                    "dataset2": stats2[key],
                    "difference": diff,
                    "percent_change": (diff / stats1[key] * 100) if stats1[key] != 0 else 0,
                }

    return comparison


def generate_quality_report(quality_data: Dict[str, Any], data_type: str = "fastq",
                          output_path: Optional[str | Path] = None) -> str:
    """Generate a comprehensive quality control report.

    Args:
        quality_data: Quality metrics data
        data_type: Type of biological data
        output_path: Optional path to save the report

    Returns:
        Formatted quality report as string
    """
    report_lines = []
    report_lines.append("=" * 60)
    report_lines.append(f"QUALITY CONTROL REPORT - {data_type.upper()}")
    report_lines.append("=" * 60)
    report_lines.append("")

    # Overall quality score
    quality_score = calculate_quality_score(quality_data, data_type)
    report_lines.append(f"Overall Quality Score: {quality_score['overall_score']:.1f}/100")
    report_lines.append(f"Grade: {quality_score['grade']}")
    report_lines.append("")

    # Data integrity
    integrity = calculate_data_integrity_score(quality_data, data_type)
    report_lines.append(f"Data Integrity Score: {integrity['integrity_score']:.1f}/100")
    report_lines.append(f"Passed Checks: {integrity['passed_checks']}/{integrity['total_checks']}")
    report_lines.append("")

    # Component breakdown
    if "components" in quality_score:
        report_lines.append("Quality Components:")
        for name, component in quality_score["components"].items():
            score = component["score"]
            weight = component["weight"] * 100
            details = component["details"]
            report_lines.append(f"  {name}: {score:.1f}/{weight:.0f}% - {details}")
        report_lines.append("")

    # Key metrics
    if "basic_statistics" in quality_data:
        stats = quality_data["basic_statistics"]
        report_lines.append("Basic Statistics:")
        report_lines.append(f"  Total Reads: {stats.get('total_reads', 'N/A'):,}")
        report_lines.append(f"  Mean Quality: {stats.get('mean_quality', 'N/A'):.1f}")
        report_lines.append(f"  Mean Length: {stats.get('mean_length', 'N/A'):.1f}")
        report_lines.append(f"  Mean GC Content: {stats.get('mean_gc', 'N/A'):.1f}%")
        report_lines.append("")

    # Recommendations
    report_lines.append("Recommendations:")
    if quality_score["overall_score"] >= 90:
        report_lines.append("  ✓ Data quality is excellent. Proceed with analysis.")
    elif quality_score["overall_score"] >= 80:
        report_lines.append("  ✓ Data quality is good. Minor improvements may be beneficial.")
    elif quality_score["overall_score"] >= 70:
        report_lines.append("  ⚠ Data quality is acceptable but could be improved.")
    elif quality_score["overall_score"] >= 60:
        report_lines.append("  ⚠ Data quality needs improvement before analysis.")
    else:
        report_lines.append("  ✗ Data quality is poor. Significant issues need to be addressed.")

    report = "\n".join(report_lines)

    if output_path:
        output_path = Path(output_path)
        with open(output_path, 'w') as f:
            f.write(report)
        logger.info(f"Quality report saved to {output_path}")

    return report


def batch_quality_analysis(file_paths: List[str | Path], data_type: str = "fastq",
                          n_reads: Optional[int] = None) -> Dict[str, Any]:
    """Perform quality analysis on multiple files.

    Args:
        file_paths: List of file paths to analyze
        data_type: Type of biological data
        n_reads: Number of reads to sample per file

    Returns:
        Dictionary with batch analysis results
    """
    results = {}

    for file_path in file_paths:
        try:
            logger.info(f"Analyzing {file_path}")
            if data_type == "fastq":
                from .fastq import analyze_fastq_quality
                quality_data = analyze_fastq_quality(file_path, n_reads)
            else:
                # Placeholder for other data types
                quality_data = {}

            quality_score = calculate_quality_score(quality_data, data_type)

            results[str(file_path)] = {
                "quality_data": quality_data,
                "quality_score": quality_score,
                "status": "success",
            }

        except Exception as e:
            logger.error(f"Failed to analyze {file_path}: {e}")
            results[str(file_path)] = {
                "error": str(e),
                "status": "failed",
            }

    # Summary statistics
    successful_analyses = [r for r in results.values() if r["status"] == "success"]
    if successful_analyses:
        scores = [r["quality_score"]["overall_score"] for r in successful_analyses]
        summary = {
            "total_files": len(file_paths),
            "successful_analyses": len(successful_analyses),
            "failed_analyses": len(file_paths) - len(successful_analyses),
            "mean_score": statistics.mean(scores),
            "median_score": statistics.median(scores),
            "min_score": min(scores),
            "max_score": max(scores),
        }
    else:
        summary = {
            "total_files": len(file_paths),
            "successful_analyses": 0,
            "failed_analyses": len(file_paths),
        }

    return {
        "results": results,
        "summary": summary,
    }


def calculate_complexity_metrics(sequences: List[str]) -> Dict[str, Any]:
    """Calculate sequence complexity metrics.

    Args:
        sequences: List of DNA/RNA/protein sequences

    Returns:
        Dictionary with complexity metrics
    """
    if not sequences:
        return {"error": "No sequences provided"}

    metrics = {
        "total_sequences": len(sequences),
        "average_length": sum(len(seq) for seq in sequences) / len(sequences),
    }

    # Linguistic complexity (Shannon entropy of k-mers)
    k = 2  # Use dinucleotides for complexity
    all_kmers = []
    for seq in sequences:
        if len(seq) >= k:
            kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
            all_kmers.extend(kmers)

    if all_kmers:
        from collections import Counter
        kmer_counts = Counter(all_kmers)
        total_kmers = len(all_kmers)

        # Shannon entropy
        entropy = 0
        for count in kmer_counts.values():
            p = count / total_kmers
            entropy -= p * math.log2(p)

        # Normalized entropy (0-1 scale)
        max_entropy = math.log2(len(kmer_counts))
        normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0

        metrics["kmer_entropy"] = entropy
        metrics["normalized_kmer_entropy"] = normalized_entropy
        metrics["unique_kmers"] = len(kmer_counts)
        metrics["total_kmers"] = total_kmers

        # Complexity score (combination of uniqueness and entropy)
        uniqueness = len(kmer_counts) / total_kmers
        metrics["complexity_score"] = (uniqueness + normalized_entropy) / 2

    # Sequence diversity (proportion of unique sequences)
    unique_sequences = len(set(sequences))
    metrics["sequence_diversity"] = unique_sequences / len(sequences)

    # GC content complexity (deviation from 0.5)
    gc_contents = []
    for seq in sequences:
        gc_count = seq.upper().count('G') + seq.upper().count('C')
        gc_content = gc_count / len(seq) if seq else 0
        gc_contents.append(gc_content)

    if gc_contents:
        avg_gc = sum(gc_contents) / len(gc_contents)
        gc_complexity = 1 - abs(avg_gc - 0.5) * 2  # Higher when closer to 0.5
        metrics["gc_complexity"] = gc_complexity

    return metrics

