"""Information theory workflow functions.

This module provides high-level workflow functions for batch processing
and comprehensive information-theoretic analysis.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from metainformant.core import io, paths
from metainformant.information.analysis import (
    analyze_sequence_information,
    compare_sequences_information,
    information_profile,
)
from metainformant.information.syntactic import (
    mutual_information,
    shannon_entropy_from_counts,
)


def batch_entropy_analysis(
    sequences: list[str],
    k: int = 1,
    output_dir: Path | None = None
) -> dict[str, Any]:
    """Perform entropy analysis on multiple sequences.
    
    Args:
        sequences: List of sequences to analyze
        k: K-mer size for analysis
        output_dir: Optional output directory for saving results
        
    Returns:
        Dictionary with analysis results for each sequence
    """
    results: dict[str, Any] = {
        "sequences": [],
        "summary": {},
    }
    
    entropies = []
    for i, seq in enumerate(sequences):
        analysis = analyze_sequence_information(seq, k_values=[k])
        if k in analysis["kmer_analyses"]:
            entropy = analysis["kmer_analyses"][k]["entropy"]
            entropies.append(entropy)
            results["sequences"].append({
                "index": i,
                "length": len(seq),
                "entropy": entropy,
                "unique_kmers": analysis["kmer_analyses"][k]["unique_kmers"],
            })
    
    # Summary statistics
    if entropies:
        results["summary"] = {
            "mean_entropy": float(np.mean(entropies)),
            "std_entropy": float(np.std(entropies)),
            "min_entropy": float(np.min(entropies)),
            "max_entropy": float(np.max(entropies)),
            "total_sequences": len(sequences),
        }
    
    # Save results if output directory provided
    if output_dir:
        output_dir = Path(output_dir)
        paths.ensure_directory(output_dir)
        output_file = output_dir / "batch_entropy_analysis.json"
        io.dump_json(results, output_file)
    
    return results


def information_workflow(
    sequences: list[str],
    k_values: list[int] | None = None,
    output_dir: Path | None = None
) -> dict[str, Any]:
    """Complete information-theoretic workflow for sequence analysis.
    
    Args:
        sequences: List of sequences to analyze
        k_values: List of k-mer sizes to analyze
        output_dir: Optional output directory for saving results
        
    Returns:
        Comprehensive workflow results
    """
    if k_values is None:
        k_values = [1, 2, 3]
    
    results: dict[str, Any] = {
        "profiles": {},
        "comparisons": {},
        "summary": {},
    }
    
    # Information profiles for each k
    for k in k_values:
        profile = information_profile(sequences, k=k)
        results["profiles"][k] = profile
    
    # Pairwise comparisons
    if len(sequences) >= 2:
        comparisons = []
        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                comparison = compare_sequences_information(
                    sequences[i], sequences[j], k=k_values[0]
                )
                comparisons.append({
                    "seq1_index": i,
                    "seq2_index": j,
                    "comparison": comparison,
                })
        results["comparisons"] = comparisons
    
    # Summary
    if results["profiles"]:
        first_k = k_values[0]
        if first_k in results["profiles"]:
            profile = results["profiles"][first_k]
            results["summary"] = {
                "total_sequences": len(sequences),
                "k_values": k_values,
                "mean_entropy": profile.get("entropy", 0.0),
                "unique_kmers": profile.get("unique_kmers", 0),
            }
    
    # Save results if output directory provided
    if output_dir:
        output_dir = Path(output_dir)
        paths.ensure_directory(output_dir)
        output_file = output_dir / "information_workflow_results.json"
        io.dump_json(results, output_file)
    
    return results


def compare_datasets(
    dataset1: list[str] | dict[str, Any],
    dataset2: list[str] | dict[str, Any],
    k: int = 1,
    method: str = "entropy"
) -> dict[str, Any]:
    """Compare two datasets using information measures.
    
    Args:
        dataset1: First dataset (sequences or data dictionary)
        dataset2: Second dataset (sequences or data dictionary)
        k: K-mer size for sequence analysis
        method: Comparison method ("entropy", "mutual_information", "kl_divergence")
        
    Returns:
        Comparison results
    """
    # Convert to sequences if needed
    if isinstance(dataset1, dict):
        seqs1 = dataset1.get("sequences", [])
    else:
        seqs1 = dataset1
    
    if isinstance(dataset2, dict):
        seqs2 = dataset2.get("sequences", [])
    else:
        seq2 = dataset2
    
    results: dict[str, Any] = {
        "method": method,
        "k": k,
    }
    
    if method == "entropy":
        # Compare entropy distributions
        profile1 = information_profile(seqs1, k=k)
        profile2 = information_profile(seqs2, k=k)
        
        results["dataset1_entropy"] = profile1["entropy"]
        results["dataset2_entropy"] = profile2["entropy"]
        results["entropy_difference"] = abs(profile1["entropy"] - profile2["entropy"])
    
    elif method == "mutual_information":
        # Calculate MI between datasets
        # Flatten sequences into single strings for comparison
        flat1 = "".join(seqs1)
        flat2 = "".join(seqs2)
        
        # Sample k-mers from each
        kmers1 = [flat1[i:i+k] for i in range(len(flat1) - k + 1)]
        kmers2 = [flat2[i:i+k] for i in range(min(len(flat2) - k + 1, len(kmers1)))]
        
        if len(kmers1) == len(kmers2):
            mi = mutual_information(kmers1, kmers2)
            results["mutual_information"] = mi
    
    elif method == "kl_divergence":
        # Calculate KL divergence between k-mer distributions
        from metainformant.information.syntactic import kl_divergence
        
        profile1 = information_profile(seqs1, k=k)
        profile2 = information_profile(seqs2, k=k)
        
        # Get all k-mers
        all_kmers = set(profile1["kmer_frequencies"].keys()) | set(
            profile2["kmer_frequencies"].keys()
        )
        
        # Create probability distributions
        total1 = sum(profile1["kmer_frequencies"].values())
        total2 = sum(profile2["kmer_frequencies"].values())
        
        p1 = [
            profile1["kmer_frequencies"].get(kmer, 0) / total1 if total1 > 0 else 0
            for kmer in sorted(all_kmers)
        ]
        p2 = [
            profile2["kmer_frequencies"].get(kmer, 0) / total2 if total2 > 0 else 0
            for kmer in sorted(all_kmers)
        ]
        
        kl = kl_divergence(p1, p2)
        results["kl_divergence"] = kl
    
    return results


def information_report(
    data: dict[str, Any],
    output_path: Path | str,
    format: str = "json"
) -> dict[str, Any]:
    """Generate comprehensive information-theoretic report.
    
    Args:
        data: Analysis results dictionary
        output_path: Path to save report
        format: Report format ("json", "markdown", "text")
        
    Returns:
        Report metadata
    """
    output_path = Path(output_path)
    paths.ensure_directory(output_path.parent)
    
    if format == "json":
        io.dump_json(data, output_path)
        return {"output_path": str(output_path), "format": format}
    
    elif format == "markdown":
        # Generate markdown report
        markdown_lines = ["# Information Theory Analysis Report\n"]
        
        if "summary" in data:
            markdown_lines.append("## Summary\n")
            for key, value in data["summary"].items():
                markdown_lines.append(f"- **{key}**: {value}\n")
        
        if "profiles" in data:
            markdown_lines.append("\n## Information Profiles\n")
            for k, profile in data["profiles"].items():
                markdown_lines.append(f"### K-mer size {k}\n")
                markdown_lines.append(f"- Entropy: {profile.get('entropy', 0):.4f} bits\n")
                markdown_lines.append(f"- Unique k-mers: {profile.get('unique_kmers', 0)}\n")
        
        if "comparisons" in data:
            markdown_lines.append("\n## Sequence Comparisons\n")
            for comp in data["comparisons"]:
                markdown_lines.append(
                    f"### Sequence {comp['seq1_index']} vs {comp['seq2_index']}\n"
                )
                comp_data = comp.get("comparison", {})
                markdown_lines.append(f"- KL Divergence: {comp_data.get('kl_divergence', 0):.4f}\n")
                markdown_lines.append(
                    f"- Mutual Information: {comp_data.get('mutual_information', 0):.4f}\n"
                )
        
        output_path.write_text("".join(markdown_lines))
        return {"output_path": str(output_path), "format": format}
    
    elif format == "text":
        # Generate text report
        text_lines = ["Information Theory Analysis Report\n"]
        text_lines.append("=" * 50 + "\n\n")
        
        if "summary" in data:
            text_lines.append("Summary:\n")
            for key, value in data["summary"].items():
                text_lines.append(f"  {key}: {value}\n")
        
        output_path.write_text("".join(text_lines))
        return {"output_path": str(output_path), "format": format}
    
    else:
        raise ValueError(f"Unknown format: {format}")

