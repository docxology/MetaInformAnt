"""High-level information analysis workflows for METAINFORMANT.

This module provides comprehensive workflows that combine multiple
information-theoretic analyses into end-to-end pipelines for
biological sequence and data analysis.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np

from metainformant.core import io, logging
from metainformant.information.metrics import analysis, continuous, estimation, semantic, syntactic

logger = logging.get_logger(__name__)


def batch_entropy_analysis(
    sequences: List[str],
    k: int = 1,
    output_dir: Optional[Union[str, Path]] = None,
    methods: Optional[List[str]] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Batch entropy analysis for multiple sequences.

    Args:
        sequences: List of biological sequences
        k: k-mer size for analysis
        output_dir: Directory to save results
        methods: List of entropy estimation methods
        **kwargs: Additional parameters

    Returns:
        Dictionary with batch analysis results
    """
    if methods is None:
        methods = ["plugin", "miller_madow"]

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    results = {
        "batch_info": {
            "n_sequences": len(sequences),
            "sequence_lengths": [len(seq) for seq in sequences],
            "k": k,
            "methods": methods,
        },
        "sequence_results": [],
        "summary": {},
    }

    start_time = time.time()

    for i, sequence in enumerate(sequences):
        seq_result = {
            "sequence_index": i,
            "sequence_length": len(sequence),
            "methods": {},
        }

        # Analyze with each method
        for method in methods:
            try:
                # Calculate entropy for k-mers
                kmers = [sequence[j : j + k] for j in range(len(sequence) - k + 1)]
                from collections import Counter

                kmer_counts = Counter(kmers)

                entropy_val = estimation.entropy_estimator(kmer_counts, method=method, **kwargs)

                seq_result["methods"][method] = {
                    "entropy": entropy_val,
                    "n_unique_kmers": len(kmer_counts),
                    "n_total_kmers": len(kmers),
                }

            except Exception as e:
                logger.warning(f"Entropy calculation failed for sequence {i}, method {method}: {e}")
                seq_result["methods"][method] = {"error": str(e)}

        results["sequence_results"].append(seq_result)

    # Calculate summary statistics
    all_entropies = {}
    for method in methods:
        method_entropies = [
            seq_res["methods"].get(method, {}).get("entropy", None) for seq_res in results["sequence_results"]
        ]
        method_entropies = [e for e in method_entropies if e is not None]

        if method_entropies:
            all_entropies[method] = {
                "mean": float(np.mean(method_entropies)),
                "std": float(np.std(method_entropies)),
                "min": float(np.min(method_entropies)),
                "max": float(np.max(method_entropies)),
                "n_successful": len(method_entropies),
            }

    results["summary"] = {
        "entropy_statistics": all_entropies,
        "processing_time": time.time() - start_time,
        "success_rate": len([r for r in results["sequence_results"] if all(m in r["methods"] for m in methods)])
        / len(sequences),
    }

    # Save results if output directory specified
    if output_dir:
        output_file = output_dir / "batch_entropy_analysis.json"
        io.dump_json(results, output_file)
        logger.info(f"Saved batch entropy analysis to {output_file}")

    logger.info(
        f"Completed batch entropy analysis: {len(sequences)} sequences, "
        f"{len(methods)} methods, {results['summary']['success_rate']:.1%} success rate"
    )

    return results


def information_workflow(
    sequences: List[str],
    k_values: Optional[List[int]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    include_profiles: bool = True,
    include_signatures: bool = False,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Comprehensive information workflow for sequence analysis.

    Args:
        sequences: List of biological sequences
        k_values: List of k-mer sizes to analyze
        output_dir: Directory to save results
        include_profiles: Whether to calculate information profiles
        include_signatures: Whether to calculate information signatures
        **kwargs: Additional parameters

    Returns:
        Dictionary with comprehensive workflow results
    """
    if k_values is None:
        k_values = [1, 2, 3]

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    workflow_results = {
        "workflow_info": {
            "n_sequences": len(sequences),
            "k_values": k_values,
            "include_profiles": include_profiles,
            "include_signatures": include_signatures,
        },
        "sequence_analyses": [],
        "aggregate_results": {},
    }

    start_time = time.time()

    # Analyze each sequence
    for i, sequence in enumerate(sequences):
        seq_analysis = analysis.analyze_sequence_information(sequence, k_values=k_values, **kwargs)
        seq_analysis["sequence_index"] = i

        workflow_results["sequence_analyses"].append(seq_analysis)

    # Aggregate results
    if len(sequences) > 1:
        workflow_results["aggregate_results"] = _aggregate_sequence_analyses(workflow_results["sequence_analyses"])

    # Calculate information profiles if requested
    if include_profiles and len(sequences) > 1:
        try:
            # Assume sequences are aligned for profile calculation
            profile_results = {}
            for k in k_values:
                if all(len(seq) >= k for seq in sequences):
                    profile = analysis.information_profile(sequences, k=k, **kwargs)
                    profile_results[f"k{k}"] = profile

            workflow_results["information_profiles"] = profile_results

        except Exception as e:
            logger.warning(f"Information profile calculation failed: {e}")
            workflow_results["information_profiles"] = {"error": str(e)}

    # Calculate information signatures if requested
    if include_signatures:
        try:
            # Convert sequences to numerical representation for signature analysis
            numerical_sequences = []
            for seq in sequences:
                # Simple numerical encoding (can be improved)
                char_to_num = {char: i for i, char in enumerate(set(seq))}
                numerical = [char_to_num[char] for char in seq]
                numerical_sequences.append(numerical)

            signature = analysis.information_signature(np.array(numerical_sequences).T, **kwargs)
            workflow_results["information_signature"] = signature

        except Exception as e:
            logger.warning(f"Information signature calculation failed: {e}")
            workflow_results["information_signature"] = {"error": str(e)}

    workflow_results["processing_time"] = time.time() - start_time
    workflow_results["workflow_status"] = "completed"

    # Save results
    if output_dir:
        output_file = output_dir / "information_workflow_results.json"
        io.dump_json(workflow_results, output_file)
        logger.info(f"Saved information workflow results to {output_file}")

    logger.info(
        f"Completed information workflow: {len(sequences)} sequences, "
        f"{len(k_values)} k-values, {workflow_results['processing_time']:.2f}s"
    )

    return workflow_results


def _aggregate_sequence_analyses(analyses: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Aggregate results from multiple sequence analyses."""
    if not analyses:
        return {}

    aggregate = {
        "sequence_types": {},
        "k_mer_statistics": {},
        "entropy_statistics": {},
        "complexity_statistics": {},
    }

    # Aggregate sequence types
    seq_types = [analysis.get("sequence_type", "unknown") for analysis in analyses]
    from collections import Counter

    aggregate["sequence_types"] = dict(Counter(seq_types))

    # Aggregate k-mer statistics
    all_k_stats = {}
    for analysis in analyses:
        for k_key, k_stats in analysis.get("k_mer_analysis", {}).items():
            if k_key not in all_k_stats:
                all_k_stats[k_key] = []
            all_k_stats[k_key].append(k_stats)

    for k_key, stats_list in all_k_stats.items():
        if stats_list:
            entropies = [s.get("entropy", 0) for s in stats_list if s.get("entropy") is not None]
            if entropies:
                aggregate["k_mer_statistics"][k_key] = {
                    "mean_entropy": float(np.mean(entropies)),
                    "std_entropy": float(np.std(entropies)),
                    "min_entropy": float(np.min(entropies)),
                    "max_entropy": float(np.max(entropies)),
                    "n_sequences": len(entropies),
                }

    # Aggregate entropy statistics
    entropies = []
    for analysis in analyses:
        entropy_data = analysis.get("methods", {}).get("entropy", {})
        if "shannon_entropy" in entropy_data:
            entropies.append(entropy_data["shannon_entropy"])

    if entropies:
        aggregate["entropy_statistics"] = {
            "mean": float(np.mean(entropies)),
            "std": float(np.std(entropies)),
            "min": float(np.min(entropies)),
            "max": float(np.max(entropies)),
            "n_sequences": len(entropies),
        }

    return aggregate


def compare_datasets(
    dataset1: List[str],
    dataset2: List[str],
    k: int = 1,
    output_dir: Optional[Union[str, Path]] = None,
    comparison_methods: Optional[List[str]] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Compare information content between two datasets.

    Args:
        dataset1: First dataset (list of sequences)
        dataset2: Second dataset (list of sequences)
        k: k-mer size for analysis
        output_dir: Directory to save results
        comparison_methods: List of comparison methods
        **kwargs: Additional parameters

    Returns:
        Dictionary with dataset comparison results
    """
    if comparison_methods is None:
        comparison_methods = ["entropy", "diversity", "similarity"]

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    comparison_results = {
        "dataset_info": {
            "dataset1_size": len(dataset1),
            "dataset2_size": len(dataset2),
            "k": k,
            "comparison_methods": comparison_methods,
        },
        "comparisons": {},
    }

    start_time = time.time()

    # Analyze each dataset
    dataset1_analysis = batch_entropy_analysis(dataset1, k=k, **kwargs)
    dataset2_analysis = batch_entropy_analysis(dataset2, k=k, **kwargs)

    comparison_results["dataset1_analysis"] = dataset1_analysis
    comparison_results["dataset2_analysis"] = dataset2_analysis

    # Perform comparisons
    for method in comparison_methods:
        try:
            if method == "entropy":
                comparison_results["comparisons"]["entropy"] = _compare_entropy_distributions(
                    dataset1_analysis, dataset2_analysis
                )

            elif method == "diversity":
                comparison_results["comparisons"]["diversity"] = _compare_sequence_diversity(dataset1, dataset2, k=k)

            elif method == "similarity":
                comparison_results["comparisons"]["similarity"] = _compare_sequence_similarity(
                    dataset1[: min(10, len(dataset1))], dataset2[: min(10, len(dataset2))], k=k  # Sample for efficiency
                )

        except Exception as e:
            logger.warning(f"Comparison method {method} failed: {e}")
            comparison_results["comparisons"][method] = {"error": str(e)}

    comparison_results["processing_time"] = time.time() - start_time

    # Save results
    if output_dir:
        output_file = output_dir / "dataset_comparison_results.json"
        io.dump_json(comparison_results, output_file)
        logger.info(f"Saved dataset comparison results to {output_file}")

    logger.info(
        f"Completed dataset comparison: {len(dataset1)} vs {len(dataset2)} sequences, "
        f"{len(comparison_methods)} methods"
    )

    return comparison_results


def _compare_entropy_distributions(analysis1: Dict, analysis2: Dict) -> Dict[str, Any]:
    """Compare entropy distributions between two datasets."""
    entropies1 = []
    entropies2 = []

    # Extract entropy values
    for seq_result in analysis1.get("sequence_results", []):
        for method_results in seq_result.get("methods", {}).values():
            if "entropy" in method_results:
                entropies1.append(method_results["entropy"])

    for seq_result in analysis2.get("sequence_results", []):
        for method_results in seq_result.get("methods", {}).values():
            if "entropy" in method_results:
                entropies2.append(method_results["entropy"])

    if not entropies1 or not entropies2:
        return {"error": "No entropy values found"}

    # Statistical comparison
    from scipy import stats

    comparison = {
        "dataset1": {
            "n_values": len(entropies1),
            "mean": float(np.mean(entropies1)),
            "std": float(np.std(entropies1)),
        },
        "dataset2": {
            "n_values": len(entropies2),
            "mean": float(np.mean(entropies2)),
            "std": float(np.std(entropies2)),
        },
        "difference": {
            "mean_diff": float(np.mean(entropies1) - np.mean(entropies2)),
            "std_diff": float(np.sqrt(np.var(entropies1) / len(entropies1) + np.var(entropies2) / len(entropies2))),
        },
    }

    # Statistical test
    try:
        t_stat, p_value = stats.ttest_ind(entropies1, entropies2)
        comparison["statistical_test"] = {
            "test": "t-test",
            "t_statistic": float(t_stat),
            "p_value": float(p_value),
            "significant": p_value < 0.05,
        }
    except Exception as e:
        comparison["statistical_test"] = {"error": str(e)}

    return comparison


def _compare_sequence_diversity(seqs1: List[str], seqs2: List[str], k: int = 1) -> Dict[str, Any]:
    """Compare sequence diversity between two datasets."""
    from collections import Counter

    # Calculate k-mer diversity for each dataset
    kmers1 = []
    for seq in seqs1:
        kmers1.extend([seq[i : i + k] for i in range(len(seq) - k + 1)])

    kmers2 = []
    for seq in seqs2:
        kmers2.extend([seq[i : i + k] for i in range(len(seq) - k + 1)])

    diversity1 = len(set(kmers1)) / len(kmers1) if kmers1 else 0
    diversity2 = len(set(kmers2)) / len(kmers2) if kmers2 else 0

    # Calculate shared k-mers
    set1 = set(kmers1)
    set2 = set(kmers2)
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    jaccard = intersection / union if union > 0 else 0

    return {
        "diversity1": diversity1,
        "diversity2": diversity2,
        "diversity_ratio": diversity1 / diversity2 if diversity2 > 0 else float("inf"),
        "jaccard_similarity": jaccard,
        "shared_kmers": intersection,
        "unique_kmers_1": len(set1 - set2),
        "unique_kmers_2": len(set2 - set1),
    }


def _compare_sequence_similarity(seqs1: List[str], seqs2: List[str], k: int = 1) -> Dict[str, Any]:
    """Compare sequence similarity between datasets."""
    similarities = []

    # Calculate pairwise similarities (sample for efficiency)
    for seq1 in seqs1[:5]:  # Limit for computational efficiency
        for seq2 in seqs2[:5]:
            if len(seq1) == len(seq2):  # Only compare same-length sequences
                sim = analysis.compare_sequences_information(seq1, seq2, k=k, method="mutual_info")
                similarities.append(sim.get("mean_mi", 0))

    return {
        "mean_similarity": float(np.mean(similarities)) if similarities else 0,
        "std_similarity": float(np.std(similarities)) if similarities else 0,
        "n_comparisons": len(similarities),
        "similarity_range": [
            float(np.min(similarities)) if similarities else 0,
            float(np.max(similarities)) if similarities else 0,
        ],
    }


def information_report(
    results: Dict[str, Any], output_path: Optional[Union[str, Path]] = None, format: str = "markdown"
) -> None:
    """Generate a comprehensive report from information analysis results.

    Args:
        results: Results dictionary from information analysis
        output_path: Path to save report
        format: Report format ('markdown', 'json', 'text')
    """
    if format == "markdown":
        report = _generate_markdown_report(results)
    elif format == "json":
        report = results  # Already JSON-serializable
    elif format == "text":
        report = _generate_text_report(results)
    else:
        raise ValueError(f"Unknown format: {format}")

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if format == "json":
            io.dump_json(report, output_path)
        else:
            with open(output_path, "w") as f:
                f.write(report)

        logger.info(f"Generated information report: {output_path}")

    else:
        # Print to console
        if format == "json":
            print(io.dumps_json(report))
        else:
            print(report)


def _generate_markdown_report(results: Dict[str, Any]) -> str:
    """Generate markdown report from results."""
    lines = ["# Information Analysis Report\n"]

    # Summary section
    if "workflow_info" in results:
        info = results["workflow_info"]
        lines.append("## Summary\n")
        lines.append(f"- **Sequences analyzed**: {info.get('n_sequences', 'N/A')}")
        lines.append(f"- **K-mer sizes**: {info.get('k_values', 'N/A')}")
        lines.append(f"- **Processing time**: {results.get('processing_time', 'N/A'):.2f}s")
        lines.append("")

    # Results section
    if "aggregate_results" in results:
        agg = results["aggregate_results"]
        lines.append("## Aggregate Results\n")

        if "entropy_statistics" in agg:
            lines.append("### Entropy Statistics\n")
            ent_stats = agg["entropy_statistics"]
            if isinstance(ent_stats, dict) and "mean" in ent_stats:
                lines.append(
                    f"- **Shannon entropy**: {ent_stats['mean']:.3f} ± {ent_stats.get('std', 0):.3f}"
                    f" (range: {ent_stats.get('min', 0):.3f} – {ent_stats.get('max', 0):.3f},"
                    f" n={ent_stats.get('n_sequences', 'N/A')})"
                )
            lines.append("")

    # Detailed results
    if "sequence_analyses" in results and results["sequence_analyses"]:
        lines.append("## Sequence Details\n")
        lines.append("| Sequence | Length | Shannon Entropy | Type |")
        lines.append("|----------|--------|----------------|------|")

        for analysis in results["sequence_analyses"][:10]:  # Limit for readability
            seq_idx = analysis.get("sequence_index", "N/A")
            length = analysis.get("sequence_length", "N/A")
            entropy = analysis.get("methods", {}).get("entropy", {}).get("shannon_entropy", "N/A")
            seq_type = analysis.get("sequence_type", "N/A")

            if isinstance(entropy, (int, float)):
                entropy = f"{entropy:.3f}"

            lines.append(f"| {seq_idx} | {length} | {entropy} | {seq_type} |")

        if len(results["sequence_analyses"]) > 10:
            lines.append(f"\n... and {len(results['sequence_analyses']) - 10} more sequences")
        lines.append("")

    return "\n".join(lines)


def _generate_text_report(results: Dict[str, Any]) -> str:
    """Generate plain text report from results."""
    lines = ["INFORMATION ANALYSIS REPORT\n"]
    lines.append("=" * 50)

    # Summary
    if "workflow_info" in results:
        info = results["workflow_info"]
        lines.append(f"Sequences analyzed: {info.get('n_sequences', 'N/A')}")
        lines.append(f"K-mer sizes: {info.get('k_values', 'N/A')}")
        lines.append(f"Processing time: {results.get('processing_time', 'N/A'):.2f}s")
        lines.append("")

    # Results summary
    if "aggregate_results" in results:
        agg = results["aggregate_results"]
        lines.append("AGGREGATE RESULTS")
        lines.append("-" * 20)

        if "entropy_statistics" in agg:
            lines.append("Entropy Statistics:")
            for method, stats in agg["entropy_statistics"].items():
                lines.append(f"  {method}: {stats['mean']:.3f} ± {stats['std']:.3f}")
            lines.append("")

    return "\n".join(lines)
